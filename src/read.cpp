#include "read.h"
#include "constants.h"
#include "functions.h"

#include <mkl_scalapack.h>   // numroc_, descinit_
#include <mkl_blacs.h>       // blacs_gridinfo, blacs_gridexit_
#include <algorithm>         // std::max
#include <vector>
#include <istream>
#include <fstream>
#include <iostream>          // For std::cout, std::endl
#include <cstdlib>           // calloc, exit
#include <cstring>           // for 'strcpy', 'strtok'
#include <cmath>             // for std::fabs



int read_individual_ids(
		MKL_INT ictxt,
		const std::string& id_filename,
		std::vector<std::string>& ordered_ids
		) {
	std::ifstream id_file(id_filename);
	if (!id_file) {
		std::cerr << "ERROR: Cannot open individual ID file: " 
			<< id_filename << std::endl;
		blacs_gridexit_(&ictxt);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	ordered_ids.clear();
	std::string line;
	while (std::getline(id_file, line)) {
		// trim whitespace
		line.erase(0, line.find_first_not_of(" \t\n\r"));
		line.erase(line.find_last_not_of(" \t\n\r") + 1);
		if (line.empty()) continue;
		ordered_ids.push_back(line);
	}
	id_file.close();

	if (ordered_ids.empty()) {
		std::cerr << "ERROR: No valid individual IDs found in " 
			<< id_filename << std::endl;
		blacs_gridexit_(&ictxt);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	std::cout << "Read " << ordered_ids.size() 
		<< " individual IDs from " << id_filename << std::endl;
	return 0;
}


int read_phenotype_data(
		const Options &options,
		MKL_INT ictxt,
		const std::vector<std::string>& individual_ids,
		std::vector<int>& rec2ind,
		int &num_rec,
		std::vector<double>& phenotypes
		) {
	std::ifstream input_phen(options.pheno_file);
	if (!input_phen) {
		std::cerr << "ERROR: Cannot open phenotype file: " << options.pheno_file << std::endl;
		blacs_gridexit_(&ictxt);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	phenotypes.clear();
	num_rec = 0;
	std::string line;
	int pheno_column=options.trait_col - 1;
	int line_num = 0; 

	while (std::getline(input_phen, line)) {
		line_num++;
		// TODO: add csv support, comma separator
		std::vector<std::string> words = str_split(" \r\t\n\x0b", line);

		if (words.empty()) continue; // Skip empty lines

		// Basic check for enough columns
		// Ensure trait_col is within bounds (0-based)
		if (words.size() <= static_cast<size_t>(pheno_column)) {
			std::cerr << "WARNING: Skipping line " << line_num << " due to insufficient columns for trait." << std::endl;
			continue;
		}
		// Also check ID column (assuming index 0)
		if (words.size() < 1) {
			std::cerr << "WARNING: Skipping line " << line_num << " due to missing individual ID column." << std::endl;
			continue;
		}

		// Get individual ID (assuming it's the first column, index 0)
		const std::string& ind_id = words[0];

		// Print name of column
		if (ind_id == "id" || ind_id == "ID") {
			std::cout << "name of phenotype column " << options.trait_col << ": " << words[pheno_column] << std::endl;
			continue;
		}

		double y_val;
		try {
			y_val = std::stod(words[pheno_column]); // trait_col is already 0-based
		} catch (const std::invalid_argument& ia) {
			std::cerr << "WARNING: Invalid numeric value for trait on line " << line_num << ". Skipping." << std::endl;
			continue;
		} catch (const std::out_of_range& oor) {
			std::cerr << "WARNING: Trait value out of range on line " << line_num << ". Skipping." << std::endl;
			continue;
		}


		// Check for missing phenotype value
		// TODO: change to string check
		if (std::fabs(y_val - options.missing_phenotype) < 1e-9) { // Use tolerance for float comparison
			std::cout << "missing at: " << line_num << "," << y_val << std::endl;
			continue; // Skip record with missing phenotype
		}

		// Find the individual's row index by searching individual_ids[]
		auto it = std::find(individual_ids.begin(), individual_ids.end(), ind_id);
		if (it == individual_ids.end()) {
			// pheno file ID not in genotype (.id.txt list)
			std::cout << "id: "
				<< ind_id << ","
				<< "not in grm ids (genotyped)"
				<< std::endl;
			continue;
		}
		int ind_col_idx = int(std::distance(individual_ids.begin(), it));

		// Store the valid record information
		//std::cout << "valid record (line,ind,,num_rec,idx,val): "
		//	<< line_num << ","
		//	<< ind_id << ","
		//	<< num_rec << ","
		//	<< ind_col_idx << ","
		//	<< y_val << std::endl;
		rec2ind.push_back(ind_col_idx);
		phenotypes.push_back(y_val);
		num_rec++;  // Increment count of valid records
	}

	input_phen.close();
    std::cout << "num_rec: " << num_rec << std::endl;

	if (rec2ind.empty()) {
		std::cerr << "ERROR: No valid phenotype records found after filtering and matching IDs." << std::endl;
		blacs_gridexit_(&ictxt);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	return 0; 
}

void load_grms(
		const Options &options,
		MKL_INT ictxt,
		MKL_INT &mkl_num_ind,
		std::map<std::string,int> &var_idx,
		std::map<int,std::string> &idx_var,
		std::map<std::string,double*> &Glocal,
		MKL_INT descG[],
		std::map<std::string,double> &h2,
		MKL_INT mem_tot
		) {
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// Build the map of file extensions internally
	std::map<std::string,std::string> var_extensions;
	var_extensions["a"]  = ".g.A";
	var_extensions["d"]  = ".g.D";
	var_extensions["h"]  = ".g.AH";
	var_extensions["aa"] = ".g.AA";
	var_extensions["ad"] = ".g.AD";
	var_extensions["dd"] = ".g.DD";
	var_extensions["aa-inter"] = ".g.AA-inter";
	var_extensions["aa-intra"] = ".g.AA-intra";
	var_extensions["ad-inter"] = ".g.AD-inter";
	var_extensions["ad-intra"] = ".g.AD-intra";
	var_extensions["dd-inter"] = ".g.DD-inter";
	var_extensions["dd-intra"] = ".g.DD-intra";
	var_extensions["aaa"] = ".g.AAA";
	var_extensions["aad"] = ".g.AAD";
	var_extensions["add"] = ".g.ADD";
	var_extensions["ddd"] = ".g.DDD";
	var_extensions["1"]  = ".g.1";
	var_extensions["2"]  = ".g.2";

	MKL_INT nprow, npcol, myrow, mycol;
	blacs_gridinfo(&ictxt, &nprow, &npcol, &myrow, &mycol);


	if (myrank == 0)
		std::cout << "loading project: " << options.load_name << std::endl;

	int idx = 0;
	for (const auto &kv : options.variances) {
		const std::string &name = kv.first;
		double var = kv.second;
		if (name == "e") continue;

		if (myrank == 0)
			std::cout << "var " << name << " and idx " << var << std::endl;

		const std::string &ext = var_extensions.at(name);
		std::string load_file = options.load_name + ext;
		if (myrank == 0)
			std::cout << "reading file: " << load_file << std::endl;

		verify_and_set_mat_size(load_file, myrank, mkl_num_ind);

		MKL_INT nr = std::max(IONE_,
				numroc_(&mkl_num_ind, &options.mb, &myrow, &IZERO_, &options.nprow));
		MKL_INT nc = std::max(IONE_,
				numroc_(&mkl_num_ind, &options.nb, &mycol, &IZERO_, &options.npcol));
		MKL_INT mat_size = nr * nc;
		MKL_INT ldd = std::max(IONE_, nr);

		MKL_INT info;
		descinit_(descG,
				&mkl_num_ind, &mkl_num_ind,
				&options.mb, &options.nb,
				&IZERO_, &IZERO_,
				&ictxt, &ldd, &info);
		if (info != 0) {
			if (myrank == 0)
				std::cerr << "Error in descinit for " << name
					<< ", info=" << info << std::endl;
			blacs_gridexit_(&ictxt);
			exit(EXIT_FAILURE);
		}

		double *Gptr = static_cast<double*>(calloc(mat_size, sizeof(double)));
		if (!Gptr) {
			if (myrank == 0)
				std::cerr << "ERROR: allocation failed for GRM '" << name << "'" << std::endl;
			blacs_gridexit_(&ictxt);
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
		mem_tot += mat_size * sizeof(double);
		read_binmat(load_file.c_str(), Gptr, descG, mat_size);

		var_idx[name] = idx;
		idx_var[idx]  = name;
		Glocal[name]   = Gptr;
		h2[name]      = 0.0;
		++idx;
	}
}

int fillBuff( std::istream & infile, std::vector <char> & buff ) {
	infile.read( &buff[0], buff.size() );
	return infile.gcount();
}

int countLines( const std::vector <char> & buff, int bs ) {
	int newlines = 0;
	const char * p = &buff[0];
	for ( int i = 0; i < bs; i++ ) {
		if ( p[i] == '\n' ) {
			newlines++;
		}
	}
	return newlines;
}

void read_file_size(const char* filename, int &r, int &c){

	std::ifstream infile(filename, std::ios::in);
	if (infile.fail()){
		std::cout << "ERROR: can't open file: " << filename << std::endl;
		// TODO: throw proper error and close
	}
	std::string line;
	std::getline(infile, line);
	std::vector<std::string> words = str_split(" \r\t\n\x0b", line);
	words.erase(words.begin());
	c = words.size();
	r = 0;
	const int buf_size = 1024 * 1024;
	std::vector <char> buff(buf_size);
	while( int bs = fillBuff( infile, buff ) ) {
		r += countLines( buff, bs );
	}

	//std::cout << filename << " dimensions " << r << "," << c << std::endl;
	infile.close();
}


void verify_and_set_mat_size(const std::string& load_file,
		int                myrank_mpi,
		MKL_INT&           mkl_num_ind) {
	// read dimensions
	MKL_INT n_rows = 0, n_cols = 0;
	read_binmat_size(load_file.c_str(), n_rows, n_cols);

	if (myrank_mpi == 0)
		std::cout << load_file << " -> r,c: " << n_rows << "," << n_cols << '\n';

	bool bad_size   = (n_rows == 0 || n_cols == 0);
	bool not_square = (n_rows != n_cols);

	if (bad_size || not_square) {
		if (myrank_mpi == 0)
			std::cerr << "ERROR: " << load_file
				<< " has invalid size (" << n_rows << " x " << n_cols
				<< "). Must be non-zero and square.\n";

		std::cerr.flush();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	// First file sets the global size; later files must match
	if (mkl_num_ind == 0) {
		mkl_num_ind = n_rows;
		if (myrank_mpi == 0)
			std::cout << "Matrix size set to " << mkl_num_ind << std::endl;
	} else if (n_rows != mkl_num_ind) {
		if (myrank_mpi == 0)
			std::cerr << "ERROR: " << load_file << " dimensions (" << n_rows
				<< ") do not match previously established size ("
				<< mkl_num_ind << ").\n";

		std::cerr.flush();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
}


void read_binmat_size(const char* filename, MKL_INT& r, MKL_INT& c){
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	if (!in) { std::cerr << "ERROR: cannot open " << filename << std::endl; return; }
	in.read((char*) (&r),sizeof(MKL_INT));
	in.read((char*) (&c),sizeof(MKL_INT));
	in.close();
}


void read_binmat(const char *filename, double *mat, MKL_INT *descmat, MKL_INT nelem){
	int myrank_mpi, nprocs_mpi, ierr;
	const int ndims = 2;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
	MPI_File fh;
	MPI_Status status;
	MPI_Datatype MPI_darray;
	MPI_Aint lower_bound, darray_extent;
	// skip s for row,col info
	MPI_Offset skip = SKIP_BYTES_;
	//MKL_INT ictxt = descmat[1];
	MKL_INT ictxt = static_cast<MKL_INT>(descmat[1]);
	MKL_INT nelements = nelem;
	//blacs context
	MKL_INT nprow, npcol, myrow, mycol;
	blacs_gridinfo(&ictxt, &nprow, &npcol, &myrow, &mycol ); 

	int dims[ndims] = {(int) descmat[2], (int) descmat[3]};
	int distribs[ndims] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
	int dargs[ndims] = {(int) descmat[4],(int) descmat[5]};
	int proc_dims[ndims] = {static_cast<int>(nprow), static_cast<int>(npcol)};

	if (myrank_mpi == 0) {
		std::cout << "skip: " << skip << std::endl;
		std::cout << "dims: " << dims[0] << "," << dims[1] << std::endl;
		std::cout << "dargs: " << dargs[0] << "," << dargs[1] << std::endl;
	}

	MPI_Type_create_darray(
			nprocs_mpi, // size of process group (positive integer)
			myrank_mpi, // rank in process group (non-negative integer)
			ndims, //   number of array dimensions as well as process grid dimensions (positive integer)
			dims, // number of elements of type oldtype in each dimension of global array (array of positive integers)
			distribs, // distribution of array in each dimension (array of state)
			dargs, // distribution argument in each dimension (array of positive integers)
			proc_dims, // size of process grid in each dimension (array of positive integers)
			//MPI_ORDER_C, // array storage order flag (state)
			MPI_ORDER_FORTRAN, // array storage order flag (state)
			MPI_DOUBLE, // old datatype (handle)
			&MPI_darray // new datatype (handle)
			);
	MPI_Type_commit(&MPI_darray);
	MPI_Type_get_extent(MPI_darray, &lower_bound, &darray_extent);

	if (myrank_mpi == 0) {
		std::cout << "nelements " << nelements << std::endl;
		std::cout << "lower bound " << lower_bound << std::endl;
		std::cout << "darray extent " << darray_extent << std::endl;
	}

	if (myrank_mpi == 0) std::cout << "reading file: " << std::endl;
	ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	if(ierr != MPI_SUCCESS) {
		if (myrank_mpi == 0) std::cerr << "Error: Could not open file " << filename << std::endl;
		return;  
	} else {
		ierr = MPI_File_set_view(fh, skip, MPI_DOUBLE, MPI_darray, "native", MPI_INFO_NULL);
		if(ierr != MPI_SUCCESS) {
			if (myrank_mpi == 0) std::cerr << "Error: Could not set file view for " << filename << std::endl;
			MPI_File_close(&fh);  
			return;
		} else {
			if (myrank_mpi == 0) std::cout << "set view succeeded, skipping bytes " << skip << std::endl;
			ierr = MPI_File_read_all(fh, mat, nelements, MPI_DOUBLE, &status);
			if(ierr != MPI_SUCCESS) {
				if (myrank_mpi == 0) std::cerr << "Error: MPI_File_read_all failed for " << filename << std::endl;
				MPI_File_close(&fh);
				return;
			}

		}

	}
	if (myrank_mpi == 0) std::cout << std::endl;
	MPI_File_close(&fh);
}


// split strings based on delimiters
std::vector<std::string> str_split(const std::string &delimiter, const std::string &str)
{
	std::vector<std::string> words;

	int str_len = str.length();
	int del_len = delimiter.length();
	if (str_len==0 || del_len==0){
		return words; //no change
	}

	char *line = new char [str_len+1];
	strcpy(line, str.c_str());
	char *sep  = new char [del_len+1];
	strcpy(sep, delimiter.c_str());

	std::string tmp_str;
	char   *token = strtok(line,sep);
	while(token!=NULL){
		tmp_str = token;
		words.push_back(tmp_str);
		tmp_str.clear();
		token = strtok(NULL,sep);
	}
	delete[] line;
	delete[] sep;
	return words;
}

std::vector<int> str_split_int(const std::string &delimiter, const std::string &str)
{
	std::vector<int> words;

	int str_len = str.length();
	int del_len = delimiter.length();
	if (str_len==0 || del_len==0){
		return words;//no change
	}

	char *line = new char [str_len+1];
	strcpy(line, str.c_str());
	char *sep  = new char [del_len+1];
	strcpy(sep, delimiter.c_str());

	std::string tmp_str;
	char   *token = strtok(line,sep);
	while(token!=NULL){
		tmp_str = token;
		words.push_back(stoi(tmp_str));
		tmp_str.clear();
		token = strtok(NULL,sep);
	}
	delete[] line;
	delete[] sep;
	//delete token;
	return words;
}
