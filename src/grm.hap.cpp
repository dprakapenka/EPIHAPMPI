#include "constants.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>
#include <cassert>
#include "grm.hpp"
#include "mkl.h"
#include "mkl_scalapack.h"
#include "mkl_trans.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"
#include "mpi.h"

extern "C" void pdlaprnt_(MKL_INT*, MKL_INT*, double*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, char*, int const&, double *);

// check that blocks are not larger than matrix/proc
void check_block_size(int& b, int n, int np) {
    double per_proc = n/np;
    int limit = std::trunc(per_proc);
    if (b > limit) {
        b = limit;
        std::cout << "changed number of blocks to " << limit << std::endl;
    }

}

std::vector<int> str_split_int(const std::string &delimiter,
        const std::string &str)
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
    delete[] token;
    return words;
}

std::vector<std::string> str_split(const std::string &delimiter,
        const std::string &str)
{
    std::vector<std::string> words;

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
        words.push_back(tmp_str);
        tmp_str.clear();
        token = strtok(NULL,sep);
    }
    delete[] line;
    delete[] sep;
    delete[] token;
    return words;
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
        MPI_Finalize();
        exit(0);
    }
    std::string line;
    std::getline(infile, line);
    std::vector<std::string> words = str_split(" \r\t\n\x0b", line);
    words.erase(words.begin());
    c = words.size();
    if ((c % 2) != 0) {
        std::cout << "ERROR: haplotype file " << filename 
            << " has odd number of columns " << c << std::endl;
        throw;
    }

    r = 0;
    const int buf_size = 1024 * 1024;
    std::vector <char> buff(buf_size);
    while( int bs = fillBuff( infile, buff ) ) {
        r += countLines( buff, bs );
    }

    //std::cout << filename << " dimensions " << r << "," << c << std::endl;
    infile.close();
}

// read mat, slow 
void read_binmat_slow(const char *filename, double *mat, MKL_INT *descmat){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    MKL_INT rows;
    MKL_INT cols;
    in.read((char*) (&rows),sizeof(MKL_INT));
    in.read((char*) (&cols),sizeof(MKL_INT));
    double *rowbuff;
    rowbuff = (double *)malloc(rows*sizeof(double));
    in.read( (char *) rowbuff, rows*sizeof(double));

    for (MKL_INT i=1; i<=rows; i++){
        //std::cout << "i,j (" << i << "," << j << "): "<< rowbuff[i-1] << std::endl;
    }
    delete(rowbuff);
}

// write the size of matrix in first 16 bytes
void write_binmat_size_16(const char* filename, MKL_INT& r, MKL_INT& c){
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    out.write((char*) (&r),sizeof(MKL_INT));
    out.write((char*) (&c),sizeof(MKL_INT));
    //std::cout << filename << " dimensions " << r << "," << c << std::endl;
    out.close();
}

// skip
void write_binmat(const char *filename, double *mat, MKL_INT *descmat, int s, MKL_INT nelem){
    int myrank_mpi, nprocs_mpi, ierr;
    const int ndims = 2;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
    MPI_File fh;
    MPI_Status status;
    MPI_Datatype MPI_darray;
    MPI_Aint lower_bound, darray_extent;
    MPI_Offset skip = s;
    MKL_INT ictxt = descmat[1];
    MKL_INT nelements = nelem;
    MKL_INT nprow, npcol, myrow, mycol;
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    int dims[ndims] = {(int) descmat[2], (int) descmat[3]};
    int distribs[ndims] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
    int dargs[ndims] = {(int) descmat[4],(int) descmat[5]};
    int proc_dims[ndims] = {(int) nprow, (int) npcol};

    if (myrank_mpi == 0) {
        std::cout << "skip: " << skip << std::endl;
        std::cout << "dims: " << dims[0] << "," << dims[1] << std::endl;
        std::cout << "dargs: " << dargs[0] << "," << dargs[1] << std::endl;
        std::cout << "proc dims: " << proc_dims[0] << "," << proc_dims[1] << std::endl;
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

    //std::cout << "rank/nprocs " << myrank_mpi << "/" << nprocs_mpi <<  " nelements " << nelements  << " skip " << skip << std::endl;
    if (myrank_mpi == 0) {
        //    std::cout << "darray size: " << darray_size << std::endl;
        std::cout << "nelements " << nelements << std::endl;
        std::cout << "lower bound " << lower_bound << std::endl;
        std::cout << "darray extent " << darray_extent << std::endl;
    }

    if (myrank_mpi == 0) std::cout << "writing file: " << std::endl;
    ierr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if(ierr != MPI_SUCCESS) {
        if (myrank_mpi == 0) std::cout << "file open failed!" << std::endl;

    } else {
        ierr = MPI_File_set_view(fh, skip, MPI_DOUBLE, MPI_darray, "native", MPI_INFO_NULL);
        if(ierr != MPI_SUCCESS) {
            if (myrank_mpi == 0) std::cout << "file set view failed!" << std::endl;
        } else {
            if (myrank_mpi == 0) std::cout << "set view succeeded, skipping bytes " << skip << std::endl;
            ierr = MPI_File_write_all(fh, mat, nelements, MPI_DOUBLE, &status);
            if(ierr != MPI_SUCCESS) {
                if (myrank_mpi == 0) std::cout << "write_all failed!" << std::endl;
            } else {
                if (myrank_mpi == 0) std::cout << "write ierr " << ierr << std::endl;
            }
            //if (myrank_mpi == 0) std::cout << "status" << status << std::endl;

        }

    }
    if (myrank_mpi == 0) std::cout << std::endl;
    MPI_File_close(&fh);
}

int main(int argc, char **argv) {

    // Initialize MPI
    int myrank_mpi, nprocs_mpi;
    int num_threads;
    MKL_INT info;
    //unsigned long long int mem_try=0;
    //unsigned long long int  mat_size=0;
    MKL_INT mem_try=0;
    MKL_INT mat_size=0;
    double mem_tot=0.0;
    double mem_tot_gb=0.0;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
    MPI_File fh;
    MPI_Status status;
    // set error return
    MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_RETURN );

    // auto calculate as scquare as possible grid
    const int ndims=2;
    int dim_arr[ndims] = { 0 };
    info = MPI_Dims_create(nprocs_mpi, ndims, dim_arr);
    // processor rows
    MKL_INT nprow = dim_arr[0];
    // processor cols
    MKL_INT npcol = dim_arr[1];
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "suggested *square processor grid: " << nprow << "," << npcol << std::endl;

    //options
    int opt;
    std::string chr_filename;
    std::vector<std::string> chr_filenames;
    // row and col block size
    // TODO check size of smallest dim and set num_blocks to that size

    // helpers 
    int itmp;
    double dtmp;

    // matrices
    int *blocks_local;
    int *buffer_local;

    // timing
    double time_start;
    double time_start_read;
    double time_end_read;
    double time_start_stat;
    double time_end_stat;
    double time_start_dist;
    double time_end_dist;
    double time_start_mult;
    double time_end_mult;
    double time_start_write;
    double time_end_write;

    int skip = 16;

    std::string out_filename = "grm";
    MKL_INT mb = 2;
    MKL_INT nb = 2;


    // get num mkl threads
    num_threads = mkl_get_max_threads();


    if (myrank_mpi == MPI_ROOT_PROC_) printf("Usage: ./grm CHANGE\n");

    while((opt = getopt(argc, argv, ":r:c:n:o:f")) != -1)
    {
        switch(opt)
        {

            case 'n': // block size
                //std::cout << "block: " << optarg << std::endl;
                if (myrank_mpi == MPI_ROOT_PROC_) printf("block size: %i\n", atoi(optarg));
                mb = nb = atoi(optarg);
                break;
            case 'r': // processor rows
                if (myrank_mpi == MPI_ROOT_PROC_) printf("override num proc rows: %i\n", atoi(optarg));
                nprow = atoi(optarg);
                break;
            case 'c': // processor cols
                if (myrank_mpi == MPI_ROOT_PROC_) printf("override num proc cols: %i\n", atoi(optarg));
                npcol = atoi(optarg);
                break;
            case 'o': //input data (load)
                if (myrank_mpi == MPI_ROOT_PROC_) printf("proj_name: %s\n", optarg);
                out_filename = optarg;
                break;
            case 'f':
                if (myrank_mpi == MPI_ROOT_PROC_) printf("f flag set: %c\n", opt);
                break;
            case ':':
                if (myrank_mpi == MPI_ROOT_PROC_) printf("option needs a value\n");
                break;
            case '?':
                if (myrank_mpi == MPI_ROOT_PROC_) printf("unknown option: %c\n", optopt);
                break;
        }
    }

    for(; optind < argc; optind++) {
        chr_filenames.push_back(argv[optind]);
    }

    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "number of chromosomes: " << chr_filenames.size() << std::endl;
        std::cout << "mpi_procs: " << nprocs_mpi << " nprow: " << nprow << " npcol: " << npcol << std::endl;
        std::cout << "num_threads: " << num_threads << std::endl;
    }

    time_start = MPI_Wtime();
    //TODO: CHECK that columns is not less than number of procs
    // read number of inds
    int r = 0;
    int c = 0;
    //chr_filename = "hap_geno";
    read_file_size(chr_filenames[0].c_str(), r, c);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "found individuals: " << r << std::endl;
    }

    double *toprint;
    toprint = (double *)calloc(nb,sizeof(double)) ;

    // set up G matrix
    MKL_INT mkl_num_ind = r;

    // initialize BLACS
    MKL_INT ictxt, myid, myrow, mycol, nprocs;
    blacs_pinfo_(&myid, &nprocs) ; // BLACS rank and world size
    blacs_get_(&IZERO_, &IZERO_, &ictxt ); // -> Create context
    blacs_gridinit_(&ictxt, &CHAR_LAYOUT_, &nprow, &npcol ); // Context -> Initialize the grid
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    //gather row and column of myid
    MKL_INT myrowcol [2] = {myrow, mycol};
    MKL_INT rowcol[nprocs][2];
    MPI_Allgather(&myrowcol,
            2,
            MPI_LONG_LONG,
            &rowcol[0],
            2,
            MPI_LONG_LONG,
            MPI_COMM_WORLD);

    // size of the local matrices
    MKL_INT nr_ind    = std::max(IONE_, numroc_(&mkl_num_ind, &mb, &myrow, &IZERO_, &nprow)); // My proc -> row of local A
    MKL_INT nc_ind    = std::max(IONE_, numroc_(&mkl_num_ind, &nb, &mycol, &IZERO_, &npcol)); // My proc -> col of local A
    MKL_INT lldG = std::max(IONE_,nr_ind);

    // initialize local matrices
    double *G;
    mat_size = nr_ind * nc_ind;
    MKL_INT nelements = mat_size;
    mem_try = sizeof(double) * mat_size;

    if (myid == MPI_ROOT_PROC_) {
        std::cout << "nr_ind " << nr_ind << " * nc_ind " << nc_ind << " = "
            << mat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
        std::cout << "mkl_nr_ind " << nr_ind << " mkl_nc_ind " << nc_ind << " mkl_num_ind " << mkl_num_ind << std::endl;
        std::cout << "trying to allocate G local: " << mem_try/1000000.0 << "MB" << std::endl;
    }
    //G = (double *)mkl_malloc(mat_size,  64) ;
    G = (double *)calloc(mat_size, sizeof(double)) ;
    if (G==NULL){
        std::cout << "ERROR: memory allocation for G on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
        exit(0);
    } else {
        if (myid == MPI_ROOT_PROC_) std::cout << "...success" << std::endl;
        mem_tot += mem_try;
    }
    // matrix descriptor
    MKL_INT descG[DESC_LEN_];
    descinit_( descG,  &mkl_num_ind, &mkl_num_ind, &mb, &nb, &IZERO_, &IZERO_, &ictxt, &lldG, &info);
    if(info != 0) {
        printf("Error in descinit, info = %d\n", int(info));
    }

    //iterate over chrs
    for (int chr = 0; chr < chr_filenames.size(); chr++) { 

        chr_filename = chr_filenames[chr];
        if (myrank_mpi == MPI_ROOT_PROC_) {
            std::cout << "loading chr: " << chr_filename << std::endl;
        }

        // read number of inds and cols
        //chr_filename = "hap_geno";
        read_file_size(chr_filename.c_str(), r, c);
        if (myid == MPI_ROOT_PROC_) {
            std::cout << "r: " << r << " c: " << c << std::endl;
        }
        if (mkl_num_ind != (MKL_INT) r) {
            std::cout << "ERROR: mismatch in number of inds in " << chr_filename << std::endl;
        }

        int using_mpi_tasks = nprocs_mpi;
        if (c % 2 != 0) {
            std::cout << "ERROR: number of columns in hap file must be divisible by 2" << std::endl;
        }
        else if ((c/2) < nprocs_mpi) {
            using_mpi_tasks = c/2;
            if (myid == MPI_ROOT_PROC_) {
                std::cout << "warning: number of tasks is more than number of hap blocks in chr," << 
                    "using " << using_mpi_tasks << " tasks for reading and block statistics" << std::endl;
            }
        }

        int columns [nprocs_mpi];
        memset(columns, 0, nprocs_mpi*sizeof(int));
        int displ [nprocs_mpi];
        int remainder = 0;
        remainder = (c / 2) % using_mpi_tasks;
        int max_cols = 0;
        int block_size = 0;

        for ( int i = using_mpi_tasks-1; i >=0; --i) {
            columns[i] = (c /2) / using_mpi_tasks;
            if (remainder > 0) {
                columns[i]++;
                remainder--;
            }
            columns[i]*=2;
        }


        int sum = 0;
        for ( int i = 0; i < nprocs_mpi; i++) {
            displ[i] = sum;
            sum += columns[i];
            if (columns[i] > max_cols)
            {
                max_cols = columns[i];
            }
        }


        /*
           for (int mrank=0; mrank < nprocs_mpi; mrank++){
           std::cout << "rank: " << mrank << " columns: " << columns[mrank] << " displacements: " << displ[mrank] << " max cols: " << max_cols << std::endl;
           }
           */

        //allocate local blocks
        block_size = r * columns[myid];
        //std::vector<int> my_block ;
        mem_try = sizeof(int) * block_size;
        for (int mrank=0; mrank < nprocs_mpi; mrank++){
            if (myid == mrank) {
                std::cout << "proc " << mrank << " blocks_local: " << block_size << " * " << sizeof(int) << " " << mem_try << " bytes" << std::endl;
            }
        }
        if (myid == MPI_ROOT_PROC_) std::cout << "trying to allocate blocks_local: " << mem_try/1000000 << "MB" << std::endl;
        //blocks_local = (int *)mkl_malloc(block_size,  64) ;
        blocks_local = (int *)calloc(block_size, sizeof(int)) ;
        if (blocks_local==NULL){
            std::cout << "Error of memory allocation blocks_local on proc " << myid << std::endl;
            MPI_Finalize();
            exit(0);
        } else {
            //std::cout << "blocks_local allocated on " << myid << std::endl;
            mem_tot += mem_try;
        }


        //allocate local buffers
        //block_size = max_cols;
        block_size = columns[myid];
        //std::vector<int> my_block ;
        mem_try = sizeof(int) * block_size;
        //buffer_local = (int *)mkl_malloc(block_size,  64) ;
        buffer_local = (int *)calloc(block_size, sizeof(int)) ;
        if (buffer_local==NULL){
            std::cout << "Error of memory allocation buffer_local on proc " << myid << std::endl;
            MPI_Finalize();
            exit(0);
        }
        /*
           else {
           std::cout << "buffer_local allocated on " << myid << std::endl;
           mem_tot += mem_try;
           }
           */


        std::string ind;
        std::ifstream infile(chr_filename, std::ios::in);
        std::string line;
        std::getline(infile, line);
        std::vector<std::string> header = str_split(" \r\t\n\x0b", line);
        std::vector<int> words;
        //TODO: read x lines at a time then send
        time_start_read = MPI_Wtime();


        for (int row = 0; row < r; row++) {
            if (myid == MPI_ROOT_PROC_) {
                std::getline(infile, line);
                words = str_split_int(" \r\t\n\x0b", line);
                words.erase(words.begin());
                //memcpy(&blocks_local[0*c],&words[0],sizeof(int)*c);
                //for (int i=0; i < c; i++) std::cout << blocks_local[i] << ",";
            }
            /*
            //for (auto i: words) std::cout << i << ",";
            //std::cout << std::endl;
            if (myid == MPI_ROOT_PROC_) {
            for (int i=0; i < c; i++) std::cout << blocks_local[i] << ",";
            std::cout << std::endl;
            }
            */


            //MPI_Scatterv(blocks_local, columns, displ, MPI_INT, buffer_local, max_cols, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
            MPI_Scatterv(&words[0], columns, displ, MPI_INT, buffer_local, max_cols, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
            memcpy(&blocks_local[row*columns[myid]],buffer_local,sizeof(int)*columns[myid]);
            /*
               for ( int j = 0; j < nprocs_mpi; j++) {
               if (myid == j) {
               for (int i=0; i < columns[j]; i++) std::cout << j << "_buffer_local__ " << buffer_local[i] << ",";
               std::cout << std::endl;
               }
               }
               */


        }
        infile.close();
        free(buffer_local);

        time_end_read = MPI_Wtime(); // end of read data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished reading chr, time %fs.\n", myid, time_end_read-time_start_read);

        //hap_block_het_geno  hb_het_geno;


        time_start_stat = MPI_Wtime();
        int il0;
        int il1;
        int n_col = columns[myid]/2;
        int max_n_col = max_cols/2;
        int tmp_row = 2;
        hb_allele.resize(n_col);
        hb_het_geno.resize(n_col);
        //std::vector<int> hap_max_key(n_col);
        hap_max_key.resize(n_col);

        for (int col=0; col < n_col; col++) {
            for (int row = 0; row < r; row++){
                //index of allele 0
                il0 = (row *columns[myid])+(col*2);
                //index of allele 1
                il1 = (row *columns[myid])+(col*2)+1;
                //std::cout << myid << "...." <<  blocks_local[il0] << ", " <<  blocks_local[il1] << " | ";


                int geno_tmp_0 = blocks_local[il0];
                int geno_tmp_1 = blocks_local[il1];
                allele_list::const_iterator got;
                if (geno_tmp_0 != 0) {
                    got = hb_allele[col].find(geno_tmp_0);
                    if (got == hb_allele[col].end()) {
                        hb_allele[col][geno_tmp_0].count = 1.0;
                    }
                    else {
                        hb_allele[col][geno_tmp_0].count += 1.0;
                    }
                }
                if (geno_tmp_1 != 0) {
                    got = hb_allele[col].find(geno_tmp_1);
                    if (got == hb_allele[col].end()) {
                        hb_allele[col][geno_tmp_1].count = 1.0;
                    }
                    else {
                        hb_allele[col][geno_tmp_1].count += 1.0;
                    }
                }
                std::array<int, 2> het_geno;
                het_geno_list::const_iterator got_het_geno;
                if (geno_tmp_0 != 0 && geno_tmp_1 != 0 && geno_tmp_0 != geno_tmp_1) {
                    if (geno_tmp_0 < geno_tmp_1) {
                        het_geno[0] = geno_tmp_0;
                        het_geno[1] = geno_tmp_1;
                    }
                    else if (geno_tmp_0 > geno_tmp_1) {
                        het_geno[0] = geno_tmp_1;
                        het_geno[1] = geno_tmp_0;
                    }
                    got_het_geno = hb_het_geno[col].find(het_geno);
                    if (got_het_geno == hb_het_geno[col].end()) {
                        hb_het_geno[col][het_geno].count = 1.0;
                    }
                    else {
                        hb_het_geno[col][het_geno].count += 1.0;
                    }
                }

            }
        }


        //hap_max_key.resize(n_col);

        for (int i = 0; i < n_col; i++)
        {
            hb_allele[i].erase(0);
            double sum = 0.0;
            int max_tmp = hb_allele[i].begin()->first;
            for (auto it = hb_allele[i].begin(); it != hb_allele[i].end(); ++it)
            {
                if ((hb_allele[i][max_tmp]).count < (it->second).count) {
                    max_tmp = it->first;
                }
                sum += (it->second).count;
            }
            //std::cout << "haplotype block " << i << " max_key : "<< max_tmp << std::endl;
            hap_max_key[i] = max_tmp;
            int count = 0;
            for (auto it = hb_allele[i].begin(); it != hb_allele[i].end(); ++it)
            {
                (it->second).count /= sum;
                if (it->first != hap_max_key[i])
                {
                    (it->second).idx = count;
                    count++;
                }
                else
                {
                    (it->second).idx = -1;
                }
                /*std::cout << it->first << ":"
                  << it->second.idx << "-"
                  << it->second.count << std::endl; */
            }

            count = 0;
            for (auto it = hb_het_geno[i].begin(); it != hb_het_geno[i].end(); ++it)
            {
                //(it->second).count /= sum;
                (it->second).idx = count;
                count++;
            }
            //hap_block_pos(0,i) = hap_count;
            //hap_count += (hb_allele[i].size()-1);
            //hap_block_pos(1,i) = hap_count-1;
            //std::cout << "hap_block_pos " << std::endl << hap_block_pos.col(i) << endl;
        }


        //MKL_INT col_hap_count_tot = 0;
        MKL_INT col_hap_count[max_n_col][nprocs];
        MKL_INT col_hap_disp[max_n_col][nprocs];
        MKL_INT col_hap_count_total[max_n_col];
        MKL_INT task_hap_count[max_n_col];

        //dont need
        MKL_INT task_hap_displ[max_n_col];
        MKL_INT all_task_hap_counts[nprocs];

        MKL_INT hap_count = 0;
        MKL_INT w_hap_idx = 0;
        double *W_hap;
        double *W_tmp;
        MKL_INT descW_tmp[nprocs][DESC_LEN_];
        MKL_INT descW_hap[DESC_LEN_];
        MKL_INT descW_me[DESC_LEN_];
        MKL_INT descempty[DESC_LEN_];
        memset(descempty, 0, DESC_LEN_*sizeof(MKL_INT));
        // gather all hap_counts for each block in task
        for (int col=0; col < max_n_col; col++)
        {
            if (col >= n_col){
                task_hap_count[col] = 0;
            }
            else
            {
                //task_hap_displ[col] = task_hap_count_tot;
                task_hap_count[col] = hb_allele[col].size() - 1;
            }
            MPI_Allgather(&task_hap_count[col],
                    1,
                    MPI_LONG_LONG,
                    &col_hap_count[col],
                    1,
                    MPI_LONG_LONG,
                    MPI_COMM_WORLD);

        }
        for (int col=0; col < max_n_col; col++) {
            col_hap_count_total[col] = 0;
            for (int me=0; me < nprocs; me++) {
                col_hap_disp[col][me] = col_hap_count_total[col];
                col_hap_count_total[col] += col_hap_count[col][me];
            }
        }

        time_end_stat = MPI_Wtime(); // end of stat data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished calculating statistics chr, time %fs.\n", myid, time_end_stat-time_start_stat);

        time_start_dist = MPI_Wtime(); // start of dist data


        for (int col=0; col < max_n_col; col++) {

            //MKL_INT tmpcol;
            // if there are more columns and if the haplotype block has more than one "genotype"
            // CHANGED: added hb_allele size check - TEST_THIS
            if (col < n_col  && hb_allele[col].size() > 1) {
                hap_count = hb_allele[col].size() - 1;
                if (hap_count == 0) {
                    std::cout << "hap_count is ZERO" 
                        << "i_col: " << col
                        << " rank: " << myid
                        << " hap_count " << hap_count
                        << " mat_size: " << mat_size << std::endl;


                }
                mat_size = mkl_num_ind * hap_count;
                W_tmp = (double *)calloc(mat_size, sizeof(double));
                if (W_tmp==NULL){
                    std::cout << "ERROR: memory allocation for W_hap on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
                    MPI_Finalize();
                    exit(0);
                } 

                for (MKL_INT row = 0; row < mkl_num_ind; row++){
                    //index of allele 0
                    il0 = (row *columns[myid])+(col*2);
                    //index of allele 1
                    il1 = (row *columns[myid])+(col*2)+1;
                    int geno_tmp_0 = blocks_local[il0];
                    int geno_tmp_1 = blocks_local[il1];
                    //MKL_INT whapidx = row*hap_count + tmpcol;
                    //MKL_INT whapidx = row*mkl_num_ind + tmpcol;
                    /*
                       if (myid == MPI_ROOT_PROC_) {
                       std::cout << "id: " << myid << " col: " << col << " allele[" << il0 << "," << il1 << "]: " << geno_tmp_0 << "," << geno_tmp_1 << std::endl;
                       }
                       */
                    if (geno_tmp_0 == 0 || geno_tmp_1 == 0) {
                        for (int k = 0; k < hb_allele[col].size() - 1; k++) {
                            int tmpcol = k;
                            W_tmp[tmpcol*mkl_num_ind + row] = 0.0;
                        }
                    }
                    else {
                        for (auto it = hb_allele[col].begin(); it != hb_allele[col].end(); ++it) {
                            if (it->second.idx >= 0) {
                                int tmpcol = (it->second.idx);
                                W_tmp[tmpcol*mkl_num_ind + row] = (it->second.count) * 2.0;
                            }
                        }
                        if (geno_tmp_0 == hap_max_key[col] && geno_tmp_1 == hap_max_key[col]) {
                        }
                        else if (geno_tmp_0 != hap_max_key[col] && geno_tmp_1 == hap_max_key[col]) {
                            int tmpcol = hb_allele[col][geno_tmp_0].idx;
                            W_tmp[tmpcol*mkl_num_ind + row] -= 1.0;
                        }
                        else if (geno_tmp_0 == hap_max_key[col] && geno_tmp_1 != hap_max_key[col]) {
                            int tmpcol = hb_allele[col][geno_tmp_1].idx;
                            W_tmp[tmpcol*mkl_num_ind + row] -= 1.0;
                        }
                        else if (geno_tmp_0 != hap_max_key[col] && geno_tmp_1 != hap_max_key[col]) {
                            int tmpcol = hb_allele[col][geno_tmp_0].idx;
                            W_tmp[tmpcol*mkl_num_ind + row] -= 1.0;
                            tmpcol = hb_allele[col][geno_tmp_1].idx;
                            W_tmp[tmpcol*mkl_num_ind + row] -= 1.0;
                        }

                        /*
                           if (row*mkl_num_ind + tmpcol >= mat_size){
                           std::cout << "id: " << myid << " row: " << row << " col: " << col << " tmpcol: " << tmpcol <<
                           " idx: " << row*hap_count+tmpcol << " matsize: " << mat_size << " idx2: " << row*mkl_num_ind + tmpcol << std::endl;
                           std::cout << "iam: " << myid << " FUCKED" << std::endl;
                           }
                           if (myid == MPI_ROOT_PROC_) {
                           std::cout << "id: " << myid << " col: " << col << " w_tmp[" << row << "," << tmpcol << "]: " << W_tmp[row*mkl_num_ind + tmpcol] << std::endl;
                           }
                           */

                    }

                } // end of create W_tmp for block

                // print wtmp
                /*
               if (myid == MPI_ROOT_PROC_)
                {
                   std::cout << "id: " << myid << " hapcount: " << hap_count << " col:" << col << std::endl;
                    for (MKL_INT row = 0; row < mkl_num_ind; row++){
                        for (int h=0; h<hap_count; h++)
                        {
                           std::cout << W_tmp[row+mkl_num_ind*h] << ", " ;
                        }
                        std::cout << std::endl;
                    }
                }
                */



                // describe matrix of the block on each this task
                descinit_( descW_me,
                        &mkl_num_ind, &hap_count,
                        &mkl_num_ind, &hap_count,
                        &myrow, &mycol,
                        &ictxt,
                        &mkl_num_ind,
                        &info);
                if(info != 0) {
                    //error on 2snp, hap_count = 0
                    std::cout << "ERROR has hap_count: " << hap_count << std::endl;
                    std::cout << "ERROR in descW_me descriptor: " << info << std::endl;
                }

            } // end if has more cols
            else {
                //set descriptor to 0
                //hap_count = 1;
                //std::cout << "iter_col: " << col << " rank: " << myid << " hap_count " << hap_count << std::endl;
                W_tmp = (double *)calloc(1, sizeof(double));
                memcpy(descW_me, descempty, sizeof(descW_me));
            }

            // gather all descriptors
            MPI_Allgather(&descW_me,
                    DESC_LEN_,
                    MPI_LONG_LONG,
                    &descW_tmp,
                    DESC_LEN_,
                    MPI_LONG_LONG,
                    MPI_COMM_WORLD);

            // size of the local matrices
            MKL_INT nc_w_hap    = std::max(IONE_, numroc_(&col_hap_count_total[col], &nb, &mycol, &IZERO_, &npcol)); // My proc -> col of local W
            MKL_INT lldW_hap = std::max(IONE_,nr_ind);

            // initialize local matrices
            mat_size = nr_ind * nc_w_hap;
            mem_try = sizeof(double) * mat_size;

            if (myid == MPI_ROOT_PROC_) {
                //std::cout << "w_cols " << col_hap_count_total[col] << std::endl;
                //std::cout << "nr_ind " << nr_ind << " * nc_w_hap " << nc_w_hap << " = "
                //    << mat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
                //std::cout << "trying to allocate W local: " << mem_try/1000000.0 << "MB" << std::endl;
            }

            W_hap = (double *)calloc(mat_size, sizeof(double)) ;
            if (W_hap==NULL){
                std::cout << "ERROR: memory allocation for W_hap on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
                MPI_Finalize();
                exit(0);
            } else {
                mem_tot += mem_try;
            }
            // matrix descriptor
            descinit_( descW_hap,  &mkl_num_ind, &col_hap_count_total[col], &mb, &nb, &IZERO_, &IZERO_, &ictxt, &lldW_hap, &info);
            if(info != 0) {
                printf("Error in descinit W_hap, info = %d\n", int(info));
            }


            /*
               for (int col=0; col < max_n_col; col++) {
               for (int me=0; me < nprocs; me++) {
               std::cout << " id: " << myid << ":: col: " << col << " on proc: " <<  me <<
               " col_hap_count: " <<  col_hap_count[col][me]  <<
               " col_hap_disp: " <<  col_hap_disp[col][me]  <<
               " col_hap_count: " <<  col_hap_count_total[col]  << std::endl;
               }
               }

               for (int me=0; me < nprocs; me++) {
               std::cout << "task " << myid << " thinks task " << me << " has desc ";
               for (int ass=0; ass < DESC_LEN_; ass++) {
               std::cout << descW_tmp[me][ass] << ", " ;
               }
               std::cout << std::endl;
               }
               */

            // distribute to W_hap for this block
            for (int me=0; me < nprocs; me++) {
                // tasks with no column do not participate
                if (descW_tmp[me][3] == 0) {
                    //std::cout << "EMPTY DESCRIPTOR ON : " << me << std::endl;
                }
                else {
                    w_hap_idx = col_hap_disp[col][me]+1;
                    // std::cout << "task " << myid << " distributing mat from " << me << 
                    //     " has big index: " << w_hap_idx << " of size: " << col_hap_count[col][me] << std::endl;
                    pdgeadd_(&CHAR_NOTRANS_,
                            &mkl_num_ind, &col_hap_count[col][me],
                            &DONE_,
                            W_tmp,   &IONE_, &IONE_, descW_tmp[me],
                            &DZERO_,
                            W_hap, &IONE_, &w_hap_idx, descW_hap);
                }
            }

            /*
               if (myid==MPI_ROOT_PROC_) {
               std::cout << myid << " W_tmp: " <<  W_tmp[0] << ", " << W_tmp[1] << " | " <<
               "W_hap: " <<  W_hap[0] << ", " << W_hap[1] << std::endl;

               }
               */

            free(W_tmp);

            /*
               if (col == 0) {
               std::string tmp_filename = "wmat";
               tmp_filename = tmp_filename + std::to_string(col);
               write_binmat_size_16(tmp_filename.c_str(), mkl_num_ind, col_hap_count_total[col]);
               write_binmat(tmp_filename.c_str(), W_hap, descW_hap, skip, mat_size);
               }
               */
            //pdlaprnt_(&mkl_num_ind, &col_hap_count_total[col], W_hap, &IONE_, &IONE_, descW_hap, &IZERO_, &IZERO_, "W", 6, toprint);

            //time_start_mult = MPI_Wtime(); // start of mult
            //multiply while adding. G was initiated to 0 on first chr
            pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
                    &mkl_num_ind, &mkl_num_ind, &col_hap_count_total[col],
                    &DONE_,
                    W_hap,   &IONE_, &IONE_, descW_hap,
                    W_hap,   &IONE_, &IONE_, descW_hap,
                    &DONE_,
                    G, &IONE_, &IONE_, descG);

            //time_end_mult = MPI_Wtime(); // end of mult
            //printf("[proc_%lli] finished multiplying and adding WW^t to G, time %fs.\n", myid, time_end_mult-time_start_mult);
            free(W_hap);

            /*
               if (myid==MPI_ROOT_PROC_) {
               for (int iii=0;iii<5;iii++){
               std::cout << myid << " G: " <<  G[iii] << ", " ;
               }
               std::cout << std::endl;
               for (int iii=0;iii<5;iii++){
               std::cout << myid << " G: " <<  G[mkl_num_ind*iii] << ", " ;
               }
               std::cout << std::endl;
               }
               */

        } //end hap column iterations
        time_end_dist = MPI_Wtime(); // end of dist data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished multiplication and distribution, time %fs.\n", myid, time_end_dist-time_start_dist);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished chr %s, time %fs.\n", myid, chr_filename.c_str(), time_end_dist-time_start_read);
    } // end of chr iteration

    free(blocks_local);

    /*
       if (myid==MPI_ROOT_PROC_) {
       for (int iii=0;iii<5;iii++){
       std::cout << myid << " G: " <<  G[iii] << ", " ;
       }
       std::cout << std::endl;
       }
       */
    double trace = pdlatra_(&mkl_num_ind, G, &IONE_, &IONE_, descG);
    trace = trace / mkl_num_ind;
    if (myid == 0) std::cout << "trace: " << trace << std::endl;
    for (int element = 0; element < nelements; element++)
    {
        G[element] = G[element] / trace ;

    }
    if (myid == 0) std::cout << "divide by trace completed " << std::endl;
    //pdlaprnt_(&mkl_num_ind, &mkl_num_ind, G, &IONE_, &IONE_, descG, &IZERO_, &IZERO_, "G", 6, toprint);

    time_start_write = MPI_Wtime(); // start of write
    out_filename = out_filename + ".g.AH";
    if (myid == 0) {
        std::cout << "writing binary grm to " << out_filename << std::endl;
    }
    write_binmat_size_16(out_filename.c_str(), mkl_num_ind, mkl_num_ind);
    write_binmat(out_filename.c_str(), G, descG, skip, nelements);

    time_end_write = MPI_Wtime(); // end of write
    if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing to binary file, time %fs.\n", myid, time_end_write - time_start_write);
    if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
    printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);

    free(G);

    //time_end = MPI_Wtime();
    //printf("[proc_%lli] iteration done, time %fs.\n", myid, time_end-time_start);

    mem_tot_gb = mem_tot/1000000000;
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank_mpi == 0) {
        printf("mem per task: %f\n", mem_tot_gb);
        printf("time: %f", time_end_write-time_start);
    }


    MPI_Finalize();
    return EXIT_SUCCESS;
    return 0;

}
