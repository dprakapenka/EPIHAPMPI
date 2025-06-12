#include "reml.h"
#include "gblup.h"
#include "options.h"
#include "functions.h"
#include "read.h"
#include "constants.h"
#include "print.h"

#include <iostream>
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <sys/time.h>
#include <algorithm>
#include <cassert>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>  // for memset, memcpy, etc.
#include <numeric>  // for std::accumulate
#include "mkl.h"
#include "mkl_scalapack.h"
#include "mkl_trans.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"
#include "mpi.h"


int main(int argc, char **argv) {

    // Initialize MPI
    const MKL_INT DESC_LEN_ = 9;
    int myrank_mpi, nprocs_mpi;
    MKL_INT info;

    MKL_INT mem_try=0;
    MKL_INT mat_size=0;
    MKL_INT mem_tot=0.0;

    int MPI_ROOT_PROC_ = 0;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
    //MPI_File fh;
    MPI_Status status;
    // set error return
    MPI_Comm_set_errhandler( MPI_COMM_WORLD, MPI_ERRORS_RETURN );

    MPI_Datatype MPI_MKL_INT;
    MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(MKL_INT), &MPI_MKL_INT);

    // timing 
    double time_start;
    time_start = MPI_Wtime();

    //options
    Options options;

    // Compute close to square BLACS 2D grid
    std::tie(options.nprow, options.npcol) = compute_blacs_grid(nprocs_mpi);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "Suggested *square processor grid: " << options.nprow << "," << options.npcol << std::endl;
    }

    // TODO check size of smallest dim and set num_blocks to that size

    double dtmp;

    MKL_INT lwork;

    std::map<std::string,double>::iterator var_iter;
    std::map<std::string,double*> aiv;

    // get options
    parse_options(argc, argv, myrank_mpi, options);


    // print options
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "parsed user options" << std::endl;

    }

    // copy variances to tmp_var
    std::map<std::string,double> tmp_var(options.variances);
    std::map<std::string,double> delta_var(options.variances);
    std::map<std::string,double> delta_ai(options.variances);

    // initialize BLACS grid
    MKL_INT iam, nprocs;
    MKL_INT ictxt, myrow, mycol;
    blacs_pinfo(&iam, &nprocs) ; // BLACS rank and world size
    blacs_get(&IZERO_, &IZERO_, &ictxt ); // -> Create context
    blacs_gridinit(&ictxt, &CHAR_LAYOUT_, &options.nprow, &options.npcol ); // Context -> Initialize the grid
    blacs_gridinfo(&ictxt, &options.nprow, &options.npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    // DEBUG variable to print distributed matrix
    double *toprint;
    toprint = (double *)calloc(options.nb,sizeof(double)) ;
    update_mem(toprint, "toprint", mem_try, mem_tot, ictxt);

    // TODO: check options.mb and options.nb, should be no bigger than inds
    if (myrank_mpi == MPI_ROOT_PROC_){
        std::cout << "blacs grid:" << std::endl
            << "iam: " << iam << std::endl
            << " nprocs: " << nprocs << std::endl
            << " mb: " << options.mb << std::endl
            << " nb: " << options.nb << std::endl
            << " nprow: " << options.nprow << std::endl
            << " npcol: " << options.npcol << std::endl;

    }

    // print options
    if (myrank_mpi == MPI_ROOT_PROC_) {
        print_options(options);

    }


    MKL_INT nr_ind;
    MKL_INT nc_ind;
    MKL_INT lddG;
    MKL_INT descG[DESC_LEN_];

    std::map<std::string,int>    var_idx;
    std::map<int,std::string>    idx_var;
    std::map<std::string,double*>        Glocal;
    std::map<std::string,double> h2;
    std::map<std::string,double> h2_new;


    // read G matrices
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "loading project: " << options.load_name << std::endl;

    // TODO: check mkl_num_ind from id file in each grm
    std::vector<std::string> individual_ids;
    MKL_INT mkl_num_ind = 0; // Will be determined by reading IDs


    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::string id_filename = options.load_name + ".id.txt";
        std::cout << "Reading individual IDs from: " << id_filename << std::endl;
        int id_read_status = read_individual_ids(ictxt, id_filename, individual_ids);
        if (id_read_status != 0) {
            std::cerr << "ERROR: Failed to read individual IDs." << std::endl;
            // Signal other processes to abort before calling MPI_Abort
            mkl_num_ind = -1; // Use a sentinel value
        } else {
            mkl_num_ind = static_cast<MKL_INT>(individual_ids.size());
        }
    }
    // Broadcast the number of individuals (or error signal)
    MPI_Bcast(&mkl_num_ind, 1, MPI_MKL_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
    if (mkl_num_ind <= 0) {
        if (myrank_mpi != MPI_ROOT_PROC_) {
            std::cerr << "Rank " << myrank_mpi << " aborting due to error on root during ID read." << std::endl;
        }
        MPI_Finalize();
        return 1; // Exit after error
    }
    // for now leave on root

    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "Total individuals from ID file: " << mkl_num_ind << std::endl;

    // Load all G matrices into the map
    load_grms(
            options,
            ictxt,
            mkl_num_ind,
            var_idx,
            idx_var,
            Glocal,
            descG,
            h2,
            mem_tot
            );

    // size of the local matrices G above. TODO: pass them instead?
    nr_ind    = std::max(IONE_, numroc_(&mkl_num_ind, &options.mb, &myrow, &IZERO_, &options.nprow)); // My proc -> row of local A
    nc_ind    = std::max(IONE_, numroc_(&mkl_num_ind, &options.nb, &mycol, &IZERO_, &options.npcol)); // My proc -> col of local A
    lddG = std::max(IONE_,nr_ind);


    //pdlaprnt_(&mkl_num_ind, &mkl_num_ind, Gh, &IONE_, &IONE_, descG, &IZERO_, &IZERO_, "Gh", 6, toprint);


    //load Y matrix ====================================

    //temporary until covars are implemented
    MKL_INT mkl_Ycol = 1;

    // --- Read Phenotype Data on Root ---
    int num_rec = 0;
    //TODO: remove
    std::vector<int>      rec2ind;
    std::vector<double> Y_values;           // Store actual Y values

    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "Reading phenotype data from: " << options.pheno_file << std::endl;
        // Pass the id_to_index map created earlier

        int pheno_read_status = read_phenotype_data(
                options,
                ictxt,
                individual_ids,
                rec2ind,
                num_rec,
                Y_values
                );
        std::cout << "nummmmmreeeeeec " << num_rec <<"." << std::endl;

        if (pheno_read_status != 0) {
            std::cerr << "ERROR: Failed to process phenotype data. Aborting." << std::endl;
            num_rec = -1; // Signal error
        } else {
            std::cout << "Found " << num_rec << " valid phenotype records matching IDs." << std::endl;
        }
    }

    // Broadcast num_rec 
    MPI_Bcast(&num_rec, 1, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
    if (num_rec <= 0) {
        if (myrank_mpi != MPI_ROOT_PROC_) {
            std::cerr << "Rank " << myrank_mpi << " aborting due to error on root during phenotype read." << std::endl;
        }
        // Cleanup GRMs
        for(auto const& [key, val] : Glocal) { free(val); }
        MPI_Finalize();
        return 1;
    }
    // ---

	// Broadcast rec2ind
	if (myrank_mpi != MPI_ROOT_PROC_)
		rec2ind.resize(num_rec);
	MPI_Bcast(rec2ind.data(), num_rec, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);

    MKL_INT mkl_num_rec = 0;
	mkl_num_rec = (MKL_INT) num_rec;


    MKL_INT nr_num_rec  = std::max(IONE_, numroc_(&mkl_num_rec, &options.mb, &myrow, &IZERO_, &options.nprow)); // My proc -> row of local A
    MKL_INT nc_num_rec  = std::max(IONE_, numroc_(&mkl_num_rec, &options.nb, &mycol, &IZERO_, &options.npcol)); // My proc -> col of local A
    MKL_INT nr_one      = std::max(IONE_, numroc_(&IONE_, &options.mb, &myrow, &IZERO_, &options.nprow)); // My proc -> row of local A
    MKL_INT nc_one      = std::max(IONE_, numroc_(&IONE_, &IONE_, &mycol, &IZERO_, &options.npcol)); // My proc -> col of local A
    MKL_INT lddY = std::max(IONE_,nr_num_rec);
    MKL_INT lddYt = std::max(IONE_,nr_one);


    // initialize local matrices
    double *Ylocal;
    // calloc initializes to 0, malloc doesn't initialize
    // in bytes
    mat_size = nr_num_rec*nc_one;
    mem_try = sizeof(double) * mat_size;
    Ylocal = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(Ylocal, "Y", mem_try, mem_tot, ictxt);

    // matrix descriptor Y
    MKL_INT descY[DESC_LEN_];
    descinit_( descY,  &mkl_num_rec, &mkl_Ycol, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddY, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }

    // matrix descriptor Yt
    MKL_INT descYt[DESC_LEN_];
    //descinit_( descYt,  &mkl_Ycol, &mkl_num_rec, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddYt, &info);
    descinit_( descYt,  &mkl_Ycol, &mkl_num_rec, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &IONE_, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }


    // --- Distribute Y ---
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "Distributing Y matrix..." << std::endl;
    for(MKL_INT i=0; i < mkl_num_rec; ++i) {
        MKL_INT global_row = i + 1; // 1-based index for ScaLAPACK
        MKL_INT global_col = 1;
        double val_to_set = (myrank_mpi == MPI_ROOT_PROC_) ? Y_values[i] : 0.0;
        pdelset_(Ylocal, &global_row, &global_col, descY, &val_to_set);
    }
    blacs_barrier_(&ictxt, &CHAR_ACHAR_); // Ensure all sets are processed
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "Finished distributing Y matrix." << std::endl;
        Y_values.clear(); Y_values.shrink_to_fit(); // Clear root's copy
    }

    //pdlaprnt_(&mkl_num_rec, &mkl_Ycol, Ylocal, &IONE_, &IONE_, descY, &IZERO_, &IZERO_, "Y", 6, toprint);

    //make X matrix
    //temporary until covars are implemented
    MKL_INT mkl_num_c = 1;
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "num_rec,num_c: " << mkl_num_rec << "," << mkl_num_c << std::endl;

    // size of the local matrices
    MKL_INT nr_num_c    = std::max(IONE_, numroc_(&mkl_num_c, &options.mb, &myrow, &IZERO_, &options.nprow)); // My proc -> row of local A
    MKL_INT nc_num_c    = std::max(IONE_, numroc_(&mkl_num_c, &options.nb, &mycol, &IZERO_, &options.npcol)); // My proc -> col of local A

    // initialize local matrices
    double *X;
    // calloc initializes to 0, malloc doesn't initialize in bytes
    mat_size = nr_num_rec*nc_num_c;
    mem_try = sizeof(double) * mat_size;
    X = (double *)calloc(mat_size,sizeof(double)) ;
    // temporary unti covars
    for (int x=0; x < mat_size; x++)
    {
        X[x] = 1.0;
    }
    update_mem(X, "X", mem_try, mem_tot, ictxt);
    //
    // matrix descriptor
    MKL_INT lddX = std::max(IONE_, nr_num_rec);
    MKL_INT lddXt = std::max(IONE_, nr_num_c);
    MKL_INT descX[DESC_LEN_];
    MKL_INT descXt[DESC_LEN_];
    descinit_( descX,  &mkl_num_rec, &mkl_num_c, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddX, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }
    descinit_( descXt,  &mkl_num_c, &mkl_num_rec, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddXt, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }
    // read in the local matrices

    // initialize local matrices
    double *Z;
    // calloc initializes to 0, malloc doesn't initialize in bytes
    mat_size = nr_num_rec*nc_ind;
    mem_try = sizeof(double) * mat_size;
    Z = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(Z, "Z", mem_try, mem_tot, ictxt);

    // matrix descriptor
	MKL_INT lldZ   = std::max(IONE_, nr_num_rec);
    MKL_INT descZ[DESC_LEN_];
    descinit_( descZ,  &mkl_num_rec, &mkl_num_ind, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lldZ, &info);
    MKL_INT descZt[DESC_LEN_];
    descinit_( descZt,  &mkl_num_ind, &mkl_num_rec, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddG, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }

    for (int r = 0; r < num_rec; ++r) {
        // global 1-based indices
        MKL_INT I = static_cast<MKL_INT>(r + 1);          // records
        MKL_INT J = static_cast<MKL_INT>(rec2ind[r] + 1); // grm inds

        // which process owns (I,J)?
        MKL_INT owner_row = indxg2p_(&I, &options.mb, &myrow, &IZERO_, &options.nprow);
        MKL_INT owner_col = indxg2p_(&J, &options.nb, &mycol, &IZERO_, &options.npcol);

        if (owner_row == myrow && owner_col == mycol) {
            // compute local 1-based indices
            MKL_INT li = indxg2l_(&I, &options.mb, &myrow, &IZERO_, &options.nprow);
            MKL_INT lj = indxg2l_(&J, &options.nb, &mycol, &IZERO_, &options.npcol);
            Z[(lj - 1) * lldZ + (li - 1)] = 1.0;
        }
    }

	//blacs_barrier_(&ictxt, &CHAR_ACHAR_); // Ensure all sets are processed
	if (myrank_mpi == MPI_ROOT_PROC_) {
		std::cout << "Finished distributing Z matrix." << std::endl;
		//pheno_mapping.clear(); pheno_mapping.shrink_to_fit(); // Clear root's copy
	}

    // temporary Z storage
    double *tmpZ;
    // calloc initializes to 0, malloc doesn't initialize
    // in bytes
    mat_size = nr_num_rec*nc_ind;
    mem_try = sizeof(double) * mat_size;
    tmpZ = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(tmpZ, "tmpZ", mem_try, mem_tot, ictxt);

    std::cout << "Initialized tmpZ  matrix." << std::endl;

    // tmp storage V, IV, ident, mkl_num_rec * mkl_num_rec or V
    double *V;
    mat_size = nr_num_rec*nc_num_rec;
    mem_try = sizeof(double) * mat_size;
    V = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(V, "V", mem_try, mem_tot, ictxt);

    double *IV;
    IV = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(IV, "IV", mem_try, mem_tot, ictxt);

    double *ident;
    ident = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(ident, "ident", mem_try, mem_tot, ictxt);

    std::cout << "Initialized ident  matrix." << std::endl;

    MKL_INT lddV = std::max(IONE_,nr_num_rec);
    // matrix descriptor
    MKL_INT descV[DESC_LEN_];
    descinit_( descV,  &mkl_num_rec, &mkl_num_rec, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddV, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }
    double *PY;
    mat_size = nr_one*nc_num_rec;
    mem_try = sizeof(double) * mat_size;
    PY = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(PY, "PY", mem_try, mem_tot, ictxt);

    double *aia;
    mat_size = nr_one*nc_num_rec;
    mem_try = sizeof(double) * mat_size;
    aia = (double *)calloc(mat_size,sizeof(double));
    update_mem(aia, "aia", mem_try, mem_tot, ictxt);

    // set up map of pointers to distributed vectors aiv
    for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
    {
        // added the next two line, check TODO
        mat_size = nr_one*nc_num_rec;
        mem_try = sizeof(double) * mat_size;
        aiv[var_iter->first] = (double *)calloc(mat_size,sizeof(double));
        //if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "\t" << var_iter->first << ": " << aiv[var_iter->first] << std::endl; 
        update_mem(aiv[var_iter->first], "aiv", mem_try, mem_tot, ictxt);
    }

    std::cout << "got pointers for aiv" << std::endl;

    // initialize IVX mkl_num_rec*num_c
    double *IVX;
    // calloc initializes to 0, malloc doesn't initialize
    // in bytes
    mat_size = nr_num_rec*nc_num_c;
    mem_try = sizeof(double) * mat_size;
    IVX = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(IVX, "IVX", mem_try, mem_tot, ictxt);

    // matrix descriptor
    MKL_INT lddicc = std::max(IONE_,nr_num_c);
    MKL_INT descicc[DESC_LEN_];
    descinit_( descicc,  &mkl_num_c, &mkl_num_c, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &lddicc, &info);
    if(info != 0) {
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ERROR in descinit" << std::endl;
        blacs_gridexit_(&ictxt);
        MPI_Finalize();
    }
    double *IXIVX;
    double *icc;
    double *icc2;
    double *s;
    double *tmpwork;
    double *work;

    // initialize ixivx num_c*num_c
    // calloc initializes to 0, malloc doesn't initialize
    // in bytes
    mat_size = nr_num_c*nc_num_c;
    mem_try = sizeof(double) * mat_size;
    IXIVX = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(IXIVX, "IXIVX", mem_try, mem_tot, ictxt);

    std::cout << "initialized ixivx" << std::endl;

    if (mkl_num_c != IONE_) {
        std::cout << "Number of levels is not 1! Trying generalized inverse of xivx" << std::endl;
        // initialize icc num_c*num_c
        // calloc initializes to 0, malloc doesn't initialize
        // in bytes
        mat_size = nr_num_c*nc_num_c;
        mem_try = sizeof(double) * mat_size;
        icc = (double *)calloc(mat_size,sizeof(double)) ;
        update_mem(icc, "icc", mem_try, mem_tot, ictxt);

        //set ixivx to identity
        for (MKL_INT i=1; i<=mkl_num_c;i++){
            pdelset_(icc, &i, &i, descicc, &DONE_); 
        }

        //pdlaprnt_(&mkl_num_c, &mkl_num_c, IXIVX, &IONE_, &IONE_, descicc, &IZERO_, &IZERO_, "IXIVX", 6, toprint);

        // initialize icc2 num_c*num_c
        // calloc initializes to 0, malloc doesn't initialize
        // in bytes
        mat_size = nr_num_c*nc_num_c;
        mem_try = sizeof(double) * mat_size;
        icc2 = (double *)calloc(mat_size,sizeof(double)) ;
        update_mem(icc2, "icc2", mem_try, mem_tot, ictxt);

        // initialize s of size num_c
        // calloc initializes to 0, malloc doesn't initialize
        // in bytes
        mem_try = sizeof(double)*mkl_num_c;
        s = (double *)calloc(mkl_num_c,sizeof(double)) ;
        update_mem(s, "s", mem_try, mem_tot, ictxt);

        //SVD get work sizes
        tmpwork = (double *)calloc(10,sizeof(double)) ;
        //work = (double *)calloc(10,sizeof(double)) ;
        lwork = -1;
        //if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "starting pdgesvd" << std::endl;
        //pdgesvd_(&CHAR_VCHAR, &CHAR_VCHAR, &mkl_num_c, &mkl_num_c, icc, &IONE_, &IONE_, descicc, s, icc2, &IONE_, &IONE_, descicc, IXIVX, &IONE_, &IONE_, descicc, &dtmp, &lwork, &info);
        //if(info != MPI_SUCCESS) {
        //std::cout << "pdgesvd failed " << info << std::endl;
        //}
        pdgels_(&CHAR_NOTRANS_,
                &mkl_num_c, &mkl_num_c, &mkl_num_c,
                IXIVX, &IONE_, &IONE_, descicc,
                icc, &IONE_, &IONE_, descicc,
                tmpwork, &lwork, &info);
        if(info != MPI_SUCCESS) {
            std::cout << "pdgels failed " << info << std::endl;
        }
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "optimal lwork: " << tmpwork[0] << std::endl;
        lwork = round(tmpwork[0]);
        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "setting lwork: " << lwork << std::endl;
        //lwork = round(dtmp);
        // calloc initializes to 0, malloc doesn't initialize
        // in bytes
        mem_try = lwork*sizeof(double);
        work = (double *)calloc(lwork,sizeof(double)) ;
        update_mem(work, "work", mem_try, mem_tot, ictxt);
    }

    std::cout << "after mkl_num_c check" << std::endl;
    // initialize ipiv
    MKL_INT *ipiv;
    int ipiv_size = nr_num_rec + options.nb;
    mem_try = ipiv_size*sizeof(MKL_INT);
    ipiv = (MKL_INT *)calloc(ipiv_size,sizeof(MKL_INT)) ;
    update_mem(ipiv, "ipiv", mem_try, mem_tot, ictxt);

    std::cout << "initialized ipiv" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    MKL_INT var_size = options.variances.size();
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "number of variances: " << var_size << std::endl;
    int ai_size = var_size * (var_size + 1) / 2;
    if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "ai size: " << var_size << std::endl;

    double *IAI;
    // calloc initializes to 0, malloc doesn't initialize
    // in bytes
    mat_size = var_size * var_size;
    mem_try = sizeof(double) * mat_size;
    IAI = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(IAI, "IAI", mem_try, mem_tot, ictxt);

    double *ai;
    mat_size = ai_size;
    mem_try = sizeof(double) * mat_size;
    ai = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(ai, "ai", mem_try, mem_tot, ictxt);

    double *dltmp;
    mat_size = var_size;
    mem_try = sizeof(double) * mat_size;
    dltmp = (double *)calloc(mat_size,sizeof(double)) ;
    update_mem(dltmp, "dltmp", mem_try, mem_tot, ictxt);

    // single point descriptor
    MKL_INT descone[DESC_LEN_];
    //descinit_( descone,  &IONE_, &IONE_, &IONE_, &IONE_, &IZERO_, &IZERO_, &ictxt, &IONE_, &info);
    descinit_( descone,  &IONE_, &IONE_, &options.mb, &options.nb, &IZERO_, &IZERO_, &ictxt, &IONE_, &info);

    //pdlaprnt_(&mkl_num_ind, &mkl_num_rec, Z, &IONE_, &IONE_, descZ, &IZERO_, &IZERO_, "Z", 6, toprint);
    //pdlaprnt_(&mkl_num_ind, &mkl_num_c, X, &IONE_, &IONE_, descX, &IZERO_, &IZERO_, "X", 6, toprint);
    // Exit and finalize

    // iterations
    double time_start_iter;
    double time_end_iter;
    double trace;
    double mem_tot_gb;

    mem_tot_gb = mem_tot/1000000000;
    //MPI_Barrier(MPI_COMM_WORLD);
    //std::cout << "total memory allocated: " << mem_tot_gb << "GB" << std::endl;
    //printf("[proc_%i] allocated total %fGB\n", myrank_mpi, mem_tot_gb);

    // TODO: make average memory over all tasks?
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << std::endl
            << "memory allocated for root task: " << mem_tot_gb << "GB" << std::endl;
    }

    for (int iter=1; iter<=options.iter_n; iter++) {
        // Timing
        time_start_iter = MPI_Wtime();
        if (myrank_mpi == MPI_ROOT_PROC_) {
            std::cout << std::endl << "iteration " << iter << " starts -----------------------------" << std::endl;
            std::cout << std::endl;
        }


        // --------start get V

        // zero out V
		std::memset(V, 0, sizeof(double)*nr_num_rec*nc_num_rec);

        // get V matrix with (new) variances
        // calculate V = var_a * Z * S * Zt for var_
        // Glocal is symmetric, transposed V ok
        for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
        {
            std::string current_name = var_iter->first;
            double current_var = var_iter->second;
            // skip residual since no grm
            if (current_name == "e") continue;
            if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "calculating V = var * Z * S * Zt for v" << current_name << std::endl;

            //V =  var_a * Z * S * Z'
            pdsymm_(&CHAR_RIGHT_, &CHAR_LOWER_,
                    &mkl_num_rec, &mkl_num_ind,
                    &current_var,
                    Glocal[current_name], &IONE_, &IONE_, descG,
                    Z, &IONE_, &IONE_, descZ,
                    &DZERO_,
                    tmpZ, &IONE_, &IONE_, descZ);

            pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
                    &mkl_num_rec, &mkl_num_rec, &mkl_num_ind,
                    &DONE_,
                    tmpZ, &IONE_, &IONE_, descZ,
                    Z,    &IONE_, &IONE_, descZ,
                    &DONE_,
                    V, &IONE_, &IONE_, descV);
        } // end of variances loop

        // set iv to 0 because we want to use it as identity matrix next
        std::memset(IV, 0.0, nr_num_rec*nc_num_rec*sizeof(double));

        // add var_e to diagonal of V
		// TODO: use local matrices - faster
        for (MKL_INT i=1; i<=mkl_num_rec;i++){

            pdelget_( &CHAR_TRANS_, &CHAR_TRANS_, &dtmp, V, &i, &i, descV);
            dtmp += options.variances["e"];
            pdelset_(V, &i, &i, descV, &dtmp); 
            // set IV to identity matrix, set
            pdelset_(IV, &i, &i, descV, &DONE_); 
        }
        // --------end get V

        // end of get V matrix


        // --------start get P
        pdpotrf_(&CHAR_LOWER_, &mkl_num_rec, V, &IONE_, &IONE_, descV, &info);
        if (myrank_mpi == MPI_ROOT_PROC_) printf("pdpotrf done %i\n", (int)info);
        pdlacpy(&CHAR_TRANS_, &mkl_num_rec, &mkl_num_c, X, &IONE_, &IONE_, descX, IVX, &IONE_, &IONE_, descX);
        pdpotrs_(&CHAR_LOWER_, &mkl_num_rec, &mkl_num_c, V, &IONE_, &IONE_, descV, IVX, &IONE_, &IONE_, descX, &info);
        if (myrank_mpi == MPI_ROOT_PROC_) printf("symm solve ivx done %i\n", (int)info);
        //pdlaprnt_(&mkl_num_rec, &mkl_num_rec, IV, &IONE_, &IONE_, descV, &IZERO_, &IZERO_, "Vfactorized", 6, toprint);

        // figure out how to copy symm to full then use pdpotri, faster.
        //pdpotri_(&CHAR_LOWER_, &mkl_num_rec, IV, &IONE_, &IONE_, descV, &info);
        pdpotrs_(&CHAR_LOWER_, &mkl_num_rec, &mkl_num_rec, V, &IONE_, &IONE_, descV, IV, &IONE_, &IONE_, descV, &info);
        if (myrank_mpi == MPI_ROOT_PROC_) printf("inverse V done %i\n", (int)info);


        //pdlaprnt_(&mkl_num_rec, &mkl_num_c, IVX, &IONE_, &IONE_, descX, &IZERO_, &IZERO_, "IVX", 6, toprint);
        //pdlaprnt_(&mkl_num_rec, &mkl_num_rec, ident, &IONE_, &IONE_, descV, &IZERO_, &IZERO_, "ident", 6, toprint);

        // XIVX
        pdgemm_(&CHAR_TRANS_, &CHAR_NOTRANS_,
                &mkl_num_c, &mkl_num_c, &mkl_num_rec,
                &DONE_,
                X, &IONE_, &IONE_, descX,
                IVX, &IONE_, &IONE_, descX,
                &DZERO_,
                IXIVX, &IONE_, &IONE_, descicc);
        if (myrank_mpi == MPI_ROOT_PROC_) printf("multiply xivx done %i\n", (int)info);

        //pdlaprnt_(&mkl_num_c, &mkl_num_c, IXIVX, &IONE_, &IONE_, descicc, &IZERO_, &IZERO_, "XIVX", 6, toprint);

        // IXIVX=inverse(XIVX)
        // if no fixed effects or covariables then just use devision
        if (mkl_num_c != IONE_) {
            // maybe do pdgetrf/pdgetrs
            pdgels_(&CHAR_NOTRANS_,
                    &mkl_num_c, &mkl_num_c, &mkl_num_c,
                    IXIVX, &IONE_, &IONE_, descicc,
                    icc, &IONE_, &IONE_, descicc,
                    work, &lwork, &info);
            if(info != MPI_SUCCESS) {
                std::cout << "pdgels failed " << info << std::endl;
            }

        } else {
            IXIVX[0] = 1/IXIVX[0];
        }

        //pdlaprnt_(&mkl_num_c, &mkl_num_c, IXIVX, &IONE_, &IONE_, descicc, &IZERO_, &IZERO_, "IXIVX", 6, toprint);
        //pdlaprnt_(&mkl_num_c, &mkl_num_c, icc, &IONE_, &IONE_, descicc, &IZERO_, &IZERO_, "icc", 6, toprint);


        //need to use v for answer, bigger than required but already exists
        //(ixivx)*xiv
        pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
                &mkl_num_c, &mkl_num_rec, &mkl_num_c,
                &DONE_,
                IXIVX, &IONE_, &IONE_, descicc,
                IVX, &IONE_, &IONE_, descX,
                &DZERO_,
                V, &IONE_, &IONE_, descXt);

        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "multiply to get P "  << std::endl;
        // IV is P
        pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                &mkl_num_rec, &mkl_num_rec, &mkl_num_c,
                &DNEGONE_,
                IVX, &IONE_, &IONE_, descX,
                V, &IONE_, &IONE_, descXt,
                &DONE_,
                IV, &IONE_, &IONE_, descV);
        // --------end get P

        //pdlaprnt_(&mkl_num_rec, &mkl_num_rec, IV, &IONE_, &IONE_, descV, &IZERO_, &IZERO_, "P", 6, toprint);


        // TODO: PY is not used anymore?
        // PY REMOVE, calculate with y after to speed up but more memory
        pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                &mkl_num_rec, &IONE_, &mkl_num_rec,
                &DONE_,
                IV, &IONE_, &IONE_, descV,
                Ylocal, &IONE_, &IONE_, descY,
                &DZERO_,
                aiv["e"], &IONE_, &IONE_, descY);
        //pdlaprnt_(&mkl_num_rec, &IONE_, Y, &IONE_, &IONE_, descY, &IZERO_, &IZERO_, "Y", 6, toprint);
        //pdlaprnt_(&mkl_num_rec, &IONE_, aiv["e"], &IONE_, &IONE_, descY, &IZERO_, &IZERO_, "PY aiv["e"]", 6, toprint);

        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "getting new variances "  << std::endl;

        //----------------------------------------------------------
        // get varience
        // Z*G
        for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
        {
            std::string current_name = var_iter->first;
            double current_var = var_iter->second;
            if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "calculating new start variance for var_" << current_name << std::endl;
            // skip residual
            if (current_name == "e") continue;

            pdsymm_(&CHAR_RIGHT_, &CHAR_LOWER_,
                    &mkl_num_rec, &mkl_num_ind,
                    &DONE_,
                    Glocal[current_name], &IONE_, &IONE_, descG,
                    Z, &IONE_, &IONE_, descZ,
                    &DZERO_,
                    tmpZ, &IONE_, &IONE_, descZ);
            std::cout << "pdsymm" << current_name << std::endl;
            // Z*G*Z
            pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
                    &mkl_num_rec, &mkl_num_rec, &mkl_num_ind,
                    &DONE_,
                    tmpZ, &IONE_, &IONE_, descZ,
                    Z, &IONE_, &IONE_, descZ,
                    &DZERO_,
                    ident, &IONE_, &IONE_, descV);
            std::cout << "pdgemm1" << current_name << std::endl;


            pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                    &mkl_num_rec, &mkl_num_rec, &mkl_num_rec,
                    &DONE_,
                    IV, &IONE_, &IONE_, descV,
                    ident, &IONE_, &IONE_, descV,
                    &DZERO_,
                    V, &IONE_, &IONE_, descV);
            std::cout << "pdgemm2" << current_name << std::endl;


            // Y'*P*Z*G*Z'
            pdgemm_(&CHAR_TRANS_, &CHAR_NOTRANS_,
                    &IONE_, &mkl_num_rec, &mkl_num_rec,
                    &DONE_,
                    aiv["e"], &IONE_, &IONE_, descY,
                    ident, &IONE_, &IONE_, descV,
                    &DZERO_,
                    aiv[current_name], &IONE_, &IONE_, descYt);
            //pdlaprnt_(&IONE_, &mkl_num_rec, aiv[current_name], &IONE_, &IONE_, descYt, &IZERO_, &IZERO_, "aiv_", 6, toprint);
            std::cout << "ypzgz" << current_name << std::endl;

            // get numerator
            pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                    &IONE_, &IONE_, &mkl_num_rec,
                    &DONE_,
                    aiv[current_name], &IONE_, &IONE_, descYt,
                    aiv["e"], &IONE_, &IONE_, descY,
                    &DZERO_,
                    &dtmp, &IONE_, &IONE_, descone);
            std::cout << "numerator" << current_name << std::endl;

            MPI_Bcast(&dtmp, IONE_, MPI_DOUBLE, MPI_ROOT_PROC_, MPI_COMM_WORLD);

            trace = pdlatra_(&mkl_num_rec, V, &IONE_, &IONE_, descV);
            //if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "trace: " << trace  << std::endl;

            tmp_var[current_name] = (current_var * dtmp) / trace;
            delta_var[current_name] = fabs(tmp_var[current_name] - options.variances[current_name]);

            // for AI, removed for now
            delta_ai[current_name] = (dtmp - trace)/2.0;

        }


        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "calculating new variance for var_e" << std::endl;
        trace = pdlatra_(&mkl_num_rec, IV, &IONE_, &IONE_, descV);
        //if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "trace e: " << trace  << std::endl;

        pdgemm_(&CHAR_TRANS_, &CHAR_NOTRANS_,
                &IONE_, &IONE_, &mkl_num_rec,
                &DONE_,
                aiv["e"], &IONE_, &IONE_, descY,
                aiv["e"], &IONE_, &IONE_, descY,
                &DZERO_,
                &dtmp, &IONE_, &IONE_, descone);

        MPI_Bcast(&dtmp, IONE_, MPI_DOUBLE, MPI_ROOT_PROC_, MPI_COMM_WORLD);

        tmp_var["e"] = options.variances["e"] * dtmp / trace;
        delta_ai["e"] = (dtmp - trace)/2.0;
        delta_var["e"] = fabs(tmp_var["e"] - options.variances["e"]);

        //----------------------------------------------------------


        //----------------------------------------------------------
        // AI calculations
        int ai_idx = 0;
        double tmpdouble=0;
        if (iter >= options.ai_start && options.ai_start > 0) {
            if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "trying to use AI" << std::endl;
            for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
            {

                //char thisvar = var_iter->first;
                if (myrank_mpi == MPI_ROOT_PROC_)
                {
                    dltmp[var_idx[var_iter->first]]=delta_ai[var_iter->first];
                }

                if (var_iter->first != "e") {
                    pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                            &IONE_, &mkl_num_rec, &mkl_num_rec,
                            //&IONE_, &IONE_, &mkl_num_rec,
                            &DONE_,
                            aiv[var_iter->first], &IONE_, &IONE_, descYt,
                            //aiv[var_iter->first], &IONE_, &IONE_, descYt,
                            IV, &IONE_, &IONE_, descV,
                            &DZERO_,
                            PY, &IONE_, &IONE_, descYt);
                    //pdlaprnt_(&IONE_, &mkl_num_rec, aiv[var_iter->first], &IONE_, &IONE_, descYt, &IZERO_, &IZERO_, "aiv", 6, toprint);
                    //pdlaprnt_(&IONE_, &mkl_num_rec, PY, &IONE_, &IONE_, descYt, &IZERO_, &IZERO_, "PY", 6, toprint);
                }
                else {
                    pdgemm_(&CHAR_TRANS_, &CHAR_NOTRANS_,
                            &IONE_, &mkl_num_rec, &mkl_num_rec,
                            //&IONE_, &IONE_, &mkl_num_rec,
                            &DONE_,
                            aiv[var_iter->first], &IONE_, &IONE_, descY,
                            //aiv[var_iter->first], &IONE_, &IONE_, descYt,
                            IV, &IONE_, &IONE_, descV,
                            &DZERO_,
                            PY, &IONE_, &IONE_, descYt);
                }

                for ( int i=0; i <= var_idx[var_iter->first]; i++)
                {
                    //pddot(&mkl_num_rec, &tmpdouble, PY, &IONE_, &IONE_, descYt, &IONE_, aiv["a"], &IONE_, &IONE_, descYt, &IONE_);
                    if (i != var_idx["e"]) {
                        pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
                                &IONE_, &IONE_, &mkl_num_rec,
                                &DONE_,
                                PY, &IONE_, &IONE_, descYt,
                                aiv[idx_var[i]], &IONE_, &IONE_, descYt,
                                &DZERO_,
                                &tmpdouble, &IONE_, &IONE_, descone);
                        //pddot(&mkl_num_rec, &tmpdouble, PY, &IONE_, &IONE_, descYt, &IONE_, aiv[idx_var[i]], &IONE_, &IONE_, descYt, &IONE_);
                    }
                    else {
                        pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                                &IONE_, &IONE_, &mkl_num_rec,
                                &DONE_,
                                PY, &IONE_, &IONE_, descYt,
                                aiv[idx_var[i]], &IONE_, &IONE_, descY,
                                &DZERO_,
                                &tmpdouble, &IONE_, &IONE_, descone);
                    }
                    if (myrank_mpi == MPI_ROOT_PROC_)
                    {
                        ai[ai_idx] = tmpdouble / 2;
                    }
                    ai_idx++;

                }
            }

            if (myrank_mpi == MPI_ROOT_PROC_)
            {
                //tmpdouble = 0;
                dppsv(&CHAR_UPPER_, &var_size, &IONE_, ai, dltmp, &var_size, &info);
                std::cout << "iai info: " << info << std::endl;
                for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++) {
                    dltmp[var_idx[var_iter->first]] += var_iter->second;
                }
                tmpdouble = dltmp[0];
                for (int i=0; i<var_size;i++) {
                    if (dltmp[i] < tmpdouble) tmpdouble = dltmp[i];
                }
            }

            MPI_Bcast(&tmpdouble, IONE_, MPI_DOUBLE, MPI_ROOT_PROC_, MPI_COMM_WORLD);
            //std::cout << "proc_" << myrank_mpi << " minimum_ai: " << tmpdouble << std::endl;

            if (tmpdouble >= options.tolerance) {
                if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "using AI-REML" << std::endl;
                MPI_Bcast(&dltmp[0], var_size, MPI_DOUBLE, MPI_ROOT_PROC_, MPI_COMM_WORLD);

                /*
                   for (int i=0; i<var_size;i++) {
                   std::cout << "proc_" << myrank_mpi <<  " dltmp_" << i << ": " << dltmp[i] << std::endl;
                   }
                   */

                for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++) {
                    tmp_var[var_iter->first] = dltmp[var_idx[var_iter->first]];

                    if (myrank_mpi == MPI_ROOT_PROC_) {
                        std::cout << "\tai_variance_" << var_iter->first << ": " 
                            << std::scientific << std::setprecision(10) << tmp_var[var_iter->first]
                            << std::endl;
                    }
                }
                if (myrank_mpi == MPI_ROOT_PROC_) std::cout  << std::endl;

            }
            else {
                if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "using EM-REML" << std::endl;
            }
        }
        else {
            if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "using EM-REML" << std::endl;
        }
        //----------------------------------------------------------

        // get maximum diff
        double max_diff = 0.0;

        for (const auto & [name, diff] : delta_var) {

            if (myrank_mpi == MPI_ROOT_PROC_) {
            print_variances(name, tmp_var.at(name), diff);
            }
            max_diff = std::max(max_diff, diff);
        }


        //set new variances to options
        for (var_iter = tmp_var.begin(); var_iter != tmp_var.end(); var_iter++)
        {
            options.variances[var_iter->first] = var_iter->second;
        }

        if (myrank_mpi == MPI_ROOT_PROC_) std::cout << std::endl;


        // sum of variances
        double var_sum = std::accumulate(
                options.variances.begin(),
                options.variances.end(),
                DZERO_,
                [](double sum, const std::pair<const std::string,double>& kv) {
                return sum + kv.second;
                }
                );

        h2_new["e"] = 0;

        for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
        {
            if (var_iter->first != "e") {
                h2_new[var_iter->first] = var_iter->second / var_sum;
                h2_new["e"] += h2_new[var_iter->first];
            }
        }

        double h2_diff;
        double max_h2_diff = 0;

        for (const auto& [name, _] : options.variances) {
            // compute the difference
            h2_diff = std::fabs(h2_new[name] - h2[name]);
            max_h2_diff = std::max(max_h2_diff, h2_diff);

            if (myrank_mpi == MPI_ROOT_PROC_) {
                print_heritabilities(name, h2_new[name], h2_diff);
            }

            // update for next iteration
            h2[name] = h2_new[name];
        }

        //if (myrank_mpi == MPI_ROOT_PROC_) std::cout << "max heritability diff: " << max_h2_diff << std::endl;

        // if max_diff below heritability tolerance, exit loop
        if (max_h2_diff <= options.htolerance && options.htolerance > 0)
        {
            time_end_iter = MPI_Wtime();
            if (myrank_mpi == MPI_ROOT_PROC_) {
                print_end_iteration(iter, time_end_iter - time_start_iter);
                print_reached_tolerance("heritability", options.htolerance, iter);
            }
            break;
        }

        // if max_diff below variance tolerance, exit loop
        if (max_diff <= options.tolerance)
        {
            time_end_iter = MPI_Wtime();
            if (myrank_mpi == MPI_ROOT_PROC_) {
                print_end_iteration(iter, time_end_iter - time_start_iter);
                print_reached_tolerance("variance", options.tolerance, iter);
            }
            break;
        }


        time_end_iter = MPI_Wtime();
        //printf("[proc_%i] iteration done, time %fs.\n", myrank_mpi, time_end_iter-time_start_iter);
        if (myrank_mpi == MPI_ROOT_PROC_) {
            print_end_iteration(iter, time_end_iter - time_start_iter);
        }



    } //end of iterations

    double var_sum = std::accumulate(
            options.variances.begin(),
            options.variances.end(),
            0.0,
            [](double sum, const auto& kv) { return sum + kv.second; }
            );


    h2_new["e"] = 0;

    for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
    {
        if (var_iter->first != "e") {
            h2_new[var_iter->first] = var_iter->second / var_sum;
            h2_new["e"] += h2_new[var_iter->first];
        }
    }

    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "\n----final variances:" << std::endl;
        for (const auto& [name, variance] : options.variances) {
            print_variances(name, variance);
        }
        std::cout << "\n----final heritabilities:" << std::endl;
        for (const auto& [name, variance] : options.variances) {
            print_heritabilities(name, h2_new[name]);
        }

        // print greml file
        save_reml_results(
                options,
                h2_new);
    }

    blacs_barrier_(&ictxt, &CHAR_ACHAR_); // Ensure all sets are processed
   //  container for the gathered GEBV vectors
   std::map<std::string, std::vector<double>> u_effects_global;

	int ve_i = (int)var_idx.size();
    calculate_gblup(
        options,
        ictxt,
        mkl_num_ind,
        Z,     descZ,
        PY,    descYt,
        Glocal, descG,
        u_effects_global
    );

    if (iam == MPI_ROOT_PROC_) {
        print_gblup(options, individual_ids, rec2ind, u_effects_global);
    }

    // aiv["e"] holds the final P*Y
    // IXIVX (or icc) holds the final inv(X'*IV*X)
    // Ensure 'invXtIVX_ptr' points to the correct matrix (IXIVX or icc)
//    double* invXtIVX_ptr = (mkl_num_c == IONE_) ? IXIVX : icc; // Get the correct ptr
//
//    calculate_gblup(
//            ictxt, myrank_mpi, nprocs_mpi, options.nprow, options.npcol, myrow, mycol,
//            options, // Contains final variances
//            mkl_num_rec, mkl_num_ind, mkl_num_c,
//            individual_ids, // Make sure this is available
//            rec2ind,
//            Ylocal, descY,
//            X, descX,
//            Z, descZ,
//            Glocal, descG,
//            IV, descV,
//            aiv["e"],      // Pass P*Y (stored in aiv["e"])
//            invXtIVX_ptr, descicc, // Pass inv(X'*IV*X)
//            V, descXt,      // Pass X'*IV (stored in V)
//            var_idx, idx_var
//            );
//

    double time_end = MPI_Wtime();

    if (myrank_mpi == MPI_ROOT_PROC_) {
        printf("\nmemory per task: %f\n", mem_tot_gb);
        printf("total running time: %f\n", time_end-time_start);
    }

    // free aiv
    for (var_iter = options.variances.begin(); var_iter != options.variances.end(); var_iter++)
    {
        free(aiv[var_iter->first]);
    }

    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "freed Gs" << std::endl;
    }


    MPI_Type_free(&MPI_MKL_INT);

    free(X);
    free(Z);
    free(tmpZ);
    free(Ylocal);
    free(V);
    free(IV);
    free(ident);
    free(PY);
    free(aia);
    free(IVX);
    free(ipiv);

    if (mkl_num_c != IONE_) {
        free(icc);
        free(icc2);
        free(IXIVX);
        free(s);
        free(work);
    }

    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "freed the rest" << std::endl;
    }

    // Exit and finalize
    // calls mpi finalize internally
    blacs_gridexit_(&ictxt);
    //Cblacs_gridexit(ictxt);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "exited blac grid" << std::endl;
    }
    blacs_exit_(&IZERO_);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "exited blacs" << std::endl;
    }
    return 0;
    //return EXIT_SUCCESS;


}
