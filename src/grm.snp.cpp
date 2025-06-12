#include "constants.h"
#include "grm.hpp"
#include "read.h"
#include "compute.hpp"

#include <cmath>
#include <float.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>
#include <cassert>
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
    int limit = trunc(per_proc);
    if (b > limit) {
        b = limit;
        std::cout << "changed number of blocks to " << limit << std::endl;
    }

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
    int proc_dims[ndims] = {(int) nprow,(int)  npcol};

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


void update_ga(double* g_local, 
        const int* columns, const int* displ, 
        MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp, 
        int* blocks_local, double* gstat_local, 
        MKL_INT* descG, 
        MKL_INT ictxt, 
        MKL_INT mkl_num_snp, MKL_INT mb, MKL_INT nb, 
        MKL_INT lldG, 
        MPI_Comm mpi_comm) {

    MKL_INT myid, nprocs, myrow, mycol, nprow, npcol;
    blacs_pinfo_(&myid, &nprocs) ; // BLACS rank and world size
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    double *WA;
    double *W_tmp;
    MKL_INT descW_tmp[nprocs][DESC_LEN_];
    MKL_INT descWA[DESC_LEN_];
    MKL_INT descW_me[DESC_LEN_];
    MKL_INT descempty[DESC_LEN_];
    MKL_INT info;
    MKL_INT mem_try=0;

    int r = (int) mkl_num_ind;


    //TODO: can get this from descG, check
    //MKL_INT ictxt = descG[1];

    memset(descempty, 0, DESC_LEN_*sizeof(MKL_INT));

    MKL_INT mat_size = mkl_num_ind * columns[myid];
    W_tmp = (double *)calloc(mat_size, sizeof(double));

    // switch to col major for scalapack on w_tmp?
    for(int j=0;j<columns[myid];j++){
        for(int i=0;i<r;i++){
            if(blocks_local[i * columns[myid] + j]==0.0)      W_tmp[j * mkl_num_ind + i] = gstat_local[2 * columns[myid] + j];
            else if(blocks_local[i * columns[myid] + j]==1.0) W_tmp[j * mkl_num_ind + i] = gstat_local[3 * columns[myid] + j];
            else if(blocks_local[i * columns[myid] + j]==2.0) W_tmp[j * mkl_num_ind + i] = gstat_local[4 * columns[myid] + j];
            else                                              W_tmp[j * mkl_num_ind + i] = 0.0;
        }
    }


    // w hap gather
    MKL_INT  n_col = columns[myid];

    //MKL_INT tmpcol;
    if (n_col > 0) {
        descinit_( descW_me,
                &mkl_num_ind, &n_col,
                &mkl_num_ind, &n_col,
                &myrow, &mycol,
                &ictxt,
                &mkl_num_ind,
                &info);
        if(info != 0) {
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

    // initialize local matrices
    mat_size = nr_ind * nc_snp;
    mem_try = sizeof(double) * mat_size;

    if (myid == MPI_ROOT_PROC_) {
        //std::cout << "w_cols " << col_hap_count_total[col] << std::endl;
        //std::cout << "nr_ind " << nr_ind << " * nc_snp " << nc_snp << " = "
        //    << mat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
        //std::cout << "trying to allocate W local: " << mem_try/1000000.0 << "MB" << std::endl;
    }

    WA = (double *)calloc(mat_size, sizeof(double)) ;
    if (WA==NULL){
        std::cout << "ERROR: memory allocation for WA on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
        MPI_Finalize();
        exit(0);
    } else {
        //TODO: mem_tot is not useful here since we delete wa, change
        //mem_tot += mem_try;
    }
    // matrix descriptor
    descinit_( descWA,  &mkl_num_ind, &mkl_num_snp, &mb, &nb, &IZERO_, &IZERO_, &ictxt, &lldG, &info);
    if(info != 0) {
        printf("Error in descinit WA, info = %d\n", int(info));
    }


    // distribute to WA
    for (int me=0; me < nprocs; me++) {
        // tasks with no column do not participate
        if (descW_tmp[me][3] == 0) {
            //std::cout << "EMPTY DESCRIPTOR ON : " << me << std::endl;
        }
        else {
            MKL_INT WA_idx = displ[me]+1;
            MKL_INT WA_cols = columns[me];
            // std::cout << "task " << myid << " distributing mat from " << me << 
            //     " has big index: " << WA_idx << " of size: " << col_hap_count[col][me] << std::endl;
            pdgeadd_(&CHAR_NOTRANS_,
                    &mkl_num_ind, &WA_cols,
                    &DONE_,
                    W_tmp,   &IONE_, &IONE_, descW_tmp[me],
                    &DZERO_,
                    WA, &IONE_, &WA_idx, descWA);
        }
    }

    free(W_tmp);

    // multiply WW'
    pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
            &mkl_num_ind, &mkl_num_ind, &mkl_num_snp,
            &DONE_,
            WA,   &IONE_, &IONE_, descWA,
            WA,   &IONE_, &IONE_, descWA,
            &DONE_,
            g_local, &IONE_, &IONE_, descG);

    free(WA);
}

void update_gd(double* g_local, 
        const int* columns, const int* displ, 
        MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp, 
        int* blocks_local, double* gstat_local, 
        MKL_INT* descG, 
        MKL_INT ictxt, 
        MKL_INT mkl_num_snp, MKL_INT mb, MKL_INT nb, 
        MKL_INT lldG, 
        MPI_Comm mpi_comm) {

    MKL_INT myid, nprocs, myrow, mycol, nprow, npcol;
    blacs_pinfo_(&myid, &nprocs) ; // BLACS rank and world size
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    double *WA;
    double *W_tmp;
    MKL_INT descW_tmp[nprocs][DESC_LEN_];
    MKL_INT descWA[DESC_LEN_];
    MKL_INT descW_me[DESC_LEN_];
    MKL_INT descempty[DESC_LEN_];
    MKL_INT info;
    MKL_INT mem_try=0;

    int r = (int) mkl_num_ind;

    //TODO: can get this from descG, check
    //MKL_INT ictxt = descG[1];

    memset(descempty, 0, DESC_LEN_*sizeof(MKL_INT));

    MKL_INT mat_size = mkl_num_ind * columns[myid];
    W_tmp = (double *)calloc(mat_size, sizeof(double));

    // switch to col major for scalapack on w_tmp?
    for(int j=0;j<columns[myid];j++){
        for(int i=0;i<r;i++){
            if(blocks_local[i * columns[myid] + j]==0.0)      W_tmp[j * mkl_num_ind + i] = gstat_local[5 * columns[myid] + j];
            else if(blocks_local[i * columns[myid] + j]==1.0) W_tmp[j * mkl_num_ind + i] = gstat_local[6 * columns[myid] + j];
            else if(blocks_local[i * columns[myid] + j]==2.0) W_tmp[j * mkl_num_ind + i] = gstat_local[7 * columns[myid] + j];
            else                                              W_tmp[j * mkl_num_ind + i] = 0.0;
        }
    }


    // w hap gather
    MKL_INT  n_col = columns[myid];

    //MKL_INT tmpcol;
    if (n_col > 0) {
        descinit_( descW_me,
                &mkl_num_ind, &n_col,
                &mkl_num_ind, &n_col,
                &myrow, &mycol,
                &ictxt,
                &mkl_num_ind,
                &info);
        if(info != 0) {
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

    // initialize local matrices
    mat_size = nr_ind * nc_snp;
    mem_try = sizeof(double) * mat_size;

    if (myid == MPI_ROOT_PROC_) {
        //std::cout << "w_cols " << col_hap_count_total[col] << std::endl;
        //std::cout << "nr_ind " << nr_ind << " * nc_snp " << nc_snp << " = "
        //    << mat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
        //std::cout << "trying to allocate W local: " << mem_try/1000000.0 << "MB" << std::endl;
    }

    WA = (double *)calloc(mat_size, sizeof(double)) ;
    if (WA==NULL){
        std::cout << "ERROR: memory allocation for WA on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
        MPI_Finalize();
        exit(0);
    } else {
        //TODO: mem_tot is not useful here since we delete wa, change
        //mem_tot += mem_try;
    }
    // matrix descriptor
    descinit_( descWA,  &mkl_num_ind, &mkl_num_snp, &mb, &nb, &IZERO_, &IZERO_, &ictxt, &lldG, &info);
    if(info != 0) {
        printf("Error in descinit WA, info = %d\n", int(info));
    }

    // distribute to WA
    for (int me=0; me < nprocs; me++) {
        // tasks with no column do not participate
        if (descW_tmp[me][3] == 0) {
            //std::cout << "EMPTY DESCRIPTOR ON : " << me << std::endl;
        }
        else {
            MKL_INT WA_idx = displ[me]+1;
            MKL_INT WA_cols = columns[me];
            // std::cout << "task " << myid << " distributing mat from " << me << 
            //     " has big index: " << WA_idx << " of size: " << col_hap_count[col][me] << std::endl;
            pdgeadd_(&CHAR_NOTRANS_,
                    &mkl_num_ind, &WA_cols,
                    &DONE_,
                    W_tmp,   &IONE_, &IONE_, descW_tmp[me],
                    &DZERO_,
                    WA, &IONE_, &WA_idx, descWA);
        }
    }

    free(W_tmp);

    // multiply WW'
    pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
            &mkl_num_ind, &mkl_num_ind, &mkl_num_snp,
            &DONE_,
            WA,   &IONE_, &IONE_, descWA,
            WA,   &IONE_, &IONE_, descWA,
            &DONE_,
            g_local, &IONE_, &IONE_, descG);

    free(WA);
}

int main(int argc, char **argv) {

    // Initialize MPI
    int myrank_mpi, nprocs_mpi;
    int num_threads;
    MKL_INT info;
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

    // temp vars
    int itmp;
    double dtmp;

    // matrices
    int *blocks_local;
    int *buffer_local;
    double *gstat_local;

    // timing
    double time_start;
    double time_end;
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

    //the size of row and column ints in the saved binary
    int skip = 16;

    //defaults
    std::string out_filename = "grm";
    std::string final_filename;
    MKL_INT mb = 2;
    MKL_INT nb = 2;
    bool needCopy = true;
    bool secondOrder = true;
    bool thirdOrder = true;
    bool needAdditive = true;
    bool needDominance = true;


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
            case 'o': //output data (save)
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
    read_file_size(chr_filenames[0].c_str(), r, c);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "found individuals: " << r << std::endl;
    }

    double *toprint;
    toprint = (double *)calloc(nb,sizeof(double)) ;

    // set up G matrix size
    MKL_INT mkl_num_ind = r;
    MKL_INT mkl_num_snp = c;

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
    MKL_INT nc_snp    = std::max(IONE_, numroc_(&mkl_num_snp, &nb, &mycol, &IZERO_, &npcol)); // My proc -> col of local A
    MKL_INT lldG = std::max(IONE_,nr_ind);


    // initialize local matrices
    MKL_INT g_local_size = nr_ind * nc_ind;
    //int g_local_size_int = (int)g_local_size;
    mem_try = sizeof(double) * g_local_size;

    if (myid == MPI_ROOT_PROC_) {
        std::cout << "---nr_ind " << nr_ind << " * nc_ind " << nc_ind << " = "
            << " mklint(" << g_local_size << ")" << " * " << sizeof(double)
            << " " << mem_try << " bytes" << std::endl;
        std::cout << "mkl_nr_ind " << nr_ind << " mkl_nc_ind " << nc_ind << " mkl_num_ind " << mkl_num_ind << std::endl;
        std::cout << "trying to allocate ga local: " << mem_try/1000000.0 << "MB" << std::endl;
    }

    // ga_local
    double *ga_local;
    if (needAdditive) {
        // switch to mkl malloc aligned allocation
        //ga_local = (double *)mkl_malloc(g_local_size,  64) ;
        ga_local = (double *)calloc(g_local_size, sizeof(double)) ;
        if (ga_local==NULL){
            std::cout << "ERROR: memory allocation for ga_local on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
            MPI_Finalize();
            exit(0);
        } else {
            if (myid == MPI_ROOT_PROC_) std::cout << "...success" << std::endl;
            mem_tot += mem_try;
        }
    }
    double *gd_local;
    if (needDominance) {
        // switch to mkl malloc aligned allocation
        //gd_local = (double *)mkl_malloc(g_local_size,  64) ;
        gd_local = (double *)calloc(g_local_size, sizeof(double)) ;
        if (gd_local==NULL){
            std::cout << "ERROR: memory allocation for gd_local on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
            MPI_Finalize();
            exit(0);
        } else {
            if (myid == MPI_ROOT_PROC_) std::cout << "...success" << std::endl;
            mem_tot += mem_try;
        }
    }

    // matrix descriptor for G
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
        read_file_size(chr_filename.c_str(), r, c);
        if (myid == MPI_ROOT_PROC_) {
            std::cout << "r: " << r << " c: " << c << std::endl;
        }
        if (mkl_num_ind != (MKL_INT) r) {
            std::cout << "ERROR: mismatch in number of inds in " << chr_filename << std::endl;
        }

        int using_mpi_tasks = nprocs_mpi;
        if (c < nprocs_mpi) {
            using_mpi_tasks = c;
            if (myid == MPI_ROOT_PROC_) {
                std::cout << "warning: number of tasks is more than number of snps in chr," << 
                    "using " << using_mpi_tasks << " tasks for reading and snp statistics" << std::endl;
            }
        }

        int columns [nprocs_mpi];
        memset(columns, 0, nprocs_mpi*sizeof(int));
        int displ [nprocs_mpi];
        int remainder = 0;
        remainder = c % using_mpi_tasks;
        int max_cols = 0;
        // CHANGED to MKL_INT
        //int block_size = 0;
        MKL_INT block_size = 0;

        for ( int i = using_mpi_tasks-1; i >=0; --i) {
            columns[i] = c / using_mpi_tasks;
            if (remainder > 0) {
                columns[i]++;
                remainder--;
            }
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



        //allocate geno_stat as gstat
        int gstat_size = 12 * columns[myid];
        //std::vector<int> my_block ;
        mem_try = sizeof(double) * gstat_size;
        for (int mrank=0; mrank < nprocs_mpi; mrank++){
            if (myid == mrank) {
                std::cout << "proc " << mrank << " gstat_local: " << gstat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
            }
        }
        if (myid == MPI_ROOT_PROC_) std::cout << "trying to allocate gstat_local: " << mem_try/1000000 << "MB" << std::endl;
        //gstat_local = (int *)mkl_malloc(gstat_size,  64) ;
        gstat_local = (double *)calloc(gstat_size, sizeof(double)) ;
        if (gstat_local==NULL){
            std::cout << "Error of memory allocation gstat_local on proc " << myid << std::endl;
            MPI_Finalize();
            exit(0);
        } else {
            //std::cout << "gstat_local allocated on " << myid << std::endl;
            mem_tot += mem_try;
        }

        //allocate local blocks
        //block size may need MKL_INT
        block_size = r * columns[myid];
        //std::vector<int> my_block ;
        mem_try = sizeof(int) * block_size;
        for (int mrank=0; mrank < nprocs_mpi; mrank++){
            if (myid == mrank) {
                std::cout << "proc " << mrank << " block_local: " << block_size << " * " << sizeof(int) << " " << mem_try << " bytes" << std::endl;
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


            //MPI_Scatterv(blocks_local, columns, displ, MPI_INT, buffer_local, max_cols, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
            MPI_Scatterv(&words[0], columns, displ, MPI_INT, buffer_local, max_cols, MPI_INT, MPI_ROOT_PROC_, MPI_COMM_WORLD);
            memcpy(&blocks_local[row*columns[myid]],buffer_local,sizeof(int)*columns[myid]);

            for(int col=0; col<columns[myid]; col++){
                //double geno_tmp = (double)blocks_local[row*columns[myid]+col];
                double geno_tmp = 1.0 * blocks_local[row*columns[myid]+col];
                //if (myid == MPI_ROOT_PROC_) std::cout << "geno_tmp:::" << geno_tmp << std::endl;
                if(geno_tmp==0.0){
                    gstat_local[0 * columns[myid] + col] += 1.0;
                    gstat_local[1 * columns[myid] + col] += geno_tmp;
                    gstat_local[9 * columns[myid] + col] += 1.0;
                }
                else if(geno_tmp==1.0){
                    gstat_local[0 * columns[myid] + col] += 1.0;
                    gstat_local[1 * columns[myid] + col] += geno_tmp;
                    gstat_local[8 * columns[myid] + col] += 1.0;
                }
                else if(geno_tmp==2.0){
                    gstat_local[0  * columns[myid] + col] += 1.0;
                    gstat_local[1  * columns[myid] + col] += geno_tmp;
                    gstat_local[10 * columns[myid] + col] += 1.0;
                }
                else{
                    //if (myid == MPI_ROOT_PROC_) std::cout << "geno_tmp[" << row << "," << col << "]: " << geno_tmp << " " << blocks_local[row*columns[MPI_ROOT_PROC_] + col] << std:: endl;
                    gstat_local[11*columns[myid] + col] += 1.0;
                }
            }



        }
        infile.close();
        free(buffer_local);

        time_end_read = MPI_Wtime(); // end of read data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished reading chr, time %fs.\n", myid, time_end_read-time_start_read);
        time_start_stat = MPI_Wtime();

        for(int j=0;j<columns[myid];j++){
            gstat_local[1 * columns[myid] + j] = gstat_local[1 * columns[myid] + j]/(2*gstat_local[0 * columns[myid] + j]);
            double freq = gstat_local[1 * columns[myid] + j];
            gstat_local[2 * columns[myid] + j] = 2*freq;
            gstat_local[3 * columns[myid] + j] = 2*freq-1;
            gstat_local[4 * columns[myid] + j] = 2*freq-2;
            gstat_local[5 * columns[myid] + j] = -2*freq*freq;
            gstat_local[6 * columns[myid] + j] = 2*freq*(1-freq);
            gstat_local[7 * columns[myid] + j] = -2*(1-freq)*(1-freq);
            gstat_local[9 * columns[myid] + j] /= gstat_local[0 * columns[myid] + j];
            gstat_local[8 * columns[myid] + j] /= gstat_local[0 * columns[myid] + j];
            gstat_local[10 * columns[myid] + j] /= gstat_local[0 * columns[myid] + j];
        }



        time_end_stat = MPI_Wtime(); // end of stat data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished calculating statistics chr, time %fs.\n", myid, time_end_stat-time_start_stat);




        // compute W, multiply and add to G
        time_start_dist = MPI_Wtime(); // start of dist data
        if (needAdditive) {
            update_ga(ga_local,
                    columns, displ,
                    mkl_num_ind, nr_ind, nc_snp, 
                    blocks_local, gstat_local, 
                    descG, 
                    ictxt, 
                    mkl_num_snp, mb, nb, 
                    lldG,
                    MPI_COMM_WORLD);
        }
        time_end_dist = MPI_Wtime(); // end of dist data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished WAT, time %fs.\n", myid, time_end_dist-time_start_dist);
        time_start_dist = MPI_Wtime(); // start of dist data
        if (needDominance) {
            update_gd(gd_local,
                    columns, displ,
                    mkl_num_ind, nr_ind, nc_snp, 
                    blocks_local, gstat_local, 
                    descG, 
                    ictxt, 
                    mkl_num_snp, mb, nb, 
                    lldG,
                    MPI_COMM_WORLD);
        }


        free(blocks_local);
        free(gstat_local);

        time_end_dist = MPI_Wtime(); // end of dist data
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished WDT, time %fs.\n", myid, time_end_dist-time_start_dist);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished chr %s, time %fs.\n", myid, chr_filename.c_str(), time_end_dist-time_start_read);

    } // end of chr iteration


    double trace;

    double *out_local;
    if (needCopy) {
        // make output matrix
        // switch to mkl malloc aligned allocation
        //out_local = (double *)mkl_malloc(g_local_size,  64) ;
        out_local = (double *)calloc(g_local_size, sizeof(double)) ;
        if (out_local==NULL){
            std::cout << "ERROR: memory allocation for out_local on proc " << myid << ": [" << myrow << "," << mycol << "]" << std::endl;
            exit(0);
        } else {
            if (myid == MPI_ROOT_PROC_) std::cout << "...success" << std::endl;
            mem_tot += mem_try;
        }
    }

    if (secondOrder) {

        // ----------------------------------------------------------------------------------output AD
        // hadamard product
#pragma omp parallel for
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = ga_local[element] * gd_local[element] ;

        }

        // mean of diag
        trace = pdlatra_(&mkl_num_ind, out_local, &IONE_, &IONE_, descG);
        trace = trace / mkl_num_ind;
        if (myid == 0) std::cout << "trace: " << trace << std::endl;
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = out_local[element] / trace ;

        }
        if (myid == 0) std::cout << "divide by trace completed " << std::endl;

        time_start_write = MPI_Wtime(); // start of write
        final_filename = out_filename + ".g.AD";
        if (myid == 0) {
            std::cout << "writing binary grm to " << final_filename << std::endl;
        }
        write_binmat_size_16(final_filename.c_str(), mkl_num_ind, mkl_num_ind);
        write_binmat(final_filename.c_str(), out_local, descG, skip, g_local_size);

        time_end_write = MPI_Wtime(); // end of write
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing %s, time %fs.\n", myid,final_filename.c_str(), time_end_write - time_start_write);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);

        // ----------------------------------------------------------------------------------output AA
        // hadamard product
#pragma omp parallel for
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            //out_local[element] = std::pow(ga_local[element],2);
            out_local[element] = ga_local[element] * ga_local[element];

        }

        // mean of diag
        trace = pdlatra_(&mkl_num_ind, out_local, &IONE_, &IONE_, descG);
        trace = trace / mkl_num_ind;
        if (myid == 0) std::cout << "trace: " << trace << std::endl;
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = out_local[element] / trace ;

        }
        if (myid == 0) std::cout << "divide by trace completed " << std::endl;

        time_start_write = MPI_Wtime(); // start of write
        final_filename = out_filename + ".g.AA";
        if (myid == 0) {
            std::cout << "writing binary grm to " << final_filename << std::endl;
        }
        write_binmat_size_16(final_filename.c_str(), mkl_num_ind, mkl_num_ind);
        write_binmat(final_filename.c_str(), out_local, descG, skip, g_local_size);

        time_end_write = MPI_Wtime(); // end of write
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing %s, time %fs.\n", myid,final_filename.c_str(), time_end_write - time_start_write);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);

        // ----------------------------------------------------------------------------------output DD
        // hadamard product
#pragma omp parallel for
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            //out_local[element] = std::pow(gd_local[element],2);
            out_local[element] = gd_local[element] * gd_local[element];

        }

        // mean of diag
        trace = pdlatra_(&mkl_num_ind, out_local, &IONE_, &IONE_, descG);
        trace = trace / mkl_num_ind;
        if (myid == 0) std::cout << "trace: " << trace << std::endl;
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = out_local[element] / trace ;

        }
        if (myid == 0) std::cout << "divide by trace completed " << std::endl;

        time_start_write = MPI_Wtime(); // start of write
        final_filename = out_filename + ".g.DD";
        if (myid == 0) {
            std::cout << "writing binary grm to " << final_filename << std::endl;
        }
        write_binmat_size_16(final_filename.c_str(), mkl_num_ind, mkl_num_ind);
        write_binmat(final_filename.c_str(), out_local, descG, skip, g_local_size);

        time_end_write = MPI_Wtime(); // end of write
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing %s, time %fs.\n", myid,final_filename.c_str(), time_end_write - time_start_write);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        
    }

    //TODO: if only need a or d don't copy
    //TODO: if only need a or d don't copy
    if (needAdditive) {
        // ----------------------------------------------------------------------------------output A
        // mean of diag
        trace = pdlatra_(&mkl_num_ind, ga_local, &IONE_, &IONE_, descG);
        trace = trace / mkl_num_ind;
        if (myid == 0) std::cout << "trace: " << trace << std::endl;
#pragma omp parallel for
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = ga_local[element] / trace ;

        }
        if (myid == 0) std::cout << "divide by trace completed " << std::endl;
        //pdlaprnt_(&mkl_num_ind, &mkl_num_ind, G, &IONE_, &IONE_, descG, &IZERO_, &IZERO_, "G", 6, toprint);

        time_start_write = MPI_Wtime(); // start of write
        final_filename = out_filename + ".g.A";
        if (myid == 0) {
            std::cout << "writing binary grm to " << final_filename << std::endl;
        }
        write_binmat_size_16(final_filename.c_str(), mkl_num_ind, mkl_num_ind);
        write_binmat(final_filename.c_str(), out_local, descG, skip, g_local_size);

        time_end_write = MPI_Wtime(); // end of write
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing %s, time %fs.\n", myid,final_filename.c_str(), time_end_write - time_start_write);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
    }
    if (needDominance) {
        // ----------------------------------------------------------------------------------output D
        // mean of diag
        trace = pdlatra_(&mkl_num_ind, gd_local, &IONE_, &IONE_, descG);
        trace = trace / mkl_num_ind;
        if (myid == 0) std::cout << "trace: " << trace << std::endl;
#pragma omp parallel for
        for (MKL_INT element = 0; element < g_local_size; element++)
        {
            out_local[element] = gd_local[element] / trace ;

        }
        if (myid == 0) std::cout << "divide by trace completed " << std::endl;
        //pdlaprnt_(&mkl_num_ind, &mkl_num_ind, G, &IONE_, &IONE_, descG, &IZERO_, &IZERO_, "G", 6, toprint);

        time_start_write = MPI_Wtime(); // start of write
        final_filename = out_filename + ".g.D";
        if (myid == 0) {
            std::cout << "writing binary grm to " << final_filename << std::endl;
        }
        write_binmat_size_16(final_filename.c_str(), mkl_num_ind, mkl_num_ind);
        write_binmat(final_filename.c_str(), out_local, descG, skip, g_local_size);

        time_end_write = MPI_Wtime(); // end of write
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] finished writing %s, time %fs.\n", myid,final_filename.c_str(), time_end_write - time_start_write);
        if (myid == MPI_ROOT_PROC_) printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
        printf("[proc_%lli] running time %fs.\n", myid, time_end_write - time_start);
    }


    free(ga_local);
    free(gd_local);
    free(out_local);

    //TODO: fix timing
    time_end = MPI_Wtime();
    printf("[proc_%lli] iteration done, time %fs.\n", myid, time_end-time_start);

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
