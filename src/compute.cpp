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

// Helper function to initialize local_w_segment_col_major and its descriptor
static void initialize_local_w_segment_and_descriptor(
    MKL_INT mkl_num_ind, MKL_INT n_col_local, int num_individuals_int,
    int* blocks_local, double* gstat_local,
    int gstat_offset_0, int gstat_offset_1, int gstat_offset_2,
    MKL_INT context_from_descG,
    MKL_INT myid, MKL_INT myrow, MKL_INT mycol,
    const MKL_INT empty_descriptor[DESC_LEN_],
    double*& out_local_w_segment_col_major, // Output: allocated segment
    MKL_INT out_local_w_segment_descriptor[DESC_LEN_] // Output: descriptor for the segment
) {
    MKL_INT info;

    if (n_col_local > 0) {
        MKL_INT mat_size = mkl_num_ind * n_col_local;
        out_local_w_segment_col_major = (double *)calloc(mat_size, sizeof(double));
        if (out_local_w_segment_col_major == NULL) {
            std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): memory allocation for local_w_segment_col_major failed in helper." << std::endl;
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < n_col_local; j++) {
            for (int i = 0; i < num_individuals_int; i++) {
                if (blocks_local[i * n_col_local + j] == 0.0)      out_local_w_segment_col_major[j * mkl_num_ind + i] = gstat_local[gstat_offset_0 * n_col_local + j];
                else if (blocks_local[i * n_col_local + j] == 1.0) out_local_w_segment_col_major[j * mkl_num_ind + i] = gstat_local[gstat_offset_1 * n_col_local + j];
                else if (blocks_local[i * n_col_local + j] == 2.0) out_local_w_segment_col_major[j * mkl_num_ind + i] = gstat_local[gstat_offset_2 * n_col_local + j];
                else                                              out_local_w_segment_col_major[j * mkl_num_ind + i] = 0.0;
            }
        }

        descinit_(out_local_w_segment_descriptor,
                  &mkl_num_ind, &n_col_local,
                  &mkl_num_ind, &n_col_local,
                  &myrow, &mycol,
                  &context_from_descG,
                  &mkl_num_ind,
                  &info);
        if (info != 0) {
            std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): local_w_segment_descriptor initialization failed in helper, info = " << info << std::endl;
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    } else { // n_col_local <= 0
        out_local_w_segment_col_major = (double *)calloc(1, sizeof(double));
        if (out_local_w_segment_col_major == NULL) {
            std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): memory allocation for local_w_segment_col_major (n_col_local <=0) failed in helper." << std::endl;
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        memcpy(out_local_w_segment_descriptor, empty_descriptor, DESC_LEN_ * sizeof(MKL_INT));
    }
}

// Helper function to distribute local data segments to the global gathered matrix
static void distribute_local_data_to_global_matrix(
    MKL_INT nprocs, const int* displ, const int* columns, MKL_INT mkl_num_ind,
    const MKL_INT all_procs_local_w_segment_descriptors[][DESC_LEN_],
    const double* local_w_segment_col_major,
    double* gathered_w_matrix_dist,
    const MKL_INT gathered_w_matrix_descriptor_dist[DESC_LEN_]
) {
    for (int me = 0; me < nprocs; me++) {
        if (all_procs_local_w_segment_descriptors[me][3] == 0) { // DNUMCOL is at index 3
            // Skip processes with no columns
        } else {
            MKL_INT WA_idx = displ[me] + 1;
            MKL_INT WA_cols = columns[me];
            pdgeadd_(&CHAR_NOTRANS_,
                     &mkl_num_ind, &WA_cols,
                     &DONE_,
                     local_w_segment_col_major, &IONE_, &IONE_, all_procs_local_w_segment_descriptors[me],
                     &DZERO_,
                     gathered_w_matrix_dist, &IONE_, &WA_idx, gathered_w_matrix_descriptor_dist);
        }
    }
}


static void update_g_internal(
    double* g_local,
    const int* columns, const int* displ,
    MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp,
    int* blocks_local, double* gstat_local,
    MKL_INT* descG,
    // MKL_INT ictxt, // Removed, will use descG[1]
    MKL_INT mkl_num_snp, MKL_INT mb, MKL_INT nb,
    MKL_INT lldG,
    MPI_Comm mpi_comm,
    int gstat_offset_0,
    int gstat_offset_1,
    int gstat_offset_2) {

    MKL_INT myid, nprocs, myrow, mycol, nprow, npcol;
    MKL_INT context_from_descG = descG[1]; // Get context from descG
    blacs_pinfo_(&myid, &nprocs) ; // BLACS rank and world size
    blacs_gridinfo_(&context_from_descG, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    double *gathered_w_matrix_dist; // Renamed from WA
    double *local_w_segment_col_major; // Renamed from W_tmp
    MKL_INT all_procs_local_w_segment_descriptors[nprocs][DESC_LEN_]; // Renamed from descW_tmp
    MKL_INT gathered_w_matrix_descriptor_dist[DESC_LEN_]; // Renamed from descWA
    MKL_INT local_w_segment_descriptor[DESC_LEN_]; // Renamed from descW_me
    MKL_INT empty_descriptor[DESC_LEN_]; // Renamed from descempty
    MKL_INT info;
    MKL_INT mem_try=0;

    int num_individuals_int = (int) mkl_num_ind; // Renamed from r

    // The ictxt is now derived from descG[1] via context_from_descG
    // No longer need MKL_INT ictxt = descG[1];

    memset(empty_descriptor, 0, DESC_LEN_*sizeof(MKL_INT)); // Renamed descempty

    // Step 2 & 3: Preparation of local_w_segment_col_major and its descriptor
    MKL_INT n_col_local = columns[myid];
    initialize_local_w_segment_and_descriptor(
        mkl_num_ind, n_col_local, num_individuals_int,
        blocks_local, gstat_local,
        gstat_offset_0, gstat_offset_1, gstat_offset_2,
        context_from_descG, myid, myrow, mycol,
        empty_descriptor,
        local_w_segment_col_major, // Output
        local_w_segment_descriptor // Output
    );

    // Step 4: gather all descriptors
    int mpi_error_code = MPI_Allgather(local_w_segment_descriptor,
            DESC_LEN_,
            MPI_LONG_LONG,
            &all_procs_local_w_segment_descriptors, // Renamed descW_tmp
            DESC_LEN_,
            MPI_LONG_LONG,
            MPI_COMM_WORLD);
    if (mpi_error_code != MPI_SUCCESS) {
        std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): MPI_Allgather failed with error code " << mpi_error_code << std::endl;
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Step 5: initialize local matrices (for gathered_w_matrix_dist)
    MKL_INT mat_size = nr_ind * nc_snp; // mat_size reused, but context is clear
    mem_try = sizeof(double) * mat_size;

    if (myid == MPI_ROOT_PROC_) {
        //std::cout << "w_cols " << col_hap_count_total[col] << std::endl;
        //std::cout << "nr_ind " << nr_ind << " * nc_snp " << nc_snp << " = "
        //    << mat_size << " * " << sizeof(double) << " " << mem_try << " bytes" << std::endl;
        //std::cout << "trying to allocate W local: " << mem_try/1000000.0 << "MB" << std::endl;
    }

    gathered_w_matrix_dist = (double *)calloc(mat_size, sizeof(double)) ; // Renamed WA
    if (gathered_w_matrix_dist==NULL){ // Renamed WA
        std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): memory allocation for gathered_w_matrix_dist failed." << std::endl; // Renamed WA
        MPI_Finalize();
        exit(EXIT_FAILURE);
    } else {
        // The mem_tot related TODO and commented line are removed as mem_tot is not used.
    }
    // matrix descriptor
    descinit_( gathered_w_matrix_descriptor_dist,  &mkl_num_ind, &mkl_num_snp, &mb, &nb, &IZERO_, &IZERO_, &context_from_descG, &lldG, &info); // Renamed descWA
    if(info != 0) {
        std::cerr << "ERROR on proc " << myid << " (" << myrow << "," << mycol << "): gathered_w_matrix_descriptor_dist initialization failed, info = " << info << std::endl;
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }


    // Step 6: distribute local data to the global gathered matrix
    distribute_local_data_to_global_matrix(
        nprocs, displ, columns, mkl_num_ind,
        all_procs_local_w_segment_descriptors,
        local_w_segment_col_major,
        gathered_w_matrix_dist,
        gathered_w_matrix_descriptor_dist
    );

    free(local_w_segment_col_major); 

    // Step 7: multiply WW'
    pdgemm_(&CHAR_NOTRANS_, &CHAR_TRANS_,
            &mkl_num_ind, &mkl_num_ind, &mkl_num_snp,
            &DONE_,
            gathered_w_matrix_dist,   &IONE_, &IONE_, gathered_w_matrix_descriptor_dist, // Renamed WA, descWA
            gathered_w_matrix_dist,   &IONE_, &IONE_, gathered_w_matrix_descriptor_dist, // Renamed WA, descWA
            &DONE_,
            g_local, &IONE_, &IONE_, descG);

    free(gathered_w_matrix_dist); // Renamed WA
}

void update_gd(double* g_local, 
        const int* columns, const int* displ, 
        MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp, 
        int* blocks_local, double* gstat_local, 
        MKL_INT* descG, 
        // MKL_INT ictxt, // Removed
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
