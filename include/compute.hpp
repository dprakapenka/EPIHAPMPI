#ifndef COMPUTE_HPP
#define COMPUTE_HPP

#include <mpi.h>
#include "mkl.h"

// Updates the G matrix component for additive genetic effects.
void update_ga(double* g_local,
        const int* columns, const int* displ,
        MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp,
        int* blocks_local, double* gstat_local,
        MKL_INT* descG,
        MKL_INT mkl_num_snp, MKL_INT mb, MKL_INT nb,
        MKL_INT lldG,
        MPI_Comm mpi_comm);

// Updates the G matrix component for dominance genetic effects.
void update_gd(double* g_local,
        const int* columns, const int* displ,
        MKL_INT mkl_num_ind, MKL_INT nr_ind, MKL_INT nc_snp,
        int* blocks_local, double* gstat_local,
        MKL_INT* descG,
        MKL_INT mkl_num_snp, MKL_INT mb, MKL_INT nb,
        MKL_INT lldG,
        MPI_Comm mpi_comm);

#endif // COMPUTE_HPP
