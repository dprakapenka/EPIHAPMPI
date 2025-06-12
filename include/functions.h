#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "errors.h"

#include <vector>
#include <string>
#include <map>
#include <utility>  // For std::pair
#include "mkl.h"
#include "mpi.h"

//
void check_block_size(
        int& b,
        int n,
        int np
        );


// suggests the most square blacs grid give nprocs_mpi tasks
std::pair<MKL_INT, MKL_INT> compute_blacs_grid(
        int nprocs_mpi
        );

// If matrix==nullptr, calls allocation_error.
// otherwise adds memory required in bytes to mem_tot.
void update_mem(void*            matrix,
                const std::string& name,
                MKL_INT          mem_req,
                MKL_INT          mem_tot,
                MKL_INT          ictxt);

#endif // FUNCTIONS_H
