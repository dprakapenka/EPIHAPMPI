#include "functions.h"  // Function declarations for utility functions
#include "errors.h"
#include "constants.h"

// std C++ headers
#include <vector>      // Used for handling dynamic arrays
#include <string>      // String handling and parsing
#include <iostream>    // Input/output operations
#include <fstream>     // File input/output handling
#include <map>        // Data storage using key-value pairs

// C std Library functions
#include <cmath>       // for 'trunc'

// External libraries for matrix operations
#include "mkl.h"
#include "mkl_blacs.h"
#include "mpi.h"

// check that blocks are not larger than matrix/proc
void check_block_size(int& b, int n, int np) {
    double per_proc = n/np;
    int limit = trunc(per_proc);
    if (b > limit) {
        b = limit;
        std::cout << "changed number of blocks to " << limit << std::endl;
    }
}

// Function to compute the most square BLACS 2D grid
std::pair<MKL_INT, MKL_INT> compute_blacs_grid(int nprocs_mpi) {
    const int ndims = 2;
    int dim_arr[ndims] = { 0 };
    int info = MPI_Dims_create(nprocs_mpi, ndims, dim_arr);
    if (info != MPI_SUCCESS) {
        std::cerr << "MPI_Dims_create failed" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, info);
    }
    return {dim_arr[0], dim_arr[1]};
}

void update_mem(void*            matrix,
                const std::string& name,
                MKL_INT          mem_req,
                MKL_INT          mem_tot,
                MKL_INT          ictxt)
{
    if (matrix == nullptr) {
        allocation_error(name, mem_req, ictxt);
    }
    mem_tot += mem_req;
}
