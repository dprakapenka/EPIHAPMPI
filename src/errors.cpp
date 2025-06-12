#include "errors.h"
#include "constants.h"    // for blacs_gridexit_
#include <iostream>
#include <cstdlib>        // for std::exit
#include "mpi.h"
#include "mkl_blacs.h"

// print info on allocation error
// TODO: add memory information
void allocation_error(const std::string& name, MKL_INT mem_req, MKL_INT ictxt)
{
    // get MPI rank
    int iam;
    MPI_Comm_rank(MPI_COMM_WORLD, &iam);

    // get BLACS grid coords
    MKL_INT nprow, npcol, myrow, mycol;
    blacs_gridinfo(&ictxt, &nprow, &npcol, &myrow, &mycol);

    std::cout 
        << "ERROR: failed memory allocation of "
        << mem_req
        << "ERROR: memory allocation for "
        << name
        << " on proc " << iam
        << ": [" << myrow
        << "," << mycol << "]"
        << std::endl;

    blacs_gridexit_(&ictxt);
    blacs_gridexit_(&ictxt);
    blacs_exit_(&IZERO_);
    MPI_Finalize();
    std::exit(1);
}
