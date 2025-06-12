#ifndef REML_H
#define REML_H

#include <map>
#include <vector>
#include "functions.h"
#include "mkl.h"
#include "mpi.h"

extern "C" {
    void pdelget_(char*, char*, double*, double*, MKL_INT*, MKL_INT*, MKL_INT*);
    void pdelset_(double*, MKL_INT*, MKL_INT*, MKL_INT*, double*);
    void pdlaprnt_(MKL_INT*, MKL_INT*, double*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, char*, MKL_INT*, double *);
    
    MKL_INT indxg2p_(MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*);
    MKL_INT indxg2l_(MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*);
    //int indxg2p_(MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*);
    //int indxg2p_(int const& glob, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
    int indxl2g_(int const& loc, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
    //int indxg2l_(int const& glob, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
    
    //int infog2l_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
}

#endif // REML_H
