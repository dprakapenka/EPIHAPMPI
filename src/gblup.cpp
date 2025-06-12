#include "gblup.h"
#include "options.h"
#include "constants.h"
#include "print.h"
#include <algorithm>
#include "mkl.h"
#include "mkl_scalapack.h"
#include "mkl_trans.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"
#include "mpi.h"

void calculate_gblup(
    const Options &options,
    MKL_INT        ictxt,
    MKL_INT        mkl_num_ind,  // # individuals = Z.rows()
    // Data buffers from main():
    double        *Z,             const MKL_INT *descZ,
    double        *PY,            const MKL_INT *descY,
    std::map<std::string,double*> &Glocal, const MKL_INT *descG,
    // Outputs:
    std::map<std::string, std::vector<double>> &u_effects_global
) {
    // BLACS grid info
    MKL_INT iam, nprocs;
    MKL_INT nprow, npcol, myrow, mycol;
    blacs_pinfo(&iam, &nprocs);
    blacs_gridinfo(&ictxt, &nprow, &npcol, &myrow, &mycol);

    if (iam == MPI_ROOT_PROC_) {
        fprintf(stderr, "\n========= GBLUP begins =========\n\n");
        // print final variances
        for (auto &kv : options.variances) {
            if (kv.first == "e") break;
            printf("var_%s = %.6e\n", kv.first.c_str(), kv.second);
        }
    }

    // Z' * PY  distributed ZPY (mkl_num_ind x 1)
    double *ZPY;
    MKL_INT descZPY[DESC_LEN_], info;
    // allocate local buffer for ZPY
    // compute local block size
    MKL_INT nr = std::max< MKL_INT >(1, numroc_(&mkl_num_ind, &options.mb, &myrow, &IZERO_, &nprow));
    MKL_INT nc = std::max< MKL_INT >(1, numroc_(&mkl_num_ind, &options.nb, &mycol, &IZERO_, &npcol));
    ZPY = (double*)calloc(nr*nc, sizeof(double));
    descinit_(descZPY, &mkl_num_ind, &IONE_, &options.mb, &IONE_, &IZERO_, &IZERO_, &ictxt, &nr, &info);

    std::cout << "allocated zpy" << std::endl;

    pdgemm_(&CHAR_TRANS_,&CHAR_NOTRANS_,
            &mkl_num_ind, &IONE_, &mkl_num_ind,
            &DONE_,
            Z,   &IONE_, &IONE_, descZ,
            PY,  &IONE_, &IONE_, descY,
            &DZERO_,
            ZPY, &IONE_, &IONE_, descZPY);

    std::cout << "pdgemm zpy" << std::endl;

    // For each component k in options.variances (skip residual):
    u_effects_global.clear();
    for (auto &kv : options.variances) {
        const std::string &name = kv.first;
        double var_k = kv.second;
        if (name == "e") break;

        printf("var_%s = %.6e\n", kv.first.c_str(), kv.second);

        // u_local = var_k * G_k * ZPY to local u of length mkl_num_ind
        double *u_local;
        MKL_INT descU[DESC_LEN_];
        u_local = (double*)calloc(nr*nc, sizeof(double));
        descinit_(descU, &mkl_num_ind, &IONE_, &options.mb, &IONE_, &IZERO_, &IZERO_, &ictxt, &nr, &info);
        pdgemm_(&CHAR_NOTRANS_, &CHAR_NOTRANS_,
                &mkl_num_ind, &IONE_, &mkl_num_ind,
                &var_k,
                Glocal[name], &IONE_, &IONE_, descG,
                ZPY,          &IONE_, &IONE_, descZPY,
                &DZERO_,
                u_local,      &IONE_, &IONE_, descU);
        std::cout << "pdgemm u_local" << std::endl;

        // Gather u_local to u_global via pdgemr2d into contiguous on root
        std::vector<double> u_global;
        {
            // rep buffer on root
            double *Urep = nullptr;
            MKL_INT descR[DESC_LEN_];
            if (iam == MPI_ROOT_PROC_) {
                Urep = (double*)malloc(sizeof(double)*mkl_num_ind);
                descinit_(descR, &mkl_num_ind, &IONE_, &mkl_num_ind, &IONE_, &IZERO_, &IZERO_, &ictxt, &mkl_num_ind, &info);
            } else {
                descinit_(descR, &mkl_num_ind, &IONE_, &mkl_num_ind, &IONE_, &IZERO_, &IZERO_, &ictxt, &IONE_, &info);
            }
            pdgemr2d_(&mkl_num_ind, &IONE_,
                     u_local, &IONE_, &IONE_, descU,
                     Urep,     &IONE_, &IONE_, descR,
                     &ictxt);
            if (iam == MPI_ROOT_PROC_) {
                u_global.assign(Urep, Urep + mkl_num_ind);
                free(Urep);
            }
        }
        u_effects_global[name] = std::move(u_global);

        free(u_local);
    }

    // cleanup
    free(ZPY);
}
