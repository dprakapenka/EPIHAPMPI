#include "read.h"
#include "options.h"
#include "constants.h"
#include "functions.h"
#include "mkl.h"
#include "mkl_scalapack.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"
#include "mpi.h"
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>

// Compute and print the numerical rank of a symmetric distributed matrix A
double compute_distributed_rank(
        MKL_INT N,
        double* A,
        MKL_INT descA[],
        MKL_INT ictxt,
        const Options& options)
{
    // blacs indices
    const MKL_INT IA = 1, JA = 1;

    // Workspace query
    MKL_INT lwork = -1, info;
    double work_query = 0;
    std::vector<double> W(N);

    // query optimal workspace
    // use lower triangle
    const char JOBZ = 'N';  // compute eigenvalues only
    pdsyev(&JOBZ, &CHAR_LOWER_,
            &N,
            A, &IA, &JA, descA,
            W.data(),
            nullptr, &IA, &JA, descA,
            &work_query, &lwork,
            &info);

    lwork  = static_cast<MKL_INT>(work_query);
    std::vector<double> WORK(lwork);

    // compute eigenvalues
    pdsyev(&JOBZ, &CHAR_LOWER_,
            &N, A, &IA, &JA, descA,
            W.data(),
            nullptr, &IA, &JA, descA,
            WORK.data(), &lwork,
            &info);
    if (info != 0) MPI_Abort(MPI_COMM_WORLD, info);

    // Count eigenvalues above tolerance
    double eps = std::numeric_limits<double>::epsilon();
    double lambda_max = W.back();
    double tol = eps * static_cast<double>(N) * lambda_max;
    int count = 0;
    for (MKL_INT i = 0; i < N; ++i) {
        if (std::abs(W[i]) > tol) ++count;
    }

    return count;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int myrank_mpi, nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    Options options;
    parse_options(argc, argv, myrank_mpi, options);

    // Compute close to square BLACS 2D grid
    std::tie(options.nprow, options.npcol) = compute_blacs_grid(nprocs_mpi);
    if (myrank_mpi == MPI_ROOT_PROC_) {
        std::cout << "Suggested *square processor grid: " << options.nprow << "," << options.npcol << std::endl;
    }

    // Initialize BLACS context
    MKL_INT ictxt;
    blacs_get(&IZERO_, &IZERO_, &ictxt);
    // For simplicity, use a square grid
    blacs_gridinit(&ictxt, &CHAR_LAYOUT_, &options.nprow, &options.npcol);


    MKL_INT block_size = 128;

    // Build map of possible GRM extensions
    std::map<std::string,std::string> var_extensions;
    var_extensions["a"]  = ".g.A";
    var_extensions["d"]  = ".g.D";
    var_extensions["h"]  = ".g.AH";
    var_extensions["aa"] = ".g.AA";
    var_extensions["ad"] = ".g.AD";
    var_extensions["dd"] = ".g.DD";
    var_extensions["aa-inter"] = ".g.AA-inter";
    var_extensions["aa-intra"] = ".g.AA-intra";
    var_extensions["ad-inter"] = ".g.AD-inter";
    var_extensions["ad-intra"] = ".g.AD-intra";
    var_extensions["dd-inter"] = ".g.DD-inter";
    var_extensions["dd-intra"] = ".g.DD-intra";
    var_extensions["aaa"] = ".g.AAA";
    var_extensions["aad"] = ".g.AAD";
    var_extensions["add"] = ".g.ADD";
    var_extensions["ddd"] = ".g.DDD";
    var_extensions["1"]  = ".g.1";
    var_extensions["2"]  = ".g.2";
    // ... add other extensions and
    // TODO:  move it to constants or init in read.cpp

    // collect (extension, rank) pairs
    //std::vector<std::pair<std::string,int>> collected_ranks;
    // collect (extension, rank, ind)
    std::vector<std::tuple<std::string, int, MKL_INT>> collected;

    // Descriptor array
    MKL_INT descG[DESC_LEN_];

    // For each variance type, try loading the corresponding GRM
    for (auto const& kv : var_extensions) {
        const std::string& name = kv.first;
        const std::string  ext  = kv.second;
        std::string filename = options.load_name + ext;

        // Check for file existence
        std::ifstream f(filename, std::ios::binary);
        if (!f.good()) {
            if (myrank_mpi == 0)
                std::cout << "Skipping missing file: " << filename << std::endl;
            continue;
        }
        f.close();

        // Verify and set global size (and ensure square)
        MKL_INT mkl_num_ind=0;
        verify_and_set_mat_size(filename, myrank_mpi, mkl_num_ind);

        // Get BLACS grid info
        MKL_INT my_row, my_col;
        blacs_gridinfo(&ictxt, &options.nprow, &options.npcol, &my_row, &my_col);

        // Set block size to bloc_size, or smaller if matrix is smaller
        options.mb = (mkl_num_ind < block_size ? mkl_num_ind : block_size);
        options.nb = (mkl_num_ind < block_size ? mkl_num_ind : block_size);

        // Initialize descriptor
        MKL_INT info;
        MKL_INT lld = std::max<MKL_INT>(1, numroc_(&mkl_num_ind, &options.mb, &my_row, &IZERO_, &options.nprow));
        descinit_(descG,
                &mkl_num_ind, &mkl_num_ind,
                &options.mb, &options.nb,
                &IZERO_, &IZERO_,
                &ictxt, &lld,
                &info);
        if (info != 0) {
            if (myrank_mpi == 0)
                std::cerr << "Error in descinit for " << name << " (info=" << info << ")" << std::endl;
            continue;
        }

        // Allocate local matrix buffer (using calloc to zero-initialize)
        MKL_INT local_rows = numroc_(&mkl_num_ind, &options.mb, &my_row, &IZERO_, &options.nprow);
        MKL_INT local_cols = numroc_(&mkl_num_ind, &options.nb, &my_col, &IZERO_, &options.npcol);
        MKL_INT local_size = std::max<MKL_INT>(1, local_rows) * std::max<MKL_INT>(1, local_cols);
        double* Gptr = static_cast<double*>(calloc(local_size, sizeof(double)));

        // Read the binary matrix into the distributed buffer
        read_binmat(filename.c_str(), Gptr, descG, local_size);

        // Compute and print rank
        int rank = compute_distributed_rank(
                mkl_num_ind,
                Gptr,
                descG,
                ictxt,
                options);
        if (myrank_mpi == 0)
            std::cout << "Rank of matrix '" << name << "' is: " << rank << std::endl;
        //collected_ranks.emplace_back(name, rank);
        collected.emplace_back(name, rank, mkl_num_ind);

        // Free the buffer for reuse
        free(Gptr);
    }

    if (myrank_mpi == 0) {
        // nothing left
        std::cout << std::endl << options.load_name << std::endl;
        std::cout << "GRM, RANK, IND" << std::endl;
         for (auto const& [varname, r, ind] : collected) {
            std::cout << varname << ", "
                << r << ", "
                << ind << std::endl;
        }
    }

    // Finalize BLACS and MPI
    blacs_gridexit(&ictxt);
    MPI_Finalize();
    return 0;
}
