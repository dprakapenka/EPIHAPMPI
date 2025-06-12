#include <iostream>
#include <mpi.h>
#include <cassert>
#include <fstream>
#include <vector>
#include <random>
#include "../include/functions.h"

// Use the Makefile-defined TESTDIR variable
#ifndef TESTDIR
#define TESTDIR "test/"
#endif

// Function to determine the best BLACS grid using MPI_Dims_create
void get_best_blacs_grid(int &nprow, int &npcol) {
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int dims[2] = {0, 0};  // MPI will determine the best layout
	MPI_Dims_create(num_procs, 2, dims);
	nprow = dims[0];
	npcol = dims[1];
}

// Function to generate a positive definite and symmetric matrix
void generate_test_matrix_file(const std::string& filename, MKL_INT rows, MKL_INT cols, std::vector<double>& original_matrix) {
	std::string file_path = std::string(TESTDIR) + filename;
	std::ofstream file(file_path, std::ios::binary);
	assert(file && "Error: Cannot create test file.");

	file.write(reinterpret_cast<const char*>(&rows), sizeof(MKL_INT));
	file.write(reinterpret_cast<const char*>(&cols), sizeof(MKL_INT));

	std::vector<double> A(rows * cols);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);

	for (MKL_INT i = 0; i < rows * cols; i++) {
		A[i] = dist(gen);
	}

	std::vector<double> M(rows * cols, 0.0);
	for (MKL_INT i = 0; i < rows; i++) {
		for (MKL_INT j = 0; j < cols; j++) {
			double sum = 0.0;
			for (MKL_INT k = 0; k < cols; k++) {
				sum += A[i * cols + k] * A[j * cols + k];
			}
			M[i * cols + j] = sum;
		}
	}

	original_matrix = M;
	file.write(reinterpret_cast<const char*>(M.data()), M.size() * sizeof(double));
	file.close();
}

// Function to validate matrix read from binary file
void validate_matrix(const std::string& filename, MKL_INT rows, MKL_INT cols, const std::vector<double>& original_matrix) {
	std::string file_path = std::string(TESTDIR) + filename;

	/*
	   int argc = 0;
	   char** argv = nullptr;
	   MPI_Init(&argc, &argv);
	   */

	int flag;
	MPI_Initialized(&flag);
	if (!flag) {
		int argc = 0;
		char** argv = nullptr;
		MPI_Init(&argc, &argv);
	}


	int myrank_mpi;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);

	int nprow, npcol;
	get_best_blacs_grid(nprow, npcol);

	int ictxt, myrow, mycol;
	blacs_get(&ictxt, 0, &ictxt);
	blacs_gridinit(&ictxt, "R", nprow, npcol);
	blacs_gridinfo(&ictxt, &nprow, &npcol, &myrow, &mycol);

	MKL_INT descG[9];
	std::vector<double> Ga(rows * cols, 0.0);
	read_binmat(file_path.c_str(), Ga.data(), descG, 2 * sizeof(MKL_INT), rows * cols, myrank_mpi, nprow * npcol);

	blacs_gridexit(&ictxt);
	MPI_Finalize();
	blacs_exit(0);
}
