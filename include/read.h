#ifndef READ_H
#define READ_H

#include "options.h"

#include <map>
#include <string>
#include <vector>
#include <istream>

#include "mkl.h"
#include "mpi.h"

/**
 * @brief Reads individual IDs from a file, preserving order and creating an index map.
 * @param ictxt is the blacs context
 * @param id_filename Path to the ID file (one ID per line).
 * @param ordered_ids Output parameter: A vector storing IDs in file order.
 * @return 0 on success, non-zero on error.
 */
int read_individual_ids(
        MKL_INT ictxt,
        const std::string& id_filename,
        std::vector<std::string>& ordered_ids
        );


// Structure to map phenotype record index to individual column index
struct PhenoRecord {
    MKL_INT record_index; // 0-based index of the valid record
    MKL_INT ind_col_index; // 0-based column index in Z/G corresponding to the individual
};

/**
 * @brief Reads the phenotype file, handling missing values and mapping records to individuals.
 * @param options Program options containing file paths, missing value indicator, etc.
 * @param id_to_index A map from individual ID (string) to its 0-based column index (must be pre-populated based on genotype order).
 * @param pheno_records Output parameter: A vector storing the mapping of valid phenotype records to individual indices.
 * @param mkl_num_rec Output parameter: the total number of valid phenotype records found.
 * @param Y_values Output parameter: A vector containing the phenotype values for valid records.
 * @return 0 on success, non-zero on error.
 */
int read_phenotype_data(
        const Options &options,
        MKL_INT ictxt,
        const std::vector<std::string>& individual_ids,
        std::vector<int>   &rec2ind,                         // output record ids to geno ind ids
        int &num_rec,                                         // Output count of valid records
        std::vector<double>& phenotypes                      // Output phenotype values
        );


/**
 * Load all G-matrices specified in options.variances (excluding 'e'), allocate
 * local block-cyclic storage, read each via read_binmat, and record pointers
 * in Glocal. Also populate var_idx, idx_var, h2, and update mem_tot and mkl_num_ind.
 * This function queries MPI rank and BLACS grid internally.
 */
void load_grms(
        const Options &options,
        MKL_INT ictxt,
        MKL_INT &mkl_num_ind, 
        std::map<std::string,int> &var_idx,
        std::map<int,std::string> &idx_var,
        std::map<std::string,double*> &Glocal,
        MKL_INT descG[],
        std::map<std::string,double> &h2,
        MKL_INT mem_tot 
        );


int fillBuff(
        std::istream & infile,
        std::vector<char> &buff
        );

int countLines(
        const std::vector<char> &buff,
        int bs
        );

void read_file_size(
        const char* filename,
        int &r,
        int &c
        );

void read_binmat_size(
        const char* filename,
        MKL_INT& r,
        MKL_INT& c
        );

void read_binmat(
        const char *filename,
        double *mat,
        MKL_INT *descmat,
        MKL_INT nelem);

/** Read matrix size from file, ensure it is non-zero, square,
 *  and consistent across all files; on first call sets
 *  number of individuals.  Aborts if any check fails.
 */
void verify_and_set_mat_size(
        const std::string& load_file,
        int                myrank_mpi,
        MKL_INT&           mkl_num_ind
        );

std::vector<std::string> str_split(
        const std::string &delimiter,
        const std::string &str
        );

std::vector<int> str_split_int(
        const std::string &delimiter,
        const std::string &str
        );

#endif // READ_H
