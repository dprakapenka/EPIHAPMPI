#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <map>
#include "mkl.h"

// Struct definition for Options
struct Options {
    int iter_n = 1;
    int num_threads = 1;
    std::string pheno_file = "phenotype.txt";
    std::string save_name = "out";
    std::string load_name = "grm";
    int trait_col = 2;
    MKL_INT mb = 128;
    MKL_INT nb = 128;
    int ai_start = 3;
    double tolerance = 1.0e-8;
    double htolerance = 1.0e-6;
    //remove, covered by variances
    double var_snp_a = 0.0;
    double var_snp_d = 0.0;
    double var_hap_a = 0.0;
    double var_snp_aa = 0.0;
    double var_snp_ad = 0.0;
    double var_snp_dd = 0.0;
    double var_snp_aaa = 0.0;
    double var_e = 0.0;
    double missing_phenotype = -9999111.0;

    std::map<std::string,double> variances;
    MKL_INT nprow;
    MKL_INT npcol;
};

// Function declarations
void parse_options(int argc, char *argv[], int myrank_mpi, Options &options);

#endif
