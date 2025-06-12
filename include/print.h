#ifndef PRINT_H
#define PRINT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include "options.h" // Include necessary headers
#include "mkl.h"     // For MKL_INT

// Print variance update line:
//   v<name>:    <var>    diff: <diff>
void print_variances(const std::string& name, double var, double diff);

// Print variance update line:
//   v<name>:    <var>
void print_variances(const std::string& name, double var);

// Print heritability (h2) update line:
//   h2_<name>:  <h2_value>   diff: <diff>
//   For residual ("e"), label "H2:" is used.
void print_heritabilities(const std::string& name, double h2_value, double diff);

// Print heritability (h2) update line:
//   h2_<name>:  <h2_value>
//   For residual ("e"), label "H2:" is used.
void print_heritabilities(const std::string& name, double h2_value);

// Print tolerance reached exit for variance or heritability
void print_update_header(const std::string& tol_type, int myrank);

// Print iteration and time at the end of iteration
void print_end_iteration(int iteration, double iteration_time);

// Print when reached either heirtability or variance tolerance
void print_reached_tolerance(const std::string& tolerance_type, double tolerance, int iteration);

// Print out options that are set
void print_options(const Options &options);

/**
 * @brief Saves the final REML variance components and heritabilities to a CSV file.
 *
 * Only the root MPI process performs the file writing.
 *
 * @param myrank_mpi            Rank of the current MPI process.
 * @param options               The program options, containing output_name.
 * @param final_variances       Map containing the final variance estimates (component name -> variance).
 * @param final_heritabilities  Map containing the final heritability estimates (component name -> h2).
 */
void save_reml_results(
    const Options& options,
    const std::map<std::string, double>& final_heritabilities
);

/**
 * @brief Saves the calculated total GEBVs to a CSV file.
 *
 * Only the root MPI process performs the file writing. It assumes the
 * total_gebv vector already exists fully formed on the root process.
 *
 * @param myrank_mpi            Rank of the current MPI process.
 * @param options               The program options, containing output_name.
 * @param individual_ids        Vector containing the IDs of the individuals corresponding to the GEBVs.
 * @param total_gebv            Vector containing the calculated total GEBV for each individual.
 * @param expected_individuals  The expected number of individuals (size of the vectors).
 */
void save_gblup_results(
    const Options& options,
    const std::vector<std::string>& individual_ids,
    const std::vector<double>& total_gebv,
    MKL_INT expected_individuals
);


/**
 * @brief Saves the calculated GBLUP effects for each component to a CSV file.
 *
 * Only the root MPI process performs the file writing. Assumes component_gebv
 * map contains fully gathered vectors on the root process.
 * Output format: ID <tab> Component1_BLUP <tab> Component2_BLUP ...
 *
 * @param options             The program options, containing output_name.
 * @param individual_ids      Vector containing the IDs of the individuals.
 * @param component_gebv      Map where key is component name (e.g., "A", "D")
 * and value is the gathered vector of GEBVs for that component.
 * @param expected_individuals The expected number of individuals (size of vectors).
 * @param component_order     Vector specifying the desired order of component columns in the output file.
 */
void save_multi_component_gblup(
    const Options& options, // Contains options.variances and save_name
    const std::vector<std::string>& individual_ids,
    const std::map<std::string, std::vector<double>>& component_gebv,
    MKL_INT expected_individuals,
    // const std::vector<std::string>& component_order, // Parameter removed
    const std::set<std::string>& validation_ids
);

void print_gblup(
    const Options& options,
    const std::vector<std::string>& individual_ids,
    const std::vector<int>& rec2ind,
    const std::map<std::string, std::vector<double>>& u_effects_global
);

#endif // PRINT_H
