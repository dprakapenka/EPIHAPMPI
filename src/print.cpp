#include "print.h"
#include "options.h"
#include "constants.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <numeric> // For std::accumulate to calculating total h2
#include <set> 


void print_variances(const std::string& name, double variance, double diff) {
	std::string tmpStrName = "v" + name + ":";
	std::cout
		<< std::setw(20) << tmpStrName
		<< std::scientific << std::setprecision(12) << variance
		<< std::setw(20) << " diff:"
		<< std::scientific << std::setprecision(12) << diff
		<< std::endl;
}


void print_variances(const std::string& name, double variance) {
	std::string tmpStrName = "v" + name + ":";
	std::cout
		<< std::setw(20) << tmpStrName
		<< std::scientific << std::setprecision(12) << variance
		<< std::endl;
}


void print_heritabilities(const std::string& name, double h2_value, double diff) {
	std::string tmpStrName;
	if (name == "e") {
		tmpStrName = "H2:";
	} else {
		tmpStrName = "h2_" + name + ":";
	}
	std::cout
		<< std::setw(20) << tmpStrName
		<< std::scientific << std::setprecision(12) << h2_value
		<< std::setw(20) << " diff:"
		<< std::scientific << std::setprecision(12) << diff
		<< std::endl;
}


void print_heritabilities(const std::string& name, double h2_value) {
	std::string tmpStrName;
	if (name == "e") {
		tmpStrName = "H2:";
	} else {
		tmpStrName = "h2_" + name + ":";
	}
	std::cout
		<< std::setw(20) << tmpStrName
		<< std::scientific << std::setprecision(12) << h2_value
		<< std::endl;
}


void print_end_iteration(int iteration, double iteration_time) {
	std::cout
		<< std::endl
		<< "iteration " << iteration
		<< " completed in "
		<< std::fixed
		<< std::setprecision(2)
		<< iteration_time
		<< " seconds" << std::endl;

}


void print_reached_tolerance(const std::string& tolerance_type, double tolerance, int iteration) {
	std::cout
		<< std::endl
		<< "reached "
		<< tolerance_type
		<< " tolerance of "
		<< std::setprecision(2)
		<< std::scientific
		<< tolerance
		<< " at iteration " << iteration << std::endl;
}


void print_options(const Options &options) {
	const std::string sep = "--------------------------------------------";

	std::cout << "\n" << sep << "\n";
	std::cout << "                  OPTIONS\n";
	std::cout << sep << "\n";

	std::cout << std::left << std::setw(28) << "Num threads:"              << options.num_threads << "\n";
	std::cout << std::left << std::setw(28) << "Max iterations:"           << options.iter_n << "\n";
	std::cout << std::left << std::setw(28) << "Variance tolerance:"       << options.tolerance << "\n";
	std::cout << std::left << std::setw(28) << "Heritability Tolerance:"   << options.htolerance << "\n";
	std::cout << std::left << std::setw(28) << "AI starts at iteration:"   << options.ai_start << "\n";

	std::cout << "\nStarting variances:\n";
	for (const auto& [var, val] : options.variances) {
		std::cout << "  " << std::setw(3) << var << " = " << val << std::endl;
	}

	std::cout << "\nPaths and Files:\n";
	std::cout << std::left << std::setw(28) << "project name (GRMs):"      << options.load_name << "\n";
	std::cout << std::left << std::setw(28) << "Output prefix:"            << options.save_name << "\n";
	std::cout << std::left << std::setw(28) << "Phenotype file:"           << options.pheno_file << "\n";
	std::cout << std::left << std::setw(28) << "Phenotype/Trait column:"   << options.trait_col << "\n";
	std::cout << std::left << std::setw(28) << "Missing value identifier:" << options.missing_phenotype << "\n";

	std::cout << sep << "\n\n";
}

void save_reml_results(
		const Options& options,
		const std::map<std::string, double>& final_heritabilities)
{

	if (options.save_name.empty()) {
		std::cerr << "WARNING: Output name prefix is empty. Cannot save REML results file." << std::endl;
		return;
	}

	std::string filename = options.save_name + "_greml.csv";
	std::ofstream outfile(filename);

	if (!outfile.is_open()) {
		std::cerr << "ERROR: Could not open file " << filename << " for writing REML results." << std::endl;
		return;
	}

	std::cout << "Writing final REML results to: " << filename << std::endl;

	// Write Header
	outfile << "Component,Variance,Heritability" << std::endl;

	// Write component data
	// Iterate through variances map, assuming heritabilities map has corresponding keys
	double total_h2 = 0.0;
	for (const auto& var_pair : options.variances) {
		const std::string& name = var_pair.first;
		double variance = var_pair.second;
		double heritability = 0.0;

		auto h2_it = final_heritabilities.find(name);
		if (h2_it != final_heritabilities.end()) {
			heritability = h2_it->second;
			if (name != "e") { // Accumulate non-residual h2 for total
				total_h2 += heritability;
			}
		} else if (name != "e") {
			std::cerr << "WARNING: Heritability not found for component '" << name << "' in save_reml_results." << std::endl;
		}

		outfile << name << ","
			<< std::scientific << std::setprecision(10) << variance << ","
			<< std::fixed << std::setprecision(10) << heritability << std::endl;
	}

	outfile << "Total_h2," << std::fixed << std::setprecision(10) << total_h2 << ",NA" << std::endl;

	outfile.close();
}


void save_gblup_results(
		const Options& options,
		const std::vector<std::string>& individual_ids,
		const std::vector<double>& total_gebv,
		MKL_INT expected_individuals)
{

	if (options.save_name.empty()) {
		std::cerr << "WARNING: Output name prefix is empty. Cannot save GBLUP results file." << std::endl;
		return;
	}

	// Validate input sizes
	if (total_gebv.empty()) {
		std::cerr << "WARNING: Total GEBV vector is empty. GBLUP output file not written." << std::endl;
		return;
	}
	if (static_cast<MKL_INT>(total_gebv.size()) != expected_individuals) {
		std::cerr << "WARNING: Total GEBV vector size (" << total_gebv.size()
			<< ") does not match expected individuals (" << expected_individuals
			<< "). GBLUP output file not written." << std::endl;
		return;
	}
	if (individual_ids.size() != static_cast<size_t>(expected_individuals)) {
		std::cerr << "WARNING: Individual ID vector size (" << individual_ids.size()
			<< ") does not match expected individuals (" << expected_individuals
			<< "). GBLUP output file may have incorrect IDs." << std::endl;
		// Decide whether to proceed or return. Proceeding for now.
	}


	std::string filename = options.save_name + "_gblup.csv";
	std::ofstream outfile(filename);

	if (!outfile.is_open()) {
		std::cerr << "ERROR: Could not open file " << filename << " for writing GBLUP results." << std::endl;
		return;
	}

	std::cout << "Writing total GEBV to: " << filename << std::endl;

	// Write Header
	outfile << "ID,TotalGEBV" << std::endl;

	// Write data
	for (MKL_INT i = 0; i < expected_individuals; ++i) {
		std::string id = (i < static_cast<MKL_INT>(individual_ids.size())) ? individual_ids[i] : "UnknownID_" + std::to_string(i);
		outfile << id << ","
			<< std::scientific << std::setprecision(10) << total_gebv[i] << std::endl;
	}

	outfile.close();
}


void save_multi_component_gblup(
    const Options& options, // Contains options.variances and save_name
    const std::vector<std::string>& individual_ids,
    const std::map<std::string, std::vector<double>>& component_gebv,
    MKL_INT expected_individuals,
    // const std::vector<std::string>& component_order, // Parameter removed
    const std::set<std::string>& validation_ids
)
{
    // Assumes function is only called by the root process

    // Define the separator character
    const char sep = ',';

    // --- Initial checks ---
    if (options.save_name.empty()) {
        std::cerr << "WARNING: Output name prefix is empty. Cannot save multi-component GBLUP results file." << std::endl;
        return;
    }
    // Check if component_gebv map is empty - might be valid if only residual variance exists
    // Check if options.variances only contains 'e' - also valid
    bool only_residual = true;
    for (const auto& pair : options.variances) {
        if (pair.first != "e") {
            only_residual = false;
            break;
        }
    }
     if (component_gebv.empty() && !only_residual) {
         std::cerr << "WARNING: GBLUP component results map is empty, but non-residual variances exist. Output file may be incomplete." << std::endl;
         // Continue for now, will likely print NAs
     }
    if (individual_ids.size() != static_cast<size_t>(expected_individuals)) {
        std::cerr << "WARNING: Individual ID vector size (" << individual_ids.size()
                  << ") does not match expected individuals (" << expected_individuals
                  << "). GBLUP output file may have incorrect row count or IDs." << std::endl;
    }
    // --- End Initial checks ---


    // --- Calculate Total GEBV vector ---
    std::vector<double> total_gebv(expected_individuals, 0.0); // Initialize with zeros
    bool total_calculation_ok = true; // Flag to track if total is valid

    // Iterate directly over options.variances to determine order and sum
    for (const auto& pair : options.variances) {
        const std::string& comp_name = pair.first;
        if (comp_name == "e") continue; // Skip residual component

        auto it = component_gebv.find(comp_name);
        if (it != component_gebv.end()) {
            // Found the component vector in the results
            const std::vector<double>& u_i_vector = it->second;
            if (static_cast<MKL_INT>(u_i_vector.size()) == expected_individuals) {
                // Add this component's values to the total
                for (MKL_INT i = 0; i < expected_individuals; ++i) {
                    total_gebv[i] += u_i_vector[i];
                }
            } else {
                // Size mismatch, total cannot be calculated correctly
                total_calculation_ok = false;
                std::cerr << "WARNING: Cannot add component '" << comp_name
                          << "' to GBLUP_G due to size mismatch (" << u_i_vector.size()
                          << " vs " << expected_individuals << ")." << std::endl;
            }
        } else {
            // Component from options.variances not found in results map
            total_calculation_ok = false;
            std::cerr << "WARNING: Cannot add component '" << comp_name
                      << "' to GBLUP_G as it was not found in results map." << std::endl;
        }
    } // End loop over options.variances for summing
    // --- END: Calculate Total GEBV vector ---


    // --- Open File ---
    std::string filename = options.save_name + "_gblup.csv";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "ERROR: Could not open file " << filename << " for writing multi-component GBLUP results." << std::endl;
        return;
    }

    std::cout << "Writing multi-component GBLUP results to: " << filename << std::endl;

    // --- Write Header ---
    outfile << "ID";
    // Iterate directly over options.variances for header order
    for (const auto& pair : options.variances) {
        const std::string& comp_name = pair.first;
        if (comp_name == "e") continue; // Skip residual component
         outfile << sep << "GBLUP_" << comp_name;
    }
    outfile << sep << "GBLUP_G"; // Use "GBLUP_G" for the total column
    outfile << sep << "Train./Valid."; // Add Train/Valid header
    outfile << std::endl; // End header line

    // --- Write Data Rows ---
    for (MKL_INT i = 0; i < expected_individuals; ++i) {
        // Get ID, handling potential size mismatch
        std::string id = (i < static_cast<MKL_INT>(individual_ids.size())) ? individual_ids[i] : "UnknownID_" + std::to_string(i);
        outfile << id;

        // Write individual component values iterating directly over options.variances
        for (const auto& pair : options.variances) {
            const std::string& comp_name = pair.first;
            if (comp_name == "e") continue; // Skip residual component

            outfile << sep; // Separator
            auto it = component_gebv.find(comp_name);
            if (it != component_gebv.end()) {
                // Found the component vector
                const std::vector<double>& u_i_vector = it->second;
                // Check vector size and index validity
                if (i < static_cast<MKL_INT>(u_i_vector.size())) {
                    outfile << std::scientific << std::setprecision(10) << u_i_vector[i];
                } else {
                    outfile << "SizeMismatch"; // Indicate error
                }
            } else {
                // Component from options.variances not found in results map
                outfile << "NA";
            }
        } // end loop through components (options.variances)

        // Write the total GEBV value ("GBLUP_G")
        outfile << sep; // Separator for the total column
        if (total_calculation_ok && i < static_cast<MKL_INT>(total_gebv.size())) { // Check flag and index
             outfile << std::scientific << std::setprecision(10) << total_gebv[i];
        } else {
             outfile << "NA"; // Indicate total could not be calculated reliably or index issue
        }

        // Write Train/Valid status
		// TODO: change to validation IDs - smaller
        outfile << sep; // Separator for the Train/Valid column
        if (validation_ids.count(id)) { // Check if the ID is in the validation set
            outfile << "V";
        } else {
            outfile << "T";
        }

        outfile << std::endl; // End data row
    } // end loop through individuals

    outfile.close();
} // End of save_multi_component_gblup function

void print_gblup(
    const Options& options,
    const std::vector<std::string>& individual_ids,
    const std::vector<int>& rec2ind,
    const std::map<std::string, std::vector<double>>& u_effects_global
) {
    // Build a training mask: T if individual had a phenotype (in rec2ind)
    std::vector<bool> is_training(individual_ids.size(), false);
    for (int r : rec2ind) {
        if (r >= 0 && r < (int)is_training.size())
            is_training[r] = true;
    }

    std::ofstream fout(options.save_name + "_gblup.csv");
    // Header
    fout << "ID";
    for (auto &kv : options.variances) {
        const std::string &name = kv.first;
        if (name == "e") break;
        fout << ",GBLUP_" << name;
    }
    fout << ",Train./Valid.\n";

    // Rows
    int n = individual_ids.size();
    for (int i = 0; i < n; ++i) {
        fout << individual_ids[i];
        for (auto &kv : options.variances) {
            const std::string &name = kv.first;
            if (name == "e") break;
            fout << ",";
            auto it = u_effects_global.find(name);
            if (it != u_effects_global.end() && i < (int)it->second.size()) {
                fout << std::setprecision(7) << it->second[i];
            } else {
                fout << "NA";
            }
        }
        fout << "," << (is_training[i] ? "T" : "V") << "\n";
    }
    fout.close();
}
