#!/bin/env python
import pandas as pd
import numpy as np
import argparse
import os
import multiprocessing
from multiprocessing import Pool

# Try importing numpy
try:
    import numpy as np
except ImportError as e:
    print("\nNumPy are not installed. Please install it.")

# Try importing matplotlib
try:
    import matplotlib.pyplot as plt
    import numpy as np
    PLOTTING_AVAILABLE = True
except ImportError as e:
    PLOTTING_AVAILABLE = False

import matplotlib.pyplot as plt
import itertools

# Toggle this to True for muted colors, False for bright vibrant colors
USE_MUTED_COLORS = True

# Define a consistent color mapping for GBLUP columns
# vibrant color scheme
GBLUP_COLOR_MAP_BRIGHT = {
    "GBLUP_A": "blue",
    "GBLUP_D": "brown",
    "GBLUP_AA": "green",
    "GBLUP_AA-intra": "purple",
    "GBLUP_AA-inter": "orange",
    "GBLUP_AD": "cyan",
    "GBLUP_DD": "pink",
    "GBLUP_AH": "gray",
    "GBLUP_G": "red",
}

# Muted color scheme for softer visualization
GBLUP_COLOR_MAP_MUTED = {
    "GBLUP_A": "#5e81ac",  # Muted Blue
    "GBLUP_D": "#d08770",  # Warm Muted Orange
    "GBLUP_AA": "#88c0d0",  # Soft Cyan
    "GBLUP_AA-intra": "#bf616a",  # Muted Red
    "GBLUP_AA-inter": "#a3be8c",  # Soft Olive Green
    "GBLUP_AD": "#ebcb8b",  # Warm Beige/Gold
    "GBLUP_DD": "#b48ead",  # Muted Purple
    "GBLUP_AH": "#7d8491",  # Cool Grayish Blue
    "GBLUP_G": "#4c566a",  # Dark Slate Gray
}


# Assign the selected color scheme
GBLUP_COLOR_MAP = GBLUP_COLOR_MAP_MUTED if USE_MUTED_COLORS else GBLUP_COLOR_MAP_BRIGHT

def check_file_validity(file_path):
    """Check if a file exists and is not empty."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False
    return True


def read_file_with_auto_delimiter(file_path):
    """Read a file and automatically detect its delimiter (comma, tab, or space)."""
    with open(file_path, 'r') as f:
        first_line = f.readline()
        if ',' in first_line:
            return pd.read_csv(file_path)
        else:
            return pd.read_csv(file_path, sep=r'\s+')


def calculate_file_correlations(gblup_path, pheno, phenotype_col, missing_value):
    """Calculate correlations for a single GBLUP file."""
    # Load GBLUP file
    gblup = read_file_with_auto_delimiter(gblup_path)

    pheno = pheno.rename(columns={'id': 'ID'})
    gblup = gblup.rename(columns={'id': 'ID'})

    # Filter out rows with missing phenotypes
    pheno_filtered = pheno[pheno[phenotype_col] != missing_value]

    # Merge the data using the identified ID column
    merged_df = pd.merge(gblup, pheno_filtered, on='ID')

    # Split into training and validation sets
    train_df = merged_df[merged_df['Train./Valid.'] == 'T']
    valid_df = merged_df[merged_df['Train./Valid.'] == 'V']

    # Extract GBLUP columns
    gblup_cols = [col for col in gblup.columns if col.startswith('GBLUP_')]

    # Compute correlations safely
    def safe_corr(x, y):
        if len(x) > 1 and len(y) > 1:
            return np.corrcoef(x, y)[0, 1]
        return np.nan  # Return NaN if correlation is undefined

    train_corr = {col: safe_corr(train_df[col], train_df[phenotype_col]) for col in gblup_cols}
    valid_corr = {col: safe_corr(valid_df[col], valid_df[phenotype_col]) for col in gblup_cols}

    return gblup_path, train_corr, valid_corr


def process_gblup_files_parallel(gblup_files, pheno, phenotype_col, missing_value, n_cores):
    """Process multiple GBLUP files in parallel."""
    args = [(file, pheno, phenotype_col, missing_value) for file in gblup_files]
    # Limit number of processes to the number of GBLUP files
    n_cores = min(n_cores, len(gblup_files))
    with Pool(processes=n_cores) as pool:
        results = pool.starmap(calculate_file_correlations, args)
    return results


def print_results(results):
    """Print individual file correlations and calculate averages."""
    if not results:
        print("\nNo valid GBLUP files to process. Skipping correlation calculations.")
        return
    all_train_corrs = {}
    all_valid_corrs = {}
    print("\nCorrelations for Individual Files:")

    for gblup_path, train_corr, valid_corr in results:
        print(f"\nFile: {gblup_path}")
        print("Training Set:")
        for col, corr in train_corr.items():
            print(f"  {col}: {corr:.4f}" if not np.isnan(corr) else f"  {col}: NaN")
            if not np.isnan(corr):
                all_train_corrs.setdefault(col, []).append(corr)
        print("Validation Set:")
        for col, corr in valid_corr.items():
            print(f"  {col}: {corr:.4f}" if not np.isnan(corr) else f"  {col}: NaN")
            if not np.isnan(corr):
                all_valid_corrs.setdefault(col, []).append(corr)

    # Calculate average correlations while excluding NaNs
    num_valid_files = len(results)
    if num_valid_files > 1:
        print(f"\nAverage correlations across {num_valid_files} runs:")
        print("  Column       ,   Mean   ,  Std Dev , Set ")
        print("  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,")
        for col, corrs in all_train_corrs.items():
            mean_corr = np.nanmean(corrs)
            std_corr = np.nanstd(corrs, ddof=1) if len(corrs) > 1 else np.nan  # Compute standard deviation
            #print(f"  {col}: mean = {mean_corr:.4f}, std = {std_corr:.4f}")
            print(f"  {col:<12} , {mean_corr:>8.4f} , {std_corr:>8.4f} ,  T  ")
        print("  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,")
        for col, corrs in all_valid_corrs.items():
            mean_corr = np.nanmean(corrs)
            std_corr = np.nanstd(corrs, ddof=1) if len(corrs) > 1 else np.nan  # Compute standard deviation
            #print(f"  {col}: mean = {mean_corr:.4f}, std = {std_corr:.4f}")
            print(f"  {col:<12} , {mean_corr:>8.4f} , {std_corr:>8.4f} ,  V  ")


def save_results_to_csv(results, csv_filename):
    """Save correlation data to a CSV file."""
    data = []
    
    # Collect all individual correlations
    for gblup_path, train_corr, valid_corr in results:
        for col, corr in train_corr.items():
            data.append([os.path.basename(gblup_path), col, corr, "T"])
        for col, corr in valid_corr.items():
            data.append([os.path.basename(gblup_path), col, corr, "V"])
    
    # Convert to DataFrame
    df = pd.DataFrame(data, columns=["File", "Column", "Correlation", "Set"])
    
    # Compute and append averages and standard deviations
    grouped = df.groupby(["Column", "Set"])["Correlation"].agg(["mean", "std"]).reset_index()
    grouped.insert(0, "File", "Average")  # Mark averages separately
    
    # Combine individual and summary results
    final_df = pd.concat([df, grouped], ignore_index=True)
    
    # Save to CSV
    final_df.to_csv(csv_filename, index=False)
    print(f"Correlation results saved to {csv_filename}")


def plot_correlation_graph(results, gblup_cols, save_plot=None, plot_set="V"):
    """Plot a graph of correlations for the chosen GBLUP columns."""
    file_numbers = []
    correlations_dict = {col: [] for col in gblup_cols}

    for idx, (gblup_path, train_corr, valid_corr) in enumerate(results, 1):
        for col in gblup_cols:
            # Choose which set to plot based on the `plot_set` argument
            corr = train_corr.get(col) if plot_set == "T" else valid_corr.get(col)
            if corr is not None and not np.isnan(corr):
                correlations_dict[col].append(corr)
            else:
                correlations_dict[col].append(np.nan)

        file_numbers.append(idx)

    # Plotting the correlations
    plt.figure(figsize=(10, 6))

    # Update plotting to use the defined colors
    for col in gblup_cols:
        if any(not np.isnan(c) for c in correlations_dict[col]):  # Only plot non-empty data
            plt.plot(file_numbers, correlations_dict[col], marker='o', linestyle='-',
                    label=col, color=GBLUP_COLOR_MAP.get(col, "black"))

    plt.title(f"Correlation of GBLUP Columns across samples ({'Training' if plot_set == 'T' else 'Validation'} Set)")
    plt.xlabel("File Number")
    plt.ylabel("Correlation with Phenotype")
    plt.xticks(file_numbers)
    plt.legend()
    plt.grid(True)

    # Save the plot if specified
    if save_plot:
        plt.savefig(save_plot)
        print(f"Plot saved to {save_plot}")
    else:
        plt.show()


def main():
    # Argument parser for command-line options
    parser = argparse.ArgumentParser(description="Calculate correlations between GBLUP values and phenotype.")
    parser.add_argument("--gblup", nargs="+", required=True, help="List of GBLUP file(s) (space-separated).")
    parser.add_argument("--pheno", required=True, help="Path to the phenotype file.")
    parser.add_argument("--phenotype_col", default="afc", help="Phenotype column name in the phenotype file.")
    parser.add_argument("--missing_value", type=float, default=-9999111, help="Value denoting missing phenotypes.")
    parser.add_argument("--n_cores", type=int, default=multiprocessing.cpu_count(),
                        help="Number of parallel processes to use (default: system's CPU count).")
    parser.add_argument("--gblup_cols", nargs="+", default=["all"],
                        help="The GBLUP columns to plot (space-separated). Use GBLUP_A, GBLUP_G, etc. or 'all' to plot all columns.")
    parser.add_argument("--save_plot", type=str, help="File to save the plot. If not specified, the plot will be shown.")
    parser.add_argument("--plot_set", choices=["T", "V"], default="V", help="Choose which set to plot: 'T' for training, 'V' for validation.")
    parser.add_argument("--save_csv", type=str, help="File to save the correlation results as a CSV.")
    args = parser.parse_args()

    # If 'all' is chosen, select all GBLUP columns from the global color mapping
    if "all" in args.gblup_cols:
        args.gblup_cols = list(GBLUP_COLOR_MAP.keys())


    # Check gblup files
    valid_gblup_files = [file for file in args.gblup if check_file_validity(file)]
    if not valid_gblup_files:
        print("No valid GBLUP files provided. Exiting.")
        return

    # Load phenotype
    pheno = read_file_with_auto_delimiter(args.pheno)

    # Process GBLUP files in parallel
    results = process_gblup_files_parallel(valid_gblup_files, pheno, args.phenotype_col, args.missing_value, args.n_cores)


    # Print results
    print_results(results)

    # If CSV output is requested, save results and exit
    if args.save_csv:
        save_results_to_csv(results, args.save_csv)

    # Plot the correlation graph for the chosen GBLUP columns, only if plotting libraries are available
    if PLOTTING_AVAILABLE:
        plot_correlation_graph(results, args.gblup_cols, args.save_plot, args.plot_set)
    else:
        print("\nMatplotlib library is not installed. Skipping plot functionality.")


if __name__ == "__main__":
    main()
