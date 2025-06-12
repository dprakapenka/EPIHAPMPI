#!/usr/bin/env python3
import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator  # automatic ticker
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FuncFormatter
import numpy as np

def get_dynamic_top_margin(title_fontsize=16, base_margin=0.98):
    """
    Dynamically adjusts the top margin based on the title font size.
    - base_margin (default=0.97) sets the margin for a normal font size (e.g., 16).
    Returns a rect[3] value to use in plt.tight_layout(rect=[0, 0, 1, dynamic_top]).
    """
    default_fontsize = 18  # Baseline fontsize
    scaling_factor = 0.002  # Controls how much space increases per font size unit
    adjusted_margin = base_margin - (title_fontsize - default_fontsize) * scaling_factor
    
    # Ensure top margin is within reasonable limits
    return max(0.92, min(adjusted_margin, 0.99))

def parse_epi_log(file_path):
    """
    Parses the epi.log file to extract iteration data while ensuring the first iteration is included even if some values are missing.
    """
    iter_header_re = re.compile(r'Iteration:\s*(\d+)', re.IGNORECASE)
    var_re = re.compile(r'\b(V[A-Z]+)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)')
    herit_re = re.compile(r'\b(h2_[A-Z]+)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)')
    reml_re = re.compile(r'^(EM-REML|AI-REML)', re.IGNORECASE)
    # don't include total iteration time
    time_re = re.compile(r'The time \(seconds\) cost for this iteration is\s*:\s*(\d+)', re.IGNORECASE)
    #time_re = re.compile(r'The time.*?(\d+)', re.IGNORECASE)
    
    iterations = []
    current_iter = {}
    in_iterations = False  
    reml_method = None
    
    with open(file_path, 'r') as f:
        for line in f:
            # Identify a new iteration
            iter_match = iter_header_re.search(line)
            if iter_match:
                if in_iterations and current_iter:  
                    iterations.append(current_iter)
                in_iterations = True  
                current_iter = {"Iteration": int(iter_match.group(1))}  
            
            if not in_iterations:
                continue  

            # Detect REML method
            reml_match = reml_re.search(line)
            if reml_match:
                reml_method = reml_match.group(1)
                current_iter["Method"] = reml_method
            
            # Search for variance components
            for match in var_re.finditer(line):
                current_iter[match.group(1)] = float(match.group(2))
            
            # Search for heritability components
            for match in herit_re.finditer(line):
                current_iter[match.group(1)] = float(match.group(2))
            
            # Search for iteration time
            time_match = time_re.search(line)
            if time_match:
                current_iter["Time"] = float(time_match.group(1))
    
    if current_iter:
        iterations.append(current_iter)
    
    return pd.DataFrame(iterations)

def identify_large_variances(df, threshold=5):
    """
    Identifies variance columns that are significantly larger at each iteration.
    """
    variance_cols = [col for col in df.columns if col.startswith('V')]
    large_variances = set()
    
    for _, row in df.iterrows():
        median_var = row[variance_cols].median(skipna=True)  # Compute median per iteration
        for col in variance_cols:
            if row[col] > threshold * median_var:
                large_variances.add(col)
    
    return list(large_variances)

def scale_selected(df, threshold=5):
    """
    Scales selected variance columns to the next whole number greater than the max of the remaining variances.
    """
    large_variances = identify_large_variances(df, threshold)
    if not large_variances:
        return df, {}  # No scaling needed
    
    variance_cols = [col for col in df.columns if col.startswith('V')]
    non_scaled_cols = [col for col in variance_cols if col not in large_variances]
    max_non_scaled = df[non_scaled_cols].max().max(skipna=True)
    scale_factor = np.ceil(max_non_scaled) if max_non_scaled > 0 else 1
    scale_factors = {}
    
    for col in large_variances:
        max_large = df[col].max(skipna=True)
        if max_large > 2 * max_non_scaled:  # Ensure it's at least 2x larger before scaling
            df[col] = df[col] / (max_large / scale_factor)
            scale_factors[col] = scale_factor
    
    return df, scale_factors

def label_reml_transitions(ax, df):
    """ Adds vertical lines and labels at REML method transitions. """
    transitions = df[df["Method"] != df["Method"].shift()]
    for _, row in transitions.iterrows():
        iteration = row["Iteration"]
        # Change "EM-REML" to "EM", "AI-REML" to "AI"
        method = row["Method"].replace("-REML", "")  
        ax.axvline(iteration, color='red', linestyle=':', alpha=0.7)
        ax.text(iteration, df["Time"].max(), method, color='red', fontsize=10, rotation=90, ha='right', va='top', alpha=0.8)

def shade_reml_regions(ax, df):
    """ Shades EM-REML and AI-REML regions in different colors and adds a legend. """
    previous_iteration = df["Iteration"].iloc[0]
    previous_method = df["Method"].iloc[0]

    em_color = 'lightcoral'
    ai_color = 'lightblue'

    alpha=0.2

    for i in range(1, len(df)):
        current_iteration = df["Iteration"].iloc[i]
        current_method = df["Method"].iloc[i]

        color = em_color if previous_method == "EM-REML" else ai_color
        ax.axvspan(previous_iteration, current_iteration, facecolor=color, alpha=alpha)

        previous_iteration = current_iteration
        previous_method = current_method

    # Shade the last region until the final iteration
    final_color = em_color if previous_method == "EM-REML" else ai_color
    ax.axvspan(previous_iteration, df["Iteration"].iloc[-1], facecolor=final_color, alpha=alpha)

    # return legend handles for shaded regions
    return [
        plt.Rectangle((0, 0), 1, 1, color=em_color, alpha=alpha, label="EM-REML"),
        plt.Rectangle((0, 0), 1, 1, color=ai_color, alpha=alpha, label="AI-REML"),
    ]


def plot_iteration_data(df, threshold=5, debug=False, title_fontsize=18,save=None):

    # set monospace font
    plt.rc('font', family=['DejaVu Sans Mono', 'Courier New', 'Monospace'])
    #plt.rc('font', family=['Courier New', 'DejaVu Sans Mono', 'Monospace'])

    """Creates three vertically stacked plots: iteration time, variances, and heritabilities."""
    if df.empty:
        print("No valid data found.")
        return
    
    df, scale_factors = scale_selected(df, threshold)
    if debug:
        print(df)
        print(scale_factors)
    
    print(df)

    x = df["Iteration"]
    # Find reml method change rows
    transitions = df[df["Method"] != df["Method"].shift()].copy()
    var_cols = [col for col in df.columns if col.startswith('V')]
    herit_cols = [col for col in df.columns if col.startswith('h2_')]
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [2, 2, 1]})
    fig.subplots_adjust(left=0.15)  # Increase left margin dynamically
    dynamic_top = get_dynamic_top_margin(title_fontsize)

# --- Plot heritabilities (Top Graph) ---
    ax1 = axes[0]
    for col in herit_cols:
        ax1.plot(x, df[col], marker='o', label=col)
    ax1.set_ylabel("Heritability Estimates")
    ax1.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    ax1.grid(axis='y', linestyle='--', alpha=0.3)  # Faint horizontal grid lines

# --- Plot Variances (Middle Graph) ---
    ax2 = axes[1]
# Define fixed widths
    max_col_width = max(len(col) for col in var_cols)  # Longest variance column name
    right_align_width = 8  # Total width for (x1) or (x1/N)
    for col in var_cols:
        padded_col = col.ljust(max_col_width + 2)  # Extra spacing for separation
        if col in scale_factors:
            scale_text = f"(x1/{int(scale_factors[col])})"
        else:
            scale_text = "(x1)"
		
# Right-align the scale factor label
        label = f"{padded_col}{scale_text.rjust(right_align_width)}"
        ax2.plot(x, df[col], marker='o', label=label)

    ax2.set_ylabel("Variance Estimates (Scaled)")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    ax2.grid(axis='y', linestyle='--', alpha=0.3)  # Faint horizontal grid lines

# --- Plot Iteration Time (Bottom Graph) ---
    ax3 = axes[2]
    ax3.plot(x, df["Time"], color='black', marker='x', linestyle='--', label="Time (s)")

    # average of the time per iteration and label it
    avg_time = df["Time"].mean()
    avg_line = ax3.axhline(avg_time, color='brown', linestyle='--', alpha=0.5, label="Avg Time")

# Choose method
    method_display = "shade"
# Collect existing legend handles
    handles, labels = ax3.get_legend_handles_labels()
# Apply either labels or shading
    if method_display == "shade":
        legend_handles = shade_reml_regions(ax3, df)  # Get legend handles for shading
        if legend_handles:  # Ensure shading legend exists
            #handles = list(handles) + legend_handles
            #labels = list(labels) + ["EM-REML", "AI-REML"]
            handles.extend(legend_handles)
            labels.extend(["EM-REML", "AI-REML"])
    elif method_display == "label":
        label_reml_transitions(ax3, df)


    ax3.text(x.iloc[-1], avg_time, f"{avg_time:.2f}s", color='brown', fontsize=10, va='bottom', ha='right', alpha=0.7)

    ax3.set_xlabel("Iteration")
    ax3.set_ylabel("Iteration Time (s)")

    ax3.grid(axis='y', linestyle='--', alpha=0.3)  # Faint horizontal grid lines
    #ax3.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    ax3.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.05, 1), frameon=True)

# x-axis ticks and labels:
# Ensure all iterations have tick marks as **minor ticks**
    ax3.xaxis.set_minor_locator(FixedLocator(df["Iteration"]))

# Use MaxNLocator to determine which ticks should be **major ticks** (labeled)
    major_locator = MaxNLocator(integer=True, steps=[1, 2, 5, 10], min_n_ticks=20, prune="lower")
    ax3.xaxis.set_major_locator(major_locator)

# Get selected major ticks **only once**
    major_ticks = ax3.xaxis.get_majorticklocs()

# Function to only format labels for major ticks
    ax3.xaxis.set_major_formatter(FuncFormatter(lambda value, _: str(int(value)) if value in major_ticks else ""))


# y-axis ticks and labels:
    # Automatically align y-axis labels across subplots
    fig.align_ylabels()

# title:
    plt.suptitle("Iteration Time, Variance, and Heritability Trends")

# plot layout:
	# rect: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0, 1, dynamic_top])
# plot show or save to file:
    if save:
        plt.savefig(save, dpi=300)
        print(f"Plot saved to {save}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Process and visualize iteration data from epi.log.")
    parser.add_argument("-i", "--input", default="epi.log", help="Input file (default: epi.log)")
    parser.add_argument("-t", "--threshold", type=float, default=5, help="Threshold for variance scaling (default: 5)")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    parser.add_argument("-f", "--fontsize", type=int, default=18, help="Title font size (default: 18)")
    parser.add_argument("-s", "--save", metavar="FILE", help="Save plot to file (fi1.png, plot3.svg, etc) instead of displaying")
    parser.add_argument("--title", type=str, default="Iteration Time, Variance, and Heritability Trends",
                        help="Custom title for the plot (default: 'Iteration Time, Variance, and Heritability Trends')")

    args = parser.parse_args()
    
    df = parse_epi_log(args.input)
    plot_iteration_data(df, threshold=args.threshold, debug=args.debug, title_fontsize=args.fontsize, save=args.save)


if __name__ == "__main__":
    main()
    #TODO: add option to change title
