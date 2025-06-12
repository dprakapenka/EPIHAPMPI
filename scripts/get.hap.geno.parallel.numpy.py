#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def make_arg_parser():
    parser = argparse.ArgumentParser(
        prog="get-hap-geno.py",
        description="Creates haplotype genotype files based on block info files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--hap", required=True, nargs='+', help="Path(s) to hap chr file(s) [required]")
    parser.add_argument("--hapinfo", required=True, nargs='+', help="Path(s) to hap_info file(s) [required]")
    parser.add_argument("-m", "--missing", type=str, help="Coding for missing alleles")
    parser.add_argument("--missing-value", type=str, default="-9999", help="Value for missing alleles")
    parser.add_argument("-o", "--output", default='hap_geno', help="Output folder for haplotype genotypes")
    parser.add_argument("-p", "--output-prefix", default='hap_geno', help="Output prefix for genotype files")
    parser.add_argument("--nosort", action="store_true", help="Do not sort input files by filename number")
    parser.add_argument("--header", action="store_true", help="Set if there is a header in the hap input files")
    parser.add_argument("-V", "--verbose", action="store_true", help="Verbose output")
    return parser

def get_digits(s):
    return int(''.join(filter(str.isdigit, s)))

def sort_files(tosort, nosort):
    return tosort if nosort else sorted(tosort, key=lambda a: get_digits(a))

def process_file_pair(hap_file, hap_info_file, args, outdir):
    hap_cnum = get_digits(os.path.basename(hap_file))
    hap_info_cnum = get_digits(os.path.basename(hap_info_file))
    chr_name = str(hap_cnum) if hap_cnum == hap_info_cnum else '_u_' + str(hap_cnum)

    outfile = os.path.join(outdir, f"{args.output_prefix}_{chr_name}")

    # Load blocks as numpy array
    blocks = np.loadtxt(hap_info_file, dtype=int, usecols=(1, 2))

    # Load haps as numpy array, skipping header if necessary
    haps = np.loadtxt(hap_file, dtype=str, skiprows=1 if args.header else 0)

    # Extract individuals and haplotypes
    individuals = haps[:, 0]
    haplotype_data = np.array([''.join(row[1:]) for row in haps], dtype=str)

    # Dictionary to store block allele codes
    codes = []
    for blk_start, blk_end in blocks:
        hap_dict, label = {}, 1
        h_parts = haplotype_data[:, blk_start * 2:(blk_end + 1) * 2]
        parent_haps = np.array([h_parts[:, ::2], h_parts[:, 1::2]])

        for parent in parent_haps:
            for allele in parent:
                if args.missing and args.missing in allele:
                    hap_dict[allele] = args.missing_value
                elif allele not in hap_dict:
                    hap_dict[allele] = str(label)
                    label += 1
        codes.append(hap_dict)

    # Write to output file
    with open(outfile, 'w') as o:
        header = ["ID"] + [f"hap_{chr_name}_{b}\thap_{chr_name}_{b}" for b in range(1, len(blocks) + 1)]
        print('\t'.join(header), file=o)

        for ind, hap in zip(individuals, haplotype_data):
            output = [ind]
            for c, (blk_start, blk_end) in enumerate(blocks):
                h_part = hap[blk_start * 2:(blk_end + 1) * 2]
                parents = [h_part[::2], h_part[1::2]]
                output.extend(codes[c][parent] for parent in parents)
            print('\t'.join(output), file=o)

    if args.verbose:
        print(f"Completed processing chromosome {chr_name}.")

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    hap_files = sort_files(args.hap, args.nosort)
    hap_info_files = sort_files(args.hapinfo, args.nosort)

    if len(hap_files) != len(hap_info_files):
        sys.exit("ERROR: The number of hap and hap_info files do not match.")

    outdir = os.path.join(os.getcwd(), args.output)
    os.makedirs(outdir, exist_ok=True)

    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_file_pair, hap_file, hap_info_file, args, outdir)
            for hap_file, hap_info_file in zip(hap_files, hap_info_files)
        ]
        for future in futures:
            future.result()  # Ensure completion

    if args.verbose:
        print(f"All genotype files are written to {outdir}.")
