#!/usr/bin/env python
# generates haplotype block genotypes using process pool
# by Dzianis Prakapenka
import os
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="get-hap-geno.py",
            description="Creates haplotype genotype files based on block info files",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument("--hap",
                        required=True,
                        nargs='+',
                        help="Path(s) to hap chr file(s) [required]")
    parser.add_argument("--hapinfo",
                        required=True,
                        nargs='+',
                        help="Path(s) to hap_info file(s) [required]")
    parser.add_argument("-m", "--missing",
                        type=str,
                        help="Coding for missing alleles")
    parser.add_argument("--missing-value",
                        type=str,
                        default="-9999",
                        help="Value for missing alleles")
    parser.add_argument("-o", "--output",
                        default='hap_geno',
                        help="Output folder for haplotype genotypes")
    parser.add_argument("-p", "--output-prefix",
                        default='hap_geno',
                        help="Output prefix for genotype files")
    parser.add_argument("--nosort",
                        action="store_true",
                        help="Do not sort input files by filename number")
    parser.add_argument("--header",
                        action="store_true",
                        help="Set if there is a header in the hap input files")
    parser.add_argument("-V", "--verbose",
                        action="store_true",
                        help="Verbose output")
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
    blocks, haps = [], []

    # Read hap_info file for block data
    with open(hap_info_file, 'r') as hi:
        blocks = [line.strip().split()[1:] for line in hi]

    # Read hap file for haplotypes
    with open(hap_file, 'r') as h:
        if args.header:
            next(h)
        haps = [line.strip().split() for line in h]

    codes = []
    for blk in blocks:
        hap_dict, label = {}, 1
        for hh in haps:
            h = ''.join(hh[1:]) if len(hh) > 2 else hh[1]
            h_part = h[int(blk[0])*2:int(blk[-1])*2+2]
            parents = [h_part[::2], h_part[1::2]]
            for p in parents:
                if args.missing and args.missing in p:
                    hap_dict[p] = args.missing_value
                elif p not in hap_dict:
                    hap_dict[p] = str(label)
                    label += 1
        codes.append(hap_dict)

    # Write to output file
    with open(outfile, 'w') as o:
        header = ["ID"] + [f"hap_{chr_name}_{b}\thap_{chr_name}_{b}" for b in range(1, len(blocks) + 1)]
        print('\t'.join(header), file=o)

        for hh in haps:
            ind, h = (hh[0], ''.join(hh[1:])) if len(hh) > 2 else hh
            output = [ind]
            for c, blk in enumerate(blocks):
                h_part = h[int(blk[0])*2:int(blk[-1])*2+2]
                parents = [h_part[::2], h_part[1::2]]
                output.extend(codes[c][p] for p in parents)
            print('\t'.join(output), file=o)

    if args.verbose:
        print(f"Completed processing chromosome {chr_name}.")

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    hap_files = sort_files(args.hap, args.nosort)
    hap_info_files = sort_files(args.hapinfo, args.nosort)

    # Verify matching file lengths
    if len(hap_files) != len(hap_info_files):
        sys.exit("ERROR: The number of hap and hap_info files do not match.")

    # Output directory setup
    outdir = os.path.join(os.getcwd(), args.output)
    os.makedirs(outdir, exist_ok=True)

    # Process each pair of hap and hap_info files in parallel
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_file_pair, hap_file, hap_info_file, args, outdir)
            for hap_file, hap_info_file in zip(hap_files, hap_info_files)
        ]
        for future in futures:
            future.result()  # Ensure completion

    if args.verbose:
        print(f"All genotype files are written to {outdir}.")
