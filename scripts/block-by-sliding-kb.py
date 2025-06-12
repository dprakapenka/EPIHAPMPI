#!/usr/bin/env python
# make x kb block definition file
# by Dzianis Prakapenka
from __future__ import print_function
import os
import sys
import argparse
from epihap.utils import create_common_block_parser, prepare_output_directory, read_map_file
import itertools
import math
import re
from collections import defaultdict

def make_arg_parser():
    app_name = "block-by-sliding-kb.py"
    description = "generate block definition files based number of kbs, the minimum distance per block the start of each block slides by user defined number of snps"

    # Get the common parser
    parent_parser = create_common_block_parser()

    # Create a new parser that inherits from the parent parser
    parser = argparse.ArgumentParser(
        prog=app_name,
        description=description,
        parents=[parent_parser],  # Inherit arguments from parent_parser
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True # Ensure help message is added
    )

    # Add script-specific arguments
    parser.add_argument("-s", "--size",
                        default=100,
                        type=int,
                        help="minimum size of each block in kbs")
    parser.add_argument("-t", "--step",
                        default=1,
                        type=int,
                        help="number of snps between block starting points")
    parser.add_argument("--noheader",
                        action="store_true",
                        default=False,
                        help="use if map file has no header")
    return parser


if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    #name output if default
    if args.output=='hap_info':
            args.output='hap_info_' + str(args.size) + 'skb'
    if args.verbose: print('output path:', args.output)
    cwd = os.getcwd()
    outdir=cwd + "/" + args.output
    hap_block=outdir + '/hap_block_'
    hap_info=hap_block + 'info_'
    prepare_output_directory(outdir, args.verbose)

    if args.verbose:
        print('reading map file:', args.map)

    # Read map file using the utility function
    # read_map_file returns: dict where keys are chromosome identifiers (str) and
    # values are lists of [snp_index_within_chr, snp_id, position_int].
    snp_data_from_util = read_map_file(args.map, noheader=args.noheader)

    if not snp_data_from_util:
        # read_map_file already prints warnings, so just exit.
        sys.exit(f"Error: No SNP data loaded from map file {args.map}. Exiting.")

    block_idx = 0
    # Sort chromosome keys numerically if they are strings representing numbers, otherwise lexicographically.
    # The keys from read_map_file are strings.
    sorted_chromosomes = sorted(snp_data_from_util.keys(), key=lambda k: int(k) if k.isdigit() else k)

    for c in sorted_chromosomes:
        if args.verbose: print("writing to chromosome", c)

        try:
            hb=open(hap_block+str(c),'w')
        except IOError:
            print("can't open ", hap_block+str(c))
            continue # Skip to next chromosome if file can't be opened
        
        # snp_list_for_chrom_c is a list of [snp_index_within_chr, snp_id, position_int]
        snp_list_for_chrom_c = snp_data_from_util[c]

        if not snp_list_for_chrom_c:
            if args.verbose:
                print(f"Info: No SNPs found for chromosome {c}. Skipping.", file=sys.stderr)
            hb.close()
            continue

        for i,snp1_data in enumerate(snp_list_for_chrom_c[::args.step]):
            out = [i] 
            
            # max position of snp for the block
            # snp1_data[2] is position_int
            max_pos = snp1_data[2] + (args.size*1000) 
            
            start_index_for_snp2_search = (i * args.step) + 1

            for snp2_data in snp_list_for_chrom_c[start_index_for_snp2_search:]:
                # snp2_data[2] is position_int for the potential end SNP
                if snp2_data[2] >= max_pos:
                    # snp2_data[0] is snp_idx_within_chr for the end SNP
                    end_snp_original_idx = snp2_data[0] 
                    out.append(end_snp_original_idx)
                    print("blk_" + str(block_idx), "\t".join([str(o) for o in out]), file = hb)
                    block_idx+=1
                    break # Found the end for this block, move to next snp1_data
        hb.close()
