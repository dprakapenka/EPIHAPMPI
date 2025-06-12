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
    app_name = "block-by-kb.py"
    description = "generate block definition files based on a list of positions for each block"
    
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
                        default=500,
                        type=int,
                        help="size of each block in kbp")
    return parser

def parseBlock(positions):
    #block = list(map(int, re.split("\D",positions)))
    block = re.split("\D",positions)
    return block[0], int(block[1]), int(block[2])

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    #name output if default
    if args.output=='hap_info':
            args.output='hap_info_' + str(args.size) + 'kbp'
    if args.verbose: print('output path:', args.output)
    cwd = os.getcwd()
    outdir=cwd + "/" + args.output
    hap_block=outdir + '/hap_block_'
    hap_info=hap_block + 'info_'
    prepare_output_directory(outdir, args.verbose)

    #change to kb
    args.size*=1000

    if args.verbose:
        print('reading:', args.map)
    
    # Read map file using the utility function
    # read_map_file returns: dict where keys are chromosome identifiers (str) and
    # values are lists of [snp_index_within_chr, snp_id, position_int].
    snp_map_data = read_map_file(args.map)

    if not snp_map_data:
        sys.exit(f"Error: Could not read or parse map file {args.map}. Exiting.")

    block_idx = 0
    # Sort chromosome keys numerically if they are strings representing numbers
    # The keys from read_map_file are strings as they appear in the file.
    # The original script treated them as integers.
    sorted_chromosomes = sorted(snp_map_data.keys(), key=lambda k: int(k) if k.isdigit() else k)

    for c_str in sorted_chromosomes:
        # Convert chromosome string key to int if it was originally treated as int for output file naming
        # However, it's safer to keep it as string if it might not always be numeric.
        # For file naming (hap_block+str(c)), using c_str directly is fine.
        if args.verbose: print("writing to chromosome", c_str)

        try:
            hb=open(hap_block+str(c_str),'w')
        except IOError:
            print("can't open ", hap_block+str(c_str))
            continue # Skip to next chromosome if file can't be opened
        
        blk = []
        max_pos = 0
        
        # snp_map_data[c_str] is a list of [snp_index_within_chr, snp_id, position_int]
        # Sort by position_int (item[2])
        for snp_item in sorted(snp_map_data[c_str], key=lambda item: item[2]):
            snp_idx_within_chr = snp_item[0]
            # snp_id_str = snp_item[1] # Not directly used in this logic block
            position_int = snp_item[2]

            # initialize first snp in block
            if not blk:
                blk.append(snp_idx_within_chr)
                max_pos = position_int + args.size
            # snp fits in block, save it
            elif (position_int < max_pos):
                blk.append(snp_idx_within_chr)
            # end of block, if it is not single snp write block
            elif (len(blk) > 1):
                out = [min(blk), max(blk)]
                print("blk_" + str(block_idx), "\t".join([str(o) for o in out]), file = hb)
                block_idx+=1
                blk = [snp_idx_within_chr]
                max_pos = position_int + args.size
            # throw out snp if it is the only one in the block
            else:
                if args.verbose: 
                    print("\nless than 2 snps in block")
                    print("\tlone snp at index", blk[0], "is thrown out")
                blk = [snp_idx_within_chr]
                max_pos = position_int + args.size
        
        # After iterating through all SNPs for a chromosome,
        # check if there's a pending block to write (especially the last block)
        if (len(blk) > 1):
            out = [min(blk), max(blk)]
            print("blk_" + str(block_idx), "\t".join([str(o) for o in out]), file = hb)
            block_idx+=1
        elif blk: # A single SNP was left in blk
             if args.verbose: 
                print("\nless than 2 snps in block at the end of chromosome processing")
                print("\tlone snp at index", blk[0], "is thrown out")
        hb.close()
