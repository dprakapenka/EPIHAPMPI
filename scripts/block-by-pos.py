#!/usr/bin/env python
# make x kb block definition file
# by Dzianis Prakapenka
from __future__ import print_function
import os
import sys
import argparse
from epihap.utils import prepare_output_directory, read_map_file
import itertools
import math
import re
from collections import defaultdict

def make_arg_parser():
    app_name="block-by-pos.py"
    description="generate block definition files based on a list of\
            positions for each block"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p","--pos",
            default=argparse.SUPPRESS,
            required=True,
            help="path to block position file [required]\n\t\
                    format: 'chr:begin_pos:end_pos' without header") 
    parser.add_argument("-m","--map",
            default=argparse.SUPPRESS,
            required=True,
            help="path to map file [required]\n\t\
                    format: 'SNPID Chr Pos' with header") 
    parser.add_argument("-o","--output",
            default='hap_info',
            help="output file name")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
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
            args.output='hap_info_pos'
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
    snp_map_data = read_map_file(args.map)

    if not snp_map_data:
        # read_map_file already prints warnings, so just exit.
        sys.exit(f"Error: No SNP data loaded from map file {args.map}. Exiting.")

    if args.verbose:
        print('reading block positions from:',args.pos)
    #open and read file
    try:
        g=open(args.pos, "r")
    except IOError:
        print("can't open ", args.pos)

    chromosomes = {}
    block_idx = 0
    for line in g:
        parsed_block_info = parseBlock(line)
        if not parsed_block_info:
            # parseBlock prints its own error, skip this line
            continue 
        c, blk_begin, blk_end = parsed_block_info
        
        blk = []
        
        # Get the list of SNPs for the current chromosome c
        # The keys in snp_map_data are strings, parseBlock also returns c as string.
        snp_list_for_chrom_c = snp_map_data.get(c)
        
        if not snp_list_for_chrom_c:
            if args.verbose:
                print(f"Warning: No SNPs found in map file for chromosome {c} specified in {args.pos} line: {line.strip()}", file=sys.stderr)
            continue # Skip to the next line in the pos file

        for snp_item in snp_list_for_chrom_c:
            # snp_item is [snp_index_within_chr, snp_id_str, position_int]
            position_int = snp_item[2]
            snp_index_within_chr = snp_item[0]
            
            if blk_begin <= position_int <= blk_end:
                blk.append(snp_index_within_chr)
        
        if len(blk) > 1:
            # Store the min and max SNP index for this block
            chromosomes.setdefault(c, {})[block_idx] = [min(blk), max(blk)] # Already sorted by append, just need min/max
            block_idx += 1
        else:
            if args.verbose or len(blk) == 1 : # Print warning if a block from pos file results in <2 SNPs
                print(f"Warning: Less than 2 SNPs found for block defined by '{line.strip()}' in chromosome {c}. This block will be skipped.", file=sys.stderr)

    # Sort chromosome keys for consistent output order if they represent numbers
    # Original script did not explicitly sort chromosome keys before writing output files.
    # Adding sorting for consistency.
    sorted_chromosome_keys = sorted(chromosomes.keys(), key=lambda k: int(k) if k.isdigit() else k)

    for c in sorted_chromosome_keys:
        if args.verbose: print("writing to chromosome", c)

        try:
            hb=open(hap_block+str(c),'w')
        except IOError:
            print("can't open ", hap_block+str(c))

        for b in chromosomes[c]:
            print("blk_"+str(b), "\t".join([str(s) for s in chromosomes[c][b]]), file=hb)

        hb.close()
