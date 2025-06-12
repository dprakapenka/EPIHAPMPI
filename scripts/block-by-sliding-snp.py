#!/usr/bin/env python
# make x kb block definition file
# by Dzianis Prakapenka
from __future__ import print_function
import os
import sys
import argparse
import itertools
import math
import re
from collections import defaultdict

def make_arg_parser():
    app_name="block-by-sliding-snp.py"
    description="generate block definition files based number of\
            SNPs per, the start of each block slides by user defined\
            number of SNPs"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s","--size",
            default=2,
            type=int,
            help="size of each block in number of SNPs")
    parser.add_argument("-t","--step",
            default=1,
            type=int,
            help="number of SNPs between block starting points")
    parser.add_argument("-m","--map",
            default=argparse.SUPPRESS,
            required=True,
            help="path to map file [required]\n\t\
                    format: 'SNPID Chr Pos' with header") 
    parser.add_argument("--noheader",
            action="store_true",
            default=False,
            help="use if map file has no header")
    parser.add_argument("-o","--output",
            default='hap_info',
            help="output folder name")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser


if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    #name output if default
    if args.output=='hap_info':
            args.output='hap_info_' + str(args.size) + 'ssnp'
    if args.verbose: print('output path:', args.output)
    cwd = os.getcwd()
    outdir=cwd + "/" + args.output
    hap_block=outdir + '/hap_block_'
    hap_info=hap_block + 'info_'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.verbose:
            print('reading:',args.map)
    #open and read file
    try:
            f=open(args.map, "r")
    except IOError:
            print("can't open ", args.map)

    if not args.noheader:
        f.readline()
    snp_idx = 0
    snp_dict = {}
    current_chr = ""
    for line in f:
        l = line.strip().split()
        if (current_chr != l[1]):
            current_chr = l[1]
            snp_idx = 0
            snp_dict[l[1]]=[]
        snp_dict[l[1]].append([snp_idx, l[0], l[2]])
        snp_idx += 1

    block_idx = 0
    for c in sorted(snp_dict.keys()):
        if args.verbose: print("writing to chromosome", c)

        try:
            hb=open(hap_block+str(c),'w')
        except IOError:
            print("can't open ", hap_block+str(c))
        blk = []
        max_pos = 0
        max_snps = len(snp_dict[c])
        for i,snp in enumerate(snp_dict[c][::args.step]):
            end_snp = i + args.size
            if end_snp > max_snps:
                #end_snp = max_snp
                # last block ends with the last snp (no blocks smaller than block size)
                break
            out = [snp[0], snp_dict[c][end_snp-1][0]]
            print("blk_" + str(block_idx), "\t".join([str(o) for o in out]), file = hb)
            block_idx+=1

        hb.close()
