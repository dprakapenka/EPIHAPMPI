#!/usr/bin/env python3
# make x snp block definition file
# by Dzianis Prakapenka
from __future__ import print_function
import sys, os
import argparse
import itertools
import math

def make_arg_parser():
    app_name="hap-block-stat"
    description="calculate statistics from hap_info files"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m","--map",
            default="map.txt",
            required=True,
            help="path to map file [required]")
    parser.add_argument("-i","--hap_info",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="path(s) to hap_info chr file(s) [required]")
    parser.add_argument("-o","--output",
            default=False,
            help="output folder")
    parser.add_argument("--nosort",
            action="store_false",
            default=True,
            help="set to not sort hap chr by number in filename")
    parser.add_argument("--noheader",
            action="store_false",
            default=True,
            help="do not output header")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

def sort_files(tosort):
    # sorts list of strings by digits in the string
    if args.nosort:
        return tosort
    else:
        return (sorted(tosort,key=lambda a: get_digits(a)))

def get_digits(s):
    # get digits from a string
    return int(''.join(list(filter(str.isdigit, s))))

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    file_list = sort_files(args.hap_info)
    if args.verbose: print('file list', file_list, len(file_list))

    allmap = []

    try:
        if args.verbose: print('reading', args.map)
        f=open(args.map,'r')
        f.readline()
        for line in f:
            allmap.append(line.strip().split())
        f.close()

    except IOError:
        print("can't open ", args.map)

    blocks = []
    for cf in file_list:
        cnum = get_digits(cf.split('/')[-1])
        map_local = [m for m in allmap if (int(m[1]) == cnum)]
        if args.verbose: print('chr ', cnum, 'len:' , len(map_local))

        try:
            if args.verbose: print('processing', cf)
            f=open(cf,'r')
            for line in f:
                l = [int(b) for b in line.strip().split()[-2:]]
                # if block exceeds size of chromosome, limit it
                if (l[1] >= len(map_local)):
                    l[1] = len(map_local) - 1
                # if the block is 0, skip it
                if l[0] == l[1]:
                    continue
                distance = int(map_local[l[1]][2]) - int(map_local[l[0]][2])
                num_in_block = l[1]-l[0]+1
                if args.verbose:
                    print(cnum, l, num_in_block, distance)
                blocks.append([num_in_block, distance])
            f.close()

        except IOError:
            print("can't open ", cf)

    num_blocks = len(blocks)
    min_snp = blocks[0][0]
    max_snp = blocks[0][0]
    avg_snp = 0
    min_dist = blocks[0][1]
    max_dist = blocks[0][1]
    avg_dist = 0

    for b in blocks:
        avg_snp = avg_snp + b[0]
        if b[0] < min_snp:
            min_snp = b[0]
        if b[0] > max_snp:
            max_snp = b[0]

        avg_dist = avg_dist + b[1]
        if b[1] < min_dist:
            min_dist = b[1]
        if b[1] > max_dist:
            max_dist = b[1]

    avg_snp = avg_snp / num_blocks
    avg_dist = avg_dist / num_blocks

    header = ["num_blocks",
            "min_snp",
            "max_snp",
            "avg_snp",
            "min_dist",
            "max_dist",
            "avg_dist"]
    result = [str(num_blocks),
            str(min_snp),
            str(max_snp),
            str(avg_snp),
            str(min_dist),
            str(max_dist),
            str(avg_dist)]

    if args.noheader: print(','.join(header))
    print(','.join(result))
    if args.output:
        try:
            if args.verbose: print('reading', args.output)
            g=open(args.output,'w')
            if args.noheader: print(','.join(header), file=g)
            print(','.join(result), file=g)
            f.close()

        except IOError:
            print("can't open ", args.output)






