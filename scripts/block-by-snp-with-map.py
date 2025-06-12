#!/usr/bin/env python
# make x snp block definition file
# by Dzianis Prakapenka
from __future__ import print_function
import sys, os
import argparse
import itertools
import math

def make_arg_parser():
    app_name="block-by-snp"
    description="make x snp block definition file using map file as input"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m","--map",
            default="map.txt",
            required=True,
            help="path to map file [required]")
    parser.add_argument("--snp",
            default=4,
            type=int,
            help="block size in snp [required]")
    parser.add_argument("-o","--output",
            default='hap_info',
            help="output folder")
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

    file_list = sort_files(args.input)

    if args.output=='hap_info':
        args.output='hap_info_snp_'+str(args.snp)
    cwd = os.getcwd()
    outdir=cwd + "/" + args.output
    hap_block=outdir + '/hap_block_'
    hap_info=hap_block + 'info_'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    sd=int(args.snp)

    if args.verbose:
        print('block size:',args.snp,'snp')
    #open and read file

    count = 1
    current_chr=0
    for cf in file_list:
        cnum = get_digits(cf.split('/')[-1])
        unknown=''
        print(cf, cnum)
        if (str(cnum).isdigit()):
             chr_name=cnum
        else:
            current_chr=current_chr+1
            chr_name='_u_'+ str(current_chr)
            print('WARNING, not sure which chromosome is', cf)
            print('assigning the name', chr_name, 'fix it before proceeding')

        try:
            f=open(cf,'r')
        except IOError:
            print("can't open ", cf)
        if args.verbose: print('processing', cf)
        if (not args.noheader):
            f.readline()
        hap = f.readline().strip().split()
        f.close()
        if len(hap)<=2:
            num_hap = len(hap[-1])/2
        else:
            num_hap = (len(hap)-1)/2
        rem=num_hap%sd
        blks=(num_hap-rem)/sd
        if args.verbose:
            print('\t', num_hap, 'haplotype snps')
            print('\t', blks, 'haplotype blocks plus', rem, 'snps')

        try:
            hb=open(hap_block+str(chr_name),'w')
        except IOError:
            print("can't open ", hap_block+str(chr_name))
        try:
            hbi=open(hap_info+str(chr_name),'w')
        except IOError:
            print("can't open ", hap_info+str(chr_name))

        for i in range(0,blks*sd,sd):
            print('blk' + str(count) + '\t',i, i-1+sd,file=hb)
            print('blk' + str(count) + '\t','\t'.join([str(ii) for ii in range(i,i+sd)]),file=hbi)
            count=count+1
        if (rem != 0):
            if args.verbose: print('\t last block ', rem)
            print('blk' + str(count) + '\t',num_hap-rem, num_hap,file=hb)
            print('blk' + str(count) + '\t','\t'.join([str(ii) for ii in range(num_hap-rem, num_hap)]),file=hbi)
            count=count+1
        if args.verbose: print('\t', count-1, 'blocks in', cf)
