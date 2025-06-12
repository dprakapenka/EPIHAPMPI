#!/usr/bin/env python3
# splits findhap haplotypes.txt output using chromosome.data by chromosomes
# by Dzianis Prakapenka
from __future__ import print_function
import sys
import os
import argparse
import collections

def make_arg_parser():
    app_name="split-by-chr.py"
    description="splits findhap genotypes.filled and haplotypes.txt \
            output using chromosome.data by chromosomes"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--hap",
            help="path to haplotypes.txt type file")
    parser.add_argument("--chr",
            default="chromosome.data",
            help="path to chromosome.data type file")
    parser.add_argument("--geno",
            default="genotypes.txt",
            help="path to genotypes.txt type file")
    parser.add_argument("--genofilled",
            default="genotypes.filled",
            help="path to genotypes.filled type file")
    parser.add_argument("--cutoff",
            type=int,
            help="(int) minimum number of snps in original genotype.txt")
    parser.add_argument("--genofolder",
            default="geno",
            help="name for the seperated geno folder")
    parser.add_argument("--hapfolder",
            default="hap",
            help="name for the seperated hap folder")
    parser.add_argument("--genoprefix",
            default="chr",
            help="prefix for the seperated geno files")
    parser.add_argument("--happrefix",
            default="chr",
            help="prefix for the seperated hap files")
    parser.add_argument("-m","--missing-value",
            default='9999',
            help="missing value symbol")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    cwd = os.getcwd()
    genodir=cwd + "/" + args.genofolder + "/"
    hapdir=cwd + "/" + args.hapfolder + "/"
    if not os.path.exists(genodir):
        os.makedirs(genodir)
    if args.hap is not None:
        if not os.path.exists(hapdir):
            os.makedirs(hapdir)

    #print(args)
    if args.verbose:
        print("\nverbose mdoe\t-----------------*\n")
        print("working dir:", cwd)
        if args.hap is not None: print("haplotype.txt file:", args.hap)
        print("chromosome file:", args.chr)
        print("genotype filled file:", args.genofilled)
        if args.cutoff:
            print("cutoff:", args.cutoff)
            print("genotype file:", args.geno)
        print("genotype directory:", genodir)
        if args.hap is not None:
        #if args.hap is not None or args.haplist is not None:
            print("haplotype directory:", hapdir)
        print()

    # open and read original snp number info
    num_snps=[]
    if args.cutoff is not None:
        if args.geno is not None:
            try:
                if args.verbose:
                    print("reading", args.geno)
                g=open(args.geno, "r")
            except IOError:
                print("can't open ", args.geno)
            for line in g:
                pair=[line.strip().split()[i] for i in [0,2]]
                if args.cutoff:
                    if int(args.cutoff) > int(pair[1]):
                        continue
                    else:
                        num_snps.append(pair[0])
                else:
                    num_snps.append(pair[0])
            g.close()
    else:
        try:
            if args.verbose:
                print("reading", args.genofilled)
            g=open(args.genofilled, "r")
        except IOError:
            print("can't open ", args.genofilled)
        num_snps=[line.strip().split()[0] for line in g]
        g.close()

    # open and read chromosome info
    try:
        if args.verbose:
            print("reading", args.chr)
        f=open(args.chr, "r")
    except IOError:
        print("can't open ", args.chr)
    # skip header
    next(f)
    # fill in [chr, num within, num overall] and collect header
    header=[]
    chr_info=[]
    for line in f:
        l=line.strip().split()[:5]
        #print(l)
        header.append([l[0],l[1],l[4]])
        chr_info.append([int(i) for i in l[1:4]])
    f.close()
    chr_list=set([c[0] for c in chr_info])
    #print(header)

    # Dictionary to hold open file handles
    geno_files = {}
    hap_files = {}

    if args.hap is not None:
        for c in chr_list:
            try:
                hap_files[c]=open(hapdir + args.happrefix + str(c), "w")
            except IOError:
                print("can't open ", hapdir + args.happrefix + str(c))
    for c in chr_list:
        try:
            geno_files[c]=open(genodir + args.genoprefix + str(c), "w")
        except IOError:
            print("can't open ", genodir + args.genoprefix + str(c))
        head = [h[0] for h in header if h[1]==str(c)]
        print("ID\t", ' '.join(head),file=geno_files[c])

    if args.verbose:
        print()
        print(len(chr_list), "chromosomes in chromosome file")
        print(len(chr_info), "snps in chromosome file")
        print()

    #open and read geno and hap files
    try:
        if args.verbose:
            print("reading", args.genofilled)
        g=open(args.genofilled, "r")
    except IOError:
        print("can't open ", args.genofilled)

    if args.hap is not None:
        try:
            if args.verbose:
                print("reading", args.hap)
            h=open(args.hap, "r")
        except IOError:
            print("can't open ", args.hap)
        if args.verbose:
            print("writing")
        for hapline in h:
            geno=next(g).strip().split()
            hap=hapline.strip().split()
            if (geno[0]==hap[0]):
                id=geno[0]
                if id in num_snps:
                    genotype=list(geno[-1])
                    haplotype=list(map(''.join, zip(*[iter(hap[-1])]*2)))
                    for c in chr_list:
                        #order=[i[-1]-1 for i in chr_info if (i[0]==c)]
                        gout=[genotype[i[-1]-1] for i in chr_info if (i[0]==c)]
                        hout=[haplotype[i[-1]-1] for i in chr_info if (i[0]==c)]
                        print(id, ' '.join(gout),file=geno_files[c])
                        print(id, ''.join(hout),file=hap_files[c])
            else:
                print("individuals", geno[0],hap[0],"don't match")
    else:
        for genoline in g:
            geno=genoline.strip().split()
            id=geno[0]
            if id in num_snps:
                genotype=list(geno[-1])
                #haplotype=list(map(''.join, zip(*[iter(hap[-1])]*2)))
                for c in chr_list:
                    gout=[genotype[i[-1]-1] for i in chr_info if (i[0]==c)]
                    print(id, ' '.join(gout),file=geno_files[c])

    for f in geno_files.values():
        f.close()

    for f in hap_files.values():
        f.close()

