#!/usr/bin/env python3.9
# filters gvchap genotype, haplotype, and chromosome.data files by maf
# by Dzianis Prakapenka
from __future__ import print_function
import sys
import os
import argparse
import random
import copy
#from operator import add

def make_arg_parser():
    app_name="create.kfold.pheno.py"
    description="create k files with N/k random rows set to missing\
            value"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pheno",
            help="path to phenotype/data file")
    parser.add_argument("-n", "--column-name",
            help="name of the column if header is there")
    parser.add_argument("-c", "--column",
            type=int,
            default=2,
            help="name for the output folder")
    parser.add_argument("-k", "--kfolds",
            type=int,
            default=10,
            help="number of k-folds")
    parser.add_argument("-p", "--percent",
            type=int,
            default=10,
            help="percent of the population to use as validation")
    parser.add_argument("-m","--missing-value",
            default="-9999111",
            help="missing value symbol")

    parser.add_argument("--noheader",
            action="store_true",
            default=False,
            help="verbose output")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

if __name__ == "__main__":

    parser = make_arg_parser()
    args = parser.parse_args()

    cwd = os.getcwd()
    split_filename=args.pheno.split(".")
    out_file = split_filename.copy()
    if args.verbose:
        print("\nverbose mdoe\t-----------------*\n")
        print("working dir:", cwd)
        print("input file:", args.pheno)
        print("k-folds:", args.kfolds)
        print("validation percentage:", args.percent)
        out_file.insert(-1,"KFOLD")
        print("out pheno:", ".".join(out_file))

    try:
        if args.verbose:
            print("reading", args.pheno)
        with open(args.pheno, "r") as f:
            if not args.noheader:
                header = f.readline().split()
                if args.column_name is not None:
                    try:
                        args.column = int(header.index(args.column_name))
                    except ValueError:
                        print("column", args.column_name, "not found")
                    except:
                        print("something wrong in header")
                
                if args.verbose:
                    print("working on column:", header[args.column])
                header.insert(args.column+1, header[args.column]+"_val")
            new_col=args.column+1
            data = []
            for line in f:
                l = line.split()
                data.append(([*l[:new_col],*l[args.column:]]))

            num_ind = len(data)
            num_val = round(num_ind * (args.percent / 100))

            if args.verbose:
                print("number of records:", num_ind)
                print("number of validation records:", num_val)

            
            for k in range(1,args.kfolds+1):
                new_data = copy.deepcopy(data)
                out_file = split_filename.copy()
                out_file.insert(-1, str(k))
                out_filename = '.'.join(out_file)
                if args.verbose:
                    print(out_filename)

                val_sample = random.sample(range(num_ind), k=num_val)
                for v in val_sample:
                    new_data[v][new_col] = args.missing_value
                try:
                    if args.verbose:
                        print("writing", args.pheno)
                    with open(out_filename, "w") as g:
                        if not args.noheader:
                            print(' '.join(header), file=g)
                        for d in new_data:
                            print(' '.join(d), file=g)
                except IOError:
                    print("can't open ", geno_file)

    except IOError:
        print("can't open ", geno_file)

    sys.exit()
