#!/usr/bin/env python
# Counts snps in chr files (columns minus id header)
# by Dzianis Prakapenka
from __future__ import print_function
import argparse

def make_arg_parser():
    app_name="count-snps.py"
    description="counts snps in chr file"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="path to input files [required]") 
    parser.add_argument("-o","--output",
            default='count-snps.log',
            help="output file name")
    parser.add_argument("-p","--prefix",
            default='geno_snp',
            help="prefix in front of output lines")
    parser.add_argument("--nosort",
	    action="store_true",
	    default=False,
	    help="do not sort by number in filename")


    return parser

def get_int(name):
    return int(''.join(list(filter(str.isdigit, name))))

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()


    if not args.nosort:
      file_list = (sorted(args.input,key=lambda a: int(''.join(list(filter(str.isdigit, a))))))
    else:
      file_list = args.input
    count = []
    for cf in file_list:
        #open and read file
        try:
            f=open(cf, "r")
        except IOError:
            print("can't open ", cf)
        count.append(len(f.readline().strip().split())-1)
        f.close()
    for cf,c in zip(file_list,count):
        print(args.prefix,c,cf)
