#!/usr/bin/env python
# by Dzianis Prakapenka

from __future__ import print_function
import argparse

def make_arg_parser():
    app_name="mege-mrk-hap"
    description="merge the marker effects and heritabilities\
                with haplotype heritabilities for graphing with SNPEVG"
    parser = argparse.ArgumentParser(prog=app_name,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-sh", "--snp_effects",
            default=argparse.SUPPRESS,
            type=str,
            #required=True,
            help="path to SNP effects and heritabilites file from GVCHAP")
    parser.add_argument("-hh", "--hap_heritabilities",
            default=argparse.SUPPRESS,
            type=str,
            required=True,
            help="path to haplotype heritabilites file from GVCHAP [required]")
    parser.add_argument("-hi", "--hapinfo",
            default=argparse.SUPPRESS,
            nargs='+',
            #type=list,
            required=True,
            help="path(s) to hap_info file(s) [required]")
    parser.add_argument("--map",
            default='map.txt',
            required=True,
            help="path to map file")
    parser.add_argument("-o","--output",
            default='combined.snpe',
            help="output folder")
    parser.add_argument("--nosort",
            action="store_false",
            default=True,
            help="set to not sort hap chr by number in filename")
    parser.add_argument("-V","--verbose",
            action="store_true",
            default=False,
            help="verbose output")
    return parser

def sort_files(tosort):
    # sorts list of strings by digits in the string
    if args.nosort:
        return (sorted(tosort,key=lambda a: get_digits(a)))
    else:
        return tosort

def get_digits(s):
    # get digits from a string
    return int(''.join(list(filter(str.isdigit, s))))

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    hap_info_files=sort_files(args.hapinfo)
    if args.verbose:
        print(len(hap_info_files), "haplotype definition files found")

    try:
        f=open(args.map, "r")
    except IOError:
        print("Can't open ", args.map)
    snp_data = [line.strip().split() for line in f]
    f.close()
    #header = snp_data[0]
    #snp_data.pop(0)
    header = snp_data.pop(0)

    if ("snp_effects" in args):
        try:
            f=open(args.snp_effects, "r")
        except IOError:
            print("Can't open ", args.snp_effects)
        header=header + f.readline().strip().split()[3:]
        for i,line in enumerate(f):
            l = line.strip().split()
            snp_data[i].extend(l[3:])
        f.close()

    try:
        f=open(args.hap_heritabilities, "r")
    except IOError:
        print("Can't open ", args.hap_heritabilities)
    hap_data = [line.strip().split()[1:] for line in f]
    f.close()
    hap_cols=len(hap_data[0])+1
    header=header + hap_data[0]
    hap_data.pop(0)

    try:
        g=open(args.output, "w")
    except IOError:
        print("Can't open ", args.output)

    print("\t".join(header),file=g)

    blk=0
    for blk_def_file in hap_info_files:
        tmp_str=blk_def_file.split("\\")
        tmp_str=tmp_str[-1].split("/")
        chrom=get_digits(tmp_str[-1])
        snp_list = [snp for snp in snp_data if snp[1] == str(chrom)]
        if args.verbose:
            print("from", blk_def_file, "assuming chr nubmer", chrom)
        try:
            print("opening ", blk_def_file)
            h=open(blk_def_file, "r")
        except IOError:
            print("Can't open ", blk_def_file)
        for i,line in enumerate(h):
            l = line.strip().split()[1:]
            for s in range(int(l[0]),int(l[-1])+1):
                print('\t'.join(snp_list[s] + hap_data[blk]), file = g)
            blk = blk + 1

        #  if there are no blocks in the chr just write the markers
        #print(i)
        #if (i==0):
        #    for snp in snp_list:
        #        print('\t'.join(snp + ["0" for s in range(1,hap_cols)]), file = g)
        h.close()
    g.close()
