#!/usr/bin/env python
# Counts haps in chr files (columns minus id header)
# by Dzianis Prakapenka
from __future__ import print_function
import argparse
from collections import Counter

def make_arg_parser():
    app_name="hap-stat.py"
    description="statistics about hap input files of gvchap"
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
            default='hap-stat.log',
            help="output file name")
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
    counters = []
    total_blocks = 0
    for cf in file_list:
        #open and read file
        try:
            f=open(cf, "r")
            num_blocks = int((len(f.readline().strip().split())-1)/2)
            #print("blocks", num_blocks)
            haps = [line.strip().split()[1:] for line in f]
            f.close()
            total_blocks += num_blocks
            a_count = 0
            for blk in zip(*haps):
                cnt = Counter()
                for h in blk:
                    cnt[h] += 1

                if a_count == 0:
                    a_count += 1
                    counters.append(cnt)
                else:
                    a_count = 0
                    counters[-1] = counters[-1] + cnt

            #for b in range(num_blocks):
        except IOError:
            print("can't open ", cf)

    if (len(counters) != total_blocks):
        print("error num block counters != num blocks", len(counters), total_blocks)

    total_haps = 0
    for c in counters:
        total_haps += len(list(c))
    haps_per_block = total_haps / float(total_blocks)

    print(total_haps, haps_per_block, total_blocks)
