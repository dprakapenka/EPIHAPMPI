#!/usr/bin/env python3

import os
import argparse
import numpy as np
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

GENO_DIR = "geno"
HAP_DIR = "hap"
MAP_FILE_TXT = "map.txt"
MAP_FILE_CSV = "map.csv"
RUNTIME_LOG = "runtime.log"
FAILED_CHROMOSOMES_FILE = "failed_chromosomes.txt"

def make_arg_parser():
    parser = argparse.ArgumentParser(description="Filter genotype and haplotype files by MAF.")
    parser.add_argument("--geno",
                        nargs='+',
                        help="Path to genotype chr files")
    parser.add_argument("--hap",
                        nargs='+',
                        help="Path to haplotype chr files")
    parser.add_argument("--chr",
                        default="chromosome.data",
                        help="Path to chromosome.data file")
    parser.add_argument("--maf",
                        type=float,
                        default=0.05,
                        help="Minimum MAF")
    parser.add_argument("--out",
                        default="out",
                        help="Output directory")
    parser.add_argument("--threads",
                        type=int,
                        help="Number of parallel processes (default: all available cores)")
    parser.add_argument("-m", "--missing-value",
                        default='9999',
                        help="Missing value symbol")
    parser.add_argument("-V", "--verbose",
                        action="store_true",
                        help="Verbose output")
    parser.add_argument("--continue-on-error",
                        action="store_true",
                        help="Continue processing other chromosomes if one fails")
    parser.add_argument("--map-format",
                        choices=["txt", "csv"],
                        default="txt",
                        help="Format for combined map output (txt or csv)")
    parser.add_argument("--version",
                        action="store_true",
                        help="Show script version and exit")
    return parser

def log(message, verbose, logfile=None):
    if verbose:
        print(message)
    if logfile:
        with open(logfile, 'a') as logf:
            logf.write(message + "\n")

def get_digits(s):
    return int(''.join(filter(str.isdigit, s)))

def prepare_output_directories(outdir, haplotypes_provided):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, GENO_DIR), exist_ok=True)
    if haplotypes_provided:
        os.makedirs(os.path.join(outdir, HAP_DIR), exist_ok=True)

def validate_file_counts(geno_files, hap_files):
    if hap_files and len(geno_files) != len(hap_files):
        raise ValueError("Number of genotype files must match number of haplotype files when both are provided.")

def read_genotype_file(filepath):
    with open(filepath) as f:
        header = f.readline().strip().split()[1:]
        snpsum = np.zeros(len(header), dtype=int)
        snpmissing = np.zeros(len(header), dtype=int)
        geno_lines = []

        for line in f:
            parts = line.strip().split()
            id_, snps = parts[0], parts[1:]
            geno_lines.append((id_, snps))

            for i, snp in enumerate(snps):
                if snp in ("0", "1", "2"):
                    snpsum[i] += int(snp)
                    snpmissing[i] += 1

    return header, geno_lines, snpsum, snpmissing

def read_chromosome_map(chr_file, chrnum):
    #TODO: store entire thing in memory, chr number is redundant
    chr_map = []
    with open(chr_file) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if int(parts[1]) == chrnum:
                # Keep only SNPID (col 0), Chromosome (col 1), and Position (col 4)
                chr_map.append((parts[0], int(parts[1]), int(parts[4])))
    return np.array(chr_map, dtype=object)

def open_output_files(out_geno_file, out_hap_file, hap_file):
    gout = open(out_geno_file, 'w')
    hout = open(out_hap_file, 'w') if hap_file else open(os.devnull, 'w')
    hap_f = open(hap_file) if hap_file else None
    return gout, hout, hap_f

def process_chromosome(chrnum, geno_file, hap_file, chr_file, outdir, maf, missing_value, verbose):
    genodir = os.path.join(outdir, GENO_DIR)
    hapdir = os.path.join(outdir, HAP_DIR)
    out_geno_file = os.path.join(genodir, f"chr{chrnum}")
    out_hap_file = os.path.join(hapdir, f"chr{chrnum}") if hap_file else None
    start_time = time.time()

    try:
        header, geno_lines, snpsum, snpmissing = read_genotype_file(geno_file)

        snpfreq = np.divide(snpsum, 2 * snpmissing, where=snpmissing > 0)
        valid_snps = np.logical_and(snpfreq >= min(maf, 1 - maf), snpfreq <= max(maf, 1 - maf))

        chr_map = read_chromosome_map(chr_file, chrnum)
        filtered_map = chr_map[valid_snps]
        #map_entries = [(snp[0], chrnum, snp[2]) for snp in filtered_map]

        gout, hout, hap_f = open_output_files(out_geno_file, out_hap_file, hap_file)

        gout.write("ID " + " ".join(filtered_map[:, 0]) + "\n")

        try:
            for i, (id_, snps) in enumerate(geno_lines):
                filtered_geno = [snps[j] for j, valid in enumerate(valid_snps) if valid]
                gout.write(f"{id_} {' '.join(filtered_geno)}\n")

                if hap_f:
                    hap_line = hap_f.readline().strip().split()
                    hap_snps = [hap_line[-1][j:j+2] for j in range(0, len(hap_line[-1]), 2)]
                    filtered_hap = ''.join([hap_snps[j] for j, valid in enumerate(valid_snps) if valid])
                    hout.write(f"{id_} {filtered_hap}\n")
        finally:
            gout.close()
            if hap_f:
                hout.close()
                hap_f.close()

        num_snps_retained = np.sum(valid_snps)
        total_snps = len(valid_snps)  # Total SNP count before filtering
        #minor_allele_freqs = np.minimum(snpfreq, 1 - snpfreq)
        #mean_maf = np.mean(minor_allele_freqs[valid_snps]) if np.any(valid_snps) else 0.0
        mean_af = np.mean(snpfreq[valid_snps]) if np.any(valid_snps) else 0.0
        duration = time.time() - start_time

        log( ( f"Chromosome {chrnum} processed - "
              f"{num_snps_retained} of {total_snps} SNPs retained - "
              f"Mean AF: {mean_af:.4f} - "
              f"Time: {duration:.2f}s"),
            verbose,
            os.path.join(outdir, RUNTIME_LOG) )


        #return map_entries
        return filtered_map

    except Exception as e:
        print(f"Error processing chromosome {chrnum}: {e}")
        raise

def write_combined_map(outdir, all_map_entries, format):
    if format == "txt":
        map_file = os.path.join(outdir, MAP_FILE_TXT)
        with open(map_file, 'w') as combined:
            combined.write("SNPID\tChr\tPos\n")
            #lambda x: (chrnum, int(pos))):
            for snpid, chrnum, pos in sorted(all_map_entries, key=lambda x: (x[1], x[2])):
                combined.write(f"{snpid}\t{chrnum}\t{pos}\n")
    elif format == "csv":
        map_file = os.path.join(outdir, MAP_FILE_CSV)
        with open(map_file, 'w') as combined:
            combined.write("SNPID,Chr,Pos\n")
            for snpid, chrnum, pos in sorted(all_map_entries, key=lambda x: (x[1], x[2])):
                combined.write(f"{snpid},{chrnum},{pos}\n")

def main():
    args = make_arg_parser().parse_args()

    if args.version:
        print("MAF Filter Script\nVersion: 1.1")
        return

    outdir = os.path.abspath(args.out if args.out != "out" else f"out.maf{args.maf}")

    start_time = datetime.datetime.now()
    log(f"Job started at {start_time}", args.verbose, os.path.join(outdir, RUNTIME_LOG))

    prepare_output_directories(outdir, bool(args.hap))
    validate_file_counts(args.geno, args.hap)

    geno_files = {get_digits(os.path.basename(f)): f for f in args.geno}
    hap_files = {get_digits(os.path.basename(f)): f for f in args.hap} if args.hap else {}

    chromosomes = sorted(geno_files.keys())
    all_map_entries = []
    failed_chromosomes = []

    num_threads = args.threads or os.cpu_count()
    log(f"Using {num_threads} threads", args.verbose, os.path.join(outdir, RUNTIME_LOG))

    with ProcessPoolExecutor(max_workers=num_threads) as executor:

        futures = {
            executor.submit(
                process_chromosome, chrnum, geno_files[chrnum], hap_files.get(chrnum),
                args.chr, outdir, args.maf, args.missing_value, args.verbose
            ): chrnum for chrnum in chromosomes
        }

        for future in (tqdm(as_completed(futures), total=len(futures)) if TQDM_AVAILABLE else as_completed(futures)):
            chrnum = futures[future]
            try:
                all_map_entries.extend(future.result())
            except Exception:
                failed_chromosomes.append(chrnum)
                if not args.continue_on_error:
                    raise

    if failed_chromosomes:
        with open(os.path.join(outdir, FAILED_CHROMOSOMES_FILE), 'w') as f:
            f.write("\n".join(map(str, failed_chromosomes)))

    write_combined_map(outdir, all_map_entries, args.map_format)

    end_time = datetime.datetime.now()
    log(f"Job finished at {end_time}", args.verbose, os.path.join(outdir, RUNTIME_LOG))
    log(f"Total runtime: {end_time - start_time}", args.verbose, os.path.join(outdir, RUNTIME_LOG))


if __name__ == '__main__':
    main()
