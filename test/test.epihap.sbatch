#!/bin/bash
#SBATCH -p atlas
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=0
#SBATCH -A dairyxbreed
#SBATCH -t 1-0
#SBATCH -o logs/epihap.out.%j.%N
#SBATCH -e logs/epihap.err.%j.%N
date
echo $SLURM_CPUS_ON_NODE
module load intel-oneapi-mkl

./run.epihap.sh
