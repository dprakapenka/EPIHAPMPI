#!/bin/bash -l
#SBATCH --job-name=test.reml.mpi
#SBATCH --account dairyxbreed
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 1-0
#SBATCH -o ./logs/reml.mpi.%J.%j.%N.out
#SBATCH -e ./logs/reml.mpi.%J.%j.%N.err

date

save_dir="saved"

if [ ! -d "${save_dir}" ]; then
    mkdir -p ${save_dir}
fi

module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load intel-oneapi-compilers

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "threads=" $MKL_NUM_THREADS
echo "threads=" $OMP_NUM_THREADS
echo "tasks="  $SLURM_NTASKS
#export MKL_VERBOSE=1

mpirun -n ${SLURM_NTASKS} ../bin/grm.snp -n 128 ../example/example.dat -o ${save_dir}/grm.mpi
../bin/grm.snp.3order ${save_dir}/grm.mpi
