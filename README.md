# EPIHAPMPI

## Introduction

EPIHAPMPI is a distributed-memory parallel computing package for genomic prediction and variance component estimation of large datasets. 

Supported models include any combination of:
- SNP additive (A) and dominance (D) effects
- Pairwise epistasis effects (AxA, AxD, DxD)
- Thirds order pairwise epistasis effects (AxAxA, AxAxD, AxDxD, DxDxD)
- Haplotype effects (HA)
- Single chromosome or region GRMs of any of the above types

Estimates variance components with GREML, calculates GBLUP and heritabilites.

## Install

For full step-by-step install instructions, see **[Installation section](INSTALL.md)**.

To quickly build the main executables:

~~~bash
# load compiler & MPI modules
# the names of the modules are system dependent
# for SciNet users:
module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load intel-oneapi-compilers

# or set the environment variables from intel
# on your system:
source /opt/intel/oneapi/setvars.sh

# build the reml and grm targets
make reml grm
~~~

## Run

For more in-depth instructions, see **[Running epihap.mpi section](manual/03.running.md)**.

~~~bash
# set threads to match your SLURM allocation
# example:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -n ${SLURM_NTASKS} ./bin/grm.snp \
                        --input genofile \
                        --out grmprefix

srun -n ${SLURM_NTASKS} ./bin/reml \
                        --grm grmprefix \
                        --pheno example.phen \
                        --pcol 2 \
                        --va 2 --vd 3.1 --vaa 1.3 --ve 1.4e1 \
                        --out project
~~~

## Dependencies

This software relies on external libraries that are **not** included in this repository:

- [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/articles/license/onemkl-license.html)

Or (coming soon)

- MPI (e.g. OpenMPI)
- BLAS/LAPACK
- [ScaLAPACK](http://www.netlib.org/scalapack/)
Blackford _et al._, *ScaLAPACK Users' Guide*, SIAM, 1997.  

Users are responsible for obtaining and complying with any licenses for these dependencies.

## License

This software is licensed under the GNU General Public License v3.0 (GPLv3). See the [LICENSE](./LICENSE) file details.

## How to Cite

If you use this software in published research, you are **required to cite** the following:

https://github.com/dprakapenka/epihap.mpi

You can also find citation metadata in the [CITATION.cff](./CITATION.cff) file.

Failure to cite this work when using it in publications may be considered a breach of academic ethics.
