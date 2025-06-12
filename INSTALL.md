# Installation

This guide shows how to compile and run the **reml** and **grm**, and other applications using the provided Makefile. GPU support is under active development and described in the final section.

---

## Prerequisites

1. **Modules / Environment**  
```bash
# if using on a cluster, load Intel oneAPI components, e.g.:
module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load intel-oneapi-compilers
# or source the environment variables from your
# installation. replace the path to setvars.sh
# if needed
source /opt/intel/oneapi/setvars.sh
```
2. **Make sure** `$MKLROOT` is defined:  
```bash
echo "$MKLROOT"
/opt/intel/oneapi/mkl/2025.0
```
3. **MPI implementation** installed and in your `PATH`.

---

## Directory Layout

```
epihap.mpi/
|-- src/                     # source code
|   |-- reml.cpp
|   |-- options.cpp
|   |-- ...
|-- obj/                     # object files
|-- bin/                     # executables built here
|-- include/                 # header files
|-- manual/                  # markdown documentation files
|-- scripts/                 # helper scripts
|-- test/                    # example test folder
|-- tests/                   # unit tests
|-- example/                 # example input files
|-- Makefile                 # makefile for compilation
```

---

## Makefile Targets

- **`make reml`**  
Builds the CPU-only version (`bin/reml`) for variance component estimation, heritabilites, and GBLUP.
- **`make grm`**  
Builds the CPU-only version (`bin/grm.snp`) for second order snp genomic relationship matrices(GRMs).
- **`make grm.hap`**  
Builds the CPU-only version (`bin/grm.hap`) for haplotype additive genomic relationship matrix (GRM).
- **`make grm.3order`**  
Builds the CPU-only version (`bin/grm.snp.3order`) for third order genomic relationship matrices (GRMs).
- **`make rank`**  
Builds the CPU-only version (`bin/rank`) calcualtes the rank of given GRMs.
- **`make read.bin`**  
Builds the CPU-only version (`bin/read.bin`) displays a portion of the binary GRM.
- **`make manual`**  
Generates a pdf version of the manual from .md and manual/\*.md files in this repository
- **`make clean`**  
Deletes object files and bianries

---

## Building

1. **Clean previous builds**  
```bash
make clean
```
2. **Build CPU-only**  
```bash
make all
# produces reml, grm.snp, grm.hap, and grm.3order executables
```

---

## Running

Invoke with MPI as usual:

```bash
# build the GRMs for the model you want to use
# replace 4 with the number of MPI tasks you have available
mpirun -np 4 ./bin/grm.snp [options]
mpirun -np 4 ./bin/grm.hap [options]
./bin/grm.snp.3order [options]

# run the prediction using GRMs
mpirun -np 4 ./bin/reml [options]
```

---

## Troubleshooting

- **Header not found**: Ensure MKL include path is in your compiler flags (the variable MKLROOT is set).
- **Library not found**: Ensure `$MKLROOT/lib` is in your linker flags.

---

## Upcoming GPU acceleration

**This section is in development.** GPU acceleration is provided via an override library that intercepts BLAS calls:

1. **Environment**  
```bash
module load cuda          # sets $CUDA_HOME
export NGPUS=1            # GPUs per MPI rank
```
2. **Directory additions**  
```
obj/gpu/                  # GPU object files, including dgemm.override.o
```
3. **Makefile changes**  
- `make reml.gpu` target compiles sources into `obj/gpu/*.o`  
														 - Intercepts `dgemm_` and `cblas_dgemm` via `src/dgemm.override.cpp`  
														 - Links with `-lcublas` from `$CUDA_HOME/lib64`  
4. **Building GPU version**  
```bash
make reml.gpu
# produces bin/reml.gpu
```
5. **Running GPU version**  
```bash
mpirun -np 4 ./bin/reml.gpu [options]
```
