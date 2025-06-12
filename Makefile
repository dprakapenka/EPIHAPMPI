# ## to set mkl variables: # set mkl env before running make: source /opt/intel/oneapi/setvars.sh
# compiler
CC = icpx
MPICC = mpiicpx
GPU_MPICC = mpiicpx
# compiler flags
CXXFLAGS = -std=c++17 -Wall -DMKL_ILP64 -Iinclude -I${MKLROOT}/include
# for include files
CXXFLAGS += -MMD -MP
# GPU nvidia
# Lets override GPU_ARCH or via env
GPU_ARCH ?= sm_80
GPU_EXTRAS = -fiopenmp -fopenmp-targets=nvptx64-nvidia-cuda --offload-arch=$(GPU_ARCH)

DEBUGFLAGS = -g

OPTFLAGS = -O3 -march=native
LTOFLAG = -flto


# linker flags
# static
MKL_LIBS_STATIC = -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a \
           ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
           ${MKLROOT}/lib/intel64/libmkl_intel_thread.a \
           ${MKLROOT}/lib/intel64/libmkl_core.a \
           ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a \
           -Wl,--end-group
# dynamic
MKL_LIBS =  -L${MKLROOT}/lib \
			-lmkl_scalapack_ilp64 \
			-lmkl_intel_ilp64 \
			-lmkl_intel_thread \
			-lmkl_core \
			-lmkl_blacs_intelmpi_ilp64

LDFLAGS = $(MKL_LIBS) -ilp64 -liomp5 -lpthread -lm -ldl

GPU_LDFLAGS = $(LDFLAGS) -lcublasXt -lcublas -lcudart -ilp64 -liomp5 -lpthread -lm -ldl

# parallel compilation
MAKEFLAGS += -j$(shell nproc)

# directories
SRCDIR = src
INCDIR = include
BINDIR = bin
OBJDIR = obj
OBJDIR_GPU = obj/gpu
TESTSDIR = tests
TESTDIR = test

# output executables
TARGET_REML = $(BINDIR)/reml
TARGET_REML_GPU = $(BINDIR)/reml.gpu
DEBUG_TARGET_REML = $(BINDIR)/reml_debug

TARGET_GRM_SNP = $(BINDIR)/grm.snp
DEBUG_TARGET_GRM_SNP = $(BINDIR)/grm.snp.debug

TARGET_GRM_HAP = $(BINDIR)/grm.hap

TARGET_GRM_SNP_3ORDER = $(BINDIR)/grm.snp.3order

TARGET_READ_BIN      = $(BINDIR)/read.bin

TARGET_RANK = $(BINDIR)/rank

SMALL_TEST_BIN = $(BINDIR)/run_small_matrix_test
LARGE_TEST_BIN = $(BINDIR)/run_large_matrix_test

# SYCL-specific flags
SYCLFLAGS = -fiopenmp -fopenmp-targets=spir64 -fsycl \
            $(MKLROOT)/lib/libmkl_sycl.a \
            -lsycl -lOpenCL

# source and obj files for grm
SRCS_GRM_SNP = $(SRCDIR)/grm.snp.cpp \
           $(SRCDIR)/read.cpp
OBJS_GRM_SNP = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS_GRM_SNP))
SRCS_GRM_HAP = $(SRCDIR)/grm.hap.cpp
OBJS_GRM_HAP = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS_GRM_HAP))

# source and obj files for reml
SRCS_REML = $(SRCDIR)/reml.cpp \
            $(SRCDIR)/gblup.cpp \
            $(SRCDIR)/options.cpp \
            $(SRCDIR)/read.cpp \
            $(SRCDIR)/print.cpp \
            $(SRCDIR)/errors.cpp \
            $(SRCDIR)/functions.cpp
OBJS_REML = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS_REML))

OBJS_REML_GPU = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR_GPU)/%.o, $(SRCS_REML))

SRCS_RANK = $(SRCDIR)/rank.cpp \
            $(SRCDIR)/options.cpp \
            $(SRCDIR)/read.cpp \
            $(SRCDIR)/functions.cpp \
            $(SRCDIR)/errors.cpp
OBJS_RANK = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS_RANK))

###
all: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
all: dirs $(TARGET_REML) $(TARGET_GRM_SNP) $(TARGET_GRM_HAP) $(TARGET_GRM_SNP_3ORDER) $(TARGET_READ_BIN) $(TARGET_RANK)

reml: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
reml: dirs $(TARGET_REML)

reml.gpu: CXXFLAGS += $(GPU_EXTRAS)
reml.gpu: dirs $(TARGET_REML_GPU)

grm: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
grm: dirs $(TARGET_GRM_SNP)

grm.snp: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
grm.snp: dirs $(TARGET_GRM_SNP)

grm.hap: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
grm.hap: dirs $(TARGET_GRM_HAP)

grm.snp.3order: dirs $(TARGET_GRM_SNP_3ORDER)

read.bin: dirs $(TARGET_READ_BIN)

rank: CXXFLAGS += $(OPTFLAGS) $(LTOFLAG)
rank: dirs $(TARGET_RANK)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: dirs $(DEBUG_TARGET_REML) $(DEBUG_TARGET_GRM_SNP) $(DEBUG_TARGET_GRM_HAP)

# Compile the final executable for reml
$(TARGET_REML): $(OBJS_REML)
	$(MPICC) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compile the final executable for reml.gpu
$(TARGET_REML_GPU): $(OBJS_REML_GPU)
	$(GPU_MPICC) $(CXXFLAGS) $^ -o $@ $(GPU_LDFLAGS)

# Compile the final executable for reml.debug
$(DEBUG_TARGET_REML): $(OBJS_REML)
	$(MPICC) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compile grm.snp
$(TARGET_GRM_SNP): $(OBJS_GRM_SNP)
	$(MPICC) $(CXXFLAGS) $(SYCLFLAGS) -o $@ $^ $(LDFLAGS)

# Compile debug version of grm.snp
$(DEBUG_TARGET_GRM_SNP): $(OBJS_GRM_SNP)
	$(MPICC) $(CXXFLAGS) $(DEBUGFLAGS) $(SYCLFLAGS) -o $@ $^ $(LDFLAGS)

# Compile grm.hap
$(TARGET_GRM_HAP): $(OBJS_GRM_HAP)
	$(MPICC) $(CXXFLAGS) $(SYCLFLAGS) -o $@ $^ $(LDFLAGS)

# Compile debug version of grm.hap
$(DEBUG_TARGET_GRM_HAP): $(OBJS_GRM_HAP)
	$(MPICC) $(CXXFLAGS) $(DEBUGFLAGS) $(SYCLFLAGS) -o $@ $^ $(LDFLAGS)

$(TARGET_GRM_SNP_3ORDER): $(SRCDIR)/grm.snp.3order.cpp
	$(CC) -std=c++17 -O3 -xHost -qopenmp $< -o $@

$(TARGET_READ_BIN): $(SRCDIR)/read.bin.cpp | dirs
	$(CC) $< -o $@

# Compile the final executable for rank
$(TARGET_RANK): $(OBJS_RANK)
	$(MPICC) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compile object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | dirs
	$(MPICC) $(CXXFLAGS) -c $< -o $@

# Compile gpu object files
$(OBJDIR_GPU)/%.o: $(SRCDIR)/%.cpp | dirs
	$(MPICC) $(CXXFLAGS) -c $< -o $@

# Ensure necessary directories exist
dirs:
	@mkdir -p $(OBJDIR) $(OBJDIR_GPU) $(BINDIR) $(TESTSDIR) $(TESTDIR)

# Compile small matrix test
$(SMALL_TEST_BIN): tests/test_small_matrix.cpp tests/test_read_binmat_shared.cpp $(OBJDIR)/functions.o | dirs
	$(MPICC) $(CXXFLAGS) -DTESTDIR=\"$(TESTDIR)\" $^ -o $@ $(LDFLAGS)

# Compile large matrix test
$(LARGE_TEST_BIN): tests/test_large_matrix.cpp tests/test_read_binmat_shared.cpp $(OBJDIR)/functions.o | dirs
	$(MPICC) $(CXXFLAGS) -DTESTDIR=\"$(TESTDIR)\" $^ -o $@ $(LDFLAGS)

# generate the manual pdf
manual: manual/manual.pdf

manual/manual.pdf: manual/gen.manual.sh \
				README.md INSTALL.md \
				$(wildcard manual/*.md)
	@echo "Building PDF manual…"
	@bash manual/gen.manual.sh

# Create test files before running tests
setup_tests: dirs
	@echo "Generating test matrices in $(TESTDIR)..."
	$(SMALL_TEST_BIN) --generate-only
	$(LARGE_TEST_BIN) --generate-only
	@echo "Test matrices created in $(TESTDIR)."

# Run tests
TEST_PROCS = 2
run_tests: $(SMALL_TEST_BIN) $(LARGE_TEST_BIN)
	@echo "Running Small Matrix Tests..."
	mpirun -np $(TEST_PROCS) $(SMALL_TEST_BIN)
	@echo "Running Large Matrix Tests..."
	mpirun -np $(TEST_PROCS) $(LARGE_TEST_BIN)
	@echo "All tests completed."

# Clean test binaries but keep test files
clean_tests:
	rm -f $(SMALL_TEST_BIN) $(LARGE_TEST_BIN)
	@echo "Test binaries removed."

# Clean object, binary, and dependency files
clean:
	rm -rf $(OBJDIR) $(OBJDIR_GPU) $(BINDIR) $(OBJS_REML:.o=.d) $(OBJS_REML_GPU:.o=.d) $(OBJS_GRM_SNP:.o=.d) $(OBJS_GRM_HAP:.o=.d)
	@echo "Build files cleaned."

# Clean everything (including test matrices)
clean_all:
	rm -rf $(OBJDIR) $(OBJDIR_GPU) $(BINDIR) $(OBJS_REML:.o=.d) $(OBJS_REML_GPU:.o=.d) $(OBJS_GRM_SNP:.o=.d) $(OBJS_GRM_HAP:.o=.d) $(TESTDIR)
	@echo "Build and test files cleaned."

# Include dependency files
-include $(OBJS_REML:.o=.d) $(OBJS_REML_GPU:.o=.d) $(OBJS_GRM_SNP:.o=.d) $(OBJS_GRM_HAP:.o=.d)

# Help target
help:
	@echo "Available targets:"
	@echo "  all          - Build all executables"
	@echo "  reml         - Build reml executable"
	@echo "  reml.gpu     - Build reml executables with gpu support"
	@echo "  grm          - Build grm snp executable"
	@echo "  grm.snp      - Build grm snp executable"
	@echo "  grm.hap      - Build grm hap executable"
	@echo "  grm.snp3     - Build grm snp third order executable"
	@echo "  debug        - Build debug version of the main target"
	@echo "  clean        - Clean all build files"
	@echo "  clean_all    - Clean all build and test files"
	@echo "  setup_tests  - Generate test matrices"
	@echo "  run_tests    - Run all tests"
	@echo "  manual       - Generate the manual PDF"

.PHONY: all debug dirs clean clean_all setup_tests run_tests manual
