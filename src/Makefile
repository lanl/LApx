# to add a
# Compiler/Linker settings
FC = mpif90
INCLUDEDIRS= -I$(FFTW_DIR)/include  -Imacro -I$(HDF5_DIR)/include -I$(F2018_ENV_DIR)/include
INCLUDELIBS= -L$(FFTW_DIR)/lib  -L$(HDF5_DIR)/lib  -L$(F2018_ENV_DIR)/lib -lfftw3 -lfftw3_mpi



## Compiler type
FOUND_COMPILER=0
IFORT=0
GFORTRAN=0
VARIABLE:= $(shell mpif90 --version | grep IFORT)
$(info mpif90 --version | grep IFORT: $(VARIABLE))
ifneq (,$(findstring IFORT,$(VARIABLE)))
 $(info INTEL compiler detected)
 IFORT=1
 FOUND_COMPILER=1
endif

ifneq ($(FOUND_COMPILER),1)
 VARIABLE:= $(shell mpif90 --version | grep "GNU Fortran")
 $(info mpif90 --version | grep gfortran: $(VARIABLE))
 ifneq (,$(findstring GNU Fortran,$(VARIABLE)))
  $(info gFortran compiler detected)
  GFORTRAN=1
  FOUND_COMPILER=1
 endif
endif

## os dependent flags
ifeq ($(OSTYPE),linux-gnu)
 # LINUX
 ifeq (IFORT,1)
  H5BLD_LDFLAGS= -free -warn all 
 endif
 ifeq (GFORTRAN,1)
  H5BLD_LDFLAGS= -Wl,-headerpad_max_install_names -Wl,-rpath,$(HDF5_DIR)/lib
 endif 
elseifeq ($(OSTYPE),darwin))
 # Mac OSX
 H5BLD_LDFLAGS= -Wl,-pie -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,$(HDF5_DIR)/lib
endif


ifeq ($(IFORT),1) 
# -stand f18                      # issue compile-time messages for nonstandard language elements.
# -O3                             # set the optimizations level
# -unroll                         # [=n] set the maximum number of times to unroll loops (no number n means automatic).
# -unroll-aggressive              # use more aggressive unrolling for certain loops.
# -diag-disable=10346             # optimization reporting will be enabled at link time when performing interprocedural optimizations.
# -diag-disable=10397             # optimization reporting will be enabled at link time when performing interprocedural optimizations.
 #-guide-vec=4                    # enable guidance for auto-vectorization, causing the compiler to generate messages suggesting ways to improve optimization (default=4, highest).
 #-parallel                       # generate multithreaded code for loops that can be safely executed in parallel. This option requires MKL libraries.
 #-qopt-subscript-in-range        # assumes there are no "large" integers being used or being computed inside loops. A "large" integer is typically > 2^31.
# -ftz                            # Flushes subnormal results to zero.
# -inline-forceinline # Instructs the compiler to force inlining of functions suggested for inlining whenever the compiler is capable doing so.
# -finline-functions  # enables function inlining for single file compilation.
# -ipo                # enable interprocedural optimization between files.
# -ip                 # determines whether additional interprocedural optimizations for single-file compilation are enabled.
FSTD = -stand f18
 PP = -fpp
 TRACEBACK = -traceback
 FFLAGS:= -D__INTEL__ -m64 -convert big_endian -extend-source -diag-disable=5268 -unroll -unroll-aggressive -diag-disable=10346 -diag-disable=10397 -ftz -inline-forceinline -finline-functions -ipo -ip 
endif

ifeq ($(GFORTRAN),1)
 FSTD = -std=f2018
 PP = -cpp
 TRACEBACK = -ftraceback
 FFLAGS:= -D__GFORTRAN__ -m64 -fconvert=big-endian -fimplicit-none -ffree-line-length-none -fall-intrinsics
endif
FFLAGS:= $(FSTD) $(PP) $(FFLAGS)
$(info FLAGS = $(FFLAGS))

H5BLD_LIBS= -lz -ldl -lm -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran
#FFLAGS = -m64 -fpp -fall-intrinsics -fno-align-commons -fimplicit-none  -ffree-line-length-none -fconvert=big-endian -fbacktrace
BLAS = -llapack
BASE_PROGRAM_NAME = LApx
RM = rm -f
# the USE calls we want to skip, seprate multiple arguments with a space
EXCLUDE_USE="omp_lib mpi mpi_f08 hdf5"

## OPTIONS:
## OPT: the level of optimization: 1, 2, or 3, default 2. Only useful for METHOD=opt
ifndef OPT
 OPT=2
endif

## METHOD: the type of compilation. Valid options are: opt and dbg. opt is default
ifndef METHOD
 METHOD=opt
endif

# define debug flags
ifeq ($(METHOD),dbg)
 ifeq ($(IFORT),1)
 #-debug full           # generate full debug information
 #-g3                   # generate full debug information
 #-O0                   # disable optimizations
 #-CB                   # Perform run-time bound-checks on array subscript and substring references (same as the -check bounds option)
 #-init:snan,arrays     # initialize arrays and scalars to NaN
 #-warn all             # enable all warning
 #-gen-interfaces       # generate interface block for each routine in the source files
 #-traceback            # trace back for debugging
 #-check all            # check all
 #-check bounds         # check array bounds
 #-fpe-all=0            # Floating-point invalid, divide-by-zero, and overflow exceptions are enabled
 #-fpe0                 # Ignore underflow (yield 0.0); Abort on other IEEE exceptions.
 #-diag-error-limit=10  # max diagnostic error count
 #-diag-disable=5268    # Extension to standard: The text exceeds right hand column allowed on the line.
 #-diag-disable=7025    # This directive is not standard Fxx.
 #-diag-disable=10346   # optimization reporting will be enabled at link time when performing interprocedural optimizations.
 #-ftrapuv              # Initializes stack local variables to an unusual value to aid error detection.
   DEBUGFLAGS = -D__DEBUG__  -traceback -check bounds -check uninit -fpe0 -warn all -diag-disable=8889 -diag-disable=5268 -diag-disable=10346 -diag-disable=702 -diag-disable=406
  endif
  ifeq ($(GFORTRAN),1)
   DEBUGFLAGS = -D__DEBUG__ -gdwarf-2 -Wall -fbounds-check -pedantic-errors -ggdb -O -Wall  -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
  endif
 OPT=0
endif

## OMP: valid options are: 0 or 1. ) is default
ifeq ($(OMP),1)
 OMPFLAGS = -fopenmp
endif

##
## EXAMPLE: make OMP=1 METHOD=dbg

##
## TARGETS
## make: build LApx with the provided options, see above
default:
 # Find all source files in all directories and subdorectories and resolve program dependencies
		$(eval ORDERED_SRCS:= $(shell python3 fortran_dependency_resolution.py --exclude-use=$(EXCLUDE_USE)))

 # actually compile
		$(FC)  $(FSTD) $(DEBUGFLAGS) -O$(OPT) $(INCLUDEDIRS) $(INCLUDELIBS) $(FFLAGS) $(BLAS) $(H5BLD_LDFLAGS) $(ORDERED_SRCS) $(H5BLD_LIBS) -o $(BASE_PROGRAM_NAME)-$(METHOD)
  # delete modules
		@$(RM) *.mod
		@$(RM) *.smod
		@$(RM) *.pyc


## make clean: remove executables and other files
clean:
		@rm -f LApx-opt
		@rm -f LApx-dbg
		@rm -rf LApx-dbg.*
		@rm -rf *.mod
		@rm -rf *.smod
		@rm -rf *.pyc
## make help: prints the help
.PHONY: help
help:
		@awk '/^##/ {$$1=""; print $$0}' Makefile
