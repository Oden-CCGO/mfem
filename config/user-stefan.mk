# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# This file describes the default MFEM build options.
#
# See the file INSTALL for description of these options.  You can
# customize them below, or copy this file to user.mk and modify it.


# Some choices below are based on the OS type:
NOTMAC := $(subst Darwin,,$(shell uname -s))

ETAGS_BIN = $(shell command -v etags 2> /dev/null)
EGREP_BIN = $(shell command -v egrep 2> /dev/null)

CXX = g++
MPICXX = mpicxx

BASE_FLAGS  = -std=c++11 -Wno-deprecated-declarations
OPTIM_FLAGS = -O3 $(BASE_FLAGS)
DEBUG_FLAGS = -g $(XCOMPILER)-Wall $(BASE_FLAGS)

# Prefixes for passing flags to the compiler and linker when using CXX or MPICXX
CXX_XCOMPILER =
CXX_XLINKER   = -Wl,

# Destination location of make install
# PREFIX = $(HOME)/mfem
PREFIX = ./mfem
# Install program
INSTALL = /usr/bin/install

STATIC = YES
SHARED = NO

# CUDA configuration options
#
# If you set MFEM_USE_ENZYME=YES, CUDA_CXX has to be configured to use cuda with
# clang as its host compiler.
CUDA_CXX = nvcc
CUDA_ARCH = sm_60
CUDA_FLAGS = -x=cu --expt-extended-lambda -arch=$(CUDA_ARCH)
# Prefixes for passing flags to the host compiler and linker when using CUDA_CXX
CUDA_XCOMPILER = -Xcompiler=
CUDA_XLINKER   = -Xlinker=

# HIP configuration options
HIP_CXX = hipcc
# The HIP_ARCH option specifies the AMD GPU processor, similar to CUDA_ARCH. For
# example: gfx600 (tahiti), gfx700 (kaveri), gfx701 (hawaii), gfx801 (carrizo),
# gfx900, gfx1010, etc.
HIP_ARCH = gfx900
HIP_FLAGS = --amdgpu-target=$(HIP_ARCH)
HIP_XCOMPILER =
HIP_XLINKER   = -Wl,

# Flags for generating dependencies.
DEP_FLAGS = -MM -MT

ifneq ($(NOTMAC),)
   AR      = ar
   ARFLAGS = crv
   RANLIB  = ranlib
   PICFLAG = $(XCOMPILER)-fPIC
   SO_EXT  = so
   SO_VER  = so.$(MFEM_VERSION_STRING)
   BUILD_SOFLAGS = -shared $(XLINKER)-soname,libmfem.$(SO_VER)
   BUILD_RPATH = $(XLINKER)-rpath,$(BUILD_REAL_DIR)
   INSTALL_SOFLAGS = $(BUILD_SOFLAGS)
   INSTALL_RPATH = $(XLINKER)-rpath,@MFEM_LIB_DIR@
else
   # Silence "has no symbols" warnings on Mac OS X
   AR      = ar
   ARFLAGS = Scrv
   RANLIB  = ranlib -no_warning_for_no_symbols
   PICFLAG = $(XCOMPILER)-fPIC
   SO_EXT  = dylib
   SO_VER  = $(MFEM_VERSION_STRING).dylib
   MAKE_SOFLAGS = $(XLINKER)-dylib,-install_name,$(1)/libmfem.$(SO_VER),\
      -compatibility_version,$(MFEM_VERSION_STRING),\
      -current_version,$(MFEM_VERSION_STRING),\
      -undefined,dynamic_lookup
   BUILD_SOFLAGS = $(subst $1 ,,$(call MAKE_SOFLAGS,$(BUILD_REAL_DIR)))
   BUILD_RPATH = $(XLINKER)-undefined,dynamic_lookup
   INSTALL_SOFLAGS = $(subst $1 ,,$(call MAKE_SOFLAGS,$(MFEM_LIB_DIR)))
   INSTALL_RPATH = $(XLINKER)-undefined,dynamic_lookup
   # Silence unused command line argument warnings when generating dependencies
   # with mpicxx and clang
   DEP_FLAGS := -Wno-unused-command-line-argument $(DEP_FLAGS)
endif

# Set CXXFLAGS to overwrite the default selection of DEBUG_FLAGS/OPTIM_FLAGS
# CXXFLAGS = -O3 -march=native

# Optional extra compile flags, in addition to CXXFLAGS:
# CPPFLAGS =

# Library configurations:
# Note: symbols of the form @VAR@ will be replaced by $(VAR) in derived
#       variables, like MFEM_FLAGS, defined in config.mk.

# Command used to launch MPI jobs
MFEM_MPIEXEC    = mpirun
MFEM_MPIEXEC_NP = -np
# Number of mpi tasks for parallel jobs
MFEM_MPI_NP = 4

# MFEM configuration options: YES/NO values, which are exported to config.mk and
# config.hpp. The values below are the defaults for generating the actual values
# in config.mk and config.hpp.

MFEM_USE_MPI           = YES
MFEM_USE_METIS         = $(MFEM_USE_MPI)
MFEM_USE_METIS_5       = YES
MFEM_DEBUG             = NO
MFEM_USE_EXCEPTIONS    = NO
MFEM_USE_ZLIB          = NO
MFEM_USE_LIBUNWIND     = NO
MFEM_USE_LAPACK        = YES
MFEM_THREAD_SAFE       = NO
MFEM_USE_OPENMP        = NO
MFEM_USE_LEGACY_OPENMP = NO
MFEM_USE_MEMALLOC      = YES
MFEM_TIMER_TYPE        = $(if $(NOTMAC),2,4)
MFEM_USE_SUNDIALS      = NO
MFEM_USE_MESQUITE      = NO
MFEM_USE_SUITESPARSE   = NO
MFEM_USE_SUPERLU       = YES
MFEM_USE_SUPERLU5      = NO
MFEM_USE_MUMPS         = YES
MFEM_USE_STRUMPACK     = NO
MFEM_USE_GINKGO        = NO
MFEM_USE_AMGX          = NO
MFEM_USE_GNUTLS        = NO
MFEM_USE_NETCDF        = NO
MFEM_USE_PETSC         = YES
MFEM_USE_SLEPC         = NO
MFEM_USE_MPFR          = NO
MFEM_USE_SIDRE         = NO
MFEM_USE_FMS           = NO
MFEM_USE_CONDUIT       = NO
MFEM_USE_PUMI          = NO
MFEM_USE_HIOP          = NO
MFEM_USE_GSLIB         = NO
MFEM_USE_CUDA          = NO
MFEM_USE_HIP           = NO
MFEM_USE_RAJA          = NO
MFEM_USE_OCCA          = NO
MFEM_USE_CEED          = NO
MFEM_USE_CALIPER       = NO
MFEM_USE_ALGOIM        = NO
MFEM_USE_UMPIRE        = NO
MFEM_USE_SIMD          = NO
MFEM_USE_ADIOS2        = NO
MFEM_USE_MKL_CPARDISO  = NO
MFEM_USE_MOONOLITH     = NO
MFEM_USE_ADFORWARD     = NO
MFEM_USE_CODIPACK      = NO
MFEM_USE_BENCHMARK     = NO
MFEM_USE_PARELAG       = NO
MFEM_USE_ENZYME        = NO

# MPI library compile and link flags
# These settings are used only when building MFEM with MPI + HIP
ifeq ($(MFEM_USE_MPI)$(MFEM_USE_HIP),YESYES)
   # We determine MPI_DIR assuming $(MPICXX) is in $(MPI_DIR)/bin
   MPI_DIR := $(patsubst %/,%,$(dir $(shell which $(MPICXX))))
   MPI_DIR := $(patsubst %/,%,$(dir $(MPI_DIR)))
   MPI_OPT = -I$(MPI_DIR)/include
   MPI_LIB = -L$(MPI_DIR)/lib $(XLINKER)-rpath,$(MPI_DIR)/lib -lmpi
endif

# ROCM/HIP directory such that ROCM/HIP libraries like rocsparse and rocrand are
# found in $(HIP_DIR)/lib, usually as links. Typically, this directory is of
# the form /opt/rocm-X.Y.Z which is called ROCM_PATH by hipconfig.
ifeq ($(MFEM_USE_HIP),YES)
   HIP_DIR := $(patsubst %/,%,$(dir $(shell which $(HIP_CXX))))
   HIP_DIR := $(patsubst %/,%,$(dir $(HIP_DIR)))
   ifeq (,$(wildcard $(HIP_DIR)/lib/librocsparse.*))
      HIP_DIR := $(shell hipconfig --rocmpath 2> /dev/null)
      ifeq (,$(wildcard $(HIP_DIR)/lib/librocsparse.*))
         $(error Unable to determine HIP_DIR. Please set it manually.)
      endif
   endif
endif

# Compile and link options for zlib.
ZLIB_DIR =
ZLIB_OPT = $(if $(ZLIB_DIR),-I$(ZLIB_DIR)/include)
ZLIB_LIB = $(if $(ZLIB_DIR),$(ZLIB_RPATH) -L$(ZLIB_DIR)/lib ,)-lz
ZLIB_RPATH = $(XLINKER)-rpath,$(ZLIB_DIR)/lib

LIBUNWIND_OPT = -g
LIBUNWIND_LIB = $(if $(NOTMAC),-lunwind -ldl,)

PETSC_ARCH := arch-osx-real
PETSC_DIR  := /Users/stefan/CCGO/LIBS/MFEM/petsc/$(PETSC_ARCH)

# HYPRE library configuration (needed to build the parallel version)
HYPRE_DIR = /Users/stefan/opt/hypre-2.25.0/src/hypre
HYPRE_OPT = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE
ifeq (YES,$(MFEM_USE_CUDA))
   # This is only necessary when hypre is built with cuda:
   HYPRE_LIB += -lcusparse -lcurand -lcublas
endif
ifeq (YES,$(MFEM_USE_HIP))
   # This is only necessary when hypre is built with hip:
   HYPRE_LIB += -L$(HIP_DIR)/lib $(XLINKER)-rpath,$(HIP_DIR)/lib\
 -lrocsparse -lrocrand
endif

# METIS library configuration
ifeq ($(MFEM_USE_SUPERLU)$(MFEM_USE_STRUMPACK)$(MFEM_USE_MUMPS),NONONO)
   ifeq ($(MFEM_USE_METIS_5),NO)
     METIS_DIR = @MFEM_DIR@/../metis-4.0
     METIS_OPT =
     METIS_LIB = -L$(METIS_DIR) -lmetis
   else
#		METIS_DIR = @MFEM_DIR@/../metis-5.0
#		METIS_OPT = -I$(METIS_DIR)/include
#		METIS_LIB = -L$(METIS_DIR)/lib -lmetis
		METIS_OPT = -I/Users/stefan/opt/parmetis-4.0.3/metis/include
		METIS_LIB = -L/Users/stefan/opt/parmetis-4.0.3/build/Darwin-x86_64/libmetis\
			-lmetis
   endif
else
   # ParMETIS: currently needed by SuperLU or STRUMPACK. We assume that METIS 5
   # (included with ParMETIS) is installed in the same location.
   # Starting with STRUMPACK v2.2.0, ParMETIS is an optional dependency while
   # METIS is still required.
#   METIS_DIR = @MFEM_DIR@/../parmetis-4.0.3
#   METIS_OPT = -I$(METIS_DIR)/include
#   METIS_LIB = -L$(METIS_DIR)/lib -lparmetis -lmetis
	METIS_OPT = -I/Users/stefan/opt/parmetis-4.0.3/include
	METIS_LIB = -L/Users/stefan/opt/parmetis-4.0.3/build/Darwin-x86_64/libparmetis\
		-lparmetis
   MFEM_USE_METIS_5 = YES
endif

# LAPACK library configuration
LAPACK_OPT =
LAPACK_LIB = $(if $(NOTMAC),-llapack -lblas,-framework Accelerate)

# OpenMP configuration
OPENMP_OPT = $(XCOMPILER)-fopenmp
OPENMP_LIB =

# Used when MFEM_TIMER_TYPE = 2
POSIX_CLOCKS_LIB = -lrt

# SuperLU library configuration
ifeq ($(MFEM_USE_SUPERLU5),YES)
   SUPERLU_DIR = @MFEM_DIR@/../SuperLU_DIST_5.1.0
   SUPERLU_OPT = -I$(SUPERLU_DIR)/include
   SUPERLU_LIB = $(XLINKER)-rpath,$(SUPERLU_DIR)/lib -L$(SUPERLU_DIR)/lib\
      -lsuperlu_dist_5.1.0
else
#   SUPERLU_DIR = @MFEM_DIR@/../SuperLU_DIST_6.3.1
#   SUPERLU_OPT = -I$(SUPERLU_DIR)/include
#   SUPERLU_LIB = $(XLINKER)-rpath,$(SUPERLU_DIR)/lib64 -L$(SUPERLU_DIR)/lib64\
#      -lsuperlu_dist -lblas
	SUPERLU_DIR = $(PETSC_DIR)
   SUPERLU_OPT = -I$(SUPERLU_DIR)/include
   SUPERLU_LIB = $(XLINKER)-rpath,$(SUPERLU_DIR)/lib -L$(SUPERLU_DIR)/lib\
      -lsuperlu_dist
endif

# SCOTCH library configuration (required by STRUMPACK <= v2.1.0, optional in
# STRUMPACK >= v2.2.0)
SCOTCH_DIR = @MFEM_DIR@/../scotch_6.0.4
SCOTCH_OPT = -I$(SCOTCH_DIR)/include
SCOTCH_LIB = -L$(SCOTCH_DIR)/lib -lptscotch -lptscotcherr -lscotch -lscotcherr\
 -lpthread

# SCALAPACK library configuration (required by STRUMPACK and MUMPS)
SCALAPACK_DIR = /usr/local
SCALAPACK_OPT = -I$(SCALAPACK_DIR)/include
SCALAPACK_LIB = -L$(SCALAPACK_DIR)/lib -lscalapack $(LAPACK_LIB)

# MPI Fortran library, needed e.g. by STRUMPACK or MUMPS
# MPICH:
# MPI_FORTRAN_LIB = -lmpifort
# OpenMPI:
MPI_FORTRAN_LIB = -lmpi_mpifh
# Additional Fortran library:
# MPI_FORTRAN_LIB += -lgfortran

# MUMPS library configuration
#MUMPS_DIR = @MFEM_DIR@/../MUMPS_5.2.0
MUMPS_DIR = $(PETSC_DIR)
MUMPS_OPT = -I$(MUMPS_DIR)/include
MUMPS_LIB = $(XLINKER)-rpath,$(MUMPS_DIR)/lib -L$(MUMPS_DIR)/lib -ldmumps\
 -lmumps_common -lpord $(SCALAPACK_LIB) $(LAPACK_LIB) $(MPI_FORTRAN_LIB)

# NetCDF library configuration
NETCDF_DIR = $(HOME)/local
HDF5_DIR   = $(HOME)/local
NETCDF_OPT = -I$(NETCDF_DIR)/include -I$(HDF5_DIR)/include $(ZLIB_OPT)
NETCDF_LIB = $(XLINKER)-rpath,$(NETCDF_DIR)/lib -L$(NETCDF_DIR)/lib\
 $(XLINKER)-rpath,$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib\
 -lnetcdf -lhdf5_hl -lhdf5 $(ZLIB_LIB)

# PETSc library configuration (version greater or equal to 3.8 or the dev branch)
PETSC_VARS := $(PETSC_DIR)/lib/petsc/conf/petscvariables
PETSC_FOUND := $(if $(wildcard $(PETSC_VARS)),YES,)
PETSC_INC_VAR = PETSC_CC_INCLUDES
PETSC_LIB_VAR = PETSC_EXTERNAL_LIB_BASIC
ifeq ($(PETSC_FOUND),YES)
   PETSC_OPT := $(shell sed -n "s/$(PETSC_INC_VAR) = *//p" $(PETSC_VARS))
   PETSC_DEP := $(shell sed -n "s/$(PETSC_LIB_VAR) = *//p" $(PETSC_VARS))
   PETSC_LIB = $(XLINKER)-rpath,$(abspath $(PETSC_DIR))/lib\
      -L$(abspath $(PETSC_DIR))/lib -lpetsc\
      $(subst $(CXX_XLINKER),$(XLINKER),$(PETSC_DEP))
endif

# MPFR library configuration
MPFR_OPT =
MPFR_LIB = -lmpfr

# If YES, enable some informational messages
VERBOSE = NO

# Optional build tag
MFEM_BUILD_TAG = $(shell uname -snm)
