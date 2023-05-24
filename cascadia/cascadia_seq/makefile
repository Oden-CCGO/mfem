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

# Use the MFEM build directory
MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
SRC = $(if $(MFEM_DIR:../..=),$(MFEM_DIR)/cascadia/cascadia_seq/,)
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_INSTALL_DIR = ../../mfem
# CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

MINIAPP = cascadia_seq

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec
.PRECIOUS:

CASCADIA_SRC = wave_map.cpp wave_obs.cpp wave_op.cpp wave_p2o.cpp wave_sol.cpp wave_vis.cpp

CASCADIA_OBJ = $(CASCADIA_SRC:.cpp=.o)

# Remove built-in rules
%: %.cpp
%.o: %.cpp

all: $(MINIAPP)

# Rules for building the miniapps
%: %.o $(CASCADIA_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) $^ -o $@ $(MFEM_LIBS)

%.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) -c $< -o $@

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build #clean-exec

clean-build:
	rm -f *.o *~ $(MINIAPP)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f  sol_p.gf sol_u.gf
	@rm -rf Paraview
	@rm -f  cascadia_seq.mesh cascadia_par.mesh
