AMREX_HOME ?= ../

EBASE := Lithium

PRECISION  = DOUBLE
PROFILE    = FALSE

# DEBUG      = TRUE
DEBUG      = FALSE

# DIM       = 2
DIM       = 3
USE_HYPRE = FALSE
USE_PETSC = FALSE

COMP	   = gnu

USE_MPI    = TRUE
# USE_OMP    = TRUE
ifeq ($(USE_HYPRE),TRUE)
	HYPRE_DIR ?= $(AMREX_HOME)/Src/LinearSolvers/hypre
endif
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bpack   := ./Make.package 
Blocs   := . 

INCLUDE_LOCATIONS += .
VPATH_LOCATIONS += .

include $(Bpack)


INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs 	:= Base Boundary AmrCore Particle LinearSolvers/C_CellMG LinearSolvers/MLMG
ifeq ($(USE_HYPRE),TRUE)
    Pdirs += Extern/HYPRE
endif

ifeq ($(USE_PETSC),TRUE)
    Pdirs += Extern/PETSc
endif
Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

all: $(executable) 
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
