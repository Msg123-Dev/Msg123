# Define the directories for the object and module files,
# the executable, and the executable name and path.
OBJDIR = ./obj_temp
PROGRAM = ../msg123
SOURCEDIR = ../src

VPATH = ${SOURCEDIR}

.SUFFIXES: .f90 .o

SHELL := /bin/bash

# Compiler: gfortran / ifort / ifx
FC = gfortran

# Mode: release or debug
MODE = release

# Parallelism: Serial / OpenMP / MPI / Hybrid
PARA = Serial

# ---- MPI flag ----
_FC_ORIG := $(FC)
ifneq ($(filter MPI Hybrid,$(PARA)),)
    USE_MPI := 1
    FC = mpif90
    CPPFLAGS = -DMPI_MSG
endif

# ---- gfortran ----
ifeq ($(_FC_ORIG),gfortran)
	ifeq ($(MODE),debug)
		FFLAGS = -cpp -Wall -g -fbacktrace -fbounds-check  -Wuninitialized -pedantic -ffpe-trap=invalid,zero,overflow
	else
		FFLAGS = -cpp -Wall -O3 -Wno-unused-function
	endif
	ifneq ($(filter OpenMP Hybrid,$(PARA)),)
		FFLAGS += -fopenmp
	endif
endif

# ---- ifort / ifx ----
ifneq ($(filter $(_FC_ORIG),ifort ifx),)
	ifeq ($(MODE),debug)
		FFLAGS = -fpp -warn all -g -traceback -check bounds -fpe0 -check all -std -ftrapuv -gen_interfaces
	else
		FFLAGS = -fpp -warn all -O3
	endif
	FFLAGS += -diag-disable=10448 -fpscomp logicals
	ifneq ($(filter OpenMP Hybrid,$(PARA)),)
		FFLAGS += -qopenmp
	endif
	ifneq ($(filter MPI OpenMP Hybrid,$(PARA)),)
		FFLAGS += -diag-disable=8889
	endif
endif

# ---- C flags ----
CC = gcc
CFLAGS = -O3 -D_UF

# ---- Object list ----
# MPI-specific objects resolve to empty string when USE_MPI is not set,
# preserving compilation order without duplicating the full list.
ifeq ($(USE_MPI),1)
    MPI_INITFIN = $(OBJDIR)/Mpi_initfin.o
    MPI_UTILITY = $(OBJDIR)/Mpi_utility.o
    MPI_READ = $(OBJDIR)/Mpi_read.o
    MPI_SET = $(OBJDIR)/Mpi_set.o
    MPI_SOLVE = $(OBJDIR)/Mpi_solve.o
    MPI_WRITE = $(OBJDIR)/Mpi_write.o
endif

OBJECTS = \
$(OBJDIR)/Kind_module.o \
$(OBJDIR)/Constval_module.o \
$(MPI_INITFIN) \
$(OBJDIR)/Utility_module.o \
$(MPI_UTILITY) \
$(OBJDIR)/Initial_module.o \
$(MPI_READ) \
$(OBJDIR)/Read_module.o \
$(OBJDIR)/Open_file.o \
$(OBJDIR)/Read_input.o \
$(MPI_SET) \
$(OBJDIR)/Check_condition.o \
$(OBJDIR)/Set_cell.o \
$(OBJDIR)/Make_cell.o \
$(OBJDIR)/Set_condition.o \
$(OBJDIR)/Assign_calculation.o \
$(OBJDIR)/Prep_calculation.o \
$(OBJDIR)/Assign_boundary.o \
$(OBJDIR)/Calc_boundary.o \
$(OBJDIR)/Set_boundary.o \
$(OBJDIR)/Calc_parameter.o \
$(OBJDIR)/Allocate_solution.o \
$(MPI_SOLVE) \
$(OBJDIR)/Calc_function.o \
$(OBJDIR)/Make_amg_matrix.o \
$(OBJDIR)/Make_linearsystem.o \
$(OBJDIR)/Check_simulation.o \
$(OBJDIR)/Calc_simulation.o \
$(OBJDIR)/Linear_solution.o \
$(OBJDIR)/Allocate_output.o \
$(MPI_WRITE) \
$(OBJDIR)/Calc_output.o \
$(OBJDIR)/Write_module.o \
$(OBJDIR)/Time_module.o \
$(OBJDIR)/Nonlinear_solution.o \
$(OBJDIR)/Write_output.o \
$(OBJDIR)/Msg.o

# ---- Build rules ----
.PHONY: all clean cleanobj

all: $(PROGRAM)

$(PROGRAM) : $(OBJECTS)
	@echo "$(FC) $(FFLAGS) -o $@ $(OBJECTS)"
	@$(FC) $(FFLAGS) -o $@ $(OBJECTS)

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	@echo "$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@"
	@$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

# Clean the object and module files and the executable
clean :
	-rm -rf $(OBJDIR)
	-rm -rf $(PROGRAM)
	-rm -rf *.mod

# Clean the object and module files
cleanobj :
	-rm -rf $(OBJDIR)
