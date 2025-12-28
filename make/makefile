# Define the directories for the object and module files,
# the executable, and the executable name and path.
OBJDIR = ./obj_temp
BINDIR = .
PROGRAM = ../msg123

SOURCEDIR = ../src

VPATH = \
${SOURCEDIR}

.SUFFIXES: .f90 .o

# Define the Fortran compile flags (gfortran or ifort)
FC = ifort

# Define the mode (release or debug)
MODE = release

# Define the parallel (MPI, OpenMP, Hybrid)
PARA = OpenMP

#for gfortran
ifeq ($(FC),gfortran)
	ifeq ($(MODE),debug)
		FFLAGS = -cpp -Wall -g -fbacktrace -fbounds-check  -Wuninitialized -pedantic -ffpe-trap=invalid,zero,overflow
	else
		FFLAGS = -cpp -Wall -O3 -fmax-stack-var-size=65535
	endif
	ifeq ($(PARA),Hybrid)
		FC = mpif90
		FFLAGS = -cpp -Wall -O3 -fopenmp
		CPPFLAGS = -DMPI_MSG
	endif
	ifeq ($(PARA),MPI)
		FC = mpif90
		CPPFLAGS = -DMPI_MSG
	endif
	ifeq ($(PARA),OpenMP)
		FFLAGS = -cpp -Wall -O3 -fopenmp
	endif
endif

#for ifortran
ifeq ($(FC),ifort)
	ifeq ($(MODE),debug)
		FFLAGS = -fpp -warn all -g -traceback -check bounds -fpe0 -check all -std -ftrapuv -gen_interfaces
	else
		FFLAGS = -fpp -warn all -O3 -heap-arrays
	endif
	ifeq ($(PARA),Hybrid)
		FC = mpif90
		FFLAGS += -qopenmp
		CPPFLAGS = -DMPI_MSG
	endif
	ifeq ($(PARA),MPI)
		FC = mpif90
		CPPFLAGS = -DMPI_MSG
	endif
	ifeq ($(PARA),OpenMP)
		FFLAGS += -qopenmp
	endif
endif

# Define the C compile flags
CC = gcc
CFLAGS = -O3 -D_UF

ifeq ($(CPPFLAGS),-DMPI_MSG)
OBJECTS = \
$(OBJDIR)/Kind_module.o \
$(OBJDIR)/Constval_module.o \
$(OBJDIR)/Mpi_initfin.o \
$(OBJDIR)/Utility_module.o \
$(OBJDIR)/Mpi_utility.o \
$(OBJDIR)/Initial_module.o \
$(OBJDIR)/Mpi_read.o \
$(OBJDIR)/Read_module.o \
$(OBJDIR)/Open_file.o \
$(OBJDIR)/Read_input.o \
$(OBJDIR)/Mpi_set.o \
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
$(OBJDIR)/Mpi_solve.o \
$(OBJDIR)/Calc_function.o \
$(OBJDIR)/Make_amg_matrix.o \
$(OBJDIR)/Make_linearsystem.o \
$(OBJDIR)/Check_simulation.o \
$(OBJDIR)/Calc_simulation.o \
$(OBJDIR)/Linear_solution.o \
$(OBJDIR)/Allocate_output.o \
$(OBJDIR)/Mpi_write.o \
$(OBJDIR)/Calc_output.o \
$(OBJDIR)/Write_module.o \
$(OBJDIR)/Time_module.o \
$(OBJDIR)/Nonlinear_solution.o \
$(OBJDIR)/Write_output.o \
$(OBJDIR)/Msg.o
else
OBJECTS = \
$(OBJDIR)/Kind_module.o \
$(OBJDIR)/Constval_module.o \
$(OBJDIR)/Utility_module.o \
$(OBJDIR)/Initial_module.o \
$(OBJDIR)/Read_module.o \
$(OBJDIR)/Open_file.o \
$(OBJDIR)/Read_input.o \
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
$(OBJDIR)/Calc_function.o \
$(OBJDIR)/Make_amg_matrix.o \
$(OBJDIR)/Make_linearsystem.o \
$(OBJDIR)/Check_simulation.o \
$(OBJDIR)/Calc_simulation.o \
$(OBJDIR)/Linear_solution.o \
$(OBJDIR)/Allocate_output.o \
$(OBJDIR)/Calc_output.o \
$(OBJDIR)/Write_module.o \
$(OBJDIR)/Time_module.o \
$(OBJDIR)/Nonlinear_solution.o \
$(OBJDIR)/Write_output.o \
$(OBJDIR)/Msg.o
endif

# Define task functions

# Create the bin directory and compile and link the program
all: makebin | $(PROGRAM)

# Make the bin directory for the executable
makebin :
	mkdir -p $(BINDIR)

# Define the objects that make up the program
$(PROGRAM) : $(OBJECTS)
	-$(FC) $(FFLAGS) -o $@ $(OBJECTS)

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

# Clean the object and module files and the executable
.PHONY : clean
clean :
	-rm -rf $(OBJDIR)
	-rm -rf $(PROGRAM)
	-rm -rf *.mod

# Clean the object and module files
.PHONY : cleanobj
cleanobj :
	-rm -rf $(OBJDIR)
