SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

COMP=mpif90

DIR_FMETIS=/home/luca/fmetis/build
DIR_METIS=/home/luca/metis-5.1.0/build/Linux-x86_64

COMP_FC_FLAGS=-O5 -g -c -cpp -ffree-form -ffree-line-length-none -fopenmp -fbounds-check
COMP_LD_FLAGS=-O5 -fopenmp

DEBUG_FC_FLAGS = -Wunused -Wconversion -fcheck=all -fcheck-array-temporaries  -ffree-line-length-none -finit-real=snan -ffpe-trap=invalid -fbacktrace -ffpe-trap=underflow -Wuninitialized #-Wall -fcheck=all -fbacktrace

METIS_FC_FLAGS=-I$(DIR_FMETIS)/include -I$(DIR_METIS)/include
METIS_LD_FLAGS=-L$(DIR_FMETIS)/lib -lfmetis -L$(DIR_METIS)/libmetis -lmetis

EXEC=Lymph3D

.PHONY: all
all: $(OBJS) $(EXEC)

$(EXEC): $(OBJS)
	$(COMP) -o $@  $(OBJS) $(COMP_LD_FLAGS) $(METIS_LD_FLAGS)
#rm -f *.o *.mod

$(OBJS):%.o: %.f90
	$(COMP) $(METIS_FC_FLAGS) $(DEBUG_FC_FLAGS) -c $< -o $@

# ubuntu package dependencies created with makedepf90
include .depend

# Poly_setup_MPI:%.o: %.f90
# matrix_free.o: matrix_free.f90 Poly_setup_MPI.o problem_data_and_properties.o
#Lymph3D:%.o: %.f90 Poly_data.o Poly_setup_MPI.o problem_data_and_properties.o Poly_mesh.o post_processing.o SET_PETSC_SYSTEM.o matrix_free.o global_parameters.o

# $(OBJS):%.o: %.f90 assemble_element.o post_processing.o SET_PETSC_SYSTEM.o MOD_VTK.o Poly_ref_mappings.o basis_function.o Poly_setup_MPI.o Poly_mesh.o Poly_data.o Poly_global.o vet_mat_operations.o problem_data_and_properties.o MOD_MPI.o matrix_free.o global_parameters.o
# 	$(COMP) $(PETSC_FC_FLAGS) $(COMP_FC_FLAGS) $(METIS_FC_FLAGS) $(DEBUG_FC_FLAGS) $< -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod Lymph3D

depend: $(SRCS)
	makedepf90 *.f90 > .depend