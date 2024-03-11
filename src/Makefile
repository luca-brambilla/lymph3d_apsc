SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

COMP=mpif90

COMP_FC_FLAGS=-O5 -g -c -cpp -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -DPETSC_AVOID_MPIF_H
COMP_LD_FLAGS=-O5  -fopenmp

METIS_FC_FLAGS=-I/mnt/c/Users/fulvi/Desktop/Thesis/LYMPH3D_elasticity/fmetis/build/include
METIS_LD_FLAGS=-L/mnt/c/Users/fulvi/Desktop/Thesis/LYMPH3D_elasticity/fmetis/build/lib -lfmetis -lmetis

#METIS_FC_FLAGS=-I/u/gorlezza/LYMPH3D_elasticity/fmetis/build/include
#METIS_LD_FLAGS=-L/u/gorlezza/LYMPH3D_elasticity/fmetis/build/lib -lfmetis -L/u/sw/pkgs/toolchains/gcc-glibc/9/pkgs/metis/5/lib -lmetis

PETSC_FC_FLAGS=-I/mnt/c/Users/fulvi/Desktop/Thesis/petsc/arch-linux-c-debug/include -I/mnt/c/Users/fulvi/Desktop/Thesis/petsc/include
PETSC_LD_FLAGS=-L/mnt/c/Users/fulvi/Desktop/Thesis/petsc/arch-linux-c-debug/lib -Wl,-rpath=/mnt/c/Users/fulvi/Desktop/Thesis/petsc/arch-linux-c-debug/lib -lpetsc

#PETSC_LD_FLAGS=-L/usr/lib/petsc/lib -lpetsc
#PETSC_FC_FLAGS=-I/usr/lib/petsc/include

#PETSC_LD_FLAGS=$(shell pkg-config --libs PETSc)
#PETSC_FC_FLAGS=$(shell pkg-config --cflags PETSc)

EXEC=Lymph3D

.PHONY: all
all: $(OBJS) $(EXEC)

$(EXEC): $(OBJS)
	$(COMP) -o $@  $(OBJS) $(PETSC_LD_FLAGS) $(COMP_LD_FLAGS) $(METIS_LD_FLAGS) 

$(OBJS):%.o: %.f90 assemble_local.o post_processing.o SET_PETSC_SYSTEM.o MOD_VTK.o Poly_ref_mappings.o basis_function.o Poly_setup_MPI.o Poly_mesh.o Poly_data.o Poly_global.o vet_mat_operations.o problem_data_and_properties.o
	$(COMP) $(PETSC_FC_FLAGS) $(COMP_FC_FLAGS) $(METIS_FC_FLAGS) $< -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod Lymph3D
