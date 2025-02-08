SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

COMP=mpif90

DIR_FMETIS=/home/user/lymph3d_apsc/fmetis/build/
DIR_METIS=/home/user/lymph3d_apsc/metis-5.1.0/build/Linux-x86_64

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

$(OBJS):%.o: %.f90
	$(COMP) $(METIS_FC_FLAGS) $(DEBUG_FC_FLAGS) -c $< -o $@

# ubuntu package dependencies created with makedepf90
include .depend

depend: $(SRCS)
	makedepf90 *.f90 > .depend

.PHONY: clean
clean:
	-rm -f *.o *.mod Lymph3D
