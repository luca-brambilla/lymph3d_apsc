# LYMPH3D Library

LYMPH3D is a library written in Fortran which applies the PolyDG method to solve linear elastodynamic problem under both Dirichlet and Neumann boundary conditions.
The library numerically computes volume and surface integrals using Gaussian quadrature formulas, in particular by means of Gauss-Legendre nodes and weights. Given that the computational mesh consists of an agglomeration of tetrahedral and triangular elements, the implementation necessitates appropriate mappings from the reference square to the reference triangle and from the reference cube to the reference tetrahedron.

The solution is computed in a matrix-free paradigm.

## Installation

### Prerequisites

- Fortran compiler `gfortran`
- `make` and `cmake`
- `METIS` (for mesh agglomeration)
- `fMETIS` (METIS Fortran interface)
- `MPI`   (parallelization)
- `makedepf90` (automatically update dependencies in `src` directory)
- `doxygen` (for documentation generation)
- `Paraview` (for solution visualization)

Install dependencies:
```
apt-get update
apt-get install -y python3 make cmake g++ gfortran openmpi-bin openmpi-common libopenmpi-dev makedepf90
```

For a fast installation of METIS and fMETIS execute the `INSTALL.sh` script.

For manual installation follow the instructions in `metis-5.1.0` and [`fMETIS`](https://github.com/ivan-pi/fmetis) directories in the repo and install [MPI](https://www.open-mpi.org/software/ompi/v4.1/) (for message passing) on your machine.


- During the installation of METIS and fMETIS, make sure to configure them with `INT=32` and `REAL=64` settings. For METIS, the header file has already been modified. For fMETIS consider the following `cmake .. -DMETIS_LIB="path/to/libmetis.a" -DREAL=64` in the `build` directory.


### Compilation

Change the variables in the Makefile to the correct directories of METIS and fMETIS: `DIR_METIS`, `DIR_FMETIS`

Compile the program in the folder `src` with the command:
```
make
```


## User-Guide
Then you can perform the following steps in a specific test folder (e.g. `tests/analytical/tet_125`):
1. Modify the file `Poly.input` to change the mesh file to read
2. Modify the file `test1.mate` to modify the total degree of basis function
3. Modify the file `problem_data_and_properties.f90` to modify the forcing term
4. Modify the file `problem_data_and_properties.f90` to change alpha (the coefficient in the definition of the penalty function) or theta to change the method
5. In case modification to the source code, recompile with
```sh
make -C ../../..
```
6. Run the parallel program with the following command:
```sh
mpirun -np 4 ../../../Lymph3D
```
where the number of processes is set to 4 in this case.


## Structure

The library includes the following files (notice that they are listed in in the order they are invoked from the main file rather than in alphabetical order):

* Lymph3D.f90: is the main file, which is structured as follows:
    - it initializes the MPI environment with the subroutine INITIALIZATION;
    - it calls the subroutine READ_INPUT_FILES to read the input file and the mesh file;
    - it calls the subroutine MAKE_PARTITION_AND_MPI_FILES to make the partition into processes and to store the mesh information among processes;
    - it defines the matrices and vectors and assembles them by calling MAKE_MATRICES_FREE and MAKE_RHS_FREE.
    - it exchanges data across processes and advances in time inside TIME_STEP_MATRIX_FREE;
    - it scatters the solution so that it can be exported;
    - it computes the modal solution for post-processing purposes;
    - it computes the errors in $L^2$ norm and DG norm with the subroutines COMPUTE_ERROR_L2_FREE and COMPUTE_ERROR_DG_FREE;
    - it exports the solution with the subroutine EXPORT_SOLUTION.

Documentation is available through doxygen.