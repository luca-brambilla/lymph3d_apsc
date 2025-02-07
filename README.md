# LYMPH3D Library

LYMPH3D is a library written in Fortran which applies the PolyDG method to solve the linear elasticity problem under both Dirichlet and Neumann boundary conditions.
The library numerically computes volume and surface integrals using Gaussian quadrature formulas, in particular by means of Gauss-Legendre nodes and weights. Given that the computational mesh consists of an agglomeration of tetrahedral and triangular elements, the implementation necessitates appropriate mappings from the reference square to the reference triangle and from the reference cube to the reference tetrahedron.

## Installation

### Prerequisites

- Fortran compiler `gfortran`
- `METIS`
- `fMETIS`
- `MPI`
- `doxygen` (for documentation generation)
- `Paraview` (for solution visualization)

For a fast installation execute the `INSTALL.sh` script

Ensure that you have installed `METIS` from `metis-5.1.0.tar.gz` in the repo, [fMETIS](https://github.com/ivan-pi/fmetis) (for mesh agglomeration) and [MPI](https://www.open-mpi.org/software/ompi/v4.1/) (for message passing) on your computer.

- During the installation of METIS and fMETIS, make sure to configure them with `INT=32` and `REAL=64` settings.

- To install `METIS`
### Compilation
Compile the program in the folder `src` with the command:
```
make
```


## User-Guide
Then you can perform the following steps:
1. Modify the file `Poly.input` to change the mesh file to read
2. Modify the file `test1.mate` to modify the total degree of basis function
3. Modify the file `problem_data_and_properties.f90` to modify the forcing term
4. Modify the file `problem_data_and_properties.f90` to change alpha (the coefficient in the definition of the penalty function) or theta to change the method
5. In case modification to the source code, recompile with
```sh
make
```
6. Run the program in a specific test folder (`TEST_TET` or `TEST_POLY`) with the following command:
```sh
mpirun -np 2 ../Lymph3D
```
where the number of processors is set to two in this case.

## Notes
- Remember to remove the file convergence_test.vtk before performing a new convergence test!

## Structure

The library includes the following files (notice that they are listed in in the order they are invoked from the main file rather than in alphabetical order):

* Lymph3D.f90: is the main file, which is structured as follows:
    - it initializes the PETSc and MPI environments with the subroutine INITIALIZATION;
    - it calls the subroutine READ_INPUT_FILES to read the input file and the mesh file;
    - it calls the subroutine MAKE_PARTITION_AND_MPI_FILES to make the partition into processes and to store the mesh information among processes;
    - it defines the PETSc matrices and vectors by calling SET_PETSC_MATRIX and SET_PETSC_VECTOR and assembles them by calling MAKE_MATRICES and MAKE_RHS.
    - it sets the algebraic solvers with the subroutine SOLVER_SETTINGS;
    - it solves the linear system by calling the PETSc function KSPSolve;
    - it scatters the solution so that it can be exported;
    - it computes the modal solution for post-processing purposes;
    - it computes the errors in L^2 norm and DG norm with the subroutines COMPUTE_ERROR_L2 and COMPUTE_ERROR_DG;
    - it exports the solution with the subroutine EXPORT_SOLUTION.

* problem_data_and_properties.f90: stores the forcing term, boundary data and exact solution. Additionally, it contains a subroutine called set_properties, responsible for configuring the properties of the problem (penalization coefficient, IP method, reaction coefficient).

* Poly_setup_MPI.f90: contains a subroutine with the same name that stores the subroutine INITIALIZATION which performs the initialization of the MPI and PETSc envinronments.

* Poly_mesh.f90: contains the definition of the structs Element, Polyhedron and Mesh_Structure.
It also contains the following subroutines:
  - allocate_Mesh_Structure to allocate the Mesh Structure, namely it sets the number of elements and allocates the connettivity matrices and the vector elem_in_poly;
  - print_Dime_Mesh_Structure to print the Mesh Structure, namely it prints the number of nodes, hexahedra, tetrahedra, prysmas, quadrilateral faces and triangular faces;

* Poly_data.f90: contains the definition of the struct Data_Structure, which stores the material file parameters such as the density and Lamé parameters, and subroutines related to print, allocation and setting of default values to the parameters.

* Poly_global.f90: contains the subroutine calc_time which takes as input the wall-time in seconds time_in_seconds and compute the corresponding time in hours, minutes and seconds. It also contains the following modules:
    - Poly_global, which stores some useful variables as head_file (the name of the input file stored in the variable), mate_file (the name of the file containing information on the materials), grid_file (the name of the file containing the mesh), and some variables related to the output and the measuring of the computational time;
    - Poly_exit_codes, Poly_fail_codes and Poly_default_codes, which are related to error and fail codes if something goes wrong reading the .input file and the .mate file;
    - qsort, which stores the qsort algorithm for sorting elements of an array;
    - local_search, which contains the subroutine GET_EL_LOC_FROM_EL_GLO that we need in order to find the local index of an element of the mesh;
    - find_poly, which contains the function FIND_TET_IN_POLY, needed to find all the tetrahedra contained in a single polyhedra.

* Poly_readfile.f90: contains the subroutines that actually reads all the input files line by line.

* READ_INPUT_FILES.f90: contains the subroutine READ_INPUT_FILES which is called directly from the main file Lymph3D and it calls all the modules we saw above in order to set the header and the material files and to define and initialize the struct Mesh_Structure.

* vet_mat_operations.f90: contains a module with the same name that stores various subroutines and functions for the computation of the determinant of a 3x3 matrix, the inverse of a 3x3 matrix and the flipping of a vector.

* MAKE_PARTITION_AND_MPI_FILES.f90: contains a subroutine with the same name that performs the partition of the mesh into the different processors and stores the local properties of the mesh into the correspondent fields of the struct Mesh_Structure. Furthermore, this subroutine generates output files in the folders FILES_MPI.
It calls also the following subroutines:
  - MESH_PARTITIONING: performs a contigous partition of the mesh into different processors using METIS and writes the mpi file elem4proc.mpi;
  - WRITE_PARTITION: writes the connettivity of tetrahedra and triangles for each processor in the mpi files con_tet.mpi and con_tri.mpi;
  - MESH_AGGLOMERATION: generates the polyhedral mesh by agglomerating the tetrahedral mesh read from the mesh file and writes the mpi file elem_in_poly.mpi;
  - CREATE_GLOBAL_POLY_MAP: allocates and stores the field elem_in_poly using the information collected in elem_in_poly.mpi;
  - FIND_POS_LOC_NODE: finds the id of the node of the mesh starting from the vertex of an element;
  - CREATE_LOCAL2GLOBAL_MAP: creates the local-to-global map dof numbering Dof_glo stored in the struct Element;
  - CREATE_LOCAL_MESH: reads the local properties of the tetrahedral mesh from the mpi files and stores them into the corresponding fields of the struct Mesh_Structure;
  - CREATE_VERT_LIST: stores the coordinates of the vertices of the tetrahedra locally;
  - CREATE_POLY_LIST: stores the map poly_loc2glo to transition from local polyhedra to the global ones and initializes the other fields of the structure Polyhedron;
  - CREATE_NORMAL_FACE: computes the coordinates of the normal vector for each face and determines the area associated with each face;
  - CREATE_BBOX_EL: performs the computation of the coordinates of the bounding box for each polyhedron;
  - CREATE_NEIGH_EL_TRIA: stores the information of the neighbouring tetrahedra and polyhedra by reading data from con_tri.mpi;
  - WRITE_MESH_INFO: writes the mesh properties for the elements contained in each processor in the mpi files mesh.mpi;
  - WRITE_MESH_VISUALIZATION_VTK: it is defined in MOD_VTK.f90 (see below).

* MESH_CORRECTION.f90: ensures that each polyhedron in the mesh contains at least one tetrahedron.

* Poly_ref_mappings.f90: contains the module Poly_ref_mappings which defines different mappings that are then composed to project the quadrature nodes onto the physical frame. It includes the following subroutines:
  - jacobians, which performs the computation of the matrices of the reference map Fk from the reference tetrahedron (0,0,0), (1,0,0), (0,1,0), (0,0,1) to the physical tetrahedron, the corresponding determinant of the Jacobian Jdet and its inverse Jinv;
  - tria2tetfaces_maps, which computes the maps from the two-dimensional reference triangle (0,0), (1,0), (0,1) to the faces of the three-dimensional reference tetrahedron and their inverses.
  - quadrature_map_2D, which computes the map from the reference square (-1,1)^2 to the reference triangle (0,0), (1,0), (0,1);
  - quadrature_map_3D which computes the map from the reference cube (-1,1)^3 to the reference tetrahedron (0,0,0), (1,0,0), (0,1,0), (0,0,1);

* basis_function.f90: contains a module with the same name which stores the following subroutines:
  - blist: returns the list of the degrees of monomials of the Np basis functions up to a total degree p;
  - quadrature: computes Gauss-Legendre quadrature nodes and weights over the reference cube (-1,1)^3 and reference square (-1,1)^2;
  - LegendreP_notscaled: evaluates the non-scaled Legendre Polynomial Ln(x) in one dimension on the interval int, which corresponds to the edge of the bounding box in one particular direction, at points x, of order given by blist considering that the total degree must be p;
  - LegendreP: evaluates the scaled Legendre Polynomial Ln(x) in one dimension on the interval int, which corresponds to the edge of the bounding box in one particular direction, at points x, of order given by blist considering that the total degree must be p. It makes use of the subroutine LegendreP_notscaled defined before;
  - GradLegendreP: evaluates the derivative L′n(x) of the scaled Legendre Polynomial Ln(x) in one dimension on the interval int, which corresponds to the edge of the bounding box in one particular direction, at points x, of order given by blist considering that the total degree must be p. It makes use of the subroutine LegendreP_notscaled defined before;
  - basis: evaluates the basis functions and their partial derivatives at the three-dimensional quadrature nodes for every element;
  - basis_boundary: evaluates the basis functions and their partial derivatives, for every face F of the two neighbouring tetrahedra E1 and E2, at the two-dimensional quadrature nodes;

* SET_PETSC_SYSTEM.f90: contains a module with the same name that contains the subroutines SET_PETSC_VECTOR and SET_PETSC_MATRIX which create and initialize vector and matrices within the PETSc environment, respectively.

* assemble_element.f90: assembles the local matrices and the local rhs vector by iterating over either the three-dimensional quadrature nodes or the two-dimensional quadrature nodes. In particular, this module contains the following 5 subroutines:
  - MAKE_STIFFNESS_VOLUME: assembles V_loc, the local term of the stiffness matrix approximating the integral on the tetrahedron;
  - MAKE_STIFFNESS_FACE: assembles I_loc, S_loc, the terms of the stiffness matrix approximating the integrals on the faces of the tetrahedron, and IN_loc and SN_loc on the faces of the neighbouring tetrahedra;
  - MAKE_RHS_VOLUME: assembles rhs_tet_loc, the local rhs term approximating the integral on the tetrahedron;
  - MAKE_RHS_FACE: assembles rhs_face_bd_loc, the local rhs term that approximates the integral on the faces of the tetrahedron that are boundary faces;
  - MAKE_MASS_VOLUME: assembles of M_loc, the local term of the mass matrix that approximates the integral on the tetrahedron.

* MAKE_MATRICES.f90: contains a subroutine with the same name which assembles the global stiffness matrix pestc_stiff, mass matrix petsc_mass and DG matrix mat_dg by performing a loop on the elements.

* MAKE_RHS.f90: contains a subroutine with the same name which assembles the global rhs vector petsc_rhs by performing a loop on the elements.

* SOLVER_SETTINGS.f90: contains a subroutine with the same name that simultaneously sets the solver (direct or iterative) and, optionally, the preconditioner. It also creates the object KSP for solving the linear system using the PETSc library.

* COMPUTE_MODAL_COEFFICIENTS.f90: contains a subroutine with the same name that computes the modal coefficients of the exact solution. These coefficients are subsequently utilized in the main file to compute the modal solution.

* MOD_VTK.f90: contains a module with the same name that stores different subroutines for the creation of vtk files. Specifically, the file MOD_VTK.f90 contains the following subroutines:
  - WRITE_ERRORS: writes the degree of the basis functions, grid size and values of $L^2$ and DG errors wo separate files located in two different folders within the POST-PROC folder, namely CONVERGENCE_TEST and ERRORS. In the former case, these values are appended to the existing content of the file, facilitating their use in convergence tests. In the latter case, the existing content of the file is replaced by the new values;
  - WRITE_MESH_VISUALIZATION: writes the partition of the mesh into different processors in mesh_partition.mpi and the agglomeration of the mesh in mesh_agglomeration.mpi;
  - WRITE_SOLUTION_VISUALIZATION: writes the values of the numerical solution at the vertices of the tetrahedra in num_sol.mpi;
For WRITE_MESH_VISUALIZATION and WRITE_SOLUTION_VISUALIZATION, the corresponding vtk files generated can be visualized using software such as [Paraview](https://www.paraview.org/).

* solution_processing.f90: contains a module with the same name that performs the post-processing. It stores the following subroutines:
  - EXPORT_SOLUTION: evaluates the solution function at the vertices of the tetrahedra and then calls the subroutine WRITE_SOLUTION_VISUALIZATION to write these values to a vtk file;
  - COMPUTE_ERROR_L2: computes the error in the $L^2$ norm, utilizing both the (exact) modal solution and the numerical solution.
  - COMPUTE_ERROR_DG: computes the error in the DG norm, again utilizing both the (exact) modal solution and the numerical solution.

In conclusion, alongside the main file, a folder named convergence test verify the convergence of the method, with respect to both the grid size h and the degree p of the basis functions. This verification is achieved by extracting error values from the vtk file generated by the subroutine WRITE ERRORS. The folder then computes the convergence rates and produces plots of the errors in different png files using [Gnuplot](http://www.gnuplot.info/)