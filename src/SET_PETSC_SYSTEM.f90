!> setup and allocate PETSc matrices and vectors
module SET_PETSC_SYSTEM

#include<petsc/finclude/petscksp.h>

    use petscksp
    use Poly_setup_mpi
    use Poly_global

    implicit none

    contains

    ! Create and initialize a vector in PETSc
    subroutine SET_PETSC_VECTOR(petsc_vec, local_dof, global_dof)

        Vec :: petsc_vec

        integer(kind=4) :: global_dof, local_dof

        PetscCall(VecCreate(PETSC_COMM_WORLD, petsc_vec, mpi_ierr))
        PetscCall(VecSetSizes(petsc_vec, local_dof, global_dof, mpi_ierr))
        PetscCall(VecSetFromOptions(petsc_vec,mpi_ierr))

    end subroutine SET_PETSC_VECTOR

    ! Create and initialize a matrix in PETSc
    subroutine SET_PETSC_MATRIX(petsc_mat, local_dof, global_dof)

        Mat :: petsc_mat

        !integer(kind=4) :: Np
        integer(kind=4) :: global_dof, local_dof
        integer (kind=4) :: nrows

        nrows = 3000; ! 1000
        if (global_dof <= 3000 ) nrows = global_dof
        ! nrows = global_dof / mpi_np / 100

        ! 3 for 3D, 4 for neighbor tetrahedra + 1 element itself
        !nrows = Np*3*5

        ! Create sparse matrix stiff
        PetscCall(MatCreate(PETSC_COMM_WORLD, petsc_mat, mpi_ierr))

        ! Set matrix stiff in AIJ (compressed sparse row) format
        ! Let PETSc automatically decide how to distribute the matrix among processes with PETSC_DECIDE
        ! (otherwise set local_dof by yourself, it is better for parallel computation)
        PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,local_dof,local_dof,global_dof,global_dof,nrows,PETSC_NULL_INTEGER_ARRAY,nrows,PETSC_NULL_INTEGER_ARRAY,petsc_mat,mpi_ierr))

        ! Allows us to configure various options for the matrix through command-line arguments or a configuration file
        PetscCall(MatSetFromOptions(petsc_mat, mpi_ierr))

    end subroutine SET_PETSC_MATRIX

end module SET_PETSC_SYSTEM