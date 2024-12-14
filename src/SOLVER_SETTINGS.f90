!> Set PETSc solver settings: type of solver (direct, iterative), algorithm and preconditioner
subroutine SOLVER_SETTINGS(stiff, ksp, pc)

#include<petsc/finclude/petscksp.h>

    use petscksp
    use Poly_setup_mpi
    use Poly_global
    use problem_data_and_properties

    implicit none

!-------------------------------------------------------------------------------

    ! PETSC

    !> FEM Matrices for acoustic materials
    Mat :: stiff ! stiffness matrices
    !Mat :: RR

    !> Further PetSc variables
    !PetscBool :: flg
    PetscReal :: rel_tolerance
    PetscInt :: max_krylov_it

    !> Index variables
    !PetscInt icntl, ival

    !> PETSc Solver context
    KSP :: ksp ! Krylov solver context
    PC :: pc ! Preconditioner object for the solver

    !> Flag variables
    integer (kind=4) :: solver_type, preconditioner_type

    real(kind=8) :: alpha, theta, c

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! LINEAR SOLVER SETTINGS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Set up the solver for A ------------------------------------------------

    !> 0 :: Direct solver
    !> 1 :: Iterative solver
    solver_type = 1

    !> Create solver object ksp
    if (IS_MatrixFree .eqv. .true.) then
        !! CREATE ONLY ONCE OUTSIDE
        !PetscCall(KSPCreate(PETSC_COMM_SELF, ksp, mpi_ierr))
    else
        PetscCall(KSPCreate(PETSC_COMM_WORLD, ksp, mpi_ierr))
    endif

    !> Set the linear operators for the KSP object
    !> Choose A itself as preconditioner as default
    PetscCall(KSPSetOperators(ksp, stiff, stiff, mpi_ierr))

    !> Choose the concrete solver
    if (solver_type == 0) then !> Direct solver

        ! Set the type of the KSP solver
        ! KSPPREONLY: apply a preconditioner to the linear system but not to perform any iterative solving
        PetscCall(KSPSetType(ksp, KSPPREONLY, mpi_ierr))
        !> Use a specific preconditioner, that is actually a direct solver
        !> 0 :: LU
        !> 1 :: Cholesky
        preconditioner_type = 0

    elseif (solver_type == 1) then !> Iterative solver

        call set_properties(alpha, theta, c)

        if(theta == -1) then
            !> KSPCG: Conjugate Gradient method
            PetscCall(KSPSetType(ksp, KSPCG, mpi_ierr))
        else
            !> KSPGMRES: GMRES method
            PetscCall(KSPSetType(ksp, KSPGMRES, mpi_ierr))
        endif

    endif

    !> Specify the iterative solver settings
    if (solver_type > 0) then
        !> Use a specific preconditioner, that is NOT a direct solver
        !> 2 :: ILU
        !> 3 :: Incomplete Cholesky
        !> 4 :: Jacobi
        !> 5 :: SOR
        preconditioner_type = 5

        !> Iterative solver settings
        rel_tolerance = 1.d-16
        max_krylov_it = 100 ! 100

        !> Set the convergence tolerances and iteration limits for the KSP solver
        !> Absolute tolerance and divergence tolerance are set to the default value
        PetscCall(KSPSetTolerances(ksp, rel_tolerance, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, max_krylov_it, mpi_ierr))
        PetscCall(KSPSetFromOptions(ksp, mpi_ierr))

        !> Set the flag that indicates whether the initial guess for the solver should be considered nonzero
        !> In this case, it is set to PETSC_TRUE, so that the initial guess is nonzero
        !> In particular, use the previous solution as initial guess for the next solve step
        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, mpi_ierr))
    endif

    ! ------------------------------------------------------------------------------

    !> Set up preconditioner -------------------------------------------------------
    !> Gets chosen by the linear solver

    !> Get the preconditioner associated with the KSP solver
    !> pc is a pointer to a variable (of type PC) that will store the preconditioner associated with the KSP object after the function call
    PetscCall(KSPGetPC(ksp, pc, mpi_ierr))

    !> Choose the concrete preconditioner
    if (preconditioner_type == 0) then
        !> Use LU decomposition as a preconditioner
        PetscCall(PCSetType(pc, PCLU, mpi_ierr))
        ! PetscCall(PCSetType(pc, PCCHOLESKY, mpi_ierr))

        !> Use MUMPS to perform the LU-decomposition
        ! PetscCall(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,mpi_ierr))
        ! PetscCall(PCFactorSetUpMatSolverType(pc,mpi_ierr))
        ! PetscCall(PCFactorGetMatrix(pc,RR,mpi_ierr))
        ! PetscCall(KSPSetFromOptions(ksp,mpi_ierr))

        ! PetscCall(PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS, mpi_ierr))
        ! PetscCall(PCFactorSetUpMatSolverPackage(pc, mpi_ierr))
        ! PetscCall(PCFactorGetMatrix(pc, RR, mpi_ierr))

        !> MUMPS settings:
        !> Index    Value of the index
        ! icntl = 14
        ! ival  = 70
        ! PetscCall(MatMumpsSetIcntl(RR, icntl, ival, mpi_ierr))

    elseif (preconditioner_type == 1) then
        !> Use Cholesky factorization as a preconditioner
        PetscCall(PCSetType(pc, PCCHOLESKY, mpi_ierr))

    elseif (preconditioner_type == 2) then
        !> Use ILU as a preconditioner
        PetscCall(PCSetType(pc, PCILU, mpi_ierr))

    elseif (preconditioner_type == 3) then
        !> Use Incomplete Cholesky factorization as a preconditioner
        PetscCall(PCSetType(pc, PCICC, mpi_ierr))

    elseif (preconditioner_type == 4) then
        !> Use Jacobi preconditioner
        PetscCall(PCSetType(pc, PCJACOBI, mpi_ierr))

    elseif (preconditioner_type == 5) then
        !> Use SOR method as a preconditioner
        PetscCall(PCSetType(pc, PCSOR, mpi_ierr))
    endif

end subroutine SOLVER_SETTINGS