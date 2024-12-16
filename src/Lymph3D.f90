!   Author: Lorenzo Gorlezza
!           Luca Brambilla
!   This file is part of the library LYMPH3D

!> @brief Lymph3D (Discontinuous Galerkin methods on polyhedral meshes for PDE problems)

! Here starts the code Lymph3D

program Lymph3D

#include<petsc/finclude/petscksp.h>

    use petscksp
    use petscmat
    use mpi
    use Poly_setup_mpi
    use problem_data_and_properties
    use Poly_global
    use Poly_data
    use Poly_mesh
    use solution_processing
    use SET_PETSC_SYSTEM
    use matrix_free
    use exchange_data
    !use MOD_MPI_CUSTOM
    use global_parameters
    use utilities

    implicit none

    ! read additional flags
    ! integer :: iarg
    ! character(len=256) :: arg

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEFINITION OF PETSC VARIABLES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    Mat :: petsc_stiff, petsc_mass, mat_dg, petsc_mass_modal ! stiffness, mass and DG matrices, reconstruction matrix
    Vec :: petsc_sol, petsc_rhs ! solution and rhs vectors for the algebraic system
    Vec :: petsc_uex, petsc_modal_coeff_uex ! modal (exact) solution and its coefficients

    Vec :: petsc_v_tmp, petsc_f ! temp vectors to store solution for sums
    Vec :: petsc_u0, petsc_v0 ! initial condition read from data

    !PetscViewer viewer ! abstract PETSc object for displaying PETSc objects and their data

    KSP :: ksp, ksp2, ksp3 ! linear system solvers
    PC :: pc, pc2, pc3 ! preconditioners

    logical :: IsPoly
    integer(kind=4) :: local_dof, global_dof
    integer(kind=4) :: Np, p, Npoly
    integer(kind=4), dimension(:), allocatable :: petsc_num
    integer(kind=4), dimension(:), allocatable :: nnod_num
    real(kind=8) :: err_L2, err_L2_mpi
    real(kind=8) :: err_DG, err_DG_mpi
    real(kind=8) :: hmax, hmax_mpi
    real(kind=8), pointer :: sol_ptr(:)
    !real(kind=8), dimension(:), allocatable :: u_loc, u_glo
    real(kind=8), dimension(:,:), allocatable :: u
    integer(kind=4), dimension(:), allocatable :: gathered_sizes, displacements

    !real(kind=8) :: E, nu
    integer(kind=4) :: i, mat_id, ie_loc

    type(Data_Structure) :: PolyData
    type(Mesh_Structure) :: PolyMesh

    integer(kind=4) :: num_dt = 1     ! number of iteration
    real(kind=8) :: t             ! time variable

    real(kind=8), parameter :: one = 1.0


    ! each local has a matrix
    real(kind=8), dimension(:,:,:), allocatable :: M_loc
    real(kind=8), dimension(:,:,:), allocatable :: M_modal_loc
    ! real(kind=8), dimension(:,PolyMesh%num_elem_loc,DIM,DIM,Np,Np), allocatable, intent(out) :: K_loc
    real(kind=8), dimension(:,:), allocatable :: v0_loc
    real(kind=8), dimension(:,:), allocatable :: u0_loc
    real(kind=8), dimension(:,:), allocatable :: un_loc
    real(kind=8), dimension(:,:), allocatable :: usol_loc
    real(kind=8), dimension(:,:), allocatable :: rhs_loc


    real(kind=8), dimension(:,:,:,:), allocatable :: K_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: A_dg_loc

    ! type(KRowArray), dimension(:), allocatable :: K_loc
    ! type(KRowArray), dimension(:), allocatable :: A_dg_loc

    ! integer(kind=4), pointer :: internal_neigh(:)
    ! integer(kind=4), allocatable :: tmp_pointer(:)
    integer(kind=4) :: n_neigh !, max_faces

    ! IsTime_dependent = .false.

    type(PetscMatStruct), dimension(:,:), allocatable:: massa
    type(PetscMatStruct), dimension(:,:), allocatable:: massa_modale
    real(kind=8), dimension(:), allocatable :: tmp
    Mat :: petsc_m_tmp

    !integer(kind=4) :: j,k,m,n,row,col
    integer(kind=4) :: row

    real(kind=8), dimension(:), allocatable :: prova_in, prova_out
    integer(kind=4) :: tmp_size, unit_print
    type(ScatteredArray), dimension(:,:), allocatable :: send_data, recv_data
    
    integer(kind=4) :: E1, E2, iface, ie_neigh_loc, j
    logical :: is_E2_local
    real(kind=8), dimension(:), pointer :: v_ptr
    real(kind=8), dimension(:,:), pointer :: m1_ptr, m2_ptr
    
    real(kind=8), dimension(:), allocatable :: x, b
    real(kind=8), dimension(:,:), allocatable :: A, R

    ! read parameter
    ! iarg = getarg(1,arg)
    ! open(unit=10, file=arg, status="new")
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    allocate(mpi_stat(MPI_STATUS_SIZE))

    call INITIALIZATION()

    ! call MPI_OP_CREATE(MPI_ZERO_OVERWRITE, .TRUE., MPI_ZERO_OVERWRITE_OP, mpi_user_reduction_error)
    ! call MPI_OP_CREATE(MPI_OVERWRITE_BY_NEW, .false., MPI_OVERWRITE_BY_NEW_OP, mpi_user_reduction_error)

    start = MPI_WTIME()

    if(mpi_id == 0) then
        write(*,'(A)')''
        write(*,'(A)')'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>'
        write(*,'(A)')'<                                                     >'
        write(*,'(A)')'<                        Lymph3D                      >'
        write(*,'(A)')'<                                                     >'
        write(*,'(A)')'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>'
        write(*,'(A)')''
    endif

    ! check FILES_MPI directory
    if (mpi_id == 0) call CHECK_MPI_FILES

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     READ INPUT FILES AND ALLOCATE VARIABLES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call READ_INPUT_FILES(PolyData,PolyMesh)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'--------------------Type of problem--------------------'
        write(*,'(A)')
        write(*,'(A,L)')'Matrix free: ',IS_MatrixFree
        write(*,'(A,L)')'Time dependent problem: ', IsTime_dependent
        if (IsTime_dependent .eqv. .false.) then
            write(*,'(A)')'Stationary problem'
        else
            write(*,'(A,F8.5)')'start time [s] : ', start_time
            write(*,'(A,F8.5)')'end time   [s] : ', stop_time
            write(*,'(A,F8.5)')'timestep   [s] : ', time_step
        endif
        write(*,'(A)')
    endif

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CHOOSE TO SOLVE WITH TETRAHEDRA
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    IsPoly = .false.
    Npoly = PolyMesh%num_poly

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CHOOSE TO SOLVE WITH POLYHEDRA
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! set a number of polyhedra which may be different from the one read from the mesh file
    ! Npoly=900

    if (Npoly /= PolyMesh%num_elem) IsPoly = .true.

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     PARTITION OF THE GRID AND GENERATION OF LOCAL CONNECTIVITY
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (mpi_id == 0)  &
        write(*,'(A)') '---------------------Partitioning----------------------'

    call MAKE_PARTITION_AND_MPI_FILES(PolyData, PolyMesh, Npoly)

    ! if (mpi_id==0) then
    !     print *, PolyMesh%elem_glo2loc
    ! endif

    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call FLUSH

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    allocate(A(Np,Np), x(Np), b(Np))

    ! Initialize random seed
    call random_seed()
    ! Generate a symmetric random matrix
    do i = 1, Np
       do j = i, Np
          call random_number(A(i, j))
          A(j, i) = A(i, j)  ! Symmetric property
       end do
    end do
    ! Add a value to the diagonal to make it diagonally dominant
    do i = 1, Np
       A(i, i) = A(i, i) + Np  ! Increase diagonal for positive definiteness
    end do

    x = 1.0
    b = matmul(A, x)

    print *, 'factor matrix'
    R = cholesky(A, Np)
    print *, 'solve system'
    x = solve_LU(transpose(R),R,b,Np)

    print *, x

    call STOP_LYMPH3D

    tmp_size = sum(PolyMesh%num_elem_inter_comm(mpi_id+1,:))
    allocate(prova_in(PolyMesh%num_elem_loc*DIM*Np))
    allocate(prova_out(tmp_size*DIM*Np))

    prova_in=0
    do i=1,PolyMesh%num_elem_loc*DIM*Np
        prova_in(i) = 1000000*mpi_id + i
    enddo

    !print *, 'proc:', mpi_id, 'data out: ', prova_in

    call MPI_EXCHANGE_ALLOCATE(PolyMesh, send_data, recv_data)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    call MPI_EXCHANGE_DOF(PolyMesh, prova_in, prova_out, send_data, recv_data)

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call MPI_EXCHANGE_DEALLOCATE(PolyMesh, send_data, recv_data)
    print *, 'end exchange'
    !if (mpi_id == 0) print *, 'proc:', mpi_id, 'data out: ', prova_out

    !call STOP_LYMPH3D

    call WRITE_MESH_VISUALIZATION_VTK(PolyMesh%num_elem_loc, PolyMesh, mpi_id)

    !! MAY VARY FROM ELEMENT TO ELEMENT !!!
    Np = PolyMesh%Elem_loc(1)%NDof_elem ! number of degrees of freedom of each element
    p = PolyMesh%Elem_loc(1)%Degree ! order of basis functions
    global_dof = DIM*PolyMesh%num_poly*Np;
    if (mpi_id == 0 ) then
        print *, 'Global number of degrees of freedom: ', global_dof
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    call FLUSH
    print *, 'Process ID:', mpi_id
    write(*,'(A28,I10,A28,I10)') 'Global number of polyhedra:', PolyMesh%num_poly, &
        'Local number of polyhedra:', PolyMesh%num_poly_loc
    print *, 'Done reading mesh'

    Np = PolyMesh%Elem_loc(1)%NDof_elem ! number of degrees of freedom of each element
    p = PolyMesh%Elem_loc(1)%Degree ! order of basis functions

    global_dof = DIM*PolyMesh%num_poly*Np;
    local_dof = DIM*PolyMesh%num_poly_loc*Np;
    PRINT *, 'Local number of degrees of freedom: ', local_dof
    call FLUSH

    ! PETSc numbering starting from 0 not 1
    allocate(petsc_num(global_dof))
    do i=1,global_dof
        petsc_num(i) = i-1
    end do

    ! E = 32e9
    ! nu = 0.2

    ! set the value of rho, lambda and mu
    do mat_id=1,PolyData%nmat
        ! PolyData%prop_mat(mat_id,1) = 2400
        ! PolyData%prop_mat(mat_id,2) = E*nu / ((1+nu)*(1-2*nu))
        ! PolyData%prop_mat(mat_id,3) = E / (2 * (1+nu))
        PolyData%prop_mat(mat_id,1) = 1
        PolyData%prop_mat(mat_id,2) = 1
        PolyData%prop_mat(mat_id,3) = 1
    enddo

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     SET PETSC VECTORS AND MATRICES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'--------------------Compute solution-------------------'
    endif
    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    !IS_MatrixFree = .true.
    if (IS_MatrixFree .eqv. .true.) then
        print *, 'matrix-free - set matrices and vectors'
        ! allocate(internal_neigh(PolyMesh%num_elem_loc))
        ! allocate(A_dg_loc(PolyMesh%num_elem_loc))
        ! allocate(K_loc(PolyMesh%num_elem_loc))

        n_neigh = PolyMesh%Elem_loc(1)%num_faces
        !! WASTE OF MEMORY... MAKE SCATTERED SIZE VECTOR?
        allocate( K_loc(PolyMesh%num_elem_loc, n_neigh+1, DIM*Np, DIM*Np) )
        allocate( A_dg_loc(PolyMesh%num_elem_loc, n_neigh+1, DIM*Np, DIM*Np) )
        allocate( rhs_loc(PolyMesh%num_elem_loc, DIM*Np) )

        ! allocate the struct containing PETSc Mat for the mass data matrix-free form
        allocate(massa(PolyMesh%num_elem_loc,DIM))
        allocate(massa_modale(PolyMesh%num_elem_loc,DIM))
        call SET_PETSC_MASS_MATRIX_FREE(PolyMesh%num_elem_loc, Np, massa)
        call SET_PETSC_MASS_MATRIX_FREE(PolyMesh%num_elem_loc, Np, massa_modale)

        ! create local vector to each proces for matrix vector multiplicaton
        ! mass matrix-free temporary vector and linear system solution
        PetscCallA(VecCreate(PETSC_COMM_SELF, petsc_v_tmp, mpi_ierr))
        PetscCallA(VecSetSizes(petsc_v_tmp, Np, Np, mpi_ierr))
        PetscCallA(VecSetFromOptions(petsc_v_tmp, mpi_ierr))

        PetscCallA(VecCreate(PETSC_COMM_SELF, petsc_sol, mpi_ierr))
        PetscCallA(VecSetSizes(petsc_sol, Np, Np, mpi_ierr))
        PetscCallA(VecSetFromOptions(petsc_sol, mpi_ierr))

        ! create local matrix to each process for matrix vector multiplication
        ! mass matrix-free temporary matrix

        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, PETSC_NULL_SCALAR_ARRAY, petsc_m_tmp, mpi_ierr))

        ! PetscCallA(MatCreate(PETSC_COMM_SELF, petsc_m_tmp, mpi_ierr))
        ! PetscCallA(MatSetSizes(petsc_m_tmp, Np, Np, Np, Np, mpi_ierr))
        ! PetscCallA(MatSetFromOptions(petsc_m_tmp, mpi_ierr))
        ! PetscCallA(MatSetUp(petsc_m_tmp, mpi_ierr)) !! what?
        ! PetscCallA(MatAssemblyBegin(petsc_m_tmp,MAT_FINAL_ASSEMBLY,mpi_ierr))
        ! PetscCallA(MatAssemblyEnd(petsc_m_tmp,MAT_FINAL_ASSEMBLY,mpi_ierr))

    else
        call FLUSH
        print *, 'SET PETSC MATRICES AND VECTORS'

        call SET_PETSC_MATRIX(petsc_stiff, local_dof, global_dof)
        call SET_PETSC_MATRIX(petsc_mass, local_dof, global_dof)
        call SET_PETSC_MATRIX(mat_dg, local_dof, global_dof)
        call SET_PETSC_MATRIX(petsc_mass_modal, local_dof, global_dof)

        call SET_PETSC_VECTOR(petsc_rhs, local_dof, global_dof)
        call SET_PETSC_VECTOR(petsc_sol, local_dof, global_dof)
        call SET_PETSC_VECTOR(petsc_uex, local_dof, global_dof)

        call SET_PETSC_VECTOR(petsc_v_tmp, local_dof, global_dof)
        call SET_PETSC_VECTOR(petsc_f, local_dof, global_dof)
        call SET_PETSC_VECTOR(petsc_u0, local_dof, global_dof)
        call SET_PETSC_VECTOR(petsc_v0, local_dof, global_dof)
    end if

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ASSEMBLE MATRICES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (IS_MatrixFree .eqv. .true.) then
        print *, 'assemble matrix-free matrices'

        ! make all matrices
        call MAKE_MATRICES_FREE(PolyMesh, PolyData, global_dof, Np, K_loc, A_dg_loc, massa, massa_modale, n_neigh)

    else
        print *, 'ASSEMBLE PETSC MATRICES'

        call MAKE_MATRICES(PolyMesh, PolyData, petsc_num, global_dof, Np, petsc_stiff, petsc_mass, mat_dg, petsc_mass_modal)

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat_stiff.m',viewer,mpi_ierr))
        ! call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, mpi_ierr)
        ! PetscCallA(MatView(petsc_stiff,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat_mass.m',viewer,mpi_ierr))
        ! call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, mpi_ierr)
        ! PetscCallA(MatView(petsc_mass,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat_dg.m',viewer,mpi_ierr))
        ! call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, mpi_ierr)
        ! PetscCallA(MatView(mat_dg,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'mat_mass_modal.m',viewer,mpi_ierr))
        ! call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, mpi_ierr)
        ! PetscCallA(MatView(petsc_mass_modal,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
    end if

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ASSEMBLE RIGHT HAND SIDE
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (IS_MatrixFree .eqv. .true.) then

        print *, 'assemble matrix-free RHS'
        call MAKE_RHS_FREE(PolyMesh, PolyData, global_dof, Np, rhs_loc)

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call flush
        ! print *, mpi_id, rhs_loc(1,:)
        call flush
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    else
        print *, 'ASSEMBLE RHS'

        call MAKE_RHS(PolyMesh, PolyData, petsc_num, global_dof, Np, petsc_rhs)

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_rhs',viewer,mpi_ierr))
        ! PetscCallA(VecView(petsc_rhs,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
    end if

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     SETTING SOLVERS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (IS_MatrixFree .eqv. .true.) then

        print *, 'SET LOCAL matrix-free solvers'

        PetscCallA(KSPCreate(PETSC_COMM_SELF, ksp, mpi_ierr))
        call SOLVER_SETTINGS(petsc_m_tmp, ksp, pc)
    else
        print *, 'SETTING SOLVERS'

        call SOLVER_SETTINGS(petsc_stiff, ksp, pc)
        call SOLVER_SETTINGS(petsc_mass, ksp2, pc2)
        call SOLVER_SETTINGS(petsc_mass_modal, ksp3, pc3)
    end if

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CALLING SOLVER
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! allocate(internal_neigh(10))
    ! do i=1,10
    !     allocate(tmp_pointer(2))
    !     internal_neigh(i) => tmp_pointer
    !     deallocate(tmp_pointer)

    !     if (.not. allocated(internal_neigh(i))) then
    !         print *, "Error: internal_neigh(", i, ")%values is not allocated"
    !         stop
    !     else
    !         print *, size(internal_neigh(i))
    !     endif
    ! end do
    ! n_neigh = size(internal_neigh)

    IsSave_output = .true.

    !! CHECK IF START AT num_dt=0
    t = 0.0
    num_dt = 1

    !stop_time = SQRT2 / 4.0 + 1.0 * SQRT2
    time_step = 0.001
    num_dt_mon = 20

    !stop_time = SQRT2 / 4.0 + 9.0 * SQRT2 ! 10 peaks
    stop_time = SQRT2 / 4.0 + 0.5 * SQRT2
    !stop_time = 0.1
    
    dt2 = time_step*time_step
    half_dt2 = 0.5*dt2

    ! --------------------- MATRIX FREE -----------------------
    if (IS_MatrixFree .eqv. .true.) then

        print *, 'matrix free solver'
        if(mpi_id == 0) print *, "                   TIME LOOP START                   "

        ! allocate nnod_num, gathered_sizes, displacements, u
        call PREPROCESS_SOLUTION_MATRIX_FREE(PolyMesh, local_dof, nnod_num, gathered_sizes, displacements, u)

        if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", 0, " Time: ", t

        if (mpi_id==0) print *, "Assemble initial conditions"
        ! initial conditions
        allocate(u0_loc(PolyMesh%num_poly_loc, DIM*Np))
        allocate(un_loc(PolyMesh%num_poly_loc, DIM*Np))
        allocate(v0_loc(PolyMesh%num_poly_loc, DIM*Np))
        allocate(tmp(DIM*Np))

        call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, massa_modale, ic_displacement, u0_loc)
        call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, massa_modale, ic_velocity, v0_loc)

        ! SAVE SOLUTION
        if (IsSave_output .eqv. .true.) then
            if (mpi_id==0) print *, '--- IC displacement ---'
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            call POST_PROCESS_MATRIX_FREE(PolyMesh, u0_loc, u, gathered_sizes, displacements)
            call EXPORT_SOLUTION(PolyMesh, u, IsPoly, 0)
        endif
        if (IsSave_output .eqv. .true.) then
            if (mpi_id==0) print *, '---   IC velocity   ---'
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            call POST_PROCESS_MATRIX_FREE(PolyMesh, v0_loc, u, gathered_sizes, displacements)
            call EXPORT_SOLUTION(PolyMesh, u, IsPoly)

            ! if (mpi_id==0) call SAVE_VECTOR(v0_loc,PolyMesh%num_poly_loc,'v0_pre.txt')
            ! if (mpi_id==0) call SAVE_VECTOR(u,Np,'v0_post.txt')

        endif

        ! call SAVE_MATRIX_PETSC(massa(1,1)%data, Np, Np, 'massa.txt')
        ! call SAVE_MATRIX_PETSC(massa_modale(1,1)%data, Np, Np, 'massa_modale.txt')
        ! if (mpi_id==0) call SAVE_MATRIX_F90(K_loc(1,:,:,:), DIM*Np, 'rigidezza.txt')

        ! u_1 = M^-1(dt^2/2 * f_0 - dt^2/2*A*u_0) + u0 + dt*v_0
        !! refactor matrices
        !! select correct u0_loc

        call FLUSH
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        if(mpi_id == 0) print *, ""
        if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", num_dt, " Time: ", t
        call FLUSH
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

        !! computation on tetra or on poly??? solution dof on poly

        call FIRST_TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, massa, ksp, rhs_loc, u0_loc, v0_loc, un_loc, petsc_m_tmp, petsc_v_tmp, petsc_sol)

        if (IsSave_output .eqv. .true.) then
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
            call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
        endif

        num_dt = num_dt + 1
        t = t + time_step

        tmatrix = 0.0
        tsolve = 0.0
        tstiffness = 0.0
        tvector = 0.0
        tsset = 0.0

        ! loop start
        do while (t <= stop_time)
            if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", num_dt, " Time: ", t
            
            call TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, massa, ksp, rhs_loc, u0_loc, un_loc, v0_loc)
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

            ! update solution
            u0_loc = un_loc
            un_loc = v0_loc

                ! SAVE SOLUTION
            if ( (IsSave_output .eqv. .true.) .and. (mod(num_dt, num_dt_mon) == 0) ) then
                call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
                call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
                call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
            endif

            ! update time
            num_dt = num_dt + 1
            t = t + time_step
        end do

        call calc_time(time_hour, time_min, time_sec, int(tstiffness))
        print *, 'tstiff  = ', time_min,' m ' , time_sec,' s'
        call calc_time(time_hour, time_min, time_sec, int(tsolve))
        print *, 'tsolve  = ', time_min,' m ' , time_sec,' s'
        call calc_time(time_hour, time_min, time_sec, int(tmatrix))
        print *, 'tmatrix = ' , time_min,' m ' , time_sec,' s'
        call calc_time(time_hour, time_min, time_sec, int(tvector))
        print *, 'tvector = ' , time_min,' m ' , time_sec,' s'
        call calc_time(time_hour, time_min, time_sec, int(tsset))
        print *, 'tsset   = ' , time_min,' m ' , time_sec,' s'

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

        if (IsSave_output .eqv. .true.) then
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
            call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
        endif

        !! STOP
        call STOP_LYMPH3D

    ! ----------------------------- PETSc -------------------------------------
    else

        ! allocate nnod_num, gathered_sizes, displacements, u
        call PREPROCESS_SOLUTION(PolyMesh, local_dof, nnod_num, gathered_sizes, displacements, u)

        if (IsTime_dependent .eqv. .false.) then
            ! Au=f
            PetscCallA(KSPSolve(ksp, petsc_rhs, petsc_sol, mpi_ierr))

        else
            if(mpi_id == 0) print *, "                   TIME LOOP START                   "

            if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", 0, " Time: ", t

            print *, "Assemble initial conditions"

            ! vectors initial conditions
            ! compute modal coefficients of initial conditions
            call COMPUTE_MODAL_COEFFICIENTS_GEN(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_u0, ic_displacement)

            ! SAVE SOLUTION
            if (IsSave_output .eqv. .true.) then
                if (mpi_id==0) print *, 'displacement'
                call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
                call POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_u0, sol_ptr, u, nnod_num, gathered_sizes, displacements)
                call EXPORT_SOLUTION(PolyMesh, u, IsPoly, 0)
            endif

            call COMPUTE_MODAL_COEFFICIENTS_GEN(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_v0, ic_velocity)
            ! SAVE SOLUTION
            if (IsSave_output .eqv. .true.) then
                if (mpi_id==0) print *, 'velocity'
                call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
                call POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_v0, sol_ptr, u, nnod_num, gathered_sizes, displacements)
                call EXPORT_SOLUTION(PolyMesh, u, IsPoly)


                open(newunit=unit_print, action='WRITE', file='v0_petsc_post.txt', &
                form='FORMATTED', status='replace')
                do i=1,Np
                    write(unit_print, *) u(i,:)
                enddo
                close(unit=unit_print)
            endif
            PetscCallA(KSPSolve(ksp3, petsc_u0, petsc_v_tmp, mpi_ierr))
            PetscCallA(VecCopy(petsc_v_tmp, petsc_u0, mpi_ierr))
            PetscCallA(KSPSolve(ksp3, petsc_v0, petsc_v_tmp, mpi_ierr))
            PetscCallA(VecCopy(petsc_v_tmp, petsc_v0, mpi_ierr))

            ! print to file
            ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_u0',viewer,mpi_ierr))
            ! PetscCallA(VecView(petsc_u0,viewer,mpi_ierr))
            ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
            ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_v0',viewer,mpi_ierr))
            ! PetscCallA(VecView(petsc_v0,viewer,mpi_ierr))
            ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

            if(mpi_id == 0) print *, ""
            if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", num_dt, " Time: ", t
            ! Mu_1 = (M-dt^2/2*A)u_0 + dt*M*v_0 + dt^2/2 * f_0

            ! petsc_v_tmp = -0.5 dt^2 A u0
            PetscCallA(MatMult(petsc_stiff, petsc_u0, petsc_v_tmp, mpi_ierr))
            PetscCallA(VecScale(petsc_v_tmp, -half_dt2, mpi_ierr))
            ! petsc_f = M u0
            PetscCallA(MatMult(petsc_mass, petsc_u0, petsc_f, mpi_ierr))
            ! petsc_f = (M-dt^2/2*A)u_0
            PetscCallA(VecAXPY(petsc_f, one, petsc_v_tmp, mpi_ierr))
            ! petsc_v_tmp = M v0
            PetscCallA(MatMult(petsc_mass, petsc_v0, petsc_v_tmp, mpi_ierr))
            ! petsc_f = (M-dt^2/2*A)u_0 + dt*M*v_0
            PetscCallA(VecAXPY(petsc_f, time_step, petsc_v_tmp, mpi_ierr))

            ! petsc_v_tmp = dt^2/2 * f_0(x)*f'_0(t)
            ! PetscCallA(VecCopy(petsc_rhs, petsc_v_tmp, mpi_ierr))
            ! PetscCallA(VecScale(petsc_v_tmp, half_dt2*time_function(time), mpi_ierr))

            !!! check if the same f'(t) applies for both forcing and BC
            ! petsc_f = [(M-dt^2/2*A)u_0 + dt*M*v_0] + dt^2/2 * f_0(x)*f'_0(t)
            PetscCallA(VecAXPY(petsc_f, half_dt2*time_function(t), petsc_rhs, mpi_ierr))
            ! M u1 = F
            PetscCallA(KSPSolve(ksp2, petsc_f, petsc_sol, mpi_ierr))

            ! SAVE SOLUTION
            ! if (IsSave_output .eqv. .true.) then
            !     call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            !     call POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_sol, sol_ptr, u, nnod_num, gathered_sizes, displacements)
            !     call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
            ! endif

            num_dt = num_dt + 1
            t = t + time_step

            ! loop
            ! Mu_{n+1) = F_n = (2M-dt^2*A)u_n - Mu_{n-1} + dt^2 * f_n
            do while (t <= stop_time)

                if(mpi_id == 0) print *, ""
                if(mpi_id == 0) write(*,'(A,I10,A,F8.5)') "Iteration: ", num_dt, " Time: ", t

                !!! petsc_u0 -> u_{n-1}
                !!! petsc_sol -> u_n
                !!! at the end write on petsc_sol for u_{n+1}

                ! petsc_v_tmp = dt^2 A u_n
                PetscCallA(MatMult(petsc_stiff, petsc_sol, petsc_v_tmp, mpi_ierr))
                PetscCallA(VecScale(petsc_v_tmp, -dt2, mpi_ierr))
                ! petsc_f = 2M u_n
                PetscCallA(MatMult(petsc_mass, petsc_sol, petsc_f, mpi_ierr))
                PetscCallA(VecScale(petsc_f, 2.0*one, mpi_ierr))
                ! petsc_f = (2M-dt^2*A)u_n
                PetscCallA(VecAXPY(petsc_f, one, petsc_v_tmp, mpi_ierr))
                ! petsc_v_tmp = M u_{n-1}
                PetscCallA(MatMult(petsc_mass, petsc_u0, petsc_v_tmp, mpi_ierr))
                ! petsc_f = (M-dt^2/2*A)u_n - M*u_{n-1}
                PetscCallA(VecAXPY(petsc_f, -one, petsc_v_tmp, mpi_ierr))

                ! Compute new RHS - keep initial RHS, muliply by time function
                !!! check if the same f'(t) applies for both forcing and BC

                ! petsc_rhs = dt^2 * f_n(x)*f'_n(t)
                ! F = petsc_f = [(M-dt^2/2*A)u_n + dt*M*u_{n-1}] + dt^2 * f_n(x)*f'_n(t)
                PetscCallA(VecAXPY(petsc_f, dt2*time_function(t), petsc_rhs, mpi_ierr))

                ! copy old solution before overwriting solution
                PetscCallA(VecCopy(petsc_sol, petsc_u0, mpi_ierr))

                ! solve linear system M u_{n+1} = F
                PetscCallA(KSPSolve(ksp2, petsc_f, petsc_sol, mpi_ierr))

                ! SAVE SOLUTION
                if ( (IsSave_output .eqv. .true.) .and. (mod(num_dt, num_dt_mon) == 0) ) then
                    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
                    call POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_sol, sol_ptr, u, nnod_num, gathered_sizes, displacements)
                    call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
                endif

                num_dt = num_dt + 1
                t = t + time_step

            end do

            if(mpi_id == 0) print *, "Solution end time: ", t

        end if

        ! final timestep save solution
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_sol, sol_ptr, u,  nnod_num, gathered_sizes, displacements)
        call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
    endif

    ! deallocate solution for post-processing
    deallocate(u)
    deallocate(nnod_num)

    if(mpi_np > 1) then
        deallocate(displacements)
        deallocate(gathered_sizes)
    endif

    !call PVD_SETUP(num_dt_mon, 0, num_dt)

    ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'petsc_sol',viewer,mpi_ierr))
    ! PetscCallA(VecView(petsc_sol,viewer,mpi_ierr))
    ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     STORE LOCAL NUMERATION TO RECONSTRUCT THE SOLUTION
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! allocate(nnod_num(local_dof))

    ! call CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     SCATTER PETSC SOLUTION AND STORE IN A FORTRAN ARRAY
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! allocate(u_loc(local_dof))
    ! allocate(u_glo(global_dof))

    ! if(mpi_np > 1) then

    !       allocate(gathered_sizes(mpi_np))

    !       call MPI_AllGather(local_dof, 1, MPI_INTEGER, gathered_sizes, 1, &
    !                   MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !       allocate(displacements(mpi_np))
    !       displacements(1) = 0
    !       do i = 2, mpi_np
    !             displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
    !       end do

    ! endif

    ! print *, 'SCATTER SOLUTION'
    ! PetscCallA(VecGetArrayF90(petsc_sol, sol_ptr, mpi_ierr))
    ! u_loc(1:local_dof) = sol_ptr
    ! if(mpi_np == 1) then
    !       u_glo = u_loc
    ! else
    !       call MPI_ALLGATHERV(u_loc, local_dof, MPI_DOUBLE_PRECISION, &
    !                   u_glo, gathered_sizes, displacements, MPI_DOUBLE_PRECISION, &
    !                   MPI_COMM_WORLD, mpi_ierr)
    ! endif
    ! deallocate(u_loc)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     RECONSTRUCT SOLUTION MATRIX FOR POST-PROCESSING
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! allocate(u(Np, DIM*PolyMesh%num_poly))
    ! u = RESHAPE(u_glo, (/Np, DIM*PolyMesh%num_poly /))
    ! deallocate(u_glo)

    ! print *,'Done with the solution'

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     COMPUTE MODAL SOLUTION
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !!! NOT WORKING WITH PROCESSES
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'---------------Compare with exact solution-------------'
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if (IS_MatrixFree .eqv. .true.) then
        print *, 'matrix-free exact solution'
    else
        print *, 'Computing modal coefficients...'
        call COMPUTE_MODAL_COEFFICIENTS(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_modal_coeff_uex)
        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'petsc_modal_coeff_uex',viewer,mpi_ierr))
        ! PetscCallA(VecView(petsc_modal_coeff_uex,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        print *, 'Calling solver for modal solution...'
        PetscCallA(KSPSolve(ksp3, petsc_modal_coeff_uex, petsc_uex, mpi_ierr))

        ! check time dependence
        if (IsTime_dependent .eqv. .true.) then
            !! WHAT TIME IS THE CORRECT?
            !t=t-time_step
            if (mpi_id == 0) print *, 'Scale modal solution by time function at t = ', t
            PetscCallA(VecScale(petsc_uex, time_function(t), mpi_ierr))
        endif

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_uex',viewer,mpi_ierr))
        ! PetscCallA(VecView(petsc_uex,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
    endif

    print *, 'Done with modal solutions'

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DESTROYING KSP SOLVERS AND VECTORS ON THE RIGHT HAND SIDE OF THE SYSTEMS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (IS_MatrixFree .eqv. .true.) then
        print *, 'matrix free destroy elements'
    else
        PetscCallA(KSPDestroy(ksp, mpi_ierr))
        PetscCallA(KSPDestroy(ksp2, mpi_ierr))
        PetscCallA(KSPDestroy(ksp3, mpi_ierr))
        PetscCallA(VecDestroy(petsc_rhs, mpi_ierr))
        PetscCallA(VecDestroy(petsc_modal_coeff_uex, mpi_ierr))
    endif

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     POST-PROCESSING: COMPUTING THE ERRORS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    err_L2 = 0.0
    err_DG = 0.0

    if (IS_MatrixFree .eqv. .true.) then
        print *, 'matrix free solver'
    else
        if(mpi_id == 0) print *,'Computing the errors...'

        call COMPUTE_ERROR_L2(petsc_mass_modal, petsc_sol, petsc_uex, global_dof, err_L2_mpi, local_dof)
        call COMPUTE_ERROR_DG(mat_dg, petsc_sol, petsc_uex, global_dof, err_DG_mpi, local_dof)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    call MPI_REDUCE(err_L2_mpi, err_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, MPI_COMM_WORLD, mpi_ierr)

    call MPI_REDUCE(err_DG_mpi, err_DG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id == 0) then
        err_L2 = sqrt(err_L2)
        err_DG = sqrt(err_DG)
        print *,'Done with the errors'
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id == 0) then
        print *,'ERROR IN NORM L2: ', err_L2
        print *,'ERROR IN NORM DG: ', err_DG
    endif

    hmax_mpi = compute_hmax(PolyMesh)

    call MPI_REDUCE(hmax_mpi, hmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                    0, MPI_COMM_WORLD, mpi_ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id == 0) print *, 'GRID SIZE: ', hmax
    if (mpi_id == 0) call WRITE_ERRORS(p, err_DG, err_L2, hmax, PolyMesh, IsPoly)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEALLOCATING PETSC MATRICES AND VECTORS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (IS_MatrixFree .eqv. .true.) then
        print *, 'matrix free destroy vector and matrices'

        deallocate(M_loc, M_modal_loc, u0_loc, un_loc, usol_loc, rhs_loc)

        do i=1,PolyMesh%num_elem_loc
            !deallocate(K_loc(i)%values, A_dg_loc(i)%values, internal_neigh(i)%values)
        end do

        !deallocate(K_loc, A_dg_loc, internal_neigh)

    else
        PetscCallA(MatDestroy(petsc_stiff, mpi_ierr))
        PetscCallA(MatDestroy(petsc_mass, mpi_ierr))
        PetscCallA(MatDestroy(mat_dg, mpi_ierr))
        PetscCallA(VecDestroy(petsc_sol,mpi_ierr))
        PetscCallA(VecDestroy(petsc_uex,mpi_ierr))
        !!! DEALLOCATE NEW OBJECTS
        PetscCallA(MatDestroy(petsc_mass_modal, mpi_ierr))
        PetscCallA(VecDestroy(petsc_f, mpi_ierr))
        PetscCallA(VecDestroy(petsc_v_tmp,mpi_ierr))
        PetscCallA(VecDestroy(petsc_u0,mpi_ierr))
        PetscCallA(VecDestroy(petsc_v0,mpi_ierr))
    endif

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     POST-PROCESSING: EXPORTING THE SOLUTION
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! call EXPORT_SOLUTION(PolyMesh, u, IsPoly, mpi_id, num_dt)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    END SETUP
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    finish = MPI_WTIME()
    call calc_time(time_hour, time_min, time_sec, int(finish-start))

    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'-------------------------------------------------------'
        write(*,'(A,I2,A,I2,A,I2,A)') &
                'Set-up time = ', time_hour,' h ' , time_min,' m ' , time_sec,' s'
        write(*,'(A)')'-------------------------------------------------------'
        write(*,'(A)')
    endif

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     FINALIZE MPI AND PETSC
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call PetscFinalize(mpi_ierr)
    call MPI_FINALIZE(mpi_ierr)

end program Lymph3D