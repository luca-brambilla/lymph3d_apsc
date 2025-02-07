! This file is part of the library LYMPH3D

!> @brief Lymph3D (Discontinuous Galerkin methods on polyhedral meshes for PDE problems)

! Here starts the code Lymph3D

program Lymph3D

    use mpi
    use Poly_setup_mpi
    use problem_data_and_properties
    use Poly_global
    use Poly_data
    use Poly_mesh
    use solution_processing
    use matrix_free
    use exchange_data
    use global_parameters
    use utilities
    use checks

    use find_tet

    implicit none

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEFINITION OF  VARIABLES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    logical :: IsPoly
    integer(kind=4) :: local_dof, global_dof
    integer(kind=4) :: Np, p, Npoly
    integer(kind=4), dimension(:), allocatable :: nnod_num
    real(kind=8) :: err_L2
    real(kind=8) :: err_DG
    real(kind=8) :: hmax, hmax_mpi
    real(kind=8), dimension(:,:), allocatable :: u
    integer(kind=4), dimension(:), allocatable :: gathered_sizes, displacements

    integer(kind=4) :: i, mat_id, ie_loc

    type(Data_Structure) :: PolyData
    type(Mesh_Structure) :: PolyMesh

    integer(kind=4) :: num_dt = 1     ! number of iteration
    real(kind=8) :: t             ! time variable

    real(kind=8), parameter :: one = 1.0

    ! each local has a matrix

    real(kind=8), dimension(:,:), allocatable :: v0_loc
    real(kind=8), dimension(:,:), allocatable :: u0_loc
    real(kind=8), dimension(:,:), allocatable :: un_loc
    real(kind=8), dimension(:,:), allocatable :: rhs_stat_loc
    real(kind=8), dimension(:,:), allocatable :: rhs_dyn_loc

    real(kind=8), dimension(:,:), allocatable :: un_mpi, uex_mpi

    real(kind=8), dimension(:,:,:,:), allocatable :: M_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: M_modal_loc
    !real(kind=8), dimension(:,:,:), allocatable :: D_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: K_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: A_dg_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: R_M_loc
    real(kind=8), dimension(:,:,:,:), allocatable :: R_M_modal_loc
    real(kind=8), dimension(:,:), pointer :: M_tmp

    integer(kind=4) :: n_neigh

    type(ScatteredArray), dimension(:,:), allocatable :: send_data, recv_data

    real(kind=8) :: t1, t2

    integer(kind=4) :: num_inter_loc!, col, col_mpi, istart, iestart

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    allocate(mpi_stat(MPI_STATUS_SIZE))

    call INITIALIZATION()

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

    ! compute h_max
    hmax_mpi = compute_hmax(PolyMesh)
    call MPI_REDUCE(hmax_mpi, hmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                    0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if (mpi_id == 0) print *, 'GRID SIZE: ', hmax
    if (mpi_id == 0) print *, 'TIME STEP: ', time_step

    !! CHECK IF IT WORKS
    IS_failCFL = .true.
    call COMPUTE_CFL(hmax, time_step)

    call LYMPH3D_BARRIER

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    !call CHECK_MPI_EXCHANGE(PolyMesh, Np)

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

    ! E = 32e9
    ! nu = 0.2

    ! set the value of rho, lambda and mu
    do mat_id=1,PolyData%nmat
        ! PolyData%prop_mat(mat_id,1) = 2400
        ! PolyData%prop_mat(mat_id,2) = E*nu / ((1+nu)*(1-2*nu))
        ! PolyData%prop_mat(mat_id,3) = E / (2 * (1+nu))
        PolyData%prop_mat(mat_id,1) = 10
        PolyData%prop_mat(mat_id,2) = 1
        PolyData%prop_mat(mat_id,3) = 1
    enddo

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     SET VECTORS AND MATRICES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'--------------------Compute solution-------------------'
    endif
    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    print *, 'matrix-free - set matrices and vectors'

    n_neigh = PolyMesh%Elem_loc(1)%num_faces

    allocate( K_loc(PolyMesh%num_elem_loc, n_neigh+1, DIM*Np, DIM*Np) )
    allocate( A_dg_loc(PolyMesh%num_elem_loc, n_neigh+1, DIM*Np, DIM*Np) )

    !allocate( D_loc(PolyMesh%num_elem_loc, DIM*Np, DIM*Np) )

    allocate( M_loc(PolyMesh%num_elem_loc, DIM, Np, Np) )
    allocate( M_modal_loc(PolyMesh%num_elem_loc, DIM, Np, Np) )
    allocate( rhs_stat_loc(PolyMesh%num_elem_loc, DIM*Np) )
    allocate( rhs_dyn_loc(PolyMesh%num_elem_loc, DIM*Np) )

    allocate( R_M_loc(PolyMesh%num_elem_loc, DIM, Np, Np) )
    allocate( R_M_modal_loc(PolyMesh%num_elem_loc, DIM, Np, Np) )

    allocate( M_tmp(Np, Np) )

    ! initial conditions
    allocate(u0_loc(PolyMesh%num_poly_loc, DIM*Np))
    allocate(un_loc(PolyMesh%num_poly_loc, DIM*Np))
    allocate(v0_loc(PolyMesh%num_poly_loc, DIM*Np))

    num_inter_loc = PolyMesh%num_elem_inter_vec(mpi_id+1)

    allocate(un_mpi(num_inter_loc, DIM*Np))
    allocate(uex_mpi(num_inter_loc, DIM*Np))

    u0_loc = 0.0d0
    un_loc = 0.0d0
    v0_loc = 0.0d0
    un_mpi = 0.0d0
    uex_mpi= 0.0d0

    if (mpi_np>1) call MPI_EXCHANGE_ALLOCATE(PolyMesh, send_data, recv_data)


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ASSEMBLE MATRICES
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    print *, 'assemble matrix-free matrices'

    ! make all matrices
    call MAKE_MATRICES_FREE(PolyMesh, PolyData, PolyMesh%num_elem_loc, Np, K_loc, A_dg_loc, M_loc, M_modal_loc, n_neigh)

    do ie_loc=1,PolyMesh%num_elem_loc
        do i=1,DIM
            M_tmp = M_loc(ie_loc,i,:,:)
            R_M_loc(ie_loc,i,:,:) = cholesky(M_tmp,Np)
            M_tmp = M_modal_loc(ie_loc,i,:,:)
            R_M_modal_loc(ie_loc,i,:,:) = cholesky(M_tmp,Np)
        enddo
    enddo

    deallocate(M_tmp)


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ASSEMBLE RIGHT HAND SIDE
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    print *, 'assemble matrix-free RHS'
    call MAKE_RHS_FREE(PolyMesh, PolyData, PolyMesh%num_elem_loc, Np, rhs_stat_loc, f_null, gd_null, gn_null)

    call MAKE_RHS_FREE(PolyMesh, PolyData, PolyMesh%num_elem_loc, Np, rhs_dyn_loc, f_time, gd, gn)

    rhs_dyn_loc = 0.0d0

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CALLING SOLVER
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    IsSave_output = .true.

    t = 0.0d0
    num_dt = 0

    !stop_time = SQRT2 / 4.0 + 1.0 * SQRT2
    time_step = 1.0d-3
    num_dt_mon = 20

    !stop_time = SQRT2 / 4.0 + 9.0 * SQRT2 ! 10 peaks
    stop_time = SQRT2 / 4.0 + 1.0 * SQRT2
    !stop_time = 0.101
    !stop_time = 0.001
    
    dt2 = time_step*time_step
    half_dt2 = 0.5d0*dt2

    if(mpi_id == 0) print *, "                   TIME LOOP START                   "

    ! allocate nnod_num, gathered_sizes, displacements, u
    call PREPROCESS_SOLUTION_MATRIX_FREE(PolyMesh, local_dof, nnod_num, gathered_sizes, displacements, u)


    if(mpi_id == 0) write(*,'(A,I10,A,F14.5)') "Iteration: ", 0, " Time: ", t


    if (mpi_id==0) print *, "Assemble initial conditions"
    call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, ic_displacement, u0_loc)
    ! call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, ic_velocity, v0_loc)
    call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, ic_displacement, v0_loc)

    ! SAVE IC
    if (IsSave_output .eqv. .true.) then
        t1 = MPI_WTIME()
        if (mpi_id==0) print *, '--- IC displacement ---'
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call POST_PROCESS_MATRIX_FREE(PolyMesh, u0_loc, u, gathered_sizes, displacements)
        call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)

        if (mpi_id==0) print *, '---   IC velocity   ---'
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call POST_PROCESS_MATRIX_FREE(PolyMesh, v0_loc, u, gathered_sizes, displacements)
        call EXPORT_SOLUTION(PolyMesh, u, IsPoly)

        t2 = MPI_WTIME()
        tp_export = tp_export + t2 - t1
    endif
    
    num_dt = num_dt + 1

    ! u_1 = M^-1(dt^2/2 * f_0 - dt^2/2*A*u_0) + u0 + dt*v_0

    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if(mpi_id == 0) print *, ""
    if(mpi_id == 0) write(*,'(A,I10,A,F14.5)') "Iteration: ", num_dt, " Time: ", t+time_step
    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! update interface solution for u^{0}
    if (mpi_np>1) then
        call MPI_EXCHANGE_DOF(PolyMesh, u0_loc, un_mpi, send_data, recv_data)
    endif

    call LYMPH3D_BARRIER

    ! compute u^{1}
    call FIRST_TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, R_M_loc, rhs_stat_loc, rhs_dyn_loc, u0_loc, v0_loc, un_loc, un_mpi)

    ! SAVE FIRST ITERATION
    ! if (IsSave_output .eqv. .true.) then
    !     t1 = MPI_WTIME()
    !     call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
    !     call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
    !     t2 = MPI_WTIME()
    !     tp_export = tp_export + t2 - t1
    ! endif

    ! update time
    num_dt = num_dt + 1
    t = t + time_step

    ! update interface solution for u^{1}
    if (mpi_np>1) then
        call MPI_EXCHANGE_DOF(PolyMesh, un_loc, un_mpi, send_data, recv_data)
    endif

    ! loop start
    do while (t <= stop_time)
        if(mpi_id == 0) write(*,'(A,I10,A,F14.6)') "Iteration: ", num_dt, " Time: ", t+time_step
        
        ! v0_loc is used for u^{n+1}
        call TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, R_M_loc, rhs_stat_loc, rhs_dyn_loc, u0_loc, un_loc, v0_loc, un_mpi)

        ! update solution
        u0_loc = un_loc
        un_loc = v0_loc

        ! SAVE SOLUTION
        if ( (IsSave_output .eqv. .true.) .and. (mod(num_dt, num_dt_mon) == 0) ) then
            t1 = MPI_WTIME()
            call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
            call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
            t2 = MPI_WTIME()
            tp_export = tp_export + t2 - t1
        endif

        ! update interface solution for u^{n}
        if (mpi_np>1) then
            call MPI_EXCHANGE_DOF(PolyMesh, un_loc, un_mpi, send_data, recv_data)
        endif

        ! update time
        num_dt = num_dt + 1
        t = t + time_step

        call LYMPH3D_BARRIER
    end do

    num_dt = num_dt - 1

    ! LAST ITERATION
    ! if (IsSave_output .eqv. .true.) then
    !     t1 = MPI_WTIME()
    !     call POST_PROCESS_MATRIX_FREE(PolyMesh, un_loc, u, gathered_sizes, displacements)
    !     call EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
    !     t2 = MPI_WTIME()
    !     tp_export = tp_export + t2 - t1
    ! endif

    ! deallocate solution for post-processing
    deallocate(u)
    deallocate(nnod_num)

    if(mpi_np > 1) then
        deallocate(displacements)
        deallocate(gathered_sizes)
    endif

    if (mpi_id==0) print *, 'Writing Ensight .case file...'
    if (mpi_id==0) call ENSIGHT_WRITE_CASE('MONITORS/', num_dt, 'DISPLACEMENT')

    !call PVD_SETUP(num_dt_mon, 0, num_dt)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     COMPUTE MODAL SOLUTION
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'---------------Compare with exact solution-------------'
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    t1 = MPI_WTIME()
    ! u0_loc for exact solution - time dependent component
    call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, uex, u0_loc)
    ! v0_loc for exact solution - static component
    call COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, uex_stat, v0_loc)

    ! check time dependence
    if (IsTime_dependent .eqv. .true.) then
        if (mpi_id == 0) write(*, '(A,F14.5)') 'Scale modal solution by time function at t = ', t
        u0_loc = v0_loc*0.0d0 + u0_loc * time_function(t)
    endif
    t2 = MPI_WTIME()
    tp_exact = t2-t1

    if (mpi_np>1) then
        call MPI_EXCHANGE_DOF(PolyMesh, u0_loc, uex_mpi, send_data, recv_data)
    endif

    print *, 'Done with modal solutions'

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     POST-PROCESSING: COMPUTING THE ERRORS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    err_L2 = 0.0d0
    err_DG = 0.0d0

    t1 = MPI_WTIME()
        print *, 'Computing the errors...'
        ! u0_loc for exact solution
        call COMPUTE_ERROR_L2_MATRIX_FREE(PolyMesh, Np, M_modal_loc, un_loc, u0_loc, err_L2)
        call COMPUTE_ERROR_DG_MATRIX_FREE(PolyMesh, Np, A_dg_loc, un_loc, u0_loc, un_mpi, uex_mpi, err_DG)

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, err_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, err_DG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, MPI_COMM_WORLD, mpi_ierr)
    else
        call MPI_REDUCE(err_L2, err_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(err_DG, err_DG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, MPI_COMM_WORLD, mpi_ierr)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! print final errors
    if (mpi_id == 0) then
        err_L2 = dsqrt(err_L2)
        err_DG = dsqrt(err_DG)
        print *, 'Done with the errors'

        print *, 'GRID SIZE: ', hmax
        print *, 'ERROR IN NORM L2: ', err_L2
        print *, 'ERROR IN NORM DG: ', err_DG
    endif

    t2 = MPI_WTIME()
    tp_error = t2 - t1

    !if (mpi_id == 0) call WRITE_ERRORS(p, err_DG, err_L2, hmax, PolyMesh, IsPoly)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEALLOCATING MATRICES AND VECTORS
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    print *, 'Destroy vector and matrices'

    deallocate(u0_loc, un_loc, v0_loc, rhs_stat_loc, rhs_dyn_loc)

    do i=1,PolyMesh%num_elem_loc
        !deallocate(K_loc(i)%values, A_dg_loc(i)%values, internal_neigh(i)%values)
    end do

    deallocate(K_loc, A_dg_loc)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    END SETUP
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call PRINT_PROFILING

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
!     FINALIZE MPI
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call MPI_FINALIZE(mpi_ierr)

end program Lymph3D