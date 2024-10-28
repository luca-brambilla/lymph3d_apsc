!    Author: Lorenzo Gorlezza
!    This file is part of the library LYMPH3D

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
      use post_processing
      use SET_PETSC_SYSTEM

      use MOD_MPI_CUSTOM

      implicit none

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEFINITION OF PETSC VARIABLES
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      Mat :: petsc_stiff, petsc_mass, mat_dg, petsc_mass_modal ! stiffness, mass and DG matrices, reconstruction matrix
      Vec :: petsc_sol, petsc_rhs ! solution and rhs vectors for the algebraic system
      Vec :: petsc_uex, petsc_modal_coeff_uex ! modal (exact) solution and its coefficients

      Vec :: petsc_tmpv, petsc_f ! temp vectors to store solution for sums
      Vec :: petsc_u0, petsc_v0 ! initial condition read from data

      PetscViewer viewer ! abstract PETSc object for displaying PETSc objects and their data

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
      real(kind=8), dimension(:), allocatable :: u_loc, u_glo
      real(kind=8), dimension(:,:), allocatable :: u
      integer(kind=4), dimension(:), allocatable :: gathered_sizes, displacements

      real(kind=8) :: E, nu
      integer(kind=4) :: i, mat_id

      type(Data_Structure) :: PolyData
      type(Mesh_Structure) :: PolyMesh 

      integer(kind=4) :: num_dt = 1     ! number of iteration
      real(kind=8) :: t             ! time variable
      
      real(kind=8) dt2 ! dt^2
      real(kind=8) half_dt2 ! dt^2 / 2
      real(kind=8), parameter :: one = 1.0
      ! IsTime_dependent = .false.

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      
      allocate(mpi_stat(MPI_STATUS_SIZE))

      call INITIALIZATION()

      call MPI_OP_CREATE(MPI_ZERO_OVERWRITE, .TRUE., MPI_ZERO_OVERWRITE_OP, mpi_user_reduction_error)
      call MPI_OP_CREATE(MPI_OVERWRITE_BY_NEW, .FALSE., MPI_OVERWRITE_BY_NEW_OP, mpi_user_reduction_error)

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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     READ INPUT FILES AND ALLOCATE VARIABLES
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call READ_INPUT_FILES(PolyData,PolyMesh)
      if(mpi_id == 0) then
            write(*,'(A)') 
            write(*,'(A)')'--------------------Type of problem--------------------'
            write(*,'(A)')   
            if (IsTime_dependent .eqv. .false.) then
                  write(*,'(A)')'Stationary problem'
                  write(*,'(A,L)')'variable: ',IsTime_dependent
            else
                  write(*,'(A)')'Time dependent problem'
                  write(*,'(A,L)')'variable      : ',IsTime_dependent
                  write(*,'(A,E8.3)')'start time [s] : ', start_time
                  write(*,'(A,E8.3)')'end time   [s] : ', stop_time
                  write(*,'(A,E8.3)')'timestep   [s] : ', time_step
            endif
            write(*,'(A)')
      endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CHOOSE TO SOLVE WITH TETRAHEDRA
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IsPoly = .false.
      Npoly = PolyMesh%num_poly

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CHOOSE TO SOLVE WITH POLYHEDRA
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      ! set a number of polyhedra which may be different from the one read from the mesh file
      ! Npoly=900

      if (Npoly /= PolyMesh%num_elem) IsPoly = .true.
     
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     PARTITION OF THE GRID AND GENERATION OF LOCAL CONNECTIVITY 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call MAKE_PARTITION_AND_MPI_FILES(PolyData, PolyMesh, Npoly)
      write(*,'(A28,I5,A28,I5)') 'Global number of polyhedra:', PolyMesh%num_poly, &
            'Local number of polyhedra:', PolyMesh%num_poly_loc
      print *, 'Done reading mesh'

      Np = PolyMesh%Elem_loc(1)%NDof_loc ! number of degrees of freedom of each element
      p = PolyMesh%Elem_loc(1)%Degree ! order of basis functions
      
      global_dof = 3*PolyMesh%num_poly*Np;
      local_dof = 3*PolyMesh%num_poly_loc*Np;
      PRINT *, 'Local number of degrees of freedom: ', local_dof

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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     SET PETSC VECTORS AND MATRICES 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !!! NOT WORKING PRINT WITH PROCESSES
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
      if (mpi_id .eq. 0) then
            write(*,'(A)') 
            write(*,'(A)')'--------------------Compute solution-------------------'
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      print *, 'SET PETSC MATRICES AND VECTORS'

      call SET_PETSC_MATRIX(petsc_stiff, local_dof, global_dof)
      call SET_PETSC_MATRIX(petsc_mass, local_dof, global_dof)
      call SET_PETSC_MATRIX(mat_dg, local_dof, global_dof)
      call SET_PETSC_MATRIX(petsc_mass_modal, local_dof, global_dof)

      call SET_PETSC_VECTOR(petsc_rhs, local_dof, global_dof)
      call SET_PETSC_VECTOR(petsc_sol, local_dof, global_dof)
      call SET_PETSC_VECTOR(petsc_uex, local_dof, global_dof)
      
      call SET_PETSC_VECTOR(petsc_tmpv, local_dof, global_dof)
      call SET_PETSC_VECTOR(petsc_f, local_dof, global_dof)
      call SET_PETSC_VECTOR(petsc_u0, local_dof, global_dof)
      call SET_PETSC_VECTOR(petsc_v0, local_dof, global_dof)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     ASSEMBLE MATRICES 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  

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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>       
!     ASSEMBLE RIGHT HAND SIDE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      print *, 'ASSEMBLE RHS'
      
      call MAKE_RHS(PolyMesh, PolyData, petsc_num, global_dof, Np, petsc_rhs)

      ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_rhs',viewer,mpi_ierr))
      ! PetscCallA(VecView(petsc_rhs,viewer,mpi_ierr))
      ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
!     SETTING SOLVERS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      print *, 'SETTING SOLVERS'

      call SOLVER_SETTINGS(petsc_stiff, ksp, pc)
      call SOLVER_SETTINGS(petsc_mass, ksp2, pc2)
      call SOLVER_SETTINGS(petsc_mass_modal, ksp3, pc3)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     CALLING SOLVER
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
      print *, 'CALLING SOLVER'
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
      
      if (IsTime_dependent .eqv. .false.) then
            ! Au=f
            PetscCallA(KSPSolve(ksp, petsc_rhs, petsc_sol, mpi_ierr))

      else
            if(mpi_id == 0) print *, "                   TIME LOOP START                   "
            ! initialization for time problem
            t = 0.0
            num_dt = 1 !!! compute number of iteration if start not t=0?
            stop_time = SQRT2/4.0 + 9*SQRT2 ! 10 peaks

            if(mpi_id == 0) print *, "Iteration: ", 0, " Time: ", t
            print *, "Assemble initial conditions"

            !! MODAL OR SOLUTION???
            ! vectors initial conditions

            ! integral
            ! call MAKE_VECTOR(PolyMesh, petsc_num, global_dof, Np, petsc_u0, ic_displacement)
            ! call MAKE_VECTOR(PolyMesh, petsc_num, global_dof, Np, petsc_v0, ic_velocity)
            ! PRINT ?
            ! compute modal coefficients
            ! PetscCallA(KSPSolve(ksp3, petsc_u0, petsc_tmpv, mpi_ierr))
            ! PetscCallA(VecCopy(petsc_tmpv, petsc_u0, mpi_ierr))
            ! PetscCallA(KSPSolve(ksp3, petsc_v0, petsc_tmpv, mpi_ierr))
            ! PetscCallA(VecCopy(petsc_tmpv, petsc_v0, mpi_ierr))

            ! modal coefficients of initial conditions 
            call COMPUTE_MODAL_COEFFICIENTS_GEN(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_u0, ic_displacement)
            call COMPUTE_MODAL_COEFFICIENTS_GEN(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_v0, ic_velocity)
            PetscCallA(KSPSolve(ksp3, petsc_u0, petsc_tmpv, mpi_ierr))
            PetscCallA(VecCopy(petsc_tmpv, petsc_u0, mpi_ierr))
            PetscCallA(KSPSolve(ksp3, petsc_v0, petsc_tmpv, mpi_ierr))
            PetscCallA(VecCopy(petsc_tmpv, petsc_v0, mpi_ierr))
            
            ! print to file
            PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_u0',viewer,mpi_ierr))
            PetscCallA(VecView(petsc_u0,viewer,mpi_ierr))
            PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
            PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_v0',viewer,mpi_ierr))
            PetscCallA(VecView(petsc_v0,viewer,mpi_ierr))
            PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
            
            ! compute solution with modal coeff
            ! PetscCallA(MatMult(petsc_mass_modal, petsc_u0, petsc_tmpv, mpi_ierr))
            ! PetscCallA(VecCopy(petsc_tmpv, petsc_u0, mpi_ierr))
            ! PetscCallA(MatMult(petsc_mass_modal, petsc_v0, petsc_tmpv, mpi_ierr))
            ! PetscCallA(VecCopy(petsc_tmpv, petsc_v0, mpi_ierr))

            ! print to file
            ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_u0_mod',viewer,mpi_ierr))
            ! PetscCallA(VecView(petsc_u0,viewer,mpi_ierr))
            ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
            ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_v0_mod',viewer,mpi_ierr))
            ! PetscCallA(VecView(petsc_v0,viewer,mpi_ierr))
            ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

            ! petsc_u0 = 
            ! petsc_v0 = 
            if(mpi_id == 0) print *, ""
            if(mpi_id == 0) print *, "Iteration: ", num_dt, " Time: ", t

            ! Mu_1 = (M-dt^2/2*A)u_0 + dt*M*v_0 + dt^2/2 * f_0

            dt2 = time_step * time_step
            half_dt2 = 0.5 * dt2
            
            ! petsc_tmpv = -0.5 dt^2 A u0
            PetscCallA(MatMult(petsc_stiff, petsc_u0, petsc_tmpv, mpi_ierr))
            PetscCallA(VecScale(petsc_tmpv, -half_dt2, mpi_ierr))
            ! petsc_f = M u0
            PetscCallA(MatMult(petsc_mass, petsc_u0, petsc_f, mpi_ierr))
            ! petsc_f = (M-dt^2/2*A)u_0
            PetscCallA(VecAXPY(petsc_f, one, petsc_tmpv, mpi_ierr))
            ! petsc_tmpv = M v0
            PetscCallA(MatMult(petsc_mass, petsc_v0, petsc_tmpv, mpi_ierr))
            ! petsc_f = (M-dt^2/2*A)u_0 + dt*M*v_0
            PetscCallA(VecAXPY(petsc_f, time_step, petsc_tmpv, mpi_ierr))

            ! petsc_tmpv = dt^2/2 * f_0(x)*f'_0(t)
            ! PetscCallA(VecCopy(petsc_rhs, petsc_tmpv, mpi_ierr))
            ! PetscCallA(VecScale(petsc_tmpv, half_dt2*time_function(time), mpi_ierr))

            !!! check if the same f'(t) applies for both forcing and BC
            ! petsc_f = [(M-dt^2/2*A)u_0 + dt*M*v_0] + dt^2/2 * f_0(x)*f'_0(t)
            PetscCallA(VecAXPY(petsc_f, half_dt2*time_function(t), petsc_rhs, mpi_ierr))
            ! M u1 = F
            PetscCallA(KSPSolve(ksp2, petsc_f, petsc_sol, mpi_ierr))

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            !!! SAVE SOLUTION
            ! ...

            ! copy old solution
            PetscCallA(VecCopy(petsc_sol, petsc_u0, mpi_ierr))

            num_dt = num_dt + 1
            t = t + time_step

            ! loop
            ! Mu_{n+1) = F_n = (2M-dt^2*A)u_n - Mu_{n-1} + dt^2 * f_n
            do while (t <= stop_time)

                  if(mpi_id == 0) print *, ""
                  if(mpi_id == 0) print *, "Iteration: ", num_dt, " Time: ", t

                  ! petsc_tmpv = dt^2 A u_n
                  PetscCallA(MatMult(petsc_stiff, petsc_sol, petsc_tmpv, mpi_ierr))
                  PetscCallA(VecScale(petsc_tmpv, -dt2, mpi_ierr))
                  ! petsc_f = 2M u_n
                  PetscCallA(MatMult(petsc_mass, petsc_sol, petsc_f, mpi_ierr))
                  PetscCallA(VecScale(petsc_f, 2.0*one, mpi_ierr))
                  ! petsc_f = (2M-dt^2*A)u_n
                  PetscCallA(VecAXPY(petsc_f, one, petsc_tmpv, mpi_ierr))
                  ! petsc_tmpv = M u_{n-1}
                  PetscCallA(MatMult(petsc_mass, petsc_u0, petsc_tmpv, mpi_ierr))
                  ! petsc_f = (M-dt^2/2*A)u_n - M*u_{n-1}
                  PetscCallA(VecAXPY(petsc_f, -one, petsc_tmpv, mpi_ierr))

                  ! Compute new RHS - keep initial RHS, muliply by time function
                  !!! check if the same f'(t) applies for both forcing and BC

                  ! petsc_tmpv = dt^2 * f_n(x)*f'_n(t)
                  ! PetscCallA(VecCopy(petsc_rhs, petsc_tmpv, mpi_ierr))
                  ! PetscCallA(VecScale(petsc_tmpv, time_function(t)*dt2, mpi_ierr))

                  ! F = petsc_f = [(M-dt^2/2*A)u_n + dt*M*u_{n-1}] + dt^2 * f_n(x)*f'_n(t)
                  PetscCallA(VecAXPY(petsc_f, dt2*time_function(t), petsc_rhs, mpi_ierr))
                  ! solve linear system M u1 = F
                  PetscCallA(KSPSolve(ksp2, petsc_f, petsc_sol, mpi_ierr))

                  call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
                  !!! SAVE SOLUTION
                  ! ...

                  ! copy old solution
                  PetscCallA(VecCopy(petsc_sol, petsc_u0, mpi_ierr))

                  num_dt = num_dt + 1
                  t = t + time_step
            end do

      end if

      ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'petsc_sol',viewer,mpi_ierr))
      ! PetscCallA(VecView(petsc_sol,viewer,mpi_ierr))
      ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     STORE LOCAL NUMERATION TO RECONSTRUCT THE SOLUTION
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>       

      allocate(nnod_num(local_dof))

      call CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     SCATTER PETSC SOLUTION AND STORE IN A FORTRAN ARRAY
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(u_loc(local_dof))
      allocate(u_glo(global_dof))

      if(mpi_np > 1) then

            allocate(gathered_sizes(mpi_np))
      
            call MPI_AllGather(local_dof, 1, MPI_INTEGER, gathered_sizes, 1, & 
                        MPI_INTEGER, MPI_COMM_WORLD, ierr)
            
            allocate(displacements(mpi_np))
            displacements(1) = 0
            do i = 2, mpi_np
                  displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
            end do

      endif

      print *, 'SCATTER SOLUTION'
      PetscCallA(VecGetArrayF90(petsc_sol, sol_ptr, mpi_ierr))
      u_loc(1:local_dof) = sol_ptr
      if(mpi_np == 1) then
            u_glo = u_loc
      else
            call MPI_ALLGATHERV(u_loc, local_dof, MPI_DOUBLE_PRECISION, &
                        u_glo, gathered_sizes, displacements, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, mpi_ierr)
      endif
      deallocate(u_loc)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     RECONSTRUCT SOLUTION MATRIX FOR POST-PROCESSING
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
      allocate(u(Np, 3*PolyMesh%num_poly))    
      u = RESHAPE(u_glo, (/Np, 3*PolyMesh%num_poly /))
      deallocate(u_glo)

      print *,'Done with the solution'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     COMPUTE MODAL SOLUTION
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !!! NOT WORKING WITH PROCESSES
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
      if (mpi_id == 0) then
            write(*,'(A)')
            write(*,'(A)')'---------------Compare with exact solution-------------'
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      print *, 'Computing modal coefficients...'
      call COMPUTE_MODAL_COEFFICIENTS(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_modal_coeff_uex)
      ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'petsc_modal_coeff_uex',viewer,mpi_ierr))
      ! PetscCallA(VecView(petsc_modal_coeff_uex,viewer,mpi_ierr))
      ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
      
      print *, 'Calling solver for modal solution...'
      PetscCallA(KSPSolve(ksp3, petsc_modal_coeff_uex, petsc_uex, mpi_ierr))

      ! check time dependence
      if (IsTime_dependent .eqv. .true.) PetscCallA(VecScale(petsc_uex, 1.0/time_function(t), mpi_ierr))

      PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'vec_uex',viewer,mpi_ierr))
      PetscCallA(VecView(petsc_uex,viewer,mpi_ierr))
      PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
      
      print *, 'Done with modal solutions'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     DESTROYING KSP SOLVERS AND VECTORS ON THE RIGHT HAND SIDE OF THE SYSTEMS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  

      PetscCallA(KSPDestroy(ksp, mpi_ierr))
      PetscCallA(KSPDestroy(ksp2, mpi_ierr))
      PetscCallA(KSPDestroy(ksp3, mpi_ierr))
      PetscCallA(VecDestroy(petsc_rhs, mpi_ierr))
      PetscCallA(VecDestroy(petsc_modal_coeff_uex, mpi_ierr))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     POST-PROCESSING: COMPUTING THE ERRORS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      err_L2 = 0.0
      err_DG = 0.0
      
      if(mpi_id == 0) print *,'Computing the errors...'
      
      !! WHICH MASS MATRIX?
      call COMPUTE_ERROR_L2(petsc_mass_modal, petsc_sol, petsc_uex, global_dof, err_L2_mpi, local_dof)
      call COMPUTE_ERROR_DG(mat_dg, petsc_sol, petsc_uex, global_dof, err_DG_mpi, local_dof)

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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     DEALLOCATING PETSC MATRICES AND VECTORS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      PetscCallA(MatDestroy(petsc_stiff, mpi_ierr))
      PetscCallA(MatDestroy(petsc_mass, mpi_ierr))
      PetscCallA(MatDestroy(mat_dg, mpi_ierr))
      PetscCallA(VecDestroy(petsc_sol,mpi_ierr))
      PetscCallA(VecDestroy(petsc_uex,mpi_ierr))
      !!! DEALLOCATE NEW OBJECTS
      PetscCallA(MatDestroy(petsc_mass_modal, mpi_ierr))
      PetscCallA(VecDestroy(petsc_f, mpi_ierr))
      PetscCallA(VecDestroy(petsc_tmpv,mpi_ierr))
      PetscCallA(VecDestroy(petsc_u0,mpi_ierr))
      PetscCallA(VecDestroy(petsc_v0,mpi_ierr))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     POST-PROCESSING: EXPORTING THE SOLUTION
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   

      call EXPORT_SOLUTION(PolyMesh, u, IsPoly, mpi_id, num_dt)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    END SETUP 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       
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
 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!     FINALIZE MPI AND PETSC
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call PetscFinalize(mpi_ierr)
      call MPI_FINALIZE(mpi_ierr)

      end program Lymph3D