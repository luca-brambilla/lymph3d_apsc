module solution_processing

#include<petsc/finclude/petscksp.h>

    use petscksp
    use Poly_setup_mpi
    use Poly_global
    use problem_data_and_properties
    use Poly_mesh
    use Poly_ref_mappings
    use basis_function
    use export_file_formats
    use SET_PETSC_SYSTEM
    use mesh_partition_and_mpi_files
    
    use global_parameters

    implicit none

    contains

!> allocate global solution and prepare size and displacement vectors for MPI
subroutine PREPROCESS_SOLUTION(PolyMesh, local_dof, nnod_num, gathered_sizes, displacements, u)


    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    integer(kind=4), intent(in) :: local_dof        !< number of local dof
    ! integer(kind=4), intent(in) :: mpi_id
    integer(kind=4), dimension(:), allocatable, intent(out) :: nnod_num
    integer(kind=4), dimension(:), allocatable, intent(out) :: gathered_sizes !< MPI gather size
    integer(kind=4), dimension(:), allocatable, intent(out) :: displacements !< MPI gather displacements
    real(kind=8), dimension(:,:), allocatable, intent(out) :: u !< solution (Np, 3*num_poly)
    integer(kind=4) :: i, Np

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    ! Allocate solution for post-processing
    allocate(u(Np, DIM * PolyMesh%num_poly))

    ! STORE LOCAL NUMERATION TO RECONSTRUCT THE SOLUTION
    allocate(nnod_num(local_dof))
    call CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)

    if(mpi_np > 1) then

        allocate(gathered_sizes(mpi_np))

        call MPI_AllGather(local_dof, 1, MPI_INTEGER, gathered_sizes, 1, &
                    MPI_INTEGER, MPI_COMM_WORLD, ierr)

        allocate(displacements(mpi_np))
        displacements(1) = 0
        do i = 2, mpi_np
                displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
        end do
        
        ! print *, 'displacements', displacements
        ! print *, 'gathered_sizes', gathered_sizes
    endif

end subroutine PREPROCESS_SOLUTION

!> gather global solution from PETSc vector into a Fortran vector
subroutine POST_PROCESS(PolyMesh, local_dof, global_dof, petsc_sol, sol_ptr, u, nnod_num, gathered_sizes, displacements)

    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh
    real(kind=8), dimension(:,:), intent(inout) :: u        !< global solution
    integer(kind=4), intent(in) :: local_dof                !< number of local dofs for the process
    integer(kind=4), intent(in) :: global_dof               !< number of global dofs
    real(kind=8), pointer, intent(inout) :: sol_ptr(:)      !< fortran pointer

    !! NOT USED??
    integer(kind=4), dimension(:), intent(in) :: nnod_num
    real(kind=8), dimension(:), allocatable :: u_loc, u_glo
    integer(kind=4), dimension(:), intent(in) :: gathered_sizes, displacements
    integer(kind=4) :: Np

    type(tVec), intent(in) :: petsc_sol    !< local solutions
    ! assuming same degree everywhere
    Np = PolyMesh%Elem_loc(1)%NDof_elem ! ndof local per dimension (for 3D we need 3*Np)

    ! SCATTER PETSC SOLUTION AND STORE IN A FORTRAN ARRAY
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    allocate(u_loc(local_dof))
    allocate(u_glo(global_dof))

    ! print *, 'SCATTER SOLUTION'
    PetscCallA(VecGetArrayReadF90(petsc_sol, sol_ptr, mpi_ierr))
    u_loc(1:local_dof) = sol_ptr
    PetscCall(VecRestoreArrayReadF90(petsc_sol,sol_ptr,mpi_ierr))

    if(mpi_np == 1) then
        u_glo = u_loc
    else
        call MPI_ALLGATHERV(u_loc, local_dof, MPI_DOUBLE_PRECISION, &
                    u_glo, gathered_sizes, displacements, MPI_DOUBLE_PRECISION, &
                    MPI_COMM_WORLD, mpi_ierr)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    deallocate(u_loc)


    ! RECONSTRUCT SOLUTION MATRIX FOR POST-PROCESSING
    ! allocate(u(Np, 3*PolyMesh%num_poly))
    u = RESHAPE(u_glo, (/Np, DIM * PolyMesh%num_poly /))

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    deallocate(u_glo)

    ! print *,'Done with the solution'

end subroutine POST_PROCESS


!> matrix free - gather fortran vector solution from local to global.
!> Dof are split not in ascending order, but gathered as blocks in order of process
! subroutine GATHER_SOLUTION(PolyMesh, local_dof, global_dof, u_loc, mpi_id, u_glo)

!     type(Mesh_Structure), intent(in) :: PolyMesh    !< Mesh structure
!     integer(kind=4), intent(in) :: local_dof        !< Number of local dof for the process
!     integer(kind=4), intent(in) :: global_dof       !< Number of global dof
!     integer(kind=4), intent(in) :: mpi_id           !< Process ID
!     real(kind=8), dimension(:,:), intent(out) :: u_glo  !< Vector of global solution
!     real(kind=8), dimension(:,:), intent(in) :: u_loc   !< Vector of local solution for each process

!     integer(kind=4), dimension(:), allocatable :: nnod_num
!     integer(kind=4), dimension(:), allocatable :: gathered_sizes, displacements
!     integer(kind=4) :: i, Np

!     Np = PolyMesh%Elem_loc(1)%NDof_elem

!     ! STORE LOCAL NUMERATION TO RECONSTRUCT THE SOLUTION
!     allocate(nnod_num(local_dof))
!     call CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)

!     ! SCATTER PETSC SOLUTION AND STORE IN A FORTRAN ARRAY
!     call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

!     allocate(u_glo(PolyMesh%num_elem,DIM*Np))

!     !! send and receive matrices NOT VECTORS?

!     ! contruct gathered_sized for MPI_ALLGATHERV
!     if(mpi_np > 1) then

!         allocate(gathered_sizes(mpi_np))

!         call MPI_AllGather(local_dof, 1, MPI_INTEGER, gathered_sizes, 1, &
!                     MPI_INTEGER, MPI_COMM_WORLD, ierr)

!         allocate(displacements(mpi_np))
!         displacements(1) = 0
!         do i = 2, mpi_np
!                 displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
!         end do

!     endif

!     ! print *, 'SCATTER SOLUTION'

!     if(mpi_np == 1) then
!         u_glo = u_loc
!     else
!         call MPI_ALLGATHERV(u_loc, local_dof, MPI_DOUBLE_PRECISION, &
!                     u_glo, gathered_sizes, displacements, MPI_DOUBLE_PRECISION, &
!                     MPI_COMM_WORLD, mpi_ierr)
!     endif

!     ! RECONSTRUCT SOLUTION MATRIX FOR POST-PROCESSING
!     allocate(u(Np, DIM*PolyMesh%num_poly))
!     u = RESHAPE(u_glo, (/Np, DIM*PolyMesh%num_poly /))
!     deallocate(u_glo)

!     ! print *,'Done with the solution'

! end subroutine GATHER_SOLUTION

!> Evaluate nodal values of the solution and store them in various formats
!> Input global modal solution, compute local vertex solution
subroutine EXPORT_SOLUTION(PolyMesh, u, IsPoly, num_dt)
        
        !! COMPUTE POINTS ONLY ONCE?
        !! HARD CODED NUMBERS
        !! WHY COMPUTE LOCAL AGAIN AFTER GATHER??? POLYGON SPLIT IN DIFFERENT PROCESSORS?

        use local_search ! see Poly_global.f90
        use problem_data_and_properties
        use mpi
        use Poly_setup_MPI

        implicit none

        integer(kind=4), intent(in), optional :: num_dt     !< number of timesteps

        type(Mesh_Structure), intent(inout) :: PolyMesh     !< mesh
        real(kind=8), dimension(:,:), intent(in) :: u       !< global solution
        logical, intent(in) :: IsPoly                       !< boolean for polyhedra
        integer(kind=4) :: p, Np, Npoly
        integer(kind=4), dimension(:,:), allocatable :: blist
        real(kind=8), dimension(:), allocatable :: x_p, y_p, z_p
        real(kind=8), dimension(:), allocatable :: valx, valy, valz
        real(kind=8), dimension(:), allocatable :: dvalx, dvaly, dvalz
        real(kind=8), dimension(2) :: intx, inty, intz
        real(kind=8), dimension(:,:), allocatable :: temp
        real(kind=8), dimension(:,:,:), allocatable :: u_nod_vet    !< local vertex solution
        !real(kind=8), dimension(3) :: points
        integer(kind=4) :: ie_loc, ie_glob, ipoly_loc, ipoly_glob, ivert, id_node, i, j, k, index

        Np = PolyMesh%Elem_loc(1)%NDof_elem
        Npoly = PolyMesh%num_poly
        p = PolyMesh%Elem_loc(1)%Degree

        if (mpi_id==0) print *, 'Saving the solution in a file ...'
        
        ! list of the degrees of monomials of the Np basis functions up to order p
        ! (see basis_functions.f90)
        allocate(blist(Np,3))
        call basis_list(blist,p,Np)
        
        allocate(u_nod_vet(DIM,4,PolyMesh%num_elem_loc))
        allocate(temp(4,PolyMesh%num_elem_loc))

        u_nod_vet = 0.0

        do k=1,DIM

            ! Computing nodal values of the solution (on each element)
            do ie_loc = 1,PolyMesh%num_elem_loc

                allocate(x_p(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(y_p(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(z_p(PolyMesh%Elem_loc(ie_loc)%num_vert))

                do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert
                
                    ! see MAKE_PARTITION_AND_MPI_FILES.f90
                    call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                        PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)      
                    
                    x_p(ivert)=PolyMesh%coord_x(id_node)
                    y_p(ivert)=PolyMesh%coord_y(id_node)
                    z_p(ivert)=PolyMesh%coord_z(id_node)

                enddo

                ! find the polyhedron ipoly_glob that contains the tetrahedron ie_loc
                ie_glob = PolyMesh%elem_loc2glo(ie_loc)
                ipoly_glob = PolyMesh%elem_in_poly(ie_glob)

                ! see module local_search in Poly_global.f90
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                                            PolyMesh%num_poly_loc, &
                                            ipoly_glob,ipoly_loc)

                do j=1,2
                    intx(j)=PolyMesh%Poly(ipoly_loc)%b_box(1,j)
                    inty(j)=PolyMesh%Poly(ipoly_loc)%b_box(2,j)
                    intz(j)=PolyMesh%Poly(ipoly_loc)%b_box(3,j)
                end do

                allocate(valx(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(valy(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(valz(PolyMesh%Elem_loc(ie_loc)%num_vert))

                allocate(dvalx(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(dvaly(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(dvalz(PolyMesh%Elem_loc(ie_loc)%num_vert))

                index = (k-1)*Npoly + ipoly_glob

                ! evaluation of L_blist(i,j), j=1,2,3, at the vertices of the given tetrahedra
                ! and evaluation of nodal values of the solution
                do ivert=1,PolyMesh%Elem_loc(ie_loc)%num_vert
                    temp(ivert,ie_loc)=0.0
                    do i = 1,Np
                        
                        call LegendreP(valx, x_p, blist(i,1), intx, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvalx, x_p, blist(i,1), intx, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call LegendreP(valy, y_p, blist(i,2), inty, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvaly, y_p, blist(i,2), inty, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call LegendreP(valz, z_p, blist(i,3), intz, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvalz, z_p, blist(i,3), intz, PolyMesh%Elem_loc(ie_loc)%num_vert)

                        temp(ivert,ie_loc) = temp(ivert,ie_loc) + u(i,index)*valx(ivert)*valy(ivert)*valz(ivert) 

                    end do

                end do

                deallocate(x_p,y_p,z_p)
                deallocate(valx,valy,valz)
                deallocate(dvalx,dvaly,dvalz)

            enddo

            u_nod_vet(k,:,:) = temp

        enddo
        
        deallocate(temp)
        deallocate(blist)

        if (present(num_dt)) then
            call WRITE_SOLUTION(PolyMesh%num_elem_loc, PolyMesh, u_nod_vet, IsPoly, num_dt) ! see export_file_formats.f90
        else
            call WRITE_SOLUTION(PolyMesh%num_elem_loc, PolyMesh, u_nod_vet, IsPoly)
        endif

        deallocate(u_nod_vet)
        if (mpi_id==0) print *,'Done exporting solution'

    end subroutine EXPORT_SOLUTION

!> Compute L2 norm square of the error (u - u_ex) per processor
subroutine COMPUTE_ERROR_L2(mass, petsc_sol, petsc_uex, global_dof, err_L2_mpi, local_dof)

    use Poly_setup_mpi, only : mpi_ierr

    implicit none

    Mat :: mass
    Vec :: petsc_sol, petsc_uex, temp, error, e_L2
    PetscScalar coeff

    !PetscViewer viewer

    integer(kind=4) :: global_dof, local_dof
    real(kind=8), pointer, dimension(:) :: error_L2_pointer
    real(kind=8), dimension(local_dof) :: error_L2
    real(kind=8) :: err_L2_mpi

    coeff = -1

    call SET_PETSC_VECTOR(e_L2, local_dof, global_dof)
    call SET_PETSC_VECTOR(temp, local_dof, global_dof)
    call SET_PETSC_VECTOR(error, local_dof, global_dof)

    PetscCall(VecWAXPY(error, coeff, petsc_sol, petsc_uex, mpi_ierr))

    ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'error',viewer,mpi_ierr))
    ! PetscCallA(VecView(error,viewer,mpi_ierr))
    ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

    PetscCall(MatMult(mass, error, temp, mpi_ierr))

    ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'temp',viewer,mpi_ierr))
    ! PetscCallA(VecView(temp,viewer,mpi_ierr))
    ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

    PetscCall(VecPointwiseMult(e_L2, error, temp, mpi_ierr))
    PetscCallA(VecGetArrayReadF90(e_L2, error_L2_pointer, mpi_ierr))
    error_L2(1:local_dof) = error_L2_pointer
    PetscCall(VecRestoreArrayReadF90(e_L2, error_L2_pointer,mpi_ierr))
    err_L2_mpi = sum(error_L2)

    ! print *, "errL2 ^2: ", err_L2_mpi

end subroutine COMPUTE_ERROR_L2

!> Compute DG norm square of the error (u - u_ex) per processor
subroutine COMPUTE_ERROR_DG(mat_dg, petsc_sol, petsc_uex, global_dof, err_DG_mpi, local_dof)

    use Poly_setup_mpi, only : mpi_ierr

    implicit none

    type(tMat), intent(in) :: mat_dg       !< PETSc DG matrix
    type(tVec), intent(in) :: petsc_sol    !< PETSc computed solution vector
    type(tVec), intent(in) :: petsc_uex    !< PETSc exact solution vector
    Vec :: temp, error, e_DG
    PetscScalar :: coeff

    !PetscViewer viewer

    integer(kind=4), intent(in) :: global_dof   !< global number of dof
    integer(kind=4), intent(in) :: local_dof    !< local number of dof
    real(kind=8), pointer, dimension(:) :: error_DG_pointer
    real(kind=8), dimension(local_dof) :: error_DG
    real(kind=8) :: err_DG_mpi

    coeff = -1

    call SET_PETSC_VECTOR(e_DG, local_dof, global_dof)
    call SET_PETSC_VECTOR(temp, local_dof, global_dof)
    call SET_PETSC_VECTOR(error, local_dof, global_dof)

    PetscCall(VecWAXPY(error, coeff, petsc_sol, petsc_uex, mpi_ierr))

    ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Error_poly_10',viewer,mpi_ierr))
    ! PetscCallA(VecView(error,viewer,mpi_ierr))
    ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

    PetscCall(MatMult(mat_dg, error, temp, mpi_ierr))

    ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Temp_poly_10',viewer,mpi_ierr))
    ! PetscCallA(VecView(temp,viewer,mpi_ierr))
    ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

    PetscCall(VecPointwiseMult(e_DG, error, temp, mpi_ierr))
    PetscCallA(VecGetArrayReadF90(e_DG, error_DG_pointer, mpi_ierr))
    error_DG(1:local_dof) = error_DG_pointer
    PetscCall(VecRestoreArrayReadF90(e_DG, error_DG_pointer,mpi_ierr))
    err_DG_mpi = sum(error_DG)

    ! print *, "errDG ^2: ", err_DG_mpi

end subroutine COMPUTE_ERROR_DG

!> Compute the maximum value of the diameter of the elements of the mesh per processor
function compute_hmax(PolyMesh) result(hmax)

    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    real(kind=8) :: hmax                            !< maximun size
    integer(kind=4) :: ipoly_loc
    real(kind=8) :: h

    hmax = PolyMesh%Poly(1)%hk

    do ipoly_loc=2,PolyMesh%num_poly_loc

        h = PolyMesh%Poly(ipoly_loc)%hk

        if(h > hmax) hmax = h

    enddo

end function compute_hmax

end module solution_processing