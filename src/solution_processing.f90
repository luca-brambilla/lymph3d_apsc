!> pre-process, post-process, export solution and compute errors
module solution_processing

    use Poly_setup_mpi
    use Poly_global
    use problem_data_and_properties
    use Poly_mesh
    use Poly_ref_mappings
    use basis_function
    use export_file_formats
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
                    MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)

        allocate(displacements(mpi_np))
        displacements(1) = 0
        do i = 2, mpi_np
                displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
        end do
        
        ! print *, 'displacements', displacements
        ! print *, 'gathered_sizes', gathered_sizes
    endif

end subroutine PREPROCESS_SOLUTION

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