module find_tet

    use Poly_mesh
    use global_parameters

    implicit none

    type :: MpiPair
        real(kind=8) :: val
        integer(kind=4) :: rank
    end type MpiPair

    contains

    ! compute TET baricentric coordinates, only used for sign, not value
    ! not divided by 6
    function COMPUTE_BARICENTRIC_COORD(v1, v2, v3, v4) result(vol)

        implicit none
        real(kind=8), intent(in) :: v1(3), v2(3), v3(3), v4(3)
        real(kind=8) :: vol
        real(kind=8) :: a(3), b(3), c(3), cross(3)

        ! Compute vectors a, b, c
        a = v2 - v1
        b = v3 - v1
        c = v4 - v1

        ! Compute the cross product b x c
        cross(1) = b(2) * c(3) - b(3) * c(2)
        cross(2) = b(3) * c(1) - b(1) * c(3)
        cross(3) = b(1) * c(2) - b(2) * c(1)

        ! Compute the scalar triple product a . (b x c)
        vol = dot_product(a, cross)

        ! Take the absolute value and divide by 6
        ! vol = abs(vol) / 6.0
    end function COMPUTE_BARICENTRIC_COORD

    ! if the sign of all baricentric coordinates are negative
    ! the point belongs to the tetrahedron
    ! function FIND_ELEM_FROM_POINT(PolyMesh,point)result(id_tet)

    !     implicit none

    !     type(Mesh_Structure) :: PolyMesh
    !     real(kind=8), dimension(DIM) :: point
    !     real(kind=8) :: vol1, vol2, vol3, vol4

    !     real(kind=8) :: v(DIM,NVERT_TET)

    !     integer(kind=4) :: ie_loc, id_tet, ivert_loc

    !     id_tet = -1

    !     do ie_loc=1,PolyMesh%num_elem_loc

    !         ! get coordinates of tetrahedron vertices
    !         do ivert_loc=1,NVERT_TET
    !             !ivert_glo = PolyMesh%Elem_loc(ie_loc)%vert(ivert_loc)
    !             v(1,ivert_loc)=PolyMesh%coord_x(ie_loc)
    !             v(2,ivert_loc)=PolyMesh%coord_y(ie_loc)
    !             v(3,ivert_loc)=PolyMesh%coord_z(ie_loc)
    !         enddo

    !         vol1 = COMPUTE_BARICENTRIC_COORD(point,v(:,2),v(:,3),v(:,4))
    !         vol2 = COMPUTE_BARICENTRIC_COORD(v(:,1),point,v(:,3),v(:,4))
    !         vol3 = COMPUTE_BARICENTRIC_COORD(v(:,1),v(:,2),point,v(:,4))
    !         vol4 = COMPUTE_BARICENTRIC_COORD(v(:,1),v(:,2),v(:,3),point)

    !         if (vol1>=0.0d0 .and. vol2>=0.0d0 .and. vol3>=0.0d0 .and. vol4>=0.0d0) then
    !             id_tet = ie_loc
    !             exit
    !         endif

    !     enddo

    ! end function

    !> given a point, find tet to which it belongs to
    subroutine FIND_ELEM_FROM_POINT(PolyMesh, point, id_tet, rank_id)

        use Poly_setup_MPI

        implicit none

        type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
        real(kind=8), dimension(DIM), intent(in) :: point   !< point
        real(kind=8), dimension(DIM) :: baryc
        ! real(kind=8), dimension(NVERT_TET) :: x,y,z, xmin,ymin,zmin
        integer(kind=4), intent(out) :: id_tet  !< local tet id
        integer(kind=4), intent(out) :: rank_id !< rank of process containing tet
        integer(kind=4) :: ie_loc
        real(kind=8) :: d2, dmin

        type(MpiPair) :: local, global

        dmin = 1.0d15
        id_tet = -1

        do ie_loc=1,PolyMesh%num_elem_loc
            baryc = COMPUTE_BARYCENTER(PolyMesh,ie_loc)
            d2 = POINT_DISTANCE_SQUARED(point,baryc)
            if (d2<dmin) then
                dmin=d2
                id_tet = ie_loc
                ! xmin=x
                ! ymin=y
                ! zmin=z
            endif
        enddo

        local%val = dmin
        local%rank = mpi_id

        ! pair of double and integer: distance and process
        call MPI_ALLREDUCE(local, global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD, ierr)

        rank_id = global%rank

        !if (mpi_id==rank_id) write(*,'(A,I0,A,I0)') 'mpi_id: ', mpi_id, ' tet_id: ', id_tet

    end subroutine FIND_ELEM_FROM_POINT

    subroutine TEST_TET(PolyMesh)
        use Poly_setup_MPI

        implicit none
        type(Mesh_Structure), intent(in) :: PolyMesh
        integer(kind=4) :: rank_id, id_tet

        call FIND_ELEM_FROM_POINT(PolyMesh, (/ 0.7d0, 0.4d0, 0.6d0 /),id_tet,rank_id)

        print *, mpi_id, rank_id, id_tet
    end subroutine TEST_TET

    !> compute the square of the distance between two points
    function POINT_DISTANCE_SQUARED(p1,p2) result(d2)
        implicit none

        real(kind=8), dimension(DIM) :: p1 !< first point coordinates
        real(kind=8), dimension(DIM) :: p2 !< second point coordinates
        real(kind=8) :: d2  !< square of the distance

        d2 = (p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2

    end function POINT_DISTANCE_SQUARED

    !> compute element barycenter
    function COMPUTE_BARYCENTER(PolyMesh,ie_loc) result(point)
        use mesh_partition_and_mpi_files
        implicit none

        type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
        real(kind=8), dimension(DIM) :: point   !< barycenter coordinates
        integer(kind=4) :: ie_loc   !< local element id
        real(kind=8), dimension(PolyMesh%Elem_loc(ie_loc)%num_vert) :: x,y,z
        integer(kind=4) :: ivert, id_node


        ! computation of the coordinates of the tetrahedron
        do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert

            ! see MAKE_PARTITION_AND_MPI_FILES.f90
            call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                    PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)

            x(ivert)=PolyMesh%coord_x(id_node)
            y(ivert)=PolyMesh%coord_y(id_node)
            z(ivert)=PolyMesh%coord_z(id_node)

        enddo

        point(1) = sum(x)
        point(2) = sum(y)
        point(3) = sum(z)

        point = point/PolyMesh%Elem_loc(ie_loc)%num_vert

    endfunction COMPUTE_BARYCENTER


    function volume(PolyMesh,ie_loc) result(vol)
        use mesh_partition_and_mpi_files
        use vet_mat_operations
        implicit none

        type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
        real(kind=8), dimension(DIM+1,DIM+1) :: mat   !< barycenter coordinates
        integer(kind=4) :: ie_loc   !< local element id
        integer(kind=4) :: ivert, id_node
        real(kind=8) :: vol

        ! computation of the coordinates of the tetrahedron
        do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert

            ! see MAKE_PARTITION_AND_MPI_FILES.f90
            call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                    PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)

            mat(ivert,1) = PolyMesh%coord_x(id_node)
            mat(ivert,2) = PolyMesh%coord_y(id_node)
            mat(ivert,3) = PolyMesh%coord_z(id_node)
            mat(ivert,4) = 1.0d0
        enddo

        vol = abs(det4(mat)) / 6.0d0

    end function volume

end module find_tet