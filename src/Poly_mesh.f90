!    Author: Nicoletta De Giosa
!    This file is part of the library LYMPH3D

!> @brief Module containing the definition of the structs Element, Polyhedron and Mesh_Structure. It also contains utilites
module Poly_mesh

    use Poly_global
    use Poly_exit_codes, only: EXIT_NO_ELEMENTS

    implicit none

    !> Properties of each element
    type Element

        character(len=3) :: el_type     !< Type of the element
        integer(kind=4)  :: mat_prop    !< ID for heterogeneous materials
        integer(kind=4)  :: num_vert    !< Number of vertices of the element
        integer(kind=4)  :: num_faces   !< Number of faces of the element
        integer(kind=4)  :: Degree      !< Local degree of the basis function for the element
        integer(kind=4)  :: NDof_loc    !< Local number of degrees of freedom
        integer(kind=4), dimension(:),   pointer :: vert    !< Indexes of the vertices of the element
        integer(kind=4), dimension(:,:), pointer :: faces   !< Indexes of the vertices for every face of the element

        !> Properties of the neighbor elements (6 columns):
        !> 0 - mpi_proc: processor containing the neighbor element data
        !> 1 - mat_id: ID of the neighbor material or tag of boundary face
        !> 2 - el_id: global element across face ID with -1 for Dirichlet boundary, -2 for Neumann boundary
        !> 3 - face_id: global face ID? !! DO NOT KNOW
        !> 4 - poly_id: numbering for processors dof ID? !! DO NOT KNOW - FUTURE CHECK WITH POLYGONS
        !> 5 - tag: tag from .mate file, numbering of BC with 0 if internal face
        integer(kind=4), dimension(:,:), pointer :: neigh_el

        integer(kind=4), dimension(:), pointer :: Dof_glo       !< Mapping to global degrees of freedom
        real(kind=8), dimension(:,:), pointer :: normal     !< Coordinates of the normal to each face
        real(kind=8), dimension(:), pointer :: area         !< Area of each face
        integer(kind=4), dimension(:), pointer :: flag      !!<

    end type Element

    !> Properties of each polyhedron
    type Polyhedron

        integer(kind=4),dimension(:),allocatable :: tet_in_poly     !< Number of tetrahedra contained in the polyhedron
        integer(kind=4) :: num_tet_in_poly      !< Global indexes of the tetrahedra contained in the polyhedron
        real(kind=8), dimension(3,2) :: b_box   !< Coordinates of two diametrically opposite points of the bounding box of the polyhedron
        real(kind=8) :: hk  !< Diameter of the polyhedron
        real(kind=8), dimension(:,:,:), allocatable :: neigh_bbox !< Coordinates of two diametrically opposite points of the bounding box of the neighbouring polyhedra
        real(kind=8), dimension(:),allocatable :: neigh_hk !< Diameters of the neighbouring polyhedra

    end type Polyhedron

    !> Properties of the mesh
    type Mesh_Structure
    !*******************************************************************************
    ! Mesh file parameters -
    !*******************************************************************************
        integer(kind=4) :: num_node         !< Total number of vertices
        integer(kind=4) :: num_poly         !< Total number of polyhedra
        integer(kind=4) :: num_hex          !< Total number of hexahedra
        integer(kind=4) :: num_tet          !< Total number of tetrahedra
        integer(kind=4) :: num_prysm        !< Total number of prysms
        integer(kind=4) :: num_quad         !< Total number of quadrilateral faces
        integer(kind=4) :: num_tria         !< Total number of triangular faces
        integer(kind=4) :: id_node          !< 
        integer(kind=4) :: num_elem         !< 

        integer(kind=4), dimension(:,:), allocatable :: con_hex         !< Connettivity matrix of the hexahedra
        integer(kind=4), dimension(:,:), allocatable :: con_quad        !< Connettivity matrix of the quadrilaterals on the boundary
        integer(kind=4), dimension(:,:), allocatable :: con_tet         !< Connettivity matrix of the tetrahedra
        integer(kind=4), dimension(:,:), allocatable :: con_tria        !< Connettivity matrix of the triangles on the boundary
        integer(kind=4), dimension(:,:), allocatable :: con_prysm       !< Connettivity matrix of the prysms

        integer(kind=4), dimension(:), allocatable :: part_elem
        integer(kind=4), dimension(:), allocatable :: elem_in_poly
        integer(kind=4), dimension(:), allocatable :: elem_in_poly_loc

        real(kind=8), dimension(:), allocatable :: coord_x  !< Coordinates of the vertices in x
        real(kind=8), dimension(:), allocatable :: coord_y  !< Coordinates of the vertices in y
        real(kind=8), dimension(:), allocatable :: coord_z  !< Coordinates of the vertices in z

        integer(kind=4) :: num_elem_loc     !< Local number of elements
        integer(kind=4) :: num_node_loc     !< Local number of vertices
        integer (kind=4) :: num_poly_loc        !< Local number of polyhedra
        integer(kind=4), dimension(:), allocatable :: elem_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the elements
        integer(kind=4), dimension(:), allocatable :: node_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the nodes
        integer(kind=4), dimension(:), allocatable :: poly_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the polyhedra

        type(Element), dimension(:), pointer :: Elem_loc !< For each element it stores all the properties

        type(Polyhedron),dimension(:),pointer :: Poly  !< For each polyhedron it stores all the properties

    end type Mesh_Structure

    contains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Print properties of the mesh
    subroutine print_Dime_Mesh_Structure(Struct)

        implicit none

        type(Mesh_Structure), intent(inout) :: Struct

        write(*,'(A,I8)')'#Nodes          : ', Struct%num_node
        write(*,'(A,I8)')'#Hexahedra      : ', Struct%num_hex
        write(*,'(A,I8)')'#Tetrahedra     : ', Struct%num_tet
        write(*,'(A,I8)')'#Prysms         : ', Struct%num_prysm
        write(*,'(A,I8)')'#Quad faces     : ', Struct%num_quad
        write(*,'(A,I8)')'#Tria faces     : ', Struct%num_tria


    end subroutine print_Dime_Mesh_Structure


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Allocate memory for mesh properties
    subroutine allocate_Mesh_Structure(Struct)

        implicit none

        type(Mesh_Structure), intent(inout) :: Struct

        Struct%num_elem = Struct%num_hex + Struct%num_tet + Struct%num_prysm

        if (Struct%num_elem  == 0) then
            write(*,'(A)')'ERROR! NUM OF ELEMENTS = 0'
            call EXIT(EXIT_NO_ELEMENTS)
        endif
        if (Struct%num_hex > 0) then
            allocate (Struct%con_hex(Struct%num_hex,9))
            Struct%con_hex = 0;
        endif
        if (Struct%num_tet > 0) then
            allocate (Struct%con_tet(Struct%num_tet,5))
            Struct%con_tet = 0;
        endif
        if (Struct%num_prysm > 0) then
            allocate (Struct%con_prysm(Struct%num_prysm,6))
            Struct%con_prysm = 0;
        endif

        !write(*,*) Struct%num_hex, Struct%num_tet, Struct%num_prysm, Struct%con_tet
        !read(*,*)

        if (Struct%num_quad > 0) allocate (Struct%con_quad(Struct%num_quad,5))
        if (Struct%num_tria > 0) allocate (Struct%con_tria(Struct%num_tria,5)) !it is modified so that it contains the tag of the boundary

        !Remove this for big simulation
        !allocate(Struct%vert_x(Struct%num_node), &
        !         Struct%vert_y(Struct%num_node), &
        !         Struct%vert_z(Struct%num_node))
        !allocate(Struct%part_elem(Struct%num_elem))
        allocate(Struct%elem_in_poly(Struct%num_elem))

    end subroutine allocate_Mesh_Structure


        !subroutine allocate_Poly_in_Mesh_Structure(Struct)
        !  implicit none

          !type(Mesh_Structure), intent(inout) :: Struct
          !integer(kind=4) :: ipoly

          !allocate(Struct%Poly(Struct%num_poly))

          !do ipoly=1,Struct%num_poly
          !  Struct%Poly(ipoly)%hk=0.0
          !  Struct%Poly(ipoly)%b_box(1,:)=[0.0,0.0]
          !  Struct%Poly(ipoly)%b_box(2,:)=[0.0,0.0]
          !  Struct%Poly(ipoly)%b_box(3,:)=[0.0,0.0]
          !end do

        !end subroutine allocate_Poly_in_Mesh_Structure


        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! IT WORKS ONLY IN SERIAL

    !> Save to a different file the mesh for each processor
    subroutine print_Local_Mesh_Structure_VTK(Struct)

        implicit none

        type(Mesh_Structure), intent(inout) :: Struct
        integer(kind=4) :: i

        open(50,file='mesh_visualization/mesh_partition.vtk')

        write(50,'(A)') '# vtk DataFile Version 2.0'
        write(50,'(A)') 'Comment'
        write(50,'(A)') 'ASCII'
        write(50,'(A)') 'DATASET UNSTRUCTURED_GRID'
        write(50,*)     'POINTS  ', Struct%num_node, '  float'

        do i = 1, Struct%num_node
            write(50,*) Struct%coord_x(i), Struct%coord_y(i), Struct%coord_z(i)
        enddo

        write(50,*) 'CELLS ', Struct%num_elem, 9*Struct%num_hex &
                            + 5*Struct%num_tet + 6*Struct%num_prysm

        do i = 1, Struct%num_hex
            write(50,*) 8, Struct%con_hex(i,2)-1, Struct%con_hex(i,3)-1, &
                        Struct%con_hex(i,4)-1, Struct%con_hex(i,5)-1, &
                        Struct%con_hex(i,6)-1, Struct%con_hex(i,7)-1, &
                        Struct%con_hex(i,8)-1, Struct%con_hex(i,9)-1
        enddo
        do i = 1, Struct%num_tet
            write(50,*) 4, Struct%con_tet(i,2)-1, Struct%con_tet(i,3)-1, &
                        Struct%con_tet(i,4)-1, Struct%con_tet(i,5)-1
        enddo
        do i = 1, Struct%num_prysm
            write(50,*) 5, Struct%con_prysm(i,2)-1, Struct%con_prysm(i,3)-1, &
                        Struct%con_prysm(i,4)-1, Struct%con_prysm(i,5)-1, &
                        Struct%con_prysm(i,6)-1
        enddo

        ! Vtk Cell type file formats
        write(50,*) 'CELL_TYPES ', Struct%num_elem
        do i = 1, Struct%num_hex
            write(50,*) 12  ! hexahedra
        enddo
        do i = 1, Struct%num_tet
            write(50,*) 10  ! tetrahedra
        enddo
        do i = 1, Struct%num_prysm
            write(50,*) 14 ! pyramids? !!
        enddo

        write(50,*) 'CELL_DATA ', Struct%num_elem

        write(50,*) 'SCALARS  mpi_id int 1'
        !write(50,*) 'SCALARS  poly_id int 1'
        write(50,*) 'LOOKUP_TABLE default'

        do i = 1, Struct%num_elem
            write(50,*) Struct%part_elem(i)
            !write(50,*) Struct%elem_in_poly(i)
        enddo

        close(50);

    end subroutine print_Local_Mesh_Structure_VTK

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module Poly_mesh