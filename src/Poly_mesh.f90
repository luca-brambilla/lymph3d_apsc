!    Author: Nicoletta De Giosa
!    This file is part of the library LYMPH3D

!> @brief Module containing the definition of the structs Element, Polyhedron and Mesh_Structure. It also contains utilites
module Poly_mesh

    use Poly_global
    use Poly_exit_codes, only: EXIT_NO_ELEMENTS

    implicit none

    !> Properties of each element - tetrahedron (or hexahedron?)
    type Element

        character(len=3) :: el_type     !< Type of the element
        integer(kind=4)  :: mat_prop    !< ID for heterogeneous materials
        integer(kind=4)  :: num_vert    !< Number of vertices of the element
        integer(kind=4)  :: num_faces   !< Number of faces of the element
        integer(kind=4)  :: Degree      !< Degree of the basis function for the element
        integer(kind=4)  :: NDof_elem    !< Number of degrees of freedom per dimension
        integer(kind=4), dimension(:),   pointer :: vert    !< Indexes of the vertices of the element
        integer(kind=4), dimension(:,:), pointer :: faces   !< Indexes of the vertices for every face of the element

        !> Properties of the neighbor elements (6 columns)
        !> - `0` - mpi_proc: processor containing the neighbor element data
        !> - `1` - mat_id: ID of the neighbor material or tag of boundary face
        !> - `2` - el_id: global element across face ID with `-1` for Dirichlet boundary, `-2` for Neumann boundary @n
        !> - `3` - face_id: ID of neighbor local face ID (1-4 tetrahedron, 1-6 hexahedron) @n
        !> - `4` - poly_id: global ID of the polyhedron the element belongs to @n
        !> - `5` - tag: tag from `.mate` file, numbering of BC with `0` if internal face
        integer(kind=4), dimension(:,:), pointer :: neigh_el

        integer(kind=4), dimension(:), pointer :: Dof_glo       !< Mapping to global degrees of freedom
        real(kind=8), dimension(:,:), pointer :: normal     !< Coordinates of the normal to each face
        real(kind=8), dimension(:), pointer :: area         !< Area of each face
        integer(kind=4), dimension(:), pointer :: flag      !! not used?

    end type Element

    !> Properties of each polyhedron
    type Polyhedron

        integer(kind=4),dimension(:),allocatable :: tet_in_poly     !< Number of tetrahedra contained in the polyhedron
        integer(kind=4) :: num_tet_in_poly      !< Global indexes of the tetrahedra contained in the polyhedron
        real(kind=8), dimension(3,2) :: b_box   !< Coordinates of two diametrically opposite points of the bounding box of the polyhedron
        real(kind=8) :: hk  !< Diameter of the polyhedron
        real(kind=8), dimension(:,:,:), allocatable :: neigh_bbox !< Coordinates of two diametrically opposite points of the bounding box of the neighbouring polyhedra
        real(kind=8), dimension(:), allocatable :: neigh_hk !< Diameters of the neighbouring polyhedra

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
        integer(kind=4) :: id_node          !! not used?
        integer(kind=4) :: num_elem         !< Sum of number of tetrahedra, hexahedra and prysms

        integer(kind=4), dimension(:,:), allocatable :: con_hex         !< Connettivity matrix of the hexahedra
        integer(kind=4), dimension(:,:), allocatable :: con_quad        !< Connettivity matrix of the quadrilaterals on the boundary
        integer(kind=4), dimension(:,:), allocatable :: con_tet         !< Connettivity matrix of the tetrahedra
        integer(kind=4), dimension(:,:), allocatable :: con_tria        !< Connettivity matrix of the triangles on the boundary
        integer(kind=4), dimension(:,:), allocatable :: con_prysm       !< Connettivity matrix of the prysms

        integer(kind=4), dimension(:), allocatable :: part_elem         !< For each tetrahedron, it stores the ID of the process to which it belongs to after the partition
        integer(kind=4), dimension(:), allocatable :: elem_in_poly      !< For each tetrahedron, it stores the index of the polyhedron to which it belongs to
        integer(kind=4), dimension(:), allocatable :: elem_in_poly_loc  !< For each local tetrahedron, it stores the global index of the polyhedron to which it belongs to

        real(kind=8), dimension(:), allocatable :: coord_x  !< Coordinates of the vertices in x
        real(kind=8), dimension(:), allocatable :: coord_y  !< Coordinates of the vertices in y
        real(kind=8), dimension(:), allocatable :: coord_z  !< Coordinates of the vertices in z

        integer(kind=4) :: num_elem_loc     !< Local number of elements
        integer(kind=4) :: num_node_loc     !< Local number of vertices
        integer (kind=4) :: num_poly_loc        !< Local number of polyhedra
        integer(kind=4), dimension(:), allocatable :: elem_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the elements
        integer(kind=4), dimension(:), allocatable :: elem_glo2loc  !< Global-to-local maps to transition from the global enumeration to the local one for the elements
        integer(kind=4), dimension(:), allocatable :: node_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the nodes
        integer(kind=4), dimension(:), allocatable :: poly_loc2glo  !< Local-to-global maps to transition from the local enumeration to the global one for the polyhedra

        !! change names num_elem_inter_comm and num_elem_inter
        integer(kind=4), dimension(:,:), allocatable :: num_elem_inter_comm !< Matrix containing the number of interface elements per process for send and receive communication
        integer(kind=4) :: num_elem_inter !< Total number of interface elements
        integer(kind=4), dimension(:,:), allocatable :: inter_disp !< indices for elem_inter setting the initial position for data of each process
        integer(kind=4), dimension(:), allocatable :: elem_inter_glo !< interface element global ID per process. IDs stored by receive process, then for each stored by send process.
        integer(kind=4), dimension(:), allocatable :: elem_inter_loc !< interface element local ID per process. IDs stored by receive process, then for each stored by send process.

        type(Element), dimension(:), pointer :: Elem_loc !< For each element it stores all the properties

        type(Polyhedron),dimension(:),pointer :: Poly  !< For each polyhedron it stores all the properties

    end type Mesh_Structure

contains

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> Print properties of the mesh
subroutine print_Dime_Mesh_Structure(PolyMesh)

    implicit none

    type(Mesh_Structure), intent(inout) :: PolyMesh

    write(*,'(A,I8)')'#Nodes          : ', PolyMesh%num_node
    write(*,'(A,I8)')'#Hexahedra      : ', PolyMesh%num_hex
    write(*,'(A,I8)')'#Tetrahedra     : ', PolyMesh%num_tet
    write(*,'(A,I8)')'#Prysms         : ', PolyMesh%num_prysm
    write(*,'(A,I8)')'#Quad faces     : ', PolyMesh%num_quad
    write(*,'(A,I8)')'#Tria faces     : ', PolyMesh%num_tria


end subroutine print_Dime_Mesh_Structure


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> Allocate memory for mesh properties
subroutine allocate_Mesh_Structure(PolyMesh)

    implicit none

    type(Mesh_Structure), intent(inout) :: PolyMesh

    PolyMesh%num_elem = PolyMesh%num_hex + PolyMesh%num_tet + PolyMesh%num_prysm

    if (PolyMesh%num_elem  == 0) then
        write(*,'(A)')'ERROR! NUM OF ELEMENTS = 0'
        call EXIT(EXIT_NO_ELEMENTS)
    endif
    if (PolyMesh%num_hex > 0) then
        allocate (PolyMesh%con_hex(PolyMesh%num_hex,9))
        PolyMesh%con_hex = 0;
    endif
    if (PolyMesh%num_tet > 0) then
        allocate (PolyMesh%con_tet(PolyMesh%num_tet,5))
        PolyMesh%con_tet = 0;
    endif
    if (PolyMesh%num_prysm > 0) then
        allocate (PolyMesh%con_prysm(PolyMesh%num_prysm,6))
        PolyMesh%con_prysm = 0;
    endif

    !write(*,*) PolyMesh%num_hex, PolyMesh%num_tet, PolyMesh%num_prysm, PolyMesh%con_tet
    !read(*,*)

    if (PolyMesh%num_quad > 0) allocate (PolyMesh%con_quad(PolyMesh%num_quad,5))
    if (PolyMesh%num_tria > 0) allocate (PolyMesh%con_tria(PolyMesh%num_tria,5)) !it is modified so that it contains the tag of the boundary

    !Remove this for big simulation
    !allocate(PolyMesh%vert_x(PolyMesh%num_node), &
    !         PolyMesh%vert_y(PolyMesh%num_node), &
    !         PolyMesh%vert_z(PolyMesh%num_node))
    !allocate(PolyMesh%part_elem(PolyMesh%num_elem))
    allocate(PolyMesh%elem_in_poly(PolyMesh%num_elem))
    allocate(PolyMesh%elem_glo2loc(PolyMesh%num_elem))

end subroutine allocate_Mesh_Structure


        !subroutine allocate_Poly_in_Mesh_Structure(PolyMesh)
        !  implicit none

          !type(Mesh_Structure), intent(inout) :: PolyMesh
          !integer(kind=4) :: ipoly

          !allocate(PolyMesh%Poly(PolyMesh%num_poly))

          !do ipoly=1,PolyMesh%num_poly
          !  PolyMesh%Poly(ipoly)%hk=0.0
          !  PolyMesh%Poly(ipoly)%b_box(1,:)=[0.0,0.0]
          !  PolyMesh%Poly(ipoly)%b_box(2,:)=[0.0,0.0]
          !  PolyMesh%Poly(ipoly)%b_box(3,:)=[0.0,0.0]
          !end do

        !end subroutine allocate_Poly_in_Mesh_Structure


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! IT WORKS ONLY IN SERIAL

!> Save to a different file the mesh for each processor
subroutine print_Local_Mesh_Structure_VTK(PolyMesh)

    implicit none

    type(Mesh_Structure), intent(inout) :: PolyMesh
    integer(kind=4) :: i

    open(50,file='mesh_visualization/mesh_partition.vtk')

    write(50,'(A)') '# vtk DataFile Version 2.0'
    write(50,'(A)') 'Comment'
    write(50,'(A)') 'ASCII'
    write(50,'(A)') 'DATASET UNSTRUCTURED_GRID'
    write(50,*)     'POINTS  ', PolyMesh%num_node, '  float'

    do i = 1, PolyMesh%num_node
        write(50,*) PolyMesh%coord_x(i), PolyMesh%coord_y(i), PolyMesh%coord_z(i)
    enddo

    write(50,*) 'CELLS ', PolyMesh%num_elem, 9*PolyMesh%num_hex &
                        + 5*PolyMesh%num_tet + 6*PolyMesh%num_prysm

    do i = 1, PolyMesh%num_hex
        write(50,*) 8, PolyMesh%con_hex(i,2)-1, PolyMesh%con_hex(i,3)-1, &
                    PolyMesh%con_hex(i,4)-1, PolyMesh%con_hex(i,5)-1, &
                    PolyMesh%con_hex(i,6)-1, PolyMesh%con_hex(i,7)-1, &
                    PolyMesh%con_hex(i,8)-1, PolyMesh%con_hex(i,9)-1
    enddo
    do i = 1, PolyMesh%num_tet
        write(50,*) 4, PolyMesh%con_tet(i,2)-1, PolyMesh%con_tet(i,3)-1, &
                    PolyMesh%con_tet(i,4)-1, PolyMesh%con_tet(i,5)-1
    enddo
    do i = 1, PolyMesh%num_prysm
        write(50,*) 5, PolyMesh%con_prysm(i,2)-1, PolyMesh%con_prysm(i,3)-1, &
                    PolyMesh%con_prysm(i,4)-1, PolyMesh%con_prysm(i,5)-1, &
                    PolyMesh%con_prysm(i,6)-1
    enddo

    ! Vtk Cell type file formats
    write(50,*) 'CELL_TYPES ', PolyMesh%num_elem
    do i = 1, PolyMesh%num_hex
        write(50,*) 12  ! hexahedra
    enddo
    do i = 1, PolyMesh%num_tet
        write(50,*) 10  ! tetrahedra
    enddo
    do i = 1, PolyMesh%num_prysm
        write(50,*) 14 ! pyramids? !!
    enddo

    write(50,*) 'CELL_DATA ', PolyMesh%num_elem

    write(50,*) 'SCALARS  mpi_id int 1'
    !write(50,*) 'SCALARS  poly_id int 1'
    write(50,*) 'LOOKUP_TABLE default'

    do i = 1, PolyMesh%num_elem
        write(50,*) PolyMesh%part_elem(i)
        !write(50,*) PolyMesh%elem_in_poly(i)
    enddo

    close(50);

end subroutine print_Local_Mesh_Structure_VTK

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module Poly_mesh