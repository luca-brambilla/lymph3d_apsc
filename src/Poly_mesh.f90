!    Author: Nicoletta De Giosa
!    This file is part of the library LYMPH3D
 
     module Poly_mesh
     
     use Poly_global
     use Poly_exit_codes, only: EXIT_NO_ELEMENTS
     
     implicit none 
     
     
     type Element
     
        character(len=3) :: el_type
        integer(kind=4)  :: mat_prop
        integer(kind=4)  :: num_vert, num_faces
        integer(kind=4)  :: Degree, NDof_loc
        integer(kind=4), dimension(:),   pointer :: vert 
        integer(kind=4), dimension(:,:), pointer :: faces
        integer(kind=4), dimension(:,:), pointer :: neigh_el
        integer(kind=4), dimension(:), pointer :: Dof_glo
        
        real(kind=8), dimension(:,:), pointer :: normal
        real(kind=8), dimension(:), pointer :: area
        integer(kind=4), dimension(:), pointer :: flag

     end type Element

     type Polyhedron

        integer(kind=4),dimension(:),allocatable :: tet_in_poly
        integer(kind=4) :: num_tet_in_poly
        real(kind=8), dimension(3,2) :: b_box
        real(kind=8) :: hk
        real(kind=8), dimension(:,:,:), allocatable :: neigh_bbox
        real(kind=8), dimension(:),allocatable :: neigh_hk

     end type Polyhedron

     type Mesh_Structure 
     !*******************************************************************************
     ! Mesh file parameters - 
     !*******************************************************************************
        
         integer(kind=4) :: num_node, num_poly,num_hex, num_tet, num_prysm, &
                            num_quad, num_tria, id_node, num_elem
                      
         integer(kind=4), dimension(:,:), allocatable :: con_hex, con_quad, &
                                                         con_tet, con_tria, &
                                                         con_prysm
         
         integer(kind=4), dimension(:), allocatable :: part_elem
         integer(kind=4), dimension(:), allocatable :: elem_in_poly
         integer(kind=4), dimension(:), allocatable :: elem_in_poly_loc
         
         real(kind=8), dimension(:), allocatable :: coord_x, coord_y, coord_z
         
         integer(kind=4) :: num_elem_loc, num_node_loc
         integer (kind=4) :: num_poly_loc
         integer(kind=4), dimension(:), allocatable :: elem_loc2glo
         integer(kind=4), dimension(:), allocatable :: node_loc2glo
         integer(kind=4), dimension(:), allocatable :: poly_loc2glo

         type(Element), dimension(:), pointer :: Elem_loc

         type(Polyhedron),dimension(:),pointer :: Poly

       end type Mesh_Structure
    
     contains 
     
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
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


        write(50,*) 'CELL_TYPES ', Struct%num_elem
        do i = 1, Struct%num_hex
          write(50,*) 12
        enddo
        do i = 1, Struct%num_tet
          write(50,*) 10
        enddo
        do i = 1, Struct%num_prysm
          write(50,*) 14
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

