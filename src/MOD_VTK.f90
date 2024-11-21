!> Module containing type definitions and routines for real-time VTK output during computations for individual material blocks
module MOD_VTK

    use Poly_mesh

    implicit none

    contains

    !>
    !! change types?
    !! hardcoded 3 for 3D and 4 for tetrahedra vertices
    !! hardcoded VTK type for tetrahedra 10
subroutine VTK_WRITE_SOLUTION(filename, xx, yy, zz, nnode_per_el, n_elem, n_elem_tot, u_name, u, PolyMesh)

    use mpi
    use Poly_setup_MPI

    implicit none

    ! Input arguments
    character(len=*), intent(in) :: filename                !< string file name
    integer*4, intent(in) :: n_elem                     !< local number of dofs for the process
    integer*4, intent(in) :: nnode_per_el               !< number of nodes per element
    real*8, dimension(4,n_elem), intent(in) :: xx       !< x coordinates of 4 vertices of tetrahedron
    real*8, dimension(4,n_elem), intent(in) ::yy        !< y coordinates of 4 vertices of tetrahedron
    real*8, dimension(4,n_elem), intent(in) ::zz       !< z coordinates of 4 vertices of tetrahedron
    integer*4, dimension(n_elem_tot,4) :: tnew
    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh

    ! Optional input arguments
    character(len=*), intent(in), optional :: u_name            !< solution name
    real*8, dimension(3,4,n_elem), intent(inout), optional :: u !< solution

    ! Internal variables
    integer*4 :: VTK_file_unit
    integer*4 :: i,j,k

    ! added
    integer*4, intent(in) :: n_elem_tot

    integer*4 :: tmp_nelem
    real*8, dimension(:,:,:), allocatable :: tmp_solution
    real*8, dimension(:,:), allocatable :: tmp_xx, tmp_yy, tmp_zz
    integer*4, dimension(:), allocatable :: tmp_poly
    !integer(kind=4), intent(in) :: mpi_id

    ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
    real*8, parameter :: cap = 1E-40

        ! call MPI_REDUCE(err_L2_mpi, err_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        ! 0, MPI_COMM_WORLD, mpi_ierr)

    ! do i=1,n_elem
    !     do j=1,4
    !         tnew(i,j)=(i-1)*4+j-1
    !     end do
    ! end do

    if (mpi_id /= 0) then
        ! tag=1 before MPI_COMM_WORLD
        ! send size
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send coordinates
        call MPI_Send(xx, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send solution
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u, n_elem*4*3, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send polyhedra
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(PolyMesh%elem_in_poly_loc, n_elem, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)

    elseif (mpi_id == 0) then

        open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', status='replace')

        write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
        write(VTK_file_unit,'(A)')'VTKFile'
        write(VTK_file_unit,'(A)')'ASCII'
        write(VTK_file_unit,'(A)')
        write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

        ! ******************
        ! POINTS
        ! ******************
        write(VTK_file_unit,'(A7,I12,A7)')'POINTS ', n_elem_tot*4,' double';
        POINT_LOOP: do i=1,n_elem
        !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            do j=1,4
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,i),yy(j,i),zz(j,i)
            end do
        end do POINT_LOOP

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive coordinates - different sizes
            allocate(tmp_xx(4,tmp_nelem))
            allocate(tmp_yy(4,tmp_nelem))
            allocate(tmp_zz(4,tmp_nelem))
            call MPI_Recv(tmp_xx, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_yy, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_zz, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write to file
            POINT_LOOP_k: do i=1,tmp_nelem
                do j=1,4
                    write(VTK_file_unit,'(3F16.8,1X)') tmp_xx(j,i),tmp_yy(j,i),tmp_zz(j,i)
                end do
            end do POINT_LOOP_k
            ! deallocate coordinates, different sizes
            deallocate(tmp_xx, tmp_yy, tmp_zz)
        enddo
        write(VTK_file_unit,'(A)')

        ! **************
        ! CELLS
        ! **************
        write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem_tot, n_elem_tot*(nnode_per_el+1);
        ELEM_LOOP: do i=1,n_elem_tot
            ! ordering
            !! 4 hardcoded
            do j=1,4
                tnew(i,j)=(i-1)*4+j-1
            end do
            ! write ordering
            write(VTK_file_unit,'(I12,8I12)')nnode_per_el, tnew(i,:);
            !write(VTK_file_unit,'(I12,8I12)')nnode_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
            !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1;
        end do ELEM_LOOP
        write(VTK_file_unit,'(A)')

        ! ******************
        ! CELL TYPES
        ! ******************
        !! hardcoded tetrahedra
        write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem_tot;
            ELEM_TYPE_LOOP: do i=1,n_elem_tot
                write(VTK_file_unit,'(I2)')10;
            end do ELEM_TYPE_LOOP
        write(VTK_file_unit,'(A)')

        ! *********************
        ! SOLUTION
        ! *********************
        if (present(u_name) .and. present(u)) then
            ! process 0 data
            write(VTK_file_unit,'(A11,I12)')'POINT_DATA ', n_elem_tot*4;
            write(VTK_file_unit,'(A8,A12,A7)')'VECTORS ', u_name, ' double'
            !write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'
            VECTOR_FIELD_LOOP: do i=1,n_elem
                do j=1,4
                    write(VTK_file_unit,'(3F16.8)') u(1,j,i), u(2,j,i), u(3,j,i)
                    !write(VTK_file_unit,'(3F16.8)') u(1:3,j,i)
                enddo
            end do VECTOR_FIELD_LOOP

            ! mpi loop
            do k=1,mpi_np-1
                ! receive size
                call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
                ! allocate and receive solution - different sizes
                allocate(tmp_solution(3,4,tmp_nelem))
                call MPI_Recv(tmp_solution, tmp_nelem*4*3, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
                ! write solution
                VECTOR_FIELD_LOOP_k: do i=1,tmp_nelem
                    do j=1,4
                        write(VTK_file_unit,'(3F16.8)') tmp_solution(1,j,i), tmp_solution(2,j,i), tmp_solution(3,j,i)
                    enddo
                end do VECTOR_FIELD_LOOP_k
                ! deallocate solution - different sizes
                deallocate(tmp_solution)
            enddo
            write(VTK_file_unit,'(A)')
        endif

        ! **********************
        ! POLYHEDRA
        ! **********************
        write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem_tot;
        write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','POLYHEDRA', ' int 1'
        write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'

        POLY_LOOP1: do i=1,n_elem
            write(VTK_file_unit,'(I10)') PolyMesh%elem_in_poly_loc(i)
        end do POLY_LOOP1

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive polyhedra - different sizes
            allocate(tmp_poly(tmp_nelem))
            call MPI_Recv(tmp_poly, tmp_nelem, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write polyhedra
            POLY_LOOP_k: do i=1,tmp_nelem
                write(VTK_file_unit,'(I10)') tmp_poly(i)
            end do POLY_LOOP_k
            ! deallocate polyhedra - different sizes
            deallocate(tmp_poly)
        enddo

        write(VTK_file_unit,'(A)')

        close(unit=VTK_file_unit)
    endif

end subroutine VTK_WRITE_SOLUTION

!< Initial lines for PVD file
subroutine PVD_SETUP(num_dt_mon, dt_start, num_dt)

    implicit none

    ! input variables
    integer*4, intent(in) :: num_dt_mon        !< number of steps to save solution
    integer*4, intent(in) :: dt_start          !< start counter for solution
    integer*4, intent(in) :: num_dt            !< Total number of timesteps (last timestep)

    ! internal
    integer*4 :: current_dt
    integer*4 :: PVD_file_unit

    ! Open or create the .pvd file
    open(newunit=PVD_file_unit, file='MONITORS/solution.pvd', action='WRITE', status='replace', form='FORMATTED')

    ! Write XML header
    write(PVD_file_unit, '(A)') '<?xml version="1.0"?>'
    write(PVD_file_unit, '(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
    write(PVD_file_unit, '(A)') '  <Collection>'

    ! Write references to each solution file
    current_dt = dt_start
    write(PVD_file_unit, '(A, I0, A, A)') '    <DataSet timestep="', current_dt, '" part="0" file="mesh.vtu"/>'
    do while (current_dt < num_dt)
        write(PVD_file_unit, '(A, I0, A, I0, A)') '    <DataSet timestep="', current_dt, '" part="1" file="solution_', current_dt, '.vtu"/>'
        current_dt = current_dt + num_dt_mon
    end do
    write(PVD_file_unit, '(A, I0, A, I0, A)') '    <DataSet timestep="', num_dt, '" part="1" file="solution_', num_dt, '.vtu"/>'

    ! Close XML
    write(PVD_file_unit, '(A)') '  </Collection>'
    write(PVD_file_unit, '(A)') '</VTKFile>'

    close(PVD_file_unit)

end subroutine PVD_SETUP



!> Write mesh and solution into 2 separate VTU files
subroutine VTU_WRITE_SOLUTION(filename, xx, yy, zz, nnode_per_el, n_elem, n_elem_tot, u_name, u, PolyMesh, num_dt)

    use mpi
    use Poly_setup_MPI

    implicit none

    ! Input arguments
    character(len=*), intent(in) :: filename                !< string file name
    integer*4, intent(in) :: n_elem                     !< local number of dofs for the process
    integer*4, intent(in) :: nnode_per_el               !< number of nodes per element
    real*8, dimension(4,n_elem), intent(in) :: xx       !< x coordinates of 4 vertices of tetrahedron
    real*8, dimension(4,n_elem), intent(in) ::yy        !< y coordinates of 4 vertices of tetrahedron
    real*8, dimension(4,n_elem), intent(in) ::zz       !< z coordinates of 4 vertices of tetrahedron
    integer*4, dimension(n_elem_tot,4) :: tnew
    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh
    integer*4, intent(in) :: num_dt                            !< number of timestep

    ! Optional input arguments
    character(len=*), intent(in), optional :: u_name            !< solution name
    real*8, dimension(3,4,n_elem), intent(inout), optional :: u !< solution

    ! Internal variables
    integer*4 :: VTU_file_unit
    integer*4 :: i,j,k, node_id

    ! added
    integer*4, intent(in) :: n_elem_tot

    integer*4 :: tmp_nelem
    real*8, dimension(:,:,:), allocatable :: tmp_solution
    real*8, dimension(:,:), allocatable :: tmp_xx, tmp_yy, tmp_zz
    integer*4, dimension(:), allocatable :: tmp_poly
    !integer(kind=4), intent(in) :: mpi_id

    character(len=256) :: mesh_filename, solution_filename
    integer*4 :: global_offset, elem_offset, conn_offset

    ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
    real*8, parameter :: cap = 1E-40

        ! call MPI_REDUCE(err_L2_mpi, err_L2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        ! 0, MPI_COMM_WORLD, mpi_ierr)

    ! do i=1,n_elem
    !     do j=1,4
    !         tnew(i,j)=(i-1)*4+j-1
    !     end do
    ! end do

    if (mpi_id /= 0) then
        ! tag=1 before MPI_COMM_WORLD
        ! send size
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send coordinates
        call MPI_Send(xx, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, n_elem*4, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send solution
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u, n_elem*4*3, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send polyhedra
        ! call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! call MPI_Send(PolyMesh%elem_in_poly_loc, n_elem, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)

    elseif (mpi_id == 0) then

        ! Define filenames
        write(mesh_filename, '("MONITORS/mesh.vtu")')
        write(solution_filename, '("MONITORS/solution_", I0, ".vtu")') num_dt

        ! ***************************************
        ! MESH
        ! ***************************************
        open(newunit=VTU_file_unit, action='WRITE', file=mesh_filename, &
            form='FORMATTED', status='replace')

        ! XML Header
        write(VTU_file_unit, '(A)') '<?xml version="1.0"?>'
        write(VTU_file_unit, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(VTU_file_unit, '(A)') '  <UnstructuredGrid>'
        write(VTU_file_unit, '(A27,I0,A17,I0,A2)') '    <Piece NumberOfPoints="', n_elem_tot*4, '" NumberOfCells="', n_elem_tot, '">'

        ! ******
        ! POINTS Section
        write(VTU_file_unit, '(A)') '    <Points>'
        write(VTU_file_unit, '(A)') '      <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
        do i = 1, n_elem
            do j = 1, 4
                write(VTU_file_unit, '(3(F16.8,1X))') xx(j, i), yy(j, i), zz(j, i)
            end do
        end do

        ! mpi loop
        do k=1,mpi_np-1

            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

            ! allocate and receive coordinates - different sizes
            allocate(tmp_xx(4,tmp_nelem))
            allocate(tmp_yy(4,tmp_nelem))
            allocate(tmp_zz(4,tmp_nelem))
            call MPI_Recv(tmp_xx, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_yy, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_zz, tmp_nelem*4, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

            ! write to file from each process
            POINT_LOOP_k: do i=1,tmp_nelem
                do j=1,4
                    write(VTU_file_unit,'(3(F16.8,1X))') tmp_xx(j,i),tmp_yy(j,i),tmp_zz(j,i)
                end do
            end do POINT_LOOP_k

            ! deallocate coordinates, different sizes
            deallocate(tmp_xx, tmp_yy, tmp_zz)
        enddo

        write(VTU_file_unit, '(A)') '      </DataArray>'
        write(VTU_file_unit, '(A)') '    </Points>'

        ! **************
        ! CELLS
        write(VTU_file_unit, '(A)') '    <Cells>'

        ! connectivity - Node indices for each cell
        write(VTU_file_unit, '(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
        elem_offset = 0
        do i = 1, n_elem_tot
            do j = 1, 4
                node_id = (i-1)*4 + (j-1)
                write(VTU_file_unit, '(I12)', advance="no") node_id
                if (j < 4) write(VTU_file_unit, '(A)', advance="no") " "
            end do
            write(VTU_file_unit, '(A)')
        end do
        write(VTU_file_unit, '(A)') '      </DataArray>'

        ! offset - Cumulative indices in the connectivity array
        write(VTU_file_unit, '(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        conn_offset = 0
        do i = 1, n_elem_tot
            conn_offset = conn_offset + 4
            write(VTU_file_unit, '(I12)') conn_offset
        end do
        write(VTU_file_unit, '(A)') '      </DataArray>'

        ! element type
        write(VTU_file_unit, '(A)') '      <DataArray type="UInt8" Name="types" format="ascii">'
        do i = 1, n_elem_tot
            write(VTU_file_unit, '(I12)') 10  ! VTK_TETRA (code 10)
        end do
        write(VTU_file_unit, '(A)') '      </DataArray>'

        write(VTU_file_unit, '(A)') '    </Cells>'

        write(VTU_file_unit, '(A)') '    </Piece>'
        write(VTU_file_unit, '(A)') '  </UnstructuredGrid>'
        write(VTU_file_unit, '(A)') '</VTKFile>'
        close(VTU_file_unit)

        ! *********************
        ! SOLUTION
        ! *********************
        if (present(u)) then
            open(newunit=VTU_file_unit, action='WRITE', file=solution_filename, &
                form='FORMATTED', status='replace')

            ! XML Header
            write(VTU_file_unit, '(A)') '<?xml version="1.0"?>'
            write(VTU_file_unit, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
            write(VTU_file_unit, '(A)') '  <UnstructuredGrid>'
            write(VTU_file_unit, '(A27,I0,A2)') '    <Piece NumberOfPoints="', n_elem_tot*4, '">'

            ! POINT_DATA Section
            write(VTU_file_unit, '(A)') '    <PointData>'
            write(VTU_file_unit, '(A)') '      <DataArray type="Float64" Name="' // trim(u_name) // '" NumberOfComponents="3" format="ascii">'
            do i = 1, n_elem
                do j = 1, 4
                    write(VTU_file_unit, '(3(F16.8,1X))') u(1, j, i), u(2, j, i), u(3, j, i)
                end do
            end do

            ! mpi loop
            do k=1,mpi_np-1

                ! receive size
                call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

                ! allocate and receive solution - different sizes
                allocate(tmp_solution(3,4,tmp_nelem))
                call MPI_Recv(tmp_solution, tmp_nelem*4*3, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

                ! write solution
                VECTOR_FIELD_LOOP_k: do i=1,tmp_nelem
                    do j=1,4
                        write(VTU_file_unit,'(3(F16.8,1X))') tmp_solution(1,j,i), tmp_solution(2,j,i), tmp_solution(3,j,i)
                    enddo
                end do VECTOR_FIELD_LOOP_k
                ! deallocate solution - different sizes
                deallocate(tmp_solution)
            enddo

            ! POINT_DATA section end
            write(VTU_file_unit, '(A)') '      </DataArray>'
            write(VTU_file_unit, '(A)') '    </PointData>'

            ! end file
            write(VTU_file_unit, '(A)') '    </Piece>'
            write(VTU_file_unit, '(A)') '  </UnstructuredGrid>'
            write(VTU_file_unit, '(A)') '</VTKFile>'
            close(VTU_file_unit)
        end if

        ! **********************
        ! POLYHEDRA
        ! **********************
        ! write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem_tot;
        ! write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','POLYHEDRA', ' int 1'
        ! write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'

        ! POLY_LOOP1: do i=1,n_elem
        !     write(VTK_file_unit,'(I6)') PolyMesh%elem_in_poly_loc(i)
        ! end do POLY_LOOP1

        ! ! mpi loop
        ! do k=1,mpi_np-1
        !     ! receive size
        !     call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
        !     ! allocate and receive polyhedra - different sizes
        !     allocate(tmp_poly(tmp_nelem))
        !     call MPI_Recv(tmp_poly, tmp_nelem, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
        !     ! write polyhedra
        !     POLY_LOOP_k: do i=1,tmp_nelem
        !         write(VTK_file_unit,'(I6)') tmp_poly(i)
        !     end do POLY_LOOP_k
        !     ! deallocate polyhedra - different sizes
        !     deallocate(tmp_poly)
        ! enddo

        ! write(VTK_file_unit,'(A)')

    endif

end subroutine VTU_WRITE_SOLUTION

    subroutine VTK_WRITE_MESH_PARTITION(filename, xx, yy, zz, nnode_per_el, n_elem, PolyMesh)

      implicit none

      ! Input arguments
      character(len=*), intent(in) :: filename
      integer*4, intent(in) :: n_elem,nnode_per_el
      real*8, dimension(4,n_elem), intent(in) :: xx,yy,zz
      integer*4,dimension(n_elem,4) :: tnew
      type(Mesh_Structure), intent(in) :: PolyMesh

      ! Internal variables
      integer*4 :: VTK_file_unit
      integer*4 :: i,j

      ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
      real*8, parameter :: cap = 1E-40

      do i=1,n_elem
        do j=1,4
          tnew(i,j)=(i-1)*4+j-1
        end do
      end do

      open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', status='replace')

            write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
            write(VTK_file_unit,'(A)')'VTKFile'
            write(VTK_file_unit,'(A)')'ASCII'
            write(VTK_file_unit,'(A)')
            write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

      write(VTK_file_unit,'(A7,I12,A7)')'POINTS ',n_elem*4,' double';
            POINT_LOOP: do i=1,n_elem
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              do j=1,4
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,i),yy(j,i),zz(j,i)
              end do
            end do POINT_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem, n_elem*(nnode_per_el+1);
            ELEM_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I12,8I12)')nnode_per_el, tnew(i,:);
              !write(VTK_file_unit,'(I12,8I12)')nnode_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
              !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1;
            end do ELEM_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem;
            ELEM_TYPE_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I2)')10;
            end do ELEM_TYPE_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem;
      write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','mpi_id', ' int 1'
      write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'
            POLY_LOOP2: do i=1,n_elem
              write(VTK_file_unit,'(I6)') PolyMesh%part_elem(PolyMesh%elem_loc2glo(i))
            end do POLY_LOOP2

      write(VTK_file_unit,'(A)')

      close(unit=VTK_file_unit)

    end subroutine VTK_WRITE_MESH_PARTITION

    subroutine VTK_WRITE_MESH_AGGLOMERATION(filename, xx, yy, zz, nnode_per_el, n_elem, PolyMesh)

      implicit none

      ! Input arguments
      character(len=*), intent(in) :: filename
      integer*4, intent(in) :: n_elem,nnode_per_el
      real*8, dimension(4,n_elem), intent(in) :: xx,yy,zz
      integer*4,dimension(n_elem,4) :: tnew
      type(Mesh_Structure), intent(in) :: PolyMesh
      ! integer*4,dimension(n_elem) :: E2P

      ! Internal variables
      integer*4 :: VTK_file_unit
      integer*4 :: i,j

      ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
      real*8, parameter :: cap = 1E-40

      do i=1,n_elem
        do j=1,4
          tnew(i,j)=(i-1)*4+j-1
        end do
      end do

      open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', status='replace')

            write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
            write(VTK_file_unit,'(A)')'VTKFile'
            write(VTK_file_unit,'(A)')'ASCII'
            write(VTK_file_unit,'(A)')
            write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

      write(VTK_file_unit,'(A7,I12,A7)')'POINTS ',n_elem*4,' double';
            POINT_LOOP: do i=1,n_elem
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              do j=1,4
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,i),yy(j,i),zz(j,i)
              end do
            end do POINT_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem, n_elem*(nnode_per_el+1);
            ELEM_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I12,8I12)')nnode_per_el, tnew(i,:);
              !write(VTK_file_unit,'(I12,8I12)')nnode_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
              !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1;
            end do ELEM_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem;
            ELEM_TYPE_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I2)')10;
            end do ELEM_TYPE_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem;
      write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','poly_id', ' int 1'
      write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'
            POLY_LOOP3: do i=1,n_elem
              write(VTK_file_unit,'(I6)') PolyMesh%elem_in_poly_loc(i)
            end do POLY_LOOP3

      write(VTK_file_unit,'(A)')

      close(unit=VTK_file_unit)

    end subroutine VTK_WRITE_MESH_AGGLOMERATION

    !> Store the errors, the degree and hmax in an appropriate file
    subroutine VTK_WRITE_ERR(filename, p, err_L2, err_DG, hh)

        implicit none

        ! Input arguments
        character(len=*), intent(in) :: filename
        integer (kind=4) :: p
        real (kind=8) :: err_L2
        real (kind=8) :: err_DG
        real (kind=8) :: hh

        ! Internal variables
        integer*4 :: VTK_file_unit
        !integer*4 :: i,j

        open(newunit=VTK_file_unit, action='WRITE', file=filename, &
              form='FORMATTED', status='replace')

        write(VTK_file_unit,'(A)')'degree';
        write(VTK_file_unit,'(I1)') p;

        write(VTK_file_unit,'(A)')'err_L2';
              !ERRL2_LOOP: do i=1,nit
              !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              write(VTK_file_unit,'(1F16.8)') err_L2
              !end do ERRL2_LOOP

        write(VTK_file_unit,'(A)')'err_DG';
              !ERRH1_LOOP: do i=1,nit
              !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              write(VTK_file_unit,'(1F16.8)') err_DG
              !end do ERRH1_LOOP

        write(VTK_file_unit,'(A)')'h';
              !ELEM_LOOP: do i=1,nit
                !write(VTK_file_unit,'(I12,8I12)')nnode_per_el, t(i,:)-1;
              write(VTK_file_unit,'(1F16.8)') hh
              !end do ELEM_LOOP

        close(unit=VTK_file_unit)

    end subroutine

    !> Store the errors, the degree and hmax in an appropriate file for the convergence test
    !> It is similar to VTK_WRITE_ERR, but here we append the data instead of replacing
    subroutine VTK_WRITE_CONV(filename, p, err_L2, err_DG, hh)

      implicit none

      ! Input arguments
      character(len=*), intent(in) :: filename
      integer (kind=4) :: p
      real (kind=8) :: err_L2
      real (kind=8) :: err_DG
      real (kind=8) :: hh

      ! Internal variables
      integer*4 :: VTK_file_unit
      !integer*4 :: i,j

      open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', position='append', status='unknown')

      write(VTK_file_unit,'(A)')'degree';
      write(VTK_file_unit,'(I1)') p;

      write(VTK_file_unit,'(A)')'err_L2';
            !ERRL2_LOOP: do i=1,nit
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            write(VTK_file_unit,'(1F16.8)') err_L2
            !end do ERRL2_LOOP

      write(VTK_file_unit,'(A)')'err_DG';
            !ERRH1_LOOP: do i=1,nit
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            write(VTK_file_unit,'(1F16.8)') err_DG
            !end do ERRH1_LOOP

      write(VTK_file_unit,'(A)')'h';
            !ELEM_LOOP: do i=1,nit
              !write(VTK_file_unit,'(I12,8I12)')nnode_per_el, t(i,:)-1;
            write(VTK_file_unit,'(1F16.8)') hh
            !end do ELEM_LOOP
      write(VTK_file_unit,'(A)') ' ';

      close(unit=VTK_file_unit)

    end subroutine

    !> Store the numerical solution in an appropriate file
    !! HARDCODED 4 FOR 4 VERTICES OF TET
    !! MERGE SUBROUTINES...
    subroutine WRITE_SOLUTION_VTK(n_elem, PolyMesh, u, IsPoly, num_dt)

        use problem_data_and_properties
        use mpi
        use Poly_setup_MPI

        implicit none

        integer(kind=4), intent(in), optional :: num_dt       !< number of timesteps
        type(Mesh_Structure), intent(inout) :: PolyMesh       !< mesh
        integer(kind=4), intent(in) :: n_elem                  !< local number of elements
        real(kind=8), dimension(3,4,n_elem), intent(inout) :: u !< 3D, 4 vertices, n_elem
        logical, intent(in) :: IsPoly
        real(kind=8), dimension(4,n_elem) :: xx, yy, zz
        character(len=80) :: vtk_filename_num!, vtk_filename_exact
        integer(kind=4) :: ie_loc,ivert,id_node

        do ie_loc = 1,n_elem

            do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert

                call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                    PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)

                xx(ivert,ie_loc)=PolyMesh%coord_x(id_node)
                yy(ivert,ie_loc)=PolyMesh%coord_y(id_node)
                zz(ivert,ie_loc)=PolyMesh%coord_z(id_node)

            enddo

        enddo

        ! WRITE VTK-FILE ----------------------------------------------------
        if(IsPoly) then

            ! write(vtk_filename, '(A,I0,A)') 'MONITORS/sol_poly_', PolyMesh%num_poly, '.vtk'
            vtk_filename_num = 'MONITORS/sol_poly_000000.vtk'
            if (mpi_id < 10) then
                write(vtk_filename_num(24:24),'(i1)') mpi_id
            elseif (mpi_id < 100) then
                write(vtk_filename_num(23:24),'(i2)') mpi_id
            elseif (mpi_id < 1000) then
                write(vtk_filename_num(22:24),'(i3)') mpi_id
            elseif (mpi_id < 10000) then
                write(vtk_filename_num(21:24),'(i4)') mpi_id
            elseif (mpi_id < 100000) then
                write(vtk_filename_num(20:24),'(i5)') mpi_id
            elseif (mpi_id < 1000000) then
                write(vtk_filename_num(19:24),'(i6)') mpi_id
            endif

        else

            ! write(vtk_filename, '(A,I0,A)') 'MONITORS/sol_tet_', PolyMesh%num_tet, '.vtk'
            vtk_filename_num = 'MONITORS/sol_tet_000000.vtk'

            ! num_dt
            if (num_dt < 10) then
                write(vtk_filename_num(23:23),'(i1)') num_dt
            elseif (num_dt < 100) then
                write(vtk_filename_num(22:23),'(i2)') num_dt
            elseif (num_dt < 1000) then
                write(vtk_filename_num(21:23),'(i3)') num_dt
            elseif (num_dt < 10000) then
                write(vtk_filename_num(20:23),'(i4)') num_dt
            elseif (num_dt < 100000) then
                write(vtk_filename_num(19:23),'(i5)') num_dt
            elseif (num_dt < 1000000) then
                write(vtk_filename_num(18:23),'(i6)') num_dt
            endif

            ! mpi id
            ! if (mpi_id < 10) then
            !   write(vtk_filename_num(30:30),'(i1)') mpi_id
            ! elseif (mpi_id < 100) then
            !   write(vtk_filename_num(29:30),'(i2)') mpi_id
            ! elseif (mpi_id < 1000) then
            !   write(vtk_filename_num(28:30),'(i3)') mpi_id
            ! elseif (mpi_id < 10000) then
            !   write(vtk_filename_num(27:30),'(i4)') mpi_id
            ! elseif (mpi_id < 100000) then
            !   write(vtk_filename_num(26:30),'(i5)') mpi_id
            ! elseif (mpi_id < 1000000) then
            !   write(vtk_filename_num(25:30),'(i6)') mpi_id
            ! endif

        endif

        call VTK_WRITE_SOLUTION(vtk_filename_num, xx,yy,zz, 4, n_elem, PolyMesh%num_elem, 'solution', u, PolyMesh)

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call VTU_WRITE_SOLUTION(vtk_filename_num, xx,yy,zz, 4, n_elem, PolyMesh%num_elem, 'solution', u, PolyMesh, num_dt)

    end subroutine WRITE_SOLUTION_VTK

    !> Write .vtk files for the visualization of the partition and agglomeration of the mesh
    subroutine WRITE_MESH_VISUALIZATION_VTK(n_elem, PolyMesh, mpi_id)

      implicit none

      type(Mesh_Structure), intent(inout) :: PolyMesh
      integer(kind=4), intent(in) :: n_elem
      integer(kind=4), intent(in) :: mpi_id
      real(kind=8), dimension(4,n_elem) :: xx, yy, zz
      character(len=80) :: vtk_filename_partition, vtk_filename_agglomeration
      integer(kind=4) :: ie_loc,ivert,id_node

      do ie_loc = 1,n_elem

        do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert

            call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                  PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)

            xx(ivert,ie_loc)=PolyMesh%coord_x(id_node)
            yy(ivert,ie_loc)=PolyMesh%coord_y(id_node)
            zz(ivert,ie_loc)=PolyMesh%coord_z(id_node)

        enddo

      enddo

      ! WRITE VTK-FILE ----------------------------------------------------

        vtk_filename_partition = 'mesh_visualization/mesh_partition_000000.vtk'
        if (mpi_id < 10) then
          write(vtk_filename_partition(40:40),'(i1)') mpi_id
        elseif (mpi_id < 100) then
          write(vtk_filename_partition(39:40),'(i2)') mpi_id
        elseif (mpi_id < 1000) then
          write(vtk_filename_partition(38:40),'(i3)') mpi_id
        elseif (mpi_id < 10000) then
          write(vtk_filename_partition(37:40),'(i4)') mpi_id
        elseif (mpi_id < 100000) then
          write(vtk_filename_partition(36:40),'(i5)') mpi_id
        elseif (mpi_id < 1000000) then
          write(vtk_filename_partition(35:40),'(i6)') mpi_id
        endif

        vtk_filename_agglomeration = 'mesh_visualization/mesh_agglomeration_000000.vtk'
        if (mpi_id < 10) then
          write(vtk_filename_agglomeration(44:44),'(i1)') mpi_id
        elseif (mpi_id < 100) then
          write(vtk_filename_agglomeration(43:44),'(i2)') mpi_id
        elseif (mpi_id < 1000) then
          write(vtk_filename_agglomeration(42:44),'(i3)') mpi_id
        elseif (mpi_id < 10000) then
          write(vtk_filename_agglomeration(41:44),'(i4)') mpi_id
        elseif (mpi_id < 100000) then
          write(vtk_filename_agglomeration(40:44),'(i5)') mpi_id
        elseif (mpi_id < 1000000) then
          write(vtk_filename_agglomeration(39:44),'(i6)') mpi_id
        endif

      call VTK_WRITE_MESH_PARTITION(vtk_filename_partition, xx,yy,zz, 4, n_elem, PolyMesh)
      call VTK_WRITE_MESH_AGGLOMERATION(vtk_filename_agglomeration, xx,yy,zz, 4, n_elem, PolyMesh)

      return

    end subroutine WRITE_MESH_VISUALIZATION_VTK

    !> Associate files storing the errors in an appropriate name
    !> and call the functions which stores the data
    subroutine WRITE_ERRORS(p, err_DG, err_L2, hh, PolyMesh, IsPoly)

      use problem_data_and_properties

      implicit none

      type(Mesh_Structure), intent(in) :: PolyMesh
      integer(kind=4), intent(in) :: p
      real(kind=8), intent(in) :: err_L2
      real(kind=8), intent(in) :: err_DG
      real(kind=8), intent(in) :: hh
      logical, intent(in) :: IsPoly

      character(len=80) :: filename_err, filename_conv
      character(len=80) :: file_poly, file_tet
      real(kind=8) :: alpha, theta, c

      call set_properties(alpha, theta, c)

      ! WRITE VTK-FILE ----------------------------------------------------
      if(IsPoly) then

        write(file_poly, '(A,I0,A)') 'ERRORS_POLY/Errors_poly_', PolyMesh%num_poly, '.vtk'

        if(c /= 0) then
          filename_err = 'POST-PROC/DIFF_REAC_EQ/' // file_poly
          filename_conv = 'POST-PROC/DIFF_REAC_EQ/CONVERGENCE_TEST_POLY/Convergence_test.vtk'
        else
          filename_err = 'POST-PROC/DIFF_EQ/' // file_poly
          filename_conv = 'POST-PROC/DIFF_EQ/CONVERGENCE_TEST_POLY/Convergence_test.vtk'
        endif

      else

        write(file_tet, '(A,I0,A)') 'ERRORS_TET/Errors_tet_', PolyMesh%num_tet, '.vtk'

        if(c /= 0) then
          filename_err = 'POST-PROC/DIFF_REAC_EQ/' // file_tet
          filename_conv = 'POST-PROC/DIFF_REAC_EQ/CONVERGENCE_TEST_TET/Convergence_test.vtk'
        else
          filename_err = 'POST-PROC/DIFF_EQ/' // file_tet
          filename_conv = 'POST-PROC/DIFF_EQ/CONVERGENCE_TEST_TET/Convergence_test.vtk'
        endif

      endif

      call VTK_WRITE_ERR(filename_err, p, err_L2, err_DG, hh)
      call VTK_WRITE_CONV(filename_conv, p, err_L2, err_DG, hh)

      return

    end subroutine WRITE_ERRORS

end module MOD_VTK