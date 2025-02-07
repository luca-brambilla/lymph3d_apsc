!> Module containing type definitions and routines for real-time VTK output during computations for individual material blocks
module export_file_formats

    use Poly_mesh
    use global_parameters

    implicit none

    contains

    !! change types?
    !! hardcoded 3 for 3D and 4 for tetrahedra vertices
    !! hardcoded VTK type for tetrahedra 10

    !> Export solution in .vtk format
subroutine VTK_WRITE_SOLUTION(filename, xx, yy, zz, nvert_per_el, n_elem, n_elem_tot, u_name, u, PolyMesh)

    use mpi
    use Poly_setup_MPI

    implicit none

    ! Input arguments
    character(len=*), intent(in) :: filename                !< string file name
    integer*4, intent(in) :: n_elem                     !< local number elements for the process
    integer*4, intent(in) :: nvert_per_el               !< number of vertices per element
    real*8, dimension(nvert_per_el, n_elem), intent(in) :: xx       !< x coordinates of 4 vertices of tetrahedron
    real*8, dimension(nvert_per_el, n_elem), intent(in) ::yy        !< y coordinates of 4 vertices of tetrahedron
    real*8, dimension(nvert_per_el, n_elem), intent(in) ::zz       !< z coordinates of 4 vertices of tetrahedron
    integer*4, dimension(n_elem_tot,4) :: tnew
    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh

    ! Optional input arguments
    character(len=*), intent(in), optional :: u_name            !< solution name
    real*8, dimension(DIM, nvert_per_el, n_elem), intent(inout), optional :: u !< solution

    ! Internal variables
    integer*4 :: VTK_file_unit
    integer*4 :: j,k,d,ie

    ! added
    integer*4, intent(in) :: n_elem_tot

    integer*4 :: tmp_nelem
    real*8, dimension(:,:,:), allocatable :: tmp_solution
    real*8, dimension(:,:), allocatable :: tmp_xx, tmp_yy, tmp_zz
    integer*4, dimension(:), allocatable :: tmp_poly
    !integer(kind=4), intent(in) :: mpi_id

    ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
    real*8, parameter :: cap = 1d-40

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
        call MPI_Send(xx, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send solution
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u, n_elem * nvert_per_el * DIM, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
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
        write(VTK_file_unit,'(A7,I12,A7)')'POINTS ', n_elem_tot*nvert_per_el,' double'
        POINT_LOOP: do ie=1,n_elem
        !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            do j=1,nvert_per_el
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,ie),yy(j,ie),zz(j,ie)
            end do
        end do POINT_LOOP

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive coordinates - different sizes
            allocate(tmp_xx(nvert_per_el, tmp_nelem))
            allocate(tmp_yy(nvert_per_el, tmp_nelem))
            allocate(tmp_zz(nvert_per_el, tmp_nelem))
            call MPI_Recv(tmp_xx, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_yy, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_zz, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write to file
            POINT_LOOP_k: do ie=1,tmp_nelem
                do j=1,nvert_per_el
                    write(VTK_file_unit,'(3F16.8,1X)') tmp_xx(j,ie),tmp_yy(j,ie),tmp_zz(j,ie)
                end do
            end do POINT_LOOP_k
            ! deallocate coordinates, different sizes
            deallocate(tmp_xx, tmp_yy, tmp_zz)
        enddo
        write(VTK_file_unit,'(A)')

        ! **************
        ! CELLS
        ! **************
        write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem_tot, n_elem_tot*(nvert_per_el+1)
        ELEM_LOOP: do ie=1,n_elem_tot
            ! ordering
            do j=1,nvert_per_el
                tnew(ie,j) = (ie-1)*nvert_per_el +j-1
            end do
            ! write ordering
            write(VTK_file_unit,'(I12,8I12)')nvert_per_el, tnew(ie,:)
            !write(VTK_file_unit,'(I12,8I12)')nvert_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
            !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1
        end do ELEM_LOOP
        write(VTK_file_unit,'(A)')

        ! ******************
        ! CELL TYPES
        ! ******************
        !! hardcoded tetrahedra
        write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem_tot
            ELEM_TYPE_LOOP: do ie=1,n_elem_tot
                write(VTK_file_unit,'(I2)')10
            end do ELEM_TYPE_LOOP
        write(VTK_file_unit,'(A)')

        ! *********************
        ! SOLUTION
        ! *********************
        if (present(u_name) .and. present(u)) then
            ! process 0 data
            write(VTK_file_unit,'(A11,I12)')'POINT_DATA ', n_elem_tot * nvert_per_el
            write(VTK_file_unit,'(A8,A12,A7)')'VECTORS ', u_name, ' double'
            !write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'
            VECTOR_FIELD_LOOP: do ie=1,n_elem
                do j=1,nvert_per_el
                    write(VTK_file_unit,'(3E20.10)') u(1,j,ie), u(2,j,ie), u(3,j,ie)
                    !write(VTK_file_unit,'(3F16.8)') u(1:3,j,i)
                enddo
            end do VECTOR_FIELD_LOOP

            ! mpi loop
            do k=1,mpi_np-1
                ! receive size
                call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
                ! allocate and receive solution - different sizes
                allocate(tmp_solution(DIM, nvert_per_el, tmp_nelem))
                call MPI_Recv(tmp_solution, tmp_nelem*nvert_per_el*DIM, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
                ! write solution
                VECTOR_FIELD_LOOP_k: do ie=1,tmp_nelem
                    do j=1,nvert_per_el
                        write(VTK_file_unit,'(3E20.10)') ( tmp_solution(d,j,ie), d=1,DIM) !, tmp_solution(2,j,i), tmp_solution(3,j,i)
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
        write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem_tot
        write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','POLYHEDRA', ' int 1'
        write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'

        POLY_LOOP1: do ie=1,n_elem
            write(VTK_file_unit,'(I10)') PolyMesh%elem_in_poly_loc(ie)
        end do POLY_LOOP1

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive polyhedra - different sizes
            allocate(tmp_poly(tmp_nelem))
            call MPI_Recv(tmp_poly, tmp_nelem, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write polyhedra
            POLY_LOOP_k: do ie=1,tmp_nelem
                write(VTK_file_unit,'(I10)') tmp_poly(ie)
            end do POLY_LOOP_k
            ! deallocate polyhedra - different sizes
            deallocate(tmp_poly)
        enddo

        write(VTK_file_unit,'(A)')

        close(unit=VTK_file_unit)
    endif

end subroutine VTK_WRITE_SOLUTION

!> Export solution in .pdv format (joint with VTU_WRITE_SOLUTION)
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
!> Export solution in .vtu format (joint with PVD_SETUP)
subroutine VTU_WRITE_SOLUTION(filename, xx, yy, zz, nvert_per_el, n_elem, n_elem_tot, u_name, u, PolyMesh, num_dt)

    use mpi
    use Poly_setup_MPI

    implicit none

    ! Input arguments
    character(len=*), intent(in) :: filename                !< string file name
    integer*4, intent(in) :: n_elem                     !< local number of dofs for the process
    integer*4, intent(in) :: nvert_per_el               !< number of vertices per element
    real*8, dimension(nvert_per_el,n_elem), intent(in) :: xx       !< x coordinates of 4 vertices of tetrahedron
    real*8, dimension(nvert_per_el,n_elem), intent(in) ::yy        !< y coordinates of 4 vertices of tetrahedron
    real*8, dimension(nvert_per_el,n_elem), intent(in) ::zz       !< z coordinates of 4 vertices of tetrahedron
    !integer*4, dimension(n_elem_tot,4) :: tnew
    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh
    integer*4, intent(in) :: num_dt                            !< number of timestep

    ! Optional input arguments
    character(len=*), intent(in), optional :: u_name            !< solution name
    real*8, dimension(DIM,nvert_per_el,n_elem), intent(inout), optional :: u !< solution

    ! Internal variables
    integer*4 :: VTU_file_unit
    integer*4 :: i,j,k,d, node_id

    ! added
    integer*4, intent(in) :: n_elem_tot

    integer*4 :: tmp_nelem
    real*8, dimension(:,:,:), allocatable :: tmp_solution
    real*8, dimension(:,:), allocatable :: tmp_xx, tmp_yy, tmp_zz
    !integer*4, dimension(:), allocatable :: tmp_poly
    !integer(kind=4), intent(in) :: mpi_id

    character(len=256) :: mesh_filename, solution_filename
    integer*4 :: elem_offset, conn_offset

    ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
    real*8, parameter :: cap = 1d-40

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
        call MPI_Send(xx, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        ! send solution
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u, n_elem*nvert_per_el*DIM, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
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
        write(VTU_file_unit, '(A27,I0,A17,I0,A2)') '    <Piece NumberOfPoints="', n_elem_tot*nvert_per_el, '" NumberOfCells="', n_elem_tot, '">'

        ! ******
        ! POINTS Section
        write(VTU_file_unit, '(A)') '    <Points>'
        write(VTU_file_unit, '(A)') '      <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
        do i = 1, n_elem
            do j = 1, nvert_per_el
                write(VTU_file_unit, '(3(F16.8,1X))') xx(j, i), yy(j, i), zz(j, i)
            end do
        end do

        ! mpi loop
        do k=1,mpi_np-1

            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

            ! allocate and receive coordinates - different sizes
            allocate(tmp_xx(nvert_per_el,tmp_nelem))
            allocate(tmp_yy(nvert_per_el,tmp_nelem))
            allocate(tmp_zz(nvert_per_el,tmp_nelem))
            call MPI_Recv(tmp_xx, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_yy, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            call MPI_Recv(tmp_zz, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

            ! write to file from each process
            POINT_LOOP_k: do i=1,tmp_nelem
                do j=1,nvert_per_el
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
            do j = 1, nvert_per_el
                node_id = (i-1)*nvert_per_el + (j-1)
                write(VTU_file_unit, '(I12)', advance="no") node_id
                if (j < nvert_per_el) write(VTU_file_unit, '(A)', advance="no") " "
            end do
            write(VTU_file_unit, '(A)')
        end do
        write(VTU_file_unit, '(A)') '      </DataArray>'

        ! offset - Cumulative indices in the connectivity array
        write(VTU_file_unit, '(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
        conn_offset = 0
        do i = 1, n_elem_tot
            conn_offset = conn_offset + nvert_per_el
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
            write(VTU_file_unit, '(A27,I0,A2)') '    <Piece NumberOfPoints="', n_elem_tot*nvert_per_el, '">'

            ! POINT_DATA Section
            write(VTU_file_unit, '(A)') '    <PointData>'
            write(VTU_file_unit, '(A)') '      <DataArray type="Float64" Name="' // trim(u_name) // '" NumberOfComponents="3" format="ascii">'
            do i = 1, n_elem
                do j = 1, nvert_per_el
                    write(VTU_file_unit, '(3(F16.8,1X))') ( u(d, j, i), d=1,DIM )!, u(2, j, i), u(3, j, i)
                end do
            end do

            ! mpi loop
            do k=1,mpi_np-1

                ! receive size
                call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

                ! allocate and receive solution - different sizes
                allocate(tmp_solution(DIM,nvert_per_el,tmp_nelem))
                call MPI_Recv(tmp_solution, tmp_nelem*nvert_per_el*DIM, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)

                ! write solution
                VECTOR_FIELD_LOOP_k: do i=1,tmp_nelem
                    do j=1,nvert_per_el
                        write(VTU_file_unit,'(3(F16.8,1X))') ( tmp_solution(d,j,i), d=1,DIM )!, tmp_solution(2,j,i), tmp_solution(3,j,i)
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
        ! write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem_tot
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

!> receive data and print to file
subroutine RECEIVE_PRINT_FILE(file_unit, fmt, n_elem, nvert_per_el, xx, ii)

    use Poly_setup_MPI
    implicit none

    real(kind=8), dimension(:,:), optional, intent(in) :: xx
    integer(kind=4), dimension(:), optional, intent(in) :: ii
    integer(kind=4), intent(in) :: n_elem
    integer(kind=4), intent(in) :: nvert_per_el
    integer(kind=4), intent(in) :: file_unit

    real(kind=8), dimension(:,:), allocatable :: tmp_xx
    integer(kind=4), dimension(:), allocatable :: tmp_ii
    integer(kind=4) :: ie,j,tmp_nelem,k
    character(len=*), intent(in) :: fmt

    ! float
    if (fmt=='(E12.5)' ) then

        ! local
        do ie=1,n_elem
            do j=1,nvert_per_el
                write(file_unit, fmt) xx(j,ie)
            enddo
        enddo

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive data - different sizes
            allocate(tmp_xx(nvert_per_el, tmp_nelem))
            call MPI_Recv(tmp_xx, tmp_nelem*nvert_per_el, MPI_DOUBLE_PRECISION, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write to file
            do ie=1,tmp_nelem
                do j=1,nvert_per_el
                    write(file_unit, fmt) tmp_xx(j,ie)
                enddo
            enddo
            ! deallocate data, different sizes
            deallocate(tmp_xx)
        enddo

    ! integer
    elseif (fmt=='(I10)' ) then

        ! local
        do ie=1,n_elem
            write(file_unit, fmt) ii(ie)
        enddo

        ! mpi loop
        do k=1,mpi_np-1
            ! receive size
            call MPI_Recv(tmp_nelem, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! allocate and receive data - different sizes
            allocate(tmp_ii(tmp_nelem))
            call MPI_Recv(tmp_ii, tmp_nelem, MPI_INTEGER, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            ! write to file
            do ie=1,tmp_nelem
                write(file_unit, fmt) tmp_ii(ie)
            enddo
            ! deallocate data, different sizes
            deallocate(tmp_ii)
        enddo

    else
        print *, 'wrong format'
    endif

end subroutine RECEIVE_PRINT_FILE

!> writes solution file in ensight format
!> saves geometry and polyhedra into 2 separate files
!> can save geometry only once
!> FORTRAN ASCII floats E12.5, integers I10 one per line
subroutine ENSIGHT_WRITE_MESH(base_filename, xx, yy, zz, nvert_per_el, n_elem, n_elem_tot, PolyMesh)

    use mpi
    use Poly_setup_MPI

    !TODO save correct polyhedra and node numbering from input file
    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh

    ! Input arguments
    integer*4, intent(in) :: nvert_per_el             !< number of vertices per element
    integer*4, intent(in) :: n_elem,n_elem_tot                   !< Number of elements in this part
    real(8), intent(in) :: xx(nvert_per_el, n_elem)            !< x coordinates of 4 vertices of tetrahedra
    real(8), intent(in) :: yy(nvert_per_el, n_elem)            !< y coordinates of 4 vertices of tetrahedra
    real(8), intent(in) :: zz(nvert_per_el, n_elem)            !< z coordinates of 4 vertices of tetrahedra

    character(len=*), intent(in) :: base_filename

    ! Internal variables
    integer*4 :: i, j
    integer*4 :: geo_unit, poly_unit
    character(len=256) :: geo_filename, poly_filename

    if (mpi_id /= 0) then
        ! tag=1 before MPI_COMM_WORLD
        ! send size and coordinates
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(xx, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, n_elem*nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)

        ! poly
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(PolyMesh%elem_in_poly_loc, n_elem, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)

    elseif (mpi_id == 0) then

        write(geo_filename, '(A, "mesh.geo")') trim(base_filename)
        write(poly_filename, '(A, "poly.sca")') trim(base_filename)

        ! **********************
        ! Write Geometry File
        ! **********************
        open(newunit=geo_unit, file=geo_filename, action='WRITE', status='REPLACE', form='FORMATTED')

        ! EnSight header
        write(geo_unit, '(A)') 'Ensight Gold'
        write(geo_unit, '(A)') 'Unstructured mesh part'
        write(geo_unit, '(A)') 'node id assign'
        write(geo_unit, '(A)') 'element id assign'

        ! Write part header
        write(geo_unit, '(A)') 'part'
        write(geo_unit, '(I10)') 1
        write(geo_unit, '(A)') 'Static mesh file'

        ! Write coordinates
        write(geo_unit, '(A)') 'coordinates'
        write(geo_unit, '(I10)') n_elem_tot * nvert_per_el
        ! numbering first
        do i = 1, n_elem_tot * nvert_per_el
            write(geo_unit, '(I10)') i
        end do

        ! coordinates
        call RECEIVE_PRINT_FILE(geo_unit, '(E12.5)', n_elem, nvert_per_el, xx=xx)
        call RECEIVE_PRINT_FILE(geo_unit, '(E12.5)', n_elem, nvert_per_el, xx=yy)
        call RECEIVE_PRINT_FILE(geo_unit, '(E12.5)', n_elem, nvert_per_el, xx=zz)

        ! Write elements
        write(geo_unit, '(A)') 'tetra4' ! Use tetrahedral representation
        write(geo_unit, '(I10)') n_elem_tot
        do i = 1, n_elem_tot
            write(geo_unit, '(I10)') i
        end do
        do i = 1, n_elem_tot
            write(geo_unit, '(4I10)') ( (i - 1) * nvert_per_el + j, j = 1, nvert_per_el )
        end do

        close(geo_unit)


        open(newunit=poly_unit, file=poly_filename, action='WRITE', status='REPLACE', form='FORMATTED')

        ! Write solution header
        write(poly_unit, '(A)') 'scalar per element'
        write(poly_unit, '(A)') 'part'
        write(poly_unit, '(I10)') 1
        write(poly_unit, '(A)') 'tetra4'

        ! Write solution values
        call RECEIVE_PRINT_FILE(poly_unit, '(I10)', n_elem, nvert_per_el, ii=PolyMesh%elem_in_poly_loc)

        close(poly_unit)

    endif

end subroutine ENSIGHT_WRITE_MESH

subroutine ENSIGHT_WRITE_BOUNDARY(base_filename, nvert_per_el, n_elem, PolyMesh)

    use mpi
    use Poly_setup_MPI
    use mesh_partition_and_mpi_files

    !TODO save correct polyhedra and node numbering from input file
    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh

    ! Input arguments
    integer*4, intent(in) :: nvert_per_el             !< number of vertices per element
    integer*4, intent(in) :: n_elem                   !< Number of elements in this part

    character(len=*), intent(in) :: base_filename

    ! Internal variables
    integer*4 :: i, j
    integer*4 :: bd_unit
    character(len=256) :: bd_filename

    integer(kind=4) :: nbd, nbd_tot, E2, ie_loc, iface, mat_id, ivert, vert_id
    integer(kind=4), dimension(:), allocatable :: bd_faces
    real(kind=8), dimension(:,:), allocatable :: xx, yy, zz

    nbd = 0

    ! count boundaries
    do ie_loc = 1,n_elem
        do iface=1,PolyMesh%Elem_loc(1)%num_faces

            E2 = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,2)
            if (E2<0) then
                nbd = nbd +1
            endif

        enddo
    enddo

    allocate(bd_faces(nbd))
    allocate(xx(nvert_per_el,nbd), yy(nvert_per_el,nbd), zz(nvert_per_el,nbd))


    ! store boundary face id
    nbd = 0
    do ie_loc = 1,n_elem
        do iface=1,PolyMesh%Elem_loc(1)%num_faces

            E2 = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,2)
            if (E2<0) then
                ! counter
                nbd = nbd+1

                mat_id = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1)
                bd_faces(nbd) = mat_id

                ! row on connectivity
                !face_id = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,3)
                !row_vert = (ie_loc-1)*4 + face_id

                do ivert=1,nvert_per_el

                    call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                    PolyMesh%Elem_loc(ie_loc)%vert(ivert),vert_id)

                    ! global vertex id
                    ! vert_id = PolyMesh%con_tria(row_vert,3+ivert)
    
                    ! vertex coordinates
                    xx(ivert,nbd) = PolyMesh%coord_x(vert_id)
                    yy(ivert,nbd) = PolyMesh%coord_y(vert_id)
                    zz(ivert,nbd) = PolyMesh%coord_z(vert_id)
                enddo
            endif

        enddo
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call MPI_REDUCE(nbd, nbd_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id /= 0) then
        ! tag=1 before MPI_COMM_WORLD
        ! send size and coordinates
        call MPI_Send(nbd, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(xx, nvert_per_el*nbd, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(nbd, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(yy, nvert_per_el*nbd, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(nbd, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(zz, nvert_per_el*nbd, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)

        ! poly
        call MPI_Send(nbd, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(bd_faces, nbd, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)

    elseif (mpi_id == 0) then

        write(bd_filename, '(A, "boundaries.case")') trim(base_filename)

        ! **********************
        ! Write Geometry File
        ! **********************
        open(newunit=bd_unit, file=bd_filename, action='WRITE', status='REPLACE', form='FORMATTED')

        ! EnSight header
        write(bd_unit, '(A)') 'Ensight Gold'
        write(bd_unit, '(A)') 'Unstructured mesh part'
        write(bd_unit, '(A)') 'node id assign'
        write(bd_unit, '(A)') 'element id assign'

        ! Write part header
        write(bd_unit, '(A)') 'part'
        write(bd_unit, '(I10)') 1
        write(bd_unit, '(A)') 'boundary faces'

        ! Write coordinates
        write(bd_unit, '(A)') 'coordinates'
        write(bd_unit, '(I10)') nbd_tot * nvert_per_el
        ! numbering first
        do i = 1, nbd_tot * nvert_per_el
            write(bd_unit, '(I10)') i
        end do

        ! coordinates
        call RECEIVE_PRINT_FILE(bd_unit, '(E12.5)', nbd, NVERT_TRIA, xx=xx)
        call RECEIVE_PRINT_FILE(bd_unit, '(E12.5)', nbd, NVERT_TRIA, xx=yy)
        call RECEIVE_PRINT_FILE(bd_unit, '(E12.5)', nbd, NVERT_TRIA, xx=zz)

        ! Write elements
        write(bd_unit, '(A)') 'tria3' ! Use tetrahedral representation
        write(bd_unit, '(I10)') nbd_tot
        do i = 1, nbd_tot
            write(bd_unit, '(I10)') i
        end do
        do i = 1, nbd_tot
            write(bd_unit, '(3I10)') ( (i - 1) * nvert_per_el + j, j = 1, NVERT_TRIA )
        end do

        write(bd_unit, '(A)') "scalar per element"
        write(bd_unit, '(A)') "part"
        write(bd_unit, '(I10)') 1
        write(bd_unit, '(A)') 'tria3'

        ! Write solution values
        call RECEIVE_PRINT_FILE(bd_unit, '(I10)', nbd, NVERT_TRIA, ii=bd_faces)

        close(bd_unit)

    endif

end subroutine ENSIGHT_WRITE_BOUNDARY

!> writes solution file in ensight format
!> FORTRAN ASCII floats E12.5, integers I10 one per line
subroutine ENSIGHT_WRITE_SOLUTION(base_filename, nvert_per_el, n_elem, n_elem_tot, u_name, u, num_dt, PolyMesh)

    use mpi
    use Poly_setup_MPI

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh

    ! Input arguments
    integer*4, intent(in) :: nvert_per_el             !< number of vertices per element
    integer*4, intent(in) :: n_elem,n_elem_tot                   !< Number of elements in this part
    character(len=*), intent(in), optional :: u_name !< Solution name (e.g., "Velocity")
    real(8), dimension(DIM, nvert_per_el, n_elem), intent(in), optional :: u !< Solution values (3D vectors per node)
    integer*4, intent(in) :: num_dt                            !< number of timestep
    character(len=*), intent(in) :: base_filename

    real(8), dimension(:,:), pointer :: u_tmp

    ! Internal variables
    integer*4 :: sol_unit
    character(len=256) :: sol_filename, tmp_filename

    allocate(u_tmp(nvert_per_el, n_elem))

    if (mpi_id /= 0) then

        ! send size each time because receive from all processes one dimension at a time
        u_tmp = u(1,:,:)
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u_tmp, n_elem * nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)

        u_tmp = u(2,:,:)
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u_tmp, n_elem * nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)

        u_tmp = u(3,:,:)
        call MPI_Send(n_elem, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_ierr)
        call MPI_Send(u_tmp, n_elem * nvert_per_el, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpi_ierr)


    elseif (mpi_id == 0) then

        tmp_filename = 'sol_tet_000000.vec'

        ! num_dt
        if (num_dt < 10) then
            write(tmp_filename(14:14),'(i1)') num_dt
        elseif (num_dt < 100) then
            write(tmp_filename(13:14),'(i2)') num_dt
        elseif (num_dt < 1000) then
            write(tmp_filename(12:14),'(i3)') num_dt
        elseif (num_dt < 10000) then
            write(tmp_filename(11:14),'(i4)') num_dt
        elseif (num_dt < 100000) then
            write(tmp_filename(10:14),'(i5)') num_dt
        elseif (num_dt < 1000000) then
            write(tmp_filename(9:14),'(i6)') num_dt
        endif

        write(sol_filename, '(A, A)') trim(base_filename), trim(tmp_filename)

        ! **********************
        ! Write Solution File
        ! **********************
        if (present(u_name) .and. present(u)) then
            open(newunit=sol_unit, file=sol_filename, action='WRITE', status='REPLACE', form='FORMATTED')

            ! Write solution header
            write(sol_unit, '(A)') 'vector per node'
            write(sol_unit, '(A)') 'part'
            write(sol_unit, '(I10)') 1
            write(sol_unit, '(A)') 'coordinates'

            ! Write solution values
            call RECEIVE_PRINT_FILE(sol_unit, '(E12.5)', n_elem, nvert_per_el, xx=u(1,:,:))
            call RECEIVE_PRINT_FILE(sol_unit, '(E12.5)', n_elem, nvert_per_el, xx=u(2,:,:))
            call RECEIVE_PRINT_FILE(sol_unit, '(E12.5)', n_elem, nvert_per_el, xx=u(3,:,:))

            close(sol_unit)
        end if

    endif

end subroutine ENSIGHT_WRITE_SOLUTION

!> write solution in ensight format
!> saves the case file to reference the solution at each timestep, with the same mesh
subroutine ENSIGHT_WRITE_CASE(base_filename, num_dt, solution_name)

    use mpi
    use Poly_setup_MPI
    use Poly_global

    implicit none

    ! Input arguments
    character(len=*), intent(in) :: base_filename  !< Base filename (without extension)
    integer*4, intent(in) :: num_dt
    character(len=*), intent(in), optional :: solution_name !< Solution variable name

    ! Internal variables
    character(len=256) :: case_filename!, tmp_name
    integer*4 :: case_file_unit, iostat, i, ndt

    ! Generate case filename
    write(case_filename, '(A, "solution.case")') trim(base_filename)

    open(newunit=case_file_unit, file=case_filename, action='WRITE', status='REPLACE', form='FORMATTED', iostat=iostat)

    if (iostat /= 0) then
        write(*,*) 'Error opening/writing file for process:', mpi_id
        stop
    end if

    write(case_file_unit, '(A)') 'FORMAT'
    write(case_file_unit, '(A)') 'type: ensight gold'

    ! Write geometry section
    write(case_file_unit, '(A)') 'GEOMETRY'
    write(case_file_unit, '(A)') 'model: mesh.geo'

    ! Write solution section if needed
    if (present(solution_name)) then
        write(case_file_unit, '(A)') 'VARIABLE'
        write(case_file_unit, '(A, I10, A)') 'vector per node: ', 1, ' u sol_tet_******.vec'
        write(case_file_unit, '(A, I10, A)') 'scalar per element: ', 1, ' poly poly.sca'

        ndt = floor(stop_time / time_step)/num_dt_mon+1

        ! Write the TIME section
        write(case_file_unit, '(A)') "TIME"
        write(case_file_unit, '(A, I10)') "time set: ", 1
        write(case_file_unit, '(A, I10)') "number of steps: ", ndt
        write(case_file_unit, '(A, I10)') "filename start number: ", 0
        write(case_file_unit, '(A, I10)') "filename increment: ", num_dt_mon

        ! Write time values in a single line
        write(case_file_unit, '(A)', advance='yes') "time values:"
        do i = 0, ndt-1
            if (mod(i,10)==0 .and. i/=0) write(case_file_unit, '(A)', advance='yes') ''
            if (i == ndt) then
                write(case_file_unit, '(F10.4)', advance='yes') i*num_dt_mon*time_step  ! Last value, finish the line
            else
                write(case_file_unit, '(F10.4)', advance='no') i*num_dt_mon*time_step  ! Continue on the same line
            end if
        end do

    end if

    close(case_file_unit)
end subroutine ENSIGHT_WRITE_CASE

!> write mesh partition files, only geometric data
    subroutine VTK_WRITE_MESH_PARTITION(filename, xx, yy, zz, nvert_per_el, n_elem, PolyMesh)

      implicit none

      ! Input arguments
      character(len=*), intent(in) :: filename
      integer*4, intent(in) :: n_elem,nvert_per_el
      real*8, dimension(nvert_per_el,n_elem), intent(in) :: xx,yy,zz
      integer*4,dimension(n_elem,nvert_per_el) :: tnew
      type(Mesh_Structure), intent(in) :: PolyMesh

      ! Internal variables
      integer*4 :: VTK_file_unit
      integer*4 :: i,j

      ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
      real*8, parameter :: cap = 1d-40

      do i=1,n_elem
        do j=1,nvert_per_el
          tnew(i,j)=(i-1)*nvert_per_el+j-1
        end do
      end do

      open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', status='replace')

            write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
            write(VTK_file_unit,'(A)')'VTKFile'
            write(VTK_file_unit,'(A)')'ASCII'
            write(VTK_file_unit,'(A)')
            write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

      write(VTK_file_unit,'(A7,I12,A7)')'POINTS ',n_elem*nvert_per_el,' double'
            POINT_LOOP: do i=1,n_elem
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              do j=1,nvert_per_el
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,i),yy(j,i),zz(j,i)
              end do
            end do POINT_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem, n_elem*(nvert_per_el+1)
            ELEM_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I12,8I12)')nvert_per_el, tnew(i,:)
              !write(VTK_file_unit,'(I12,8I12)')nvert_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
              !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1
            end do ELEM_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem
            ELEM_TYPE_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I2)')10
            end do ELEM_TYPE_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem
      write(VTK_file_unit,'(A8,A12,A9)')'SCALARS ','mpi_id', ' int 1'
      write(VTK_file_unit,'(A20)')'LOOKUP_TABLE default'
            POLY_LOOP2: do i=1,n_elem
              write(VTK_file_unit,'(I6)') PolyMesh%part_elem(PolyMesh%elem_loc2glo(i))
            end do POLY_LOOP2

      write(VTK_file_unit,'(A)')

      close(unit=VTK_file_unit)

    end subroutine VTK_WRITE_MESH_PARTITION

    subroutine VTK_WRITE_MESH_AGGLOMERATION(filename, xx, yy, zz, nvert_per_el, n_elem, PolyMesh)

      implicit none

      ! Input arguments
      character(len=*), intent(in) :: filename
      integer*4, intent(in) :: n_elem,nvert_per_el
      real*8, dimension(nvert_per_el,n_elem), intent(in) :: xx,yy,zz
      integer*4,dimension(n_elem,nvert_per_el) :: tnew
      type(Mesh_Structure), intent(in) :: PolyMesh
      ! integer*4,dimension(n_elem) :: E2P

      ! Internal variables
      integer*4 :: VTK_file_unit
      integer*4 :: i,j

      ! Capping parameter (VTK format problems with e.g. 1E-300 --> set to zero)
      real*8, parameter :: cap = 1d-40

      do i=1,n_elem
        do j=1,nvert_per_el
          tnew(i,j)=(i-1)*nvert_per_el+j-1
        end do
      end do

      open(newunit=VTK_file_unit, action='WRITE', file=filename, &
            form='FORMATTED', status='replace')

            write(VTK_file_unit,'(A)')'# vtk DataFile Version 3.0'
            write(VTK_file_unit,'(A)')'VTKFile'
            write(VTK_file_unit,'(A)')'ASCII'
            write(VTK_file_unit,'(A)')
            write(VTK_file_unit,'(A)')'DATASET UNSTRUCTURED_GRID'

      write(VTK_file_unit,'(A7,I12,A7)')'POINTS ',n_elem*nvert_per_el,' double'
            POINT_LOOP: do i=1,n_elem
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              do j=1,nvert_per_el
                write(VTK_file_unit,'(3F16.8,1X)') xx(j,i),yy(j,i),zz(j,i)
              end do
            end do POINT_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A6,I12,I12)')'CELLS ', n_elem, n_elem*(nvert_per_el+1)
            ELEM_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I12,8I12)')nvert_per_el, tnew(i,:)
              !write(VTK_file_unit,'(I12,8I12)')nvert_per_el, PolyMesh%con_tet(i,2)-1,PolyMesh%con_tet(i,3)-1,&
              !                                                PolyMesh%con_tet(i,4)-1,PolyMesh%con_tet(i,5)-1
            end do ELEM_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_TYPES ', n_elem
            ELEM_TYPE_LOOP: do i=1,n_elem
              write(VTK_file_unit,'(I2)')10
            end do ELEM_TYPE_LOOP
      write(VTK_file_unit,'(A)')

      write(VTK_file_unit,'(A11,I12)')'CELL_DATA ', n_elem
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

        write(VTK_file_unit,'(A)')'degree'
        write(VTK_file_unit,'(I1)') p

        write(VTK_file_unit,'(A)')'err_L2'
              !ERRL2_LOOP: do i=1,nit
              !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              write(VTK_file_unit,'(1F16.8)') err_L2
              !end do ERRL2_LOOP

        write(VTK_file_unit,'(A)')'err_DG'
              !ERRH1_LOOP: do i=1,nit
              !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
              write(VTK_file_unit,'(1F16.8)') err_DG
              !end do ERRH1_LOOP

        write(VTK_file_unit,'(A)')'h'
              !ELEM_LOOP: do i=1,nit
                !write(VTK_file_unit,'(I12,8I12)')nvert_per_el, t(i,:)-1
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

      write(VTK_file_unit,'(A)')'degree'
      write(VTK_file_unit,'(I1)') p

      write(VTK_file_unit,'(A)')'err_L2'
            !ERRL2_LOOP: do i=1,nit
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            write(VTK_file_unit,'(1F16.8)') err_L2
            !end do ERRL2_LOOP

      write(VTK_file_unit,'(A)')'err_DG'
            !ERRH1_LOOP: do i=1,nit
            !write(VTK_file_unit,'(3E16.8,1X)') xx(1:4,i),yy(1:4,i),zz(1:4,i)
            write(VTK_file_unit,'(1F16.8)') err_DG
            !end do ERRH1_LOOP

      write(VTK_file_unit,'(A)')'h'
            !ELEM_LOOP: do i=1,nit
              !write(VTK_file_unit,'(I12,8I12)')nvert_per_el, t(i,:)-1
            write(VTK_file_unit,'(1F16.8)') hh
            !end do ELEM_LOOP
      write(VTK_file_unit,'(A)') ' '

      close(unit=VTK_file_unit)

    end subroutine

    !> Store the numerical solution in an appropriate file
    !! HARDCODED 4 FOR 4 VERTICES OF TET
    !! MERGE SUBROUTINES...
    subroutine WRITE_SOLUTION(n_elem, PolyMesh, u, IsPoly, num_dt)

        use problem_data_and_properties
        use mpi
        use Poly_setup_MPI
        use mesh_partition_and_mpi_files

        implicit none

        integer(kind=4), intent(in), optional :: num_dt       !< number of timesteps
        type(Mesh_Structure), intent(inout) :: PolyMesh       !< mesh
        integer(kind=4), intent(in) :: n_elem                  !< local number of elements
        real(kind=8), dimension(DIM, NVERT_TET, n_elem), intent(inout) :: u !< 3D, 4 vertices, n_elem
        logical, intent(in) :: IsPoly
        real(kind=8), dimension(NVERT_TET, n_elem) :: xx, yy, zz
        character(len=80) :: vtk_filename_num!, vtk_filename_exact
        integer(kind=4) :: ie_loc,ivert,id_node!,i
        ! integer(kind=4) :: start_node(mpi_np), start_elem(mpi_np), gathered_sizes(mpi_np), start_solution(mpi_np)
        integer(kind=4) :: nvert_per_el

        nvert_per_el = PolyMesh%Elem_loc(1)%num_vert

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
        if (present(num_dt)) then
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

        else
            if(.not. IsPoly) vtk_filename_num = 'MONITORS/sol_tet_000000.vtk'

            write(vtk_filename_num(23:23),'(A)') 'V'
            call VTK_WRITE_SOLUTION(vtk_filename_num, xx,yy,zz, nvert_per_el, n_elem, PolyMesh%num_elem, 'solution', u, PolyMesh)
        endif

        !if (mpi_id==0) print *, 'Writing .vtk file...'
        !call VTK_WRITE_SOLUTION(vtk_filename_num, xx,yy,zz, nvert_per_el, n_elem, PolyMesh%num_elem, 'solution', u, PolyMesh)

        if (present(num_dt)) then
            !call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            !if (mpi_id==0) print *, 'Writing .vtu file...'

            !call VTU_WRITE_SOLUTION(vtk_filename_num, xx,yy,zz, nvert_per_el, n_elem, PolyMesh%num_elem, 'solution', u, PolyMesh, num_dt)

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
            ! start_node(1) = 0
            ! start_elem(1) = 0
            ! start_solution(1) = 0
            ! only in parallel
            ! if(mpi_np > 1) then
            !     call MPI_AllGather(PolyMesh%num_elem_loc, 1, MPI_INTEGER, gathered_sizes, 1, &
            !     MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
            !     start_elem = gathered_sizes
            !     ! mpi process from 0
            !     do i = mpi_np,2,-1
            !     start_elem(i) = start_elem(i-1)
            !     enddo
            !     start_elem(1) = 0
            !     do i = 2,mpi_np
            !         start_node(i) = start_node(i-1) + nvert_per_el*start_elem(i)
            !     enddo
            !     start_solution(1) = 0
            !     do i = 2,mpi_np
            !         start_solution(i) = start_solution(i-1) + DIM*nvert_per_el*start_elem(i)
            !     enddo
            !     do i = 2,mpi_np
            !         start_elem(i) = start_elem(i-1) + start_elem(i)
            !     enddo
            ! endif

            if (num_dt==0) then
                if (mpi_id==0) print *, 'Writing Ensight mesh.geo and poly.sca file...'
                call ENSIGHT_WRITE_MESH('MONITORS/', xx,yy,zz, nvert_per_el, n_elem, PolyMesh%num_elem, PolyMesh)

                if (mpi_id==0) print *, 'Writing Ensight boundaries.case file...'
                call ENSIGHT_WRITE_BOUNDARY('MONITORS/', NVERT_TRIA, n_elem, PolyMesh)
            endif

            if (mpi_id==0) print *, 'Writing Ensight .vec file...'
            call ENSIGHT_WRITE_SOLUTION('MONITORS/', nvert_per_el, n_elem, PolyMesh%num_elem, 'DISPLACEMENT', u, num_dt, PolyMesh=PolyMesh)

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

        endif

    end subroutine WRITE_SOLUTION

    !> Write .vtk files for the visualization of the partition and agglomeration of the mesh
    subroutine WRITE_MESH_VISUALIZATION_VTK(n_elem, PolyMesh, mpi_id)

        use mesh_partition_and_mpi_files

      implicit none

      type(Mesh_Structure), intent(inout) :: PolyMesh
      integer(kind=4), intent(in) :: n_elem
      integer(kind=4), intent(in) :: mpi_id
      real(kind=8), dimension(4,n_elem) :: xx, yy, zz
      character(len=80) :: vtk_filename_partition, vtk_filename_agglomeration
      integer(kind=4) :: ie_loc,ivert,id_node
        integer(kind=4) :: nvert_per_el

      nvert_per_el = PolyMesh%Elem_loc(1)%num_vert

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

      call VTK_WRITE_MESH_PARTITION(vtk_filename_partition, xx,yy,zz, nvert_per_el, n_elem, PolyMesh)
      call VTK_WRITE_MESH_AGGLOMERATION(vtk_filename_agglomeration, xx,yy,zz, nvert_per_el, n_elem, PolyMesh)

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

end module export_file_formats