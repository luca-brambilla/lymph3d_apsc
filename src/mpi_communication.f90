!> Each process sends and receive a different amount of data from all other processes
subroutine MPI_SEND_ELEM(PolyMesh, send_dof, recv_dof)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    !integer(kind=4), dimension(mpi_np,mpi_np), intent(in) :: buffer_sizes !< sizes of the buffers to be sent
    integer(kind=4), dimension(:), intent(in) :: send_dof !< data to be sent
    integer(kind=4), dimension(:), intent(out) :: recv_dof !< data to be reveived

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: size
    integer(kind=4) :: Np

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    do id_send=1,mpi_np
        do id_recv=1,mpi_np
            size = PolyMesh%num_elem_inter_comm(id_send,id_recv)*Np
            if (id_send /= id_recv) then
                call MPI_SEND(send_dof, size, MPI_INTEGER, id_recv, 0, MPI_COMM_WORLD, mpi_ierr)
            end if
        end do
    end do

    do id_send=1,mpi_np
        do id_recv=1,mpi_np
            size = PolyMesh%num_elem_inter_comm(id_send,id_recv)*Np
            if (id_send /= id_recv) then
                call MPI_RECV(recv_dof, size, MPI_INTEGER, id_send, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
            end if
        end do
    end do



end subroutine MPI_SEND_ELEM

!> Each process sends and receive a different amount of data from all other processes
subroutine MPI_PROVA(PolyMesh, input_sol, output_sol)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    !integer(kind=4), dimension(mpi_np,mpi_np), intent(in) :: buffer_sizes !< sizes of the buffers to be sent
    integer(kind=4), dimension(:) :: input_sol !< old solution
    integer(kind=4), dimension(:) :: output_sol !< new solution

    integer(kind=4), dimension(:), allocatable :: send_dof !< data to be sent
    integer(kind=4), dimension(:), allocatable :: recv_dof !< data to be reveived

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem, n_dof
    integer(kind=4) :: Np, ie_loc, ie_glob, k
    integer(kind=4) :: row_loc, row_glo

    Np = PolyMesh%Elem_loc(1)%NDof_elem !! same for all elements per dimension

    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            ! multiply by number of dof per element and dimension 3D
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            if (n_elem /= 0 .and. mpi_id == id_send-1 ) then
                ! fill the send buffer with correct dofs
                do i=1,DIM
                    do k=1,n_elem
                        ie_glob = PolyMesh%elem_inter(PolyMesh%inter_disp(id_recv,id_send)+k)
                        ie_loc = PolyMesh%elem_glo2loc(ie_glob)
                        !! rows are wrong
                        row_glo = PolyMesh%num_elem*(i-1)*DIM + (ie_loc-1)*Np
                        row_sol = PolyMesh%num_elem_loc*(i-1)*DIM + (k-1)*Np
                        send_dof(row_loc+1:row_loc+Np) = input_sol(row_sol+1:row_sol+Np)
                    enddo
                enddo
                call MPI_SEND(send_dof, n_dof, MPI_INTEGER, id_recv-1, 0, MPI_COMM_WORLD, mpi_ierr)
            end if
        end do
    end do





end subroutine MPI_PROVA