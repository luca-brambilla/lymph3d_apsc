!   Author: Luca Brambilla
!   This file is part of the library LYMPH3D

!> module to exchage data with MPI. Communication between interfaces.
module exchange_data

contains

!> allocate send and receive buffer for interface communication
subroutine MPI_EXCHANGE_ALLOCATE(PolyMesh, send_data, recv_data)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    type(ScatteredArray), dimension(:,:), allocatable, intent(out) :: send_data !< variable size send buffers
    type(ScatteredArray), dimension(:,:), allocatable, intent(out) :: recv_data !< variable size receive buffers

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem
    integer(kind=4) :: n_dof
    integer(kind=4) :: Np

    allocate(send_data(mpi_np,mpi_np))
    allocate(recv_data(mpi_np,mpi_np))

    ! number of dof on one element per direction
    Np = PolyMesh%Elem_loc(1)%NDof_elem

    print *, 'allocate MPI interface buffers'

    ! allocate
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            if (n_elem /= 0 .and. mpi_id == id_send-1 .and. id_recv/=id_send) then
                allocate(send_data(id_recv,id_send)%data(n_dof))
                send_data(id_recv,id_send)%data = 0
                !print *, 'sending buffers -', id_send, id_recv, n_dof
            end if
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then
                allocate(recv_data(id_recv,id_send)%data(n_dof))
                recv_data(id_recv,id_send)%data = 0
                !print *, 'receive buffers -', id_send, id_recv, n_dof
            endif
        end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

end subroutine MPI_EXCHANGE_ALLOCATE

!> deallocate send and receive buffer for interface communication
subroutine MPI_EXCHANGE_DEALLOCATE(PolyMesh, send_data, recv_data)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    type(ScatteredArray), dimension(:, :), allocatable, intent(inout) :: send_data !< variable size send buffers
    type(ScatteredArray), dimension(:, :), allocatable, intent(inout) :: recv_data !< variable size receive buffers

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem

    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call FLUSH
    print *, 'deallocate MPI interface buffers'
    ! deallocate
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            if (n_elem /= 0 .and. mpi_id == id_send-1 .and. id_recv/=id_send) then
                deallocate(send_data(id_recv,id_send)%data)
            end if
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then
                deallocate(recv_data(id_recv,id_send)%data)
            endif
        end do
    end do

    deallocate(send_data, recv_data)

end subroutine MPI_EXCHANGE_DEALLOCATE

!> Each process sends and receive a different amount of data from all other processes
!> All-to-all implementation with asyncronous send and receive
subroutine MPI_EXCHANGE_DOF(PolyMesh, input_sol, output_sol, send_data, recv_data)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    real(kind=8), dimension(:), intent(in) :: input_sol   !< full input solution
    real(kind=8), dimension(:), intent(out) :: output_sol !< output solution containing interface data only


    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem, n_dof
    integer(kind=4) :: Np, ie_loc, ie_glob, k, i, ireq
    integer(kind=4) :: row_sol, row_send, row_recv

    type(ScatteredArray), dimension(mpi_np, mpi_np), intent(inout) :: send_data, recv_data
    integer, allocatable :: requests(:), statuses(:,:)  ! For tracking operations

    allocate(requests(2 * mpi_np), statuses(2 * mpi_np, MPI_STATUS_SIZE))

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! send
    ireq=1
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            ! multiply by number of dof per element and dimension 3D
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            ! check if data needs to be sent
            if (n_elem /= 0 .and. mpi_id == id_send-1 .and. id_recv/=id_send) then

                ! fill the send buffer with correct dofs
                do i=1,DIM
                    do k=1,n_elem
                        ! get global index
                        ie_glob = PolyMesh%elem_inter_glo(PolyMesh%inter_disp(id_recv,id_send)+k)
                        ! find local index
                        ie_loc = PolyMesh%elem_glo2loc(ie_glob)
                        ! ie_glob = PolyMesh%elem_inter_loc(PolyMesh%inter_disp(id_recv,id_send)+k)

                        ! insert in buffer
                        row_send = (i-1)*n_elem*Np + (k-1)*Np
                        row_sol = (i-1)*PolyMesh%num_elem_loc*Np + (ie_loc-1)*Np
                        send_data(id_recv,id_send)%data(row_send+1:row_send+Np) = input_sol(row_sol+1:row_sol+Np)
                    enddo
                enddo

                !call flush
                !print *, 'id', mpi_id, 'n_dof', n_dof, 'data send', send_data(id_recv,id_send)%data
                !print *, 'sending   - send: ', id_send-1, 'receive: ', id_recv-1, 'n_dof', n_dof

                ! all-to-all asyncronous communication
                call MPI_ISEND(send_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_recv-1, 0, MPI_COMM_WORLD, requests(ireq), mpi_ierr)

                ireq = ireq + 2 ! request for send and request for receive
            end if
        end do
    end do

    !call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    !print *, '$$$$$$$$$$$ END SENDING DATA $$$$$$$$$$$$'

    ! receive
    ireq = 2
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            !call FLUSH
            ! multiply by number of dof per element and dimension 3D
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then
                row_recv = sum(PolyMesh%num_elem_inter_comm(id_recv,1:id_send-1))*DIM*Np
                !print *, PolyMesh%num_elem_inter_comm(id_recv, 1:id_send)
                !print *, 'sum:', row_recv/(DIM*Np)

                ! all-to-all asyncronous communication
                call MPI_IRECV(recv_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_send-1, 0, MPI_COMM_WORLD, requests(ireq), mpi_ierr)

                ireq = ireq + 2 ! request for send and request for receive

                !call flush
                !print *, 'receiving - send: ', id_send-1, 'receive: ', id_recv-1, '- sum:', row_recv, ' - data:', PolyMesh%num_elem_inter_comm(id_recv, 1:id_send-1)
                ! print *, 'data receive', recv_data(id_recv,id_send)%data
                output_sol(row_recv+1:row_recv+n_dof) = recv_data(id_recv,id_send)%data
            end if
        end do
    end do

    call MPI_WAITALL(ireq-2, requests, statuses, mpi_ierr)
    !print *, '$$$$$$$$$$$$ END RECEIVING DATA $$$$$$$$$$$$', mpi_id

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

end subroutine MPI_EXCHANGE_DOF

end module exchange_data