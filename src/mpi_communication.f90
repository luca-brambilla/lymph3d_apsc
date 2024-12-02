!> module to exchage data with MPI. Communication between interfaces.
module exchange_data

contains

!> Each process sends and receive a different amount of data from all other processes
subroutine MPI_EXCHANGE_DOF(PolyMesh, input_sol, output_sol, comm_dof)

    use mpi
    use Poly_setup_mpi
    use Poly_mesh
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh !< mesh
    !integer(kind=4), dimension(mpi_np,mpi_np), intent(in) :: buffer_sizes !< sizes of the buffers to be sent
    integer(kind=4), intent(in) :: comm_dof
    real(kind=8), dimension(:), intent(in) :: input_sol
    real(kind=8), dimension(comm_dof), intent(out) :: output_sol
    ! real(kind=8), dimension(comm_dof) :: send_dof !< data to be sent
    ! real(kind=8), dimension(comm_dof) :: recv_dof !< data to be reveived

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem, n_dof
    integer(kind=4) :: Np, ie_loc, ie_glob, k, i
    integer(kind=4) :: row_sol, row_send, row_recv

    type(ScatteredArray), dimension(mpi_np, mpi_np) :: send_data, recv_data

    Np = PolyMesh%Elem_loc(1)%NDof_elem
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! allocate
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            if (n_elem /= 0 .and. mpi_id == id_send-1 .and. id_recv/=id_send) then
                allocate(send_data(id_recv,id_send)%data(n_dof))
                send_data(id_recv,id_send)%data = 0
                print *, 'send_data', id_send-1, id_recv-1
            end if
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then
                allocate(recv_data(id_recv,id_send)%data(n_dof))
                recv_data(id_recv,id_send)%data = 0
                print *, 'recv_data', id_send-1, id_recv-1
            endif
        end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! send
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            ! multiply by number of dof per element and dimension 3D
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
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
                        !print *, id_send-1, id_recv-1, '||', row_send+1, row_send+Np, row_sol+1, row_sol+Np
                        send_data(id_recv,id_send)%data(row_send+1:row_send+Np) = input_sol(row_sol+1:row_sol+Np)
                    enddo
                enddo
                call flush
                ! print *, 'data send', send_data(id_recv,id_send)%data
                print *, 'sending   - send: ', id_send-1, 'receive: ', id_recv-1
                !print *, '****************************'
                call MPI_SEND(send_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_recv-1, 0, MPI_COMM_WORLD, mpi_ierr)
            end if
        end do
    end do

    call FLUSH
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    !print *, '$$$$$$$$$$$ END SENDING DATA $$$$$$$$$$$$'

    ! receive
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            call FLUSH
            ! multiply by number of dof per element and dimension 3D
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            n_dof = n_elem * Np * DIM
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then
                row_recv = sum(PolyMesh%num_elem_inter_comm(id_recv,1:id_send-1))*DIM*Np
                !print *, PolyMesh%num_elem_inter_comm(id_recv, 1:id_send)
                !print *, 'sum:', row_recv/(DIM*Np)
                if (.not. allocated(recv_data(id_recv,id_send)%data)) then
                    print *, "Error: recv_data is not allocated on rank ", id_send-1, id_recv-1
                end if
                call MPI_RECV(recv_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_send-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_ierr)
                call flush
                print *, 'receiving - send: ', id_send-1, 'receive: ', id_recv-1, '- sum:', row_recv, ' - data:', PolyMesh%num_elem_inter_comm(id_recv, 1:id_send-1)
                ! print *, 'data receive', recv_data(id_recv,id_send)%data
                !print *, '**************************'
                output_sol(row_recv+1:row_recv+n_dof) = recv_data(id_recv,id_send)%data
            end if
        end do
    end do

    print *, '$$$$$$$$$$$$ END RECEIVING DATA $$$$$$$$$$$$', mpi_id


    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! deallocate
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            if (n_elem /= 0 .and. mpi_id == id_send-1 ) then
                deallocate(send_data(id_recv,id_send)%data)
            end if
            if (n_elem /= 0 .and. mpi_id == id_recv-1 ) then
                deallocate(recv_data(id_recv,id_send)%data)
            endif
        end do
    end do

end subroutine MPI_EXCHANGE_DOF

end module exchange_data