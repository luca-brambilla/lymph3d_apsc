!   Author: Luca Brambilla
!   This file is part of the library LYMPH3D

!> module to exchage data with MPI. Communication between interfaces.
module exchange_data

contains

subroutine CHECK_MPI_EXCHANGE(PolyMesh, Np)

    use utilities
    use Poly_setup_MPI
    use global_parameters
    use Poly_mesh

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    integer(kind=4), intent(in) :: Np               !< number of element dofs per direction
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np) :: check_in
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np) :: check_out
    type(ScatteredArray), dimension(:,:), allocatable:: send_data !< variable size send buffers
    type(ScatteredArray), dimension(:,:), allocatable :: recv_data !< variable size receive buffers

    integer(kind=4) :: ie_loc,i

    if (mpi_id == 0) print *, 'CHECK MPI EXCHANGE'

    check_in=0.0d0
    do ie_loc=1,PolyMesh%num_elem_loc
        do i=1,DIM*Np
            check_in(ie_loc,i) = 1000000*mpi_id + (ie_loc-1)*DIM*Np + i
        enddo
    enddo

    call LYMPH3D_BARRIER

    if (mpi_id == 0) then
        print *, "prova_in row-by-row:"
        do ie_loc = 1, PolyMesh%num_elem_loc  ! Loop over rows
            print *, check_in(ie_loc, :)
        end do
    endif

    call MPI_EXCHANGE_ALLOCATE(PolyMesh, send_data, recv_data)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    call MPI_EXCHANGE_DOF(PolyMesh, check_in, check_out, send_data, recv_data)

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call MPI_EXCHANGE_DEALLOCATE(PolyMesh, send_data, recv_data)

    if (mpi_id==0) then
        print *, "check_out row-by-row:"
        do i = 1, PolyMesh%num_elem_inter_vec(mpi_id+1)  ! Loop over rows
            print *, mpi_id, check_out(i, :)
        end do
    endif

end subroutine CHECK_MPI_EXCHANGE

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
    real(kind=8), dimension(:,:), intent(in) :: input_sol   !< full input solution
    real(kind=8), dimension(:,:), intent(out) :: output_sol !< output solution containing interface data only

    integer(kind=4) :: id_send
    integer(kind=4) :: id_recv
    integer(kind=4) :: n_elem, n_dof
    integer(kind=4) :: Np, ie_glob, k, ireq, ie_send_loc, ie_recv_loc, el_sum
    integer(kind=4) :: row_send, row_recv

    type(ScatteredArray), dimension(mpi_np, mpi_np), intent(inout) :: send_data, recv_data
    integer(kind=4), allocatable :: requests(:), statuses(:,:)  ! For tracking operations

    ! real(kind=8) :: factor, ifactor

    ! factor = 1.0d0
    ! ifactor = 1.0d0

    allocate(requests(2 * mpi_np), statuses(2 * mpi_np, MPI_STATUS_SIZE))

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! send start request index
    ireq=1
    do id_recv=1,mpi_np
        do id_send=1,mpi_np

            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            ! check if data needs to be sent
            if (n_elem /= 0 .and. mpi_id == id_send-1 .and. id_recv/=id_send) then

                ! prepare send buffer
                do k=1,n_elem
                    ! get global index
                    ie_glob = PolyMesh%elem_inter_glo(PolyMesh%inter_disp(id_recv,id_send)+k)
                    ! find local index
                    ie_send_loc = PolyMesh%elem_glo2loc(ie_glob)

                    ! insert in buffer
                    row_send = (k-1)*Np*DIM
                    send_data(id_recv,id_send)%data(row_send+1:row_send+DIM*Np) = input_sol(ie_send_loc,:) ! * factor
                enddo

                ! multiply by number of dof per element and dimension 3D
                n_dof = n_elem * Np * DIM

                ! all-to-all asyncronous communication
                call MPI_ISEND(send_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_recv-1, 0, MPI_COMM_WORLD, requests(ireq), mpi_ierr)

                ! update request for send
                ireq = ireq + 2
            end if
        end do
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! receive start request index
    ireq = 0
    do id_recv=1,mpi_np
        do id_send=1,mpi_np
            n_elem = PolyMesh%num_elem_inter_comm(id_recv,id_send)
            ! check if data needs to be received
            if (n_elem /= 0 .and. mpi_id == id_recv-1 .and. id_recv/=id_send) then

                ! update request for receive
                ireq = ireq + 2

                ! multiply by number of dof per element and dimension 3D
                n_dof = n_elem * Np * DIM
                ! all-to-all asyncronous communication
                call MPI_IRECV(recv_data(id_recv,id_send)%data, n_dof, MPI_DOUBLE_PRECISION, id_send-1, 0, MPI_COMM_WORLD, requests(ireq), mpi_ierr)

                ! save receive buffer in local output solution
                el_sum = sum(PolyMesh%num_elem_inter_comm(id_recv,1:id_send-1))
                do k=1,n_elem
                    ie_recv_loc = el_sum + k
                    row_recv = (k-1)*Np*DIM
                    output_sol(ie_recv_loc,:) = recv_data(id_recv,id_send)%data(row_recv+1:row_recv+DIM*Np) !* ifactor
                enddo

            end if
        end do
    end do

    call MPI_WAITALL(ireq, requests, statuses, mpi_ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

end subroutine MPI_EXCHANGE_DOF

end module exchange_data