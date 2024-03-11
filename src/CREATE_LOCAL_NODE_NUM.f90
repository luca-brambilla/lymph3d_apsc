! Store local numeration to reconstruct the numerical solution
subroutine CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)
    
    use mpi
    use Poly_setup_MPI
    use Poly_mesh
   
    implicit none

    integer(kind=4) :: local_dof, local_dof_send, last_elem_send
    integer(kind=4) :: ip, i
    integer(kind=4), dimension(local_dof) :: nnod_num

    local_dof_send = 0;

    do ip=1,mpi_np

          if (mpi_id == ip-1) then

            do i=1,local_dof
                  nnod_num(i) = local_dof_send+i;
            enddo

            last_elem_send = nnod_num(i-1);

          endif

          if (mpi_id == ip-1) local_dof_send = last_elem_send;

          call MPI_BCAST(local_dof_send, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

          call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
          
    enddo

end subroutine CREATE_LOCAL_NODE_NUM