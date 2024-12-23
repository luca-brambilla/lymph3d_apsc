!> Utilities to save data to file, check FILES_MPI folder
module utilities

#include<petsc/finclude/petscmat.h>
use Poly_setup_MPI

implicit none
contains

!> Stop the simulation and print the duration
subroutine LYMPH3D_STOP

    use Poly_global

    implicit none

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    finish = MPI_WTIME()
    call calc_time(time_hour, time_min, time_sec, int(finish-start))

    if (mpi_id == 0) then
        write(*,'(A)')'--------------- SIMULATION WAS STOPPED ----------------'
        write(*,'(A)')
        if (IS_MatrixFree .eqv. .true.) then
            write(*,'(A)')'Matrix free'
        else
            write(*,'(A)')'PETSc full'
        endif
        call PRINT_PROFILING
        write(*,'(A)')
        write(*,'(A)')'-------------------------------------------------------'
        write(*,'(A,I2,A,I2,A,I2,A)') &
                'Simulation time = ', time_hour,' h ' , time_min,' m ' , time_sec,' s'
        write(*,'(A)')'-------------------------------------------------------'
        write(*,'(A)')
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    stop
    call PetscFinalize(mpi_ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
    call MPI_FINALIZE(mpi_ierr)
    stop
    call exit(0)

end subroutine LYMPH3D_STOP

!> save full PETSc matrix to a file, row by row
subroutine SAVE_MATRIX_PETSC(matrix, nrows, ncols, filename)
    use petscmat

    implicit none

    Mat :: matrix
    integer(kind=4), intent(in) :: nrows
    integer(kind=4), intent(in) :: ncols
    character(len=*), intent(in) :: filename

    integer(kind=4) :: i, unit_print
    PetscInt, dimension(:), allocatable :: cols
    PetscInt :: row(1)
    PetscScalar, dimension(:), allocatable :: values

    allocate(cols(ncols))
    allocate(values(ncols))
    cols = [(i,i=0,ncols-1)]

    open(newunit=unit_print, action='WRITE', file=filename, &
    form='FORMATTED', status='replace')
    do i=0,nrows-1
        row(1) = i
        call MatGetValues(matrix, 1, row, ncols, cols, values, ierr)
        write(unit_print, *) values
    enddo
    close(unit=unit_print)

    deallocate(cols, values)

end subroutine SAVE_MATRIX_PETSC

!> save Fortan matrices of an element to a file, ordered first by neighbour then by row
subroutine SAVE_MATRIX_F90(matrix, nrows, filename)

    implicit none

    real(kind=8), dimension(:,:,:), intent(in) :: matrix
    integer(kind=4), intent(in) :: nrows
    character(len=*), intent(in) :: filename

    integer(kind=4) :: i, unit_print, ie

    open(newunit=unit_print, action='WRITE', file=filename, &
    form='FORMATTED', status='replace')
    do ie=1,5
        write(unit_print, *) ie
        do i=1,nrows
            write(unit_print, *) matrix(ie,i,:)
        enddo
    enddo
    close(unit=unit_print)

end subroutine SAVE_MATRIX_F90

!> save solution vector, element by element
subroutine SAVE_VECTOR(vector, nrows, filename)
    implicit none

    integer(kind=4), intent(in) :: nrows
    real(kind=8), dimension(:,:), intent(in) :: vector
    character(len=*), intent(in) :: filename

    integer(kind=4) :: i, unit_print

    open(newunit=unit_print, action='WRITE', file=filename, &
    form='FORMATTED', status='replace')
    do i=1,nrows
        write(unit_print, *) vector(i,:)
    enddo
    close(unit=unit_print)

end subroutine SAVE_VECTOR

!> Check files in the folder FILES_MPI and delete if the program was previously
!> ran with a different number of processes
subroutine CHECK_MPI_FILES

    implicit none

    logical :: file_exists
    character(len=100) :: filename

    ! checking mesh files
    ! A required file is missing; delete all files
    write(filename, '("FILES_MPI/mesh_", i6.6, ".mpi")') mpi_np-1
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
        write(*,'(A,I0)') "Program previously ran with fewer processes than: ", mpi_np
        call DELETE_ALL_FILES
        return
    end if

    ! Check for any extra files
    write(filename, '("FILES_MPI/mesh_", i6.6, ".mpi")') mpi_np
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
        write(*,'(A,I0)') "Program previously ran with more processes than: ", mpi_np
        call DELETE_ALL_FILES
    end if

end subroutine CHECK_MPI_FILES

!> Delete all files in the FILES_MPI directory
subroutine DELETE_ALL_FILES

    implicit none

    character(len=200) :: command

    ! System command to delete all files
    command = "rm -f FILES_MPI/*"
    call EXECUTE_COMMAND_LINE(command, wait=.true., exitstat=ierr)

    if (ierr /= 0) then
        print *, "Error deleting files. Exit code:", ierr
    else
        print *, "All files deleted in FILES_MPI directory."
    end if
end subroutine DELETE_ALL_FILES

!> Alias for MPI barrier
subroutine LYMPH3D_BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
end subroutine

subroutine PRINT_PROFILING
    use Poly_setup_MPI
    use Poly_global
    implicit none

    call LYMPH3D_BARRIER

    if (mpi_id == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, tp_copy_matrix, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_linear_system, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_KU, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_copy_vector, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_system_setup, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_setup_K, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_setup_M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_setup_RHS, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_partition, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_export, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_exact, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(MPI_IN_PLACE, tp_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
    else
        call MPI_REDUCE(tp_copy_matrix, tp_copy_matrix, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_linear_system, tp_linear_system, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_KU, tp_KU, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_copy_vector, tp_copy_vector, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_system_setup, tp_system_setup, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_setup_K, tp_setup_K, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_setup_M, tp_setup_M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_setup_RHS, tp_setup_RHS, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_partition, tp_partition, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_export, tp_export, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_exact, tp_exact, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_REDUCE(tp_error, tp_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, mpi_ierr)
    endif

    if (mpi_id == 0) then
        print *,''
        print *,'------------- PROFILING -------------'
        print *,''
        ! partition
        call calc_time(time_hour, time_min, time_sec, int(tp_partition))
        print *, 'tp_partition     = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! assemble stiffness
        call calc_time(time_hour, time_min, time_sec, int(tp_setup_K))
        print *, 'tp_setup_K       = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! assemble mass
        call calc_time(time_hour, time_min, time_sec, int(tp_setup_M))
        print *, 'tp_setup_M       = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! assemble rhs
        call calc_time(time_hour, time_min, time_sec, int(tp_setup_RHS))
        print *, 'tp_setup_RHS     = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! product K*U
        call calc_time(time_hour, time_min, time_sec, int(tp_KU))
        print *, 'tp_KU            = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! solve mass linear system
        call calc_time(time_hour, time_min, time_sec, int(tp_linear_system))
        print *, 'tp_linear_system = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! assign pointer to temp matrix
        call calc_time(time_hour, time_min, time_sec, int(tp_copy_matrix))
        print *, 'tp_copy_matrix   = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! assign pointer to temp vector
        call calc_time(time_hour, time_min, time_sec, int(tp_copy_vector))
        print *, 'tp_copy_vector   = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! setup Krylov solver
        call calc_time(time_hour, time_min, time_sec, int(tp_system_setup))
        print *, 'tp_system_setup  = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! export
        call calc_time(time_hour, time_min, time_sec, int(tp_export))
        print *, 'tp_export        = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! exact solution
        call calc_time(time_hour, time_min, time_sec, int(tp_exact))
        print *, 'tp_exact         = ', time_hour,' h ', time_min,' m ', time_sec,' s'
        ! compute error
        call calc_time(time_hour, time_min, time_sec, int(tp_exact))
        print *, 'tp_exact         = ', time_hour,' h ', time_min,' m ', time_sec,' s'
    endif
    
end subroutine

!> print all neighbor bounding boxes for each polyhedron on the process
subroutine PRINT_NEIGHBOR_BBOX(PolyMesh)
    use Poly_mesh
    use local_search

    implicit none

    type(Mesh_Structure), intent(inout) :: PolyMesh
    integer(kind=4) :: E1, E2, ipoly_glob, ipoly_loc, ie_glob, iface, ipoly2_glob, ipoly2_loc, proc_id

    do E1=1,PolyMesh%num_elem_loc
        ie_glob = PolyMesh%elem_loc2glo(E1)
        print *, '------ ie_loc', E1, ' - ie_glob:', ie_glob, '------'
        face_loop: do iface=1,4
            ! neighbor
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)
            proc_id = PolyMesh%Elem_loc(E1)%neigh_el(iface,0)
            if (E2<0) then
                write(*,'(A,I3,A,I3)') 'neighbor:', E2
                cycle face_loop
            endif

            ipoly2_glob = PolyMesh%elem_in_poly(E2)
            if (proc_id==mpi_id) then
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                PolyMesh%num_poly_loc, ipoly2_glob, ipoly2_loc)

                write(*, '(A,I3,A,I3,A,6F8.3)') 'neighbor:', ipoly2_glob, ' process ', proc_id, ' neigh_bbox: ', PolyMesh%Poly(ipoly2_loc)%b_box
            else

                ipoly_glob = PolyMesh%elem_in_poly(ie_glob)
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                PolyMesh%num_poly_loc, ipoly_glob, ipoly_loc)

                write(*, '(A,I3,A,I3,A,6F8.3)') 'neighbor:', ipoly2_glob, ' process ', proc_id, ' neigh_bbox: ', PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface,:,:)

            endif
        enddo face_loop
    enddo

end subroutine PRINT_NEIGHBOR_BBOX

end module utilities