module utilities

#include<petsc/finclude/petscmat.h>
use Poly_setup_MPI

implicit none
contains

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

end module utilities