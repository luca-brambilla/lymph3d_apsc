module utilities

#include<petsc/finclude/petscmat.h>
use Poly_setup_MPI

implicit none
contains

!> save full PETSc matrix to a file, row by row
subroutine save_matrix(matrix, nrows, ncols, filename)
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

end subroutine

!> save solution vector, element by element
subroutine save_vector(vector, nrows, filename)
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

end subroutine

!> Check files in the folder FILES_MPI and delete if the program was previously
!> ran with a different number of processes
subroutine check_mpi_files

    implicit none

    logical :: file_exists
    character(len=100) :: filename

    ! checking mesh files
    ! A required file is missing; delete all files
    write(filename, '("FILES_MPI/mesh_", i6.6, ".mpi")') mpi_np-1
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
        print *, "Program previously ran with fewer processes than", mpi_np
        call delete_all_files()
        return
    end if

    ! Check for any extra files
    write(filename, '("FILES_MPI/mesh_", i6.6, ".mpi")') mpi_np
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
        print *, "Program previously ran with more processes than", mpi_np
        call delete_all_files()
    end if

end subroutine check_mpi_files

!> Delete all files in the FILES_MPI directory
subroutine delete_all_files()

    implicit none

    character(len=200) :: command

    ! System command to delete all files
    command = "rm -f FILES_MPI/*"
    call execute_command_line(command, wait=.true., exitstat=ierr)

    if (ierr /= 0) then
        print *, "Error deleting files. Exit code:", ierr
    else
        print *, "All files deleted in FILES_MPI directory."
    end if
end subroutine delete_all_files


end module utilities