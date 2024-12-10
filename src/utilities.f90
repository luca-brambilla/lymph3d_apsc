module utilities

#include<petsc/finclude/petscmat.h>

implicit none
contains

!> save full PETSc matrix to a file, row by row
subroutine save_matrix(matrix, nrows, ncols, filename)
    use petscmat
    use Poly_setup_MPI

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

subroutine save_vector
end subroutine

end module utilities