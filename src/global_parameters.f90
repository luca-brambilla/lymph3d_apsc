!> Global constants
module global_parameters

    implicit none

    type :: ScatteredArray
        real(kind=8), dimension(:), allocatable :: data
    end type ScatteredArray

    integer(kind=4), parameter :: DIM = 3               !< dimension of the problem
    real(kind=8), parameter :: SQRT2 = sqrt(2.0d0)      !< parameter for @f$ \sqrt{2} @f$
    real(kind=8), parameter :: PI = 4.d0*datan(1.0d0)   !< parameter for @f$ \pi @f$
    real(kind=8), parameter :: TOL = 1.0d-40            !< tolerance for small numbers

contains

    subroutine STOP_LYMPH3D
        use mpi
        use Poly_setup_MPI

        implicit none

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        call PetscFinalize(mpi_ierr)
        call MPI_FINALIZE(mpi_ierr)
        stop
    end subroutine STOP_LYMPH3D

end module global_parameters