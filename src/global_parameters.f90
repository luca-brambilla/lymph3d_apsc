!> Global constants
module global_parameters

    implicit none

    type :: ScatteredArray
        real(kind=8), dimension(:), allocatable :: data
    end type ScatteredArray

    integer(kind=4), parameter :: DIM = 3               !< dimension of the problem
    real(kind=8), parameter :: SQRT2 = sqrt(2.0d0)      !< parameter for @f$ \sqrt{2} @f$
    real(kind=8), parameter :: PI = 4.d0*datan(1.0d0)   !< parameter for @f$ \pi @f$

end module global_parameters