!> @brief Module with various subroutines and functions for
!> the computation of the determinant of a 3x3 matrix, the inverse of a 3x3
!> matrix and the  flipping of a vector.
!> Cholesky decomposition and direct solver for linear systems.
module vet_mat_operations

    implicit none

    contains

    ! Inverse of a 3×3 matrix A
    function matinv3(A) result(B)

        implicit none

        real(kind=8), dimension (3,3) :: A
        real(kind=8), dimension (3,3) :: B

        real(kind=8) :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1 / det3(A)

        ! Calculate the inverse of the matrix
        B(1,1) = + detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = - detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = + detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = - detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = + detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = - detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = + detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = - detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = + detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    end function matinv3

    !> Determinant of a 3x3 matrix A
    function det3(A) result(determinant)

        implicit none

        real(kind=8), dimension (3,3) :: A
        real(kind=8) :: determinant

        determinant = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2) &
                     -A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1) &
                     +A(1,3)*A(2,1)*A(3,2)-A(2,2)*A(1,3)*A(3,1)

    end function det3

    !> Determinant of a 4x4 matrix A
    function det4(A) result(determinant)
        implicit none

        real(kind=8), intent(in) :: A(4, 4)
        real(kind=8) :: determinant
        determinant = A(1,1) * (A(2,2) * (A(3,3) * A(4,4) - A(3,4) * A(4,3)) - &
                          A(2,3) * (A(3,2) * A(4,4) - A(3,4) * A(4,2)) + &
                          A(2,4) * (A(3,2) * A(4,3) - A(3,3) * A(4,2))) - &
              A(1,2) * (A(2,1) * (A(3,3) * A(4,4) - A(3,4) * A(4,3)) - &
                          A(2,3) * (A(3,1) * A(4,4) - A(3,4) * A(4,1)) + &
                          A(2,4) * (A(3,1) * A(4,3) - A(3,3) * A(4,1))) + &
              A(1,3) * (A(2,1) * (A(3,2) * A(4,4) - A(3,4) * A(4,2)) - &
                          A(2,2) * (A(3,1) * A(4,4) - A(3,4) * A(4,1)) + &
                          A(2,4) * (A(3,1) * A(4,2) - A(3,2) * A(4,1))) - &
              A(1,4) * (A(2,1) * (A(3,2) * A(4,3) - A(3,3) * A(4,2)) - &
                          A(2,2) * (A(3,1) * A(4,3) - A(3,3) * A(4,1)) + &
                          A(2,3) * (A(3,1) * A(4,2) - A(3,2) * A(4,1)))
    end function det4

    ! Write vector vec in the inverse way
    subroutine flip_vector(vec,n,flip_vec)

        real(kind=8), dimension(n) :: vec, flip_vec
        integer(kind=4) :: i, n

        do i = 1, n
          flip_vec(i) = vec(n-i+1)
        end do

    end subroutine flip_vector

    !> Cholesky decomposition for SPD matrix A=R^T R, with R upper triangular
    function cholesky(A, dim) result(R)

        integer(kind=4), intent(in) :: dim                  !< dimension of square matrix
        real(kind=8), intent(in) :: A(dim, dim)            !< symmetric positive definite matrix
        real(kind=8) :: R(dim, dim)                        !< upper triangular matrix
        integer(kind=4) :: i, j, k
        real(kind=8) :: sum
    
        ! Initialize R to zero
        R = 0.0
    
        ! Loop over rows
        do j = 1, dim
    
            ! Compute diagonal term
            sum = A(j, j)
            do k = 1, j - 1
                sum = sum - R(k, j)**2
            end do
            if (sum <= 0.0) then
                print *, "Matrix is not positive definite!"
                stop
            end if
            R(j, j) = dsqrt(sum)
    
            ! Compute off-diagonal terms
            do i = j + 1, dim
                sum = A(j, i)
                do k = 1, j - 1
                    sum = sum - R(k, j) * R(k, i)
                end do
                R(j, i) = sum / R(j, j)
            end do
    
        end do
    
    end function cholesky
    
    !> direct solver with LU factorization
    function solve_LU(L, U, b, dim) result(x)

        integer(kind=4), intent(in) :: dim
        real(kind=8), dimension(dim, dim), intent(in) :: L
        real(kind=8), dimension(dim, dim), intent(in) :: U
        real(kind=8), dimension(dim), intent(in) :: b
        real(kind=8), dimension(dim) :: x
        real(kind=8), dimension(dim) :: y
        integer(kind=4) :: i, j
    
        ! Forward substitution: Ly = b
        y(1) = b(1) / L(1, 1)
        do i = 2, dim
            y(i) = b(i)
            do j = 1, i-1
                y(i) = y(i) - L(i, j) * y(j)
            end do
            y(i) = y(i) / L(i, i)
        end do
    
        ! Backward substitution: Ux = y
        x(dim) = y(dim) / U(dim, dim)
        do i = dim-1, 1, -1
            x(i) = y(i)
            do j = i+1, dim
                x(i) = x(i) - U(i, j) * x(j)
            end do
            x(i) = x(i) / U(i, i)
        end do
    
    end function solve_LU
    

end module vet_mat_operations