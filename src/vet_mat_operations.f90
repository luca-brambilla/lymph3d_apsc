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
        
    ! Determinant of a 3x3 matrix A
    function det3(A) result(determinant)
    
        implicit none
    
        real(kind=8), dimension (3,3) :: A
        real(kind=8) :: determinant
        
        determinant = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2) &
                     -A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1) &
                     +A(1,3)*A(2,1)*A(3,2)-A(2,2)*A(1,3)*A(3,1)
        
    end function det3

    ! Write vector vec in the inverse way
    subroutine flip_vector(vec,n,flip_vec)

        real(kind=8), dimension(n) :: vec, flip_vec
        integer(kind=4) :: i, n
      
        do i = 1, n
          flip_vec(i) = vec(n-i+1)
        end do
      
    end subroutine flip_vector
    
    ! Cholesky decomposition for SPD matrix A=R^T R, with R upper triangular

    function cholesky(mat,Np) result(res)

        integer(kind=4) :: Np
        real(kind=8), dimension(Np,Np) :: mat
        real(kind=8), dimension(Np,Np) :: res
        integer(kind=4) :: i,j,k
        res = 0.0

        res(1,1) = sqrt(mat(1,1))
        do j=2,Np

            ! extra-diagonal terms
            do i=1,j-1
                res(i,j) = mat(i,j)
                do k=1,i-1
                    res(i,j) = res(i,j) - res(k,i)*res(k,j)
                end do
                res(i,j) = res(i,j) / sqrt(res(i,i))
            end do

            ! diagonal terms
            res(1,1) = mat(j,j)
            do k=1,j-1
                res(j,j)=res(j,j) - res(k,j)*res(k,j)
            end do
            res(j,j) = sqrt(res(j,j))

        end do

    end function cholesky

    function linear_system(Np, A, b, type) result(x)

        integer(kind=4) :: Np
        integer(kind=4) :: type
        real(kind=8), dimension(Np,Np) :: A
        real(kind=8), dimension(Np,Np) :: R
        real(kind=8), dimension(Np) :: b
        real(kind=8), dimension(Np) :: x
        real(kind=8), dimension(Np) :: y
        integer(kind=4) :: i,j

        ! Cholesky
        if(type == 0) then
            R = cholesky(A,Np)

            ! Ly = b, forward substitution
            ! L = R^T
            y(1) = b(1)/R(1,1)
            do i=2,Np
                y(i) = b(i)
                do j=1,Np-1
                    y(i) = y(i) - R(j,i) * b(j)
                end do
                y(i) = y(i) / R(i,i)
            end do

            ! Ux = y, backward substitution
            ! U = R
            x(Np) = y(Np)/R(Np,Np)
            do i=Np-1,1
                x(i) = y(i)
                do j=2,Np
                    x(i) = x(i) - R(i,j) * y(j)
                end do
                x(i) = x(i) / R(i,i)
            end do
        end if
    end function linear_system


end module vet_mat_operations