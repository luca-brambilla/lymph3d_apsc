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
    
end module vet_mat_operations