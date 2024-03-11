module Poly_ref_mappings
   
    use vet_mat_operations

    implicit none

    contains
    
! Function that computes for a given tetrahedral element of coordinates x, y and z:
! - the map Fk from the reference tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1) to the physical tetrahedron; 
!       note: Fk(x_hat) = Jacobian*x_hat + c,      where Jacobian=Fk(:,1:3), c=Fk(:,4)
! - the determinant of the jacobian Jdet;
! - the inverse of the jacobian Jinv.
subroutine jacobians(x, y, z, Fk, Jinv, Jdet)

    real(kind=8), dimension(4), intent(in) :: x, y, z
    real(kind=8), dimension(3,4), intent(out) :: Fk
    real(kind=8), dimension(3,3), intent(out) :: Jinv
    real(kind=8), intent(out) :: Jdet

    real(kind=8), dimension (3,3) :: Jacobian
    real(kind=8), dimension(4,4) :: Mat
    integer(kind=4) :: i,j,l
    
    do i=1,4
        do j=1,4 

            Mat(i,j) = 0

            if (j==1) then
                if (i==4) then
                    Mat(i,j) = 1
                else
                    Mat(i,j) = -1
                end if
            end if

            if (j>1 .and. j==i+1) then
                Mat(i,j) = 1
            end if

        end do
    end do

    do j=1,4
        Fk(1,j)=0.0
        Fk(2,j)=0.0
        Fk(3,j)=0.0
        do l=1,4
            Fk(1,j) = Fk(1,j) + Mat(j,l)*x(l)
            Fk(2,j) = Fk(2,j) + Mat(j,l)*y(l)
            Fk(3,j) = Fk(3,j) + Mat(j,l)*z(l)
        end do
    end do

    do j=1,3
        do l=1,3
            Jacobian(j,l) = Fk(j,l)
        end do
    end do
    
    Jdet = det3(Jacobian) ! see vet_mat_operations.f90
    Jinv = matinv3(Jacobian) ! see vet_mat_operations.f90
    
end subroutine jacobians

! Maps from the 2D reference triangle to the faces of the 3D reference tetrahedron and their inverses
subroutine tria2tetfaces_maps(node_maps, node_maps_inv)

    real(kind=8), dimension(4,4,4), intent(out) :: node_maps
    real(kind=8), dimension(2,3,4), intent(out) :: node_maps_inv

    integer(kind=4) :: i, j, k

    do k=1,4 ! index of the face of the reference tetrahedron
        do i=1,4
            do j=1,4
                node_maps(i,j,k) = 0.0
            end do
        end do
        do i=1,2
            do j=1,3
                node_maps_inv(i,j,k) = 0.0
            end do
        end do
    end do
    
    ! k==1 
    node_maps(1,2,1) = 1.0
    node_maps(2,1,1) = 1.0
    node_maps(4,4,1) = 1.0

    node_maps_inv(1,2,1) = 1.0
    node_maps_inv(2,1,1) = 1.0

    ! k==2
    node_maps(1,1,2) = 1.0
    node_maps(3,2,2) = 1.0
    node_maps(4,4,2) = 1.0
    
    node_maps_inv(1,1,2) = 1.0
    node_maps_inv(2,3,2) = 1.0

    ! k==3
    node_maps(1,1,3) = 1.0
    node_maps(2,2,3) = 1.0
    node_maps(3,3,3) = 1.0
    node_maps(4,4,3) = 1.0
    node_maps(3,1,3) = -1.0
    node_maps(3,2,3) = -1.0
    
    node_maps_inv(1,1,3) = 1.0
    node_maps_inv(2,2,3) = 1.0

    ! k==4
    node_maps(2,2,4) = 1.0
    node_maps(3,1,4) = 1.0
    node_maps(4,4,4) = 1.0

    node_maps_inv(1,3,4) = 1.0
    node_maps_inv(2,2,4) = 1.0

end subroutine tria2tetfaces_maps

! map from the reference cube to the reference tetrahedron
subroutine mapping_quadrature_3D(nod3, wei3, nq3, nodtet3, weitet3)

    real(kind=8), dimension(4, nq3) :: nod3
    real(kind=8), dimension(nq3) :: wei3
    integer(kind=4), intent(in) :: nq3
    real(kind=8), dimension(4, nq3), intent(out) :: nodtet3
    real(kind=8), dimension(nq3), intent(out) :: weitet3

    real(kind=8), dimension(nq3) :: Jdet_tet

    ! first compute GL nodes and weights on the cube (0,1)^3
    nod3(1:3,:) = 0.5 * nod3(1:3,:) + 0.5
    wei3 = 0.125 * wei3

    nodtet3 = 0.0

    nodtet3(1,:) = nod3(1,:)
    nodtet3(2,:) = (1 - nod3(1,:)) * nod3(2,:)
    nodtet3(3,:) = (1 - nod3(1,:)) * (1 - nod3(2,:)) * nod3(3,:)
    nodtet3(4,:) = nod3(4,:)

    Jdet_tet = (1 - nod3(2,:)) * (1 - nod3(1,:))**2

    weitet3 = abs(Jdet_tet) * wei3

end subroutine

! map from the reference square to the reference triangle
subroutine mapping_quadrature_2D(nod2, wei2, nq2, nodtria2, weitria2)

    real(kind=8), dimension(4, nq2), intent(in) :: nod2
    real(kind=8), dimension(nq2), intent(in) :: wei2
    integer(kind=4), intent(in) :: nq2
    real(kind=8), dimension(4, nq2), intent(out) :: nodtria2
    real(kind=8), dimension(nq2), intent(out) :: weitria2

    real(kind=8), dimension(nq2) :: Jdet_tria

    nodtria2 = 0.0

    nodtria2(1,:) = 0.5 * (1 + nod2(1,:))
    nodtria2(2,:) = 0.25 * (1 - nod2(1,:)) * (1 + nod2(2,:))
    nodtria2(3,:) = nod2(3,:)
    nodtria2(4,:) = nod2(4,:)

    Jdet_tria = 0.125 * (1 - nod2(1,:))

    weitria2 = abs(Jdet_tria) * wei2

end subroutine

end module Poly_ref_mappings