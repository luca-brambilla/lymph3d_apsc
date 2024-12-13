!> Mappings from and to reference elements
module Poly_ref_mappings

    use vet_mat_operations
    use global_parameters

    implicit none

    contains

!! NO NEED TO COMPUTE AT EACH CALL ?????
!! Jinv NEVER USED???
!! Jdet ONLY USED FOR MODAL COEFFICIENTS
!> Compute the map from reference to physical tetrahedron, given tetrahedral element coordinates x, y, z:
!> \f[ \mathbf{x} = F_K(\widehat{\mathbf{x}}) = J_K \widehat{\mathbf{x}} + \mathbf{c} \f]
!> then compute the inverse and the determinant of the Jacobian \f(J_K\f).
!> Coordinates of the reference tetrahedron \f( (0,0,0),(1,0,0),(0,1,0),(0,0,1) \f).
subroutine jacobians(x, y, z, Fk, Jinv, Jdet)

    real(kind=8), dimension(NVERT_TET), intent(in) :: x     !< physical tetrahedron x coordinates
    real(kind=8), dimension(NVERT_TET), intent(in) :: y     !< physical tetrahedron y coordinates
    real(kind=8), dimension(NVERT_TET), intent(in) :: z     !< physical tetrahedron z coordinates
    real(kind=8), dimension(DIM,DIM+1), intent(out) :: Fk   !< coefficients for transformation from reference to physical tetrahedron.
    !< `Fk(:,1:DIM)` is the Jacobian \f(J_K\f) and `Fk(:,DIM+1)` is \f( \mathbf{c} \f)
    real(kind=8), dimension(DIM,DIM), intent(out) :: Jinv !< inverse of the Jacobian
    real(kind=8), intent(out) :: Jdet   !< determinant of the Jacobian

    real(kind=8), dimension (DIM,DIM) :: jacobian
    real(kind=8), dimension(NVERT_TET,NVERT_TET) :: mat
    integer(kind=4) :: i,j,l

    ! sign to add/subtract vertex coordinates in the transformation
    do i=1,4
        do j=1,4

            mat(i,j) = 0

            if (j==1) then
                if (i==4) then
                    mat(i,j) = 1
                else
                    mat(i,j) = -1
                end if
            end if

            if (j>1 .and. j==i+1) then
                mat(i,j) = 1
            end if

        end do
    end do

    ! compute coefficients for the tranformation
    do j=1,4
        Fk(1,j)=0.0
        Fk(2,j)=0.0
        Fk(3,j)=0.0
        do l=1,4
            Fk(1,j) = Fk(1,j) + mat(j,l)*x(l)
            Fk(2,j) = Fk(2,j) + mat(j,l)*y(l)
            Fk(3,j) = Fk(3,j) + mat(j,l)*z(l)
        end do
    end do

    ! compute inverse and determinant of the Jacobian
    do j=1,DIM
        do l=1,DIM
            jacobian(j,l) = Fk(j,l)
        end do
    end do

    Jdet = det3(jacobian)       ! see vet_mat_operations.f90
    Jinv = matinv3(jacobian)    ! see vet_mat_operations.f90

end subroutine jacobians

!! INVERSE MAP NEVER USED??? IS IT EVEN CORRECT? DIMENSIONS?
!> Maps from the 2D reference triangle to the faces of the 3D reference tetrahedron and their inverses
subroutine tria2tetfaces_maps(node_maps, node_maps_inv)

    real(kind=8), dimension(NVERT_TET,NVERT_TET,NVERT_TET), intent(out) :: node_maps    !< mapping from reference triangle to tetrahedron face
    real(kind=8), dimension(2,3,4), intent(out) :: node_maps_inv !< mapping from tetrahedron face to reference triangle

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

!> map from the reference cube to the reference tetrahedron
subroutine mapping_quadrature_3D(nod3, wei3, nq3, nodtet3, weitet3)

    real(kind=8), dimension(4, nq3) :: nod3 !< 3D quadrature nodes on reference cube
    real(kind=8), dimension(nq3) :: wei3    !< 3D quadrature weights on reference cube
    integer(kind=4), intent(in) :: nq3  !< number of 3D quadrature nodes on reference cube
    real(kind=8), dimension(4, nq3), intent(out) :: nodtet3 !< 3D quadrature nodes on reference tetrahedron
    real(kind=8), dimension(nq3), intent(out) :: weitet3    !< 3D quadrature weights on reference tetrahedron

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

!> map from the reference square to the reference triangle
subroutine mapping_quadrature_2D(nod2, wei2, nq2, nodtria2, weitria2)

    real(kind=8), dimension(4, nq2), intent(in) :: nod2 !< 2D quadrature nodes on reference square
    real(kind=8), dimension(nq2), intent(in) :: wei2    !< 2D quadrature weights on reference square
    integer(kind=4), intent(in) :: nq2                  !< number of 2D quadrature nodes on reference square
    real(kind=8), dimension(4, nq2), intent(out) :: nodtria2    !< 2D quadrature nodes on reference triangle
    real(kind=8), dimension(nq2), intent(out) :: weitria2       !< 2D quadrature weights on reference triangle

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