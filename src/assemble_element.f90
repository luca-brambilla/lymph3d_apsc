!> Assemble contributions for element matrices and element vectors
module assemble_element

    use problem_data_and_properties
    use Poly_mesh
    use Poly_ref_mappings
    use basis_function
    use global_parameters

    implicit none

    contains

!! OLD INDICES i,j,m,n
!! TODO change indices
!! TODO check hardcoded
!! change from LOCAL to ELEMENT
!! REMOVE SUMS FROM DOCUMENTATION???

!> Assemble the element stiffness matrix `stiff_tet_vol` approximating the volume integral over the tetrahedral element
!> \f[ [V_{K}]_{ij} = \int_K  \boldsymbol{\sigma}(\boldsymbol{\varphi}_{j,K}) : \boldsymbol{\varepsilon}(\boldsymbol{\varphi}_{i,K}) \f]
subroutine MAKE_STIFFNESS_VOLUME(Np, Jdet, weitet3, nq3, lambda, mu, dphi, stiff_tet_vol)

    ! dphi is provided by the subroutine basis in basis_function.f90
    ! weitet3 is provided by the subroutine mapping_quadrature_3D in Poly_ref_mappings.f90
    ! nq3 is provided by the subroutine quadrature in basis_function.f90
    ! Jdet is provided by the subroutine jacobians in Poly_ref_mappings.f90

    integer(kind=4), intent(in) :: nq3  !< number of 3D quadrature nodes
    integer(kind=4), intent(in) :: Np   !< number of element dof per dimension
    real(kind=8), intent(in) :: Jdet    !< transormation determinant
    real(kind=8), intent(in) :: lambda  !< Lamé 1st parameter
    real(kind=8), intent(in) :: mu      !< Lamé 2nd parameter
    real(kind=8), dimension(nq3), intent(in) :: weitet3 !< tetrahedron 3D quadrature weights
    real(kind=8), dimension(3,Np,nq3), intent(in) :: dphi   !< basis function derivative evaluations
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: stiff_tet_vol    !< element stiffness matrix volumetric contribution V

    integer(kind=4) :: q, i, j, m, n
    real(kind=8), dimension(DIM,DIM,Np,Np) :: temp

    stiff_tet_vol = 0.0

    ! loop on 3D quadrature points
    do q = 1,nq3

        temp = 0.0

        do m=1,Np
            do n=1,Np

                temp(1,1,m,n) = (lambda + 2*mu)*dphi(1,m,q)*dphi(1,n,q) + &
                                mu*dphi(2,m,q)*dphi(2,n,q) + mu*dphi(3,m,q)*dphi(3,n,q)
                temp(2,2,m,n) = (lambda + 2*mu)*dphi(2,m,q)*dphi(2,n,q) + &
                                mu*dphi(1,m,q)*dphi(1,n,q) + mu*dphi(3,m,q)*dphi(3,n,q)
                temp(3,3,m,n) = (lambda + 2*mu)*dphi(3,m,q)*dphi(3,n,q) + &
                                mu*dphi(1,m,q)*dphi(1,n,q) + mu*dphi(2,m,q)*dphi(2,n,q)
                do i=1,2
                    do j=i+1,3
                        temp(i,j,m,n) = lambda*dphi(i,m,q)*dphi(j,n,q) + mu*dphi(j,m,q)*dphi(i,n,q)
                        temp(j,i,m,n) = lambda*dphi(j,m,q)*dphi(i,n,q) + mu*dphi(i,m,q)*dphi(j,n,q)
                    enddo
                enddo

                do i=1,DIM
                    do j=1,DIM
                        stiff_tet_vol(i,j,m,n) = stiff_tet_vol(i,j,m,n) + weitet3(q)*abs(Jdet)*temp(i,j,m,n)
                    enddo
                enddo

            end do
        end do

    enddo

end subroutine MAKE_STIFFNESS_VOLUME

!> Assemble the mass matrix `mass_tet_vol` approximating the volume integral over the tetrahedral element
!> \f[ [M_{K}]_{ij} = \int_K \boldsymbol{\phi}_{j,K} \cdot \boldsymbol{\phi}_{i,K} \f]
subroutine MAKE_MASS_VOLUME(Np, Jdet, weitet3, nq3, phi, rho, mass_tet_vol)

    ! phi is provided by the subroutine basis in basis_function.f90
    ! weitet3 is provided by the subroutine mapping_quadrature_3D in Poly_ref_mappings.f90
    ! nq3 is provided by the subroutine quadrature in basis_function.f90
    ! Jdet is provided by the subroutine jacobians in Poly_ref_mappings.f90

    integer(kind=4), intent(in) :: nq3      !< number of 3D quadrature nodes
    integer(kind=4), intent(in) :: Np       !< number of element dofs per dimension
    real(kind=8), intent(in) :: Jdet        !< transformation determinant
    real(kind=8), dimension(nq3), intent(in) :: weitet3 !< tetrahedron 3D quadrature weights
    real(kind=8), dimension(Np,nq3), intent(in) :: phi  !< element basis function
    real(kind=8), intent(in) :: rho         !< element density
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: mass_tet_vol !< element mass matrix M

    integer(kind=4) :: q, i, m, n

    mass_tet_vol = 0.0

    ! loop on 3D quadrature points
    do q = 1,nq3
        do m=1,Np
            do n=1,Np

                do i = 1,DIM
                    mass_tet_vol(i,i,m,n) = mass_tet_vol(i,i,m,n) + rho * weitet3(q)*abs(Jdet)*phi(m,q)*phi(n,q)
                enddo

            end do
        end do
    enddo

end subroutine MAKE_MASS_VOLUME

!> Assemble the term `rhs_tet_vol` approximating the volume integral over the tetrahedral element
!> \f[ [F_{K}]_{i} = \int_K \boldsymbol{f} \cdot \boldsymbol{\varphi}_{i,K} \f]
subroutine MAKE_RHS_VOLUME(Np, Fk, Jdet, nodtet3, weitet3, nq3, lambda, mu, phi, rhs_tet_vol, rho)

    ! phi is provided by the subroutine basis in basis_function.f90
    ! weitet3 is provided by the subroutine mapping_quadrature_3D in Poly_ref_mappings.f90
    ! nq3 is provided by the subroutine quadrature in basis_function.f90
    ! Fk and Jdet are provided by the subroutine jacobians in Poly_ref_mappings.f90

    ! density is 0.0 in static case
    real(kind=8), intent(in) :: rho !< element density

    integer(kind=4), intent(in) :: nq3 !< number of 3D quadrature nodes
    integer(kind=4), intent(in) :: Np !< number of element dofs per dimension
    real(kind=8), intent(in) :: Jdet !< transformation determinant
    real(kind=8), intent(in) :: lambda      !< Lamé 1st parameter
    real(kind=8), intent(in) :: mu          !< Lamé 2nd parameter
    real(kind=8), dimension(4,nq3), intent(in) :: nodtet3 !< tetrahedron quadrature nodes
    real(kind=8), dimension(nq3), intent(in) :: weitet3 !< tetrahedron 3D quadrature weights
    real(kind=8), dimension(Np,nq3), intent(in) :: phi  !< element basis function
    real(kind=8), dimension(DIM,DIM+1), intent(in) :: Fk  !< coefficients for map from reference to physical tetrahedron
    ! rhs_tet_vol dim (3,Np)
    real(kind=8), dimension(:,:), intent(out) :: rhs_tet_vol !< element forcing term volume contribution

    integer(kind=4) :: q, i, j, k, m
    real(kind=8), dimension(DIM) :: points, forc_term

    rhs_tet_vol = 0.0

    ! time dependence
    ! loop on 3D quadrature nodes
    do q = 1,nq3

        do m=1,Np

            ! map the quadrature nodes from the reference tetrahedron to the physical tetrahedron
            do j=1,DIM
                points(j)=0.0
                do k=1,4
                    points(j) = points(j) + Fk(j,k)*nodtet3(k,q)
                end do
            end do

            ! f_time is provided by problem_data_and_properties.f90
            !! considers correct density if dynamic problem, otherwise rho=0.0
            ! forc_term = f(lambda,mu,points)
            forc_term = f_time(lambda,mu,points,rho)

            do i=1,DIM
                rhs_tet_vol(i,m) = rhs_tet_vol(i,m) + abs(Jdet)*weitet3(q)*forc_term(i)*phi(m,q)
            enddo

        end do

    end do


end subroutine MAKE_RHS_VOLUME

!> Assemble the element rhs term `rhs_tet_face` approximating the integral on the boundary faces of the tetrahedron (boundary conditions)
!> \f[ [F_{\partial K}]_i = \sum_{F \in {\mathcal{F}_h^{N}}|_{K} } \int_F \boldsymbol{g}_N \cdot \boldsymbol{\varphi}_{i,K} +  \theta \sum_{F \in \mathcal{F}_h^D |_{K} } \int_F \{\boldsymbol{\sigma}(\boldsymbol{g}_D) \} \cdot [\![ \boldsymbol{\varphi}_{i,K} ]\!] + \sum_{F\in \mathcal{F}_h^D |_{K} } \int_F \eta [\![ \boldsymbol{g}_D ]\!] \cdot [\![ \boldsymbol{\varphi}_{i,K} ]\!]  \f]
subroutine MAKE_RHS_FACE(theta, alpha, p, Np, e, E2, hk_1, hk_2, normal, area, Fk, nodtria2, weitria2, nq2, lambda, mu, node_maps, &
                            phi_b, grad_b, space_fun_tag, rhs_tet_face)

    ! theta and alpha are provided by the subroutine set_properties in problem_data_and_properties.f90
    ! phi_b and grad_b are provided by the subroutine basis_boundary in basis_functions.f90
    ! weitria2 is provided by the subroutine mapping_quadrature_2D in Poly_ref_mappings.f90
    ! nq2 is provided by the subroutine quadrature in basis_function.f90
    ! Fk is provided by by the subroutine jacobians in Poly_ref_mappings.f90


    integer(kind=4), intent(in) :: nq2  !< number of 2D quadrature nodes
    integer(kind=4), intent(in) :: Np   !< number of element dofs per dimension
    integer(kind=4), intent(in) :: p    !< element polynomial degree
    integer(kind=4), intent(in) :: e    !< current face
    integer(kind=4), intent(in) :: E2   !< E- neighbour element ID
    real(kind=8), intent(in) :: theta   !< penalty method (SIP, NIP, IIP)
    real(kind=8), intent(in) :: alpha   !< penalty constant
    real(kind=8), intent(in) :: hk_1    !< E+ element diameter
    real(kind=8), intent(in) :: hk_2    !< E- neighbour element diameter
    real(kind=8), intent(in) :: area    !< element area
    real(kind=8), intent(in) :: lambda  !< Lamé 1st parameter
    real(kind=8), intent(in) :: mu      !< Lamé 2nd parameter
    real(kind=8), dimension(4,4,4), intent(in) :: node_maps !< map from reference triangle to tetrahedral face
    real(kind=8), dimension(4,nq2), intent(in) :: nodtria2  !< 2D quatrature nodes for the reference triangle
    real(kind=8), dimension(nq2), intent(in) :: weitria2  !< 2D quatrature weights for the reference triangle
    real(kind=8), dimension(DIM), intent(in) :: normal  !< face normal vector
    real(kind=8), dimension(Np,nq2,2), intent(in) :: phi_b !< basis on the boundary
    real(kind=8), dimension(3,Np,nq2,2), intent(in) :: grad_b !< basis gradient on the boundary
    real(kind=8), dimension(DIM,DIM+1), intent(in) :: Fk !< coefficients for map from reference to physical tetrahedron
    integer(kind=4), intent(in) :: space_fun_tag    !< tag corresponding to specific function for boundary conditions, see problem_and_data_properties.f90
    real(kind=8), dimension(DIM,Np), intent(out) :: rhs_tet_face !< surface integral contributions to element rhs

    integer(kind=4) :: q, i, j, k, t, m
    real(kind=8) :: sigma, D_bar
    real(kind=8), dimension(2) :: val
    real(kind=8), dimension(DIM) :: points, diri_data, neum_data
    real(kind=8), dimension(DIM,Np) :: temp, temp2

    D_bar = lambda + 2*mu ! harmonic average of lambda+2*mu

    ! evaluation of the penalization function
    if(E2 == -1) then
        sigma = alpha*(p**2) / hk_1  * D_bar
    endif
    if(E2 /= -1 .and. E2 /= -2) then
        val(1) = hk_1
        val(2) = hk_2
        sigma = alpha*(p**2) / minval(val) * D_bar
    end if

    rhs_tet_face = 0.0

    ! loop on 2D quadrature nodes
    do q = 1,nq2

        temp = 0.0
        temp2 = 0.0

        ! if the condition is satisfied, then e is a Dirichlet boundary face
        if (E2 == -1) then

            do m=1,Np

                ! map from the 2D quadrature nodes to the 3D points of the physical tetrahedron
                do j=1,DIM
                    points(j) = 0.0
                    do k=1,4
                        do t=1,4
                            points(j) = points(j) + Fk(j,k)*node_maps(k,t,e)*nodtria2(t,q)
                        end do
                    end do
                end do

                diri_data = gd(points,space_fun_tag)

                temp(1,m) = (lambda + 2*mu)*grad_b(1,m,q,1)*diri_data(1)*normal(1) + &
                            mu*grad_b(2,m,q,1)*diri_data(1)*normal(2) + mu*grad_b(3,m,q,1)*diri_data(1)*normal(3) + &
                            lambda*grad_b(1,m,q,1)*diri_data(2)*normal(2) + mu*grad_b(2,m,q,1)*diri_data(2)*normal(1) + &
                            lambda*grad_b(1,m,q,1)*diri_data(3)*normal(3) + mu*grad_b(3,m,q,1)*diri_data(3)*normal(1)
                temp(2,m) = (lambda + 2*mu)*grad_b(2,m,q,1)*diri_data(2)*normal(2) + &
                            mu*grad_b(1,m,q,1)*diri_data(2)*normal(1) + mu*grad_b(3,m,q,1)*diri_data(2)*normal(3) + &
                            lambda*grad_b(2,m,q,1)*diri_data(1)*normal(1) + mu*grad_b(1,m,q,1)*diri_data(1)*normal(2) + &
                            lambda*grad_b(2,m,q,1)*diri_data(3)*normal(3) + mu*grad_b(3,m,q,1)*diri_data(3)*normal(2)
                temp(3,m) = (lambda + 2*mu)*grad_b(3,m,q,1)*diri_data(3)*normal(3) + &
                            mu*grad_b(1,m,q,1)*diri_data(3)*normal(1) + mu*grad_b(2,m,q,1)*diri_data(3)*normal(2) + &
                            lambda*grad_b(3,m,q,1)*diri_data(1)*normal(1) + mu*grad_b(1,m,q,1)*diri_data(1)*normal(3) + &
                            lambda*grad_b(3,m,q,1)*diri_data(2)*normal(2) + mu*grad_b(2,m,q,1)*diri_data(2)*normal(3)

                temp2(1,m) = phi_b(m,q,1)*diri_data(1) * (abs(normal(1))**2 + 0.5*abs(normal(2))**2 + 0.5*abs(normal(3))**2) + &
                                    0.5*phi_b(m,q,1)*diri_data(2)*normal(1)*normal(2) + 0.5*phi_b(m,q,1)*diri_data(3)*normal(1)*normal(3)
                temp2(2,m) = phi_b(m,q,1)*diri_data(2) * (abs(normal(2))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(3))**2) + &
                                    0.5*phi_b(m,q,1)*diri_data(1)*normal(1)*normal(2) + 0.5*phi_b(m,q,1)*diri_data(3)*normal(2)*normal(3)
                temp2(3,m) = phi_b(m,q,1)*diri_data(3) * (abs(normal(3))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(2))**2) + &
                                    0.5*phi_b(m,q,1)*diri_data(1)*normal(1)*normal(3) + 0.5*phi_b(m,q,1)*diri_data(2)*normal(2)*normal(3)

                do i=1,DIM
                    rhs_tet_face(i,m) = rhs_tet_face(i,m) &
                                    + theta*weitria2(q)*temp(i,m)*area &
                                    + sigma*weitria2(q)*temp2(i,m)*area
                enddo

            end do

        endif

        ! if the condition is satisfied, then e is a Neumann boundary face
        if (E2 == -2) then

            do m=1,Np

                ! map from the 2D quadrature nodes to the 3D points of the physical tetrahedron
                do j=1,DIM
                    points(j) = 0.0
                    do k=1,4
                        do t=1,4
                            points(j) = points(j) + Fk(j,k)*node_maps(k,t,e)*nodtria2(t,q)
                        end do
                    end do
                end do

                neum_data = gn(lambda,mu,normal,points,space_fun_tag)

                do i=1,DIM
                    rhs_tet_face(i,m) = rhs_tet_face(i,m) &
                                    + weitria2(q)*phi_b(m,q,1)*neum_data(i)*area
                enddo

            end do

        endif

    end do

end subroutine MAKE_RHS_FACE

!> Assemble the terms of the stiffness matrix  `S_E1` and `I_E1` approximating the integrals on the faces of the tetrahedron
!> and the ones `S_E2` and `I_E2` approximating the integrals on the faces of the neighbouring tetrahedron E-
!> \f[ [I_K]_{ij} = \sum_{F\in (F_h^I \cup F_h^D)|_{K} } \int_F  \{ \boldsymbol{\sigma}(\boldsymbol{\varphi}_{j,K}) \} \cdot [\![ \boldsymbol{\varphi}_{i,K} ]\!] \quad  [S_K]_{ij} =  \sum_{F\in (F_h^I \cup F_h^D) |_{K} } \int_F \eta [\![ \boldsymbol{\varphi}_{j,K} ]\!] \cdot [\![ \boldsymbol{\varphi}_{i,K} ]\!] \f]
!> components
!> \f[ [S_{K+}]_{ij} = \sum_{F\in F_h^I|_{K} } \int_F \eta (\boldsymbol{\varphi}_{j,K}^+ \textbf{n}^+):(\boldsymbol{\varphi}_{i,K}^+ \textbf{n}^+ ) + \sum_{F\in F_h^D|_{K} } \int_F \eta (\boldsymbol{\varphi}_{j,K}^+ \textbf{n}^+):(\boldsymbol{\varphi}_{i,K}^+ \textbf{n}^+) \f]
!> \f[ [S_{K-}]_{ij} = \sum_{F\in F_h^I|_{K} } \int_F \eta (\boldsymbol{\varphi}_{j,K}^+ \textbf{n}^+):( \boldsymbol{\varphi}_{i,K}^- \textbf{n}^-) \f]
!> \f[ [I_{K+}]_{ij} = \sum_{F\in F_h^I|_{K} } \int_F (\boldsymbol{\varphi}_{j,K}^+ \textbf{n}^+ ) : \frac 12 \boldsymbol{\sigma}(\boldsymbol{\varphi}_{i,K}^+) + \sum_{F\in F_h^D|_{K} } \int_F (\boldsymbol{\varphi}_{j,K}^+ \textbf{n}^+ ) : \frac 12 \boldsymbol{\sigma}(\boldsymbol{\varphi}_{i,K}^+) \f]
!> \f[ [I_{K-}]_{ij} = \sum_{F\in F_h^I|_{K} } \int_F ( \boldsymbol{\varphi}_{j,K}^- \textbf{n}^- ) : \frac 12 \boldsymbol{\sigma}(\boldsymbol{\varphi}_{i,K}^+) \f]

subroutine MAKE_STIFFNESS_FACE(alpha, p, Np, E2, hk_1, hk_2, normal, area, weitria2, nq2, lambda, mu, &
                            phi_b, grad_b, S_E1, I_E1, I_E2, S_E2)

    ! theta and alpha are provided by the subroutine set_properties in problem_data_and_properties.f90
    ! phi_b and grad_b are provided by the subroutine basis_boundary in basis_functions.f90
    ! weitria2 is provided by the subroutine mapping_quadrature_2D in Poly_ref_mappings.f90
    ! nq2 is provided by the subroutine quadrature in basis_function.f90

    integer(kind=4), intent(in) :: nq2  !< number of 2D quadrature nodes
    integer(kind=4), intent(in) :: Np   !< element dof per dimension
    integer(kind=4), intent(in) :: p    !< element polynomial degree
    integer(kind=4), intent(in) :: E2   !< E- neighbour element ID
    real(kind=8), intent(in) :: alpha   !< penaly constant
    real(kind=8), intent(in) :: hk_1    !< E+ diameter
    real(kind=8), intent(in) :: hk_2    !< E- diameter
    real(kind=8), intent(in) :: area    !< element area
    real(kind=8), intent(in) :: lambda  !< Lamé 1st parameter
    real(kind=8), intent(in) :: mu      !< Lamé 2nd parameter
    real(kind=8), dimension(nq2), intent(in) :: weitria2    !< 2d weights
    real(kind=8), dimension(DIM), intent(in) :: normal      !< normal vector
    real(kind=8), dimension(Np,nq2,2), intent(in) :: phi_b  !< basis on the boundary
    real(kind=8), dimension(3,Np,nq2,2), intent(in) :: grad_b   !< basis gradient on the boundary
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: S_E1    !< element stiffness matrix stabilization S, E+ contribution
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: I_E1         !< element stiffness matrix interior flux I, E+ contribution
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: I_E2        !< element stiffness matrix interior flux I, E- contribution
    real(kind=8), dimension(DIM,DIM,Np,Np), intent(out) :: S_E2        !< element stiffness matrix stabilization S, E- contribution

    real(kind=8), dimension(DIM,DIM,Np,Np) :: temp, temp1, temp2
    real(kind=8), dimension(2) :: val
    real(kind=8) :: sigma, D_bar
    integer(kind=4) :: q, i, j, m, n

    D_bar = lambda + 2*mu ! harmonic average of lambda+2*mu

    ! evaluation of the penalization function
    if (E2 == -1) then
        sigma = alpha*(p**2) / hk_1 * D_bar
    endif
    if(E2 /= -1 .and. E2 /= -2) then
        val(1) = hk_1
        val(2) = hk_2
        sigma = alpha*(p**2) / minval(val) * D_bar
    end if

    S_E1 = 0.0
    I_E1 = 0.0
    I_E2 = 0.0
    S_E2 = 0.0

    ! check if the actual face is not a Neumann boundary face
    if(E2 /= -2) then

        temp = 0.0

        ! loop on 2D quadrature nodes
        do q = 1,nq2

            do m=1,Np
                do n=1,Np

                    temp(1,1,m,n) = phi_b(m,q,1)*phi_b(n,q,1) * (abs(normal(1))**2 + 0.5*abs(normal(2))**2 + 0.5*abs(normal(3))**2)
                    temp(2,2,m,n) = phi_b(m,q,1)*phi_b(n,q,1) * (abs(normal(2))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(3))**2)
                    temp(3,3,m,n) = phi_b(m,q,1)*phi_b(n,q,1) * (abs(normal(3))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(2))**2)
                    do i=1,2
                        do j=i+1,3
                            temp(i,j,m,n) = 0.5*phi_b(m,q,1)*phi_b(n,q,1)*normal(i)*normal(j)
                            temp(j,i,m,n) = 0.5*phi_b(m,q,1)*phi_b(n,q,1)*normal(i)*normal(j)
                        enddo
                    enddo

                    do i=1,DIM
                        do j=1,DIM
                            S_E1(i,j,m,n) = S_E1(i,j,m,n) + sigma*weitria2(q)*temp(i,j,m,n)*area
                        enddo
                    enddo

                enddo
            enddo

            ! if the condition is satisfied, then e is a Dirichlet boundary face
            if (E2 == -1) then

                temp = 0.0

                do m=1,Np
                    do n=1,Np

                        temp(1,1,m,n) = (lambda + 2*mu)*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + mu*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3)
                        temp(2,2,m,n) = (lambda + 2*mu)*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + &
                                        mu*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1) + mu*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3)
                        temp(3,3,m,n) = (lambda + 2*mu)*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + mu*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1)
                        do i=1,2
                            do j=i+1,3
                                temp(i,j,m,n) = lambda*grad_b(i,m,q,1)*phi_b(n,q,1)*normal(j) + mu*grad_b(j,m,q,1)*phi_b(n,q,1)*normal(i)
                                temp(j,i,m,n) = lambda*grad_b(j,m,q,1)*phi_b(n,q,1)*normal(i) + mu*grad_b(i,m,q,1)*phi_b(n,q,1)*normal(j)
                            enddo
                        enddo

                        do i=1,DIM
                            do j=1,DIM
                                I_E1(i,j,m,n) = I_E1(i,j,m,n) + weitria2(q)*area*temp(i,j,m,n)
                            enddo
                        enddo

                    enddo
                enddo

            else

                temp = 0.0
                temp2 = 0.0
                temp1 = 0.0

                do m=1,Np
                    do n=1,Np

                        temp(1,1,m,n) = (lambda + 2*mu)*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + mu*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3)
                        temp(2,2,m,n) = (lambda + 2*mu)*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + &
                                        mu*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1) + mu*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3)
                        temp(3,3,m,n) = (lambda + 2*mu)*grad_b(3,m,q,1)*phi_b(n,q,1)*normal(3) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,1)*normal(2) + mu*grad_b(1,m,q,1)*phi_b(n,q,1)*normal(1)
                        do i=1,2
                            do j=i+1,3
                                temp(i,j,m,n) = lambda*grad_b(i,m,q,1)*phi_b(n,q,1)*normal(j) + mu*grad_b(j,m,q,1)*phi_b(n,q,1)*normal(i)
                                temp(j,i,m,n) = lambda*grad_b(j,m,q,1)*phi_b(n,q,1)*normal(i) + mu*grad_b(i,m,q,1)*phi_b(n,q,1)*normal(j)
                            enddo
                        enddo
                        do i=1,DIM
                            do j=1,DIM
                                I_E1(i,j,m,n) = I_E1(i,j,m,n) + 0.5*weitria2(q)*area*temp(i,j,m,n)
                            enddo
                        enddo

                        temp1(1,1,m,n) = (lambda + 2*mu)*grad_b(1,m,q,1)*phi_b(n,q,2)*normal(1) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,2)*normal(2) + mu*grad_b(3,m,q,1)*phi_b(n,q,2)*normal(3)
                        temp1(2,2,m,n) = (lambda + 2*mu)*grad_b(2,m,q,1)*phi_b(n,q,2)*normal(2) + &
                                        mu*grad_b(1,m,q,1)*phi_b(n,q,2)*normal(1) + mu*grad_b(3,m,q,1)*phi_b(n,q,2)*normal(3)
                        temp1(3,3,m,n) = (lambda + 2*mu)*grad_b(3,m,q,1)*phi_b(n,q,2)*normal(3) + &
                                        mu*grad_b(2,m,q,1)*phi_b(n,q,2)*normal(2) + mu*grad_b(1,m,q,1)*phi_b(n,q,2)*normal(1)
                        do i=1,2
                            do j=i+1,3
                                temp1(i,j,m,n) = lambda*grad_b(i,m,q,1)*phi_b(n,q,2)*normal(j) + mu*grad_b(j,m,q,1)*phi_b(n,q,2)*normal(i)
                                temp1(j,i,m,n) = lambda*grad_b(j,m,q,1)*phi_b(n,q,2)*normal(i) + mu*grad_b(i,m,q,1)*phi_b(n,q,2)*normal(j)
                            enddo
                        enddo
                        do i=1,DIM
                            do j=1,DIM
                                I_E2(i,j,m,n) = I_E2(i,j,m,n) - 0.5*weitria2(q)*area*temp1(i,j,m,n)
                            enddo
                        enddo

                        temp2(1,1,m,n) = phi_b(m,q,1)*phi_b(n,q,2) * (abs(normal(1))**2 + 0.5*abs(normal(2))**2 + 0.5*abs(normal(3))**2)
                        temp2(2,2,m,n) = phi_b(m,q,1)*phi_b(n,q,2) * (abs(normal(2))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(3))**2)
                        temp2(3,3,m,n) = phi_b(m,q,1)*phi_b(n,q,2) * (abs(normal(3))**2 + 0.5*abs(normal(1))**2 + 0.5*abs(normal(2))**2)
                        do i=1,2
                            do j=i+1,3
                                temp2(i,j,m,n) = 0.5*phi_b(m,q,1)*phi_b(n,q,2)*normal(i)*normal(j)
                                temp2(j,i,m,n) = 0.5*phi_b(m,q,1)*phi_b(n,q,2)*normal(i)*normal(j)
                            enddo
                        enddo

                        do i=1,DIM
                            do j=1,DIM
                                S_E2(i,j,m,n) = S_E2(i,j,m,n) - sigma*weitria2(q)*temp2(i,j,m,n)*area
                            enddo
                        enddo

                    enddo
                enddo

            endif

        enddo

    endif

end subroutine MAKE_STIFFNESS_FACE

!> @brief Assemble the local vector term vec_loc approximating the integral on the tetrahedron
subroutine MAKE_VECTOR_TET(Np, Fk, Jdet, nodtet3, weitet3, nq3, phi, vec_loc, f_analytic)

    ! phi is provided by the subroutine basis in basis_function.f90
    ! weitet3 is provided by the subroutine mapping_quadrature_3D in Poly_ref_mappings.f90
    ! nq3 is provided by the subroutine quadrature in basis_function.f90
    ! Fk and Jdet are provided by the subroutine jacobians in Poly_ref_mappings.f90

    ! pass function as argument
    interface
        function f_analytic(point) result(res)
            use global_parameters
            real(kind=8), dimension(DIM) :: point, res
        end function f_analytic
    end interface

    integer(kind=4), intent(in) :: nq3
    integer(kind=4), intent(in) :: Np
    real(kind=8), intent(in) :: Jdet
    real(kind=8), dimension(4,nq3), intent(in) :: nodtet3
    real(kind=8), dimension(nq3), intent(in) :: weitet3 !< tetrahedron 3D quadrature weights
    real(kind=8), dimension(Np,nq3), intent(in) :: phi
    real(kind=8), dimension(DIM,DIM+1), intent(in) :: Fk !< coefficients for map from reference to physical tetrahedron
    real(kind=8), dimension(DIM,Np), intent(out) :: vec_loc

    integer(kind=4) :: q, i, j, k, m
    real(kind=8), dimension(DIM) :: points, eval

    vec_loc = 0.0

    ! loop on 3D quadrature nodes
    do q = 1,nq3

        do m=1,Np

            ! map the quadrature nodes from the reference tetrahedron to the physical tetrahedron
            do j=1,DIM
                points(j)=0.0
                do k=1,4
                    points(j) = points(j) + Fk(j,k)*nodtet3(k,q)
                end do
            end do

            ! f_time is provided by problem_data_and_properties.f90
            eval = f_analytic(points)

            do i=1,DIM
                vec_loc(i,m) = vec_loc(i,m) + abs(Jdet)*weitet3(q)*eval(i)*phi(m,q)
            enddo

        end do

    end do

end subroutine MAKE_VECTOR_TET

end module assemble_element