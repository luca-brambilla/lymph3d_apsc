module basis_function

    use Poly_mesh

    implicit none

    contains

    ! Create a list blist of the degrees of monomials of the Np basis functions up to a total degree p
    subroutine basis_list(blist, p, Np)

        implicit none

        integer(kind=4), intent(in) :: p, Np
        integer(kind=4), dimension(Np,3), intent(out) :: blist

        integer(kind=4) :: ii
        integer(kind=4) :: q1, q2, q3

        ii = 1
        q1 = p

        blist = 0

        do while (q1 >= 0)
            q2 = p - q1
            do while (q2 >= 0)
                q3 = p - q1 - q2
                do while (q3 >= 0)
                    blist(ii,1) = q1
                    blist(ii,2) = q2
                    blist(ii,3) = q3
                    q3 = q3 - 1
                    ii = ii + 1
                end do
                q2 = q2 - 1
            end do
            q1 = q1 - 1
        end do

    end subroutine basis_list

    ! Compute Gauss-Legendre quadrature nodes and weights on the square (-1,1)^2 and cube (-1,1)^3
    subroutine quadrature(nod2, wei2, nod3, wei3, p, nq3, nq2)

        implicit none

        integer(kind=4), intent(in) :: p
        integer(kind=4), intent(out) :: nq3, nq2
        real(kind=8), dimension(:,:), allocatable, intent(out) :: nod2, nod3
        real(kind=8), dimension(:), allocatable, intent(out) :: wei2, wei3

        real(kind=8), dimension(:), allocatable :: nod1, wei1
        integer(kind=4) :: nqn ! number of quadrature nodes
        integer(kind=4) :: i, j, k, counter

        nqn = p+1 ! number of quadrature nodes in each direction

        nq3 = nqn**3
        nq2 = nqn**2

        ALLOCATE (nod1(nqn))
        ALLOCATE (wei1(nqn))
        ALLOCATE (nod2(4, nq2))
        ALLOCATE (wei2(nq2))
        ALLOCATE (nod3(4, nq3))
        ALLOCATE (wei3(nq3))

        ! Construction of GL nodes and weights in 1D on the interval (-1,1)
        ! They are hard-coded until degree 7
        if (nqn == 1) then
        ! Degree of exactness 1
            nod1 = [0.0]
            wei1 = [2.0]
        end if
        if (nqn == 2) then
        ! Degree of exactness 3
            nod1 = [sqrt(1.0/3.0), - sqrt(1.0/3.0)]
            wei1 = [1.0, 1.0]
        end if
        if (nqn == 3) then
        ! Degree of exactness 5
            nod1 = [0.0, sqrt(3.0/5.0), - sqrt(3.0/5.0)]
            wei1 = [8.0/9.0, 5.0/9.0, 5.0/9.0]
        endif
        if (nqn == 4) then
        ! Degree of exactness 7
            nod1 = [0.3399810435848563, -0.3399810435848563, 0.8611363115940526, -0.8611363115940526]
            wei1 = [(18+sqrt(30.0))/36.0, (18+sqrt(30.0))/36.0, (18-sqrt(30.0))/36.0, (18-sqrt(30.0))/36.0]
        end if
        if (nqn == 5) then
        ! Degree of exactness 9
            nod1 = [0.0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640]
            wei1 = [128.0/225.0, (322+13*sqrt(70.0))/900.0, (322+13*sqrt(70.0))/900.0, (322-13*sqrt(70.0))/900.0, (322-13*sqrt(70.0))/900.0]
        end if
        if (nqn == 6) then
        ! Degree of exactness 11
            nod1 = [0.6612093864662645, -0.6612093864662645, 0.2386191860831969, -0.2386191860831969, 0.9324695142031521, -0.9324695142031521]
            wei1 = [0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704]
        end if
        if (nqn == 7) then
        ! Degree of exactness 13
            nod1 = [0.0, 0.4058451513773972, -0.4058451513773972, 0.7415311855993945, -0.7415311855993945, 0.9491079123427585, -0.9491079123427585]
            wei1 = [0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697]
        end if

        if (nqn > 7) then
            ! Degree of exactness 2*nqn-1
            call GauLeg(-1, 1, nqn, nod1, wei1) ! see subroutine below
        endif

        ! Construction of GL nodes and weights in 2D
        counter = 0
        do i=1,nqn
            do j=1,nqn
                counter = counter + 1
                nod2(1,counter) = nod1(i)
                nod2(2,counter) = nod1(j)
                wei2(counter) = wei1(i) * wei1(j)
            end do
        end do

        ! Construction of GL nodes and weights in 3D
        counter = 0
        do i=1,nqn
            do j=1,nqn
                do k=1,nqn
                    counter = counter + 1
                    nod3(1,counter) = nod1(i)
                    nod3(2,counter) = nod1(j)
                    nod3(3,counter) = nod1(k)
                    wei3(counter) = wei1(i) * wei1(j) * wei1(k)
                enddo
            end do
        end do

        ! Ones in the quadrature nodes are needed by the translation vector of Fk and node_maps (see basis and basis_boundary below)
        nod2(3,:) = 1.0
        nod2(4,:) = 1.0
        nod3(4,:) = 1.0

    end subroutine quadrature

    ! Compute n Gauss-Legendre quadrature nodes and weights on a given interval (a,b)
    subroutine GauLeg(a, b, n, x_GL, w_GL)

        use vet_mat_operations
        use global_parameters

        implicit none

        integer(kind=4), intent(in) :: n
        integer(kind=4), intent(in) :: a, b
        real(kind=8), dimension(n), intent(out) :: x_GL, w_GL

        integer(kind=4) :: m, j
        real(kind=8) :: Err
        real(kind=8), dimension(:), allocatable :: zeros
        real(kind=8), dimension(:), allocatable :: z, z_old, flip_z
        real(kind=8), dimension(:), allocatable :: ErrV
        real(kind=8), dimension(:), allocatable :: p_app, flip_p_app
        real(kind=8), dimension(:,:), allocatable :: p

        ! Number of the nodes to compute (before reflection)
        if (MOD(n, 2) == 0) then
          m = (n + 1) / 2 + 1
        else
          m = (n + 1) / 2
        end if

        allocate(z(m))
        allocate(flip_z(m))
        allocate(z_old(m))
        allocate(ErrV(m))
        allocate(p(m,3))
        allocate(p_app(m))
        allocate(flip_p_app(m))
        allocate(zeros(m))

        ! Computation of the function
        do j=1,m
          z(j) = cos( pi * ((j-0.25) / (n+0.5) ))
        enddo

        ! Error initialization in scalar and vectorial form
        Err = epsilon(1.0d0) + 1.0
        ErrV = Err * 1.0

        zeros = 0.0
        do while (Err > epsilon(1.0d0))
          p = 0.0
          p(:,1) = 1.0

          do j = 1, n
            p(:,2:3) = p(:,1:2)
            p(:, 1) = ((2.0 * j - 1.0) / j) * z * p(:, 2) - ((j - 1.0) / j) * p(:, 3)
          end do

          ! Update of the new step solution
          p_app = n * (z * p(:, 1) - p(:, 2)) / (z**2 - 1.0)
          z_old = z
          z = z_old - merge((p(:, 1) / p_app), zeros, ErrV > epsilon(1.0d0))

          ! Compute the new errors
          ErrV = abs(z - z_old)
          Err = maxval(ErrV)

        end do

        ! see vet_mat_operations.f90
        call flip_vector(z,m,flip_z)
        call flip_vector(p_app,m,flip_p_app)

        ! Computation of the GL nodes in the requested interval
        if (MOD(n, 2) == 0) then
          x_GL = ((b + a) / 2.0) - ((b - a) / 2.0) * (/ z(1:(m-1)), -flip_z(2:m) /)
        else
          x_GL = ((b + a) / 2.0) - ((b - a) / 2.0) * (/ z(1:(m-1)), -flip_z /)
        end if

        ! Reflection of z and p_app
        if (MOD(n, 2) == 0) then
          z = (/ z(1:(m-1)), flip_z(2:m) /)
          p_app = (/ p_app(1:(m-1)), flip_p_app(2:m) /)
        else
          z = (/ z(1:(m-1)), flip_z /)
          p_app = (/ p_app(1:(m-1)), flip_p_app /)
        end if

        ! Computation of the GL weights in the requested interval
        w_GL = (b - a) / ((1.0 - z**2) * (p_app**2))

        deallocate(z)
        deallocate(flip_z)
        deallocate(z_old)
        deallocate(ErrV)
        deallocate(p)
        deallocate(p_app)
        deallocate(flip_p_app)

    end subroutine GauLeg

    ! Evaluate the non-normalized one-dimensional Legendre Polynomial on the interval int at points x for order p
    subroutine LegendreP_nonnorm(LP_nonnorm, x, p, int, nq)

        implicit none

        integer(kind=4), intent(in) :: p, nq
        real(kind=8), dimension(2), intent(in) :: int
        real(kind=8), dimension(nq), intent(in) :: x
        real(kind=8), dimension(nq), intent(out) :: LP_nonnorm ! L_p not normalized

        integer(kind=4) :: ii, start
        real(kind=8), dimension(nq) :: xp
        real(kind=8), dimension(nq) :: LP_old, LP_oold
        real(kind=8), dimension(nq) :: LP_temp ! vector which collects the non-normalized Legendre Polynomial L_p evaluated at each point of x
        real(kind=8) :: hb, mb  ! hb: half of the length of the interval int
                                ! mb: midpoint of the interval int

        ! Initialization
        LP_temp = 0.0
        LP_nonnorm = 0.0

        hb = (int(2) - int(1)) / 2.0
        mb = (int(2) + int(1)) / 2.0

        ! affine map from the Cartesian bounding box B_E for a polygon to the reference hypercube B_hat
        xp = (x - mb) / hb

        ! Computation of Legendre Polynomials for each specific degree
        ! (they are hard-coded until the degree 7)
        if (p == 0) then
            LP_temp = 1.0
        endif
        if (p == 1) then
            LP_temp = xp
        endif
        if (p == 2) then
            LP_temp = (3*xp**2 - 1) / 2
        endif
        if (p == 3) then
            LP_temp = (5*xp**3 - 3*xp) / 2
        endif
        if (p == 4) then
            LP_temp = (35*xp**4 - 30*xp**2 + 3) / 8
        endif
        if (p == 5) then
            LP_temp = (63*xp**5 - 70*xp**3 + 15*xp) / 8
        endif
        if (p == 6) then
            LP_temp = (231*xp**6 - 315*xp**4 + 105*xp**2 - 5) / 16
        endif
        if (p == 7) then
            LP_temp = (429*xp**7 - 693*xp**5 + 315*xp**3 - 35*xp) / 16
        endif

        ! From degree 7 on, compute Legendre Polynomials by recursion
        if (p > 7) then

            ! Computation of Legendre Polynomials for degree greater than 7
            LP_old = (429*xp**7 - 693*xp**5 + 315*xp**3 - 35*xp) / 16
            LP_oold = (231*xp**6 - 315*xp**4 + 105*xp**2 - 5) / 16
            start = 7
            do ii = start, p - 1
                LP_temp = 1.0 / (ii + 1.0) * ((2*ii + 1)* xp*LP_old - ii*LP_oold)
                LP_oold = LP_old
                LP_old = LP_temp
            end do

        end if

        LP_nonnorm = LP_temp

    end subroutine LegendreP_nonnorm

    ! Evaluate the normalized one-dimensional Legendre Polynomial L_p on the interval int at points x for order p
    subroutine LegendreP(LP, x, p, int, nq)

        implicit none

        integer(kind=4), intent(in) :: p, nq
        real(kind=8), dimension(2), intent(in) :: int
        real(kind=8), dimension(nq), intent(in) :: x
        real(kind=8), dimension(nq), intent(out) :: LP ! L_p

        real(kind=8) :: hb
        real(kind=8), dimension(nq) :: LP_nonnorm ! vector which collects the non-normalized Legendre Polynomial L_p evaluated at each point of x

        ! Initialization
        LP_nonnorm = 0.0
        LP = 0.0

        hb = (int(2) - int(1)) / 2.0

        call LegendreP_nonnorm(LP_nonnorm, x, p, int, nq)

        ! Normalization to extract L_p
        LP = LP_nonnorm / sqrt(hb)

    end subroutine LegendreP

    ! Evaluate the normalized derivative L'_p of the one-dimensional Legendre Polynomial L_p on the interval int at points x for order p
    subroutine GradLegendreP(LPder, x, p, int, nq)

        implicit none

        integer(kind=4), intent(in) :: p, nq
        real(kind=8), dimension(2), intent(in) :: int
        real(kind=8), dimension(nq), intent(in) :: x
        real(kind=8), dimension(nq), intent(out) :: LPder ! L'_p

        integer(kind=4) :: ii
        real(kind=8), dimension(nq) :: xp
        real(kind=8), dimension(nq) :: LP_old
        real(kind=8), dimension(nq) :: LPder_temp ! vector which collects the Legendre Polynomial L_p evaluated at each point of x
        real(kind=8) :: hb, mb ! hb_inv: half of the length of the interval int
                               ! mb: midpoint of the interval int

        ! Initialization
        LPder_temp = 0.0
        LPder = 0.0

        hb = (int(2) - int(1)) / 2.0
        mb = (int(2) + int(1)) / 2.0

        ! affine map from the Cartesian bounding box B_E for a polygon to the reference hypercube B_hat
        xp = (x - mb) / hb

        ! Computation of the derivative Legendre Polynomials for each specific degree
        ! (they are hard-coded until the degree 7)
        if (p == 1) then
            LPder_temp = 1.0
        endif
        if (p == 2) then
            LPder_temp = 3 * xp
        endif
        if (p == 3) then
            LPder_temp = (15 * xp**2 - 3) / 2
        endif
        if (p == 4) then
            LPder_temp = (35 * xp**3 - 15 * xp) / 2
        endif
        if (p == 5) then
            LPder_temp = (315 * xp**4 - 210 * xp**2 + 15) / 8
        endif
        if (p == 6) then
            LPder_temp = (1386 * xp**5 - 1260 * xp**3 + 210 * xp - 5) / 16
        endif
        if (p == 7) then
            LPder_temp = (3003 * xp**6 - 3465 * xp**4 + 945 * xp**2 - 35) / 16
        endif

        ! From degree 7 on, compute the derivative of Legendre Polynomials by recursion
        if (p > 7) then

            ! Computation of the derivative of Legendre Polynomials for degree greater than 7
            ii = p-1

            do while (ii .ge. 0)
                call LegendreP_nonnorm(LP_old, x, ii, int, nq)
                LPder_temp = LPder_temp + (2*ii + 1) * LP_old
                ii = ii - 2
            end do

        endif

        ! Normalization to extract L'_p
        LPder = LPder_temp / ((hb) * sqrt(hb))

    end subroutine GradLegendreP

    ! Evaluate the basis functions and their partial derivatives at the 3D quadrature nodes for a given polyhedral element contained in b_box
    subroutine basis(phi, dphi, b_box, Np, blist, Fk, nodtet3, nq3)

        implicit none

        integer(kind=4), intent(in) :: nq3, Np
        real(kind=8), dimension(3,2), intent(in) :: b_box
        real(kind=8), dimension(3,4), intent(in) :: Fk
        integer(kind=4), dimension(Np,3), intent(in) :: blist
        real(kind=8), dimension(4,nq3), intent(in) :: nodtet3
        real(kind=8), dimension(Np,nq3), intent(out) :: phi
        real(kind=8), dimension(3,Np,nq3), intent(out) :: dphi

        real(kind=8), dimension(nq3) :: x_p, y_p, z_p
        real(kind=8), dimension(2) :: intx, inty, intz
        real(kind=8), dimension(3,nq3) :: pt
        real(kind=8), dimension(nq3) :: valx, valy, valz
        real(kind=8), dimension(nq3) :: dvalx, dvaly, dvalz
        integer(kind=4) :: j, q, f, l

        ! initialization of phi and dphi
        do j=1,Np
            do q=1,nq3
                phi(j,q)=0.0 ! j-th basis function evaluated in the q-th 3D quadrature node
                do l=1,3
                    dphi(l,j,q) = 0.0 ! partial derivative, w.r.t. the l-th coordinate, of the j-th basis function evaluated in the q-th 3D quadrature node
                end do
            end do
        end do

        ! map the quadrature nodes from the reference tetrahedron to the physical tetrahedron
        do j=1,3
            do q=1,nq3
                pt(j,q)=0.0
                do l=1,4
                    pt(j,q) = pt(j,q) + Fk(j,l) * nodtet3(l,q) ! j-th coordinate of the q-th 3D quadrature node mapped to the physical tetrahedron
                end do
            end do
        end do

        ! collect the coordinates of pt into different vectors
        x_p = pt(1,:)
        y_p = pt(2,:)
        z_p = pt(3,:)

        ! definition of the intervals intx, inty and intz
        intx = b_box(1,:)
        inty = b_box(2,:)
        intz = b_box(3,:)

        ! loop on the basis functions
        do f = 1,Np

            ! initialization of valx, dvalx, valy, dvaly, valz and dvalz
            valx = 0.0
            dvalx = 0.0
            valy = 0.0
            dvaly = 0.0
            valz = 0.0
            dvalz = 0.0

            ! evaluation of Legendre Polynomials L_{blist(f,j)}, j=1,2,3, and their derivatives at x_p, y_p and z_p respectively
            call LegendreP(valx, x_p, blist(f,1), intx, nq3)
            call GradLegendreP(dvalx, x_p, blist(f,1), intx, nq3)
            call LegendreP(valy, y_p, blist(f,2), inty, nq3)
            call GradLegendreP(dvaly, y_p, blist(f,2), inty, nq3)
            call LegendreP(valz, z_p, blist(f,3), intz, nq3)
            call GradLegendreP(dvalz, z_p, blist(f,3), intz, nq3)

            ! evaluation of the f-th basis function and their partial derivatives in the q-th 3D quadrature node (mapped to the physical tetrahedron)
            phi(f,:) = valx*valy*valz
            dphi(1,f,:) = dvalx*valy*valz
            dphi(2,f,:) = valx*dvaly*valz
            dphi(3,f,:) = valx*valy*dvalz

        end do

    end subroutine basis

    ! This function evaluates the basis functions for every face of two neighbouring tetrahedra with global id E1 and E2 at the 2D quadrature nodes
    ! contained respectively in bounding boxes b_box1 and b_box2
    subroutine basis_boundary(phi_b, grad_b, e_E1, E2, b_box1, b_box2, blist, Np, Fk, node_maps, nodtria2, nq2)

        integer(kind=4), intent(in) :: nq2, Np
        integer(kind=4), intent(in) :: e_E1, E2
        real(kind=8), dimension(3,2), intent(in) :: b_box1, b_box2
        integer(kind=4), dimension(Np,3), intent(in) :: blist
        real(kind=8), dimension(3,4), intent(in) :: Fk
        real(kind=8), dimension(4,nq2), intent(in) :: nodtria2
        real(kind=8), dimension(4,4,4), intent(in) :: node_maps
        real(kind=8), dimension(Np,nq2,2), intent(out) :: phi_b
        real(kind=8), dimension(3,Np,nq2,2), intent(out) :: grad_b

        real(kind=8), dimension(3,nq2) :: pt
        real(kind=8), dimension(3,4) :: temp
        real(kind=8), dimension(nq2) :: x_p, y_p, z_p
        real(kind=8), dimension(2) :: intx, inty, intz
        real(kind=8), dimension(nq2) :: valx, valy, valz
        real(kind=8), dimension(nq2) :: dvalx, dvaly, dvalz
        integer(kind=4) :: j, q, l, f

        ! initialization of phi_b and grad_b
        do j=1,Np
            do q=1,nq2
                phi_b(j,q,1)=0.0 ! j-th basis function, restricted in E1, evaluated in the q-th 2D quadrature node
                phi_b(j,q,2)=0.0
                do l=1,3
                    grad_b(l,j,q,1) = 0.0 ! partial derivative w.r.t. the l-th component of the j-th basis function, restricted in E1, evaluated in the q-th 2D quadrature node
                    grad_b(l,j,q,2) = 0.0
                end do
            end do
        end do

        ! map from the face of index e_E1 of the reference tetrahedron to the corresponding face of the physical tetrahedron
        do j=1,3
            do q=1,4
                temp(j,q)=0.0
                do l=1,4
                    temp(j,q) = temp(j,q) + Fk(j,l)*node_maps(l,q,e_E1)
                end do
            end do
        end do

        ! 2D quadrature nodes mapped to the face of the physical tetrahedron
        do j=1,3
            do q=1,nq2
                pt(j,q)=0.0
                do l=1,4
                    pt(j,q) = pt(j,q) + temp(j,l)*nodtria2(l,q)
                enddo
            enddo
        enddo

        ! collect the coordinates of pt into different vectors
        x_p = pt(1,:)
        y_p = pt(2,:)
        z_p = pt(3,:)

        ! definition of the intervals intx, inty and intz
        intx = b_box1(1,:)
        inty = b_box1(2,:)
        intz = b_box1(3,:)

        ! loop on the basis functions
        do f = 1,Np

            ! initialization of valx, dvalx, valy, dvaly, valz and dvalz
            valx = 0.0
            dvalx = 0.0
            valy = 0.0
            dvaly = 0.0
            valz = 0.0
            dvalz = 0.0

            ! evaluation of Legendre Polynomials L_{blist(f,j)}, j=1,2,3, and their derivatives at x_p, y_p and z_p respectively
            call LegendreP(valx, x_p, blist(f,1), intx, nq2)
            call GradLegendreP(dvalx, x_p, blist(f,1), intx, nq2)
            call LegendreP(valy, y_p, blist(f,2), inty, nq2)
            call GradLegendreP(dvaly, y_p, blist(f,2), inty, nq2)
            call LegendreP(valz, z_p, blist(f,3), intz, nq2)
            call GradLegendreP(dvalz, z_p, blist(f,3), intz, nq2)

            ! evaluation of the f-th basis function, restricted in E1, and their partial derivatives in the q-th 2D quadrature node (mapped to the face of the physical tetrahedron)
            ! similar to what is done in the subroutine basis
            phi_b(f,:,1) = valx*valy*valz
            grad_b(1,f,:,1) = dvalx*valy*valz
            grad_b(2,f,:,1) = valx*dvaly*valz
            grad_b(3,f,:,1) = valx*valy*dvalz

        end do

        ! if E1 shares a face with E2, then repeat the whole procedure for E2
        ! equivalently, if it is false, then e_E1 is a boundary face
        if (E2 .ne. -1 .and. E2 .ne. -2) then

            intx = b_box2(1,:)
            inty = b_box2(2,:)
            intz = b_box2(3,:)

            do f=1,Np

                valx = 0.0
                dvalx = 0.0
                valy = 0.0
                dvaly = 0.0
                valz = 0.0
                dvalz = 0.0

                call LegendreP(valx, x_p, blist(f,1), intx, nq2)
                call GradLegendreP(dvalx, x_p, blist(f,1), intx, nq2)
                call LegendreP(valy, y_p, blist(f,2), inty, nq2)
                call GradLegendreP(dvaly, y_p, blist(f,2), inty, nq2)
                call LegendreP(valz, z_p, blist(f,3), intz, nq2)
                call GradLegendreP(dvalz, z_p, blist(f,3), intz, nq2)

                phi_b(f,:,2) = valx*valy*valz
                grad_b(1,f,:,2) = dvalx*valy*valz
                grad_b(2,f,:,2) = valx*dvaly*valz
                grad_b(3,f,:,2) = valx*valy*dvalz

            end do

        end if

    end subroutine basis_boundary

end module basis_function