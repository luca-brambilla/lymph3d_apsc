module matrix_free

    #include<petsc/finclude/petscksp.h>
        
    use petscksp
    use Poly_mesh
    use mpi
    use Poly_setup_mpi
    use problem_data_and_properties
    use basis_function
    use assemble_local
    use local_search
    use Poly_ref_mappings
    use Poly_data
    use vet_mat_operations

    ! for stiffness matrix that has different sizes depending on neighbors
    type :: RowArray
        real(kind=8), allocatable :: values(:,:,:,:,:) ! Dimensions you need
    end type RowArray

    contains

subroutine MAKE_MATRICES_FREE(PolyMesh, PolyData, petsc_num, global_dof, Np, K_loc, M_loc, rhs_loc)

    implicit none

    real(kind=8) :: present = 0.0
    real(kind=8) :: tmp = 0.0

    real(kind=8) :: dt2
    PetscScalar :: val(1), val1(1), val2(1), val3(1), val4(1)
    PetscInt :: irow(1), jcol(1), irow2(1), jcol2(1), irow3(1), jcol3(1)

    type(Mesh_Structure), intent(inout) :: PolyMesh
    type(Data_Structure), intent(in) :: PolyData
    integer(kind=4), intent(in) :: Np, global_dof
    integer(kind=4), dimension(global_dof), intent(in) :: petsc_num

    integer(kind=4) :: nq3, nq2, p
    real(kind=8) :: theta, alpha, c

    real(kind=8), dimension(4,4,4) :: node_maps
    real(kind=8), dimension(2,3,4) :: node_maps_inv

    real(kind=8), dimension(:,:), ALLOCATABLE :: nod3, nodtet3
    real(kind=8), dimension(:), ALLOCATABLE :: wei3, weitet3
    real(kind=8), dimension(:,:), ALLOCATABLE :: nod2, nodtria2
    real(kind=8), dimension (:), ALLOCATABLE :: wei2, weitria2
    integer(kind=4), dimension(:,:), ALLOCATABLE :: blist

    real(kind=8), dimension(:,:), ALLOCATABLE :: phi
    real(kind=8), dimension(:,:,:), ALLOCATABLE :: dphi
    real(kind=8), dimension(:,:,:), ALLOCATABLE :: phi_b
    real(kind=8), dimension(:,:,:,:), ALLOCATABLE :: grad_b

    real(kind=8), dimension(3,4) :: Fk
    real(kind=8) :: Jdet
    real(kind=8), dimension(3,3) :: Jinv
    real(kind=8), dimension(4) :: x, y, z

    real(kind=8) :: lambda, mu, rho ! density used for dynamics
    integer(kind=4) :: mat_id

    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, ipoly_loc, ipoly_glob, ipoly2_loc, ipoly2_glob
    integer(kind=4) :: n_tet_in_poly, iface_poly
    integer(kind=4) :: Npoly
    integer(kind=4) :: i, j, m, n

    integer(kind=4) :: e, E1, E2, sides
    integer(kind=4), dimension(4) :: face_flag
    real(kind=8), dimension(3) :: nn
    integer(kind=4) :: tag

    real(kind=8), dimension(3,Np) :: rhs_tet_loc
    real(kind=8), dimension(3,Np) :: rhs_face_bd_loc
    real(kind=8), dimension(3,3,Np,Np) :: V_loc
    real(kind=8), dimension(3,3,Np,Np) :: S_loc, I_loc, IN_loc, SN_loc

    ! each local has a matrix
    real(kind=8), dimension(PolyMesh%num_elem_loc,3,3,Np,Np), intent(out) :: M_loc
    ! real(kind=8), dimension(:,PolyMesh%num_elem_loc,3,3,Np,Np), allocatable, intent(out) :: K_loc
    real(kind=8), dimension(PolyMesh%num_elem_loc,3,Np), intent(out) :: rhs_loc
    type(RowArray), dimension(PolyMesh%num_elem_loc) :: K_loc

    PRINT *, '***** MATRIX FREE SETUP *****'

    dt2 = time_step * time_step
    ! set the properties of the method (see problem_data_and_properties.f90)
    call set_properties(alpha, theta, c)
    
    ! total degree of the basis functions
    p = PolyMesh%Elem_loc(1)%Degree;
    Npoly = PolyMesh%num_poly

    ! Computation of Gauss-Legendre quadrature nodes and weights over the reference square and cube
    ! (see basis_functions.f90)
    call quadrature(nod2, wei2, nod3, wei3, p, nq3, nq2)

    ! see Poly_ref_mappings.f90
    call tria2tetfaces_maps(node_maps, node_maps_inv)

    allocate(nodtet3(4,nq3))
    allocate(nodtria2(4,nq2))
    allocate(weitet3(nq3))
    allocate(weitria2(nq2))

    ! Maps to the reference tetrahedron and reference triangle (see Poly_ref_mappings.f90)
    call mapping_quadrature_3D(nod3, wei3, nq3, nodtet3, weitet3)
    call mapping_quadrature_2D(nod2, wei2, nq2, nodtria2, weitria2)
    
    ! list of the degrees of monomials of the Np basis functions up to order p (see basis_functions.f90)
    allocate(blist(Np,3))
    call basis_list(blist, p, Np)
        
    print *,'Assembling linear system...'
    
    allocate(phi(Np,nq3))
    allocate(dphi(3,Np,nq3))
    allocate(phi_b(Np,nq2,2))
    allocate(grad_b(3,Np,nq2,2))

    ! assign 0 to density for static case
    if (IsTime_dependent .eqv. .true.) then
        print *,'RHS with additional dynamic component'
        present = 1.0
    else
        print *,'RHS with only static component'
    endif

    ! initialize output once
    M_loc = 0.0
    K_loc = 0.0
    rhs_loc = 0.0

    ! loop on the tetrahedra
    do ie_loc = 1, PolyMesh%num_elem_loc

        E1 = ie_loc
        sides = PolyMesh%Elem_loc(E1)%num_faces
        
        ! initialize K_loc
        !!! CHECK IF IT POSSIBLE TO COUNT NEIGHBOURS
        allocate(K_loc(ie_loc)%values(sides, 3,3,Np,Np))
        K_loc(ie_loc)%values = 0.0
        ! initialization of V_loc
        V_loc = 0.0
        ! initialization of the rhs term on the volume rhs_tet_loc       
        rhs_tet_loc = 0.0


        mat_id = PolyMesh%Elem_loc(ie_loc)%mat_prop

        ! take correct density only for dynamic case 
        rho = PolyData%prop_mat(mat_id,1) ! DENSITY USED FOR DYNAMICS
        lambda = PolyData%prop_mat(mat_id,2)
        mu = PolyData%prop_mat(mat_id,3)

        ! computation of the coordinates of the tetrahedron
        do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert
        
            ! see MAKE_PARTITION_AND_MPI_FILES.f90
            call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)      
            
            x(ivert)=PolyMesh%coord_x(id_node)
            y(ivert)=PolyMesh%coord_y(id_node)
            z(ivert)=PolyMesh%coord_z(id_node)

        enddo

        ! computation of the reference map Fk, the inverse Jinv and the determinant Jdet of its jacobian (see Poly_ref_mappings.f90)
        call jacobians(x, y, z, Fk, Jinv, Jdet)
        
        ! find the polyhedron ipoly_glob that contains the tetrahedron ie_loc
        ie_glob = PolyMesh%elem_loc2glo(ie_loc)
        ipoly_glob = PolyMesh%elem_in_poly(ie_glob)

        ! see subroutine local_search in Poly_global.f90
        call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                                    PolyMesh%num_poly_loc, &
                                    ipoly_glob,ipoly_loc)

        ! evaluation of the basis functions and their partial derivatives at the 3D quadrature nodes for a given polyhedral element contained in b_box
        ! (see basis_functions.f90)
        call basis(phi, dphi, PolyMesh%Poly(ipoly_loc)%b_box, Np, blist, Fk, nodtet3, nq3)
        
        ! computation of the local stiffness matrix V_loc (see assemble_local.f90)
        ! already multiply by -dt^2, no need for actual V_loc
        call MAKE_STIFF_TET_LOC(Np, Jdet, weitet3, nq3, lambda, mu, dphi, V_loc)

        ! computation of the local mass matrix M_loc (see assemble_local.f90)
        !!! i only need true mass matrix
        call MAKE_MASS_LOC(Np, Jdet, weitet3, nq3, phi, rho, M_loc(ie_loc,:,:,:,:))

        ! computation of the local forcing vector, consider rhs with density only if dynamic case.

        ! call MAKE_RHS_TET(Np, Fk, Jdet, nodtet3, weitet3, nq3, lambda, mu, phi, rhs_tet_loc, rho)
        call MAKE_RHS_TET(Np, Fk, Jdet, nodtet3, weitet3, nq3, lambda, mu, phi, rhs_loc(ie_loc,:,:), rho*present)

        !!! CAN OPTIMIZE LOOPS COMMON INDEX i
        !!! NO NEED TO ADD, JUST ASSIGN?
        ! insert the values of V_loc in the entries of the local stiffness
        ! do i=1,3
        !     do j=1,3
        !         do m=1,Np
        !             do n=1,Np

        !                 tmp = V_loc(i,j,m,n)
        !                 if (tmp .ne. 0.0) then
        !                     K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) - dt2 * tmp
        !                 endif

        !             enddo
        !         enddo
        !     enddo
        ! enddo

        ! multiply the values of M_loc by density
        ! do i=1,3
        !     do m=1,Np
        !         do n=1,Np

        !             tmp = M_loc(ie_loc,i,i,m,n)
        !             if (tmp .ne. 0.0) then
        !                 M_loc(ie_loc,i,i,m,n) = rho * tmp
        !             endif
                    
        !         enddo
        !     enddo
        ! enddo

        ! insert the values of rhs_tet_loc in the entries of the global rhs vector
        ! do i=1,3
        !     do m=1,Np

        !         tmp = rhs_tet_loc(i,m)
        !         if (tmp .ne. 0.0) then
        !             rhs_loc(ie_loc,i,m) = rhs_loc(ie_loc,i,m) + tmp
        !         endif

        !     enddo
        ! enddo

        ! this allows to assemble the local matrix correctly into the global matrices
        !beg = (ipoly_glob-1)*Np + 1
        
        ! current element E+
        E1 = ie_loc

        ! begin loop on the faces of the tetrahedron E1
        do e=1,PolyMesh%Elem_loc(E1)%num_faces

            face_flag(e) = 0

            ! initialization of the face rhs term rhs_face_bd_loc
            rhs_face_bd_loc = 0.0
            ! initialization of the face matrices I_loc, S_loc, IN_loc and SN_loc
            I_loc = 0.0
            S_loc = 0.0
            IN_loc = 0.0
            SN_loc = 0.0

            ! find the neighbouring tetrahedron E2 sharing the face e with E1
            ! E2 is element E-
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(e,2)
            tag = PolyMesh%Elem_loc(E1)%neigh_el(e,5)

            ! if e is not a boundary face, then find the polyhedron in which E2 is contained
            if (E2 /= -1 .and. E2 /= -2) then

                ipoly2_glob = PolyMesh%elem_in_poly(E2)

                ! see subroutine local_search in Poly_global.f90
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                            PolyMesh%num_poly_loc, &
                            ipoly2_glob,ipoly2_loc)
            
            endif
            
            ! if e is not a boundary face, then check if E1 and E2 belong to the same polyhedron
            if (E2 /= -1 .and. E2 /= -2) then
                if (ipoly_glob == ipoly2_glob) then
                    face_flag(e) = 1
                endif
            end if

            ! if it is true, then E2 does not belong to the same polyhedron E1 belongs to
            ! or e is a boundary face
            if (face_flag(e) == 0) then

                nn = PolyMesh%Elem_loc(E1)%normal(e,:)

                ! If it is true, then the two polyhedra do not belong to the same processor 
                ! so we have to retrieve b_box of neighbouring element from neigh_bbox
                ! and hk of neighbouring element from neigh_hk
                ! Otherwise, the two polyhedra belong to the same processor 
                ! so that b_box and hk can be easily retrieved
                ! In both cases, compute the basis functions on the faces 
                ! and the local matrices on the faces
                if (ipoly2_loc==0) then

                    n_tet_in_poly=PolyMesh%Poly(ipoly_loc)%num_tet_in_poly

                    do j=1,n_tet_in_poly

                        if (PolyMesh%Poly(ipoly_loc)%tet_in_poly(j)==ie_glob) then
                            iface_poly=PolyMesh%Elem_loc(E1)%num_faces*(j-1)+e
                        endif

                    enddo

                    ! evaluation of the basis functions for every face of two neighbouring tetrahedra E1 and E2 at the 2D quadrature nodes
                    ! contained respectively in b_box1 and b_box2 (see basis_function.f90)
                    call basis_boundary(phi_b,grad_b,e, E2, PolyMesh%Poly(ipoly_loc)%b_box,&
                                        PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,:,:),blist, Np, Fk, node_maps, nodtria2, nq2)
                    
                    ! multiply lambda and mu by -dt^2 to get directly -dt^2 K
                    call MAKE_STIFF_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly), nn, &
                                        PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)
                    ! DO NOT MULTIPLY by -dt^2
                    call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly),&
                                    nn,PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,tag,rhs_face_bd_loc)

                else

                    ! check if e is not a boundary face
                    ! otherwise take a default "neighbouring" element (its information won't be read)
                    if (E2 /= -1 .and. E2 /= -2) then

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(ipoly2_loc)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        ! multiply lambda and mu by -dt^2 to get directly -dt^2 K                        
                        call MAKE_STIFF_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly2_loc)%hk,nn, &
                                                PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)

                        ! DO NOT MULTIPLY by -dt^2
                        call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly2_loc)%hk,nn, &
                                        PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,tag,rhs_face_bd_loc)

                    else

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(1)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)
                        
                        ! multiply lambda and mu by -dt^2 to get directly -dt^2 K  
                        call MAKE_STIFF_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(1)%hk, nn, &
                                                PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)
                        ! DO NOT MULTIPLY by -dt^2
                        call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(1)%hk,nn, &
                                        PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,tag,rhs_face_bd_loc)

                    endif

                endif

                ! if e is a boundary edge, then insert the values of rhs_face_bd_loc 
                ! in the entries of the rhs vector 
                if (E2 == -1 .or. E2 == -2) then

                    do i=1,3
                        do m=1,Np
                            tmp = rhs_face_bd_loc(i,m)
                            if (tmp .ne. 0.0) then
                                rhs_loc(ie_loc,i,m) = rhs_loc(ie_loc,i,m) + tmp
                            endif
                        enddo
                    enddo

                endif

            endif

            if(E2 /= -2) then

                ! insert the values of S_loc in the entries of the global stiffness and DG matrices
                do i=1,3
                    do j=1,3
                        do m=1,Np
                            do n=1,Np
                                tmp = S_loc(i,j,m,n)
                                if (tmp .ne. 0.0) then
                                    ! Add E+ contribution
                                    K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) + tmp
                                endif    
                            enddo
                        enddo
                    enddo
                enddo

                ! insert the values of I_loc, SN_loc and IN_loc in the entries of the global stiffness matrix
                ! and insert the values of SN_loc in the entries of the global DG matrix
                if (E2 == -1) then

                    do i=1,3
                        do j=1,3
                            do m=1,Np
                                do n=1,Np
                                    ! Add E+ contribution
                                    tmp = - I_loc(j,i,n,m)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) + tmp
                                    endif

                                    tmp = I_loc(i,j,m,n)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) + tmp
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo

                else

                    do i=1,3
                        do j=1,3
                            do m=1,Np
                                do n=1,Np
                                    ! Add E+ contribution
                                    tmp = - I_loc(j,i,n,m)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) + tmp
                                    endif

                                    tmp = I_loc(i,j,m,n)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(0,i,j,m,n) = K_loc(ie_loc)%values(0,i,j,m,n) + tmp
                                    endif
                                    
                                    ! Add E- contributions
                                    tmp = - IN_loc(j,i,n,m)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(e,i,j,m,n) = K_loc(ie_loc)%values(e,i,j,m,n) + tmp
                                    endif

                                    tmp = theta * IN_loc(i,j,m,n)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(e,i,j,m,n) = K_loc(ie_loc)%values(e,i,j,m,n) + tmp
                                    endif

                                    tmp = SN_loc(i,j,m,n)
                                    if (tmp .ne. 0.0) then
                                        K_loc(ie_loc)%values(e,i,j,m,n) = K_loc(ie_loc)%values(e,i,j,m,n) + tmp
                                    endif

                                enddo
                            enddo
                        enddo
                    enddo

                endif

            endif

        enddo
        
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! print *, 'Reaction coefficient c: ', c

    ! Matrix operation petsc_stiff + c*petsc_mass
    ! PetscCallA(MatAXPY(petsc_stiff, c, petsc_mass, DIFFERENT_NONZERO_PATTERN, mpi_ierr))

    ! Matrix operation mat_dg + c*petsc_mass
    ! PetscCallA(MatAXPY(mat_dg, c, petsc_mass, DIFFERENT_NONZERO_PATTERN, mpi_ierr))

    deallocate(phi)
    deallocate(dphi,phi_b)
    deallocate(grad_b)
    deallocate(nodtet3)
    deallocate(nodtria2)
    deallocate(weitet3)
    deallocate(weitria2)

    PRINT *, 'Done with assembling local matrices for matrix free'

end subroutine MAKE_MATRICES_FREE

!> @brief solver for matrix free considering only one element with time dependence
subroutine SOLVE_MATRIX_FREE(PolyMesh, Np, t, sides, K_loc_e, M_loc_e, rhs_loc_e, u0_loc, un_loc, u_loc)

    !TODO consider different Np
    !TODO compute inverse of mass matrix only once

    !> current time
    real(kind=8), intent(in) :: t
    !> degrees of freedom for local polynomial basis
    integer(kind=4), intent(in) :: Np
    !> number of contribution from neighbours E- (excluding boundary sides)
    integer(kind=4), intent(in) :: sides
    !> mass matrix for E+
    real(kind=8), dimension(3,3,Np,Np), intent(in) :: M_loc_e
    !> stiffness matrix with contributions from E+ and neighbors E-
    real(kind=8), dimension(:,3,3,Np,Np), intent(in) :: K_loc_e
    !> forcing term vector on E+
    real(kind=8), dimension(3,Np), intent(in) :: rhs_loc_e
    !> solution from two previous timesteps u^{(n-1)} for element E+
    real(kind=8), dimension(3,Np), intent(in) :: u0_loc
    !> solution from previous timestep for u^{(n)} with contributions from E+ and neighbors E-
    real(kind=8), dimension(sides,3,Np), intent(in) :: un_loc
    !> output solution for the element E+
    real(kind=8), dimension(3,Np), intent(out) :: u_loc

    !> partial result for contribution from stiffness and forcing term
    real(kind=8), dimension(3,Np) :: b

    real(kind=8) :: dt2

    !! HOW TO EXTRACT UN_LOC and U0_LOC ???

    ! u_loc = M_loc^-1 * ( - dt2 * K_loc*un_loc  + dt2*rhs_loc*time_function(t) ) +  2*un_loc - u0_loc

    ! initialize partial result vector
    b = 0.0

    ! compute effects of stiffness - matrix free in stiffness
    ! Contribution from K+ and neighbors K-
    do k=1,sides
        do i=1,3
            do j=1,3
                do m=1,Np
                    do n=1,Np
                        ! matrix vector multiplication K(k)*u_n(k)
                        b(i,m) = b(i,m) + K_loc_e(k,i,j,m,n)*un_loc(k,j,n)
                    end do
                end do
            end do
        end do
    end do

    ! compute effects of forcing term and solutions from previous time steps
    ! matrix free in mass
    do j=1,3
        do n=1,Np
            ! Add local forcing term once
            b(j,n) =   dt2 * (- b(j,n) + rhs_loc_e(j,n)*time_function(t) )
        end do

        ! matrix free solution of 3 different blocks of mass matrix
        !! FACTORIZATION AND SUBSTITUTION REPEATED FOR EACH TIME STEP
        u_loc(j,:) = linear_system(Np,M_loc_e(j,j,:,:),b(j,:),0)

        ! add solutions from previous time steps
        ! u_loc(0,:,:) has the modal coefficient corresponding to the current element E+ only
        do n=1,Np
            u_loc(j,n) = u_loc(j,n) + 2.0*un_loc(0,j,n) - u0_loc(j,n)
        end do

    end  do

end subroutine SOLVE_MATRIX_FREE

end module matrix_free