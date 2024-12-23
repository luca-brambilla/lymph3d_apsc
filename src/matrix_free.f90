!> @brief module for matrix free solution of DG FEM for elastodynamics
module matrix_free

#include<petsc/finclude/petscksp.h>

    use petscksp
    use mpi
    use Poly_setup_mpi
    use problem_data_and_properties
    use basis_function
    use assemble_element
    use local_search
    use Poly_ref_mappings
    use Poly_data
    use Poly_mesh
    use vet_mat_operations
    use mesh_partition_and_mpi_files

    use Poly_global
    use global_parameters
    use utilities

    implicit none

    ! for stiffness matrix that has different sizes depending on neighbors
    !! POINTER???????????
    type :: KRowArray
        real(kind=8), allocatable :: values(:,:,:,:,:)
    end type KRowArray

    type :: VecRowArray
        integer(kind=4), allocatable :: values(:)
    end type VecRowArray

    type :: PetscRowMat
        Mat :: values_x
        Mat :: values_y
        Mat :: values_z
    end type PetscRowMat

    type :: PetscRowVec
        Vec :: values_x
        Vec :: values_y
        Vec :: values_z
    end type PetscRowVec

    !> typedef for PETSc Matrix - to construct arrays of matrices
    !! maybe type(tMat), dimension(:,:)
    type :: PetscMatStruct
        Mat :: data
    end type

    !> typedef for PETSc Vector - to construct arrays of vectors
    !! maybe type(tVec), dimension(:)
    type :: PetscVecStruct
        Vec :: data
    end type

    contains


!> set up 2D array of PETSc matrices
subroutine SET_PETSC_MASS_MATRIX_FREE(ne_loc, Np, M)

    integer(kind=4), intent(in) :: ne_loc   !< number of rows
    integer(kind=4), intent(in) :: Np       !< PETSc square matrix size
    !> 2D array, store each block of the diagonal in a 1D array for each element
    type(PetscMatStruct), dimension(ne_loc,DIM), intent(out) :: M

    integer(kind=4) :: ie_loc, i

    do ie_loc=1,ne_loc
        do i=1,DIM
            PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, PETSC_NULL_SCALAR_ARRAY, M(ie_loc,i)%data, mpi_ierr))

            ! PetscCall(MatCreate(PETSC_COMM_SELF, M(ie_loc,i)%data, mpi_ierr))
            ! PetscCall(MatSetSizes(M(ie_loc,i)%data, Np, Np, Np, Np, mpi_ierr))
            ! PetscCall(MatSetFromOptions(M(ie_loc,i)%data, mpi_ierr))
            ! PetscCall(MatSetUp(M(ie_loc,i)%data, mpi_ierr))
        enddo
    enddo

end subroutine SET_PETSC_MASS_MATRIX_FREE

!> set up array of PETSc vectors
subroutine SET_PETSC_VECTOR_MATRIX_FREE(ne_loc, Np, V)

    integer(kind=4), intent(in) :: ne_loc   !< number of rows
    integer(kind=4), intent(in) :: Np       !< PETSc vector size
    !> array of PETSc vectors
    type(PetscVecStruct), dimension(ne_loc), intent(out) :: V

    integer(kind=4) :: ie_loc

    do ie_loc=1,ne_loc
        PetscCall(VecCreate(PETSC_COMM_SELF, V(ie_loc)%data, mpi_ierr))
        PetscCall(VecSetSizes(V(ie_loc)%data, Np, Np, mpi_ierr))
        PetscCall(VecSetFromOptions(V(ie_loc)%data, mpi_ierr))
    enddo

end subroutine SET_PETSC_VECTOR_MATRIX_FREE

!> @brief Compute the mass, stiffness, dg, modal matrices for each element.
!> The stiffness and dg matrices for the element E+ are rectangular and contain
!> the contributions also from neighboring elements E-.
!> mass matrix is directly in PETSc for later to solve linear systems.
subroutine MAKE_MATRICES_FREE(PolyMesh, PolyData, num_elem_loc, Np, K_loc, A_dg_loc, M_loc, M_modal_loc, max_faces)

    !TODO variable number of sides, do not count boundaries
    !TODO polytopal elements, face contribution to same matrices
    !TODO deallocate allocations for matrices

    implicit none

    type(Mesh_Structure), intent(inout) :: PolyMesh !< Mesh
    type(Data_Structure), intent(in) :: PolyData    !< Data
    integer(kind=4), intent(in) :: Np               !< number of degrees of freedom of each element
    integer(kind=4), intent(in) :: num_elem_loc       !< Number of global dofs

    !type(PetscMatStruct), dimension(PolyMesh%num_elem_loc, DIM), intent(inout) :: massa
    !type(PetscMatStruct), dimension(PolyMesh%num_elem_loc, DIM), intent(inout) :: massa_modale

    integer(kind=4) :: nq3, nq2, p
    real(kind=8) :: theta, alpha, c

    real(kind=8), dimension(4,4,4) :: node_maps
    real(kind=8), dimension(2,3,4) :: node_maps_inv

    real(kind=8), dimension(:,:), allocatable :: nod3, nodtet3
    real(kind=8), dimension(:), allocatable :: wei3, weitet3
    real(kind=8), dimension(:,:), allocatable :: nod2, nodtria2
    real(kind=8), dimension (:), allocatable :: wei2, weitria2
    integer(kind=4), dimension(:,:), allocatable :: blist

    real(kind=8), dimension(:,:), allocatable :: phi
    real(kind=8), dimension(:,:,:), allocatable :: dphi
    real(kind=8), dimension(:,:,:), allocatable :: phi_b
    real(kind=8), dimension(:,:,:,:), allocatable :: grad_b

    real(kind=8), dimension(DIM,DIM+1) :: Fk    ! tranformation
    real(kind=8) :: Jdet                        ! determinant of Jacobian
    real(kind=8), dimension(3,3) :: Jinv        ! inverse of Jacobian
    real(kind=8), dimension(4) :: x, y, z

    real(kind=8) :: lambda, mu, rho ! density used for dynamics
    integer(kind=4) :: mat_id

    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, ipoly_loc, ipoly_glob, ipoly2_loc, ipoly2_glob
    integer(kind=4) :: n_tet_in_poly, iface_poly
    integer(kind=4) :: Npoly
    integer(kind=4) :: i, j, m, n

    integer(kind=4) :: iface, E1, E2!, sides
    integer(kind=4), dimension(4) :: face_flag !! WHY VECTOR????
    real(kind=8), dimension(DIM) :: nn
    integer(kind=4) :: space_fun_tag
    !integer(kind=4) :: n_neigh      !< number of neighbor internal faces
    !integer(kind=4) :: neigh_count  !< counter for inserting contributions in K_loc

    ! integer(kind=4), intent(out) :: internal_neigh(PolyMesh%num_elem_loc,:) !< internal neighbors connectivity - similar to el_neigh from PolyMesh
    integer(kind=4), intent(out) :: max_faces !< maximum number of polygon faces

    real(kind=8), dimension(DIM, DIM, Np, Np) :: V_loc
    real(kind=8), dimension(DIM, DIM, Np, Np) :: S_E1, I_E1, S_E2, I_E2, IT_E2

    ! each local has a matrix
    ! local square
    real(kind=8), dimension(DIM, DIM, Np, Np) :: mass_loc

    real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: M_loc
    real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: M_modal_loc

    ! process local matrices, rectangular for each element
    real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: K_loc
    real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: A_dg_loc

    integer(kind=4) :: row, col

    real(kind=8) :: t1, t2

    ! set the properties of the method (see problem_data_and_properties.f90)
    call set_properties(alpha, theta, c)

    ! total degree of the basis functions
    p = PolyMesh%Elem_loc(1)%Degree
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

    print *,'Assembling local matrices...'

    allocate(phi(Np,nq3))
    allocate(dphi(3,Np,nq3))
    allocate(phi_b(Np,nq2,2))
    allocate(grad_b(3,Np,nq2,2))

    ! loop over elements to find the maximum number of polygonal faces
    !! IS ALL SPLIT IN POLYGONS OR TETRAHEDRA???
    max_faces = PolyMesh%Elem_loc(1)%num_faces
    ! do ie_loc = 2, PolyMesh%num_elem_loc
    !     sides = PolyMesh%Elem_loc(ie_loc)%num_faces
    !     if (sides > max_faces) max_faces = sides
    ! end do

    ! allocate(internal_neigh(PolyMesh%num_elem_loc, max_faces))
    ! internal_neigh = 0

    ! initialize output once
    K_loc = 0.0d0
    A_dg_loc = 0.0d0
    M_loc = 0.0d0
    M_modal_loc = 0.0d0

    ! loop on the tetrahedra
    elem_loop: do ie_loc = 1, PolyMesh%num_elem_loc

        ! count neighbor element contribution only to allocate K_loc
        ! current element E+
        E1 = ie_loc
        ! sides = PolyMesh%Elem_loc(E1)%num_faces ! could vary for each element
        !n_neigh = 0     ! could vary for each element
        !neigh_count = 1

        ! begin loop on the faces of the tetrahedron E1
        ! do iface=1,sides
        !     E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)
        !     ! space_fun_tag = PolyMesh%Elem_loc(E1)%neigh_el(iface,5)
        !     ! if not boundary, then internal face
        !     if (E2 /= -1 .and. E2 /= -2) n_neigh = n_neigh + 1
        ! end do

        ! initialize K_loc
        ! allocate(K_loc(ie_loc)%values(n_neigh, 3,3,Np,Np))
        ! K_loc(ie_loc)%values = 0.0

        ! allocate(A_dg_loc(ie_loc)%values(n_neigh, 3,3,Np,Np))
        ! A_dg_loc(ie_loc)%values = 0.0

        ! initialization of V_loc
        V_loc = 0.0d0
        mass_loc = 0.0d0

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

        ! computation of the local mass matrix M_loc (see assemble_element.f90)
        !!! i only need true mass matrix
        !call MAKE_MASS_VOLUME(Np, Jdet, weitet3, nq3, phi, rho, M_loc(ie_loc,:,:,:,:))

        t1 = MPI_WTIME()
        ! for petsc
        call MAKE_MASS_VOLUME(Np, Jdet, weitet3, nq3, phi, 1.0d0, mass_loc)

        ! PETSc populate matrix
        ! insert the values of M_loc in the 3 blocks of the mass matrix
        do i=1,DIM

            M_modal_loc(ie_loc,i,:,:) = mass_loc(i,i,:,:)
            M_loc(ie_loc,i,:,:) = mass_loc(i,i,:,:)*rho

        enddo
        t2 = MPI_WTIME()
        tp_setup_M = tp_setup_M + t2 - t1

        t1 = MPI_WTIME()

        ! computation of the local stiffness matrix V_loc (see assemble_element.f90)
        call MAKE_STIFFNESS_VOLUME(Np, Jdet, weitet3, nq3, lambda, mu, dphi, V_loc)

        ! current element E+
        E1 = ie_loc

        ! insert values of V_loc into stiffness matrix
        do i=1,DIM
            row = (i-1)*Np
            do j = 1,DIM
                col = (j-1)*Np
                do m = 1,Np
                    do n = 1,Np
                        K_loc(ie_loc, 1, row+m, col+n) = V_loc(i,j,m,n)
                        A_dg_loc(ie_loc, 1, row+m, col+n) = V_loc(i,j,m,n)
                    end do
                end do
            end do
        end do

        ! stiffness and rhs
        ! begin loop on the faces of the tetrahedron E1
        face_loop: do iface=1,PolyMesh%Elem_loc(E1)%num_faces

            face_flag(iface) = 0

            ! initialization of the face matrices I_E1, S_E1, I_E2 and S_E2
            I_E1 = 0.0d0
            S_E1 = 0.0d0
            I_E2 = 0.0d0
            S_E2 = 0.0d0
            IT_E2 = 0.0d0

            ! find the neighbouring tetrahedron E2 sharing the face iface with E1
            ! E2 is element E-
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)
            space_fun_tag = PolyMesh%Elem_loc(E1)%neigh_el(iface,5)

            ! if iface is not a boundary face
            if (E2 /= -1 .and. E2 /= -2) then

                ! find the polyhedron in which E2 is contained
                ipoly2_glob = PolyMesh%elem_in_poly(E2)

                ! see subroutine local_search in Poly_global.f90
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                            PolyMesh%num_poly_loc, &
                            ipoly2_glob,ipoly2_loc)

                ! check if E1 and E2 belong to the same polyhedron
                if (ipoly_glob == ipoly2_glob) then
                    face_flag(iface) = 1
                endif

            endif

            ! if it is true, then E2 does not belong to the same polyhedron E1 belongs to
            ! or iface is a boundary face
            if (face_flag(iface) == 0) then

                nn = PolyMesh%Elem_loc(E1)%normal(iface,:)

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
                            iface_poly=PolyMesh%Elem_loc(E1)%num_faces*(j-1)+iface
                        endif

                    enddo

                    ! evaluation of the basis functions for every face of two neighbouring tetrahedra E1 and E2 at the 2D quadrature nodes
                    ! contained respectively in b_box1 and b_box2 (see basis_function.f90)
                    call basis_boundary(phi_b,grad_b,iface, E2, PolyMesh%Poly(ipoly_loc)%b_box,&
                                        PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,:,:),blist, Np, Fk, node_maps, nodtria2, nq2)

                    call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly), nn, &
                                        PolyMesh%Elem_loc(E1)%area(iface),weitria2,nq2,lambda,mu,phi_b,grad_b,S_E1,I_E1,S_E2,I_E2,IT_E2)

                ! polyhedra on the same process
                else

                    ! check if iface is not a boundary face
                    ! otherwise take a default "neighbouring" element (its information won't be read)
                    if (E2 /= -1 .and. E2 /= -2) then

                        call basis_boundary(phi_b,grad_b,iface,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(ipoly2_loc)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly2_loc)%hk,nn, &
                                                PolyMesh%Elem_loc(E1)%area(iface),weitria2,nq2,lambda,mu,phi_b,grad_b,S_E1,I_E1,S_E2,I_E2,IT_E2)

                    else

                        call basis_boundary(phi_b,grad_b,iface,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(1)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(1)%hk, nn, &
                                                PolyMesh%Elem_loc(E1)%area(iface),weitria2,nq2,lambda,mu,phi_b,grad_b,S_E1,I_E1,S_E2,I_E2,IT_E2)

                    endif

                endif

            endif

            !! INDEX i,j for I matrix????? WHERE THETA
            ! Insert in stiffness matrix
            ! check face again and add boundary contributions
            ! If not Neumann boundary
            if(E2 /= -2) then

                ! Add E+ contribution S_E1 to stiffness and DG
                do i=1,DIM
                    do j=1,DIM
                        do m=1,Np
                            row = (i-1)*Np + m
                            do n=1,Np
                                col = (j-1)*Np + n
                                K_loc(ie_loc,1,row,col) = K_loc(ie_loc,1,row,col) + S_E1(i,j,m,n)
                                A_dg_loc(ie_loc,1,row,col) = A_dg_loc(ie_loc,1,row,col) + S_E1(i,j,m,n)
                            enddo
                        enddo
                    enddo
                enddo

                ! internal face contribution, not Dirichlet boundary
                ! Add neighbor face contribution E- (S_E2 to stiffness and DG - I_E2 to stiffness)
                ! leave matrix K_loc portion 0.0 if not Dirichlet
                ! E- contribution in position iface+1 since position 1 is for E+

                ! if Dirichlet
                if (E2 == -1) then

                    do i=1,DIM
                        do j=1,DIM
                            do m=1,Np
                                row = (i-1)*Np + m
                                do n=1,Np
                                    col = (j-1)*Np + n

                                    ! E+ on E+ -> theta*I - I^T
                                    K_loc(ie_loc,1,row,col) = K_loc(ie_loc,1,row,col) + theta*I_E1(i,j,m,n) - I_E1(j,i,n,m)
                                enddo
                            enddo
                        enddo
                    enddo
                    !internal_neigh(ie_loc,neigh_count) = E2
                    ! neigh_count = neigh_count + 1

                    ! print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, '- column: ', iface+1, ' - DIRICHLET'

                ! if internal face
                else

                    do i=1,DIM
                        do j=1,DIM
                            do m=1,Np
                                row = (i-1)*Np + m
                                do n=1,Np
                                    col = (j-1)*Np + n

                                    ! E+ on E+ -> theta*I - I^T
                                    K_loc(ie_loc,1,row,col) = K_loc(ie_loc,1,row,col) + theta*I_E1(i,j,m,n) - I_E1(j,i,n,m)

                                    ! E+ on E- -> theta*I + S
                                    K_loc(ie_loc,iface+1,row,col) = K_loc(ie_loc,iface+1,row,col) + theta*I_E2(i,j,m,n) + S_E2(i,j,m,n)

                                    A_dg_loc(ie_loc,iface+1,row,col) = A_dg_loc(ie_loc,iface+1,row,col) + S_E2(i,j,m,n)

                                    ! E- on E+ -> -I^T
                                    K_loc(ie_loc,iface+1,row,col) = K_loc(ie_loc,iface+1,row,col) - IT_E2(i,j,m,n)

                                enddo
                            enddo
                        enddo
                    enddo

                endif

            endif
            !neigh_count = neigh_count + 1

        enddo face_loop

        t2 = MPI_WTIME()
        tp_setup_K = tp_setup_K + t2 - t1

    enddo elem_loop

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! print *, 'Reaction coefficient c: ', c

    ! Matrix operation petsc_stiff + c*petsc_mass
    ! PetscCall(MatAXPY(petsc_stiff, c, petsc_mass, DIFFERENT_NONZERO_PATTERN, mpi_ierr))

    ! Matrix operation mat_dg + c*petsc_mass
    ! PetscCall(MatAXPY(mat_dg, c, petsc_mass, DIFFERENT_NONZERO_PATTERN, mpi_ierr))

    deallocate(phi)
    deallocate(dphi,phi_b)
    deallocate(grad_b)
    deallocate(nodtet3)
    deallocate(nodtria2)
    deallocate(weitet3)
    deallocate(weitria2)
    deallocate(blist)

    PRINT *, 'Done with assembling local matrices'

end subroutine MAKE_MATRICES_FREE

!> @brief Compute the RHS vector for each element.
subroutine MAKE_RHS_FREE(PolyMesh, PolyData, num_elem_loc, Np, rhs_loc, f_forcing, f_dirichlet, f_neumann)

    implicit none

    ! PASS FUNCTION AS ARGUMENT
    interface
        function f_forcing(lambda, mu, point, rho) result(res)
            use global_parameters, only: DIM
            real(kind=8) :: rho, lambda, mu
            real(kind=8), dimension(DIM) :: point, res
        end function f_forcing

        function f_dirichlet(point, space_fun_tag) result(res)
            use global_parameters, only: DIM
            real(kind=8), dimension(DIM) :: point, res
            integer(kind=4) :: space_fun_tag
        end function f_dirichlet

        function f_neumann(lambda, mu, normal, point, space_fun_tag) result(res)
            use global_parameters, only: DIM
            real(kind=8) :: lambda, mu
            real(kind=8), dimension(DIM) :: point, res, normal
            integer(kind=4) :: space_fun_tag
        end function f_neumann
    end interface


    type(Mesh_Structure), intent(inout) :: PolyMesh !< Mesh
    type(Data_Structure), intent(in) :: PolyData    !< Data
    integer(kind=4), intent(in) :: Np               !< number of degrees of freedom of each element
    integer(kind=4), intent(in) :: num_elem_loc       !< Number of global dofs

    real(kind=8) :: present = 0.0d0
    real(kind=8) :: tmp = 0.0d0

    integer(kind=4) :: nq3, nq2, p
    real(kind=8) :: theta, alpha, c

    real(kind=8), dimension(4,4,4) :: node_maps
    real(kind=8), dimension(2,3,4) :: node_maps_inv

    real(kind=8), dimension(:,:), allocatable :: nod3, nodtet3
    real(kind=8), dimension(:), allocatable :: wei3, weitet3
    real(kind=8), dimension(:,:), allocatable :: nod2, nodtria2
    real(kind=8), dimension (:), allocatable :: wei2, weitria2
    integer(kind=4), dimension(:,:), allocatable :: blist

    real(kind=8), dimension(:,:), allocatable :: phi
    real(kind=8), dimension(:,:,:), allocatable :: dphi
    real(kind=8), dimension(:,:,:), allocatable :: phi_b
    real(kind=8), dimension(:,:,:,:), allocatable :: grad_b

    real(kind=8), dimension(DIM,DIM+1) :: Fk
    real(kind=8) :: Jdet
    real(kind=8), dimension(3,3) :: Jinv
    real(kind=8), dimension(4) :: x, y, z

    real(kind=8) :: lambda, mu, rho ! density used for dynamics
    integer(kind=4) :: mat_id

    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, ipoly_loc, ipoly_glob, ipoly2_loc, ipoly2_glob
    integer(kind=4) :: n_tet_in_poly, iface_poly
    integer(kind=4) :: Npoly
    integer(kind=4) :: i, j, m

    integer(kind=4) :: iface, E1, E2!, sides
    integer(kind=4), dimension(4) :: face_flag !! WHY VECTOR???
    real(kind=8), dimension(DIM) :: nn
    integer(kind=4) :: space_fun_tag
    !integer(kind=4) :: n_neigh      !< number of neighbor internal faces
    !integer(kind=4) :: neigh_count  !< counter for inserting contributions in K_loc

    real(kind=8), dimension(DIM, Np) :: rhs_tet_loc
    real(kind=8), dimension(DIM, Np) :: rhs_face_bd_loc

    ! local vector
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM, Np) :: rhs_loc_tmp
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(inout) :: rhs_loc

    integer(kind=4) :: row
    real(kind=8) :: t1, t2

    call set_properties(alpha, theta, c)

    ! total degree of the basis functions
    p = PolyMesh%Elem_loc(1)%Degree
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

    print *,'Assembling RHS...'

    allocate(phi(Np,nq3))
    allocate(dphi(3,Np,nq3))
    allocate(phi_b(Np,nq2,2))
    allocate(grad_b(3,Np,nq2,2))

    ! assign 0 to density for static case
    if (IsTime_dependent .eqv. .true.) then
        print *,'RHS with additional dynamic component'
        present = 1.0d0
    else
        print *,'RHS with only static component'
    endif

    ! initialize output once
    rhs_loc = 0.0d0

    ! loop on the tetrahedra
    elem_loop: do ie_loc = 1, PolyMesh%num_elem_loc

        ! initialization of the rhs term on the volume rhs_tet_loc
        rhs_tet_loc = 0.0d0

        ! current element E+
        E1 = ie_loc

        mat_id = PolyMesh%Elem_loc(ie_loc)%mat_prop
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

        t1 = MPI_WTIME()

        call MAKE_RHS_VOLUME(Np, Fk, Jdet, nodtet3, weitet3, nq3, lambda, mu, phi, rhs_loc_tmp(ie_loc,:,:), rho*present, f_forcing)

        ! copy to output in correct format
        do i=1,DIM
            row = (i-1)*Np
            rhs_loc(ie_loc, row+1:row+Np) = rhs_loc_tmp(ie_loc,i,:)
        enddo

        ! current element E+
        E1 = ie_loc
        face_loop: do iface=1,PolyMesh%Elem_loc(E1)%num_faces
            face_flag(iface) = 0

            ! initialization of the face rhs term rhs_face_bd_loc
            rhs_face_bd_loc = 0.0d0

            ! find the neighbouring tetrahedron E2 sharing the face iface with E1
            ! E2 is element E-
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)
            space_fun_tag = PolyMesh%Elem_loc(E1)%neigh_el(iface,5)

            ! if iface is not a boundary face
            if (E2 /= -1 .and. E2 /= -2) then

                ! find the polyhedron in which E2 is contained
                ipoly2_glob = PolyMesh%elem_in_poly(E2)

                ! see subroutine local_search in Poly_global.f90
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                            PolyMesh%num_poly_loc, &
                            ipoly2_glob,ipoly2_loc)

                ! check if E1 and E2 belong to the same polyhedron
                if (ipoly_glob == ipoly2_glob) then
                    face_flag(iface) = 1
                endif

            endif

            ! if it is true, then E2 does not belong to the same polyhedron E1 belongs to
            ! or iface is a boundary face
            if (face_flag(iface) == 0) then

                nn = PolyMesh%Elem_loc(E1)%normal(iface,:)

                ! If it is true, then the two polyhedra do not belong to the same process
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
                            iface_poly=PolyMesh%Elem_loc(E1)%num_faces*(j-1)+iface
                        endif

                    enddo

                    ! evaluation of the basis functions for every face of two neighbouring tetrahedra E1 and E2 at the 2D quadrature nodes
                    ! contained respectively in b_box1 and b_box2 (see basis_function.f90)
                    call basis_boundary(phi_b,grad_b,iface, E2, PolyMesh%Poly(ipoly_loc)%b_box,&
                                        PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,:,:),blist, Np, Fk, node_maps, nodtria2, nq2)

                    call MAKE_RHS_FACE(theta,alpha,p,Np,iface,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly),&
                                    nn,PolyMesh%Elem_loc(E1)%area(iface),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc, f_dirichlet, f_neumann)

                else

                    ! check if iface is not a boundary face
                    ! otherwise take a default "neighbouring" element (its information won't be read)
                    if (E2 /= -1 .and. E2 /= -2) then

                        call basis_boundary(phi_b,grad_b,iface,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(ipoly2_loc)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_RHS_FACE(theta,alpha,p,Np,iface,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly2_loc)%hk,nn, &
                                        PolyMesh%Elem_loc(E1)%area(iface),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc, f_dirichlet, f_neumann)

                    else

                        call basis_boundary(phi_b,grad_b,iface,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(1)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_RHS_FACE(theta,alpha,p,Np,iface,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(1)%hk,nn, &
                                        PolyMesh%Elem_loc(E1)%area(iface),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc, f_dirichlet, f_neumann)

                    endif

                endif

                ! if iface is a boundary edge, then insert the values of rhs_face_bd_loc
                ! in the entries of the rhs vector
                if (E2 == -1 .or. E2 == -2) then
                    do i=1,DIM
                        do m=1,Np
                            row = (i-1)*Np + m
                            tmp = rhs_face_bd_loc(i,m)
                            if (tmp .ne. 0.0d0) then
                                rhs_loc(ie_loc,row) = rhs_loc(ie_loc,row) + tmp
                            endif
                        enddo
                    enddo
                endif

            endif

        enddo face_loop

        t2 = MPI_WTIME()
        tp_setup_RHS = tp_setup_RHS + t2 - t1

    enddo elem_loop

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    deallocate(phi)
    deallocate(dphi,phi_b)
    deallocate(grad_b)
    deallocate(nodtet3)
    deallocate(nodtria2)
    deallocate(weitet3)
    deallocate(weitria2)
    deallocate(blist)

    print *, 'Done with assembling local RHS for matrix free'

end subroutine MAKE_RHS_FREE

!> Matrix-free context, start from nodal solution and get modal solution coefficients.
!> Data is scattered already, performed in series by each processor, local PETSc definition of vectors and matrix to use solver.
subroutine COMPUTE_MODAL_COEFFICIENTS_FREE(PolyMesh, Np, R_M_modal_loc, f_analytic, modal_coeff)

    implicit none

    ! PASS FUNCTION AS ARGUMENT
    interface
        function f_analytic(point) result(res)
            use global_parameters, only: DIM
            real(kind=8), dimension(DIM) :: point, res
        end function f_analytic
    end interface

    !type(PetscMatStruct), dimension(:,:), intent(in):: massa_modale
    real(kind=8), dimension(:,:,:,:), intent(in) :: R_M_modal_loc

    type(Mesh_Structure) :: PolyMesh
    integer(kind=4) :: p, Np
    integer(kind=4) :: nq3, nq2, Npoly
    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, i, m
    integer(kind=4) :: ipoly_loc, ipoly_glob
    integer(kind=4), dimension(:,:), allocatable :: blist
    real(kind=8) :: Jdet
    real(kind=8), dimension(4) :: x, y, z
    real(kind=8), dimension(DIM,DIM+1) :: Fk
    real(kind=8), dimension(3,3) :: Jinv
    real(kind=8), dimension(:,:), allocatable :: nod2, nod3, nodtet3
    real(kind=8), dimension(:), allocatable :: wei2, wei3, weitet3
    real(kind=8), dimension(DIM) :: points
    real(kind=8), dimension(:,:), allocatable :: phi
    real(kind=8), dimension(:,:,:), allocatable :: dphi
    real(kind=8), dimension(:,:), allocatable, intent(out) :: modal_coeff
    real(kind=8), dimension(:,:,:), allocatable :: uex_integral
    real(kind=8), dimension(DIM) :: eval

    integer(kind=4) :: q, ii, jj
    integer(kind=4) :: row

    real(kind=8), dimension(Np) :: v_ptr
    real(kind=8), dimension(Np,Np) :: R, RT

    ! logical, intent(in) :: IsTime_dependent
    ! real(kind=8), intent(in), optional :: time

    p = PolyMesh%Elem_loc(1)%Degree !< degree of the polynomial
    Npoly = PolyMesh%num_poly       !< number of polyhedra

    ! Computation of Gauss-Legendre quadrature nodes and weights over the reference square and cube
    ! (see basis_functions.f90)
    call quadrature(nod2, wei2, nod3, wei3, p, nq3, nq2)

    allocate(nodtet3(4,nq3))
    allocate(weitet3(nq3))

    ! Maps to the reference tetrahedron (see Poly_ref_mappings.f90)
    call mapping_quadrature_3D(nod3, wei3, nq3, nodtet3, weitet3)

    ! list of the degrees of monomials of the Np basis functions up to order p (see basis_functions.f90)
    allocate(blist(Np,3))
    call basis_list(blist, p, Np)

    allocate(phi(Np,nq3))
    allocate(dphi(3,Np,nq3))

    ! solution initializations
    allocate(uex_integral(PolyMesh%num_poly_loc, DIM, Np))
    allocate(modal_coeff(PolyMesh%num_poly_loc, DIM*Np))
    uex_integral = 0.0d0
    modal_coeff = 0.0d0

    ! loop on the tetrahedra
    elem_loop: do ie_loc = 1, PolyMesh%num_elem_loc

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

        ! loop on the 3D quadrature nodes
        do q = 1,nq3
            do m=1,Np

                do ii=1,DIM
                    points(ii)=0.0d0
                    do jj=1,4
                        points(ii)=points(ii)+Fk(ii,jj)*nodtet3(jj,q)
                    enddo
                enddo

                eval = f_analytic(points)

                do i=1,DIM
                    uex_integral(ie_loc,i,m) = uex_integral(ie_loc,i,m) + abs(Jdet)*weitet3(q)*eval(i)*phi(m,q)
                enddo

            end do
        end do

        do i=1,DIM

            row = (i-1)*Np

            ! copy vector
            v_ptr = uex_integral(ie_loc,i,:)

            ! solve linear system matrix-free on block (i,i)

            R = R_M_modal_loc(ie_loc,i,:,:)
            RT = transpose(R)
            modal_coeff(ie_loc,row+1:row+Np) = solve_LU(RT, R, v_ptr, Np)

        enddo

    enddo elem_loop

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if (mpi_id==0) print *, 'Done with computing modal coefficients'

    deallocate(phi)
    deallocate(dphi)

    deallocate(nodtet3)
    deallocate(weitet3)
    deallocate(blist)

    deallocate(uex_integral)

end subroutine COMPUTE_MODAL_COEFFICIENTS_FREE

!> allocate global solution and prepare size and displacement vectors for MPI
subroutine PREPROCESS_SOLUTION_MATRIX_FREE(PolyMesh, local_dof, nnod_num, gathered_sizes, displacements, u)


    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    integer(kind=4), intent(in) :: local_dof        !< number of local dof
    ! integer(kind=4), intent(in) :: mpi_id
    integer(kind=4), dimension(:), allocatable, intent(out) :: nnod_num
    integer(kind=4), dimension(:), allocatable, intent(out) :: gathered_sizes !< MPI gather size per dimension
    integer(kind=4), dimension(:), allocatable, intent(out) :: displacements !< MPI gather displacements per dimension
    real(kind=8), dimension(:,:), allocatable, intent(out) :: u !< solution (Np, 3*num_poly)
    integer(kind=4) :: i, Np

    Np = PolyMesh%Elem_loc(1)%NDof_elem

    ! Allocate solution for post-processing
    allocate(u(Np, DIM * PolyMesh%num_poly))

    ! STORE LOCAL NUMERATION TO RECONSTRUCT THE SOLUTION
    allocate(nnod_num(local_dof))
    call CREATE_LOCAL_NODE_NUM(nnod_num, local_dof)

    if(mpi_np > 1) then

        allocate(gathered_sizes(mpi_np))

        call MPI_AllGather(Np*PolyMesh%num_poly_loc, 1, MPI_INTEGER, gathered_sizes, 1, &
                    MPI_INTEGER, MPI_COMM_WORLD, ierr)

        allocate(displacements(mpi_np))
        displacements(1) = 0
        do i = 2, mpi_np
                displacements(i) = displacements(i - 1) + gathered_sizes(i - 1)
        end do

        ! print *, 'displacements', displacements
        ! print *, 'gathered_sizes', gathered_sizes
    endif

end subroutine PREPROCESS_SOLUTION_MATRIX_FREE

!> gather global solution from local solution to compute errors and save output
subroutine POST_PROCESS_MATRIX_FREE(PolyMesh, u_loc, u_glo, gathered_sizes, displacements)

    use Poly_mesh
    use Poly_setup_MPI
    use global_parameters

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh            !< mesh
    real(kind=8), dimension(:,:), intent(in) :: u_loc       !< local solution
    real(kind=8), dimension(:,:), intent(inout) :: u_glo    !< global solution

    real(kind=8), dimension(PolyMesh%Elem_loc(1)%NDof_elem, PolyMesh%num_poly_loc) :: u_loc_t
    real(kind=8), dimension(PolyMesh%Elem_loc(1)%NDof_elem, PolyMesh%num_poly) :: u_tmp_dim
    integer(kind=4), dimension(:), intent(in) :: gathered_sizes, displacements
    integer(kind=4) :: Np, i, jcol_glo, jcol_tmp, Npoly, Npoly_loc

    call LYMPH3D_BARRIER

    ! assuming same degree everywhere
    Np = PolyMesh%Elem_loc(1)%NDof_elem ! ndof local per dimension (for 3D we need 3*Np)
    Npoly = PolyMesh%num_poly
    Npoly_loc = PolyMesh%num_poly_loc

    ! SCATTER PETSC SOLUTION AND STORE IN A FORTRAN ARRAY
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    if(mpi_np == 1) then

        ! transpose each dimension block
        do i=1,DIM
            jcol_glo = (i-1)*Npoly
            jcol_tmp = (i-1)*Np
            u_glo(:,jcol_glo+1:jcol_glo+Npoly) = transpose(u_loc(:,jcol_tmp+1:jcol_tmp+Np))
        enddo

    else

        do i=1,DIM
            jcol_tmp = (i-1)*Np
            jcol_glo = (i-1)*Npoly

            u_loc_t = transpose(u_loc(:,jcol_tmp+1:jcol_tmp+Np))

            call MPI_ALLGATHERV(u_loc_t, gathered_sizes(mpi_id+1), MPI_DOUBLE_PRECISION, u_tmp_dim, gathered_sizes, displacements, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_ierr)

            u_glo(:,jcol_glo+1:jcol_glo+Npoly) = u_tmp_dim

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
        enddo

    endif

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! print *,'Done with the solution'

end subroutine POST_PROCESS_MATRIX_FREE

!> @brief solver for matrix free considering only one element with time dependence
subroutine FIRST_TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, R_M_loc, rhs_stat_loc, rhs_dyn_loc, u0_loc, v0_loc, usol_loc, u0_mpi)

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh

    !> current time
    real(kind=8), intent(in) :: t
    !> degrees of freedom for local polynomial basis
    integer(kind=4), intent(in) :: Np
    !> mass matrix for E+
    real(kind=8), dimension(:,:,:,:), intent(in) :: R_M_loc

    !> stiffness matrix with contributions only from element E+ and neighbors E-
    real(kind=8), dimension(:,:,:,:), intent(in) :: K_loc
    !> forcing term vector on E+ - static component
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: rhs_stat_loc
    !> forcing term vector on E+ - dynamic component
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: rhs_dyn_loc

    !> IC displacement
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: u0_loc
    !> IC velocity
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: v0_loc
    !> solution of first timestep for u^{(1)}
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(out) :: usol_loc
    !> IC displacement from other processes
    real(kind=8), dimension(PolyMesh%num_elem_inter_vec(mpi_id+1), DIM*Np), intent(in) :: u0_mpi

    real(kind=8), dimension(DIM*Np) :: tmp
    integer(kind=4) :: ie_loc, E1, E2, iface, ie_neigh_loc, n_neigh

    !> partial result for contribution from stiffness and forcing term
    real(kind=8), dimension(Np) :: v_ptr

    real(kind=8), dimension(Np,Np) :: R, RT

    integer(kind=4) :: i, istart, iend
    integer(kind=4) :: row, elem_proc_id

    n_neigh = PolyMesh%Elem_loc(1)%num_faces

    !allocate(v_ptr(Np))

    elem_loop: do ie_loc = 1,PolyMesh%num_elem_loc
        E1 = ie_loc
        tmp = 0.0d0

        ! E+ contribution
        tmp = matmul(K_loc(E1,1,:,:), u0_loc(E1,:))

        ! E- contributions
        neigh_loop: do iface=1,n_neigh

            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)

            ! if boundary face, no contribution in E-, Dirichlet contribution already in E+
            if (E2 < 0) then

                !print *, '-- cycle --'
                !if (mpi_id==0) print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, ' - BOUNDARY'

                cycle neigh_loop
            endif

            elem_proc_id = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0)
            ! check if neighbor is in the same process
            if (elem_proc_id == mpi_id) then

                ie_neigh_loc = PolyMesh%elem_glo2loc(E2)

                !print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, ' at', ie_neigh_loc, ' column: ', iface+1, ' - LOCAL'

                ! use u0_loc directly
                tmp = tmp + matmul(K_loc(E1,iface+1,:,:), u0_loc(ie_neigh_loc,:))

            ! find index in u_mpi, given the global index of E-
            else
                ! start and end index on u_mpi
                istart = PolyMesh%inter_disp(mpi_id+1,1)+1
                iend = istart + PolyMesh%num_elem_inter_vec(mpi_id+1)

                ! finds first occurrence in full list
                do i = istart, iend
                    if (PolyMesh%elem_inter_glo(i) == E2) then
                        ie_neigh_loc = i
                        exit
                    endif
                enddo
                ! rescale to local vector index
                ie_neigh_loc = ie_neigh_loc - PolyMesh%inter_disp(mpi_id+1,1)

                !print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, ' at', ie_neigh_loc, ' column: ', iface+1

                ! use u_mpi directly
                tmp = tmp + matmul(K_loc(E1,iface+1,:,:), u0_mpi(ie_neigh_loc,:))
            endif

        end do neigh_loop

        !print *, 'mpi_id:', mpi_id, ' ---- end neigh loop', ie_loc, '----'

        ! add forcing term and rescale
        tmp = - tmp + rhs_stat_loc(ie_loc,:) + rhs_dyn_loc(ie_loc,:) * time_function(t)

        ! mass linear system in matrix-free
        do i=1,DIM

            v_ptr = tmp(row+1:row+Np)

            R = R_M_loc(ie_loc,i,:,:)
            RT = transpose(R)
            tmp(row+1:row+Np) = solve_LU(RT, R, v_ptr, Np)

        enddo

        ! sum initial condition contributions
        usol_loc(ie_loc,:) = half_dt2*tmp + u0_loc(ie_loc,:) + time_step*v0_loc(ie_loc,:)

    end do elem_loop

    !deallocate(v_ptr)
    call LYMPH3D_BARRIER

end subroutine FIRST_TIME_STEP_MATRIX_FREE

!> @brief solver for matrix free considering only one element with time dependence
subroutine TIME_STEP_MATRIX_FREE(PolyMesh, Np, t, K_loc, R_M_loc, rhs_stat_loc, rhs_dyn_loc, u0_loc, un_loc, usol_loc, un_mpi)

    !TODO consider different Np
    !DONE sides different for each row of K

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh

    !> current time
    real(kind=8), intent(in) :: t
    !> degrees of freedom for local polynomial basis
    integer(kind=4), intent(in) :: Np
    !> mass matrix factorization for E+
    real(kind=8), dimension(:,:,:,:), intent(in) :: R_M_loc

    !> stiffness matrix with contributions only from element E+ and neighbors E-
    real(kind=8), dimension(:,:,:,:), intent(in) :: K_loc
    !> forcing term vector on E+ - static component
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: rhs_stat_loc
    !> forcing term vector on E+ - dynamic component
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: rhs_dyn_loc

    !> solution from two previous timesteps u^{(n-1)} for element E+
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: u0_loc
    !> solution from previous timestep for u^{(n)} with contributions from E+ and neighbors E-
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(in) :: un_loc
    !> output solution u^{n+1} for the element E+
    real(kind=8), dimension(PolyMesh%num_elem_loc, DIM*Np), intent(out) :: usol_loc
    !> u^{n} solution E- from other processes
    real(kind=8), dimension(PolyMesh%num_elem_inter_vec(mpi_id+1), DIM*Np), intent(in) :: un_mpi

    real(kind=8), dimension(DIM*Np) :: tmp
    integer(kind=4) :: ie_loc, E1, E2, iface, ie_neigh_loc, n_neigh

    !> partial result for contribution from stiffness and forcing term
    real(kind=8), dimension(Np) :: v_ptr

    real(kind=8), dimension(Np,Np) :: R, RT

    integer(kind=4) :: i, istart, iend
    integer(kind=4) :: row, elem_proc_id

    real(kind=8) :: t1,t2

    ! u_loc = M_loc^-1 * ( - dt2 * K_loc*un_loc  + dt2*rhs_loc*time_function(t) ) +  2*un_loc - u0_loc

    !allocate(v_ptr(Np))

    n_neigh = PolyMesh%Elem_loc(1)%num_faces

    elem_loop: do ie_loc = 1,PolyMesh%num_elem_loc
        E1 = ie_loc
        tmp = 0.0d0

        !----
        t1 = MPI_WTIME()

        ! E+ contribution
        tmp = matmul(K_loc(E1,1,:,:), un_loc(E1,:))

        ! E- contributions
        neigh_loop: do iface=1,n_neigh

            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)

            ! if boundary face, no contribution in E-, Dirichlet contribution already in E+
            if (E2 < 0) then
                !print *, '-- cycle --'
                !if (mpi_id==0) print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, ' - BOUNDARY'
                cycle neigh_loop
            endif

            elem_proc_id = PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0)
            ! check if neighbor is in the same process
            if (elem_proc_id == mpi_id) then

                ie_neigh_loc = PolyMesh%elem_glo2loc(E2)

                !print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, 'at', ie_neigh_loc, 'column: ', iface+1, ' - LOCAL'
                ! use un_loc directly

                tmp = tmp + matmul(K_loc(E1,iface+1,:,:), un_loc(ie_neigh_loc,:))

            ! find index in u_mpi, given the global index of E-
            else
                ! start and end index on u_mpi
                istart = PolyMesh%inter_disp(mpi_id+1,1)+1
                iend = istart + PolyMesh%num_elem_inter_vec(mpi_id+1)

                ! finds first occurrence in full list
                do i = istart, iend
                    if (PolyMesh%elem_inter_glo(i) == E2) then
                        ie_neigh_loc = i
                        exit
                    endif
                enddo
                ! rescale to local vector index
                ie_neigh_loc = ie_neigh_loc - PolyMesh%inter_disp(mpi_id+1,1)

                !print *, 'mpi_id:', mpi_id, ' ie_loc:', ie_loc, ' E2:', E2, 'at', ie_neigh_loc, 'column: ', iface+1
                ! use u_mpi directly
                tmp = tmp + matmul(K_loc(E1,iface+1,:,:), un_mpi(ie_neigh_loc,:))
            endif

        end do neigh_loop

        t2 = MPI_WTIME()
        tp_KU = tp_KU + t2 - t1
        !----

        ! add element forcing term and rescale
        tmp = - tmp + rhs_stat_loc(ie_loc,:) + rhs_dyn_loc(ie_loc,:) * time_function(t)

        ! element mass linear system - matrix-free in each dimension
        do i=1,DIM

            row = (i-1)*Np

            ! copy vector
            !----
            t1 = MPI_WTIME()

            v_ptr = tmp(row+1:row+Np)

            t2 = MPI_WTIME()
            tp_copy_vector = tp_copy_vector + t2 - t1
            !----

            ! solve linear system matrix-free on block (i,i)
            !----
            t1 = MPI_WTIME()

            R = R_M_loc(ie_loc,i,:,:)
            RT = transpose(R)
            usol_loc(ie_loc, row+1:row+Np) = solve_LU(RT, R, v_ptr, Np)

            t2 = MPI_WTIME()
            tp_linear_system = tp_linear_system + t2 - t1
            !----

        enddo

        ! sum initial condition contributions
        usol_loc(ie_loc,:) = dt2*usol_loc(ie_loc,:) + 2.0d0*un_loc(ie_loc,:) - u0_loc(ie_loc,:)

    end do elem_loop

    !deallocate(v_ptr)

    call LYMPH3D_BARRIER

end subroutine TIME_STEP_MATRIX_FREE

!> Compute L2 error in matrix-free framework
subroutine COMPUTE_ERROR_L2_MATRIX_FREE(PolyMesh, Np, M_modal_loc, uh_loc, uex_loc, err_L2_loc)

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    integer(kind=4), intent(in) :: Np               !< number of element dofs per direction
    real(kind=8), dimension(:,:,:,:), intent(in):: M_modal_loc  !< modal mass matrix
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np), intent(in) :: uh_loc !< computed solution
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np), intent(in) :: uex_loc    !< exact solution
    real(kind=8), intent(out) :: err_L2_loc     !< process local sum of element L2 errors

    real(kind=8) :: err_L2_el

    real(kind=8), dimension(Np) :: du, tmp
    integer(kind=4) :: i,ie_loc,row

    err_L2_loc = 0.0d0

    do ie_loc=1,PolyMesh%num_elem_loc
        do i=1,3
            row = (i-1)*Np
            du = uh_loc(ie_loc,row+1:row+Np) - uex_loc(ie_loc,row+1:row+Np)
            tmp = matmul(M_modal_loc(ie_loc,i,:,:), du)

            err_L2_el = dot_product(du, tmp)
            err_L2_loc = err_L2_loc + err_L2_el

        enddo
    enddo

end subroutine COMPUTE_ERROR_L2_MATRIX_FREE

!> Compute DG error in matrix-free framework
subroutine COMPUTE_ERROR_DG_MATRIX_FREE(PolyMesh, Np, A_dg_loc, uh_loc, uex_loc, uh_mpi, uex_mpi, err_DG_loc)

    implicit none

    type(Mesh_Structure), intent(in) :: PolyMesh    !< mesh
    integer(kind=4), intent(in) :: Np       !< number of element dofs per direction
    real(kind=8), dimension(:,:,:,:), intent(in):: A_dg_loc !< local DG matrix
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np), intent(in) :: uh_loc !< computed solution
    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np), intent(in) :: uex_loc !< exact solution
    real(kind=8), dimension(PolyMesh%num_elem_inter_vec(mpi_id+1),DIM*Np), intent(in) :: uh_mpi !< computed solution from other processes
    real(kind=8), dimension(PolyMesh%num_elem_inter_vec(mpi_id+1),DIM*Np), intent(in) :: uex_mpi !< exact solution from other processes
    real(kind=8), intent(out) :: err_DG_loc !< process local sum of element DG errors

    real(kind=8) :: err_DG_el

    real(kind=8), dimension(PolyMesh%num_elem_loc,DIM*Np) :: du
    real(kind=8), dimension(PolyMesh%num_elem_inter_vec(mpi_id+1),DIM*Np) :: du_mpi
    real(kind=8), dimension(DIM*Np) :: tmp
    integer(kind=4) :: ie_loc,n_neigh,iface,E2,ie_neigh_loc
    integer(kind=4) :: E1,i,iend,istart
    logical :: is_E2_local

    err_DG_loc = 0.0d0
    n_neigh = PolyMesh%Elem_loc(1)%num_faces

    du = uh_loc - uex_loc
    du_mpi = uh_mpi - uex_mpi

    do ie_loc=1,PolyMesh%num_elem_loc

        E1 = ie_loc

        ! E+ contribution
        tmp = matmul(A_DG_loc(ie_loc,1,:,:), du(ie_loc,:))

        ! E- contributions
        neigh_loop: do iface=1,n_neigh

            is_E2_local = .false.
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(iface,2)

            ! if boundary face, no contribution in E-, Dirichlet contribution already in E+
            if (E2 < 0) then
                cycle neigh_loop
            endif

            ! check if neighbor is in the same process
            if (PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) == mpi_id) then
                is_E2_local = .true.
            endif

            if (is_E2_local) then
                ie_neigh_loc = PolyMesh%elem_glo2loc(E2)
                ! use un_loc directly
                tmp = tmp + matmul(A_DG_loc(E1,iface+1,:,:), du(ie_neigh_loc,:))

            ! find index in u_mpi, given the global index of E-
            else
                ! start and end index on u_mpi
                istart = PolyMesh%inter_disp(mpi_id+1,1)+1
                iend = istart + PolyMesh%num_elem_inter_vec(mpi_id+1)

                ! finds first occurrence in full list
                do i = istart, iend
                    if (PolyMesh%elem_inter_glo(i) == E2) then
                        ie_neigh_loc = i
                        exit
                    endif
                enddo
                ! rescale to local vector index
                ie_neigh_loc = ie_neigh_loc - PolyMesh%inter_disp(mpi_id+1,1)

                ! use u_mpi directly
                tmp = tmp + matmul(A_DG_loc(E1,iface+1,:,:), du_mpi(ie_neigh_loc,:))
            endif

        end do neigh_loop

        err_DG_el = dot_product(du(ie_loc,:), tmp)
        err_DG_loc = err_DG_loc + err_DG_el

    enddo

end subroutine COMPUTE_ERROR_DG_MATRIX_FREE

end module matrix_free