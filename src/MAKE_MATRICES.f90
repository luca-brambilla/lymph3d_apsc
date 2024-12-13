!> Set up PETSc matrices for stiffness, mass, DG and mass modal
subroutine MAKE_MATRICES(PolyMesh, PolyData, petsc_num, global_dof, Np, petsc_stiff, petsc_mass, mat_dg, petsc_mass_modal)

#include<petsc/finclude/petscksp.h>

    use petscksp
    use Poly_mesh
    use mpi
    use Poly_setup_mpi
    use problem_data_and_properties
    use basis_function
    use assemble_element
    use local_search
    use Poly_ref_mappings
    use Poly_data
    use mesh_partition_and_mpi_files

    implicit none

    ! petsc_stiff, petsc_mass and mat_dg are provided by SET_PETSC_MATRIX.f90

    Mat, intent(inout) :: petsc_stiff          !< stiffnes matrix
    Mat, intent(inout) :: petsc_mass           !< mass matrix
    Mat, intent(inout) :: mat_dg               !< DG matrix
    Mat, intent(inout) :: petsc_mass_modal     !< mass modal matrix
    PetscScalar :: val(1), val1(1), val2(1), val3(1), val4(1)
    PetscInt :: irow(1), jcol(1), irow2(1), jcol2(1)!, irow3(1), jcol3(1)

    type(Mesh_Structure), intent(inout) :: PolyMesh     !< Mesh
    type(Data_Structure), intent(in) :: PolyData        !< Data
    integer(kind=4), intent(in) :: Np                   !< Number of dofs for each element
    integer(kind=4), intent(in) :: global_dof           !< Number of global dofs
    integer(kind=4), dimension(global_dof), intent(in) :: petsc_num !< PETSc numberbering starting from 0 not 1

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

    real(kind=8), dimension(DIM,DIM+1) :: Fk
    real(kind=8) :: Jdet
    real(kind=8), dimension(3,3) :: Jinv
    real(kind=8), dimension(4) :: x, y, z

    real(kind=8) :: lambda, mu, rho !! DENSITY USED FOR DYNAMICS
    integer(kind=4) :: mat_id

    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, ipoly_loc, ipoly_glob, ipoly2_loc, ipoly2_glob
    integer(kind=4) :: n_tet_in_poly, iface_poly
    integer(kind=4) :: Npoly
    integer(kind=4) :: i, j, m, n
    integer(kind=4) :: beg, beg2

    integer(kind=4) :: e, E1, E2
    integer(kind=4), dimension(4) :: face_flag
    real(kind=8), dimension(3) :: nn

    real(kind=8), dimension(3,3,Np,Np) :: V_loc, M_loc
    real(kind=8), dimension(3,3,Np,Np) :: S_loc, I_loc, IN_loc, SN_loc

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

    print *, 'Assembling linear system...'

    allocate(phi(Np,nq3))
    allocate(dphi(3,Np,nq3))
    allocate(phi_b(Np,nq2,2))
    allocate(grad_b(3,Np,nq2,2))

    ! loop on the tetrahedra
    do ie_loc = 1, PolyMesh%num_elem_loc
        ! if (mod(ie_loc,50) == 0) print *, ie_loc
        ! initialization of V_loc and M_loc
        V_loc = 0.0
        M_loc = 0.0

        mat_id = PolyMesh%Elem_loc(ie_loc)%mat_prop
        rho = PolyData%prop_mat(mat_id,1) ! density used for dynamics
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

        ! computation of the local stiffness matrix V_loc (see assemble_element.f90)
        call MAKE_STIFFNESS_VOLUME(Np, Jdet, weitet3, nq3, lambda, mu, dphi, V_loc)

        ! computation of the local mass matrix M_loc (see assemble_element.f90)
        ! used rho=1.0 for modal matrix, then multiplied for true mass matrix
        call MAKE_MASS_VOLUME(Np, Jdet, weitet3, nq3, phi, 1.0d0, M_loc)

        ! this allows to assemble the local matrix correctly into the global matrices
        beg = (ipoly_glob-1)*Np + 1

        ! insert the values of V_loc in the entries of the global stiffness and DG matrices
        do i=1,3
            do j=1,3
                do m=1,Np
                    do n=1,Np

                        val(1)  = V_loc(i,j,m,n)
                        irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)
                        jcol(1) = petsc_num((j-1)*Np*Npoly + beg+n-1)

                        if (val(1) .ne. 0.0) then
                            ! set value val to the matrix petsc_stiff in the irow-th row and jcol-th column
                            PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                            PetscCall(MatSetValues(mat_dg, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                        endif

                    enddo
                enddo
            enddo
        enddo

        ! insert the values of M_loc in the entries of the global mass matrix
        do i=1,3
            do m=1,Np
                do n=1,Np

                    val(1)  = M_loc(i,i,m,n)
                    irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)
                    jcol(1) = petsc_num((i-1)*Np*Npoly + beg+n-1)

                    if (val(1) .ne. 0.0) then
                        ! set value val to the matrix petsc_mass and petsc_mass_modal in the irow-th row and jcol-th column
                        PetscCall(MatSetValues(petsc_mass_modal, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                        val(1) = val(1)*rho
                        PetscCall(MatSetValues(petsc_mass, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                    endif

                enddo
            enddo
        enddo

        E1 = ie_loc

        ! begin loop on the faces of the tetrahedron E1
        do e=1,PolyMesh%Elem_loc(E1)%num_faces

            face_flag(e) = 0

            ! initialization of the face matrices I_loc, S_loc, IN_loc and SN_loc
            I_loc = 0.0
            S_loc = 0.0
            IN_loc = 0.0
            SN_loc = 0.0

            ! find the neighbouring tetrahedron E2 sharing the face e with E1
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(e,2)

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

                    call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly), nn, &
                                        PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)

                else

                    ! check if e is not a boundary face
                    ! otherwise take a default "neighbouring" element (its information won't be read)
                    if (E2 /= -1 .and. E2 /= -2) then

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(ipoly2_loc)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(ipoly2_loc)%hk, nn, &
                                                PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)

                    else

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                                PolyMesh%Poly(1)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_STIFFNESS_FACE(alpha,p,Np,E2,PolyMesh%Poly(ipoly_loc)%hk, PolyMesh%Poly(1)%hk, nn, &
                                                PolyMesh%Elem_loc(E1)%area(e),weitria2,nq2,lambda,mu,phi_b,grad_b,S_loc,I_loc,IN_loc,SN_loc)

                    endif

                endif

            endif

            if(E2 /= -2) then

                ! insert the values of S_loc in the entries of the global stiffness and DG matrices
                do i=1,3
                    do j=1,3
                        do m=1,Np
                            do n=1,Np

                                val(1)  = S_loc(i,j,m,n)
                                irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)
                                jcol(1) = petsc_num((j-1)*Np*Npoly + beg+n-1)

                                if (val(1) .ne. 0.0) then
                                    PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                                    PetscCall(MatSetValues(mat_dg, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
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

                                    val(1) = - I_loc(j,i,n,m)
                                    val1(1) = theta * I_loc(i,j,m,n)
                                    irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)
                                    jcol(1) = petsc_num((j-1)*Np*Npoly + beg+n-1)

                                    if (val(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                                    endif

                                    if (val1(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val1, ADD_VALUES, mpi_ierr))
                                    endif

                                enddo
                            enddo
                        enddo
                    enddo

                else

                    ! this allows to assemble the local matrix correctly into the global matrices
                    beg2 = (ipoly2_glob-1)*Np + 1

                    do i=1,3
                        do j=1,3
                            do m=1,Np
                                do n=1,Np

                                    val(1) = - I_loc(j,i,n,m)
                                    val1(1) = theta * I_loc(i,j,m,n)
                                    val2(1) = - IN_loc(j,i,n,m)
                                    val3(1) = theta * IN_loc(i,j,m,n)
                                    val4(1) = SN_loc(i,j,m,n)

                                    irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)
                                    jcol(1) = petsc_num((j-1)*Np*Npoly + beg+n-1)

                                    irow2(1) = petsc_num((i-1)*Np*Npoly + beg2+m-1)
                                    jcol2(1) = petsc_num((j-1)*Np*Npoly + beg2+n-1)

                                    if (val(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val, ADD_VALUES, mpi_ierr))
                                    endif

                                    if (val1(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol, val1, ADD_VALUES, mpi_ierr))
                                    endif

                                    if (val2(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow2, 1, jcol, val2, ADD_VALUES, mpi_ierr))
                                    endif

                                    if (val3(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol2, val3, ADD_VALUES, mpi_ierr))
                                    endif

                                    if (val4(1) .ne. 0.0) then
                                        PetscCall(MatSetValues(petsc_stiff, 1, irow, 1, jcol2, val4, ADD_VALUES, mpi_ierr))
                                        PetscCall(MatSetValues(mat_dg, 1, irow, 1, jcol2, val4, ADD_VALUES, mpi_ierr))
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

    ! Assembly of stiffness matrix petsc_stiff
    PetscCall(MatAssemblyBegin(petsc_stiff, MAT_FINAL_ASSEMBLY, mpi_ierr))
    PetscCall(MatAssemblyEnd(petsc_stiff, MAT_FINAL_ASSEMBLY, mpi_ierr))

    ! Assembly of mass matrix petsc_mass
    PetscCall(MatAssemblyBegin(petsc_mass, MAT_FINAL_ASSEMBLY, mpi_ierr))
    PetscCall(MatAssemblyEnd(petsc_mass, MAT_FINAL_ASSEMBLY, mpi_ierr))

    ! Assembly of mass matrix mat_dg
    PetscCall(MatAssemblyBegin(mat_dg, MAT_FINAL_ASSEMBLY, mpi_ierr))
    PetscCall(MatAssemblyEnd(mat_dg, MAT_FINAL_ASSEMBLY, mpi_ierr))

    ! Assembly of mass matrix petsc_mass_modal
    PetscCall(MatAssemblyBegin(petsc_mass_modal, MAT_FINAL_ASSEMBLY, mpi_ierr))
    PetscCall(MatAssemblyEnd(petsc_mass_modal, MAT_FINAL_ASSEMBLY, mpi_ierr))

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

    PRINT *, 'Done with assembling linear system'

end subroutine MAKE_MATRICES