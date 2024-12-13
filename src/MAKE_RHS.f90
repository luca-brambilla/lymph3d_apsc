subroutine MAKE_RHS(PolyMesh, PolyData, petsc_num, global_dof, Np, petsc_rhs)

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
    use global_parameters

    implicit none
    
    real(kind=8) :: present = 0.0
    ! petsc_rhs is provided by SET_PETSC_VECTOR.f90

    Vec :: petsc_rhs
    PetscScalar :: val(1)
    PetscInt :: irow(1)

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

    real(kind=8) :: lambda, mu, rho !! DENSITY USED FOR DYNAMICS
    integer(kind=4) :: mat_id

    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, ipoly_loc, ipoly_glob, ipoly2_loc, ipoly2_glob
    integer(kind=4) :: n_tet_in_poly, iface_poly
    integer(kind=4) :: Npoly
    integer(kind=4) :: i, j, m
    integer(kind=4) :: beg

    integer(kind=4) :: e, E1, E2
    integer(kind=4), dimension(4) :: face_flag
    real(kind=8), dimension(3) :: nn
    integer(kind=4) :: space_fun_tag

    real(kind=8), dimension(3,Np) :: rhs_tet_loc
    real(kind=8), dimension(3,Np) :: rhs_face_bd_loc

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

    ! Map to the reference tetrahedron and triangle (see Poly_ref_mappings.f90)
    call mapping_quadrature_3D(nod3, wei3, nq3, nodtet3, weitet3)
    call mapping_quadrature_2D(nod2, wei2, nq2, nodtria2, weitria2)

    ! list of the degrees of monomials of the Np basis functions up to order p (see basis_functions.f90)
    allocate(blist(Np,3))
    call basis_list(blist, p, Np)

    print *,'Assembling rhs...'

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

    ! loop on the tetrahedra
    do ie_loc = 1, PolyMesh%num_elem_loc

        ! initialization of the rhs term on the volume rhs_tet_loc       
        rhs_tet_loc = 0.0

        mat_id = PolyMesh%Elem_loc(ie_loc)%mat_prop

        !!! check if there is a better way than this in loop over elements
        ! take correct density only for dynamic case 
        rho = PolyData%prop_mat(mat_id,1) * present !! present = 1.0 if DYNAMIC PROBLEM, 0.0 otherwise
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

        ! computation of the rhs term on the volume rhs_tet_loc (see assemble_element.f90)
        call MAKE_RHS_VOLUME(Np, Fk, Jdet, nodtet3, weitet3, nq3, lambda, mu, phi, rhs_tet_loc, rho)

        ! this allows to assemble the local vector correctly into the global vector
        beg = (ipoly_glob-1)*Np + 1
        
        ! insert the values of rhs_tet_loc in the entries of the global rhs vector
        do i=1,3
            do m=1,Np

                val(1)  = rhs_tet_loc(i,m)
                irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)

                if (val(1) .ne. 0.0) then
                    ! set value val to the vector petsc_rhs in the irow-th row
                    PetscCall(VecSetValues(petsc_rhs, 1, irow, val, ADD_VALUES, mpi_ierr))
                endif

            enddo
        enddo
        
        E1 = ie_loc

        ! begin loop on the faces of the tetrahedron E1
        do e=1,PolyMesh%Elem_loc(E1)%num_faces

            face_flag(e) = 0

            ! initialization of the face rhs term rhs_face_bd_loc
            rhs_face_bd_loc = 0.0

            ! find the neighbouring tetrahedron E2 sharing the face e with E1
            E2 = PolyMesh%Elem_loc(E1)%neigh_el(e,2)
            space_fun_tag = PolyMesh%Elem_loc(E1)%neigh_el(e,5)

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

                ! If it is true, then the two polyhedra do not belong to the same processo
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

                    call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly),&
                                        nn,PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc)

                else

                    if (E2 /= -1 .and. E2 /= -2) then

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                            PolyMesh%Poly(ipoly2_loc)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(ipoly2_loc)%hk,nn, &
                                            PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc)

                    else

                        call basis_boundary(phi_b,grad_b,e,E2,PolyMesh%Poly(ipoly_loc)%b_box,&
                                            PolyMesh%Poly(1)%b_box,blist, Np, Fk, node_maps, nodtria2, nq2)

                        call MAKE_RHS_FACE(theta,alpha,p,Np,e,E2,PolyMesh%Poly(ipoly_loc)%hk,PolyMesh%Poly(1)%hk,nn, &
                                            PolyMesh%Elem_loc(E1)%area(e),Fk,nodtria2,weitria2,nq2,lambda,mu,node_maps,phi_b,grad_b,space_fun_tag,rhs_face_bd_loc)

                    endif

                endif
                  
                ! if e is a boundary edge, then insert the values of rhs_face_bd_loc 
                ! in the entries of the global rhs vector 
                if (E2 == -1 .or. E2 == -2) then

                    do i=1,3
                        do m=1,Np

                            val(1)  = rhs_face_bd_loc(i,m)
                            irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)

                            if (val(1) .ne. 0.0) then
                                ! set value val to the vector petsc_rhs in the irow-th row
                                PetscCall(VecSetValues(petsc_rhs, 1, irow, val, ADD_VALUES, mpi_ierr))
                            endif

                        enddo
                    enddo

                endif

            endif
        
        enddo
    enddo
 
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    ! Assembly of the global rhs vector petsc_rhs
    PetscCall(VecAssemblyBegin(petsc_rhs, mpi_ierr))
    PetscCall(VecAssemblyEnd(petsc_rhs, mpi_ierr))

    deallocate(phi)
    deallocate(dphi,phi_b)
    deallocate(grad_b)
    deallocate(nodtet3)
    deallocate(nodtria2)
    deallocate(weitet3)
    deallocate(weitria2)

    deallocate(blist)

    print *, 'Done with assembling rhs'
            
end subroutine MAKE_RHS