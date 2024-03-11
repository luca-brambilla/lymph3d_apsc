subroutine COMPUTE_MODAL_COEFFICIENTS(PolyMesh, petsc_num, global_dof, local_dof, Np, petsc_modal_coeff_uex)
     
#include<petsc/finclude/petscksp.h>

    use petscksp
    use Poly_setup_mpi
    use Poly_global
    use problem_data_and_properties
    use basis_function
    use Poly_ref_mappings
    use local_search
    use SET_PETSC_SYSTEM
        
    implicit none
        
    Vec :: petsc_modal_coeff_uex
    PetscScalar :: val(1)
    PetscInt :: irow(1)

    type(Mesh_Structure) :: PolyMesh
    integer(kind=4) :: p, Np
    integer(kind=4) :: nq3, nq2, global_dof, local_dof, Npoly
    integer(kind=4) :: ie_loc, ie_glob, ivert, id_node, i, m, beg
    integer(kind=4) :: ipoly_loc, ipoly_glob
    integer(kind=4), dimension(:,:), allocatable :: blist
    real(kind=8) :: Jdet
    real(kind=8), dimension(4) :: x, y, z
    real(kind=8), dimension(3,4) :: Fk
    real(kind=8), dimension(3,3) :: Jinv
    real(kind=8), dimension(:,:), allocatable :: nod2, nod3, nodtet3
    real(kind=8), dimension(:), allocatable :: wei2, wei3, weitet3
    real(kind=8), dimension(3) :: points
    real(kind=8), dimension(:,:), allocatable :: phi
    real(kind=8), dimension(:,:,:), allocatable :: dphi
    real(kind=8), dimension(:,:), allocatable :: modal_coeff_uex
    integer(kind=4) :: petsc_num(global_dof)
    real(kind=8), dimension(3) :: exact_sol

    integer(kind=4) :: q, ii, jj

    p = PolyMesh%Elem_loc(1)%Degree
    Npoly = PolyMesh%num_poly

    call SET_PETSC_VECTOR(petsc_modal_coeff_uex, local_dof, global_dof)

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
    allocate(modal_coeff_uex(3,Np))

    ! loop on the tetrahedra
    do ie_loc = 1, PolyMesh%num_elem_loc

        modal_coeff_uex = 0.0

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

                do ii=1,3
                    points(ii)=0.0
                    do jj=1,4
                        points(ii)=points(ii)+Fk(ii,jj)*nodtet3(jj,q)
                    enddo
                enddo

                exact_sol = uex(points)
                do i=1,3
                    modal_coeff_uex(i,m) = modal_coeff_uex(i,m) + abs(Jdet)*weitet3(q)*exact_sol(i)*phi(m,q)
                enddo
                
            end do
        end do

        ! this allows to assemble the local matrix correctly into the global vector
        beg = (ipoly_glob-1)*Np + 1

        do i=1,3
            do m=1,Np

            val(1) = modal_coeff_uex(i,m)
            irow(1) = petsc_num((i-1)*Np*Npoly + beg+m-1)

            if (val(1) .ne. 0.0) then
                PetscCall(VecSetValues(petsc_modal_coeff_uex, 1, irow, val, ADD_VALUES, mpi_ierr))
            endif

            enddo
        enddo
        
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

    PetscCall(VecAssemblyBegin(petsc_modal_coeff_uex,mpi_ierr))
    PetscCall(VecAssemblyEnd(petsc_modal_coeff_uex,mpi_ierr))

    print *, 'Done with computing modal coefficients'

    deallocate(phi)
    deallocate(dphi)
    
    deallocate(nodtet3)
    deallocate(weitet3)

end subroutine COMPUTE_MODAL_COEFFICIENTS