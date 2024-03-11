module post_processing

#include<petsc/finclude/petscksp.h>
    
    use petscksp
    use Poly_setup_mpi
    use Poly_global
    use problem_data_and_properties
    use Poly_mesh
    use Poly_ref_mappings
    use basis_function
    use MOD_VTK
    use SET_PETSC_SYSTEM

    implicit none

    contains

    ! Evaluate nodal values of the solution and store them in output.vtk 
    subroutine EXPORT_SOLUTION(PolyMesh, u, IsPoly, mpi_id)
        
        use local_search ! see Poly_global.f90
        use problem_data_and_properties

        implicit none

        type(Mesh_Structure), intent(inout) :: PolyMesh
        real(kind=8), dimension(:,:), allocatable, intent(in) :: u
        integer(kind=4), intent(in) :: mpi_id
        logical, intent(in) :: IsPoly
        integer(kind=4) :: p, Np, Npoly
        integer(kind=4), dimension(:,:), allocatable :: blist
        real(kind=8), dimension(:), allocatable :: x_p, y_p, z_p
        real(kind=8), dimension(:), allocatable :: valx, valy, valz
        real(kind=8), dimension(:), allocatable :: dvalx, dvaly, dvalz
        real(kind=8), dimension(2) :: intx, inty, intz
        real(kind=8), dimension(:,:), allocatable :: temp
        real(kind=8), dimension(:,:,:), allocatable :: u_nod_vet
        real(kind=8), dimension(3) :: points
        integer(kind=4) :: ie_loc, ie_glob, ipoly_loc, ipoly_glob, ivert, id_node, i, j, k, index

        Np = PolyMesh%Elem_loc(1)%NDof_loc
        Npoly = PolyMesh%num_poly
        p = PolyMesh%Elem_loc(1)%Degree

        print *, 'Saving the solution in a file .vtk ...'
        
        ! list of the degrees of monomials of the Np basis functions up to order p
        ! (see basis_functions.f90)
        allocate(blist(Np,3))
        call basis_list(blist,p,Np)
        
        allocate(u_nod_vet(3,4,PolyMesh%num_elem_loc))
        allocate(temp(4,PolyMesh%num_elem_loc))

        u_nod_vet = 0.0

        do k=1,3

            ! Computing nodal values of the solution (on each element)
            do ie_loc = 1,PolyMesh%num_elem_loc

                allocate(x_p(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(y_p(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(z_p(PolyMesh%Elem_loc(ie_loc)%num_vert))

                do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert
                
                    ! see MAKE_PARTITION_AND_MPI_FILES.f90
                    call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                        PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node)      
                    
                    x_p(ivert)=PolyMesh%coord_x(id_node)
                    y_p(ivert)=PolyMesh%coord_y(id_node)
                    z_p(ivert)=PolyMesh%coord_z(id_node)

                enddo

                ! find the polyhedron ipoly_glob that contains the tetrahedron ie_loc
                ie_glob = PolyMesh%elem_loc2glo(ie_loc)
                ipoly_glob = PolyMesh%elem_in_poly(ie_glob)

                ! see module local_search in Poly_global.f90
                call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                                            PolyMesh%num_poly_loc, &
                                            ipoly_glob,ipoly_loc)

                do j=1,2
                    intx(j)=PolyMesh%Poly(ipoly_loc)%b_box(1,j)
                    inty(j)=PolyMesh%Poly(ipoly_loc)%b_box(2,j)
                    intz(j)=PolyMesh%Poly(ipoly_loc)%b_box(3,j)
                end do

                allocate(valx(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(valy(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(valz(PolyMesh%Elem_loc(ie_loc)%num_vert))

                allocate(dvalx(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(dvaly(PolyMesh%Elem_loc(ie_loc)%num_vert))
                allocate(dvalz(PolyMesh%Elem_loc(ie_loc)%num_vert))

                index = (k-1)*Npoly + ipoly_glob

                ! evaluation of L_blist(i,j), j=1,2,3, at the vertices of the given tetrahedra
                ! and evaluation of nodal values of the solution
                do ivert=1,PolyMesh%Elem_loc(ie_loc)%num_vert
                    temp(ivert,ie_loc)=0.0
                    do i = 1,Np
                        
                        call LegendreP(valx, x_p, blist(i,1), intx, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvalx, x_p, blist(i,1), intx, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call LegendreP(valy, y_p, blist(i,2), inty, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvaly, y_p, blist(i,2), inty, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call LegendreP(valz, z_p, blist(i,3), intz, PolyMesh%Elem_loc(ie_loc)%num_vert)
                        call GradLegendreP(dvalz, z_p, blist(i,3), intz, PolyMesh%Elem_loc(ie_loc)%num_vert)

                        temp(ivert,ie_loc) = temp(ivert,ie_loc) + u(i,index)*valx(ivert)*valy(ivert)*valz(ivert) 

                    end do

                end do

                deallocate(x_p,y_p,z_p)
                deallocate(valx,valy,valz)
                deallocate(dvalx,dvaly,dvalz)

            enddo

            u_nod_vet(k,:,:) = temp

        enddo
        
        call WRITE_SOLUTION_VTK(PolyMesh%num_elem_loc, PolyMesh, u_nod_vet, IsPoly, mpi_id) ! see MOD_VTK.f90
        print *,'Done exporting solution'

    end subroutine EXPORT_SOLUTION

    ! Compute L2 norm square of the error (u - u_ex)
    subroutine COMPUTE_ERROR_L2(mass, petsc_sol, petsc_uex, global_dof, err_L2_mpi, local_dof)

        use Poly_setup_mpi, only : mpi_ierr

        implicit none

        Mat :: mass
        Vec :: petsc_sol, petsc_uex, temp, error, e_L2
        PetscScalar coeff

        PetscViewer viewer

        integer(kind=4) :: global_dof, local_dof
        real(kind=8), pointer, dimension(:) :: error_L2_pointer
        real(kind=8), dimension(local_dof) :: error_L2
        real(kind=8) :: err_L2_mpi

        coeff = -1

        call SET_PETSC_VECTOR(e_L2, local_dof, global_dof)
        call SET_PETSC_VECTOR(temp, local_dof, global_dof)
        call SET_PETSC_VECTOR(error, local_dof, global_dof)

        PetscCall(VecWAXPY(error, coeff, petsc_sol, petsc_uex, mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'error',viewer,mpi_ierr))
        ! PetscCallA(VecView(error,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
        
        PetscCall(MatMult(mass, error, temp, mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'temp',viewer,mpi_ierr))
        ! PetscCallA(VecView(temp,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        PetscCall(VecPointwiseMult(e_L2, error, temp, mpi_ierr))
        PetscCallA(VecGetArrayF90(e_L2, error_L2_pointer, mpi_ierr))
        error_L2(1:local_dof) = error_L2_pointer
        err_L2_mpi = sum(error_L2)

        ! print *, "errL2 ^2: ", err_L2_mpi

    end subroutine COMPUTE_ERROR_L2

    ! Compute DG norm square of the error (u - u_ex)
    subroutine COMPUTE_ERROR_DG(mat_dg, petsc_sol, petsc_uex, global_dof, err_DG_mpi, local_dof)

        use Poly_setup_mpi, only : mpi_ierr

        implicit none

        Mat :: mat_dg
        Vec :: petsc_sol, petsc_uex, temp, error, e_DG
        PetscScalar coeff

        PetscViewer viewer

        integer(kind=4) :: global_dof, local_dof
        real(kind=8), pointer, dimension(:) :: error_DG_pointer
        real(kind=8), dimension(local_dof) :: error_DG
        real(kind=8) :: err_DG_mpi

        coeff = -1

        call SET_PETSC_VECTOR(e_DG, local_dof, global_dof)
        call SET_PETSC_VECTOR(temp, local_dof, global_dof)
        call SET_PETSC_VECTOR(error, local_dof, global_dof)

        PetscCall(VecWAXPY(error, coeff, petsc_sol, petsc_uex, mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Error_poly_10',viewer,mpi_ierr))
        ! PetscCallA(VecView(error,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))
        
        PetscCall(MatMult(mat_dg, error, temp, mpi_ierr))

        ! PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Temp_poly_10',viewer,mpi_ierr))
        ! PetscCallA(VecView(temp,viewer,mpi_ierr))
        ! PetscCallA(PetscViewerDestroy(viewer,mpi_ierr))

        PetscCall(VecPointwiseMult(e_DG, error, temp, mpi_ierr))
        PetscCallA(VecGetArrayF90(e_DG, error_DG_pointer, mpi_ierr))
        error_DG(1:local_dof) = error_DG_pointer
        err_DG_mpi = sum(error_DG)

        ! print *, "errDG ^2: ", err_DG_mpi

    end subroutine COMPUTE_ERROR_DG

    ! Compute the maximum value of the diameter of the elements of the mesh
    function compute_hmax(PolyMesh) result(hmax)

        type(Mesh_Structure), intent(in) :: PolyMesh
        real(kind=8) :: hmax
        integer(kind=4) :: ipoly_loc
        real(kind=8) :: h

        hmax = PolyMesh%Poly(1)%hk
        
        do ipoly_loc=2,PolyMesh%num_poly_loc

            h = PolyMesh%Poly(ipoly_loc)%hk

            if(h > hmax) hmax = h

        enddo

    end function compute_hmax

end module post_processing