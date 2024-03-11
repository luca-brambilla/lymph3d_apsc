! Correct the mesh so that each polyhedron contains at least one tetrahedron
subroutine MESH_CORRECTION(PolyMesh,mpi_np,mpi_id)
          
    use Poly_mesh
    use find_poly
    use local_search

    implicit none
              
    type(Mesh_Structure) :: PolyMesh
          
    integer(kind=4), intent(in)  :: mpi_np, mpi_id
    integer(kind=4) :: ie, ie_loc
    integer(kind=4) :: ipoly, Kbeg, Kend, i, j
    integer(kind=4) :: len_missing
    integer(kind=4), dimension(:), allocatable :: index, missing

    ! set the boundings of the loop
    Kbeg = mpi_id * PolyMesh%num_poly_loc + 1
    Kend = (mpi_id + 1) * PolyMesh%num_poly_loc
    
    ! print *,'LOC:',Kbeg,Kend
    j=1
    allocate(missing(PolyMesh%num_elem_loc))
    do ipoly=Kbeg,Kend
        ! Store the indeces of tetrahedras contained in polyhedra ipoly
        ! see subroutine find_poly in Poly_global.f90
        index=FIND_TET_IN_POLY(PolyMesh%elem_in_poly_loc,ipoly,PolyMesh%num_elem_loc)
        if (index(1) == 0) then
            missing(j)=ipoly;
            !print *,'MISSING',ipoly
            j=j+1;
        end if
    end do
    len_missing=j-1;
    Kend=Kend-len_missing;
    !K_new=Kend;
    !print *,'IN QUANTITY',len_missing

    ! Take the indices of tetrahedras contained in a given polyhedron
    ! and assign them to the polyhedra which does not contain any tetrahedra
    i=Kbeg;
    do while(len_missing>0)
      index=FIND_TET_IN_POLY(PolyMesh%elem_in_poly_loc,i,PolyMesh%num_elem)
      print *, "index: ", index
      if (size(index)>1) then
          !do j=1,size(index)
          ie=index(1)
          !print *,'MISSING',missing(len_missing),'substitute',i,'ie',ie
          call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                      PolyMesh%num_elem_loc, &
                                      ie,ie_loc)
          if (ie_loc /= 0) then                
            !print *,'local',ie_loc
            PolyMesh%elem_in_poly_loc(ie_loc)=missing(len_missing)
            !end do
            len_missing=len_missing-1;
          endif
      endif
      i=i+1;
    end do

    PolyMesh%num_poly_loc=maxval(PolyMesh%elem_in_poly_loc)
    !print *,'NEW POLY LOC',PolyMesh%num_poly_loc
          
end subroutine MESH_CORRECTION