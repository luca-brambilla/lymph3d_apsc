!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!

!TODO do not pass mpi_id and mpi_np use Poly_setup_mpi for all

!> @brief Makes Partitioning and writes files *.mpi
!! @author Ilario Mazzieri
!> @date September, 2013
!> @version 1.0

!> @note A new partition is made if elementdomain.mpi is absent in the
!! workig directory or in the folder given in SPEED.input (MPIFILE).

module mesh_partition_and_mpi_files

   implicit none

contains


!> partition of the mesh into the different processors and stores the local properties of the mesh into the correspondent fields of the struct Mesh_Structure. Furthermore, this subroutine generates output files in the folders FILES_MPI;
subroutine MAKE_PARTITION_AND_MPI_FILES(PolyData,PolyMesh,npoly)

   use mpi
   use Poly_setup_mpi
   use Poly_global
   use Poly_exit_codes
   use Poly_data
   use Poly_mesh
   !use MOD_VTK

   implicit none

   type(Data_Structure), intent(inout) :: PolyData
   type(Mesh_Structure), intent(inout) :: PolyMesh

   character(len = 70) :: part_file
   integer(kind=4), parameter :: part_unit = 400
   integer(kind=4) :: dummy, i, npoly,npoly_loc
   integer(kind=4),dimension(:,:),allocatable :: con_tet_loc

   !write(*,*) PolyMesh%num_elem, PolyMesh%num_node
   !read(*,*)

   allocate(PolyMesh%part_elem(PolyMesh%num_elem));
   PolyMesh%part_elem = mpi_id;

   if (mpi_id == 0) then
      if(len_trim(folder_mpi) /= 70) then
         part_file = folder_mpi(1:len_trim(folder_mpi)) // '/elem4proc.mpi'
      else
         part_file = 'elem4proc.mpi'
      endif

      inquire(file=part_file,exist=IS_filefound)

      if(IS_filefound .eqv. .TRUE.) then
         write(*,'(A)') 'Reading the existing partition...'
      else
         write(*,'(A)') 'Making a new partition...'
         call MESH_PARTITIONING(folder_mpi, PolyMesh%num_elem, &
                                 PolyMesh%num_node, mpi_np, &
                                 PolyMesh%num_hex, PolyMesh%con_hex, &
                                 PolyMesh%num_tet, PolyMesh%con_tet, &
                                 PolyMesh%num_prysm, PolyMesh%con_prysm)
      endif

      open(part_unit,file=part_file)
      read(part_unit,*) dummy
      do i = 1, PolyMesh%num_elem
         read(part_unit,*) dummy, PolyMesh%part_elem(i)
      enddo
      close(part_unit)

      write(*,'(A)') 'Done.'

   endif

   if(mpi_np > 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
      call MPI_BCAST(PolyMesh%part_elem,PolyMesh%num_elem, &
                     MPI_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
   endif

   !ATTENTION PROBABLY THIS HAS TO BE DONE ONLY IF ELEM4POC.MPI IS ABSENT
   !IN FILES_MPI --- CHECK LATER ---

   if (mpi_id == 0) PRINT *,'WRITE PARTITION'

   call WRITE_PARTITION(folder_mpi, PolyMesh, mpi_np, mpi_id)

   !try to agglomerate tetrahedral mesh
   allocate(con_tet_loc(PolyMesh%num_elem_loc,5))

   allocate(PolyMesh%elem_in_poly_loc(PolyMesh%num_elem_loc))

   npoly_loc=npoly/mpi_np;
   ! if(mod(npoly,mpi_np) .ne. 0) then
   !    write(*,'(A,I0,A,I0,A,I0,A)'), "Since the number of processors (", mpi_np, ") is not a multiple of the number of polyhedra (", npoly, &
   !              ") then the mesh will be agglomerated with ", npoly_loc*mpi_np, " polyhedra"
   ! endif
   PolyMesh%num_poly=npoly;
   PolyMesh%num_poly_loc=npoly_loc;

   if (npoly/=PolyMesh%num_elem) then
      print *,'WRITE MESH AGGLOMERATION'
      call MESH_AGGLOMERATION(folder_mpi, PolyMesh, PolyMesh%num_elem_loc,PolyMesh%num_node,npoly_loc, con_tet_loc,mpi_np, mpi_id)
      ! call MESH_CORRECTION(PolyMesh, mpi_np)
   else
      PolyMesh%num_poly=npoly;
      PolyMesh%num_poly_loc=PolyMesh%num_elem_loc
      call WRITE_POLY_INFO(folder_mpi,PolyMesh,mpi_id,mpi_np)
   endif

   if (mpi_id == 0) print *,'CREATE LOCAL MESH STRUCTURE'

   !creation of local mesh PolyMesh
   call CREATE_LOCAL_MESH(folder_mpi, PolyMesh, PolyData, mpi_np, mpi_id)

   if (mpi_id == 0) print *,'CREATE LOCAL VERTEX LIST'

   !creation of the local vertex list
   call CREATE_VERT_LIST(grid_file, PolyMesh, mpi_id)

   if (mpi_id == 0) print *,'CREATE POLY LIST'

   !creation of the local polyhedra list
   call CREATE_POLY_LIST(PolyMesh,mpi_id)

   !print *,PolyMesh%elem_in_poly
   if (mpi_id == 0) print *,'WRITE ELEM POLY GLOBAL'

   call CREATE_GLOBAL_POLY_MAP(folder_mpi,PolyMesh%elem_in_poly,PolyMesh%num_elem,PolyMesh,mpi_id,mpi_np)
   PolyMesh%num_poly=maxval(PolyMesh%elem_in_poly)


   if (mpi_id == 0) print *,'PRINT LOCAL MESH STRUCTURE'

   ! call print_Local_Mesh_Structure_VTK(PolyMesh)

   if(PolyMesh%num_prysm > 0) deallocate(PolyMesh%con_prysm)
   if(PolyMesh%num_tet > 0)   deallocate(PolyMesh%con_tet)
   if(PolyMesh%num_hex > 0)   deallocate(PolyMesh%con_hex)

   if (mpi_id == 0) print *,'CREATE NORMAL FACE'

   !creation of the normal vector to each face element
   call CREATE_NORMAL_FACE(PolyMesh, mpi_id)

   if (mpi_id == 0)   print *,'CREATE BBOX EL'

   call CREATE_BBOX_EL(PolyMesh)

   !creation of the bounding box (polygonal elements)faces
   !hexa-hexa / hexa-pyramids
   !call CREATE_NEIGH_EL_QUAD(folder_mpi, PolyMesh, mpi_np, mpi_id)

   if (mpi_id == 0) print *,'CREATE NEIGH EL TRIA'

   !creation of neighbouring element list for tria interfaces
   !tria-tria / tria-pyramids
   call CREATE_NEIGH_EL_TRIA(folder_mpi, PolyMesh, mpi_np, mpi_id)

   !creation of the Local/Global Dof numbering
   call CREATE_LOCAL2GLOBAL_MAP(PolyMesh, PolyData, mpi_np, mpi_id)

   if (mpi_id == 0) PRINT *,'WRITE MESH INFO'

   !optional: write infos about mesh (for debug purposes)
   call WRITE_MESH_INFO(folder_mpi, PolyMesh)

   if (mpi_id == 0) PRINT *,'WRITE INTERFACE INFO'
   call WRITE_INTERFACE_INFO(folder_mpi, PolyMesh)

end subroutine MAKE_PARTITION_AND_MPI_FILES

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> performs a contigous partition of the mesh into different processors using METIS and writes the mpi file elem4proc.mpi;
!> Mesh partitioning using METIS. Create `elem4proc.mpi` assigning an element (poly)
!> to each processor
subroutine MESH_PARTITIONING(mpi_file, n_elem, nnode, nparts, &
                                   num_hex, con_hex, &
                                   num_tet, con_tet, &
                                   num_prysm, con_prysm)

   use metis_interface, only: idx_t, real_t, METIS_SetDefaultOptions, &
                              METIS_PartMeshNodal, METIS_NOPTIONS, &
                              METIS_OPTION_NUMBERING, METIS_OPTION_CONTIG, &
                              METIS_OK, METIS_PartMeshDual

   implicit none

   integer(idx_t), intent(in) :: n_elem    !< number of elements
   integer(idx_t), intent(in) :: nnode    !< number of nodes
   integer(idx_t), parameter :: ncommon = 3

   integer(idx_t) :: eptr(n_elem+1), &
                     eind(8*num_hex+4*num_tet+5*num_prysm)   !< arrays storing mesh structure
   integer(idx_t) :: epart(n_elem), npart(nnode)       ! element and node partition vectors

   integer(idx_t) :: opts(0:METIS_NOPTIONS-1), ios, objval
   integer(idx_t), intent(in) :: nparts         !< number of parts in partition (number of processes)
   integer(idx_t), pointer :: vsize(:)  => NULL()
   integer(idx_t), pointer :: vwgt(:)   => NULL()
   integer(idx_t), pointer :: tpwgts(:) => NULL()

   character(len=70) :: mpi_file, u_name

   integer(kind=4) :: num_hex, num_tet, num_prysm
   ! integer(kind=4), dimension(num_hex,9), intent(in) :: con_hex
   ! integer(kind=4), dimension(num_tet,5), intent(in) :: con_tet
   ! integer(kind=4), dimension(num_prysm,6), intent(in) :: con_prysm
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: con_hex
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: con_tet
   integer(kind=4), dimension(:,:), allocatable, intent(in) :: con_prysm

   integer(kind=4) :: ic, ie, i
   integer(kind=4) :: u_mpi = 400

   print *, 'subroutine mesh partioning'
   if(len_trim(mpi_file) /= 70) then
      u_name = mpi_file(1:len_trim(mpi_file)) // '/elem4proc.mpi'
   else
      u_name = 'elem4proc.mpi'
   endif

   !write(*,*) n_elem
   !read(*,*)

   if (nparts == 1) then
      open(u_mpi,file = u_name)
      write(u_mpi,*) n_elem
      do i = 1, n_elem
         write(u_mpi,*) i, 0
      enddo
      close(u_mpi)

   else

      !write(*,*) num_hex, num_tet, num_prysm
      !read(*,*)
      eptr(1) = 1;
      ic = 1;
      !first hex
      do ie = 1, num_hex
         eptr(ie+1) = eptr(ie) + 8
         do i = 2, 9
            eind(ic) = con_hex(ie,i)
            ic = ic + 1
         enddo
      enddo
      !then tetra
      do ie = 1, num_tet
         eptr(ie+num_hex+1) = eptr(ie+num_hex) + 4
         do i = 2, 5
            eind(ic) = con_tet(ie,i)
            ic = ic +1
         enddo
      enddo
      !finally pyramids
      do ie = 1, num_prysm
         eptr(ie+num_hex+num_tet+1) = eptr(ie+num_hex+num_tet) + 5
         do i = 2, 6
            eind(ic) = con_prysm(ie,i)
            ic = ic +1
         enddo
      enddo

      ! write(*,*) eptr
      ! write(*,*) eind
      ! read(*,*)


      ios = METIS_SetDefaultOptions(opts)
      if (ios /= METIS_OK) then
         write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
         error stop 1
      end if
      opts(METIS_OPTION_NUMBERING) = 1    ! Fortran-style numbering
      opts(METIS_OPTION_CONTIG) = 1       ! Force contigous partitions
      ! call print_metis_options(opts)


!        ios = METIS_PartMeshNodal(n_elem,nnode,eptr,eind,nparts=nparts,options=opts, &
!                              objval=objval,epart=epart,npart=npart)

      ios = METIS_PartMeshDual(n_elem,nnode,eptr,eind, vwgt, vsize, &
                              ncommon, nparts, tpwgts, options=opts, objval=objval, &
                              epart=epart, npart=npart)


      if (ios /= METIS_OK) then
         write(*,*) "METIS_PartMeshDual failed with error: ", ios
         error stop 1
      end if

      open(u_mpi,file = u_name)
      write(u_mpi,*) n_elem
      do i = 1, n_elem
         write(u_mpi,*) i, epart(i)-1
      enddo
      close(u_mpi)

   endif

end subroutine MESH_PARTITIONING

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> generates the polyhedral mesh by agglomerating the tetrahedral mesh read from the mesh file and writes the mpi file elem_in_poly.mpi;
!>
subroutine MESH_AGGLOMERATION(mpi_file, PolyMesh,nelem_loc, nnode, nparts, con_tet,mpi_np, mpi_id)

   use metis_interface, only: idx_t, real_t, METIS_SetDefaultOptions, &
                           METIS_PartMeshNodal, METIS_NOPTIONS, &
                           METIS_OPTION_NUMBERING, METIS_OPTION_CONTIG, &
                           METIS_OK, METIS_PartMeshDual
   use Poly_mesh
   implicit none

   integer(idx_t), intent(in) :: nelem_loc    ! number of elements
   integer(idx_t), intent(in) :: nnode    ! number of nodes
   integer(idx_t), parameter :: ncommon = 3 !4

   integer(idx_t) :: eptr(nelem_loc+1), &
                     eind(4*nelem_loc)   ! arrays storing mesh structure
   integer(idx_t) :: epart(nelem_loc), npart(nnode)       ! element and node partition vectors

   integer(idx_t) :: opts(0:METIS_NOPTIONS-1), ios, objval, nparts

   integer(idx_t), pointer :: vsize(:)  => NULL()
   integer(idx_t), pointer :: vwgt(:)   => NULL()
   integer(idx_t), pointer :: tpwgts(:) => NULL()

   character(len=70) :: mpi_file, u_name
   character(len=70) :: mpi_file_tet
   integer(kind=4), intent(in)  :: mpi_np, mpi_id
   integer(kind=4)  :: unit_mpi, num_elem_list

   integer(kind=4), dimension(:), allocatable :: vect_read

   !integer(kind=4) :: num_tet
   integer(kind=4), dimension(nelem_loc,5) :: con_tet

   type(Mesh_Structure), intent(inout) :: PolyMesh

   integer(kind=4) :: ic, ie, i
   integer(kind=4) :: u_mpi = 400

   !allocate(PolyMesh%poly_loc2glo(nelem_loc))
   !PolyMesh%poly_loc2glo = 0

   mpi_file_tet = 'con_tet_000000.mpi'

   unit_mpi = 40 + mpi_id
   if (mpi_id < 10) then
      write(mpi_file_tet(14:14),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(mpi_file_tet(13:14),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(mpi_file_tet(12:14),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(mpi_file_tet(11:14),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(mpi_file_tet(10:14),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(mpi_file_tet(9:14),'(i6)') mpi_id
   endif

   mpi_file_tet = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_tet

   !then tet
   open(unit_mpi,file=mpi_file_tet)
   read(unit_mpi,*) num_elem_list
   allocate(vect_read(5))

   do ie = 1, num_elem_list/5
      read(unit_mpi,*) vect_read
      con_tet(ie,:)=vect_read(1:5)
      !print *,con_tet(ie,:)
   enddo

   deallocate(vect_read)
   close(unit_mpi)

   u_name = 'elem_in_poly_000000.mpi'

   if (mpi_id < 10) then
      write(u_name(19:19),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(u_name(13:14),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(u_name(12:14),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(u_name(11:14),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(u_name(10:14),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(u_name(9:14),'(i6)') mpi_id
   endif

   if(len_trim(mpi_file) /= 70) then
      u_name = mpi_file(1:len_trim(mpi_file)) // '/' // u_name
   else
      u_name = 'elemInpoly'
   endif

   !print *,u_name
   eptr(1) = 1;
   ic = 1;

   !then tetra
   do ie = 1, nelem_loc
      eptr(ie+1) = eptr(ie) + 4
      do i = 2, 5
         eind(ic) = con_tet(ie,i)
         ic = ic +1
      enddo
   enddo

   ! write(*,*) eptr
   ! write(*,*) eind
   ! read(*,*)

   ios = METIS_SetDefaultOptions(opts)
   if (ios /= METIS_OK) then
      write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
      error stop 1
   end if

   opts(METIS_OPTION_NUMBERING) = 1    ! Fortran-style numbering
   opts(METIS_OPTION_CONTIG) = 1       ! Force contigous partitions
   ! call print_metis_options(opts)


   !        ios = METIS_PartMeshNodal(n_elem,nnode,eptr,eind,nparts=nparts,options=opts, &
   !                              objval=objval,epart=epart,npart=npart)

   ios = METIS_PartMeshDual(nelem_loc,nnode,eptr, eind, vwgt, vsize, &
                           ncommon, nparts, tpwgts, options=opts, objval=objval, &
                           epart=epart, npart=npart)

   if (ios /= METIS_OK) then
      write(*,*) "METIS_PartMeshDual failed with error: ", ios
      error stop 1
   end if

   open(u_mpi,file = u_name)
   write(u_mpi,*) nelem_loc
   do i = 1, nelem_loc
      ! print *, "poly: ", PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(i))
      !if(PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(i)) > nparts*mpi_id &
      !      .and. PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(i)) <= nparts*(mpi_id+1)) then
         write(u_mpi,*) PolyMesh%elem_loc2glo(i), nparts*mpi_id+epart(i) ! PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(i))
         ! print *, PolyMesh%elem_loc2glo(i), PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(i))
         ! print *, ""
      !PolyMesh%elem_in_poly_loc(i)=nparts*mpi_id+epart(i)
      ! print *,nparts*mpi_id+epart(i)
      !endif
   enddo
   close(u_mpi)

end subroutine MESH_AGGLOMERATION

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> ???
!> Write element list in each polyhedron. The polygons are numbered in ascending order locally in each process, the ordering depends on the partitioning starting from the process 0.
!> Write file `elem_in_poly_X.mpi`.
subroutine WRITE_POLY_INFO(mpi_file,PolyMesh,mpi_id,mpi_np)

   use mpi
   use Poly_setup_MPI, only: mpi_ierr
   use Poly_mesh
   use find_poly

   character(len=70) :: mpi_file, u_name
   integer(kind=4), intent(in)  :: mpi_np, mpi_id
   !integer(kind=4)  :: unit_mpi
   integer(kind=4) :: num_elem_send,last_elem_send
   !integer(kind=4), dimension(:), allocatable :: index_glob
   !integer(kind=4) :: num_faces_in_poly
   integer(kind=4) :: i,ip
   integer(kind=4) :: u_mpi = 400

   type(Mesh_Structure), intent(inout) :: PolyMesh

   num_elem_send=0;

   do ip=1,mpi_np
      if (mpi_id==ip-1) then
         do i=1,PolyMesh%num_elem_loc
            !print *,num_elem_send,PolyMesh%num_elem_loc,mpi_id
            !PolyMesh%elem_in_poly_loc(i)=PolyMesh%num_elem_loc*mpi_id+i;

            ! Polygon numbered locally, ordering given by partitioning
            PolyMesh%elem_in_poly_loc(i)=num_elem_send+i;

            !PolyMesh%elem_in_poly_loc(i)=PolyMesh%elem_loc2glo(i);

         enddo
         last_elem_send=PolyMesh%elem_in_poly_loc(i-1);
         !print *,last_elem_send;
      endif
      if (mpi_id==ip-1) num_elem_send=last_elem_send;

      call MPI_BCAST(num_elem_send,1,&
                     MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

   enddo

   u_name = 'elem_in_poly_000000.mpi'

   if (mpi_id < 10) then
      write(u_name(19:19),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(u_name(13:14),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(u_name(12:14),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(u_name(11:14),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(u_name(10:14),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(u_name(9:14),'(i6)') mpi_id
   endif

   if(len_trim(mpi_file) /= 70) then
      u_name = mpi_file(1:len_trim(mpi_file)) // '/' // u_name
   else
      u_name = 'elemInpoly'
   endif

   open(u_mpi,file = u_name)
   write(u_mpi,*) PolyMesh%num_elem_loc
   do i = 1, PolyMesh%num_elem_loc
      write(u_mpi,*) i, PolyMesh%elem_in_poly_loc(i)
   enddo
   close(u_mpi)

end subroutine

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> allocates and stores the field elem_in_poly using the information collected in elem_in_poly.mpi;
!>
subroutine CREATE_GLOBAL_POLY_MAP(mpi_file,vec_glob,num_elem,PolyMesh,mpi_id,mpi_np)

   use mpi
   use Poly_setup_MPI, only: mpi_ierr
   use Poly_mesh

   implicit none

   integer(kind=4) :: num_elem_loc,num_elem,num_elem_send
   integer(kind=4) :: ip,i,unit_mpi
   integer(kind=4),dimension(:),allocatable :: vec_loc
   integer(kind=4),dimension(:),allocatable :: vec_send
   integer(kind=4),dimension(num_elem) :: vec_glob
   integer(kind=4) :: mpi_id,mpi_np
   character(len=70) :: mpi_file
   character(len=70) :: mpi_file_poly

   type(Mesh_Structure), intent(inout) :: PolyMesh

   vec_glob=0;
   num_elem_loc=PolyMesh%num_elem_loc
   allocate(vec_loc(num_elem_loc))
   vec_loc=PolyMesh%elem_in_poly_loc;

   do ip = 1, mpi_np

      allocate(vec_send(num_elem_loc))

      if(mpi_id == ip-1) vec_send = vec_loc;

      !PRINT *,vec_loc

      if(mpi_id == ip-1) num_elem_send = num_elem;

      if (mpi_id==ip-1) then
         do i=1,num_elem_loc
            vec_glob(PolyMesh%elem_loc2glo((i))) = vec_send(i)
         enddo
      endif

      call MPI_BCAST(num_elem_send, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      !print *,num_elem_send

      call MPI_BCAST(vec_glob,num_elem,&
                  MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      deallocate(vec_send)

   enddo
   call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

   if (mpi_id==0) then

      mpi_file_poly =  mpi_file(1:len_trim(mpi_file)) // '/' // 'poly_global.mpi'

      unit_mpi   = 400 + mpi_id

      open(unit_mpi,file=mpi_file_poly)
      write(unit_mpi,*) num_elem

      do i = 1, PolyMesh%num_tet
         write(unit_mpi,*) vec_glob(i)
      enddo
   endif

   !print *,vec_glob

end subroutine CREATE_GLOBAL_POLY_MAP

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> writes the connettivity of tetrahedra and triangles for each processor in the mpi files con_tet.mpi and con_tri.mpi;
!> Write partitioned connectivity to file for each element type for each processor.
!> con_tet_X.mpi columns: material ID, global element, local face number, 3 global node IDs
subroutine WRITE_PARTITION(mpi_file, PolyMesh, mpi_np, mpi_id)

   use Poly_mesh
   use qsort

   implicit none

   character(len=70), intent(in) :: mpi_file
   character(len=70) :: mpi_file_hex, mpi_file_tet, mpi_file_pry, &
                        mpi_file_quad, mpi_file_tria
   integer(kind=4), intent(in)  :: mpi_np, mpi_id
   integer(kind=4)  :: ie, unit_mpi, kiter, unit_mpi_2, unit_mpi_3, &
                        num_hex_mpi, num_tet_mpi, num_prysm_mpi, &
                        num_node_hex_mpi, num_node_tet_mpi, num_node_prysm_mpi, &
                        num_node_quad_mpi, num_node_tria_mpi
   integer(kind=4), dimension(4) :: vec_sort_quad
   integer(kind=4), dimension(3) :: vec_sort_tria

   type(Mesh_Structure), intent(inout) :: PolyMesh

   !write(*,*) PolyMesh%part_elem
   !read(*,*)

   ! first count hex
   num_hex_mpi = 0;
   do ie = 1, PolyMesh%num_hex
   !   write(*,*) ie
      if(PolyMesh%part_elem(ie) == mpi_id) &
         num_hex_mpi = num_hex_mpi + 1
   enddo
   num_node_hex_mpi = 9*num_hex_mpi

   !write(*,*)  num_hex_mpi
   !read(*,*)

   ! then count tetra
   num_tet_mpi = 0;
   do ie = 1, PolyMesh%num_tet
   !   write(*,*) ie+PolyMesh%num_hex
      if(PolyMesh%part_elem(ie+PolyMesh%num_hex) == mpi_id) &
         num_tet_mpi = num_tet_mpi + 1

   enddo
   num_node_tet_mpi = 5*num_tet_mpi

   !write(*,*)  num_tet_mpi
   !read(*,*)

   ! finally count prysm
   num_prysm_mpi = 0;
   do ie = 1, PolyMesh%num_prysm
      if(PolyMesh%part_elem(ie+PolyMesh%num_hex+PolyMesh%num_tet) == mpi_id) &
               num_prysm_mpi = num_prysm_mpi + 1
   enddo
   num_node_prysm_mpi = 6*num_prysm_mpi

   num_node_quad_mpi = 7*(num_hex_mpi*6 + num_prysm_mpi*1)
   num_node_tria_mpi = 6*(num_tet_mpi*4 + num_prysm_mpi*4)

   ! LOCAL NUMBER OF MESH ELEMENTS
   PolyMesh%num_elem_loc = num_hex_mpi + num_tet_mpi + num_prysm_mpi
   allocate(PolyMesh%elem_loc2glo(PolyMesh%num_elem_loc))
   !print *,PolyMesh%num_elem_loc
   PolyMesh%elem_loc2glo = 0

   !allocate(PolyMesh%Poly2glo(PolyMesh%num_Poly))
   !PolyMesh%Poly2glo=0

   mpi_file_hex  = 'con_hex_000000.mpi'
   mpi_file_tet  = 'con_tet_000000.mpi'
   mpi_file_pry  = 'con_pry_000000.mpi'
   mpi_file_quad = 'con_qua_000000.mpi'
   mpi_file_tria = 'con_tri_000000.mpi'

   unit_mpi   = 40 + mpi_id
   unit_mpi_2 = 4000 + mpi_id
   unit_mpi_3 = 400000 + mpi_id

   if (mpi_id < 10) then
      write(mpi_file_hex(14:14),'(i1)') mpi_id
      write(mpi_file_tet(14:14),'(i1)') mpi_id
      write(mpi_file_pry(14:14),'(i1)') mpi_id
      write(mpi_file_quad(14:14),'(i1)') mpi_id
      write(mpi_file_tria(14:14),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(mpi_file_hex(13:14),'(i2)') mpi_id
      write(mpi_file_tet(13:14),'(i2)') mpi_id
      write(mpi_file_pry(13:14),'(i2)') mpi_id
      write(mpi_file_quad(13:14),'(i2)') mpi_id
      write(mpi_file_tria(13:14),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(mpi_file_hex(12:14),'(i3)') mpi_id
      write(mpi_file_tet(12:14),'(i3)') mpi_id
      write(mpi_file_pry(12:14),'(i3)') mpi_id
      write(mpi_file_quad(12:14),'(i3)') mpi_id
      write(mpi_file_tria(12:14),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(mpi_file_hex(11:14),'(i4)') mpi_id
      write(mpi_file_tet(11:14),'(i4)') mpi_id
      write(mpi_file_pry(11:14),'(i4)') mpi_id
      write(mpi_file_quad(11:14),'(i4)') mpi_id
      write(mpi_file_tria(11:14),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(mpi_file_hex(10:14),'(i5)') mpi_id
      write(mpi_file_tet(10:14),'(i5)') mpi_id
      write(mpi_file_pry(10:14),'(i5)') mpi_id
      write(mpi_file_quad(10:14),'(i5)') mpi_id
      write(mpi_file_tria(10:14),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(mpi_file_hex(9:14),'(i6)') mpi_id
      write(mpi_file_tet(9:14),'(i6)') mpi_id
      write(mpi_file_pry(9:14),'(i6)') mpi_id
      write(mpi_file_quad(9:14),'(i6)') mpi_id
      write(mpi_file_tria(9:14),'(i6)') mpi_id
   endif

   mpi_file_hex = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_hex
   mpi_file_tet = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_tet
   mpi_file_pry = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_pry
   mpi_file_quad = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_quad
   mpi_file_tria = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_tria

   kiter = 1;

   ! Hexhedron
   open(unit_mpi,file=mpi_file_hex)
   write(unit_mpi,*) num_node_hex_mpi

   open(unit_mpi_2,file=mpi_file_quad)
   write(unit_mpi_2,*) num_node_quad_mpi

   do ie = 1, PolyMesh%num_hex
      if(PolyMesh%part_elem(ie) == mpi_id) then
         write(unit_mpi,*) PolyMesh%con_hex(ie,1:9)
         PolyMesh%elem_loc2glo(kiter) = ie
         !PolyMesh%Poly2glo(kiter) = PolyMesh%elem_in_poly(ie)
         kiter = kiter + 1

         ! Permutation of 6 nodes to get the 6 faces of the tetrahedron

         !face 1
         vec_sort_quad = [PolyMesh%con_hex(ie,2), PolyMesh%con_hex(ie,3), &
                              PolyMesh%con_hex(ie,4), PolyMesh%con_hex(ie,5)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 1, vec_sort_quad(1:4)
         !face 2
         vec_sort_quad = [PolyMesh%con_hex(ie,6), PolyMesh%con_hex(ie,7), &
                              PolyMesh%con_hex(ie,8), PolyMesh%con_hex(ie,9)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 2, vec_sort_quad(1:4)
         !face 3
         vec_sort_quad = [PolyMesh%con_hex(ie,3), PolyMesh%con_hex(ie,4), &
                              PolyMesh%con_hex(ie,8), PolyMesh%con_hex(ie,7)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 3, vec_sort_quad(1:4)
         !face 4
         vec_sort_quad = [PolyMesh%con_hex(ie,4), PolyMesh%con_hex(ie,5), &
                              PolyMesh%con_hex(ie,9), PolyMesh%con_hex(ie,8)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 4, vec_sort_quad(1:4)
         !face 5
         vec_sort_quad = [PolyMesh%con_hex(ie,5), PolyMesh%con_hex(ie,2), &
                              PolyMesh%con_hex(ie,6), PolyMesh%con_hex(ie,9)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 5, vec_sort_quad(1:4)
         !face 6
         vec_sort_quad = [PolyMesh%con_hex(ie,2), PolyMesh%con_hex(ie,3), &
                              PolyMesh%con_hex(ie,7), PolyMesh%con_hex(ie,6)]
         call QsortC(vec_sort_quad)
         write(unit_mpi_2,*) PolyMesh%con_hex(ie,1), ie, 6, vec_sort_quad(1:4)
      endif
   enddo
   close(unit_mpi)
   close(unit_mpi_2)

   ! -------
   ! Tetrahedron
   open(unit_mpi,file=mpi_file_tet)
   write(unit_mpi,*) num_node_tet_mpi

   open(unit_mpi_2,file=mpi_file_tria)
   write(unit_mpi_2,*) num_node_tria_mpi

   do ie = 1, PolyMesh%num_tet
      if(PolyMesh%part_elem(ie+PolyMesh%num_hex) == mpi_id) then
         write(unit_mpi,*) PolyMesh%con_tet(ie,1:5)
         PolyMesh%elem_loc2glo(kiter) = ie+PolyMesh%num_hex
         !PolyMesh%Poly2glo(kiter) = PolyMesh%elem_in_poly(ie+PolyMesh%num_hex)
         kiter = kiter + 1

         ! Permutation of 4 nodes to get the 4 faces of the tetrahedron

         !face 1
         vec_sort_tria = [PolyMesh%con_tet(ie,2), &
                                 PolyMesh%con_tet(ie,3), PolyMesh%con_tet(ie,4)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_2,*) PolyMesh%con_tet(ie,1), &
                              ie+PolyMesh%num_hex, 1, vec_sort_tria(1:3)
         !face 2
         vec_sort_tria = [PolyMesh%con_tet(ie,2), &
                                 PolyMesh%con_tet(ie,3), PolyMesh%con_tet(ie,5)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_2,*) PolyMesh%con_tet(ie,1), &
                              ie+PolyMesh%num_hex, 2, vec_sort_tria(1:3)
         !face 3
         vec_sort_tria = [PolyMesh%con_tet(ie,3), &
                                 PolyMesh%con_tet(ie,4), PolyMesh%con_tet(ie,5)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_2,*) PolyMesh%con_tet(ie,1), &
                              ie+PolyMesh%num_hex, 3, vec_sort_tria(1:3)
         !face 4
         vec_sort_tria = [PolyMesh%con_tet(ie,4), &
                                 PolyMesh%con_tet(ie,2), PolyMesh%con_tet(ie,5)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_2,*) PolyMesh%con_tet(ie,1), &
                              ie+PolyMesh%num_hex, 4, vec_sort_tria(1:3)

      endif
   enddo
   !print *,'kiter',kiter
   close(unit_mpi)
   close(unit_mpi_2)

   ! ---------
   ! Prysm
   open(unit_mpi,file=mpi_file_pry)
   write(unit_mpi,*) num_node_prysm_mpi

   open(unit_mpi_2,file=mpi_file_quad, position='append')
   open(unit_mpi_3,file=mpi_file_tria, position='append')

   do ie = 1, PolyMesh%num_prysm
      if(PolyMesh%part_elem(ie+PolyMesh%num_hex+PolyMesh%num_tet) == mpi_id) then
         write(unit_mpi,*) PolyMesh%con_prysm(ie,1:6)

         PolyMesh%elem_loc2glo(kiter) = ie+PolyMesh%num_hex+PolyMesh%num_tet
         !PolyMesh%Poly2glo(kiter) = PolyMesh%elem_in_poly(ie+PolyMesh%num_hex+PolyMesh%num_tet)
         kiter = kiter + 1
         !base
         vec_sort_quad = [PolyMesh%con_prysm(ie,2), PolyMesh%con_prysm(ie,3), &
                           PolyMesh%con_prysm(ie,4), PolyMesh%con_prysm(ie,5)]
         call QsortC(vec_sort_quad)

         write(unit_mpi_2,*) PolyMesh%con_prysm(ie,1), &
                           ie+PolyMesh%num_hex+PolyMesh%num_tet, 1, vec_sort_quad(1:4)
         !face 1
         vec_sort_tria = [PolyMesh%con_prysm(ie,2), PolyMesh%con_prysm(ie,3), &
                           PolyMesh%con_prysm(ie,6)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_3,*) PolyMesh%con_prysm(ie,1), &
                           ie+PolyMesh%num_hex+PolyMesh%num_tet, 2,  &
                                    vec_sort_tria(1:3)
         !face 2
         vec_sort_tria = [PolyMesh%con_prysm(ie,3), PolyMesh%con_prysm(ie,4), &
                           PolyMesh%con_prysm(ie,6)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_3,*) PolyMesh%con_prysm(ie,1), &
                           ie+PolyMesh%num_hex+PolyMesh%num_tet, 3, &
                                    vec_sort_tria(1:3)
         !face 3
         vec_sort_tria = [PolyMesh%con_prysm(ie,4), PolyMesh%con_prysm(ie,5), &
                           PolyMesh%con_prysm(ie,6)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_3,*) PolyMesh%con_prysm(ie,1), &
                           ie+PolyMesh%num_hex+PolyMesh%num_tet, 4, &
                                    vec_sort_tria(1:3)
         !face 4
         vec_sort_tria = [PolyMesh%con_prysm(ie,5), PolyMesh%con_prysm(ie,2), &
                           PolyMesh%con_prysm(ie,6)]
         call QsortC(vec_sort_tria)
         write(unit_mpi_3,*) PolyMesh%con_prysm(ie,1), &
                           ie+PolyMesh%num_hex+PolyMesh%num_tet, 5, &
                                    vec_sort_tria(1:3)
         endif
   enddo
   close(unit_mpi)
   close(unit_mpi_2)
   close(unit_mpi_3)

end subroutine WRITE_PARTITION


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> reads the local properties of the tetrahedral mesh from the mpi files and stores them into the corresponding fields of the struct Mesh_Structure;
!>
      subroutine CREATE_LOCAL_MESH(mpi_file, PolyMesh, PolyData, mpi_np, mpi_id)

      use Poly_data
      use Poly_mesh

      implicit none

      character(len=70), intent(in) :: mpi_file
      character(len=70) :: mpi_file_hex, mpi_file_tet, mpi_file_pry,u_name
      integer(kind=4), intent(in)  :: mpi_np, mpi_id
      integer(kind=4)  :: i, ie, unit_mpi, num_elem_list, kiter, num_quad_loc, &
                          num_tria_loc, ie_shift
      integer(kind=4), dimension(:), allocatable :: vect_read

      type(Mesh_Structure), intent(inout) :: PolyMesh
      type(Data_Structure), intent(inout) :: PolyData

      mpi_file_hex = 'con_hex_000000.mpi'
      mpi_file_tet = 'con_tet_000000.mpi'
      mpi_file_pry = 'con_pry_000000.mpi'

      unit_mpi = 40 + mpi_id
      if (mpi_id < 10) then
         write(mpi_file_hex(14:14),'(i1)') mpi_id
         write(mpi_file_tet(14:14),'(i1)') mpi_id
         write(mpi_file_pry(14:14),'(i1)') mpi_id
      elseif (mpi_id < 100) then
         write(mpi_file_hex(13:14),'(i2)') mpi_id
         write(mpi_file_tet(13:14),'(i2)') mpi_id
         write(mpi_file_pry(13:14),'(i2)') mpi_id
      elseif (mpi_id < 1000) then
         write(mpi_file_hex(12:14),'(i3)') mpi_id
         write(mpi_file_tet(12:14),'(i3)') mpi_id
         write(mpi_file_pry(12:14),'(i3)') mpi_id
      elseif (mpi_id < 10000) then
         write(mpi_file_hex(11:14),'(i4)') mpi_id
         write(mpi_file_tet(11:14),'(i4)') mpi_id
         write(mpi_file_pry(11:14),'(i4)') mpi_id
      elseif (mpi_id < 100000) then
         write(mpi_file_hex(10:14),'(i5)') mpi_id
         write(mpi_file_tet(10:14),'(i5)') mpi_id
         write(mpi_file_pry(10:14),'(i5)') mpi_id
      elseif (mpi_id < 1000000) then
         write(mpi_file_hex(9:14),'(i6)') mpi_id
         write(mpi_file_tet(9:14),'(i6)') mpi_id
         write(mpi_file_pry(9:14),'(i6)') mpi_id
      endif

      u_name = 'elem_in_poly_000000.mpi'

      if (mpi_id < 10) then
         write(u_name(19:19),'(i1)') mpi_id
      elseif (mpi_id < 100) then
         write(u_name(13:14),'(i2)') mpi_id
      elseif (mpi_id < 1000) then
         write(u_name(12:14),'(i3)') mpi_id
      elseif (mpi_id < 10000) then
         write(u_name(11:14),'(i4)') mpi_id
      elseif (mpi_id < 100000) then
         write(u_name(10:14),'(i5)') mpi_id
      elseif (mpi_id < 1000000) then
         write(u_name(9:14),'(i6)') mpi_id
      endif

      if(len_trim(mpi_file) /= 70) then
         u_name = mpi_file(1:len_trim(mpi_file)) // '/' // u_name
      else
         u_name = 'elemInpoly'
      endif


      mpi_file_hex = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_hex
      mpi_file_tet = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_tet
      mpi_file_pry = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_pry

      allocate(PolyMesh%Elem_loc(PolyMesh%num_elem_loc));
      !allocate(PolyMesh%Elem_loc(PolyMesh%num_Poly))

      ! first hex
      open(unit_mpi,file=mpi_file_hex)
      read(unit_mpi,*) num_elem_list
      allocate(vect_read(9))

      num_quad_loc = num_elem_list/9
      !kiter = 0
      do ie = 1, num_elem_list/9
         !do i = kiter + 1, kiter + 9
            read(unit_mpi,*) vect_read(1:9)
         ! enddo
          PolyMesh%Elem_loc(ie)%el_type  = 'HEX'
          PolyMesh%Elem_loc(ie)%mat_prop = vect_read(1)

          !Define material properties and polynomial degree element wise
          PolyMesh%Elem_loc(ie)%Degree   = PolyData%sdeg_mat(vect_read(1))
          PolyMesh%Elem_loc(ie)%NDof_loc = (PolyMesh%Elem_loc(ie)%Degree + 1) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 2) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 3) / 6;
!          PolyMesh%Elem_loc(ie)%Rho    = PolyData%prop_mat(vect_read(1),1)
!          PolyMesh%Elem_loc(ie)%Lambda = PolyData%prop_mat(vect_read(1),2)
!          PolyMesh%Elem_loc(ie)%Mu     = PolyData%prop_mat(vect_read(1),3)
!          PolyData%Elem_loc(ie)%QS     = PolyData%QS(vect_read(1))
!          PolyData%Elem_loc(ie)%QP     = PolyData%QP(vect_read(1))

          PolyMesh%Elem_loc(ie)%num_faces = 6
          PolyMesh%Elem_loc(ie)%num_vert = 8

          allocate(PolyMesh%Elem_loc(ie)%vert(8))
          PolyMesh%Elem_loc(ie)%vert(1:8) = vect_read(2:9)

          allocate(PolyMesh%Elem_loc(ie)%faces(6,4))
          PolyMesh%Elem_loc(ie)%faces(1,1:4) = &
!               [vect_read(2), vect_read(3), vect_read(4), vect_read(5)]
               [vect_read(5), vect_read(4), vect_read(3), vect_read(2)]
          PolyMesh%Elem_loc(ie)%faces(2,1:4) = &
               [vect_read(6), vect_read(7), vect_read(8), vect_read(9)]
          PolyMesh%Elem_loc(ie)%faces(3,1:4) = &
               [vect_read(3), vect_read(4), vect_read(8), vect_read(7)]
          PolyMesh%Elem_loc(ie)%faces(4,1:4) = &
               [vect_read(4), vect_read(5), vect_read(9), vect_read(8)]
          PolyMesh%Elem_loc(ie)%faces(5,1:4) = &
               [vect_read(5), vect_read(2), vect_read(6), vect_read(9)]
          PolyMesh%Elem_loc(ie)%faces(6,1:4) = &
               [vect_read(2), vect_read(3), vect_read(7), vect_read(6)]

          allocate(PolyMesh%Elem_loc(ie)%neigh_el(6,0:4))
          PolyMesh%Elem_loc(ie)%neigh_el = 0
          !write(*,*) ie, PolyMesh%Elem_loc(ie)%neigh_el
          !read(*,*)
          !kiter = kiter + 9
      enddo
      deallocate(vect_read)
      close(unit_mpi)

      !then tet
      open(unit_mpi,file=mpi_file_tet)
      read(unit_mpi,*) num_elem_list
      allocate(vect_read(5))

      num_tria_loc  = num_elem_list/5
      kiter = 1
      do ie = 1, num_elem_list/5
         !do i = kiter + 1, kiter + 5
            read(unit_mpi,*) vect_read
         ! enddo
          ie_shift = ie + num_quad_loc
          PolyMesh%Elem_loc(ie_shift)%el_type  = 'TET'
          PolyMesh%Elem_loc(ie_shift)%mat_prop = vect_read(1)

          PolyMesh%Elem_loc(ie)%Degree   = PolyData%sdeg_mat(vect_read(1))
          PolyMesh%Elem_loc(ie)%NDof_loc = (PolyMesh%Elem_loc(ie)%Degree + 1) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 2) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 3) / 6;

          PolyMesh%Elem_loc(ie_shift)%num_faces = 4
          PolyMesh%Elem_loc(ie_shift)%num_vert = 4

          allocate(PolyMesh%Elem_loc(ie_shift)%vert(4))
          PolyMesh%Elem_loc(ie_shift)%vert(1:4) = vect_read(2:5)

          allocate(PolyMesh%Elem_loc(ie_shift)%faces(4,3))
          PolyMesh%Elem_loc(ie_shift)%faces(1,1:3) = &
!               [vect_read(2), vect_read(3), vect_read(4)]
               [vect_read(4), vect_read(3), vect_read(2)]
          PolyMesh%Elem_loc(ie_shift)%faces(2,1:3) = &
               [vect_read(2), vect_read(3), vect_read(5)]
          PolyMesh%Elem_loc(ie_shift)%faces(3,1:3) = &
               [vect_read(3), vect_read(4), vect_read(5)]
          PolyMesh%Elem_loc(ie_shift)%faces(4,1:3) = &
               [vect_read(4), vect_read(2), vect_read(5)]

          allocate(PolyMesh%Elem_loc(ie_shift)%neigh_el(4,0:5))
          !allocate(PolyMesh%Elem_loc(ie_shift)%flag(4))
          PolyMesh%Elem_loc(ie_shift)%neigh_el = 0
          !PolyMesh%Elem_loc(ie_shift)%flag=-1

          !PolyMesh%Poly(i)=PolyMesh%elem_in_poly()
          !PolyMesh%Poly2glo(kiter)=PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(ie))

          !kiter = kiter + 5
      enddo
      deallocate(vect_read)
      close(unit_mpi)


      !finally prysm
      open(unit_mpi,file=mpi_file_pry)
      read(unit_mpi,*) num_elem_list
      allocate(vect_read(6))


      do ie = 1, num_elem_list/6
         !do i = kiter + 1, kiter + 6
            read(unit_mpi,*) vect_read
         !enddo
          ie_shift = ie+num_quad_loc+num_tria_loc
          PolyMesh%Elem_loc(ie_shift)%el_type  = 'PRY'
          PolyMesh%Elem_loc(ie_shift)%mat_prop = vect_read(1)

          PolyMesh%Elem_loc(ie)%Degree   = PolyData%sdeg_mat(vect_read(1))
          PolyMesh%Elem_loc(ie)%NDof_loc = (PolyMesh%Elem_loc(ie)%Degree + 1) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 2) * &
                                           (PolyMesh%Elem_loc(ie)%Degree + 3) / 6;

          PolyMesh%Elem_loc(ie_shift)%num_faces = 5
          PolyMesh%Elem_loc(ie_shift)%num_vert = 5

          allocate(PolyMesh%Elem_loc(ie_shift)%vert(5))
          PolyMesh%Elem_loc(ie_shift)%vert(1:5) = vect_read(2:6)

          allocate(PolyMesh%Elem_loc(ie_shift)%faces(5,4))
          PolyMesh%Elem_loc(ie_shift)%faces(1,1:4) = &
               [vect_read(2), vect_read(3), vect_read(4), vect_read(5)]
!               [vect_read(5), vect_read(4), vect_read(3), vect_read(2)]
          PolyMesh%Elem_loc(ie_shift)%faces(2,1:4) = &
!               [vect_read(2), vect_read(3), vect_read(6), 0]
               [vect_read(6), vect_read(3), vect_read(2), 0]
          PolyMesh%Elem_loc(ie_shift)%faces(3,1:4) = &
!               [vect_read(3), vect_read(4), vect_read(6), 0]
               [vect_read(6), vect_read(4), vect_read(3), 0]
          PolyMesh%Elem_loc(ie_shift)%faces(4,1:4) = &
!               [vect_read(4), vect_read(5), vect_read(6), 0]
               [vect_read(6), vect_read(5), vect_read(4), 0]
          PolyMesh%Elem_loc(ie_shift)%faces(5,1:4) = &
!               [vect_read(5), vect_read(2), vect_read(6), 0]
               [vect_read(6), vect_read(2), vect_read(5), 0]

          allocate(PolyMesh%Elem_loc(ie_shift)%neigh_el(5,0:4))
          PolyMesh%Elem_loc(ie_shift)%neigh_el = 0

          !kiter = kiter + 6
      enddo
      deallocate(vect_read)
      close(unit_mpi)

      unit_mpi=400;

      if (PolyMesh%num_poly/=PolyMesh%num_elem) then
         open(unit_mpi,file = u_name)
         read(unit_mpi,*) num_elem_list
         print *,num_elem_list
         allocate(vect_read(2))
         do i = 1, num_elem_list
            read(unit_mpi,*) vect_read
            PolyMesh%elem_in_poly_loc(i)=vect_read(2)
         enddo
         close(unit_mpi)
         deallocate(vect_read)
         call MESH_CORRECTION(PolyMesh,mpi_np,mpi_id)
      endif
      end subroutine CREATE_LOCAL_MESH

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> stores the coordinates of the vertices of the tetrahedra locally;
!>
      subroutine CREATE_VERT_LIST(gridfile, PolyMesh, mpi_id)

      use Poly_mesh
      use qsort

      implicit none

      character(len=70), intent(in) :: gridfile
      character(len=100)            :: inline

      integer(kind=4), intent(in) :: mpi_id
      integer(kind=4) :: ie, ivert, i,k, kiter,num_node, elem_total, id_node
      integer(kind=4) :: num_vert_with_duplicate!, num_local_vert
      integer(kind=4), dimension(:), allocatable :: vert_with_duplicate

      real(kind=8) :: xx,yy,zz

      type(Mesh_Structure), intent(inout) :: PolyMesh

      ! calculating number of local vertices with duplicates
      num_vert_with_duplicate = 0;
      do ie = 1, PolyMesh%num_elem_loc
         num_vert_with_duplicate =  num_vert_with_duplicate + &
                                           PolyMesh%Elem_loc(ie)%num_vert
      enddo

      ! storing the number of vertices with duplicates
      allocate(vert_with_duplicate(num_vert_with_duplicate))
      i = 1;
      do ie = 1, PolyMesh%num_elem_loc
         do ivert = 1, PolyMesh%Elem_loc(ie)%num_vert
           vert_with_duplicate(i) = PolyMesh%Elem_loc(ie)%vert(ivert);
           i = i + 1;
         enddo
      enddo

      !write(*,*) 'vert with duplicates'
      !write(*,*) vert_with_duplicate
      !write(*,*) 'end vert with duplicates'
      !read(*,*)

      !sorting the array and remove duplicates
      call QsortC(vert_with_duplicate)

      !write(*,*) 'sorting vert with duplicates'
      !write(*,*) vert_with_duplicate
      !write(*,*) 'end sorting vert with duplicates'
      !read(*,*)

      !counting the vertex without duplicates
      PolyMesh%num_node_loc = 1
      do i = 2, num_vert_with_duplicate
         if (vert_with_duplicate(i) /= vert_with_duplicate(i-1)) &
             PolyMesh%num_node_loc = PolyMesh%num_node_loc + 1;
      enddo
      !write(*,*) PolyMesh%num_node_loc
      !read(*,*)

      !computing global to local map
      allocate(PolyMesh%node_loc2glo(PolyMesh%num_node_loc));
      PolyMesh%node_loc2glo(1) = vert_with_duplicate(1);
      k = 2;
      do i = 2, num_vert_with_duplicate
         if (vert_with_duplicate(i) /= vert_with_duplicate(i-1)) then
             PolyMesh%node_loc2glo(k) = vert_with_duplicate(i);
             k = k + 1;
         endif
      enddo

      !write(*,*) 'local numeration'
      !write(*,*) PolyMesh%node_loc2glo
      !write(*,*) 'end local numeration'
      !read(*,*)

      allocate(PolyMesh%coord_x(PolyMesh%num_node_loc), &
               PolyMesh%coord_y(PolyMesh%num_node_loc), &
               PolyMesh%coord_z(PolyMesh%num_node_loc))

      open(40,file=gridfile)

      do
        read(40,'(A)') inline
        if (inline(1:1) /= '#') exit
      enddo

      read(inline,*) num_node, elem_total

      kiter = 1;
      do i = 1, num_node

        read(40,*) id_node, xx, yy, zz
        !print *,kiter,id_node, PolyMesh%node_loc2glo(kiter)
        if (id_node == PolyMesh%node_loc2glo(kiter)) then
            PolyMesh%coord_x(kiter) = xx;
            PolyMesh%coord_y(kiter) = yy;
            PolyMesh%coord_z(kiter) = zz;
            kiter = kiter + 1;
        endif
        if (kiter > PolyMesh%num_node_loc) exit;

      enddo

      !write(*,*) 'grid nodes'
      !do i = 1, PolyMesh%num_node_loc
      !   write(*,*) PolyMesh%coord_x(i), PolyMesh%coord_y(i), PolyMesh%coord_z(i)
      !enddo
      !write(*,*) 'end grid nodes'
      !read(*,*)


      end subroutine CREATE_VERT_LIST

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> stores the map poly_loc2glo to transition from local polyhedra to the global ones and initializes the other fields of the structure Polyhedron;
!>
   subroutine CREATE_POLY_LIST(PolyMesh, mpi_id)

      use Poly_mesh
      use qsort
      use find_poly

      implicit none

      integer(kind=4), intent(in) :: mpi_id
      integer(kind=4) :: ie_loc, i,j,k, kiter, ipoly,num_tet_in_poly,num_faces_in_poly
      integer(kind=4) :: num_poly_with_duplicate!, num_local_poly
      integer(kind=4), dimension(:), allocatable :: poly_with_duplicate,index_glob

      type(Mesh_Structure), intent(inout) :: PolyMesh

      ! calculating number of local vertices with duplicates
      num_poly_with_duplicate = 0;

      do ie_loc = 1, PolyMesh%num_elem_loc
         num_poly_with_duplicate =  num_poly_with_duplicate + 1
      enddo

      !print *,num_poly_with_duplicate

      ! storing the number of vertices with duplicates
      allocate(poly_with_duplicate(num_poly_with_duplicate))

      do ie_loc = 1, PolyMesh%num_elem_loc
         !ie_glob=PolyMesh%elem_loc2glo(ie_loc)
         ipoly=PolyMesh%elem_in_poly_loc(ie_loc)
         !print *,'tet',ie_glob,'in poly',ipoly
         poly_with_duplicate(ie_loc) = ipoly;
         !print *,mpi_id,'index',ie_loc,'poly',poly_with_duplicate(ie_loc)
      enddo

      !write(*,*) 'poly with duplicates'
      !write(*,*) poly_with_duplicate
      !write(*,*) 'end poly with duplicates'
      !read(*,*)

      !sorting the array and remove duplicates
      if (PolyMesh%num_poly /= PolyMesh%num_elem) call QsortC(poly_with_duplicate)

      !write(*,*) 'sorting poly with duplicates'
      !write(*,*) poly_with_duplicate
      !write(*,*) 'end sorting poly with duplicates'
      !read(*,*)


      !counting the vertex without duplicates
      PolyMesh%num_poly_loc = 1
      do i = 2, num_poly_with_duplicate
         if (poly_with_duplicate(i) /= poly_with_duplicate(i-1)) &
               PolyMesh%num_poly_loc = PolyMesh%num_poly_loc + 1;
      enddo

      !write(*,*) PolyMesh%num_poly_loc
      !read(*,*)

      !computing global to local map
      allocate(PolyMesh%poly_loc2glo(PolyMesh%num_poly_loc));
      PolyMesh%poly_loc2glo(1) = poly_with_duplicate(1);
      k = 2;
      do i = 2, num_poly_with_duplicate
         if (poly_with_duplicate(i) /= poly_with_duplicate(i-1)) then
               PolyMesh%poly_loc2glo(k) = poly_with_duplicate(i);
               !print *,k
               k = k + 1;
         endif
      enddo
      !poly_loc2glo

      !write(*,*) 'local numeration'
      !write(*,*) PolyMesh%poly_loc2glo

      !write(*,*) PolyMesh%elem_in_poly_loc
      !write(*,*) 'end local numeration'
      !read(*,*)

      allocate(PolyMesh%Poly(PolyMesh%num_poly_loc))

      kiter = 1;
      do ipoly = 1, PolyMesh%num_poly
         !if (ipoly==21) print *,'bbbb'
         if (ipoly == PolyMesh%poly_loc2glo(kiter)) then

            !print *,kiter, PolyMesh%poly_loc2glo(kiter)
            !if (ipoly==21) print *,'aaaaaaaa'
            index_glob=FIND_TET_IN_POLY(PolyMesh%elem_in_poly_loc,ipoly,PolyMesh%num_elem_loc)
            num_tet_in_poly=size(index_glob)
            !if (size(index_glob)==1 .and. index_glob(1)==0) print *,'AAA'
            num_faces_in_poly=num_tet_in_poly*4;
            !print *,num_tet_in_poly
            PolyMesh%Poly(kiter)%num_tet_in_poly=num_tet_in_poly
            allocate(PolyMesh%Poly(kiter)%tet_in_poly(num_tet_in_poly))
            allocate(PolyMesh%Poly(kiter)%neigh_bbox(num_faces_in_poly,3,2))
            allocate(PolyMesh%Poly(kiter)%neigh_hk(num_faces_in_poly))

            PolyMesh%Poly(kiter)%hk=0.0
            PolyMesh%Poly(kiter)%b_box(1,:)=[0.0,0.0]
            PolyMesh%Poly(kiter)%b_box(2,:)=[0.0,0.0]
            PolyMesh%Poly(kiter)%b_box(3,:)=[0.0,0.0]

            do i=1,num_faces_in_poly
               PolyMesh%Poly(kiter)%neigh_bbox(i,1,:) = [0.0,0.0]
               PolyMesh%Poly(kiter)%neigh_bbox(i,2,:) = [0.0,0.0]
               PolyMesh%Poly(kiter)%neigh_bbox(i,3,:) = [0.0,0.0]
               PolyMesh%Poly(kiter)%neigh_hk(i)=0.0
            enddo

            do j=1,num_tet_in_poly
               PolyMesh%Poly(kiter)%tet_in_poly(j)=PolyMesh%elem_loc2glo(index_glob(j))
               !print *,'tet: ',PolyMesh%Poly(kiter)%tet_in_poly(j)
            enddo

            kiter = kiter + 1;
         endif
         if (kiter > PolyMesh%num_poly_loc) exit;
      enddo

      end subroutine CREATE_POLY_LIST

   ! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> computes the coordinates of the normal vector for each face and determines the area associated with each face;
!>
subroutine CREATE_NORMAL_FACE(PolyMesh, mpi_id)

   use Poly_mesh

   implicit none

   integer(kind=4), intent(in) :: mpi_id
   integer(kind=4) :: ie, iface, id_node1, id_node2, id_node3

   integer(kind=4), dimension(3) :: ivert
   real(kind=8) :: Px, Py, Pz, Qx, Qy, Qz, n1, n2, n3, norm_n

   type(Mesh_Structure), intent(inout) :: PolyMesh

   do ie = 1, PolyMesh%num_elem_loc

      allocate(PolyMesh%Elem_loc(ie)%normal(PolyMesh%Elem_loc(ie)%num_faces,3))
      allocate(PolyMesh%Elem_loc(ie)%area(PolyMesh%Elem_loc(ie)%num_faces))

      do iface = 1, PolyMesh%Elem_loc(ie)%num_faces

         ivert = PolyMesh%Elem_loc(ie)%faces(iface,1:3)

         call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                 ivert(1),id_node1);
         call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                 ivert(2),id_node2);
         call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                 ivert(3),id_node3);


         Qx = PolyMesh%coord_x(id_node1) - PolyMesh%coord_x(id_node2);
         Qy = PolyMesh%coord_y(id_node1) - PolyMesh%coord_y(id_node2);
         Qz = PolyMesh%coord_z(id_node1) - PolyMesh%coord_z(id_node2);

         Px = PolyMesh%coord_x(id_node3) - PolyMesh%coord_x(id_node2);
         Py = PolyMesh%coord_y(id_node3) - PolyMesh%coord_y(id_node2);
         Pz = PolyMesh%coord_z(id_node3) - PolyMesh%coord_z(id_node2);

         !write(*,*) Px, Py, Pz
         !write(*,*) Qx, Qy, Qz

         n1 = Py*Qz - Pz*Qy;
         n2 = - Px*Qz + Pz*Qx;
         n3 = Px*Qy - Py*Qx;

         !write(*,*) n1, n2, n3
         !read(*,*)

         norm_n = dsqrt(n1**2+n2**2+n3**2)

         PolyMesh%Elem_loc(ie)%normal(iface,1) = n1/norm_n
         PolyMesh%Elem_loc(ie)%normal(iface,2) = n2/norm_n
         PolyMesh%Elem_loc(ie)%normal(iface,3) = n3/norm_n
         PolyMesh%Elem_loc(ie)%area(iface) = norm_n

      enddo
   enddo

end subroutine CREATE_NORMAL_FACE

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> performs the computation of the coordinates of the bounding box for each polyhedron;
!>
   subroutine CREATE_BBOX_EL(PolyMesh)

      use Poly_mesh
      use local_search

      implicit none

      integer(kind=4) :: ipoly_loc,ipoly_glob, ivert, id_node
      integer (kind=4) :: t,l,ie,ie_loc,ie_glob
      integer(kind=4) :: num_tet_in_poly,num_vert_poly
      real(kind=8) :: dist

      real(kind=8), dimension(:), allocatable :: xx_vert, yy_vert, zz_vert
      !real(kind=8), dimension(:), allocatable :: xx_send, yy_send, zz_send
      !real(kind=8), dimension(:), allocatable :: xx_recv, yy_recv, zz_recv

      type(Mesh_Structure), intent(inout) :: PolyMesh

      do ipoly_loc = 1, PolyMesh%num_poly_loc

         ipoly_glob=PolyMesh%poly_loc2glo(ipoly_loc)

         num_tet_in_poly=PolyMesh%Poly(ipoly_loc)%num_tet_in_poly;

         num_vert_poly=num_tet_in_poly*PolyMesh%Elem_loc(1)%num_vert ! number of vertices of the polyhedron ipoly_loc

         allocate(xx_vert(num_vert_poly), &
               yy_vert(num_vert_poly), &
               zz_vert(num_vert_poly))
         t=1;

         do ie=1,num_tet_in_poly

            ie_glob=PolyMesh%Poly(ipoly_loc)%tet_in_poly(ie)

            call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                       PolyMesh%num_elem_loc, &
                                       ie_glob,ie_loc)

            do ivert = 1, PolyMesh%Elem_loc(ie_loc)%num_vert

               call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                     PolyMesh%Elem_loc(ie_loc)%vert(ivert),id_node);

               xx_vert(t) = PolyMesh%coord_x(id_node)
               yy_vert(t) = PolyMesh%coord_y(id_node)
               zz_vert(t) = PolyMesh%coord_z(id_node)
               t=t+1;
            end do
         enddo

         PolyMesh%Poly(ipoly_loc)%b_box(1,1) = minval(xx_vert);
         PolyMesh%Poly(ipoly_loc)%b_box(1,2) = maxval(xx_vert);
         PolyMesh%Poly(ipoly_loc)%b_box(2,1) = minval(yy_vert);
         PolyMesh%Poly(ipoly_loc)%b_box(2,2) = maxval(yy_vert);
         PolyMesh%Poly(ipoly_loc)%b_box(3,1) = minval(zz_vert);
         PolyMesh%Poly(ipoly_loc)%b_box(3,2) = maxval(zz_vert);

         PolyMesh%Poly(ipoly_loc)%hk=0.0
         do  t = 1,num_vert_poly
            do l= t+1,num_vert_poly
               dist = SQRT((xx_vert(t)-xx_vert(l))*(xx_vert(t)-xx_vert(l))+(yy_vert(t)-yy_vert(l))*(yy_vert(t)-yy_vert(l))+(zz_vert(t)-zz_vert(l))*(zz_vert(t)-zz_vert(l)))
               if (dist > PolyMesh%Poly(ipoly_loc)%hk) then
                  PolyMesh%Poly(ipoly_loc)%hk = dist
               end if
            end do
         end do
         deallocate(xx_vert, yy_vert, zz_vert)
      enddo

      end subroutine CREATE_BBOX_EL

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> ???
      !> 
      subroutine FIND_TET_IN_PROC(n_tet,tet_in_poly,n_el,v,IsFound)

         use local_search

         implicit none

         integer(kind=4), intent(in out) :: n_tet
         integer(kind=4), intent(in out) :: n_el
         integer(kind=4), intent(in out), dimension(n_el) :: v
         integer(kind=4),intent(in out),dimension(n_tet) :: tet_in_poly
         integer(kind=4) :: i,ie_glob,ie_loc

         logical :: IsFound

         IsFound = .true.

         do i=1,n_tet
            ie_glob=tet_in_poly(i)
            call GET_EL_LOC_FROM_EL_GLO(v, n_el, &
                                       ie_glob,ie_loc)
            if (ie_loc==0) then
               isFound = .false.
            endif
         enddo

      end subroutine FIND_TET_IN_PROC

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> ???
      !>
      subroutine CREATE_NEIGH_EL_QUAD(mpifile, PolyMesh, mpi_np, mpi_id)

      use mpi
      use Poly_setup_MPI, only: mpi_ierr
      use Poly_mesh
      use qsort
      use local_search
      use Poly_exit_codes, only: EXIT_NEIGHBOUR_EL_ERROR


      implicit none

      character(len=70), intent(in) :: mpifile
      character(len=70) :: mpi_file_qua
      integer(kind=4), intent(in)  :: mpi_np, mpi_id
      integer(kind=4)  :: unit_mpi, num_elem_list, num_quad_loc, i, Row2, &
                          num_quad_send, ie, iface, mat, ie_ne, iface_ne, mat_ne, &
                          ie_loc, ie_ne_loc, ip, num_quad_send_mpi, kiter, &
                          num_quad_send_loc

      integer(kind=4), dimension(:,:), allocatable :: con_quad_loc, con_quad_recv_mpi
      integer(kind=4), dimension(:),   allocatable :: con_quad_send, &
                                                      con_quad_send_mpi


      !integer(kind=4) :: vect_sort_q(4), vect_sort_t(3)
      logical :: IsFound_int, IsFound_bnd, IsQuad

      type(Mesh_Structure), intent(inout) :: PolyMesh


      !1 - Reoder con_quad structure
      do ie = 1, PolyMesh%num_quad

         call QsortC(PolyMesh%con_quad(ie,2:5))

      enddo


      !2 - Read con_qua_000000.mpi and load faces
      mpi_file_qua = 'con_qua_000000.mpi'

      unit_mpi = 40 + mpi_id
      if (mpi_id < 10) then
         write(mpi_file_qua(14:14),'(i1)') mpi_id
      elseif (mpi_id < 100) then
         write(mpi_file_qua(13:14),'(i2)') mpi_id
      elseif (mpi_id < 1000) then
         write(mpi_file_qua(12:14),'(i3)') mpi_id
      elseif (mpi_id < 10000) then
         write(mpi_file_qua(11:14),'(i4)') mpi_id
      elseif (mpi_id < 100000) then
         write(mpi_file_qua(10:14),'(i5)') mpi_id
      elseif (mpi_id < 1000000) then
         write(mpi_file_qua(9:14),'(i6)') mpi_id
      endif

      mpi_file_qua = mpifile(1:len_trim(mpifile)) // '/' // mpi_file_qua

      open(unit_mpi,file=mpi_file_qua)
      read(unit_mpi,*) num_elem_list

      num_quad_loc =  num_elem_list/7
      num_quad_send = 0

      allocate(con_quad_loc(num_quad_loc,7))

      do i = 1, num_quad_loc
          read(unit_mpi,*) con_quad_loc(i,1:7)
      enddo

     !3 - Find neighbouring elements between elements in the same processor
     ! first internal faces, than boundary faces
     IsQuad = .true.
      do i = 1, num_quad_loc
         isFound_int = .false.
         isFound_bnd = .false.
         !if (mpi_id == 2) write(*,*) 'before', con_quad_loc(i,:)

            if (con_quad_loc(i,1) /= 0) then
                call FIND_NEIGHBOUR_EL(i,con_quad_loc,num_quad_loc,&
                                       IsFound_int,Row2,IsQuad)
                if(IsFound_int) then
                   mat   = con_quad_loc(i,1); mat_ne   = con_quad_loc(Row2,1)
                   ie    = con_quad_loc(i,2); ie_ne    = con_quad_loc(Row2,2)
                   iface = con_quad_loc(i,3); iface_ne = con_quad_loc(Row2,3)

                   call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                               PolyMesh%num_elem_loc, &
                                               ie,ie_loc)

                   call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                               PolyMesh%num_elem_loc, &
                                               ie_ne,ie_ne_loc)


                   PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = mpi_id
                   PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:3) = [mat_ne, ie_ne, iface_ne]

                   PolyMesh%Elem_loc(ie_ne_loc)%neigh_el(iface_ne,0) = mpi_id
                   PolyMesh%Elem_loc(ie_ne_loc)%neigh_el(iface_ne,1:3) = [mat, ie, iface]

                   con_quad_loc(i,:) = 0;
                   con_quad_loc(Row2,:) = 0;
                else

                   call FIND_BOUNDARY_EL(i,con_quad_loc,num_quad_loc,&
                                         PolyMesh%con_quad, PolyMesh%num_quad, &
                                         IsFound_bnd,Row2,IsQuad)

                   if(IsFound_bnd) then
                      mat   = con_quad_loc(i,1); mat_ne   = - PolyMesh%con_quad(Row2,1)
                      ie    = con_quad_loc(i,2); ie_ne    =   con_quad_loc(i,2)
                      iface = con_quad_loc(i,3); iface_ne =   con_quad_loc(i,3);

                      call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                                  PolyMesh%num_elem_loc, &
                                                  ie,ie_loc)

                      PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = mpi_id
                      PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:3) = &
                                                             [mat_ne, ie_ne, iface_ne]
                      con_quad_loc(i,:) = 0;
                      PolyMesh%con_quad(Row2,:) = 0;
                   endif
                endif
           endif

           !if (mpi_id == 2) write(*,*) 'after', IsFound_int, IsFound_bnd, con_quad_loc(i,:)

           if((IsFound_int .eqv. .false.) .and. (IsFound_bnd .eqv. .false.) &
                .and. con_quad_loc(i,4) /= 0) &
               num_quad_send =  num_quad_send + 1 !num_quad left to search
      enddo

      num_quad_send_loc = num_quad_send

 !      write(*,*) 'id', mpi_id, num_quad_send
 !    call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
!      if(mpi_id == 1) then
!        do ie = 1, PolyMesh%num_elem_loc
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(1,:)
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(2,:)
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(3,:)
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(4,:)
!           if(PolyMesh%Elem_loc(ie)%el_type == 'PRY' .or. &
!              PolyMesh%Elem_loc(ie)%el_type == 'HEX') &
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(5,:)
!           if(PolyMesh%Elem_loc(ie)%el_type == 'HEX') &
!           write(*,*) ie, 'el', PolyMesh%Elem_loc(ie)%neigh_el(6,:)
!           write(*,*) '--------------------------------'
!           !read(*,*)
!        enddo
!      endif


!      if(mpi_id == 2) then
!        do i = 1, num_quad_loc
!           write(*,*) 'con_quad', con_quad_loc(i,:)
!        enddo
!      endif

     !if(mpi_id == 2) write(*,*) '============================'

     call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

     !4 - Preparation for sending quadrilateral faces
     allocate(con_quad_send(7*num_quad_send))
     kiter = 0
     do i = 1, num_quad_loc

        if (con_quad_loc(i,1) /= 0) then
            con_quad_send(kiter+1:kiter+7) = con_quad_loc(i,1:7)
!            if(mpi_id == 0) write(*,*) 'con_quad', con_quad_loc(i,:)
            kiter=kiter+7
         endif
     enddo

!     if(mpi_id == 1) write(*,*) '======================'


     call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)


     !5 - Broadcasting elements
     do ip = 1, mpi_np

       if(mpi_id == ip-1) num_quad_send_mpi = 7*num_quad_send;
!       if(mpi_id == ip-1) write(*,*)  num_quad_send_mpi

       call MPI_BCAST(num_quad_send_mpi, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)


       !write(*,*) 'id', mpi_id, 'ip', ip-1, num_quad_send_mpi
       call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

       allocate(con_quad_send_mpi(num_quad_send_mpi))
       if(mpi_id == ip-1) con_quad_send_mpi = con_quad_send;

!       if(mpi_id == ip-1) write(*,*) mpi_id, ip-1,'con', con_quad_send_mpi

       call MPI_BCAST(con_quad_send_mpi,num_quad_send_mpi,&
                      MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      !unpaking the elements

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

!      if(mpi_id == ip-1) write(*,*) '********************************'

      allocate(con_quad_recv_mpi(num_quad_send_mpi/7,7))
      kiter = 1
      do i = 1, num_quad_send_mpi, 7
        con_quad_recv_mpi(kiter,1:7) = con_quad_send_mpi(i:i+6)

!        if(mpi_id == 1) write(*,*) num_quad_send_mpi, i,i+6,con_quad_send_mpi(i:i+6)
!        if(mpi_id == 1) write(*,*) 'recv_from', ip-1, 'con_quad_rcv', con_quad_recv_mpi(kiter,1:7)

        kiter = kiter + 1
      enddo

!      if(mpi_id == 1) write(*,*) '-------------------------------------------'

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      !6 - find neighbouring elements located in differen processors
      do i = 1, num_quad_loc
         isFound_int = .false.

         if (con_quad_loc(i,1) /= 0 .and. mpi_id /= ip-1) then

           ! if(mpi_id == 0) then
           !  do j = 1,  num_quad_send_mpi/7
           !    write(*,*) con_quad_recv_mpi(j,:)
           !  enddo
           ! endif

            call FIND_NEIGHBOURING_EL_MPI(i,con_quad_loc, num_quad_loc,&
                                            con_quad_recv_mpi, num_quad_send_mpi/7, &
                                            IsFound_int,Row2, IsQuad)

            if(IsFound_int) then
              mat   = con_quad_loc(i,1); mat_ne   = con_quad_recv_mpi(Row2,1)
              ie    = con_quad_loc(i,2); ie_ne    = con_quad_recv_mpi(Row2,2)
              iface = con_quad_loc(i,3); iface_ne = con_quad_recv_mpi(Row2,3);

              call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                          PolyMesh%num_elem_loc, &
                                          ie,ie_loc)

              PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = ip - 1;
              PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:3) = [mat_ne, ie_ne, iface_ne]

              !put zero if found
              con_quad_loc(i,:) = 0;
              con_quad_recv_mpi(Row2,:) = 0
              num_quad_send_loc = num_quad_send_loc - 1

            endif
          endif
       enddo

      !go to next processor
      deallocate(con_quad_send_mpi, con_quad_recv_mpi)
      enddo

!      if(mpi_id==0) then
!        do i = 1, num_quad_loc
!          if(con_quad_loc(i,1) /= 0)  write(*,*) 'con_quad', con_quad_loc(i,:)
!        enddo
!      endif

     call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      if(num_quad_send_loc == 0) then
           write(*,*) 'Proc ', mpi_id, ':', ' found all quad interfaces!'
      else
           write(*,*) 'Proc ', mpi_id, ':',  num_quad_send_loc, ' quad interfaces not found!'
           !call EXIT(EXIT_NEIGHBOUR_EL_ERROR)
      endif

      end subroutine CREATE_NEIGH_EL_QUAD


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> stores the information of the neighbouring tetrahedra and polyhedra by reading data from con_tri.mpi;
!>
subroutine CREATE_NEIGH_EL_TRIA(mpifile, PolyMesh, mpi_np, mpi_id)

   use mpi
   use Poly_setup_MPI, only: mpi_ierr
   use Poly_mesh
   use qsort
   use local_search
   use Poly_exit_codes, only: EXIT_NEIGHBOUR_EL_ERROR

   implicit none

   character(len=70), intent(in) :: mpifile
   character(len=70) :: mpi_file_tri
   integer(kind=4), intent(in)  :: mpi_np, mpi_id
   !integer(kind=4) :: status(MPI_STATUS_SIZE)
   integer(kind=4) :: unit_mpi, num_elem_list, num_tria_loc, i, j,Row2, &
                        num_tria_send, ie, iface, mat, ie_ne, iface_ne, mat_ne, &
                        ie_loc, ipoly_glob,ipoly2_glob,ipoly_loc,ipoly2_loc,ie_ne_loc, ip, &
                        num_tria_send_mpi, kiter,kiter2, num_tria_send_loc,&
                        iface_poly,num_tet_in_poly
   !integer (kind=4), dimension(:), allocatable :: faces_found_loc,faces_to_find_send
   !integer (kind=4) :: num_faces_send_mpi,iter_face,iter_face_2

   integer(kind=4), dimension(:,:), allocatable :: con_tria_loc, con_tria_recv_mpi
   integer(kind=4), dimension(:),   allocatable :: con_tria_send, &
                                                   con_tria_send_mpi
   !integer(kind=4), dimension(:), allocatable :: faces_to_find_send_mpi
   !integer(kind=4), dimension(1) :: one = 1

   real(kind=8), dimension(:,:), allocatable :: xx_loc
   real(kind=8), dimension(:,:), allocatable :: yy_loc
   real(kind=8), dimension(:,:), allocatable :: zz_loc
   real(kind=8), dimension(:), allocatable :: hk_send,hk_send_mpi
   real(kind=8), dimension(:), allocatable :: xx_send
   real(kind=8), dimension(:), allocatable :: yy_send
   real(kind=8), dimension(:), allocatable :: zz_send
   real(kind=8), dimension(:), allocatable :: x1_send_mpi,x2_send_mpi
   real(kind=8), dimension(:), allocatable :: y1_send_mpi,y2_send_mpi
   real(kind=8), dimension(:), allocatable :: z1_send_mpi,z2_send_mpi
   !real(kind=8), dimension(:,:), allocatable :: xx_recv_mpi,yy_recv_mpi,zz_recv_mpi

   integer(kind=4), dimension(:), allocatable :: index_send,index_send_mpi
   integer(kind=4) :: num_poly_send,num_poly_send_loc,num_poly_send_mpi,num_index_send_mpi,num_hk_send_mpi
   integer(kind=4) :: space_fun_tag

   logical :: IsFound_int, IsFound_bnd, IsQuad

   type(Mesh_Structure), intent(inout) :: PolyMesh

   ! ----------
   ! 1 - Reoder con_quad structure
   ! ----------

   do ie = 1, PolyMesh%num_tria

      call QsortC(PolyMesh%con_tria(ie,2:4))

   enddo
   ! ----------
   ! 2 - read con_tria_000000.mpi and load faces
   ! ----------

   mpi_file_tri = 'con_tri_000000.mpi'

   unit_mpi = 40 + mpi_id
   if (mpi_id < 10) then
      write(mpi_file_tri(14:14),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(mpi_file_tri(13:14),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(mpi_file_tri(12:14),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(mpi_file_tri(11:14),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(mpi_file_tri(10:14),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(mpi_file_tri(9:14),'(i6)') mpi_id
   endif

   mpi_file_tri = mpifile(1:len_trim(mpifile)) // '/' // mpi_file_tri

   open(unit_mpi,file=mpi_file_tri)
   read(unit_mpi,*) num_elem_list

   num_tria_loc =  num_elem_list/6
   num_tria_send = 0
   num_poly_send = 0

   !write(*,*) num_tria_loc

   allocate(con_tria_loc(num_tria_loc,6))

   allocate(xx_loc(PolyMesh%num_poly_loc,2))
   allocate(yy_loc(PolyMesh%num_poly_loc,2))
   allocate(zz_loc(PolyMesh%num_poly_loc,2))
   allocate(index_send(PolyMesh%num_poly_loc))
   allocate(hk_send(PolyMesh%num_poly_loc))

   do i = 1, num_tria_loc
         read(unit_mpi,*) con_tria_loc(i,1:6)
   enddo

   do i=1,PolyMesh%num_poly_loc
      xx_loc(i,1:2)=PolyMesh%Poly(i)%b_box(1,1:2)
      yy_loc(i,1:2)=PolyMesh%Poly(i)%b_box(2,1:2)
      zz_loc(i,1:2)=PolyMesh%Poly(i)%b_box(3,1:2)
      hk_send(i)=PolyMesh%Poly(i)%hk
      !if (index_loc(i)==2866) print *,'index_loc',i
   enddo

   ! ----------
   ! 3 - Find neighbouring elements between elements in the same processor
   ! first internal faces, than boundary faces
   ! ----------

   IsQuad = .false.
   tria_loop: do i = 1, num_tria_loc
      isFound_int = .false.
      isFound_bnd = .false.

         !if (mpi_id == 0) write(*,*) 'before', con_tria_loc(i,:)

         if (con_tria_loc(i,1) /= 0) then
            call FIND_NEIGHBOUR_EL(i, con_tria_loc, num_tria_loc,&
                                    IsFound_int, Row2, IsQuad)
            ! internal face
            if(IsFound_int) then
               mat   = con_tria_loc(i,1); mat_ne   = con_tria_loc(Row2,1)
               ie    = con_tria_loc(i,2); ie_ne    = con_tria_loc(Row2,2)
               iface = con_tria_loc(i,3); iface_ne = con_tria_loc(Row2,3)
!                   write(*,*) con_quad_loc(i,:)
!                   write(*,*) con_quad_loc(Row2,:)
!                   write(*,*) mat, ie, iface
!                   write(*,*) mat_ne, ie_ne, iface_ne
!                   read(*,*)
!                   write(*,*) ie, iface, PolyMesh%Elem_loc(ie)%neigh_el(iface,:)
!                   write(*,*) ie_ne, iface_ne, PolyMesh%Elem_loc(ie_ne)%neigh_el(iface_ne,:)
!                   read(*,*)

               call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                          PolyMesh%num_elem_loc, &
                                          ie,ie_loc)

               call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                          PolyMesh%num_elem_loc, &
                                          ie_ne,ie_ne_loc)

               ipoly_glob=PolyMesh%elem_in_poly(ie)
               ipoly2_glob=PolyMesh%elem_in_poly(ie_ne)

               PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = mpi_id
               PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:4) = [mat_ne, ie_ne, iface_ne,ipoly2_glob]

               PolyMesh%Elem_loc(ie_ne_loc)%neigh_el(iface_ne,0) = mpi_id
               PolyMesh%Elem_loc(ie_ne_loc)%neigh_el(iface_ne,1:4) = [mat, ie, iface,ipoly_glob]

               !write(*,*) 'el', PolyMesh%Elem_loc(ie)%neigh_el(iface,1:3)
               !write(*,*) 'ne', PolyMesh%Elem_loc(ie_ne)%neigh_el(iface_ne,1:3)
               !read(*,*)
               !PolyMesh%Elem_loc(ie_loc)%flag(iface)=1
               !PolyMesh%Elem_loc(ie_ne_loc)%flag(iface_ne)=0


               con_tria_loc(i,:) = 0;
               con_tria_loc(Row2,:) = 0;
               !print *,ipoly_loc
               !faces_found_loc(i)=1;
               !print *,'local internal'
               !print *,i,Row2
               !faces_found_loc(Row2)=1;

            ! boundary face
            else
               !if (mpi_id == 1 .and. i == 18) then
               call FIND_BOUNDARY_EL(i,con_tria_loc,num_tria_loc,&
                                    PolyMesh%con_tria, PolyMesh%num_tria, &
                                    IsFound_bnd,Row2,IsQuad)
               !write(*,*) IsFound_bnd
               !endif
                  if(IsFound_bnd) then
                  mat   = con_tria_loc(i,1); mat_ne   =   PolyMesh%con_tria(Row2,1)
                  space_fun_tag = PolyMesh%con_tria(Row2,5)
                  ie    = con_tria_loc(i,2);

                  !! PROBLEM
                  ! TODO - case with more faces on domain? case with different BC? make automatic
                  ! Dirichlet bc identifier
                  if(mat_ne == 2 .or. mat_ne == 3 .or. mat_ne == 4) ie_ne = -1
                  ! if(mat_ne == 2) ie_ne = -1

                  ! Neumann bc identifier
                  if(mat_ne == 5 .or. mat_ne == 6 .or. mat_ne == 7) ie_ne = -2
                  ! if(mat_ne == 3 .or. mat_ne == 4 .or. mat_ne == 5 .or. mat_ne == 6 .or. mat_ne == 7) ie_ne = -2

                  iface = con_tria_loc(i,3); iface_ne =   con_tria_loc(i,3);

                  call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                             PolyMesh%num_elem_loc, &
                                             ie,ie_loc)

                  ipoly_glob=PolyMesh%elem_in_poly(ie)

                  PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = mpi_id
                  PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:5) = &
                                                         [mat_ne, ie_ne, iface_ne, ipoly_glob, space_fun_tag]

                  ! print *, "neigh_el: ", PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,:)
                  ! print *, ""

                  !PolyMesh%Elem_loc(ie_loc)%flag(iface)=0;
                  con_tria_loc(i,:) = 0;
                  PolyMesh%con_tria(Row2,:) = 0;
                     !faces_found_loc(i)=1;
                  !  print *,'boundary'
                  !  print *,i,Row2
                     !faces_found_loc(Row2)=1;
               endif
            endif
         endif
         !if(mpi_id == 1) write(*,*) IsFound_int, IsFound_bnd, 'el', i, con_quad_loc(i,4:7)
         !read(*,*)

         if((IsFound_int .eqv. .false.) .and. (IsFound_bnd .eqv. .false.) &
               .and. con_tria_loc(i,4) /= 0) &
               num_tria_send =  num_tria_send + 1

   enddo tria_loop

   num_tria_send_loc = num_tria_send

   num_poly_send = PolyMesh%num_poly_loc

   num_poly_send_loc = num_poly_send

   !print *,num_tria_send

   !print *,'___________________'
   !print *,num_poly_send_loc

   !if(mpi_id == 1)  write(*,*) num_tria_send_loc
   !read(*,*)

   !callMPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

   ! ----------
   ! 4 - Preparation for sending triangular faces
   ! ----------

   allocate(con_tria_send(6*num_tria_send))
   !allocate(faces_to_find_send(num_tria_send))

   kiter = 0

   do i = 1, num_tria_loc
      if (con_tria_loc(i,1) /= 0) then
            con_tria_send(kiter+1:kiter+6) = con_tria_loc(i,1:6)
            kiter = kiter + 6
            !if(mpi_id == 1) write(*,*) mpi_id, 'ctria', con_tria_loc(i,1:6)
         endif
   enddo

   allocate(xx_send(2*num_poly_send),yy_send(2*num_poly_send),zz_send(2*num_poly_send))

   kiter2 = 0

   do i=1,PolyMesh%num_poly_loc

      xx_send(kiter2+1:kiter2+2) = xx_loc(i,1:2)
      yy_send(kiter2+1:kiter2+2)=yy_loc(i,1:2)
      zz_send(kiter2+1:kiter2+2)=zz_loc(i,1:2)
      !print *,'proc',mpi_id,'i',i,xx_send(kiter2+1:kiter2+2),yy_send(kiter2+1:kiter2+2),zz_send(kiter2+1:kiter2+2)
      index_send(i)=PolyMesh%poly_loc2glo(i)
      !print *,'index_loc',i,'index_glob',index_send(i)
      kiter2=kiter2+2
      !endif
   enddo

   !if(mpi_id ==1) write(*,*) mpi_id, con_quad_send
   !if(mpi_id ==1) write(*,*) '========================'
   call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

   ! ----------
   ! 5 - Broadcasting elements
   ! ----------

   proc_loop: do ip = 1, mpi_np

      if(mpi_id == ip-1) num_tria_send_mpi = 6*num_tria_send;
      !print *,num_tria_send_mpi
      call MPI_BCAST(num_tria_send_mpi, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      !call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(con_tria_send_mpi(num_tria_send_mpi))
      if(mpi_id == ip-1) con_tria_send_mpi = con_tria_send;
      !print *,num_tria_send_mpi
      call MPI_BCAST(con_tria_send_mpi,num_tria_send_mpi,&
                     MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(con_tria_recv_mpi(num_tria_send_mpi/6,6))

      kiter = 1

      do i = 1, num_tria_send_mpi, 6
         !print *,'i',i
         !print *,'kiter',kiter
         !print *,con_tria_send_mpi(i:i+5)
         con_tria_recv_mpi(kiter,1:6) = con_tria_send_mpi(i:i+5)
         !if(mpi_id == 1) write(*,*) i,i+6,con_quad_send_mpi(i:i+6)
         !if(mpi_id == 1) write(*,*) 'con_tria_rcv', con_tria_recv_mpi(kiter,1:6)
         kiter = kiter + 1
      enddo

      if(mpi_id==ip-1) num_index_send_mpi=PolyMesh%num_poly_loc

      call MPI_BCAST(num_index_send_mpi, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      allocate(index_send_mpi(num_index_send_mpi))

      if(mpi_id==ip-1) index_send_mpi=index_send

      call MPI_BCAST(index_send_mpi,num_index_send_mpi,&
                     MPI_INTEGER,ip-1,MPI_COMM_WORLD,mpi_ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      if(mpi_id==ip-1) num_hk_send_mpi=PolyMesh%num_poly_loc

      call MPI_BCAST(num_hk_send_mpi, 1, MPI_INTEGER, ip-1, MPI_COMM_WORLD, mpi_ierr)

      allocate(hk_send_mpi(num_hk_send_mpi))

      if(mpi_id==ip-1) hk_send_mpi=hk_send

      call MPI_BCAST(hk_send_mpi,2*num_hk_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      if (mpi_id==ip-1) num_poly_send_mpi = num_poly_send
      !print *,num_poly_send_mpi

      call MPI_BCAST(num_poly_send_mpi,1,MPI_INTEGER,ip-1,MPI_COMM_WORLD,mpi_ierr)
      !print *,num_poly_send_mpi

      allocate(x1_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1
         do i=1,2*num_poly_send_mpi,2
            x1_send_mpi(kiter2)=xx_send(i)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(x1_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(x2_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1
         do i=1,2*num_poly_send_mpi,2
            x2_send_mpi(kiter2)=xx_send(i+1)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(x2_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(y1_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1
         do i=1,2*num_poly_send_mpi,2
            y1_send_mpi(kiter2)=yy_send(i)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(y1_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(y2_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1
         do i=1,2*num_poly_send_mpi,2
            y2_send_mpi(kiter2)=yy_send(i+1)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(y2_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(z1_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1;
         do i=1,2*num_poly_send_mpi,2
            z1_send_mpi(kiter2)=zz_send(i)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(z1_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      ! unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      allocate(z2_send_mpi(num_poly_send_mpi))

      if (mpi_id==ip-1) then
         kiter2=1;

         do i=1,2*num_poly_send_mpi,2
            z2_send_mpi(kiter2)=zz_send(i+1)
            kiter2=kiter2+1;
         enddo
      endif

      call MPI_BCAST(z2_send_mpi,2*num_poly_send_mpi,&
                     MPI_REAL,ip-1,MPI_COMM_WORLD,mpi_ierr)

      !unpaking the elements
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

      !find neighbouring elements
      !iter_face=1;
      do i = 1, num_tria_loc
         !flag=0;
         !iface_poly=1;
         isFound_int = .false.
         if (con_tria_loc(i,1) /= 0 .and. mpi_id /= ip-1) then
            !print *,'PROC:',mpi_id
            call FIND_NEIGHBOURING_EL_MPI(i, con_tria_loc, num_tria_loc,&
                                          con_tria_recv_mpi, num_tria_send_mpi/6, &
                                          IsFound_int, Row2, IsQuad)

            if(IsFound_int) then
               mat   = con_tria_loc(i,1); mat_ne   = con_tria_recv_mpi(Row2,1)
               ie    = con_tria_loc(i,2); ie_ne    = con_tria_recv_mpi(Row2,2)
               iface = con_tria_loc(i,3); iface_ne = con_tria_recv_mpi(Row2,3);

               call GET_EL_LOC_FROM_EL_GLO(PolyMesh%elem_loc2glo, &
                                             PolyMesh%num_elem_loc, &
                                             ie,ie_loc)

               ipoly_glob=PolyMesh%elem_in_poly(ie)
               ipoly2_glob=PolyMesh%elem_in_poly(ie_ne)

               PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,0) = ip-1
               PolyMesh%Elem_loc(ie_loc)%neigh_el(iface,1:4) = [mat_ne, ie_ne, iface_ne,ipoly2_glob]


               if (ipoly_glob /= ipoly2_glob) then
                  !print *,'neighbor tet',ie,ie_ne

                  call GET_EL_LOC_FROM_EL_GLO(index_send_mpi,num_poly_send_mpi,ipoly2_glob,ipoly2_loc)

                  call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                                             PolyMesh%num_poly_loc, &
                                             ipoly_glob,ipoly_loc)
                  !print *,'tet',ie,' belongs to global poly',ipoly_glob,'local:',ipoly_loc,'in processor',PolyMesh%part_elem(ie),'iface',iface
                  !print *,'tet',ie_ne,'belongs to global poly',ipoly2_glob,'local:',ipoly2_loc,'in processor',PolyMesh%part_elem(ie_ne),'iface_ne',iface_ne

                  !print *,'bbbox:',x1_send_mpi(ipoly2_loc),x2_send_mpi(ipoly2_loc),&
                  !                  y1_send_mpi(ipoly2_loc),y2_send_mpi(ipoly2_loc),&
                  !                  z1_send_mpi(ipoly2_loc),z2_send_mpi(ipoly2_loc)
                  !print *,'hk',hk_send_mpi(ipoly2_loc)

                  num_tet_in_poly=PolyMesh%Poly(ipoly_loc)%num_tet_in_poly
                  !print *,num_tet_in_poly

                  do j=1,num_tet_in_poly
                     !print *,PolyMesh%Poly(ipoly_loc)%tet_in_poly(j)
                     if (PolyMesh%Poly(ipoly_loc)%tet_in_poly(j)==ie) then
                        iface_poly=PolyMesh%Elem_loc(ie_loc)%num_faces*(j-1)+iface
                        !print *,ie,iface,iface_poly
                     endif
                  enddo
                  !print *,iface_poly
                  PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,1,1:2)=[x1_send_mpi(ipoly2_loc),x2_send_mpi(ipoly2_loc)]
                  PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,2,1:2)=[y1_send_mpi(ipoly2_loc),y2_send_mpi(ipoly2_loc)]
                  PolyMesh%Poly(ipoly_loc)%neigh_bbox(iface_poly,3,1:2)=[z1_send_mpi(ipoly2_loc),z2_send_mpi(ipoly2_loc)]
                  PolyMesh%Poly(ipoly_loc)%neigh_hk(iface_poly)=hk_send_mpi(ipoly2_loc)
                  !index_loc(ipoly2_loc)=0;

                  !print *,'_______'

               endif

               con_tria_loc(i,:) = 0;
               con_tria_recv_mpi(Row2,:) = 0;
               num_tria_send_loc = num_tria_send_loc - 1
            endif
         endif
      enddo

      !go to next processor
      deallocate(con_tria_send_mpi,con_tria_recv_mpi)
      deallocate(x1_send_mpi,x2_send_mpi)
      deallocate(y1_send_mpi,y2_send_mpi)
      deallocate(z1_send_mpi,z2_send_mpi)
      !deallocate(xx_recv_mpi,yy_recv_mpi,zz_recv_mpi)
      deallocate(index_send_mpi,hk_send_mpi)

   enddo proc_loop

   !deallocate(faces_to_find_send,faces_to_find_send_mpi)

   call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
   call FLUSH
   if(num_tria_send_loc == 0) then
         write(*,*) 'Proc ', mpi_id, ':', ' found all tria interfaces!'
   else
         write(*,*) 'Proc ', mpi_id, ':',  num_tria_send_loc, ' tria interfaces not found!'
         !call EXIT(EXIT_NEIGHBOUR_EL_ERROR)
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)

end subroutine CREATE_NEIGH_EL_TRIA

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> creates the local-to-global map dof numbering Dof_glo stored in the struct Element;
!>
      subroutine CREATE_LOCAL2GLOBAL_MAP(PolyMesh, PolyData, mpi_np, mpi_id)

      use mpi
      use Poly_setup_MPI, only: mpi_ierr
      use Poly_mesh
      use Poly_data

      implicit none

      integer(kind=4), intent(in)  :: mpi_np, mpi_id
      integer(kind=4)  :: ie, iglo_el, mat_type, im, Deg, ndof_mat, ndof_shift, &
                          i, el_before_ie

      integer(kind=4), dimension(:), allocatable :: el_per_mat_loc, el_per_mat

      type(Mesh_Structure), intent(inout) :: PolyMesh
      type(Data_Structure), intent(inout) :: PolyData

      !> first compute how many elements per materials
      allocate(el_per_mat_loc(PolyData%nmat),el_per_mat(PolyData%nmat))
      el_per_mat_loc = 0; el_per_mat = 0

      do ie = 1, PolyMesh%num_elem_loc
         im = PolyMesh%Elem_loc(ie)%mat_prop
         el_per_mat(im) = el_per_mat(im) + 1;
      enddo


      call MPI_ALLREDUCE(el_per_mat_loc, el_per_mat, PolyData%nmat, MPI_INTEGER,  &
                         MPI_SUM, MPI_COMM_WORLD, mpi_ierr)


      do ie = 1, PolyMesh%num_elem_loc

         allocate(PolyMesh%Elem_loc(ie)%Dof_glo(PolyMesh%Elem_loc(ie)%NDof_loc))

         iglo_el  = PolyMesh%elem_loc2glo(ie);
         mat_type = PolyMesh%Elem_loc(ie)%mat_prop

         ndof_shift = 0
         do im = 1, mat_type-1

            Deg      = PolyData%sdeg_mat(im)
            ndof_mat = (deg+1)*(deg+2)*(deg+3)/6;

            ndof_shift = ndof_shift + ndof_mat*el_per_mat(im)

         enddo

         el_before_ie = iglo_el - sum(el_per_mat(1:mat_type-1)) - 1
         ndof_shift = ndof_shift + el_before_ie * PolyMesh%Elem_loc(ie)%NDof_loc

         do i = 1, PolyMesh%Elem_loc(ie)%NDof_loc

            PolyMesh%Elem_loc(ie)%Dof_glo(i) = ndof_shift + i

         enddo

      enddo

      deallocate(el_per_mat_loc,el_per_mat)

      end subroutine CREATE_LOCAL2GLOBAL_MAP

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> ???
      !>
      subroutine FIND_NEIGHBOUR_EL(iRow, array, num_row, &
                                   IsFound, jRow, IsQuad)


      implicit none

      integer(kind=4), intent(in) :: iRow, num_row
      integer(kind=4), dimension(num_row,*), intent(in) :: array
      logical, intent(inout) :: IsFound, IsQuad

      integer(kind=4), intent(out) :: jRow

      integer(kind=4) :: i


      IsFound = .false.
      jRow = 0

      if (IsQuad) then
        do i = 1, num_row
           if( array(i,1) /= 0 .and. i /= iRow) then
              if (array(iRow,4) == array(i,4) .and. &
                  array(iRow,5) == array(i,5) .and. &
                  array(iRow,6) == array(i,6) .and. &
                  array(iRow,7) == array(i,7)) then
                  IsFound = .true.
                  jRow = i
              endif
           endif
        enddo
      else
        do i = 1, num_row
           if( array(i,1) /= 0 .and. i /= iRow) then
              if (array(iRow,4) == array(i,4) .and. &
                  array(iRow,5) == array(i,5) .and. &
                  array(iRow,6) == array(i,6)) then
                  IsFound = .true.
                  jRow = i
              endif
           endif
        enddo
      endif


      end subroutine FIND_NEIGHBOUR_EL

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      !> ???
      !>
      subroutine FIND_BOUNDARY_EL(iRow, con_int, num_int, &
                                  con_bc, num_bc, IsFound, jRow, IsQuad)


      implicit none

      integer(kind=4), intent(in) :: iRow, num_int, num_bc
      integer(kind=4), dimension(num_int,*), intent(in) :: con_int
      integer(kind=4), dimension(num_bc,*),   intent(in) :: con_bc
      logical, intent(inout) :: IsFound, IsQuad

      integer(kind=4), intent(out) :: jRow

      integer(kind=4) :: i


      IsFound = .false.
      jRow = 0

      if(IsQuad) then
        do i = 1, num_bc
           !write(*,*) i, con_bc(i,1)
           if(con_bc(i,1) /= 0) then
              if (con_int(iRow,4) == con_bc(i,2) .and. &
                  con_int(iRow,5) == con_bc(i,3) .and. &
                  con_int(iRow,6) == con_bc(i,4) .and. &
                  con_int(iRow,7) == con_bc(i,5)) then
                  IsFound = .true.
                  jRow = i
              endif
           endif
        enddo
      else
        do i = 1, num_bc
           !write(*,*) i, con_bc(i,1)
           if(con_bc(i,1) /= 0) then
              if (con_int(iRow,4) == con_bc(i,2) .and. &
                  con_int(iRow,5) == con_bc(i,3) .and. &
                  con_int(iRow,6) == con_bc(i,4)) then
                  IsFound = .true.
                  jRow = i
              endif
           endif
        enddo

      endif


      end subroutine FIND_BOUNDARY_EL

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> ???
!>
subroutine FIND_NEIGHBOURING_EL_MPI(iRow, con_int, num_int,&
                                          con_recv, num_recv, &
                                          IsFound, jRow, IsQuad)

   implicit none

   integer(kind=4), intent(in) :: iRow, num_int, num_recv
   integer(kind=4), dimension(num_int,*),  intent(in) :: con_int
   integer(kind=4), dimension(num_recv,*), intent(in) :: con_recv
   logical, intent(inout) :: IsFound, IsQuad

   integer(kind=4), intent(out) :: jRow

   integer(kind=4) :: i

   IsFound = .false.
   jRow = 0

   if(IsQuad) then
      do i = 1, num_recv
         !write(*,*) i, con_bc(i,1)
         if(con_recv(i,1) /= 0) then
            if (con_int(iRow,4) == con_recv(i,4) .and. &
               con_int(iRow,5) == con_recv(i,5) .and. &
               con_int(iRow,6) == con_recv(i,6) .and. &
               con_int(iRow,7) == con_recv(i,7)) then
               IsFound = .true.
               jRow = i
            endif
         endif
      enddo
   else
      do i = 1, num_recv
         !write(*,*) i, con_bc(i,1)
         if(con_recv(i,1) /= 0) then
            if (con_int(iRow,4) == con_recv(i,4) .and. &
               con_int(iRow,5) == con_recv(i,5) .and. &
               con_int(iRow,6) == con_recv(i,6)) then
               IsFound = .true.
               jRow = i
            endif
         endif
      enddo
   endif

end subroutine FIND_NEIGHBOURING_EL_MPI

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> finds the id of the node of the mesh starting from the vertex of an element;
!>
subroutine FIND_POS_LOC_NODE(vect, dim_vect, is, it)

   implicit none

   integer(kind=4), intent(in) :: dim_vect, is
   integer(kind=4), dimension(dim_vect), intent(in) :: vect

   integer(kind=4), intent(out) :: it
   integer(kind=4) :: i

   it = 0

   do i = 1, dim_vect
      if (vect(i) == is ) then
         it = i;
         return
      endif
   enddo

   if (it==0) write(*,*) 'Error! Index not found in FIND_POS_LOC_NODE'

end subroutine FIND_POS_LOC_NODE


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> writes the mesh properties for the elements contained in each processor in the mpi files mesh.mpi;
!>
subroutine WRITE_MESH_INFO(mpi_file, PolyMesh)

   use Poly_setup_MPI
   use Poly_mesh
   use local_search

   implicit none

   type(Mesh_Structure), intent(inout) :: PolyMesh !< mesh
   character(len=70), intent(in) :: mpi_file    !< folder_mpi where to store files
   character(len=70) :: mpi_file_mesh = 'mesh_000000.mpi'
   !character(len=70) :: mat_file_mesh = 'matf0000000.m'
   !character(len=70) :: mat_file_mesh_1

   integer(kind=4)  :: ie, iface, unit_mpi, ivert,ipoly_loc,ipoly_glob, id_node, unit_mat

   unit_mpi = 40 + mpi_id
   unit_mat = 4000 + mpi_id

   if (mpi_id < 10) then
      write(mpi_file_mesh(11:11),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(mpi_file_mesh(10:11),'(i2)') mpi_id
   elseif (mpi_id < 1000) then
      write(mpi_file_mesh(9:11),'(i3)') mpi_id
   elseif (mpi_id < 10000) then
      write(mpi_file_mesh(8:11),'(i4)') mpi_id
   elseif (mpi_id < 100000) then
      write(mpi_file_mesh(7:11),'(i5)') mpi_id
   elseif (mpi_id < 1000000) then
      write(mpi_file_mesh(6:11),'(i6)') mpi_id
   endif

   mpi_file_mesh = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_mesh
   open(unit_mpi,file=mpi_file_mesh)

   do ie = 1, PolyMesh%num_elem_loc
      write(unit_mpi,*) &
            'Loc. Element #: ', ie, ' Glo. Element #: ', PolyMesh%elem_loc2glo(ie)

      write(unit_mpi,*) 'Global Dof: ', PolyMesh%Elem_loc(ie)%Dof_glo


      write(unit_mpi,*) 'Element type: ', PolyMesh%Elem_loc(ie)%el_type
      write(unit_mpi,*) 'Num Vert: ',     PolyMesh%Elem_loc(ie)%num_vert
      write(unit_mpi,*) 'Vertices #: ',   PolyMesh%Elem_loc(ie)%vert

      ipoly_glob=PolyMesh%elem_in_poly(PolyMesh%elem_loc2glo(ie))
      call GET_EL_LOC_FROM_EL_GLO(PolyMesh%poly_loc2glo, &
                                                PolyMesh%num_poly_loc, &
                                                ipoly_glob,ipoly_loc)

      write(unit_mpi,*) 'Belongs to Glo. Polyhedra', ipoly_glob, 'Loc. Polyhedra',ipoly_loc

      do ivert = 1, PolyMesh%Elem_loc(ie)%num_vert

         call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
                                 PolyMesh%Elem_loc(ie)%vert(ivert),id_node);

         write(unit_mpi,*) 'Vertices # ', PolyMesh%Elem_loc(ie)%vert(ivert), &
                           'of coords : ', PolyMesh%coord_x(id_node), &
                                          PolyMesh%coord_y(id_node), &
                                          PolyMesh%coord_z(id_node)
      enddo

      write(unit_mpi,*) 'BBox x_coord : ', PolyMesh%Poly(ipoly_loc)%b_box(1,:)
      write(unit_mpi,*) 'BBox y_coord : ', PolyMesh%Poly(ipoly_loc)%b_box(2,:)
      write(unit_mpi,*) 'BBox z_coord : ', PolyMesh%Poly(ipoly_loc)%b_box(3,:)
      write(unit_mpi,*) 'Diameter hk : ', PolyMesh%Poly(ipoly_loc)%hk
      write(unit_mpi,*) 'Num Faces: ',    PolyMesh%Elem_loc(ie)%num_faces

      ! loop on face for each element
      do iface = 1, PolyMesh%Elem_loc(ie)%num_faces
         write(unit_mpi,*) 'Face #', iface, ': ',  PolyMesh%Elem_loc(ie)%faces(iface,:)
         write(unit_mpi,*) 'Normal : ',    PolyMesh%Elem_loc(ie)%normal(iface,:)
         write(unit_mpi,*) 'Area : ',      PolyMesh%Elem_loc(ie)%area(iface)

      enddo

      write(unit_mpi,*) 'Neighbouring Elements: shared by (mpi-proc/mat_id/el_id/face_id/poly id/space_fun_tag)'
      write(unit_mpi,*) 'Neigh Face #1: ', PolyMesh%Elem_loc(ie)%neigh_el(1,:)
      write(unit_mpi,*) 'Neigh Face #2: ', PolyMesh%Elem_loc(ie)%neigh_el(2,:)
      write(unit_mpi,*) 'Neigh Face #3: ', PolyMesh%Elem_loc(ie)%neigh_el(3,:)
      write(unit_mpi,*) 'Neigh Face #4: ', PolyMesh%Elem_loc(ie)%neigh_el(4,:)

      !write(unit_mpi,*) 'Flag to see if the faces needs  to be integrated'
      !write(unit_mpi,*) 'Face #1: ', PolyMesh%Elem_loc(ie)%flag(1)
      !write(unit_mpi,*) 'Face #2: ', PolyMesh%Elem_loc(ie)%flag(2)
      !write(unit_mpi,*) 'Face #3: ', PolyMesh%Elem_loc(ie)%flag(3)
      !write(unit_mpi,*) 'Face #4: ', PolyMesh%Elem_loc(ie)%flag(4)


      if(PolyMesh%Elem_loc(ie)%el_type == 'PRY' .or. &
            PolyMesh%Elem_loc(ie)%el_type == 'HEX') &
      write(unit_mpi,*) 'Neigh Face #5: ',PolyMesh%Elem_loc(ie)%neigh_el(5,:)

      if(PolyMesh%Elem_loc(ie)%el_type == 'HEX') &
      write(unit_mpi,*) 'Neigh Face #6: ',PolyMesh%Elem_loc(ie)%neigh_el(6,:)
      write(unit_mpi,*) '--------------------------------'

   enddo

   close(unit_mpi)

   ! matfile output for debugging
   !-----------------------------------------------------------------------------

   !do ie = 1, PolyMesh%num_elem_loc

   !  if (ie < 10) then
   !     write(mat_file_mesh(11:11),'(i1)') ie
   !  elseif (ie < 100) then
   !     write(mat_file_mesh(10:11),'(i2)') ie
   !  elseif (ie < 1000) then
   !     write(mat_file_mesh(9:11),'(i3)') ie
   !  elseif (ie < 10000) then
   !     write(mat_file_mesh(8:11),'(i4)') ie
   !  elseif (ie < 100000) then
   !     write(mat_file_mesh(7:11),'(i5)') ie
   !  elseif (ie < 1000000) then
   !     write(mat_file_mesh(6:11),'(i6)') ie
   !  endif

   !  mat_file_mesh_1 = mpi_file(1:len_trim(mpi_file)) // '/' // mat_file_mesh
   !  open(unit_mat,file=mat_file_mesh_1)

   !  write(unit_mat,*) 'vert = [ '

   !  do ivert = 1, PolyMesh%Elem_loc(ie)%num_vert

   !    call FIND_POS_LOC_NODE(PolyMesh%node_loc2glo,PolyMesh%num_node_loc, &
   !                            PolyMesh%Elem_loc(ie)%vert(ivert),id_node);

   !    write(unit_mat,*) PolyMesh%coord_x(id_node), &
   !                     PolyMesh%coord_y(id_node), &
   !                     PolyMesh%coord_z(id_node)
   !   enddo

   !   write(unit_mat,*) ' ];'
   !   write(unit_mat,*)

   !   write(unit_mat,*) 'faces = ['
   !   do i = 1, PolyMesh%Elem_loc(ie)%num_faces
   !     write(unit_mat,*) PolyMesh%Elem_loc(ie)%faces(i,:)
   !   enddo
   !   write(unit_mat,*)  '];'

   !   write(unit_mat,*) 'normal = ['
   !   do i = 1, PolyMesh%Elem_loc(ie)%num_faces
   !     write(unit_mat,*) PolyMesh%Elem_loc(ie)%normal(i,:)
   !   enddo
   !   write(unit_mat,*)  '];'

   !   write(unit_mat,*) 'bbox = ['
   !   do i = 1,3
   !     write(unit_mat,*) PolyMesh%Elem_loc(ie)%b_box(i,:)
   !   enddo
   !   write(unit_mat,*)  '];'

   !   close(unit_mat)

   ! enddo


end subroutine WRITE_MESH_INFO

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine WRITE_INTERFACE_INFO(mpi_file, PolyMesh)

   use Poly_setup_MPI
   use Poly_mesh

   implicit none

   type(Mesh_Structure), intent(inout) :: PolyMesh !< mesh
   character(len=70), intent(in) :: mpi_file    !< folder_mpi where to store files
   character(len=70) :: mpi_file_interface = 'interface_000000.mpi'

   integer(kind=4)  :: i, ie, iface, inter_count, unit_int, elem_proc_id, j, count

   integer(kind=4) :: num_inter_loc
   integer(kind=4), dimension(mpi_np) :: num_elem_interface_loc
   integer(kind=4), dimension(mpi_np, PolyMesh%Elem_loc(1)%num_faces) :: elem_interface_tmp !< temp matrix for worst case scenario, interface on all faces
   integer(kind=4), dimension(mpi_np) :: recvcounts, displs
   integer(kind=4), dimension(mpi_np) :: inter_displs_loc    !< vector containing displacement index for each process (progressive)

   integer(kind=4), dimension(mpi_np,mpi_np) :: tmp_mat
   integer(kind=4), dimension(:), allocatable :: elem_inter_loc

   unit_int = 40000 + mpi_id

   if (mpi_id < 10) then
      write(mpi_file_interface(16:16),'(i1)') mpi_id
   elseif (mpi_id < 100) then
      write(mpi_file_interface(15:16),'(i1)') mpi_id
   elseif (mpi_id < 1000) then
      write(mpi_file_interface(14:16),'(i1)') mpi_id
   elseif (mpi_id < 10000) then
      write(mpi_file_interface(13:16),'(i1)') mpi_id
   elseif (mpi_id < 100000) then
      write(mpi_file_interface(12:16),'(i1)') mpi_id
   elseif (mpi_id < 1000000) then
      write(mpi_file_interface(11:16),'(i1)') mpi_id
   endif

   mpi_file_interface = mpi_file(1:len_trim(mpi_file)) // '/' // mpi_file_interface
   open(unit_int,file=mpi_file_interface)

   ! loop over faces of each element
   ! count local number of interfaces
   num_inter_loc = 0
   num_elem_interface_loc = 0
   elem_interface_tmp = -1
   do ie = 1, PolyMesh%num_elem_loc
      do iface = 1, PolyMesh%Elem_loc(ie)%num_faces
         elem_proc_id = PolyMesh%Elem_loc(ie)%neigh_el(iface,0)
         ! check process and store IDs if in a different one and count
         if (mpi_id /= elem_proc_id) then
            ! count total number of interfaces
            num_inter_loc = num_inter_loc + 1
            ! count number of interfaces per process
            num_elem_interface_loc(elem_proc_id+1) = num_elem_interface_loc(elem_proc_id+1) + 1
            ! insert global element ID into temporary matrix following progressive indexing
            inter_count = num_elem_interface_loc(elem_proc_id+1)
            elem_interface_tmp(elem_proc_id+1, inter_count) = PolyMesh%Elem_loc(ie)%neigh_el(iface,2)
         endif
      enddo
   enddo
   write(unit_int,*) num_inter_loc

   ! allocate after knowing the number of interfaces per processor
   allocate(elem_inter_loc(num_inter_loc))

   ! find total number of interfaces and allocate memory for global data
   call MPI_ALLREDUCE(num_inter_loc, PolyMesh%num_elem_inter, 1, MPI_INTEGER,  &
                         MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
   if (mpi_id==0) print *, 'total number of interfaces for comm:', PolyMesh%num_elem_inter

   ! allocate after knowing the total number of interfaces
   allocate(PolyMesh%elem_inter(PolyMesh%num_elem_inter))

   ! loop over processes - progressive counter increasing with each new number of elements to send per processor
   inter_count = 1
   do i=1,mpi_np
      num_inter_loc = num_elem_interface_loc(i)
      if (num_inter_loc /= 0) then
         elem_inter_loc(inter_count:inter_count+num_inter_loc-1) = elem_interface_tmp(i, 1:num_inter_loc)
         inter_count = inter_count + num_inter_loc
      endif
   end do

   ! write to file
   do ie = 1, PolyMesh%num_elem_loc
      ! loop on face for each element
      do iface = 1, PolyMesh%Elem_loc(ie)%num_faces
         ! write to file interface element IDs and process
         elem_proc_id = PolyMesh%Elem_loc(ie)%neigh_el(iface,0)
         if (mpi_id /= elem_proc_id) then
            write(unit_int,*) elem_proc_id, PolyMesh%Elem_loc(ie)%neigh_el(iface,2)
         endif
      enddo
   enddo

   close(unit_int)

   ! save the number of interface elements to send and receive for all processes
   allocate(PolyMesh%num_elem_inter_comm(mpi_np,mpi_np))
   PolyMesh%num_elem_inter_comm = 0
   recvcounts = mpi_np  ! Each process sends mpi_np elements
   displs = (/ (i*mpi_np, i=0, mpi_np-1) /)  ! Offsets for each process

   ! gather number of elements to send/receive for each processor
   ! incidence matrix
   call MPI_ALLGATHERV(num_elem_interface_loc, mpi_np, MPI_INTEGER, PolyMesh%num_elem_inter_comm, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)

   !! fortran saves column major, data sent by row
   !PolyMesh%num_elem_inter_comm = transpose(PolyMesh%num_elem_inter_comm)

   allocate(PolyMesh%inter_disp(mpi_np,mpi_np))
   ! call MPI_ALLGATHERV(inter_displs_loc, mpi_np, MPI_INTEGER, PolyMesh%inter_disp(1,1), recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)
   ! print *, PolyMesh%inter_disp

   ! progressive counter
   count = 0
   do i=1,mpi_np
      do j=1,mpi_np
         PolyMesh%inter_disp(i,j) = count
         count = count + PolyMesh%num_elem_inter_comm(i,j)
      enddo
   enddo

   ! size of each list containing the number of elements per process
   do i=1,mpi_np-1
      recvcounts(i) = PolyMesh%inter_disp(i+1,1) - PolyMesh%inter_disp(i,1)
   enddo
   recvcounts(mpi_np) = PolyMesh%num_elem_inter - PolyMesh%inter_disp(mpi_np,1)
   displs = PolyMesh%inter_disp(:,1)

   ! gather ordered (by process) lists (each of different legth) into a single array
   call MPI_ALLGATHERV(elem_inter_loc, recvcounts(mpi_id+1), MPI_INTEGER, PolyMesh%elem_inter, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, mpi_ierr)


   if (mpi_id == 0) then
      print *, "num_elem_inter_comm row-by-row:"
      do i = 1, mpi_np  ! Loop over rows
         print *, PolyMesh%num_elem_inter_comm(i, :)
      end do

      print *, "inter_disp row-by-row:"
      do i = 1, mpi_np  ! Loop over rows
         print *, PolyMesh%inter_disp(i, :)
      end do

      print *, 'elem_inter'
      print *, PolyMesh%elem_inter
   endif

end subroutine WRITE_INTERFACE_INFO

end module mesh_partition_and_mpi_files