!    Copyright (C) 2021 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of PolyWAVE.
!
!    PolyWAVE is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    PolyWAVE is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with PolyWAVE.  If not, see <http://www.gnu.org/licenses/>.
 
     module Poly_setup_MPI
#include<petsc/finclude/petscksp.h>

     use mpi
     use petscksp
     

     implicit none

     integer(kind=4), dimension (:), allocatable :: mpi_stat
     integer(kind=4) :: mpi_id, mpi_np, mpi_ierr
     logical(kind=4) :: flag
     
     PetscErrorCode :: ierr
     PetscMPIInt :: size,rank

     !integer*4 POLYSPEED_COMM
     !integer*4 POLYSPEED_TAG, POLYSPEED_TAG_MIN, POLYSPEED_TAG_MAX
     !integer*4 POLYSPEED_INTEGER, POLYSPEED_REAL, POLYSPEED_DOUBLE
     !parameter (POLYSPEED_COMM = MPI_COMM_WORLD)
     !integer(kind=4) :: mpi_comm=POLYSPEED_COMM
     !parameter (POLYSPEED_STATUS_SIZE = MPI_STATUS_SIZE)
     !parameter (POLYSPEED_TAG_MIN = 1515, POLYSPEED_TAG_MAX = 1530)
     !parameter (POLYSPEED_INTEGER = MPI_INTEGER)
     !parameter (POLYSPEED_REAL = MPI_REAL)
     !parameter (POLYSPEED_DOUBLE = MPI_DOUBLE_PRECISION)


      contains      
      
          subroutine INITIALIZATION()

          !> MPI INITIALIZATION
          
          call MPI_INIT(mpi_ierr)
          call MPI_Initialized(flag,mpi_ierr)
          call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id, mpi_ierr)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_np, mpi_ierr)

          PetscCall(PetscInitialize(ierr))
          PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,size,ierr))
          PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))

          !> Console output in case of error
          
          if (mpi_ierr.ne.0 .or. ierr .ne. 0) then
          !> Error occured on rank number: ....
               write(*,*)'MPI Initialization error - proc : ',mpi_id
          endif

          !What does it mean?
          !speed_tag = speed_tag_min
          !return
          
          end subroutine INITIALIZATION
   
     end module Poly_setup_MPI
     
     
     
