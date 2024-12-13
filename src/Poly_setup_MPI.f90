!    Copyright (C) 2021 The SPEED FOUNDATION
!    Author: Ilario Mazzieri

!> @brief Module for initialization of the MPI and PETSc envinronments with variables declaration.
module Poly_setup_MPI
#include<petsc/finclude/petscksp.h>

     use mpi
     use petscksp

     implicit none

     integer(kind=4), dimension (:), allocatable :: mpi_stat
     integer(kind=4) :: mpi_id          !< MPI process ID
     integer(kind=4) :: mpi_np          !< number of MPI processes
     integer(kind=4) :: mpi_ierr        !< MPI error code
     logical(kind=4) :: flag            !< flag is true if MPI_INIT has been called

     type(PetscErrorCode) :: ierr             !< PETSc error code
     type(PetscMPIInt) :: mpi_size            !< PETSc MPI size
     type(PetscMPIInt) :: rank                !< PETSc MPI process ID

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

!> @brief Subroutine for initialization of the MPI and PETSc envinronments.
subroutine INITIALIZATION()

     ! MPI INITIALIZATION

     call MPI_INIT(mpi_ierr)
     call MPI_Initialized(flag,mpi_ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id, mpi_ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_np, mpi_ierr)

     PetscCall(PetscInitialize(ierr))
     PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,mpi_size,ierr))
     PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))

     ! Initialize PETSc
     ! call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
     !call MPI_Comm_size(PETSC_COMM_WORLD, mpi_size, ierr)
     !call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

     ! Console output in case of error

     if (mpi_ierr.ne.0 .or. ierr .ne. 0) then
     ! Error occured on rank number: ....
          write(*,*)'MPI Initialization error - proc : ',mpi_id
     endif

     !What does it mean?
     !speed_tag = speed_tag_min
     !return

end subroutine INITIALIZATION

end module Poly_setup_MPI



