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

!> @brief Poly_global - Setup variables.
module Poly_global

    implicit none

    ! header file names
    !> Name of the header file for simulation options
    character(len = 14) :: head_file = 'Poly.input'
    !> Name of the mesh file
    character(len = 75) :: grid_file
    !> Name of the file with materials, polynomial degrees, boundary conditions
    character(len = 75) :: mate_file
    !> Name of the folder where to save MPI files
    character(len = 70) :: folder_mpi
    !> Name of the folder where to save solutions !!
    character(len = 70) :: folder_monitors
    !> Name of the folder where to read the restart !!
    character(len = 70) :: folder_restart

    !time discretization paramters
    !> timestep for time integration
    real(kind=8)        :: time_step
    !> initial time for time integration
    real(kind=8)        :: start_time
    !> stop time for time integration
    real(kind=8)        :: stop_time
    !> time for restart the simulation !!
    real(kind=8)        :: time_restart
    !> number of timestep to skip to save solution !!
    integer(kind=4)     :: num_dt_mon
    !> number of timestep for restart !!
    integer(kind=4)     :: num_dt_restart

    real(kind=8) :: half_dt2, dt2

    !dampung type: 1) Q frequency proportional/ 2) Q frequency constant
    !> damping type
    integer(kind=4) :: damping_type

    !option for output
    !> output variable !!
    integer(kind=4), dimension(6) :: opt_out_var

    !Restart and Debug
    !> logical for restart !!
    logical :: Is_Restart
    !> logical for debug !!
    logical :: Is_Debug

    !option for output file list
    !> ?!!
    real(kind=8)    :: depth_search_mon_lst
    !> ?!!
    logical         :: IS_mon_lst

    !measuring computational time
    !> wall time hours
    integer(kind=4) :: time_hour
    !> wall time start
    integer(kind=4) :: time_min
    !> wall time seconds
    integer(kind=4) :: time_sec
    !> time of the beginning of simulation
    real(kind=8)    :: start
    !> time of end of simulation
    real(kind=8)    :: finish

    ! Time profiling
    real(kind=8)    :: tp_copy_matrix = 0.d0
    real(kind=8)    :: tp_linear_system = 0.d0
    real(kind=8)    :: tp_KU = 0.d0
    real(kind=8)    :: tp_copy_vector = 0.d0
    real(kind=8)    :: tp_system_setup = 0.d0
    real(kind=8)    :: tp_setup_K = 0.d0
    real(kind=8)    :: tp_setup_M = 0.d0
    real(kind=8)    :: tp_setup_RHS = 0.d0
    real(kind=8)    :: tp_partition = 0.d0
    real(kind=8)    :: tp_export = 0.d0
    real(kind=8)    :: tp_exact = 0.d0
    real(kind=8)    :: tp_error = 0.d0


    !file found
    !> logical variable to identify if file was found
    logical :: IS_filefound

    !time dependence
    !> logical variable for time dependent problems
    logical :: IsTime_dependent

    !> logical variable for type of solver
    logical :: IS_MatrixFree

    ! save outputs
    !>  logical variable to save outputs
    logical :: IsSave_output


end module Poly_global


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    EXIT CODES
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> @brief Exit codes !! CHECK
module Poly_exit_codes

    implicit none

    !> Normal exit code
    integer, parameter :: EXIT_NORMAL         = 0
    !> Exit due to CFL (Courant–Friedrichs–Lewy) condition violation
    integer, parameter :: EXIT_CFL            = 1
    !> Exit due to numerical instability
    integer, parameter :: EXIT_INSTAB         = 2
    !> Exit due to anelastic condition
    integer, parameter :: EXIT_ANELASTIC      = 3
    !> Exit due to setup error
    integer, parameter :: EXIT_SETUP          = 4
    !> Exit due to singular matrix error
    integer, parameter :: EXIT_SINGULARMTX    = 5
    !> Exit because surface was not found
    integer, parameter :: EXIT_SURF_NOTFOUND  = 6
    !> Exit due to energy error
    integer, parameter :: EXIT_ENERGY_ERROR   = 7
    !> Exit due to syntax error
    integer, parameter :: EXIT_SYNTAX_ERROR   = 8
    !> Exit because a required file was missing
    integer, parameter :: EXIT_MISSING_FILE   = 9
    !> Exit due to root-finding error
    integer, parameter :: EXIT_ROOT           = 10
    !> Exit due to element orientation error
    integer, parameter :: EXIT_ELEM_ORIENT    = 11
    !> Exit because nodes were not found
    integer, parameter :: EXIT_NO_NODES       = 12
    !> Exit because elements were not found
    integer, parameter :: EXIT_NO_ELEMENTS    = 13
    !> Exit due to damping peak error
    integer, parameter :: EXIT_DAMPING_PEAK   = 14
    !> Exit because materials were not found
    integer, parameter :: EXIT_NO_MATERIALS   = 15
    !> Exit due to function error
    integer, parameter :: EXIT_FUNCTION_ERROR = 16
    !> Exit due to non-honoring error
    integer, parameter :: EXIT_NOTHONORING_ERROR = 17
    !> Exit due to neighboring element error
    integer, parameter :: EXIT_NEIGHBOUR_EL_ERROR = 18

end module Poly_exit_codes


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    FAIL CODES
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> Module containing specific fail codes for specific conditions and behaviors
module Poly_fail_codes

    implicit none

    !> if .TRUE., fails on negative anelastic coefficients (see: damping 2)
    logical :: IS_failoncoeffs

    !> if .TRUE., do not start the TIME_LOOP @n
    !> only setup mesh, parameters, CFL, etc. then quit
    logical :: IS_setuponly

    !> if .TRUE:, quit if CFL does not hold
    logical :: IS_failCFL

    !> if .TRUE., quit if simulation becomes unstable
    logical :: IS_instabilitycontrol

end module Poly_fail_codes

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    DEFAULT VALUES
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> @brief Module containing default error values
module Poly_default_codes

    implicit none

    ! Default values
    integer(kind=4), parameter :: damping_type_default = 1

    real(kind=8), parameter :: start_time_default = 0.d0

    !!! INCONSISTENT CAPITALIZATION
    logical, parameter :: IS_mon_lst_default = .false.
    logical, parameter :: IS_Restart_default = .false.
    logical, parameter :: IS_Debug_default = .false.
    logical, parameter :: IS_failoncoeffs_default = .false.
    logical, parameter :: IS_setuponly_default = .false.
    logical, parameter :: IS_failCFL_default = .false.
    logical, parameter :: IS_instabilitycontrol_default = .false.

    logical, parameter :: IS_timedependent_default = .false.
    logical, parameter :: IS_saveoutput_default = .false.
    logical, parameter :: IS_MatrixFree_default = .false.

end module Poly_default_codes


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    QSORT MODULE
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> @brief Module containing quick sort algorithm for vectors
module qsort

    implicit none
    public :: QsortC
    private :: Partition

    contains
    !> @brief Quick sort algorithm subroutine for vectors
    recursive subroutine QsortC(vec)

        integer(kind=4), intent(inout), dimension(:) :: vec !< vector to sort
        integer(kind=4) :: iq

        if(size(vec) > 1) then
            call Partition(vec, iq)
            call QsortC(vec(:iq-1))
            call QsortC(vec(iq:))
        endif

    end subroutine QsortC

    !> @brief Subroutine in quick sort to partition a vector
    subroutine Partition(vec, marker)

        integer(kind=4), intent(inout), dimension(:) :: vec    !< vector to partition
        integer(kind=4), intent(out) :: marker  !< pivot index
        integer(kind=4) :: i, j, temp, x

        x = vec(1)
        i = 0
        j = size(vec) + 1

        do
            j = j-1
            do
                if (vec(j) <= x) exit
                j = j-1
            end do
            i = i+1
            do
                if (vec(i) >= x) exit
                i = i+1
            end do

            if (i < j) then
                ! exchange vec(i) and vec(j)
                temp = vec(i)
                vec(i) = vec(j)
                vec(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do

    end subroutine Partition

end module qsort

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!> @brief Module to get a local element from a global element
module local_search

    implicit none

    public :: GET_EL_LOC_FROM_EL_GLO

    contains

    !! NECESSARY?
    !> @brief Subroutine to get a local element from a global element
    subroutine GET_EL_LOC_FROM_EL_GLO(v, n_el , ie, ie_loc)

        integer(kind=4), intent(inout) :: n_el      !< number of elements
        integer(kind=4), intent(inout), dimension(n_el) :: v    !< input vector
        integer(kind=4), intent(in) :: ie       !< ID of global element
        integer(kind=4), intent(out) :: ie_loc  !< ID of local element - 0 if not found on the process
        integer(kind=4) :: i

        ie_loc = 0
        do i = 1, n_el
            if (v(i) == ie) then
                ie_loc = i
                return
            endif
        enddo

    end subroutine GET_EL_LOC_FROM_EL_GLO

end module local_search

!> @brief Module to find a tetrahedron in polyhedron
module find_poly
    implicit none

    public :: FIND_TET_IN_POLY

    contains

    !> @brief Function to find a tetrahedron in polyhedron
    function FIND_TET_IN_POLY(array,val,N)result(index)
        integer(kind=4) :: N    !<
        integer(kind=4) :: val  !<
        integer(kind=4),dimension(N) :: array   !< input array
        integer(kind=4),dimension(N) :: temp1,temp2
        integer(kind=4),dimension(:), allocatable :: index !<
        integer(kind=4) i,ii
        ii=1

        do i=1,N

            if (array(i)  == val) then
                temp1(i)=1
            else
                temp1(i)=0
            end if
        end do
        do i=1,N
            if (temp1(i) == 1) then
                temp2(ii)=i
                ii=ii+1
            end if
        end do

        if (ii >1 ) then
            ALLOCATE(index(ii-1))
            do i=1,ii-1
                index(i)=temp2(i)
            end do
        else
            allocate(index(1))
            index(1)=0
        end if

    end function FIND_TET_IN_POLY

end module find_poly

! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    calc_time: description
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!> @brief Subroutine to pass from seconds to hh:mm:ss format
subroutine calc_time(time_h, time_m, time_s, time_in_seconds)

    implicit none

    integer(kind=4), intent(out) :: time_h  !< hours
    integer(kind=4), intent(out) :: time_m  !< minutes
    integer(kind=4), intent(out) :: time_s  !< seconds
    integer(kind=4), intent(in)  :: time_in_seconds !< total duration in seconds
    real(kind=8)                 :: rem_min

    time_h  = int(floor(real(time_in_seconds)/3600))
    rem_min = mod(time_in_seconds,3600)
    time_m  = int(floor(rem_min/60))
    time_s  = mod(mod(time_in_seconds,3600),60)

end subroutine calc_time

!> Global constants
module global_parameters

    implicit none

    type :: ScatteredArray
        real(kind=8), dimension(:), allocatable :: data
    end type ScatteredArray

    integer(kind=4), parameter :: DIM = 3               !< dimension of the problem
    integer(kind=4), parameter :: NVERT_TRIA = 3        !< number of vertices of a triangle
    integer(kind=4), parameter :: NVERT_QUAD = 4        !< number of vertices of a quadrilateral
    integer(kind=4), parameter :: NVERT_TET = 4         !< number of vertices of a tetrahedron
    integer(kind=4), parameter :: NVERT_HEX = 6         !< number of vertices of a hexahedron
    real(kind=8), parameter :: SQRT2 = sqrt(2.0d0)      !< parameter for @f$ \sqrt{2} @f$
    real(kind=8), parameter :: PI = 4.d0*datan(1.0d0)   !< parameter for @f$ \pi @f$
    real(kind=8), parameter :: TOL = 1.0d-40            !< tolerance for small numbers

end module global_parameters