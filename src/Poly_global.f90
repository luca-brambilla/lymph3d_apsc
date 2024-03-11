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
 
     module Poly_global

     implicit none 
     
     !header file names
     character(len = 14) :: head_file = 'Poly.input'
     character(len = 75) :: grid_file, mate_file 
     character(len = 70) :: folder_mpi, folder_monitors, folder_restart

     !time discretization paramters
     real(kind=8)        ::  time_step, start_time, stop_time, time_restart
     integer(kind=4)     ::  num_dt_mon, num_dt_restart
     
     !dampung type: 1) Q frequency proportional/ 2) Q frequency constant 
     integer(kind=4) :: damping_type
     
     !option for output
     integer(kind=4), dimension(6) :: opt_out_var
                   
     !Restart and Debug
     logical :: Is_Restart, Is_Debug 

     !option for output file list 
     real(kind=8)    :: depth_search_mon_lst 
     logical         :: IS_mon_lst          
                                  
     !measuring computational time
     integer(kind=4) :: time_hour, time_min, time_sec
     real(kind=8)    :: start, finish
     
     !file found 
     logical :: IS_filefound  
     
     
     !parameters
     real(kind=8), parameter :: pi = 4.d0*datan(1.d0)  !< p-greek
     
     end module Poly_global
     

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    EXIT CODES 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     module Poly_exit_codes

     implicit none

     integer, parameter :: EXIT_NORMAL         = 0
     integer, parameter :: EXIT_CFL            = 1
     integer, parameter :: EXIT_INSTAB         = 2
     integer, parameter :: EXIT_ANELASTIC      = 3
     integer, parameter :: EXIT_SETUP          = 4
     integer, parameter :: EXIT_SINGULARMTX    = 5
     integer, parameter :: EXIT_SURF_NOTFOUND  = 6
     integer, parameter :: EXIT_ENERGY_ERROR   = 7
     integer, parameter :: EXIT_SYNTAX_ERROR   = 8
     integer, parameter :: EXIT_MISSING_FILE   = 9
     integer, parameter :: EXIT_ROOT           = 10
     integer, parameter :: EXIT_ELEM_ORIENT    = 11
     integer, parameter :: EXIT_NO_NODES       = 12
     integer, parameter :: EXIT_NO_ELEMENTS    = 13
     integer, parameter :: EXIT_DAMPING_PEAK   = 14
     integer, parameter :: EXIT_NO_MATERIALS   = 15
     integer, parameter :: EXIT_FUNCTION_ERROR = 16
     integer, parameter :: EXIT_NOTHONORING_ERROR = 17
     integer, parameter :: EXIT_NEIGHBOUR_EL_ERROR = 18
     
     end module Poly_exit_codes


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    FAIL CODES 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
     module Poly_fail_codes
      
      implicit none 
      
      ! if .TRUE., fails on negative anelastic coefficients (see: damping 2)
      logical :: IS_failoncoeffs

      ! if .TRUE., do not start the TIME_LOOP
      ! only setup mesh, parameters, CFL, etc. then quit
      logical :: IS_setuponly

      ! if .TRUE:, quit if CFL does not hold
      logical :: IS_failCFL

      ! if .TRUE., quit if simulation becomes unstable
      logical :: IS_instabilitycontrol

     end module Poly_fail_codes

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    DEFAULT VALUES 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      module Poly_default_codes
      
      implicit none

      ! Default values
      integer(kind=4), parameter :: damping_type_default = 1
      
      real(kind=8), parameter :: start_time_default = 0.d0;
      
      logical, parameter :: IS_mon_lst_default = .FALSE.
      logical, parameter :: IS_Restart_default = .FALSE.
      logical, parameter :: IS_Debug_default = .FALSE.
      logical, parameter :: IS_failoncoeffs_default = .FALSE.
      logical, parameter :: IS_setuponly_default = .FALSE.
      logical, parameter :: IS_failCFL_default = .FALSE.
      logical, parameter :: IS_instabilitycontrol_default = .FALSE.    


      end module Poly_default_codes
     

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    QSORT MODULE 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      module qsort     

      implicit none
      public :: QsortC
      private :: Partition

      contains

        recursive subroutine QsortC(A)
        
        integer(kind=4), intent(in out), dimension(:) :: A
        integer(kind=4) :: iq

        if(size(A) > 1) then
          call Partition(A, iq)
          call QsortC(A(:iq-1))
          call QsortC(A(iq:))
        endif
        
        end subroutine QsortC

        subroutine Partition(A, marker)
  
        integer(kind=4), intent(in out), dimension(:) :: A
        integer(kind=4), intent(out) :: marker
        integer(kind=4) :: i, j, temp, x

        x = A(1)
        i = 0
        j = size(A) + 1

        do
          j = j-1
          do
            if (A(j) <= x) exit
            j = j-1
          end do
          i = i+1
          do
            if (A(i) >= x) exit
            i = i+1
          end do
         
          if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      module local_search
      
      implicit none
      
      public :: GET_EL_LOC_FROM_EL_GLO
      
      contains 
          
          subroutine GET_EL_LOC_FROM_EL_GLO(v, n_el , ie, ie_loc)
         
          integer(kind=4), intent(in out) :: n_el
          integer(kind=4), intent(in out), dimension(n_el) :: v
          integer(kind=4), intent(in) :: ie
          integer(kind=4), intent(out) :: ie_loc
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

module find_poly
  implicit none
  
  public :: FIND_TET_IN_POLY
  
  contains

  function FIND_TET_IN_POLY(array,val,N)result(index)
    integer(kind=4) :: N
    integer(kind=4) :: val
    integer(kind=4),dimension(N) :: array
    integer(kind=4),dimension(N) :: temp1,temp2
    integer(kind=4),dimension(:), ALLOCATABLE :: index
    integer(kind=4) i,ii
    ii=1

    do i=1,N

        if (array(i)  .eq. val) then
            temp1(i)=1
        else
            temp1(i)=0
        end if
    end do
    do i=1,N
        if (temp1(i) .eq. 1) then
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

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    calc_time: description 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
     subroutine calc_time(time_h, time_m, time_s, time_in_seconds)
               
        implicit none
        
        integer(kind=4), intent(out) :: time_h, time_m, time_s
        integer(kind=4), intent(in)  :: time_in_seconds
        real(kind=8)                 :: rem_min
        
        time_h  = int(floor(real(time_in_seconds)/3600));
        rem_min = mod(time_in_seconds,3600)
        time_m  = int(floor(rem_min/60))
        time_s  = mod(mod(time_in_seconds,3600),60)
         
        
     end subroutine calc_time
 
 
 
 
            

