!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario sazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.


!> @brief Reads input files and allocates memory.
!! @author Ilario Mazzieri
!> @date November, 2014
!> @version 1.0


subroutine READ_INPUT_FILES(PolyData,PolyMesh)

    !|use mpi
    use Poly_setup_mpi
    use Poly_global
    use Poly_exit_codes
    use Poly_fail_codes
    use Poly_data
    use Poly_mesh

    implicit none


    type(Data_Structure), intent(inout) :: PolyData
    type(Mesh_Structure), intent(inout) :: PolyMesh


! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Reading header file

    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A)')'------------------Reading Header File------------------'
        write(*,'(A)')
        write(*,'(A,A20)') 'Header File: ', head_file
    endif

    inquire(file=head_file,exist=IS_filefound);
    if(IS_filefound .eqv. .false.) then
        write(*,*) 'File ', head_file, ' is missing!'
        call EXIT(EXIT_MISSING_FILE)
    endif


    call READ_HEADER(head_file)

    if (mpi_id == 0) write(*,'(A)')'Read.'

    ! Instability control
    if (IS_instabilitycontrol) then
        if (mpi_id == 0) then
            write(*,'(A,L)') 'Instability control  :   ', IS_instabilitycontrol
            write(*,'(A,E8.3)') 'Instability threshold: ', 100.d0
        endif
    endif


    ! Restart control
    if (IS_restart) then
        if (mpi_id == 0) then
            write(*,'(A,L)') 'Restart Active  :   ', IS_restart
            write(*,'(A,A)') 'Backup file are saved in  : ', folder_restart
            num_dt_restart = int(time_restart/time_step)
            write(*,'(A,I10,A)')  'Backup every ',  num_dt_restart, ' time steps'
            if (time_restart > stop_time) &
                write(*,'(A)') 'No backup files since Restart time > Stop time'
        endif
    endif

!
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Reading material file


    mate_file = mate_file(1:len_trim(mate_file)) // '.mate'

    if (mpi_id == 0) then
        write(*,'(A)')
        write(*,'(A,A20)')'-----------------Reading Material File-----------------'
        write(*,'(A,A20)')'Material File    : ',mate_file
    endif

    inquire(file=mate_file,exist=IS_filefound);
    if(IS_filefound .eqv. .false.) call EXIT(EXIT_MISSING_FILE)


    call set_Data_Structure_default_value(PolyData)

    call READ_DIME_MATEFILE(mate_file,PolyData)

    if (mpi_id == 0) call print_Dime_Data_Structure(PolyData)

    call allocate_Data_Structure(PolyData)

    call READ_MATEFILE(mate_file,PolyData)

    if (mpi_id == 0) call print_Data_Structure(PolyData)

!
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Reading grid file


    grid_file = grid_file(1:len_trim(grid_file)) // '.mesh'
    if (mpi_id.eq.0)  write(*,'(A)') '-------------------Reading Grid File-------------------'
    if (mpi_id.eq.0)  write(*,'(A,A35)') 'Grid File : ',grid_file

    inquire(file=grid_file,exist=IS_filefound);
    if(IS_filefound .eqv. .false.) call EXIT(EXIT_MISSING_FILE)


    !counting hexahedras and squares
    call READ_DIME_MESHFILE(grid_file,PolyData,PolyMesh)

    if (mpi_id == 0) call print_Dime_Mesh_Structure(PolyMesh)

    call allocate_Mesh_Structure(PolyMesh)

    call READ_MESHFILE(grid_file,PolyData,PolyMesh)

    !call allocate_Poly_in_Mesh_Structure(PolyMesh)

!    if (mpi_id == 0) call print_Mesh_Structure(PolyData)


    if(mpi_id .eq. 0) write(*,'(A)') 'Read.'
    if(mpi_id .eq. 0) write(*,'(A)')

end subroutine READ_INPUT_FILES
