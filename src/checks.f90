module checks

    use mpi
    use Poly_setup_mpi
    use Poly_fail_codes
    use Poly_exit_codes
    use utilities

    implicit none

    contains

    !> check if CFL is ok
    subroutine COMPUTE_CFL(h_max, dt)

        !TODO find correct speed and constant for CFL
        implicit none

        real(kind=8), intent(in) :: h_max !< maximum grid diameter
        real(kind=8), intent(in) :: dt !< time-step

        real(kind=8) :: speed !< speed

        speed = 1.0d0

        if (mpi_id==0) then

            if (dt > h_max/speed) then
                print *, "!-------- FAIL CFL CONDITION --------!"
                if (IS_failCFL .eqv. .true.) then
                    call LYMPH3D_BARRIER
                    call PetscFinalize(mpi_ierr)
                    call MPI_FINALIZE(mpi_ierr)
                    call EXIT(EXIT_CFL)
                endif
            endif

        endif

    end subroutine COMPUTE_CFL

end module checks