module checks

    use Poly_setup_mpi
    use Poly_fail_codes
    use Poly_exit_codes

    implicit none

    contains

    !> check if CFL is ok
    subroutine COMPUTE_CFL(h_max, dt)

        !TODO find correct speed and constant for CFL
        implicit none

        real(kind=8), intent(in) :: h_max !< maximum grid diameter
        real(kind=8), intent(in) :: dt !< time-step

        real(kind=8) :: speed !< speed

        speed = 1

        if (dt > h_max/speed) then
            if(mpi_id==0) print *, "!-------- FAIL CFL CONDITION --------!"
            if (IS_failCFL .eqv. .true.) then

                call PetscFinalize(mpi_ierr)
                call MPI_FINALIZE(mpi_ierr)
                call EXIT(EXIT_CFL)
            endif

        endif

    end subroutine COMPUTE_CFL

end module checks