!    Copyright (C) 2021 The SPEED FOUNDATION
!    Author: Ilario Mazzieri

!> @brief definition of the struct `Data_Structure`, which stores the material
!> file parameters such as the density and Lamé parameters, and subroutines
!> related to print, allocation and setting of default values to the parameters.
module Poly_data

    use Poly_global
    use Poly_exit_codes, only: EXIT_NO_MATERIALS, EXIT_NOTHONORING_ERROR, &
                                EXIT_DAMPING_PEAK

    implicit none

    !> Material and loads parameters
    type Data_Structure
    !*******************************************************************************
    ! Material file parameters -
    !*******************************************************************************

        integer(kind=4) :: nmat             !< number of elastic material
        integer(kind=4) :: nmat_nle         !< number of nonlinear el. material
        integer(kind=4) :: nmat_rnd         !< number of random el. material

        integer(kind=4) :: nload_diri_el    !< number of el. Dirichlet bound. (XYZ load)
        integer(kind=4) :: nload_neum_el    !< number of el. Neumann bound. (XYZ load)
        integer(kind=4) :: nload_neuN_el    !< number of el. Neumann bound. (normal load)

        ! Number of element point load (XYZ load)
        integer(kind=4) :: nload_poiX_el    !< Number of element point load X
        integer(kind=4) :: nload_poiY_el    !< Number of element point load Y
        integer(kind=4) :: nload_poiZ_el    !< Number of element point load Z

        ! Number of element plane wave load (XYZ load)
        integer(kind=4) :: nload_plaX_el    !< Number of element plane wave load X
        integer(kind=4) :: nload_plaY_el    !< Number of element plane wave load Y
        integer(kind=4) :: nload_plaZ_el    !< Number of element plane wave load Z

        ! Number of element volume load (XYZ load)
        integer(kind=4) :: nload_forX_el    !< Number of element volumme load X
        integer(kind=4) :: nload_forY_el    !< Number of element volumme load Y
        integer(kind=4) :: nload_forZ_el    !< Number of element volumme load Z

        integer(kind=4) :: nload_abc_el     !< number of el. absorbing conditions
        integer(kind=4) :: nfunc            !< number of keyword FUNC
        integer(kind=4) :: nfunc_data       !< number of data for all FUNC keywords
        integer(kind=4) :: nload_sism_el    !< number of seismic load (kinematic source)
        integer(kind=4) :: n_case           !< CASE number for not-honoring
        integer(kind=4) :: nmat_nhe         !< material number for not-honoring enhanced
        integer(kind=4) :: srcmodflag       !<flag for srcmod - sism lines
        real(kind=8)    :: fmax             !<reference f-value for damping
        real(kind=8)    :: fpeak            !<peak frequency non linear case


        ! Integer arrays with allocatable dimensions
        integer(kind=4), dimension(:), allocatable :: sdeg_mat
        integer(kind=4), dimension(:), allocatable :: tag_mat
        integer(kind=4), dimension(:), allocatable :: sdeg_mat_nle
        integer(kind=4), dimension(:), allocatable :: tag_mat_nle
        integer(kind=4), dimension(:), allocatable :: rand_mat
        integer(kind=4), dimension(:), allocatable :: space_fun_tag_diri_el ! function in space tag for dirichlet BC in problem_data_and_properties.f90 - tag = Element%neigh_el(:,5)
        integer(kind=4), dimension(:), allocatable :: face_tag_diri_el       ! mesh face tag
        integer(kind=4), dimension(:), allocatable :: space_fun_tag_neum_el ! function in space tag for neumann BC in problem_data_and_properties.f90 - tag = Element%neigh_el(:,5)
        integer(kind=4), dimension(:), allocatable :: face_tag_neum_el       ! mesh face tag
        integer(kind=4), dimension(:), allocatable :: fun_neuN_el
        integer(kind=4), dimension(:), allocatable :: tag_neuN_el
        integer(kind=4), dimension(:), allocatable :: fun_poiX_el
        integer(kind=4), dimension(:), allocatable :: fun_poiY_el
        integer(kind=4), dimension(:), allocatable :: fun_poiZ_el
        integer(kind=4), dimension(:), allocatable :: fun_plaX_el
        integer(kind=4), dimension(:), allocatable :: tag_plaX_el
        integer(kind=4), dimension(:), allocatable :: fun_plaY_el
        integer(kind=4), dimension(:), allocatable :: tag_plaY_el
        integer(kind=4), dimension(:), allocatable :: fun_plaZ_el
        integer(kind=4), dimension(:), allocatable :: tag_plaZ_el
        integer(kind=4), dimension(:), allocatable :: fun_forX_el
        integer(kind=4), dimension(:), allocatable :: fun_forY_el
        integer(kind=4), dimension(:), allocatable :: fun_forZ_el
        integer(kind=4), dimension(:), allocatable :: tag_abc_el
        integer(kind=4), dimension(:), allocatable :: fun_sism_el
        integer(kind=4), dimension(:), allocatable :: tag_sism_el
        integer(kind=4), dimension(:), allocatable :: tag_case
        integer(kind=4), dimension(:), allocatable :: val_case
        integer(kind=4), dimension(:), allocatable :: val_nhe
        integer(kind=4), dimension(:), allocatable :: tol_nhe
        integer(kind=4), dimension(:), allocatable :: tag_func
        integer(kind=4), dimension(:), allocatable :: func_type
        integer(kind=4), dimension(:), allocatable :: func_indx

        ! Integer arrays with two-dimensional allocatable array
        integer(kind=4), dimension(:,:), allocatable :: tag_func_mat_nle

        ! Real arrays with allocatable dimensions
        real(kind=8), dimension(:), allocatable :: tol_case
        real(kind=8), dimension(:), allocatable :: func_data
        real(kind=8), dimension(:), allocatable :: QS
        real(kind=8), dimension(:), allocatable :: QP

        ! Real arrays with two-dimensional allocatable arrays
        real(kind=8), dimension(:,:), allocatable :: prop_mat
        real(kind=8), dimension(:,:), allocatable :: val_mat_nle
        real(kind=8), dimension(:,:), allocatable :: val_diri_el
        real(kind=8), dimension(:,:), allocatable :: val_neum_el
        real(kind=8), dimension(:,:), allocatable :: val_neuN_el
        real(kind=8), dimension(:,:), allocatable :: val_poiX_el
        real(kind=8), dimension(:,:), allocatable :: val_poiY_el
        real(kind=8), dimension(:,:), allocatable :: val_poiZ_el
        real(kind=8), dimension(:,:), allocatable :: val_plaX_el
        real(kind=8), dimension(:,:), allocatable :: val_plaY_el
        real(kind=8), dimension(:,:), allocatable :: val_plaZ_el
        real(kind=8), dimension(:,:), allocatable :: val_forX_el
        real(kind=8), dimension(:,:), allocatable :: val_forY_el
        real(kind=8), dimension(:,:), allocatable :: val_forZ_el
        real(kind=8), dimension(:,:), allocatable :: val_sism_el


    end type Data_Structure


    contains

    ! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Set material parameters to default values
    subroutine set_Data_Structure_default_value(PolyData)

        implicit none

        type(Data_Structure), intent(inout) :: PolyData !< materials and loads

        PolyData%nmat = 0             ! number of elastic material
        PolyData%nmat_nle = 0         ! number of nonlinear el. material
        PolyData%nmat_rnd = 0         ! number of random el. material

        PolyData%nload_diri_el = 0    ! number of el. Dirichlet bound. (XYZ load)

        PolyData%nload_neum_el = 0    ! number of el. Neumann bound. (XYZ load)

        PolyData%nload_neuN_el = 0    ! number of el. Neumann bound. (normal load)

        PolyData%nload_poiX_el = 0    ! number of el. point load (XYZ load)
        PolyData%nload_poiY_el = 0
        PolyData%nload_poiZ_el = 0

        PolyData%nload_plaX_el = 0     ! number of el. plane wave load (XYZ load)
        PolyData%nload_plaY_el = 0
        PolyData%nload_plaZ_el = 0

        PolyData%nload_forX_el = 0     ! number of el. volume load (XYZ load)
        PolyData%nload_forY_el = 0
        PolyData%nload_forZ_el = 0

        PolyData%nload_abc_el  = 0     ! number of el. absorbing conditions

        PolyData%nfunc         = 0     ! number of keyword FUNC

        PolyData%nfunc_data    = 0     ! number of data for all FUNC keywords

        PolyData%nload_sism_el = 0     ! number of seismic load (kinematic source)

        PolyData%n_case        = 0     ! CASE number for not-honoring

        PolyData%nmat_nhe      = 0     ! material number for not-honoring enhanced

        PolyData%srcmodflag    = 0     !f lag for srcmod - sism lines

    end subroutine set_Data_Structure_default_value

    ! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Print number of material parameters to terminal
    subroutine print_Dime_Data_Structure(PolyData)

        implicit none

        type(Data_Structure), intent(inout) :: PolyData !< materials and loads

        if (PolyData%nmat <= 0) then
            write(*,*)'Error ! Number of material is 0'
            call EXIT(EXIT_NO_MATERIALS)
        endif

        write(*,'(A,I8)') 'Materials        : ',PolyData%nmat

        if(PolyData%nmat_nle > 0) &
            write(*,'(A,I8)')     'Materials NL-El. : ',PolyData%nmat_nle
        if(PolyData%nmat_rnd > 0) &
            write(*,'(A,I8)')     'Materials Random : ',PolyData%nmat_rnd
        if(PolyData%nload_diri_el > 0) &
            write(*,'(A,I8)')'Dirichlet B.C.   : ',PolyData%nload_diri_el
        if(PolyData%nload_neum_el > 0) &
            write(*,'(A,I8)')'Neumann B.C.     : ',PolyData%nload_neum_el
        if(PolyData%nload_neuN_el > 0) &
            write(*,'(A,I8)')'Neumann N B.C.   : ',PolyData%nload_neuN_el
        if(PolyData%nload_poiX_el > 0) &
            write(*,'(A,I8)')'Point Loads X    : ',PolyData%nload_poiX_el
        if(PolyData%nload_poiY_el > 0) &
            write(*,'(A,I8)')'Point Loads Y    : ',PolyData%nload_poiY_el
        if(PolyData%nload_poiZ_el > 0) &
            write(*,'(A,I8)')'Point Loads Z    : ',PolyData%nload_poiZ_el
        if(PolyData%nload_plaX_el > 0) &
            write(*,'(A,I8)')'Plane Loads X    : ',PolyData%nload_plaX_el
        if(PolyData%nload_plaY_el > 0) &
            write(*,'(A,I8)')'Plane Loads Y    : ',PolyData%nload_plaY_el
        if(PolyData%nload_plaZ_el > 0) &
            write(*,'(A,I8)')'Plane Loads Z    : ',PolyData%nload_plaZ_el
        if(PolyData%nload_forX_el > 0) &
            write(*,'(A,I8)')'Force X          : ',PolyData%nload_forX_el
        if(PolyData%nload_forY_el > 0) &
            write(*,'(A,I8)')'Force Y          : ',PolyData%nload_forY_el
        if(PolyData%nload_forZ_el > 0) &
            write(*,'(A,I8)')'Force Z          : ',PolyData%nload_forZ_el
        if(PolyData%nload_abc_el > 0)  &
            write(*,'(A,I8)')'ABSO Boundaries  : ',PolyData%nload_abc_el
        if(PolyData%nfunc > 0) &
            write(*,'(A,I8)')'Functions        : ',PolyData%nfunc
        if(PolyData%srcmodflag == 1) &
            write(*,*)'Using Not Honoring Fault Method For Siesmic Sources'
        if(PolyData%nload_sism_el > 0) &
            write(*,'(A,I8)')'Seis. Mom. Load  : ',PolyData%nload_sism_el

        if (PolyData%n_case > 1) then
            write(*,'(A)')'CASE WARNING: More than one case defined,'
            write(*,'(A)')'              only the 1st case will be adopted'
        endif
        write(*,'(A,I8)') 'CASE             : ',PolyData%n_case

        if(PolyData%nmat_nhe > 0) &
            write(*,'(A,I8)')     'Not_Honoring Enhanced Blocks : ',PolyData%nmat_nhe


    end subroutine print_Dime_Data_Structure

    ! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Allocate material parameters
    subroutine allocate_Data_Structure(PolyData)

        implicit none

        type(Data_Structure), intent(inout) :: PolyData !< materials and loads
        integer(kind=4)  :: size_sism

        ! material parameters
        allocate (PolyData%prop_mat(PolyData%nmat,4), &
                    PolyData%sdeg_mat(PolyData%nmat), &
                    PolyData%tag_mat(PolyData%nmat))

        ! Quality factors for damping
        allocate(PolyData%QS(PolyData%nmat), PolyData%QP(PolyData%nmat))
        PolyData%QS = 0.d0
        PolyData%QP = 0.d0;


        ! non linear material
        if (PolyData%nmat_nle > 0) &
            allocate(PolyData%sdeg_mat_nle(PolyData%nmat_nle), &
                        PolyData%tag_func_mat_nle(PolyData%nmat_nle,1), &
                        PolyData%val_mat_nle(PolyData%nmat_nle,1), &
                        PolyData%tag_mat_nle(PolyData%nmat_nle))

        ! random material
        if (PolyData%nmat_rnd > 0) allocate(PolyData%rand_mat(PolyData%nmat_rnd))

        ! Dirichlet boundary conditions
        if (PolyData%nload_diri_el > 0) &
            allocate (PolyData%val_diri_el(PolyData%nload_diri_el,4), &
                        PolyData%space_fun_tag_diri_el(PolyData%nload_diri_el), &
                        PolyData%face_tag_diri_el(PolyData%nload_diri_el))

        ! Neumann boundary conditions
        if (PolyData%nload_neum_el > 0) &
            allocate (PolyData%val_neum_el(PolyData%nload_neum_el,4), &
                        PolyData%space_fun_tag_neum_el(PolyData%nload_neum_el), &
                        PolyData%face_tag_neum_el(PolyData%nload_neum_el))

        if (PolyData%nload_neuN_el > 0) &
            allocate (PolyData%val_neuN_el(PolyData%nload_neuN_el,4), &
                        PolyData%fun_neuN_el(PolyData%nload_neuN_el), &
                        PolyData%tag_neuN_el(PolyData%nload_neuN_el))

        ! Point Load
        if (PolyData%nload_poiX_el > 0) &
            allocate (PolyData%val_poiX_el(PolyData%nload_poiX_el,4), &
                        PolyData%fun_poiX_el(PolyData%nload_poiX_el))

        if (PolyData%nload_poiY_el > 0) &
            allocate (PolyData%val_poiY_el(PolyData%nload_poiY_el,4), &
                        PolyData%fun_poiY_el(PolyData%nload_poiY_el))

        if (PolyData%nload_poiZ_el > 0) &
            allocate (PolyData%val_poiZ_el(PolyData%nload_poiZ_el,4), &
                        PolyData%fun_poiZ_el(PolyData%nload_poiZ_el))

        ! Plane Wave
        if (PolyData%nload_plaX_el > 0) &
            allocate (PolyData%val_plaX_el(PolyData%nload_plaX_el,1), &
                        PolyData%fun_plaX_el(PolyData%nload_plaX_el), &
                        PolyData%tag_plaX_el(PolyData%nload_plaX_el))

        if (PolyData%nload_plaY_el > 0) &
            allocate (PolyData%val_plaY_el(PolyData%nload_plaY_el,1), &
                        PolyData%fun_plaY_el(PolyData%nload_plaY_el), &
                        PolyData%tag_plaY_el(PolyData%nload_plaY_el))

        if (PolyData%nload_plaZ_el > 0) &
            allocate (PolyData%val_plaZ_el(PolyData%nload_plaZ_el,1), &
                        PolyData%fun_plaZ_el(PolyData%nload_plaZ_el), &
                        PolyData%tag_plaZ_el(PolyData%nload_plaZ_el))

        ! Volume force
        if (PolyData%nload_forX_el > 0) &
            allocate (PolyData%val_forX_el(PolyData%nload_forX_el,4), &
                      PolyData%fun_forX_el(PolyData%nload_forX_el))

        if (PolyData%nload_forY_el > 0) &
            allocate (PolyData%val_forY_el(PolyData%nload_forY_el,4), &
                      PolyData%fun_forY_el(PolyData%nload_forY_el))

        if (PolyData%nload_forZ_el > 0) &
            allocate (PolyData%val_forZ_el(PolyData%nload_forZ_el,4), &
                      PolyData%fun_forZ_el(PolyData%nload_forZ_el))

        ! Absorbing boundaries
        if (PolyData%nload_abc_el > 0) &
            allocate (PolyData%tag_abc_el(PolyData%nload_abc_el))

        ! Slip mode
        if (PolyData%srcmodflag == 0) then
            size_sism = 21
            if (PolyData%nload_sism_el > 0) &
                allocate (PolyData%val_sism_el(PolyData%nload_sism_el,21), &
                            PolyData%fun_sism_el(PolyData%nload_sism_el), &
                            PolyData%tag_sism_el(PolyData%nload_sism_el))
        elseif (PolyData%srcmodflag == 1) then
            size_sism = 15
            if (PolyData%nload_sism_el > 0) &
                allocate (PolyData%val_sism_el(PolyData%nload_sism_el,15), &
                            PolyData%fun_sism_el(PolyData%nload_sism_el), &
                            PolyData%tag_sism_el(PolyData%nload_sism_el))
        endif

        ! Not honoring case
        if (PolyData%n_case > 0) &
            allocate (PolyData%val_case(PolyData%n_case), &
                        PolyData%tag_case(PolyData%n_case), &
                        PolyData%tol_case(PolyData%n_case))


        if (PolyData%n_case == 0) &
            allocate(PolyData%tag_case(1)); PolyData%tag_case(1) = 0

        ! Not honoring enhanced
        if (PolyData%nmat_nhe > 0) &
            allocate (PolyData%val_nhe(PolyData%nmat_nhe), &
                        PolyData%tol_nhe(PolyData%nmat_nhe))

        ! Functions
        if (PolyData%nfunc > 0) &
            allocate (PolyData%tag_func(PolyData%nfunc), &
                        PolyData%func_type(PolyData%nfunc), &
                        PolyData%func_indx(PolyData%nfunc +1), &
                        PolyData%func_data(PolyData%nfunc_data))


    end subroutine allocate_Data_Structure

    ! - >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !> Print to material and loading parameters terminal
    subroutine print_Data_Structure(PolyData)

        implicit none

        type(Data_Structure), intent(inout) :: PolyData   !< materials and loads
        integer(kind=4) :: i, im


        write(*,'(A)') 'Read.'
        write(*,'(A)')

        if (PolyData%n_case > 0) then
            do i = 1, PolyData%n_case
                write(*,'(A)') '------------------Not-Honoring case-------------------'
                write(*,'(A,I8)') 'CASE :', PolyData%tag_case(i)
                select case (PolyData%tag_case(i))
                    case(1); write(*,'(A)')'GRENOBLE HONORING'
                    case(2); write(*,'(A)')'GRENOBLE'
                    case(3); write(*,'(A)')'GUBBIO'
                    case(4); write(*,'(A)')'SULMONA'
                    case(5); write(*,'(A)')'VOLVI'
                    case(6); write(*,'(A)')'FRIULI'
                    case(7); write(*,'(A)')'AQUILA'
                    case(8); write(*,'(A)')'SANTIAGO'
                    case(11); write(*,'(A)')'CHRISTCHURCH NEW TOPO'
                    case(12); write(*,'(A)')'PO PLAIN (new model)'
                    case(13); write(*,'(A)')'PO PLAIN-BEDROCK (new-model)'
                    case(14); write(*,'(A)')'WELLINGTON (Benites)'
                    case(15); write(*,'(A)')'MARSICA-FUCINO'
                    case(16); write(*,'(A)')'ISTANBUL'
                    case(18); write(*,'(A)')'BEIJING-TUTORIAL'
                    case(19); write(*,'(A)')'THESSALONIKI'
                    case(20); write(*,'(A)')'ATHENS'
                    case(21); write(*,'(A)')'BEIJING'
                    case(22); write(*,'(A)')'NORCIA'
                    case(27); write(*,'(A)')'AQUILA-OB'
                    case(28); write(*,'(A)')'NORCIA-OB'
                    case(29); write(*,'(A)')'THESS-OB'
                    case(30); write(*,'(A)')'ATHENS-Parthenon'
                    case(31); write(*,'(A)')'GRONINGEN'
                    case(32); write(*,'(A)')'GRONINGEN-LAYERED MODEL'
                    case(33); write(*,'(A)')'GRONINGEN-ZE'
                    case(35); write(*,'(A)')'THESS+MYGD-FINAL'
                    case(38); write(*,'(A)')'MONTELIMAR'
                    case(40); write(*,'(A)')'KUTCH BASIN, INDIA'
                    case(46); write(*,'(A)')'KUMAMOTO, JAPAN'
                    case(60); write(*,'(A)')'JAKARTA, INDONESIA'
                    case(70); write(*,'(A)')'AQUILA MULTI-BASIN'
                    case(98,99,100); write(*,'(A)')'TEST MODE'
                    case default
                        write(*,'(A)')'ERROR.. this case was not implemented!'
                        call EXIT(EXIT_NOTHONORING_ERROR)
                end select

                write(*,'(A,I8)') 'MATERIAL TYPE    : ',PolyData%val_case(i)
                write(*,'(A,E12.4)') 'MATERIAL TOL.    : ',PolyData%tol_case(i)
                write(*,'(A)')
            enddo
        endif

        write(*,'(A)')
        do im = 1, PolyData%nmat
            write(*,'(A,I8)')    'MATERIAL#   : ',PolyData%tag_mat(im)
            write(*,'(A,I8)')    'degree      : ',PolyData%sdeg_mat(im)
            write(*,'(A,E12.4)') 'rho [kg/m3] : ',PolyData%prop_mat(im,1)
            write(*,'(A,E12.4)') 'vs  [m/s]   : ', &
                    (PolyData%prop_mat(im,3)/PolyData%prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'vp  [m/s]   : ', &
                    ((PolyData%prop_mat(im,2) + 2*PolyData%prop_mat(im,3))&
                      / PolyData%prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'zeta [1/s]  : ',PolyData%prop_mat(im,4)
            write(*,'(A,E12.4)') 'Qs [-]      : ',PolyData%QS(im)
            write(*,'(A,E12.4)') 'Qp [-]      : ',PolyData%QP(im)
            write(*,*)
        enddo

        write(*,'(A)')
        do im = 1, PolyData%nmat_nle
            write(*,'(A,I8)')    'NON LINEAR ELASTIC MATERIAL# : ', &
                PolyData%tag_mat_nle(im)
            write(*,'(A,I8)')    'Tag func nle : ',(PolyData%tag_func_mat_nle(im,1))
            write(*,'(A,E12.4)') 'Depth  [m]   : ',(PolyData%val_mat_nle(im,1))
            write(*,'(A,I8)')    'Degree       : ', PolyData%sdeg_mat_nle(im)
            write(*,'(A)')
        enddo

        if ((PolyData%fpeak == 0).and. (PolyData%nmat_nle > 0)) then
            write(*,'(A)')'ERROR: Peak frequency for damping not defined!'
            call EXIT(EXIT_DAMPING_PEAK)
        endif

    end subroutine print_Data_Structure

end module Poly_data