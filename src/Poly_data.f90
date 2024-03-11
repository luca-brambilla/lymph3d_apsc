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
 
     module Poly_data
     
     use Poly_global
     use Poly_exit_codes, only: EXIT_NO_MATERIALS, EXIT_NOTHONORING_ERROR, &
                                    EXIT_DAMPING_PEAK 
     
     implicit none 
     
     type Data_Structure 
     !*******************************************************************************
     ! Material file parameters - 
     !*******************************************************************************
        
        integer(kind=4) :: nmat             !< number of elastic material
        integer(kind=4) :: nmat_nle         !< number of nonlinear el. material
        integer(kind=4) :: nmat_rnd         !< number of random el. material     

        integer(kind=4) :: nload_diri_el
                                             !< number of el. Dirichlet bound. (XYZ load) 
        integer(kind=4) :: nload_neum_el
                                             !< number of el. Neumann bound. (XYZ load)
        integer(kind=4) :: nload_neuN_el    !< number of el. Neumann bound. (normal load)    
        integer(kind=4) :: nload_poiX_el,nload_poiY_el,nload_poiZ_el
                                             !< number of el. point load (XYZ load)
        integer(kind=4) :: nload_plaX_el,nload_plaY_el,nload_plaZ_el       
                                             !< number of el. plane wave load (XYZ load)
        integer(kind=4) :: nload_forX_el,nload_forY_el,nload_forZ_el
                                             !< number of el. volume load (XYZ load)
        integer(kind=4) :: nload_abc_el     !< number of el. absorbing conditions
        integer(kind=4) :: nfunc            !< number of keyword FUNC
        integer(kind=4) :: nfunc_data       !< number of data for all FUNC keywords
        integer(kind=4) :: nload_sism_el    !< number of seismic load (kinematic source)
        integer(kind=4) :: n_case           !< CASE number for not-honoring
        integer(kind=4) :: nmat_nhe         !< material number for not-honoring enhanced 
        integer(kind=4) :: srcmodflag       !<flag for srcmod - sism lines
        real(kind=8)    :: fmax             !<reference f-value for damping 
        real(kind=8)    :: fpeak            !<peak frequency non linear case

        integer(kind=4), dimension(:), allocatable :: sdeg_mat, tag_mat, &
                                                      sdeg_mat_nle, tag_mat_nle, &
                                                      rand_mat, &
                                                      fun_space_diri_el, tag_diri_el, &
                                                      fun_space_neum_el, tag_neum_el, &
                                                      fun_neuN_el, tag_neuN_el, &
                                                      fun_poiX_el, fun_poiY_el, &
                                                      fun_poiZ_el , &
                                                      fun_plaX_el, tag_plaX_el, &
                                                      fun_plaY_el, tag_plaY_el, &
                                                      fun_plaZ_el, tag_plaZ_el, &
                                                      fun_forX_el, fun_forY_el, &
                                                      fun_forZ_el, &
                                                      tag_abc_el, &
                                                      fun_sism_el, tag_sism_el, &
                                                      tag_case, val_case, &
                                                      val_nhe, tol_nhe, &
                                                      tag_func, func_type, func_indx
                                                      
                                                      
                                                      
                                                       
        integer(kind=4), dimension(:,:), allocatable :: tag_func_mat_nle

        real(kind=8), dimension(:), allocatable :: tol_case, func_data, &
                                                   QS, QP
        
        real(kind=8), dimension(:,:), allocatable :: prop_mat, val_mat_nle, &
                                                     val_diri_el,  &
                                                     val_neum_el, &
                                                     val_neuN_el, &
                                                     val_poiX_el, val_poiY_el, &
                                                     val_poiZ_el, val_plaX_el, &
                                                     val_plaY_el, val_plaZ_el, &
                                                     val_forX_el, val_forY_el, &
                                                     val_forZ_el, &
                                                     val_sism_el


     end type Data_Structure

    
     contains 
     
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        subroutine set_Data_Structure_default_value(Struct)
        
        implicit none

        type(Data_Structure), intent(inout) :: Struct
                
         Struct%nmat = 0             !< number of elastic material
         Struct%nmat_nle = 0         !< number of nonlinear el. material
         Struct%nmat_rnd = 0         !< number of random el. material                   

         Struct%nload_diri_el = 0    !< number of el. Dirichlet bound. (XYZ load) 
                                 
         Struct%nload_neum_el = 0    !< number of el. Neumann bound. (XYZ load)
                                 
         Struct%nload_neuN_el = 0    !< number of el. Neumann bound. (normal load) 

         Struct%nload_poiX_el = 0    !< number of el. point load (XYZ load)
         Struct%nload_poiY_el = 0
         Struct%nload_poiZ_el = 0
                                     
         Struct%nload_plaX_el = 0     !< number of el. plane wave load (XYZ load)
         Struct%nload_plaY_el = 0
         Struct%nload_plaZ_el = 0                               
                                 
         Struct%nload_forX_el = 0     !< number of el. volume load (XYZ load)
         Struct%nload_forY_el = 0
         Struct%nload_forZ_el = 0
                                      
         Struct%nload_abc_el  = 0     !< number of el. absorbing conditions
                             
         Struct%nfunc         = 0     !< number of keyword FUNC
         
         Struct%nfunc_data    = 0     !< number of data for all FUNC keywords
                            
         Struct%nload_sism_el = 0     !< number of seismic load (kinematic source) 

         Struct%n_case        = 0     !< CASE number for not-honoring
         
         Struct%nmat_nhe      = 0     !< material number for not-honoring enhanced 
         
         Struct%srcmodflag    = 0     !<flag for srcmod - sism lines        
        
        end subroutine set_Data_Structure_default_value
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
        subroutine print_Dime_Data_Structure(Struct)
        
      
        implicit none
        
        type(Data_Structure), intent(inout) :: Struct
        
        if (Struct%nmat <= 0) then
          write(*,*)'Error ! Number of material is 0'
          call EXIT(EXIT_NO_MATERIALS)
        endif

        write(*,'(A,I8)') 'Materials        : ',Struct%nmat
        
        if(Struct%nmat_nle > 0) &
            write(*,'(A,I8)')     'Materials NL-El. : ',Struct%nmat_nle        
        if(Struct%nmat_rnd > 0) &
            write(*,'(A,I8)')     'Materials Random : ',Struct%nmat_rnd        
        if(Struct%nload_diri_el > 0) &
            write(*,'(A,I8)')'Dirichlet B.C. : ',Struct%nload_diri_el
        if(Struct%nload_neum_el > 0) &
            write(*,'(A,I8)')'Neumann B.C.   : ',Struct%nload_neum_el
        if(Struct%nload_neuN_el > 0) &
            write(*,'(A,I8)')'Neumann N B.C.   : ',Struct%nload_neuN_el              
        if(Struct%nload_poiX_el > 0) &
            write(*,'(A,I8)')'Point Loads X    : ',Struct%nload_poiX_el
        if(Struct%nload_poiY_el > 0) &
            write(*,'(A,I8)')'Point Loads Y    : ',Struct%nload_poiY_el
        if(Struct%nload_poiZ_el > 0) &
            write(*,'(A,I8)')'Point Loads Z    : ',Struct%nload_poiZ_el
        if(Struct%nload_plaX_el > 0) &
            write(*,'(A,I8)')'Plane Loads X    : ',Struct%nload_plaX_el              
        if(Struct%nload_plaY_el > 0) &
            write(*,'(A,I8)')'Plane Loads Y    : ',Struct%nload_plaY_el              
        if(Struct%nload_plaZ_el > 0) &
            write(*,'(A,I8)')'Plane Loads Z    : ',Struct%nload_plaZ_el              
        if(Struct%nload_forX_el > 0) &
            write(*,'(A,I8)')'Force X          : ',Struct%nload_forX_el
        if(Struct%nload_forY_el > 0) &
            write(*,'(A,I8)')'Force Y          : ',Struct%nload_forY_el
        if(Struct%nload_forZ_el > 0) &
            write(*,'(A,I8)')'Force Z          : ',Struct%nload_forZ_el
        if(Struct%nload_abc_el > 0)  &
            write(*,'(A,I8)')'ABSO Boundaries  : ',Struct%nload_abc_el
        if(Struct%nfunc > 0) &
            write(*,'(A,I8)')'Functions        : ',Struct%nfunc
        if(Struct%srcmodflag == 1) &
            write(*,*)'Using Not Honoring Fault Method For Siesmic Sources'
        if(Struct%nload_sism_el > 0) &
            write(*,'(A,I8)')'Seis. Mom. Load  : ',Struct%nload_sism_el        

         if (Struct%n_case > 1) then
            write(*,'(A)')'CASE WARNING: More than one case defined,'        
            write(*,'(A)')'              only the 1st case will be adopted'  
         endif
         write(*,'(A,I8)') 'CASE             : ',Struct%n_case                               

         if(Struct%nmat_nhe > 0) &
            write(*,'(A,I8)')     'Not_Honoring Enhanced Blocks : ',Struct%nmat_nhe        
          
        
        end subroutine print_Dime_Data_Structure
        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
        subroutine allocate_Data_Structure(Struct)
        
              
        implicit none
        
        type(Data_Structure), intent(inout) :: Struct
        integer(kind=4)  :: size_sism
      
        ! material parameters
        allocate (Struct%prop_mat(Struct%nmat,4), &
                  Struct%sdeg_mat(Struct%nmat), &
                  Struct%tag_mat(Struct%nmat))
        
        ! Quality factors for damping
        allocate(Struct%QS(Struct%nmat), Struct%QP(Struct%nmat));
        Struct%QS = 0.d0; Struct%QP = 0.d0;


        ! non linear material
        if (Struct%nmat_nle > 0) &
            allocate(Struct%sdeg_mat_nle(Struct%nmat_nle), &
                     Struct%tag_func_mat_nle(Struct%nmat_nle,1), &
                     Struct%val_mat_nle(Struct%nmat_nle,1), &
                     Struct%tag_mat_nle(Struct%nmat_nle))
        
        ! random material             
        if (Struct%nmat_rnd > 0) allocate(Struct%rand_mat(Struct%nmat_rnd))
      
        ! Dirichlet boundary conditions
        if (Struct%nload_diri_el > 0) &
            allocate (Struct%val_diri_el(Struct%nload_diri_el,4), &
                      Struct%fun_space_diri_el(Struct%nload_diri_el), &
                      Struct%tag_diri_el(Struct%nload_diri_el)) 
      
        ! Neumann boundary conditions
        if (Struct%nload_neum_el > 0) &
            allocate (Struct%val_neum_el(Struct%nload_neum_el,4), &
                      Struct%fun_space_neum_el(Struct%nload_neum_el), &
                      Struct%tag_neum_el(Struct%nload_neum_el))
                      
        if (Struct%nload_neuN_el > 0) &
            allocate (Struct%val_neuN_el(Struct%nload_neuN_el,4), &
                      Struct%fun_neuN_el(Struct%nload_neuN_el), &
                      Struct%tag_neuN_el(Struct%nload_neuN_el))
      
        ! Point Load
        if (Struct%nload_poiX_el > 0) &
            allocate (Struct%val_poiX_el(Struct%nload_poiX_el,4), &
                      Struct%fun_poiX_el(Struct%nload_poiX_el))
                          
        if (Struct%nload_poiY_el > 0) &
            allocate (Struct%val_poiY_el(Struct%nload_poiY_el,4), &
                      Struct%fun_poiY_el(Struct%nload_poiY_el)) 
                        
        if (Struct%nload_poiZ_el > 0) &
            allocate (Struct%val_poiZ_el(Struct%nload_poiZ_el,4), &
                      Struct%fun_poiZ_el(Struct%nload_poiZ_el))

        ! Plane Wave
        if (Struct%nload_plaX_el > 0) &
            allocate (Struct%val_plaX_el(Struct%nload_plaX_el,1), &
                      Struct%fun_plaX_el(Struct%nload_plaX_el), &
                      Struct%tag_plaX_el(Struct%nload_plaX_el))
                      
        if (Struct%nload_plaY_el > 0) &
            allocate (Struct%val_plaY_el(Struct%nload_plaY_el,1), &
                      Struct%fun_plaY_el(Struct%nload_plaY_el), &
                      Struct%tag_plaY_el(Struct%nload_plaY_el))
                      
        if (Struct%nload_plaZ_el > 0) &
            allocate (Struct%val_plaZ_el(Struct%nload_plaZ_el,1), &
                      Struct%fun_plaZ_el(Struct%nload_plaZ_el), &
                      Struct%tag_plaZ_el(Struct%nload_plaZ_el))
        
        ! Volume force
        if (Struct%nload_forX_el > 0) &
            allocate (Struct%val_forX_el(Struct%nload_forX_el,4), &
                      Struct%fun_forX_el(Struct%nload_forX_el))
                      
        if (Struct%nload_forY_el > 0) &
            allocate (Struct%val_forY_el(Struct%nload_forY_el,4), &
                      Struct%fun_forY_el(Struct%nload_forY_el))
                      
        if (Struct%nload_forZ_el > 0) &
            allocate (Struct%val_forZ_el(Struct%nload_forZ_el,4), &
                      Struct%fun_forZ_el(Struct%nload_forZ_el))
     
        ! Absorbing boundaries                                           
        if (Struct%nload_abc_el > 0) &
            allocate (Struct%tag_abc_el(Struct%nload_abc_el))

        ! Slip mode
        if (Struct%srcmodflag == 0) then
          size_sism = 21
          if (Struct%nload_sism_el > 0) &
              allocate (Struct%val_sism_el(Struct%nload_sism_el,21), &
                        Struct%fun_sism_el(Struct%nload_sism_el), &
                        Struct%tag_sism_el(Struct%nload_sism_el))
        elseif (Struct%srcmodflag == 1) then
          size_sism = 15
          if (Struct%nload_sism_el > 0) &
              allocate (Struct%val_sism_el(Struct%nload_sism_el,15), &
                        Struct%fun_sism_el(Struct%nload_sism_el), &
                        Struct%tag_sism_el(Struct%nload_sism_el))
        endif

        ! Not honoring case
        if (Struct%n_case > 0) &
            allocate (Struct%val_case(Struct%n_case), &
                      Struct%tag_case(Struct%n_case), &
                      Struct%tol_case(Struct%n_case))
        
                      
        if (Struct%n_case == 0) &
            allocate(Struct%tag_case(1)); Struct%tag_case(1) = 0;
         
        ! Not honoring enhanced 
        if (Struct%nmat_nhe > 0) &
            allocate (Struct%val_nhe(Struct%nmat_nhe), &
                      Struct%tol_nhe(Struct%nmat_nhe))
      
        ! Functions 
        if (Struct%nfunc > 0) &
            allocate (Struct%tag_func(Struct%nfunc), &
                      Struct%func_type(Struct%nfunc), &
                      Struct%func_indx(Struct%nfunc +1), &
                      Struct%func_data(Struct%nfunc_data))
      

        end subroutine allocate_Data_Structure
     

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
        subroutine print_Data_Structure(Struct)
        
        
        implicit none
        
        type(Data_Structure), intent(inout) :: Struct
        integer(kind=4) :: i, im


        write(*,'(A)') 'Read.'
        write(*,'(A)')                                

        if (Struct%n_case > 0) then
          do i = 1, Struct%n_case
            write(*,'(A)') '------------------Not-Honoring case-------------------'
              write(*,'(A,I8)') 'CASE :', Struct%tag_case(i) 
            select case (Struct%tag_case(i))  
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
                                   
             write(*,'(A,I8)') 'MATERIAL TYPE    : ',Struct%val_case(i)          
             write(*,'(A,E12.4)') 'MATERIAL TOL.    : ',Struct%tol_case(i)
             write(*,'(A)')
           enddo
        endif
        
        write(*,'(A)')
        do im = 1, Struct%nmat    
            write(*,'(A,I8)')    'MATERIAL#   : ',Struct%tag_mat(im)
            write(*,'(A,I8)')    'degree      : ',Struct%sdeg_mat(im)
            write(*,'(A,E12.4)') 'rho [kg/m3] : ',Struct%prop_mat(im,1)
            write(*,'(A,E12.4)') 'vs  [m/s]   : ', &
                    (Struct%prop_mat(im,3)/Struct%prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'vp  [m/s]   : ', &
                    ((Struct%prop_mat(im,2) + 2*Struct%prop_mat(im,3))&
                      / Struct%prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'zeta [1/s]  : ',Struct%prop_mat(im,4)
            write(*,'(A,E12.4)') 'Qs [-]      : ',Struct%QS(im)
            write(*,'(A,E12.4)') 'Qp [-]      : ',Struct%QP(im)
            write(*,*)
         enddo

         write(*,'(A)') 
         do im = 1, Struct%nmat_nle
           write(*,'(A,I8)')    'NON LINEAR ELASTIC MATERIAL# : ', &
                Struct%tag_mat_nle(im)
           write(*,'(A,I8)')    'Tag func nle : ',(Struct%tag_func_mat_nle(im,1))    
           write(*,'(A,E12.4)') 'Depth  [m]   : ',(Struct%val_mat_nle(im,1))          
           write(*,'(A,I8)')    'Degree       : ', Struct%sdeg_mat_nle(im)                 
           write(*,'(A)') 
         enddo
                                                                                                           
         if ((Struct%fpeak == 0).and. (Struct%nmat_nle > 0)) then   
           write(*,'(A)')'ERROR: Peak frequency for damping not defined!'                
           call EXIT(EXIT_DAMPING_PEAK)
         endif                                                                  
                                                                                              
        end subroutine print_Data_Structure
     
     end module Poly_data  

