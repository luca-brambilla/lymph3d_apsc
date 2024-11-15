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
 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    READ_HEADER(header_file)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine READ_HEADER(header_file)

   use Poly_global, only: grid_file, mate_file, &
                              folder_mpi, folder_monitors, folder_restart, &
                              opt_out_var, damping_type, &
                              time_step, start_time, stop_time, time_restart, &
                              num_dt_mon, Is_Restart, Is_Debug, &
                              depth_search_mon_lst, IS_mon_lst, &
                              IsTime_dependent, IsSave_output, IS_MatrixFree


   use Poly_exit_codes
   use Poly_fail_codes
   use Poly_default_codes
   use Poly_setup_mpi

   implicit none
   
   character(len=14) :: header_file
   character(len=8)  :: keyword
   character(len=70) :: inline

   integer(kind=4)   :: i, status, ileft, iright, arglen, val_mon_lst
   integer(kind=4)   :: file_row = 0

   integer(kind=4) :: IS_dynamic, IS_save, IS_free


   !Setup default values 
   start_time              = start_time_default
   damping_type            = damping_type_default

   IS_mon_lst              = IS_mon_lst_default
   IS_debug                = IS_debug_default
   IS_restart              = IS_restart_default
   IS_setuponly            = IS_setuponly_default
   IS_failoncoeffs         = IS_failoncoeffs_default
   IS_failCFL              = IS_failCFL_default
   IS_instabilitycontrol   = IS_instabilitycontrol_default
   IsTime_dependent        = IS_timedependent_default
   IsSave_output           = IS_saveoutput_default
   IS_MatrixFree           = IS_MatrixFree_default

   if(mpi_id == 0) write(*,'(A)')  !! PRINT VARIABLES
   
   open(40,file=header_file)
   
   do    
      read(40,'(A)',IOSTAT = status) inline
      file_row = file_row + 1
      
      if (status.ne.0) exit
      
      ! Skip comments
      if (inline(1:1) == ' ') then
         cycle
      endif

      !!!! Parse keyword arguments
      ileft = 1
      iright = len_trim(inline)
      
      ! Compute index to first non-keyword argument
      do i = 1,iright
         if (inline(i:i) == ' ') exit
      enddo
      ileft = i + 1
      keyword = inline(1:(ileft-2))

      ! Now ileft points after the first blank
      arglen = len_trim(inline(ileft:iright))

      ! Remove this for keywords without arguments!
      if (arglen == 0) then
         write(*,'(A,I3,A,A,A)') 'FATAL in PolyWAVE.input, row', &
                  file_row, ': no argument given for command "', trim(keyword), '"'
         call EXIT(EXIT_SYNTAX_ERROR)
      endif

!         write(*,*) 'Comparing ', inline(1:(ileft-2)), ', ileft=', ileft

      if(mpi_id == 0) write(*,'(A)')keyword  !! PRINT READ VARIABLES
      
      select case (keyword)

         case('GRIDFILE')
            read(inline(ileft:iright),*) grid_file

         case('MATFILE')
            read(inline(ileft:iright),*) mate_file
         
         case('MPIFILE')
            read(inline(ileft:iright),*) folder_mpi   
         
         case('MONFILE')
            read(inline(ileft:iright),*) folder_monitors   
         
         case('BKPFILE')
            read(inline(ileft:iright),*) folder_restart   

         case('OPTIOUT')
            read(inline(ileft:iright),*) opt_out_var(1),opt_out_var(2),opt_out_var(3), &
                                       opt_out_var(4),opt_out_var(5),opt_out_var(6) 
         
         case('TIMESTEP')
            read(inline(ileft:iright),*) time_step
         
         case('STARTIME')
            read(inline(ileft:iright),*) start_time
         
         case('STOPTIME')
            read(inline(ileft:iright),*) stop_time

         case('DYNAMIC')
            read(inline(ileft:iright),*) IS_dynamic
            if (IS_dynamic /= 0) IsTime_dependent = .true.

         case('RESTART')
            read(inline(ileft:iright),*) time_restart
            IS_Restart = .true.

         case('TMONITOR')                                        
            read(inline(ileft:iright),*) num_dt_mon  
            
         case('SAVEOUT')
            read(inline(ileft:iright),*) IS_save
            if(IS_save /= 0) IsSave_output = .true.

         case('MATFREE')
            read(inline(ileft:iright),*) IS_save
            if(IS_free /= 0) IS_MatrixFree = .true.
            
         case('DAMPING')
            read(inline(ileft:iright),*) damping_type   
   
         case('MLST')                        
            read(inline(ileft:iright),*) depth_search_mon_lst, val_mon_lst
            if (val_mon_lst == 1) IS_mon_lst = .true. 
            
            
         case('FAILCFL')
            ! If specified as "FAILCFL", quit if CFL condition does not hold
            IS_failCFL = .true.

         case('FAILINST')
            ! If specified as "FAILINST", enable instability control
            IS_instabilitycontrol = .true.

         case('SETUPONL')
            ! If specified as "SETUPONL", quit before starting the time loop
            IS_setuponly = .true.

         case('FAILCOEF')
            ! If specified as "FAILCOEF", quit if any of the computed
            ! anelastic coefficients is negative [damping 2]
            IS_failoncoeffs = .true.

         case('DEBUG')
            IS_Debug = .true.
      

         ! Fail if keyword is not recognised
         case default
            write(*,'(A,I3,A,A)') 'FATAL in test1.input, row', &
                        file_row, ': unknown keyword ', keyword
            call EXIT(EXIT_SYNTAX_ERROR)

      end select  

   enddo
         
   
   close(40)
   
   if(mpi_id == 0) write(*,'(A)') !! PRINT VARIABLES

end subroutine READ_HEADER


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    READ_DIME_MATEFILE(mate_file)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


     subroutine READ_DIME_MATEFILE(filemate,PolyData)
      
     use Poly_data     
      
     implicit none
      
     character(len=50)     :: filemate, src_name                   
     character(len=100000) :: inline
     character(len=4)      :: keyword
      
     integer(kind=4) :: status
     integer(kind=4) :: tagel_func, func_typec, nfunc_datac
      
     type(Data_Structure), intent(inout) :: PolyData
                                                                                            

     open(40,file=filemate)

     do
        read(40,'(A)',IOSTAT = status) inline

     if (status.ne.0) exit

       keyword = inline(1:4)

       select case (keyword)

         case('MATE'); PolyData%nmat          = PolyData%nmat + 1
         case('MATN'); PolyData%nmat_nle      = PolyData%nmat_nle + 1     
         case('MATR'); PolyData%nmat_rnd      = PolyData%nmat_rnd + 1
         case('DIRI'); PolyData%nload_diri_el = PolyData%nload_diri_el + 1
         case('NEUM'); PolyData%nload_neum_el = PolyData%nload_neum_el + 1
         case('NEUN'); PolyData%nload_neuN_el = PolyData%nload_neuN_el + 1                 
         case('PLOX'); PolyData%nload_poiX_el = PolyData%nload_poiX_el + 1
         case('PLOY'); PolyData%nload_poiY_el = PolyData%nload_poiY_el + 1
         case('PLOZ'); PolyData%nload_poiZ_el = PolyData%nload_poiZ_el + 1
         case('PLAX'); PolyData%nload_plaX_el = PolyData%nload_plaX_el + 1                
         case('PLAY'); PolyData%nload_plaY_el = PolyData%nload_plaY_el + 1                
         case('PLAZ'); PolyData%nload_plaZ_el = PolyData%nload_plaZ_el + 1                
         case('FORX'); PolyData%nload_forX_el = PolyData%nload_forX_el + 1
         case('FORY'); PolyData%nload_forY_el = PolyData%nload_forY_el + 1
         case('FORZ'); PolyData%nload_forZ_el = PolyData%nload_forZ_el + 1
         case('ABSO'); PolyData%nload_abc_el  = PolyData%nload_abc_el + 1
         case('SISM'); PolyData%nload_sism_el = PolyData%nload_sism_el + 1                
         case('CASE'); PolyData%n_case        = PolyData%n_case + 1        
         case('NHEE'); PolyData%nmat_nhe      = PolyData%nmat_nhe + 1        
         case('SLIP')        

            ! Not-honoring Fault Plane
            read(inline(5:),*) src_name
            if (src_name.eq.'LOAD-SRCMOD2') PolyData%srcmodflag = 1; 
           
           case('FUNC')
            PolyData%nfunc = PolyData%nfunc + 1         
            read(inline(5:),*) tagel_func, func_typec
            
            select case (func_typec)
               case(0,32)
                  ! Case 32 - Ramp Source Time Function
                  PolyData%nfunc_data = PolyData%nfunc_data + 0
               case(1) 
                  ! RICKER WAVELET
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(2) 
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(3,30,31,33) 
                  ! TIME SERIES
                  ! Case 31 - Text File with Source Time Function
                  read(inline(5:),*) tagel_func, func_typec, nfunc_datac
                  PolyData%nfunc_data = PolyData%nfunc_data + 2*nfunc_datac
                
               case(4) 
                  ! DERIVATIVE OF THE RICKER WAVELET
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(5) 
                  ! DERIVATIVE OF THE GAUSSIAN WAVELET
                  PolyData%nfunc_data = PolyData%nfunc_data + 2

               case(6) 
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(7) 
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(8,9) 
                  PolyData%nfunc_data = PolyData%nfunc_data + 1

               case(12) 
                  ! SIGMOIDAL FUNC
                  PolyData%nfunc_data = PolyData%nfunc_data + 3
               case(13) 
                  ! GRENOBLE BENCHMARK
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(14) 
                  ! SCEC BENCHMARK
                  PolyData%nfunc_data = PolyData%nfunc_data + 2
               case(15) 
                  ! EXPLOSION    
                  PolyData%nfunc_data = PolyData%nfunc_data + 4
               case(50,55) 
                  ! VARIABLE TAU 
                  PolyData%nfunc_data = PolyData%nfunc_data + 2  
               case(60,62) 
                  ! LINEAR EQUIVALENT
                  read(inline(5:),*) tagel_func, func_typec, nfunc_datac               
                  PolyData%nfunc_data = PolyData%nfunc_data + 2*nfunc_datac
               case(61,63)                                      
                  read(inline(5:),*) tagel_func, func_typec, nfunc_datac               
                  PolyData%nfunc_data = PolyData%nfunc_data + 2*nfunc_datac
               case(99) 
                  ! CASHIMA   
                  PolyData%nfunc_data = PolyData%nfunc_data + 2    
               case(100)
                 PolyData%nfunc_data = PolyData%nfunc_data + 1   
               case(773) 
                  ! TIME SERIES
                  read(inline(5:),*) tagel_func, func_typec, nfunc_datac
                  PolyData%nfunc_data = PolyData%nfunc_data + nfunc_datac

            end select 
         
         end select

     enddo
      
     close(40)
      
     end subroutine READ_DIME_MATEFILE


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    READ_DIME_MATEFILE(mate_file)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


     subroutine READ_MATEFILE(filemate,PolyData)
      
     use mpi
     use Poly_setup_mpi
     use Poly_data
     use Poly_exit_codes, only: EXIT_FUNCTION_ERROR   
     use Poly_global, only: damping_type   
      
     implicit none
      
     character(len=50)     :: filemate, fileinput                   
     character(len=100000) :: inline
     character(len=4)      :: keyword
      
     integer(kind=4) :: status
     integer(kind=4) :: tagel_func, func_typec, nfunc_datac, ileft, iright
     integer(kind=4) :: im, im_nle, im_rnd, icase, &
                        idX, idY, idZ, inX, inY, inZ, inN, &
                        ipX, ipY, ipZ, iplX, iplY, iplZ, &
                        ifX, ifY, ifZ, iabc, inhee, isism, iexpl, ifunc, &
                        file_nd, dummy, i, j
     
     real(kind=8)    :: rho,vs,vp                   
      
     type(Data_Structure), intent(inout) :: PolyData

     im    = 0;  im_nle = 0;  icase = 0;  im_rnd = 0;
     idX   = 0;  idY    = 0;  idZ   = 0
     inX   = 0;  inY    = 0;  inZ   = 0;  inN    = 0;
     ipX   = 0;  ipY    = 0;  ipZ   = 0
     iplX  = 0;  iplY   = 0;  iplZ  = 0 
     ifX   = 0;  ifY    = 0;  ifZ   = 0
     iabc  = 0;  inhee  = 0;
     isism = 0;  iexpl  = 0;  
     ifunc = 0

     PolyData%fmax = 0.d0
     


     open(40,file=filemate)

     if (PolyData%nfunc > 0) PolyData%func_indx(1) = 1
      
      
     do 
         read(40,'(A)',IOSTAT = status) inline
         
         
         if (status /= 0) exit
         
         keyword = inline(1:4)
         
         ileft = 0
         iright = len(inline)
         do i = 1,iright
            if (inline(i:i) == ' ') exit
         enddo
         ileft = i
         

         select case (keyword)

         
           case('MATE')
             im = im + 1
             read(inline(ileft:iright),*) PolyData%tag_mat(im), &
                                          PolyData%sdeg_mat(im),&
                                          rho, vs, vp, & 
                                          PolyData%QS(im), &
                                          PolyData%QP(im)

             ! Constant Quality factor is divided by 2 to obtain coherent results
             ! this aspect should be investigated
             if(damping_type == 2) then 
                PolyData%QS(im) = 0.5d0*PolyData%QS(im)
                PolyData%QP(im) = 0.5d0*PolyData%QP(im)
             endif
            !rho     
            PolyData%prop_mat(im,1) = rho
            !lambda
            PolyData%prop_mat(im,2) = rho * (VP**2 - 2*VS**2) 
            !mu 
            PolyData%prop_mat(im,3) = rho * VS**2
  
           case('MATN')                                                                        
             im_nle = im_nle + 1                                                          
             read(inline(ileft:iright),*) PolyData%tag_mat_nle(im_nle), &
                                          PolyData%sdeg_mat_nle(im_nle), & 
                                          PolyData%tag_func_mat_nle(im_nle,1), &
                                          PolyData%val_mat_nle(im_nle,1) 
         
           case('MATR')                                                                        
            im_rnd = im_rnd + 1                   
            read(inline(ileft:iright),*) PolyData%rand_mat(im_rnd)                          

           case('DIRI')
            idX = idX + 1
            read(inline(ileft:iright),*) PolyData%face_tag_diri_el(idX), &
                                         PolyData%space_fun_tag_diri_el(idX), &
                                         PolyData%val_diri_el(idX,1), &
                                         PolyData%val_diri_el(idX,2), &
                                         PolyData%val_diri_el(idX,3), &
                                         PolyData%val_diri_el(idX,4)

           case('NEUM')
            inX = inX + 1
            read(inline(ileft:iright),*) PolyData%face_tag_neum_el(inX), &
                                         PolyData%space_fun_tag_neum_el(inX), &
                                         PolyData%val_neum_el(inX,1), &
                                         PolyData%val_neum_el(inX,2), &
                                         PolyData%val_neum_el(inX,3), &
                                         PolyData%val_neum_el(inX,4)

           case('NEUN')                                         
            inN = inN + 1                                                         
            read(inline(ileft:iright),*) PolyData%tag_neuN_el(inN), &
                                         PolyData%fun_neuN_el(inN), &                 
                                         PolyData%val_neuN_el(inN,1), &
                                         PolyData%val_neuN_el(inN,2), &
                                         PolyData%val_neuN_el(inN,3), &
                                         PolyData%val_neuN_el(inN,4)                 

           case('PLOX')
            ipX = ipX + 1
            read(inline(ileft:iright),*) PolyData%fun_poiX_el(ipX), &
                                         PolyData%val_poiX_el(ipX,1), &
                                         PolyData%val_poiX_el(ipX,2), &
                                         PolyData%val_poiX_el(ipX,3), &
                                         PolyData%val_poiX_el(ipX,4)

           case('PLOY')
            ipY = ipY + 1
            read(inline(ileft:iright),*) PolyData%fun_poiY_el(ipY), &
                                         PolyData%val_poiY_el(ipY,1), &
                                         PolyData%val_poiY_el(ipY,2), &
                                         PolyData%val_poiY_el(ipY,3), &
                                         PolyData%val_poiY_el(ipY,4)

           case('PLOZ')
            ipZ = ipZ + 1
            read(inline(ileft:iright),*) PolyData%fun_poiZ_el(ipZ),&
                                         PolyData%val_poiZ_el(ipZ,1), &
                                         PolyData%val_poiZ_el(ipZ,2), &
                                         PolyData%val_poiZ_el(ipZ,3), &
                                         PolyData%val_poiZ_el(ipZ,4)
 
           case('PLAX')                                        
            iplX = iplX + 1                                                        
            read(inline(ileft:iright),*) PolyData%fun_plaX_el(iplX), &                        
                                         PolyData%tag_plaX_el(iplX), &
                                         PolyData%val_plaX_el(iplX,1)  
 
           case('PLAY')                                        
            iplY = iplY + 1                                                        
            read(inline(ileft:iright),*) PolyData%fun_plaY_el(iplY), &                        
                                         PolyData%tag_plaY_el(iplY), &
                                         PolyData%val_plaY_el(iplY,1)     
 
           case('PLAZ')
            iplZ = iplZ + 1
            read(inline(ileft:iright),*) PolyData%fun_plaZ_el(iplZ), &
                                         PolyData%tag_plaZ_el(iplZ), &
                                         PolyData%val_plaZ_el(iplZ,1)

           case('FORX')
            ifX = ifX + 1
            read(inline(ileft:iright),*) PolyData%fun_forX_el(ifX), &
                                         PolyData%val_forX_el(ifX,1), &
                                         PolyData%val_forX_el(ifX,2), &
                                         PolyData%val_forX_el(ifX,3), &
                                         PolyData%val_forX_el(ifX,4)

           case('FORY')
            ifY = ifY + 1
            read(inline(ileft:iright),*) PolyData%fun_forY_el(ifY), &
                                         PolyData%val_forY_el(ifY,1), &
                                         PolyData%val_forY_el(ifY,2), &
                                         PolyData%val_forY_el(ifY,3), &
                                         PolyData%val_forY_el(ifY,4)

           case('FORZ')
            ifZ = ifZ + 1
            read(inline(ileft:iright),*) PolyData%fun_forZ_el(ifZ), &
                                         PolyData%val_forZ_el(ifZ,1), &
                                         PolyData%val_forZ_el(ifZ,2), &
                                         PolyData%val_forZ_el(ifZ,3), &
                                         PolyData%val_forZ_el(ifZ,4)

           case('ABSO')
            iabc = iabc + 1
            read(inline(ileft:iright),*) PolyData%tag_abc_el(iabc)
                                   
           case('SISM')        
              isism = isism + 1                                                                
              if (PolyData%srcmodflag == 0) then
                  read(inline(ileft:iright),*) PolyData%fun_sism_el(isism), & 
                                               PolyData%tag_sism_el(isism), &
                                               PolyData%val_sism_el(isism,1), &
                                               PolyData%val_sism_el(isism,2), & 
                                               PolyData%val_sism_el(isism,3), &
                                               PolyData%val_sism_el(isism,4), &
                                               PolyData%val_sism_el(isism,5), &            
                                               PolyData%val_sism_el(isism,6), &  
                                               PolyData%val_sism_el(isism,7), &
                                               PolyData%val_sism_el(isism,8), &                
                                               PolyData%val_sism_el(isism,9), &
                                               PolyData%val_sism_el(isism,10), &
                                               PolyData%val_sism_el(isism,11), & 
                                               PolyData%val_sism_el(isism,12), & 
                                               PolyData%val_sism_el(isism,13), & 
                                               PolyData%val_sism_el(isism,14), &         
                                               PolyData%val_sism_el(isism,15), & 
                                               PolyData%val_sism_el(isism,16), & 
                                               PolyData%val_sism_el(isism,17), &        
                                               PolyData%val_sism_el(isism,18), &
                                               PolyData%val_sism_el(isism,19), &
                                               PolyData%val_sism_el(isism,20), &        
                                               PolyData%val_sism_el(isism,21)
              elseif (PolyData%srcmodflag == 1) then
                  read(inline(ileft:iright),*) PolyData%fun_sism_el(isism), &
                                               PolyData%tag_sism_el(isism), &
                                               PolyData%val_sism_el(isism,1), &
                                               PolyData%val_sism_el(isism,2), &         
                                               PolyData%val_sism_el(isism,3), &
                                               PolyData%val_sism_el(isism,4), &
                                               PolyData%val_sism_el(isism,5), &                
                                               PolyData%val_sism_el(isism,6), &
                                               PolyData%val_sism_el(isism,7), &
                                               PolyData%val_sism_el(isism,8), &                
                                               PolyData%val_sism_el(isism,9), &
                                               PolyData%val_sism_el(isism,10), &
                                               PolyData%val_sism_el(isism,11), &  
                                               PolyData%val_sism_el(isism,12), &
                                               PolyData%val_sism_el(isism,13), &
                                               PolyData%val_sism_el(isism,14), &        
                                               PolyData%val_sism_el(isism,15)
              endif

           case('CASE')                                                        
              icase = icase + 1                                                                
              read(inline(ileft:iright),*) PolyData%tag_case(icase), &
                                           PolyData%val_case(icase), &
                                           PolyData%tol_case(icase)
                                                                           
           case('NHEE')
              inhee = inhee + 1
              PolyData%tol_nhe(inhee) = 0.0
              read(inline(ileft:iright),*) PolyData%val_nhe(inhee) 

           case('FMAX') 
              read(inline(ileft:iright),*) PolyData%fmax
              
           case('FPEK') 
              read(inline(ileft:iright),*) PolyData%fpeak        
                         
           case('FUNC')
            ifunc = ifunc + 1
            read(inline(ileft:iright),*) PolyData%tag_func(ifunc), &
                                         PolyData%func_type(ifunc)

            select case (PolyData%func_type(ifunc))
            
               case(0,32)
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) + 0 
               
               case(8,9,100)
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) + 1
                 read(inline(ileft:iright),*) dummy, dummy, &
                    (PolyData%func_data(j), j = PolyData%func_indx(ifunc), &
                                                PolyData%func_indx(ifunc +1) -1)
               case(1,2,4,5,6,7,13,14,50,55,99)
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) + 2
                 read(inline(ileft:iright),*) dummy,dummy,&
                     (PolyData%func_data(j), j = PolyData%func_indx(ifunc), &
                                                 PolyData%func_indx(ifunc +1) -1)
                                                 
               !SIGMOIDAL FUNCTION
               case(12)
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) + 3
                 read(inline(ileft:iright),*) dummy, dummy, &
                    (PolyData%func_data(j), j = PolyData%func_indx(ifunc), &
                                                PolyData%func_indx(ifunc +1) -1)

               !ERF FUNCTION
               case(15)
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) + 4
                 read(inline(ileft:iright),*) dummy, dummy, &
                    (PolyData%func_data(j), j = PolyData%func_indx(ifunc), &
                                                PolyData%func_indx(ifunc +1) -1)
               
                    
               case(3,30,31,33)
                 read(inline(ileft:iright),*) dummy, dummy, &
                                              PolyData%nfunc_data, fileinput
                                              
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc)  &
                                                 + 2*PolyData%nfunc_data

                 open(24,file = fileinput)
                 read(24,*) file_nd
                 if (PolyData%nfunc_data .ne. file_nd) then
                    write(*,*) 'Error reading function ! nfunc_data .ne. file_nd !'
                    write(*,*) 'Error reading function from file ', trim(fileinput), '!'
                    write(*,*) 'Line numbers not consistent with material file.'
                    call EXIT(EXIT_FUNCTION_ERROR)
                 endif
                 do j = 1, file_nd
                    i = PolyData%func_indx(ifunc) + 2*(j -1)
                    read(24,*) PolyData%func_data(i), PolyData%func_data(i +1)
                 enddo
                 close(24)
            
                               

               case(60,61,62,63)                                   
                 read(inline(ileft:iright),*) dummy,dummy, PolyData%nfunc_data                 
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) &
                                                + 2*PolyData%nfunc_data                  
                 read(inline(ileft:iright),*) dummy, dummy,dummy, &                 
                    (PolyData%func_data(j), j = PolyData%func_indx(ifunc), &
                                                PolyData%func_indx(ifunc +1) -1)    

                  
              case(773)
                 read(inline(ileft:iright),*) dummy, dummy, PolyData%nfunc_data, &
                                              fileinput
                 PolyData%func_indx(ifunc +1) = PolyData%func_indx(ifunc) &
                                                + PolyData%nfunc_data

                 open(24,file=fileinput)
                 read(24,*) file_nd
                 if (PolyData%nfunc_data .ne. file_nd) then
                    write(*,*) 'Error reading function ! nfunc_data .ne. file_nd !'
                    write(*,*) 'Error reading function from file ', trim(fileinput), '!'
                    write(*,*) 'Line numbers not consistent with material file.'
                    call EXIT(EXIT_FUNCTION_ERROR)
                 endif
                 do j = 1, file_nd
                    i = PolyData%func_indx(ifunc) + (j -1)
                    read(24,*) PolyData%func_data(i)
                 enddo
                 close(24)

            
            end select
                                           
         end select

     enddo
     

     close(40)


     if(PolyData%fmax == 0.d0) then
         PolyData%fmax = 3.d0; 
         if (mpi_id == 0) write(*,'(A)') 'ATTENTION: FMAX not defined!'
         if (mpi_id == 0) write(*,'(A)') 'FMAX assumed = 3'
     endif   

     ! modifying zeta value for frequency proportional damping
     if (damping_type == 1) then 
         do im = 1, PolyData%nmat
             if(PolyData%QS(im) == 0.d0) then 
                PolyData%prop_mat(im,4) = 0.d0;
             else
                PolyData%prop_mat(im,4) = pi*(PolyData%fmax)/PolyData%QS(im)
             endif   
         enddo
      else 
         PolyData%prop_mat(:,4) = 0.d0;
      endif   
      
     end subroutine READ_MATEFILE


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    READ_DIME_FILEMESH(mesh_file)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


     subroutine READ_DIME_MESHFILE(filemesh,PolyData,PolyMesh)
      
     use Poly_data     
     use Poly_mesh
      
     implicit none
      
     character(len=50)     :: filemesh                   
     character(len=100)    :: inline
     character(len=20)     :: keyword
     
     integer(kind=4)       :: i,ie, str_len, ileft, iright, trash, &
                              control, mat_code, elem_total
      
      
     type(Data_Structure), intent(in)  :: PolyData
     type(Mesh_Structure), intent(out) :: PolyMesh
     
     PolyMesh%num_hex   = 0;
     PolyMesh%num_tet   = 0;
     PolyMesh%num_prysm = 0;
     PolyMesh%num_quad = 0;
     PolyMesh%num_tria = 0;
     PolyMesh%num_node = 0;
     PolyMesh%num_poly = 0;
     
                  
     open(40,file=filemesh)
      
     do 
       read(40,'(A)') inline
       if (inline(1:1) .ne. '#') exit
     enddo
      
     read(inline,*) PolyMesh%num_node, elem_total
      
     do i = 1, PolyMesh%num_node
        read(40,'(A)') inline
     enddo
      
     do ie = 1, elem_total
       read(40,'(A)')inline
         
       str_len = len(inline)
       ileft = 0
       iright = 0 
       do i = 1, str_len
          if (inline(i:i).ge.'A') exit
       enddo
       ileft = i
       do i = ileft, str_len
          if (inline(i:i).lt.'A') exit
       enddo
       iright = i
         
       keyword = inline(ileft:iright)
       !write(*,*) keyword
       !read(*,*)
         
       read(inline(1:ileft),*) trash, mat_code
         
       if ((keyword == 'hex') .or. (keyword =='HEX')) then
         control = 0
         do i = 1, PolyData%nmat
           if (PolyData%tag_mat(i) == mat_code) control = 1
           
         enddo
            
         if (control /= 0) PolyMesh%num_hex = PolyMesh%num_hex +1

       elseif ((keyword == 'tetra') .or. (keyword =='TETRA')) then
         control = 0
         do i = 1, PolyData%nmat
           if (PolyData%tag_mat(i) == mat_code) control = 1
         enddo
            
         if (control /= 0) PolyMesh%num_tet = PolyMesh%num_tet +1

       elseif ((keyword == 'pyram') .or. (keyword =='PYRAM')) then
         control = 0
         do i = 1, PolyData%nmat
           if (PolyData%tag_mat(i) == mat_code) control = 1
         enddo
            
         if (control /= 0) PolyMesh%num_prysm = PolyMesh%num_prysm +1

         
       elseif ((keyword == 'quad') .or. (keyword == 'QUAD')) then
          control = 0
          do i = 1, PolyData%nload_diri_el
             if (PolyData%face_tag_diri_el(i) == mat_code) control = 1
          enddo
          do i = 1, PolyData%nload_neum_el
             if (PolyData%face_tag_neum_el(i) == mat_code) control = 1
          enddo
          do i = 1, PolyData%nload_neuN_el                     
             if (PolyData%tag_neuN_el(i) == mat_code) control = 1
          enddo                                                                
          do i = 1, PolyData%nload_abc_el
             if (PolyData%tag_abc_el(i) == mat_code) control = 1
          enddo

          if (control /= 0) PolyMesh%num_quad = PolyMesh%num_quad +1

       elseif ((keyword == 'tria') .or. (keyword == 'TRIA')) then
          control = 0
          do i = 1, PolyData%nload_diri_el
             if (PolyData%face_tag_diri_el(i) == mat_code) control = 1
          enddo
          do i = 1, PolyData%nload_neum_el
             if (PolyData%face_tag_neum_el(i) == mat_code) control = 1
          enddo
          do i = 1, PolyData%nload_neuN_el                     
             if (PolyData%tag_neuN_el(i) == mat_code) control = 1
          enddo                                                                
          do i = 1, PolyData%nload_abc_el
             if (PolyData%tag_abc_el(i) == mat_code) control = 1
          enddo

          if (control /= 0) PolyMesh%num_tria = PolyMesh%num_tria +1

       endif
     enddo
      
      
      
     close(40)
     end subroutine READ_DIME_MESHFILE


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    READ_MESHFILE(mesh_file)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


     subroutine READ_MESHFILE(filemesh,PolyData,PolyMesh)
      
     use Poly_data     
     use Poly_mesh
      
     implicit none
      
     character(len=50)     :: filemesh                   
     character(len=100)    :: inline
     character(len=20)     :: keyword
     
     integer(kind=4)       :: id_status, id_node, i, j, ie, str_len, ileft, iright, & 
                              trash, mat_code, ihexa, iquad, inode, control, &
                              elem_total, itetra, iprysm, itria,ipoly,check_poly, &
                              space_fun_tag
                              
     real(kind=8)          :: xx, yy, zz
      
      
     type(Data_Structure), intent(in)  :: PolyData
     type(Mesh_Structure), intent(inout) :: PolyMesh
    
    
     id_status = 0; ihexa = 0; iquad = 0; 
     itetra = 0; iprysm = 0; itria = 0; ipoly=0;
     check_poly=0;
      
      
      open(40,file=filemesh)
      
      do 
        read(40,'(A)') inline
        if (inline(1:1) /= '#') exit
      enddo
      
      read(inline,*) PolyMesh%num_node, elem_total ! elem_total for both 3D, 2D
      
      ! read nodes
      do i = 1, PolyMesh%num_node
      
        !read(40,*) id_node, PolyMesh%vert_x(i), PolyMesh%vert_y(i), PolyMesh%vert_z(i)
        read(40,*) id_node, xx, yy, zz
        if (inode /= i)  id_status = 1
        
      enddo
      
      ! read 3D elements and 2D faces
      do ie = 1, elem_total
         
         ! look for a word (alphabetic characters)
         read(40,'(A)')inline
         str_len = len(inline)
         ileft = 0
         iright = 0 
         do i = 1,str_len
            if (inline(i:i) >= 'A') exit
         enddo
         ileft = i
         do i = ileft,str_len
            if (inline(i:i) <= 'A') exit
         enddo
         iright = i
         
         keyword = inline(ileft:iright)
         
         read(inline(1:ileft),*) trash, mat_code
         
         ! 3D, 2nd column for material tag (MATE) - mate
         if ((keyword == 'hex') .or. (keyword =='HEX')) then
            control = 0
            do i = 1, PolyData%nmat
               if (PolyData%tag_mat(i) == mat_code) control = 1
            enddo
            
            if (control /= 0) then
               ihexa = ihexa + 1
               PolyMesh%con_hex(ihexa,1) = mat_code
               read(inline(iright:str_len),*)(PolyMesh%con_hex(ihexa,j),j=2,9)
            endif
        !
         elseif ((keyword == 'tetra') .or. (keyword =='TETRA')) then
            control = 0
            do i = 1, PolyData%nmat
               if (PolyData%tag_mat(i) == mat_code) control = 1
            enddo
            
            if (control /= 0) then
               itetra = itetra + 1
               !write(*,*) PolyMesh%con_tet(itetra,:)
               !read(*,*)               
               PolyMesh%con_tet(itetra,1) = mat_code
               !write(*,*) PolyMesh%con_tet(itetra,:)
               !read(*,*)
               read(inline(iright:str_len),*)(PolyMesh%con_tet(itetra,j),j=2,5)
            endif

         elseif ((keyword == 'pyram') .or. (keyword =='PYRAM')) then
            control = 0
            do i = 1, PolyData%nmat
               if (PolyData%tag_mat(i) == mat_code) control = 1
            enddo
            
            if (control /= 0) then
               iprysm = iprysm + 1
               PolyMesh%con_prysm(iprysm,1) = mat_code
               read(inline(iright:str_len),*)(PolyMesh%con_prysm(iprysm,j),j=2,6)
            endif

         ! 2D, 2nd column is the face tag - mate
         ! TODO - assign space function tag for quads and poly?
         elseif ((keyword == 'quad') .or. (keyword == 'QUAD')) then
           control = 0
           do i = 1, PolyData%nload_diri_el
              if (PolyData%face_tag_diri_el(i) == mat_code) control = 1
           enddo
           do i = 1, PolyData%nload_neum_el
              if (PolyData%face_tag_neum_el(i) == mat_code) control = 1
           enddo
           do i = 1, PolyData%nload_neuN_el                               
              if (PolyData%tag_neuN_el(i) == mat_code) control = 1          
           enddo                                                                
           do i = 1, PolyData%nload_abc_el
              if (PolyData%tag_abc_el(i) == mat_code) control = 1
           enddo

           if (control /= 0) then
             iquad = iquad + 1
             !write(*,*) iquad, PolyMesh%num_quad
             !read(*,*)
             PolyMesh%con_quad(iquad,1) = mat_code
             read(inline(iright:str_len),*)(PolyMesh%con_quad(iquad,j),j=2,5)
             !write(*,*) iquad, PolyMesh%con_quad(iquad,:)
             !read(*,*)
           endif
         
         elseif ((keyword == 'tria') .or. (keyword == 'TRIA')) then
           control = 0
           do i = 1, PolyData%nload_diri_el
              if (PolyData%face_tag_diri_el(i) == mat_code) then
               control = 1
               space_fun_tag = PolyData%space_fun_tag_diri_el(i) ! find space function tag dirichlet
              endif
           enddo
           do i = 1, PolyData%nload_neum_el
              if (PolyData%face_tag_neum_el(i) == mat_code) then
               control = 1
               space_fun_tag = PolyData%space_fun_tag_neum_el(i) ! find space function tag neumann
              endif
           enddo
           do i = 1, PolyData%nload_neuN_el                               
              if (PolyData%tag_neuN_el(i) == mat_code) control = 1          
           enddo                                                                
           do i = 1, PolyData%nload_abc_el
              if (PolyData%tag_abc_el(i) == mat_code) control = 1
           enddo

           if (control /= 0) then
             itria = itria + 1
             !write(*,*) iquad, PolyMesh%num_quad
             !read(*,*)
             PolyMesh%con_tria(itria,1) = mat_code
             read(inline(iright:str_len),*)(PolyMesh%con_tria(itria,j),j=2,4)
             PolyMesh%con_tria(itria,5) = space_fun_tag ! assign space function tag
             !write(*,*) iquad, PolyMesh%con_quad(iquad,:)
             !read(*,*)
           endif

         elseif ((keyword == 'poly') .or. (keyword =='POLY')) then
            check_poly = 1

            !do i = 1, PolyData%nmat
            !   if (PolyData%tag_mat(i) == mat_code) control = 1
            !enddo
            
            !if (control /= 0) then
            if (check_poly /= 0) then
               !check_poly=1;
               ipoly = ipoly + 1
               PolyMesh%elem_in_poly(ipoly) = mat_code
               !write(*,*) PolyMesh%con_tet(itetra,:)
               !read(*,*)
               read(inline(iright:str_len),*) PolyMesh%elem_in_poly(ipoly)             
            endif

         endif
      enddo

      ! if not a polygon, set to 0 number of elements in polygon
      if (check_poly == 0) then
         PolyMesh%num_poly=PolyMesh%num_elem
         do i=1,PolyMesh%num_poly
            !PolyMesh%elem_in_poly(i)=i
            PolyMesh%elem_in_poly(i)=0
         end do
      else
         PolyMesh%num_poly=maxval(PolyMesh%elem_in_poly)
      endif
      
      close(40)
      
      end subroutine READ_MESHFILE
























     
