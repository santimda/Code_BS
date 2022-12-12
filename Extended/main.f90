PROGRAM main
   use global
   use controls
   use initialize
   use CD, only: CD_run
   use BB, only: BB_run
   use opac_ff, only: opac_ff_run
   implicit none
   logical, parameter :: verbose=.True.
   
   interface
      subroutine calculate(shock, radiative) 
         character(LEN=*), intent(in) :: shock
         logical, intent(in) :: radiative
      end subroutine calculate
   end interface
   
   interface
      subroutine abs_maps(shock)
         character(LEN=*), intent(in) :: shock
      end subroutine abs_maps
   end interface

   print *, '\n ============================'
   print *, ' START'
   print *, '============================\n '

   ! Initialize vectors   
   if (verbose) print *, '\n Allocate vectors'    
   call initialize_allocate()

   if (BBrad_calc) call BB_run

   if (CD_calc) call CD_run

   if (opacff_calc) call opac_ff_run

   if (RS) then 
      print *, '\n----------------------------'
      print *, ' Reverse Shock'
      print *, '----------------------------'
      call calculate('RS', radiative=.false.)
      if (maps_calc) then
         if (.not. BOTH) then
            print *, '\n--------------------------------------'
            print *, ' Calculate absorption and emission maps'
            print *, '----------------------------------------'
            call abs_maps('RS')
         end if
      end if
   endif

   if (FS) then 
      print *, '\n----------------------------'
      print *, ' Forward shock'
      print *, '----------------------------'
      call calculate('FS', radiative=.true.)
      if (maps_calc) then
         if (.not. BOTH) then
            print *, '\n--------------------------------------'
            print *, ' Calculate absorption and emission maps'
            print *, '----------------------------------------'
            call abs_maps('FS')
         end if
      end if  
   endif

   if (maps_calc) then
      if (BOTH) then
         print *, '\n--------------------------------------'
         print *, ' Calculate absorption and emission maps'
         print *, '----------------------------------------'
         call abs_maps('BOTH')
      end if
   end if
      
   if (verbose) print *, '\n Deallocate vectors'    
   call initialize_deallocate()

   print *, '\n ======================================'
   print *, '         END'
   print *, '======================================'

   
 END PROGRAM main


 !============================================================================

 Subroutine calculate(shock,radiative)
   use global
   use controls
   use BB, only: BB_run
   use CD_termo, only: CD_termo_run
   use shock_thermal, only: shock_thermal_run_total,shock_thermal_run_spectrum
   use CD_diste, only: CD_diste_run
   use CD_distp, only: CD_distp_run
   use CD_radsyn, only: CD_radsyn_run
   use CD_radbr, only: CD_radbr_run
   use CD_radic, only: CD_radic_st_run,CD_radic_IR_run
   use CD_radpp, only: CD_radpp_run
   use CD_rot, only: CD_rot_run
   
   implicit none
   80 format(A9,X,Es10.2)

   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative
   real(dp) :: Lsy_one,Lbr_one,Lpp_one, LIC_one,LIC_all,LIC_int,L_X_int,L_syn_E_one
   real(dp) :: E_IC_st(mE_f),E_IC_IR(mE_f),LIC_st_tot(mE_f),LIC_IR_tot(mE_f)
   real(dp) :: E_sy(m_nu),Lsy_sh_tot(m_nu),L_syn_E_tot(m_nu)
   real(dp) :: Lsy_sh_abs(m_nu), Lsy_sh_unabs(m_nu), Lsy_unabs, Lsy_abs
   real(dp) :: E_br(mE_f),Lbr_sh_tot(mE_f)
   real(dp) :: E_pp(mE_f),Lpp_sh_tot(mE_f)
   real(dp) :: Eph_min, Eph_max
   real(dp) :: dE_sy(m_nu),dE_IC_IR(mE_f),dE_IC_st(mE_f)
   real(dp) :: u_proj,v_proj,emi,emi_abs
   real(dp) :: phi
   integer :: i_phi,i_E,i_l  

   LIC_st_tot = 0.d0
   Lbr_sh_tot = 0.d0
   Lsy_sh_abs = 0.d0
   Lsy_sh_unabs = 0.d0
   Lsy_sh_tot = 0.d0
   Lpp_sh_tot = 0.d0
 
   if (termo_calc) then
      if (shock == 'RS') then
         call CD_termo_run(shock, radiative =.false.)
      else if (shock == 'FS') then
         call CD_termo_run(shock, radiative =.true.)
      end if
   end if

   if (Xrad_calc .and. radiative) then
      call shock_thermal_run_total(shock, radiative = .true.)
      call shock_thermal_run_spectrum(shock, radiative = .true.)
   else if (Xrad_calc .and. (.not. radiative)) then
      call shock_thermal_run_total(shock, radiative = .false.)
      call shock_thermal_run_spectrum(shock, radiative = .false.)
   endif


   !===============================================================!
   !                                                               ! 
   !                      Electrons                                !
   !                                                               !
   !===============================================================!

   IF (electrons .and. (.not. radiative)) then
      print *, '\n ELECTRONS:'

      ! -----------------------------------------------------------------------------

      if (diste_calc) call CD_diste_run(shock) 

      if (radsyn_sw_calc) call CD_radsyn_run(shock)       ! Just one linear emitter      

      if (radbr_sw_calc) call CD_radbr_run(unscr,shock)  ! Just one linear emitter  

      if (radic_ir_calc) call CD_radic_ir_run(shock)      !(phi)

      ! IC-star emission, g-g and f-f absorption depend on the 3D position, so we make
      ! linear-emitters for different phi-values.      
      phi = 0.d0
      LIC_all = 0.d0
      if (rotate) then
         DO i_phi = 1, mphi
            print *, 'i(phi)=',i_phi,'/',mphi
            !     print *, 'phi(deg)=', phi*180_dp/pi

            call CD_rot_run(phi)   

            if (radic_st_calc) then
               print *, 'Calculate the IC(star) emission'
               call CD_radic_st_run(LIC_int,shock) 
               open (35,file='radic_st_'//shock//'.dat')
               do i_E = 1,mE_f    ! SED --> sum emission from different phi
                  read(35,*) E_IC_st(i_E),LIC_one
                  ! LIC_one is for linear emitters with 0 <= phi <= pi/2; I later multiply by 2
                  LIC_st_tot(i_E) = LIC_st_tot(i_E) + LIC_one
               end do
               close(35)
               LIC_all = LIC_all + LIC_int
            endif

            ! -----------------------------------------------------------------------------
            !------------------------------------------------------------------------------

            phi = phi + delta_phi           
         end do

         if (radic_ir_calc) then
            print *, 'The total L_IC_st is:'
            write(*,80) 'L_IC_st = ', LIC_all ! 
            print *, ''
         endif

         close(88) ! maps for conv
         close(89)
         close(90)

      endif

      !=============================================================
      !   
      !  Save for the total linear-emitters 
      !  
      if (escribir_e) then

         if (radic_st_calc) then       ! This has the full 3D geometry in the abs
            open (31,file='radic_st_'//shock//'_tot.dat')
            do i_E = 1,mE_f
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! (this one already considers "mirroring" for y < 0 in the calculation)
               write(31,*) E_IC_st(i_E),2.d0*LIC_st_tot(i_E)
            end do
            close(31)
         endif

         if (radic_ir_calc) then
            open (32,file='radic_ir_'//shock//'.dat')     ! One linear-emitter
            open (33,file='radic_ir_'//shock//'_tot.dat') ! The total emission IC(dust) in the SW
            do i_E = 1,mE_f
               read(32,*) E_IC_IR(i_E),LIC_one
               LIC_IR_tot(i_E) = LIC_one*mphi *4.d0
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! and "mirroring" for y < 0
               write(33,*) E_IC_IR(i_E),LIC_IR_tot(i_E) 
            end do
            close(32)
            close(33)
         endif

         !       if ( radic_ir_calc .eqv. yes ) then       ! This has the full 3D geometry in the abs
         !          open (32,file='radic_ir_tot.dat')
         !          do i_E = 1,mE_f
         !             ! Multiply by a factor 2 because of the axial symmetry in the SW
         !             write(32,*) E_IC_IR(i_E),2.d0*LIC_IR_tot(i_E)
         !          end do
         !          close(32)
         !       endif

         if (radsyn_sw_calc) then   ! NOT VERY USEFUL
            open (51,file='radsyn_'//shock//'_tot.dat')
            open (41,file='radsyn_'//shock//'.dat')
            do i_E = 1,m_nu
               read(41,*) E_sy(i_E),Lsy_one,L_syn_E_one
               Lsy_sh_tot(i_E) = Lsy_one*mphi * 4.d0
               L_syn_E_tot(i_E) = L_syn_E_one*mphi * 4.d0
               ! Multiply by a factor 4 because of the axial symmetry in the SW 
               ! and "mirroring" for y < 0
               write(51,*) E_sy(i_E),Lsy_sh_tot(i_E),L_syn_E_tot(i_E)
            end do
            close(41)
            close(51)
         endif


         if (radbr_sw_calc) then
            open (42,file='radbr_'//shock//'.dat')     ! One linear-emitter
            open (52,file='radbr_'//shock//'_tot.dat') ! The total SW
            do i_E = 1,mE_f
               read(42,*) E_br(i_E),Lbr_one
               Lbr_sh_tot(i_E) = Lbr_one*mphi *4.d0
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! and "mirroring" for y < 0
               write(52,*) E_br(i_E),Lbr_sh_tot(i_E) 
            end do
            close(42)
            close(52)
         endif


         if ( (radic_st_calc) .and. (radic_ir_calc) .and. (radsyn_sw_calc) ) then
            ! Calculates the X-ray luminosity
            L_X_int = 0.d0
            Eph_min = 0.4d3 ! 0.4 keV
            Eph_max = 4.d3 ! 4 keV
            
            open (31,file='radic_st_RS_tot.dat')
            do i_E = 1, mE_f
               read(31,*) E_IC_st(i_E), LIC_st_tot(i_E)
               dE_IC_st(i_E) = (E_IC_st(2)/E_IC_st(1) - 1.d0)
               if ( E_IC_st(i_E) >= Eph_min .and. E_IC_st(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + 2.d0*LIC_st_tot(i_E)*dE_IC_st(i_E) ! Multiply by 2 because of axial symmetry
               end if
            end do
            close(31)

            open (33,file='radic_ir_RS_tot.dat')
            do i_E = 1, mE_f
               read(33,*) E_IC_IR(i_E),LIC_IR_tot(i_E)
               dE_IC_IR(i_E) = (E_IC_IR(2)/E_IC_IR(1) - 1.d0)
               if ( E_IC_IR(i_E) >= Eph_min .and. E_IC_IR(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + LIC_IR_tot(i_E)*dE_IC_IR(i_E)
                end if
            end do
            close(33)

            open (51,file='radsyn_RS_tot.dat')
            do i_E = 1, m_nu
               read(51,*) E_sy(i_E), Lsy_sh_tot(i_E)
               dE_sy(i_E) = (E_sy(2)/E_sy(1) - 1.d0)
               if ( E_sy(i_E) >= Eph_min .and. E_sy(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + Lsy_sh_tot(i_E)*dE_sy(i_E)
               end if
            end do
            close(51)

            print *, 'Calculate the luminosity btw 0.4-4 keV'
            write(*,80) 'L_X= ', L_X_int
            print *, ''
         end if

      endif


   ENDIF

   ! -----------------------------------------------------------------------------
   !------------------------------------------------------------------------------  

   ! -----------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   !--------------------------------

   IF (protons .and. (.not. radiative)) then
      print *, '\n PROTONS:'

      if (distp_calc) call CD_distp_run(shock) 

      if (radpp_calc) call CD_radpp_run(shock) 

      !=============================================================
      !   
      !  Save for the total linear-emitters 
      !  
      if (escribir_p) then

         if (radpp_calc) then
            open (43,file='radpp_'//shock//'.dat')     ! One linear-emitter
            open (53,file='radpp_'//shock//'_tot.dat')
            do i_E = 1,mE_f
               read(43,*) E_pp(i_E),Lpp_one       
               Lpp_sh_tot(i_E) = Lpp_one*mphi *2.d0
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               write(53,*) E_pp(i_E),Lpp_sh_tot(i_E) 
            end do
            close(43)
            close(53)
         endif

      endif

   ENDIF


   !================================================================

End Subroutine Calculate



!============================================================================
!============================================================================

Subroutine abs_maps(SHOCK)
   use global
   use controls
   use CD_abssyn, only: CD_abssyn_run
   use CD_rot, only: CD_rot_run
   use SW_conv, only: SW_conv1G_run, SW_conv3G_run, SW_conv5G_run, SW_convavg_run,SW_spex_run

   implicit none

   character(LEN=*), intent(in) :: shock
   real(dp) :: Lsy_sh_abs(m_nu), Lsy_sh_unabs(m_nu), Lsy_unabs, Lsy_abs
   real(dp) :: S_th_RS, S_th_FS, S_th, S_hot
   real(dp) :: S_th_avg_RS, S_th_avg_FS, S_th_avg, S_hot_avg
   real(dp) :: emi_avg,emi_avg_abs
   real(dp) :: E_sy(m_nu)
   real(dp) :: u_proj,v_proj,emi,emi_abs
   real(dp) :: phi
   integer :: i_phi,i_E,i_l

   Lsy_sh_unabs = 0.d0
   Lsy_sh_abs = 0.d0

   if (shock /= 'BOTH') then    
      open(88,file='sync_radiomap_'//shock//'_1GHz.dat')
      open(89,file='sync_radiomap_'//shock//'_3GHz.dat')
      open(90,file='sync_radiomap_'//shock//'_5GHz.dat')
      open(91,file='sync_radiomap_avg_'//shock//'.dat')
      open(98,file='th_ff_map_'//shock//'_1G.dat')
      open(99,file='th_ff_map_'//shock//'_3G.dat')
      open(100,file='th_ff_map_'//shock//'_5G.dat')
      open(101,file='th_ff_map'//shock//'avg.dat')
   else if (shock == 'BOTH') then
      open(88,file='sync_radiomap_1GHz.dat')
      open(89,file='sync_radiomap_3GHz.dat')
      open(90,file='sync_radiomap_5GHz.dat')
      open(91,file='sync_radiomap_avg.dat')
      open(98,file='th_ff_map_1G.dat')
      open(99,file='th_ff_map_3G.dat')
      open(100,file='th_ff_map_5G.dat')
      open(101,file='th_ff_map_avg.dat')
   end if

   if ((shock == 'BOTH') .or. (shock == 'FS')) then
      open(108,file='hot_layer_1GHz.dat')
      open(109,file='hot_layer_3GHz.dat')
      open(110,file='hot_layer_5GHz.dat')
      open(111,file='hot_layer_avg.dat')
   end if

phi = 0.d0
      if (rotate) then
         DO i_phi = 1, mphi
            print *, 'i(phi)=',i_phi,'/',mphi

            call CD_rot_run(phi)   

            if (abssyn_calc) then
               print *, 'Calculate the free-free synchrotron absorption'
               call CD_abssyn_run(shock) !(phi)
               if (shock /= 'BOTH') then
                  open (45,file='radsyn_abs_'//shock//'.dat')
               else if (shock == 'BOTH') then
                  open (45,file='radsyn_abs.dat')
               end if
               do i_E = 1,m_nu    ! SED --> sum emission from different phi
                  read(45,*) E_sy(i_E),Lsy_unabs,Lsy_abs
                  Lsy_sh_unabs(i_E) = Lsy_sh_unabs(i_E) + Lsy_unabs
                  Lsy_sh_abs(i_E) = Lsy_sh_abs(i_E) + Lsy_abs
               end do
               close(45)
            endif
            ! -----------------------------------------------------------------------------
            !------------------------------------------------------------------------------ 
            if (maps_calc) then           
               print *, 'Save synthetic radio emission maps'
               open(61,file='CD_th_RS_1GHz.dat') ! 1.4 GHz
               open(62,file='CD_th_RS_3GHz.dat') ! 3.0 GHz
               open(63,file='CD_th_RS_5GHz.dat') ! 4.8 GHz
               open(64,file='CD_th_FS_1GHz.dat') ! 1.4 GHz
               open(65,file='CD_th_FS_3GHz.dat') ! 3.0 GHz
               open(66,file='CD_th_FS_5GHz.dat') ! 4.8 GHz
               open(67,file='CD_th_radiomap_avg_RS.dat') ! average ff RS
               open(68,file='CD_th_radiomap_avg_FS.dat') ! averague isothermal region
               if (shock /= 'BOTH') then
                  open(55,file='CD_sync_radiomap_'//shock//'1GHz.dat') ! 1.4 GHz sync
                  open(56,file='CD_sync_radiomap_'//shock//'3GHz.dat') ! 3.0 GHz sync
                  open(57,file='CD_sync_radiomap_'//shock//'5GHz.dat') ! 4.8 GHz sync
                  open(58,file='CD_sync_radiomap_avg'//shock//'.dat')
               else if (shock == 'BOTH') then
                  open(55,file='CD_sync_radiomap_1GHz.dat') ! 1.4 GHz
                  open(56,file='CD_sync_radiomap_3GHz.dat') ! 3.0 GHz
                  open(57,file='CD_sync_radiomap_5GHz.dat') ! 4.8 GHz
                  open(58,file='CD_sync_radiomap_avg.dat')
               end if
               if ((shock == 'BOTH') .or. (shock == 'FS'))  then
                  open(70,file='CD_hot_layer_1GHz.dat')
                  open(71,file='CD_hot_layer_3GHz.dat')
                  open(72,file='CD_hot_layer_5GHz.dat')
                  open(73,file='CD_hot_radiomap_avg.dat') ! average hot layer
               end if

               do i_l = 1,2*m_theta
                  ! 1 GHz
                  read(55,*) u_proj,v_proj,emi,emi_abs
                  read(61,*) S_th_RS
                  read(64,*) S_th_FS
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     read(70,*) S_hot
                  end if
                  S_th = S_th_FS + S_th_RS
                  write(88,*) u_proj,v_proj,emi,emi_abs
                  write(98,*) S_th
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     write(108,*) S_hot
                  end if
                  if (i_phi.ne.1) then ! .and. i_phi.ne.mphi) then
                     write(88,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
                     write(98,*) S_th
                     if ((shock == 'BOTH') .or. (shock == 'FS')) then
                        write(108,*) S_hot
                     end if
                  endif

                  ! 3 GHz
                  read(56,*) u_proj,v_proj,emi,emi_abs
                  read(62,*) S_th_RS
                  read(65,*) S_th_FS
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     read(71,*) S_hot
                  end if
                  S_th = S_th_FS + S_th_RS
                  write(89,*) u_proj,v_proj,emi,emi_abs
                  write(99,*) S_th
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     write(109,*) S_hot
                  end if
                  if (i_phi.ne.1) then !.and. i_phi.ne.mphi) then
                     write(89,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
                     write(99,*) S_th
                     if ((shock == 'BOTH') .or. (shock == 'FS')) then
                        write(109,*) S_hot
                     end if
                  endif

                  ! 5 GHz
                  read(57,*) u_proj,v_proj,emi,emi_abs
                  read(63,*) S_th_RS
                  read(66,*) S_th_FS
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     read(72,*) S_hot
                  end if
                  S_th = S_th_FS + S_th_RS
                  write(90,*) u_proj,v_proj,emi,emi_abs  
                  write(100,*) S_th
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     write(110,*) S_hot
                  end if
                  if (i_phi.ne.1) then !.and. i_phi.ne.mphi) then
                     write(90,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
                     write(100,*) S_th
                     if ((shock == 'BOTH') .or. (shock == 'FS')) then
                        write(110,*) S_hot
                     end if
                  endif

                  ! Average
                  read(58,*) u_proj,v_proj,emi_avg,emi_avg_abs
                  read(67,*) S_th_avg_RS
                  read(68,*) S_th_avg_FS
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     read(73,*) S_hot_avg
                  end if
                  S_th_avg = S_th_avg_FS + S_th_avg_RS
                  write(91,*) u_proj,v_proj,emi_avg,emi_avg_abs
                  write(101,*) S_th_avg
                  if ((shock == 'BOTH') .or. (shock == 'FS')) then
                     write(111,*) S_hot_avg
                  end if
                  if (i_phi.ne.1) then !.and. i_phi.ne.mphi) then
                     write(91,*) u_proj,-v_proj,emi_avg,emi_avg_abs
                     write(101,*) S_th_avg ! symmetry???
                     if ((shock == 'BOTH') .or. (shock == 'FS')) then
                        write(111,*) S_hot_avg
                     end if
                  endif

               end do
               close(55)
               close(56)
               close(57)
               close(58)
               close(61)
               close(62)  
               close(63)
               close(64)
               close(65)
               close(66)
               close(67)
               close(68)
               if ((shock == 'BOTH') .or. (shock == 'FS'))  then
                  close(70)
                  close(71)
                  close(72)
                  close(73)
               end if
            endif
            ! -----------------------------------------------------------------------------
            !------------------------------------------------------------------------------

            phi = phi + delta_phi           
         end do
   end if

   close(88)
   close(89)
   close(90)
   close(91)
   close(98)
   close(99)
   close(100)
   close(101)
   if ((shock == 'BOTH') .or. (shock == 'FS')) then
      close(108)
      close(109)
      close(110)
      close(111)
   end if


   if (abssyn_calc) then       ! This has the full 3D geometry in the abs
            open (51,file='radsynabs_'//shock//'_tot.dat')
            do i_E = 1,m_nu
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! (this one already considers "mirroring" for y < 0 in the calculation)
               write(51,*) E_sy(i_E),2.d0*Lsy_sh_unabs(i_E),2.d0*Lsy_sh_abs(i_E)
            end do
            close(51)
   endif



   if (convolution) then
      print *, 'Convolution of the emission map with the instrument response'
      call SW_conv1G_run(shock)
      call SW_conv3G_run(shock)
      !call SW_conv5G_run(shock)
      call SW_convavg_run(shock)
      !call SW_spex_run(shock)
   endif


   END subroutine abs_maps