PROGRAM main

   use global
   use controls
   use CD, only: CD_run
   use CD_rot, only: CD_rot_run
   use BB, only: BB_run
   use SW_conv, only: SW_conv1G_run, SW_conv3G_run, SW_conv5G_run, SW_spex_run
   implicit none
   
   interface
      subroutine calculate(shock) 
         character(LEN=*), intent(in) :: shock
      end subroutine calculate
   end interface
   

   write(*,*) '\n ============================'
   write(*,*) ' START'
   write(*,*) '============================\n '

   

   if (BBrad_calc) call BB_run

   !---------------------------------------------------------------!
   !                                                               !
   !                  Shocks                                       !
   !                                                               !
   !---------------------------------------------------------------!

   if (CD_calc) call CD_run

   


   !--------------------------------------------------------------

   write(*,*) '\n ----------------------------'
   write(*,*) ' Reverse Shock'
   write(*,*) '----------------------------'

   if (RS) call calculate('RS')

   write(*,*) '\n ----------------------------'
   write(*,*) ' Forward shock'
   write(*,*) '----------------------------'

   if (FS) call calculate('FS')

   write(*,*) '\n ======================================'
   write(*,*) '         END'
   write(*,*) '======================================'


   
 END PROGRAM main



 Subroutine calculate(shock)
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
   use opac_ff, only: opac_ff_run
   use CD_abssyn, only: CD_abssyn_run
   use CD_rot, only: CD_rot_run
   
   implicit none
   80 format(A9,X,Es10.2)

   character(LEN=*), intent(in) :: shock 
   real(dp) :: Lsy_one,Lbr_one,Lpp_one, LIC_one,LIC_all,LIC_int,L_X_int
   real(dp) :: E_IC_st(mE_f),E_IC_IR(mE_f),LIC_st_tot(mE_f),LIC_IR_tot(mE_f)
   real(dp) :: E_sy(m_nu),Lsy_sh_tot(m_nu) 
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
 
   if (termo_calc) call CD_termo_run(shock)

   if (Xrad_calc) then
      call shock_thermal_run_total
      call shock_thermal_run_spectrum
   endif


   !===============================================================!
   !                                                               ! 
   !                      Electrons                                !
   !                                                               !
   !===============================================================!

   IF (electrons) then
      write(*,*) '\n ELECTRONS:'

       ! -----------------------------------------------------------------------------
       !------------------------------------------------------------------------------
      !if (maps_calc) then      
      !   open(88,file='SW_radiomap_1GHz.dat')
      !   open(89,file='SW_radiomap_3GHz.dat')
      !   open(90,file='SW_radiomap_5GHz.dat')
      !endif
       ! -----------------------------------------------------------------------------
       !------------------------------------------------------------------------------


      if (diste_calc) call CD_diste_run(shock) 

      if (radsyn_sw_calc) call CD_radsyn_run(shock)       ! Just one linear emitter      

      if (opacff_calc) call opac_ff_run

      if (radbr_sw_calc) call CD_radbr_run(unscr,shock)  ! Just one linear emitter  

      if (radic_ir_calc) call CD_radic_ir_run(shock)      !(phi)

      ! IC-star emission, g-g and f-f absorption depend on the 3D position, so we make
      ! linear-emitters for different phi-values.      
      phi = 0.d0
      LIC_all = 0.d0
      if (rotate) then
         DO i_phi = 1, mphi
            write(*,*) 'i(phi)=',i_phi,'/',mphi
            !     write(*,*) 'phi(deg)=', phi*180_dp/pi

            call CD_rot_run(phi)   

            if (radic_st_calc) then
               write(*,*) 'Calculate the IC(star) emission'
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

            if (abssyn_calc) then
               write(*,*) 'Calculate the free-free synchrotron absorption'
               call CD_abssyn_run(shock) !(phi)
               open (45,file='radsyn_abs_'//shock//'.dat')
               do i_E = 1,m_nu    ! SED --> sum emission from different phi
                  read(45,*) E_sy(i_E),Lsy_unabs,Lsy_abs
                  Lsy_sh_unabs(i_E) = Lsy_sh_unabs(i_E) + Lsy_unabs
                  Lsy_sh_abs(i_E) = Lsy_sh_abs(i_E) + Lsy_abs
               end do
               close(45)
            endif
            ! -----------------------------------------------------------------------------
            !------------------------------------------------------------------------------ 
           ! if (maps_calc) then           
           !    write(*,*) 'Save synthetic radio emission maps'
           !    open(55,file='CD_radiomap_1GHz.dat') ! 1.4 GHz
           !    open(56,file='CD_radiomap_3GHz.dat') ! 3.0 GHz
           !    open(57,file='CD_radiomap_5GHz.dat') ! 4.8 GHz
           !    do i_l = 1,2*m_theta
           !       read(55,*) u_proj,v_proj,emi,emi_abs
           !       write(88,*) u_proj,v_proj,emi,emi_abs
           !       if (i_phi.ne.1) then ! .and. i_phi.ne.mphi) then
           !          write(88,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
           !       endif
           !       read(56,*) u_proj,v_proj,emi,emi_abs
           !       write(89,*) u_proj,v_proj,emi,emi_abs  
           !       if (i_phi.ne.1) then !.and. i_phi.ne.mphi) then
           !          write(89,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
           !       endif
           !       read(57,*) u_proj,v_proj,emi,emi_abs
           !       write(90,*) u_proj,v_proj,emi,emi_abs  
           !       if (i_phi.ne.1) then !.and. i_phi.ne.mphi) then
           !          write(90,*) u_proj,-v_proj,emi,emi_abs ! symmetry???
           !       endif
           !    end do
           !    close(55)
           !    close(56)
           !    close(57)
           ! endif
            ! -----------------------------------------------------------------------------
            !------------------------------------------------------------------------------

            phi = phi + delta_phi           
         end do

         write(*,*) 'The total L_IC_st is:'
         write(*,80) 'L_IC_st = ', LIC_all ! 

         write(*,*) ''

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
               read(41,*) E_sy(i_E),Lsy_one
               Lsy_sh_tot(i_E) = Lsy_one*mphi * 4.d0
               ! Multiply by a factor 4 because of the axial symmetry in the SW 
               ! and "mirroring" for y < 0
               write(51,*) E_sy(i_E),Lsy_sh_tot(i_E)
            end do
            close(41)
            close(51)
         endif

         if (abssyn_calc) then       ! This has the full 3D geometry in the abs
            open (51,file='radsynabs_'//shock//'_tot.dat')
            do i_E = 1,m_nu
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! (this one already considers "mirroring" for y < 0 in the calculation)
               write(51,*) E_sy(i_E),2.d0*Lsy_sh_unabs(i_E),2.d0*Lsy_sh_abs(i_E)
            end do
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
            Eph_min = 1.d3 ! 1 keV
            Eph_max = 1.d4 ! 10 keV
            
            open (31,file='radic_st_'//shock//'_tot.dat')
            do i_E = 1, mE_f
               read(31,*) E_IC_st(i_E), LIC_st_tot(i_E)
               dE_IC_st(i_E) = (E_IC_st(2)/E_IC_st(1) - 1.d0)
               if ( E_IC_st(i_E) >= Eph_min .and. E_IC_st(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + LIC_st_tot(i_E)*dE_IC_st(i_E)
               end if
            end do
            close(31)

            open (33,file='radic_ir_'//shock//'_tot.dat')
            do i_E = 1, mE_f
               read(33,*) E_IC_IR(i_E),LIC_IR_tot(i_E)
               dE_IC_IR(i_E) = (E_IC_IR(2)/E_IC_IR(1) - 1.d0)
               if ( E_IC_IR(i_E) >= Eph_min .and. E_IC_IR(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + LIC_IR_tot(i_E)*dE_IC_IR(i_E)
                end if
            end do
            close(33)

            open (51,file='radsyn_'//shock//'_tot.dat')
            do i_E = 1, m_nu
               read(51,*) E_sy(i_E), Lsy_sh_tot(i_E)
               dE_sy(i_E) = (E_sy(2)/E_sy(1) - 1.d0)
               if ( E_sy(i_E) >= Eph_min .and. E_sy(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + Lsy_sh_tot(i_E)*dE_sy(i_E)
               end if
            end do
            close(51)

            write(*,*) 'Calculate the luminosity btw 1-10 keV'
            write(*,80) 'L_X= ', L_X_int
            write(*,*) ''
         end if

      endif


   ENDIF

   ! -----------------------------------------------------------------------------
   !------------------------------------------------------------------------------  
 !  if (convolution) then
 !     write(*,*) 'Convolution of the emission map with the instrument response'
 !     call SW_conv1G_run
 !     call SW_conv3G_run  
 !     call SW_conv5G_run  
 !     call SW_spex_run  

 !  endif
   ! -----------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   !--------------------------------

   IF (protons) then
      write(*,*) '\n PROTONS:'

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



