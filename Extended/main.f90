PROGRAM main
   use global
   use controls
   use initialize
   use CD, only: CD_run
   use BB, only: BB_run
   use opac_ff, only: opac_ff_run
   use map_freqs, only: maps_init
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
   call maps_init()

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

   if (rad_alpha .and. FS) call H_alpha_map('FS')
      
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
   use radstr, only: radstr_run
   
   implicit none
   80 format(A9,X,Es10.2)

   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative
   real(dp) :: Lsy_one,Lbr_one,Lpp_one, LIC_one,LIC_all,LIC_int,L_X_int,L_syn_E_one
   real(dp) :: E_IC_st(mE_f),E_IC_IR(mE_f),LIC_st_tot(mE_f),LIC_IR_tot(mE_f)
   real(dp) :: E_sy(m_nu),Lsy_sh_tot(m_nu),L_syn_E_tot(m_nu)
   real(dp) :: Lsy_sh_abs(m_nu), Lsy_sh_unabs(m_nu), Lsy_unabs, Lsy_abs
   real(dp) :: Lff_tot(mE_ph), E_ff(mE_ph), dE_ff(mE_ph)
   real(dp) :: E_br(mE_f),Lbr_sh_tot(mE_f)
   real(dp) :: E_pp(mE_f),Lpp_sh_tot(mE_f)
   real(dp) :: Eph_min, Eph_max
   real(dp) :: dE_sy(m_nu),dE_IC_IR(mE_f),dE_IC_st(mE_f)
   real(dp) :: u_proj,v_proj,emi,emi_abs
   real(dp) :: phi
   integer :: i_phi,i_E  

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
               open(54,file='rad_IC_st_tot_'//shock//'.dat')
                  write(54,*) LIC_all
               close(54)
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

      endif

      !=============================================================
      !   
      !  Save for the total linear-emitters 
      !  
      if (escribir_e) then

         if (radic_st_calc) then       ! This has the full 3D geometry in the abs
            open (31,file='radic_st_'//system//'_'//shock//'_tot.dat')
            do i_E = 1,mE_f
               ! Multiply by a factor 2 because of the axial symmetry in the SW
               ! (this one already considers "mirroring" for y < 0 in the calculation)
               write(31,*) E_IC_st(i_E),2.d0*LIC_st_tot(i_E)
            end do
            close(31)
         endif

         if (radic_ir_calc) then
            open (32,file='radic_ir_'//shock//'.dat')     ! One linear-emitter
            open (33,file='radic_ir_'//system//'_'//shock//'_tot.dat') ! The total emission IC(dust) in the SW
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
            open (51,file='radsyn_'//system//'_'//shock//'_tot.dat')
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
            open (52,file='radbr_'//system//'_'//shock//'_tot.dat') ! The total SW
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


         if ( (radic_st_calc) .and. (radic_ir_calc) .and. (radsyn_sw_calc) .and. (Xrad_calc) ) then
            ! Calculates the X-ray luminosity
            L_X_int = 0.d0
            if (system == 'BD+43') then
               Eph_min = 0.4d3 ! 0.4 keV
               Eph_max = 4.d3 ! 4 keV
            else if (system == 'Lambda_Cep') then
               Eph_min = 2.d3 ! 2 keV
               Eph_max = 10.d3 ! 10 keV
            else
               Eph_min = 1.d3 ! 1 keV
               Eph_max = 10.d3 ! 10 keV
            end if               
            
            open (31,file='radic_st_'//system//'_RS_tot.dat')
            do i_E = 1, mE_f
               read(31,*) E_IC_st(i_E), LIC_st_tot(i_E)
               dE_IC_st(i_E) = (E_IC_st(2)/E_IC_st(1) - 1.d0)
               if ( E_IC_st(i_E) >= Eph_min .and. E_IC_st(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + 2.d0*LIC_st_tot(i_E)*dE_IC_st(i_E) ! Multiply by 2 because of axial symmetry
               end if
            end do
            close(31)

            open (33,file='radic_ir_'//system//'_RS_tot.dat')
            do i_E = 1, mE_f
               read(33,*) E_IC_IR(i_E),LIC_IR_tot(i_E)
               dE_IC_IR(i_E) = (E_IC_IR(2)/E_IC_IR(1) - 1.d0)
               if ( E_IC_IR(i_E) >= Eph_min .and. E_IC_IR(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + LIC_IR_tot(i_E)*dE_IC_IR(i_E)
                end if
            end do
            close(33)

            open (51,file='radsyn_'//system//'_RS_tot.dat')
            do i_E = 1, m_nu
               read(51,*) E_sy(i_E), Lsy_sh_tot(i_E)
               dE_sy(i_E) = (E_sy(2)/E_sy(1) - 1.d0)
               if ( E_sy(i_E) >= Eph_min .and. E_sy(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + Lsy_sh_tot(i_E)*dE_sy(i_E)
               end if
            end do
            close(51)

            open (52,file='ionized_ff_'//system//'_RS.dat')
            do i_E = 1, mE_ph
               read(52,*) E_ff(i_E), Lff_tot(i_E)
               dE_ff(i_E) = (E_ff(2)/E_ff(1) - 1.d0)
               if ( E_ff(i_E) >= Eph_min .and. E_ff(i_E) <= Eph_max ) then
                  L_X_int = L_X_int + Lff_tot(i_E)*dE_ff(i_E)
               end if
            end do
            close(52)

            if (system == 'BD+43') then
               print *, 'Calculate the luminosity btw 0.4-4 keV'
               write(*,80) 'L_X= ', L_X_int
               if (L_X_int > 10.**(30.878384)) print*, 'R-X luminosity is to high'  ! Toalá+ 2016
            else if (system == 'Lambda_Cep') then
               print *, 'Calculate the luminosity btw 2-10 keV'
               write(*,80) 'L_X= ', L_X_int
               if (L_X_int > 3.37d30) print*, 'R-X luminosity is to high'
            else
               print *, 'Calculate the luminosity btw 1-10 keV'
                write(*,80) 'L_X= ', L_X_int
            end if
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
            open (53,file='radpp_'//system//'_'//shock//'_tot.dat')
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
   
   if (stromgren .and. electrons .and. protons .and. (.not. radiative)) then
      print*, '\n Calculate the NT emission from the H II region:'
      call radstr_run(shock)
   end if
   
   !================================================================

End Subroutine Calculate

!============================================================================
!============================================================================

Subroutine H_alpha_map(shock)
   use global
   use controls
   use rad_H_alpha, only: rad_H_alpha_run, project_H_alpha_map   
   use CD_rot, only: CD_rot_run
   implicit none

   character(LEN=*), intent(in) :: shock
   real(dp) :: S_H_alpha
   real(dp) :: u_proj, v_proj
   real(dp) :: phi
   integer :: i_phi,i_l

   open(710,file='S_alpha_map_'//shock//'.dat')
   phi = 0.d0
   do i_phi = 1,mphi
      print *, 'i(phi)=',i_phi,'/',mphi
      call CD_rot_run(phi)
      call project_H_alpha_map(shock, phi)       ! This projects the H_alpha map in u,v
      open(70,file='S_alpha_proj_'//shock//'.dat')
      do i_l = 1,2*m_theta
         read(70,*) u_proj,v_proj,S_H_alpha
         write(710,*) u_proj,v_proj,S_H_alpha
         if (i_phi.ne.1) then
            write(710,*) u_proj,-v_proj,S_H_alpha
         endif
      end do
      close(70)
      phi = phi + delta_phi
   end do
   close(710)

   if (convolution_alpha) then
      call rad_H_alpha_run(shock)
   end if

   End Subroutine H_alpha_map


!============================================================================
!============================================================================

Subroutine abs_maps(shock)
   use global
   use controls
   use map_freqs, only: maps
   use CD_abssyn, only: CD_abssyn_run
   use CD_rot, only: CD_rot_run
   use SW_conv, only: SW_conv_run

   implicit none

   character(LEN=*), intent(in) :: shock
   integer, parameter :: m_l = 2*m_theta     ! CD position, both sides
   real(dp) :: Lsy_sh_abs(m_nu), Lsy_sh_unabs(m_nu), Lsy_unabs, Lsy_abs
   real(dp) :: E_sy(m_nu)
   real(dp) :: phi
   integer :: i_phi, i_E, im, n_map
   integer :: iu_syn_avg, iu_th_avg, iu_hot_avg
   integer, allocatable :: iu_syn(:), iu_th(:), iu_hot(:)   ! one output file per map frequency
   character(len=:), allocatable :: tag_cd, tag_out
   logical :: hot

   Lsy_sh_unabs = 0.d0
   Lsy_sh_abs = 0.d0

   n_map = size(maps)
   allocate( iu_syn(n_map), iu_th(n_map), iu_hot(n_map) )
   iu_hot = -1
   iu_hot_avg = -1
   hot = (shock == 'BOTH') .or. (shock == 'FS')   ! the hot layer only exists in the FS

   ! File-name tags:
   !   tag_cd : files written by CD_abssyn   ('CD_sync_radiomap_RS110MHz.dat')
   !   tag_out: files written here           ('sync_radiomap_RS_110MHz.dat')
   if (shock == 'BOTH') then
      tag_cd  = ''
      tag_out = ''
   else
      tag_cd  = shock
      tag_out = shock//'_'
   end if

   ! Output maps (all phi together): one file per selected frequency + the average
   do im = 1, n_map
      open(newunit=iu_syn(im), file='sync_radiomap_'//tag_out//trim(maps(im)%label)//'.dat')
      open(newunit=iu_th(im),  file='th_ff_map_'//tag_out//trim(maps(im)%label)//'.dat')
      if (hot) open(newunit=iu_hot(im), file='hot_layer_'//trim(maps(im)%label)//'.dat')
   end do
   open(newunit=iu_syn_avg, file='sync_radiomap_'//tag_out//'avg.dat')
   open(newunit=iu_th_avg,  file='th_ff_map_'//tag_out//'avg.dat')
   if (hot) open(newunit=iu_hot_avg, file='hot_layer_avg.dat')

   phi = 0.d0
   if (rotate) then
      DO i_phi = 1, mphi
         print *, 'i(phi)=',i_phi,'/',mphi

         call CD_rot_run(phi)   

         if (abssyn_calc) then
            print *, 'Calculate the free-free synchrotron absorption'
            call CD_abssyn_run(shock, radiative = .true.) !(phi)
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

         if (maps_calc) then           
            print *, 'Save synthetic radio emission maps'
            do im = 1, n_map
               call merge_map('CD_sync_radiomap_'//tag_cd//trim(maps(im)%label)//'.dat', &
                              'CD_th_RS_'//trim(maps(im)%label)//'.dat',                   &
                              'CD_th_FS_'//trim(maps(im)%label)//'.dat',                   &
                              'CD_hot_layer_'//trim(maps(im)%label)//'.dat',               &
                              iu_syn(im), iu_th(im), iu_hot(im))
            end do
            call merge_map('CD_sync_radiomap_avg'//tag_cd//'.dat',   &
                           'CD_th_radiomap_avg_RS.dat',              &
                           'CD_th_radiomap_avg_FS.dat',              &
                           'CD_hot_radiomap_avg.dat',                &
                           iu_syn_avg, iu_th_avg, iu_hot_avg)
         endif

         phi = phi + delta_phi           
      end do
   end if

   do im = 1, n_map
      close(iu_syn(im))
      close(iu_th(im))
      if (hot) close(iu_hot(im))
   end do
   close(iu_syn_avg)
   close(iu_th_avg)
   if (hot) close(iu_hot_avg)

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
      do im = 1, n_map
         call SW_conv_run(shock, maps(im))
      end do
      if (map_avg%on) call SW_conv_run(shock, map_avg)
   endif

   deallocate( iu_syn, iu_th, iu_hot )

contains

   !-------------------------------------------------------------------------
   ! Reads the maps of one frequency for the current phi (synchrotron from
   ! CD_abssyn, thermal RS + FS and hot layer from shock_thermal) and appends
   ! them to the output files. For i_phi > 1 it also writes the mirror image
   ! (v -> -v), line by line as before.
   !-------------------------------------------------------------------------
   subroutine merge_map(f_syn, f_th_RS, f_th_FS, f_hot, iu_s, iu_t, iu_h)
      character(len=*), intent(in) :: f_syn, f_th_RS, f_th_FS, f_hot
      integer, intent(in) :: iu_s, iu_t, iu_h
      real(dp) :: uve(4,m_l)     ! u_proj, v_proj, emi, emi_abs
      real(dp) :: S_RS(m_l), S_FS(m_l), S_th(m_l), S_hot(m_l)
      integer :: iu, il

      open(newunit=iu, file=f_syn);   read(iu,*) uve;  close(iu)
      open(newunit=iu, file=f_th_RS); read(iu,*) S_RS; close(iu)
      open(newunit=iu, file=f_th_FS); read(iu,*) S_FS; close(iu)
      S_th = S_RS + S_FS
      if (hot) then
         open(newunit=iu, file=f_hot); read(iu,*) S_hot; close(iu)
      end if

      if (i_phi == 1) then
         write(iu_s,'(4es16.8)') uve
         write(iu_t,'(es16.8)') S_th
         if (hot) write(iu_h,'(es16.8)') S_hot
      else
         write(iu_s,'(4es16.8)') (uve(:,il), uve(1,il), -uve(2,il), uve(3:4,il), il=1,m_l)
         write(iu_t,'(es16.8)') (S_th(il), S_th(il), il=1,m_l)
         if (hot) write(iu_h,'(es16.8)') (S_hot(il), S_hot(il), il=1,m_l)
      end if
   end subroutine merge_map

END subroutine abs_maps