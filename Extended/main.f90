PROGRAM main
   use global
   use controls
   use initialize
   use CD, only: CD_run
   use BB, only: BB_run
   use opac_ff, only: opac_ff_run
   use map_freqs, only: maps_init
   !$ use omp_lib
   implicit none

   ! FS: .true. = radiative (isothermal layer, no NT particles), .false. = adiabatic
   logical, parameter :: radiative_FS = .true.
   integer :: st

   interface
      subroutine calculate(shock, radiative)
         character(LEN=*), intent(in) :: shock
         logical, intent(in) :: radiative
      end subroutine calculate
      subroutine abs_maps(shock, radiative)
         character(LEN=*), intent(in) :: shock
         logical, intent(in) :: radiative
      end subroutine abs_maps
      subroutine H_alpha_map(shock)
         character(LEN=*), intent(in) :: shock
      end subroutine H_alpha_map
   end interface

   print *, '\n ============================'
   print *, ' START'
   print *, '============================\n '

   call create_output_dir()

   ! OpenMP threads: OMP_NUM_THREADS or nth
   call get_environment_variable('OMP_NUM_THREADS', status=st)
   !$ if (st /= 0) call omp_set_num_threads(nth)
   !$ if (verbose) write(*,'(a,i0,a)') ' Using OpenMP with ', omp_get_max_threads(), ' threads'

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
      if (maps_calc .and. .not. BOTH) then
         print *, '\n--------------------------------------'
         print *, ' Calculate absorption and emission maps'
         print *, '----------------------------------------'
         call abs_maps('RS', radiative_FS)
      end if
   endif

   if (FS) then
      print *, '\n----------------------------'
      print *, ' Forward shock'
      print *, '----------------------------'
      call calculate('FS', radiative=radiative_FS)
      if (maps_calc .and. .not. BOTH) then
         print *, '\n--------------------------------------'
         print *, ' Calculate absorption and emission maps'
         print *, '----------------------------------------'
         call abs_maps('FS', radiative_FS)
      end if
   endif

   if (maps_calc .and. BOTH) then
      print *, '\n--------------------------------------'
      print *, ' Calculate absorption and emission maps'
      print *, '----------------------------------------'
      call abs_maps('BOTH', radiative_FS)
   end if

   if (rad_alpha .and. FS) then
      if (radiative_FS) then
         call H_alpha_map('FS')
      else
         print *, 'WARNING (main): the H-alpha map needs a radiative FS'
      end if
   end if

   if (verbose) print *, '\n Deallocate vectors'
   call initialize_deallocate()

   print *, '\n ======================================'
   print *, '         END'
   print *, '======================================'

END PROGRAM main


!============================================================================
! Weight of each phi in the sums over phi: the mirror (z -> -z) of
! phi = 0 and pi/2 gives the same positions
!============================================================================

pure real(kind(1.d0)) function phi_weight(i_phi, mphi)
   integer, intent(in) :: i_phi, mphi
   phi_weight = 1.d0
   if (i_phi /= 1 .and. i_phi /= mphi) phi_weight = 2.d0
end function phi_weight


!============================================================================

Subroutine calculate(shock,radiative)
   use global
   use controls
   use CD_termo, only: CD_termo_run
   use shock_thermal, only: shock_thermal_run
   use CD_diste, only: CD_diste_run
   use CD_distp, only: CD_distp_run
   use CD_radsyn, only: CD_radsyn_run
   use CD_radbr, only: CD_radbr_run
   use CD_radic, only: CD_radic_st_run, CD_radic_IR_run
   use CD_radpp, only: CD_radpp_run
   use CD_rot, only: CD_rot_run
   use radstr, only: radstr_run
   implicit none

   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative
   interface
      pure real(kind(1.d0)) function phi_weight(i_phi, mphi)
         integer, intent(in) :: i_phi, mphi
      end function phi_weight
   end interface
   real(dp) :: LIC_int, LIC_all, w, phi
   integer  :: i_phi, i_E, iu, iu2
   real(dp), allocatable :: E_IC_st(:), LIC_one(:), LIC_st_tot(:), E_IC_IR(:), LIC_IR_tot(:)
   real(dp), allocatable :: E_sy(:), Lsy_sh_tot(:), L_syn_E_tot(:), E_br(:), Lbr_sh_tot(:)
   real(dp), allocatable :: E_pp(:), Lpp_sh_tot(:)

   80 format(A,X,Es10.2)

   allocate( E_IC_st(mE_f), LIC_one(mE_f), LIC_st_tot(mE_f), E_IC_IR(mE_f), LIC_IR_tot(mE_f) )
   allocate( E_sy(m_nu), Lsy_sh_tot(m_nu), L_syn_E_tot(m_nu), E_br(mE_f), Lbr_sh_tot(mE_f) )
   allocate( E_pp(mE_f), Lpp_sh_tot(mE_f) )
   LIC_st_tot = 0.d0

   if (termo_calc) call CD_termo_run(shock, radiative)

   if (Xrad_calc) call shock_thermal_run(shock, radiative)

   !===============================================================!
   !                                                               !
   !                      Electrons                                !
   !                                                               !
   !===============================================================!

   IF (electrons .and. (.not. radiative)) then
      print *, '\n ELECTRONS:'

      if (diste_calc) call CD_diste_run(shock)

      if (radsyn_sw_calc) call CD_radsyn_run(shock)       ! Just one linear emitter

      if (radbr_sw_calc) call CD_radbr_run(unscr,shock)   ! Just one linear emitter

      if (radic_ir_calc) call CD_radic_ir_run(shock)      ! Just one linear emitter

      ! IC-star emission depends on the 3D position: linear emitters for different phi
      phi = 0.d0
      LIC_all = 0.d0
      if (rotate .and. radic_st_calc) then
         DO i_phi = 1, mphi
            print *, 'i(phi)=',i_phi,'/',mphi

            call CD_rot_run(phi)

            print *, 'Calculate the IC(star) emission'
            call CD_radic_st_run(LIC_int,shock)
            open(newunit=iu, file=output_file('radic_st_'//shock//'.dat'))
            read(iu,*) (E_IC_st(i_E), LIC_one(i_E), i_E=1,mE_f)
            close(iu)

            w = phi_weight(i_phi, mphi)
            LIC_st_tot = LIC_st_tot + w*LIC_one
            LIC_all    = LIC_all    + w*LIC_int

            phi = phi + delta_phi
         end do

         write(*,80) 'L_IC_st =', LIC_all
         open(newunit=iu, file=output_file('rad_IC_st_tot_'//shock//'.dat'))
         write(iu,'(es17.8e3)') LIC_all                 ! total (whole BS)
         close(iu)
      endif

      !=============================================================
      !  Save the totals (whole BS)

      if (escribir_e) then

         if (radic_st_calc .and. rotate) then           ! this has the full 3D geometry
            open(newunit=iu, file=output_file('radic_st_'//system//'_'//shock//'_tot.dat'))
            write(iu,'(2es17.8e3)') (E_IC_st(i_E), LIC_st_tot(i_E), i_E=1,mE_f)
            close(iu)
         endif

         if (radic_ir_calc) then
            open(newunit=iu, file=output_file('radic_ir_'//shock//'.dat'))     ! One linear emitter
            read(iu,*) (E_IC_IR(i_E), LIC_IR_tot(i_E), i_E=1,mE_f)
            close(iu)
            LIC_IR_tot = LIC_IR_tot*n_emit                                    ! whole ring
            open(newunit=iu, file=output_file('radic_ir_'//system//'_'//shock//'_tot.dat'))
            write(iu,'(2es17.8e3)') (E_IC_IR(i_E), LIC_IR_tot(i_E), i_E=1,mE_f)
            close(iu)
         endif

         if (radsyn_sw_calc) then
            open(newunit=iu, file=output_file('radsyn_'//shock//'.dat'))        ! One linear emitter
            read(iu,*) (E_sy(i_E), Lsy_sh_tot(i_E), L_syn_E_tot(i_E), i_E=1,m_nu)
            close(iu)
            Lsy_sh_tot  = Lsy_sh_tot*n_emit
            L_syn_E_tot = L_syn_E_tot*n_emit
            open(newunit=iu, file=output_file('radsyn_'//system//'_'//shock//'_tot.dat'))
            write(iu,'(3es17.8e3)') (E_sy(i_E), Lsy_sh_tot(i_E), L_syn_E_tot(i_E), i_E=1,m_nu)
            close(iu)
         endif

         if (radbr_sw_calc) then
            open(newunit=iu, file=output_file('radbr_'//shock//'.dat'))         ! One linear emitter
            read(iu,*) (E_br(i_E), Lbr_sh_tot(i_E), i_E=1,mE_f)
            close(iu)
            Lbr_sh_tot = Lbr_sh_tot*n_emit
            open(newunit=iu, file=output_file('radbr_'//system//'_'//shock//'_tot.dat'))
            write(iu,'(2es17.8e3)') (E_br(i_E), Lbr_sh_tot(i_E), i_E=1,mE_f)
            close(iu)
         endif

         if (radic_st_calc .and. rotate .and. radic_ir_calc .and. radsyn_sw_calc .and. Xrad_calc) then
            call X_ray_luminosity()
         end if

      endif

   ENDIF

   !===============================================================!
   !                                                               !
   !                      Protons                                  !
   !                                                               !
   !===============================================================!

   IF (protons .and. (.not. radiative)) then
      print *, '\n PROTONS:'

      if (distp_calc) call CD_distp_run(shock)

      if (radpp_calc) call CD_radpp_run(shock)

      if (escribir_p .and. radpp_calc) then
         open(newunit=iu, file=output_file('radpp_'//shock//'.dat'))            ! One linear emitter
         read(iu,*) (E_pp(i_E), Lpp_sh_tot(i_E), i_E=1,mE_f)
         close(iu)
         Lpp_sh_tot = Lpp_sh_tot*n_emit                                        ! whole ring
         open(newunit=iu2, file=output_file('radpp_'//system//'_'//shock//'_tot.dat'))
         write(iu2,'(2es17.8e3)') (E_pp(i_E), Lpp_sh_tot(i_E), i_E=1,mE_f)
         close(iu2)
      endif

   ENDIF

   !================================================================

   if (stromgren .and. electrons .and. protons .and. (.not. radiative)) then
      print*, '\n Calculate the NT emission from the H II region:'
      call radstr_run(shock)
   end if

   deallocate( E_IC_st, LIC_one, LIC_st_tot, E_IC_IR, LIC_IR_tot )
   deallocate( E_sy, Lsy_sh_tot, L_syn_E_tot, E_br, Lbr_sh_tot, E_pp, Lpp_sh_tot )

contains

   !-------------------------------------------------------------------------
   ! X-ray luminosity: IC (star and dust), synchrotron and thermal
   !-------------------------------------------------------------------------
   subroutine X_ray_luminosity()
      real(dp) :: Eph_min, Eph_max, L_X_int, L_X_max
      real(dp), allocatable :: E_ff(:), Lff_tot(:)

      if (system == 'BD+43') then
         Eph_min = 0.4d3              ! 0.4 keV
         Eph_max = 4.d3               ! 4 keV
         L_X_max = 10.d0**30.878384d0 ! Toala+ 2016
      else if (system == 'Lambda_Cep') then
         Eph_min = 2.d3               ! 2 keV
         Eph_max = 10.d3              ! 10 keV
         L_X_max = 3.37d30
      else
         Eph_min = 1.d3               ! 1 keV
         Eph_max = 10.d3              ! 10 keV
         L_X_max = huge(1.d0)
      end if

      allocate( E_ff(m_nu), Lff_tot(m_nu) )
      open(newunit=iu, file=output_file('ionized_ff_'//system//'_'//shock//'.dat'))
      read(iu,*) (E_ff(i_E), Lff_tot(i_E), dummy, i_E=1,m_nu)
      close(iu)

      L_X_int = band(E_IC_st, LIC_st_tot, Eph_min, Eph_max) + band(E_IC_IR, LIC_IR_tot, Eph_min, Eph_max) &
              + band(E_sy, Lsy_sh_tot, Eph_min, Eph_max)    + band(E_ff, Lff_tot, Eph_min, Eph_max)

      write(*,'(a,f5.1,a,f5.1,a)') ' Luminosity between ', Eph_min/1.d3, ' and ', Eph_max/1.d3, ' keV'
      write(*,'(A,X,Es10.2)') 'L_X =', L_X_int
      if (L_X_int > L_X_max) print *, 'X-ray luminosity is too high'
      print *, ''

      deallocate( E_ff, Lff_tot )
   end subroutine X_ray_luminosity

   ! int L_E dE between emin and emax [eV]
   real(dp) function band(E, EL, emin, emax)
      real(dp), intent(in) :: E(:), EL(:), emin, emax
      band = sum(EL, mask = (E >= emin) .and. (E <= emax)) * log(E(2)/E(1))
   end function band

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
   real(dp) :: phi
   integer  :: i_phi, il, iu, iu_map
   real(dp), allocatable :: uvS(:,:)      ! u, v, S of each linear emitter

   ! phi = 0 and pi/2 are written once, the rest also mirrored (v -> -v)
   allocate( uvS(3,2*m_theta) )
   open(newunit=iu_map, file=output_file('S_alpha_map_'//shock//'.dat'))
   phi = 0.d0
   do i_phi = 1, mphi
      print *, 'i(phi)=',i_phi,'/',mphi
      call CD_rot_run(phi)
      call project_H_alpha_map(shock)       ! This projects the H_alpha map in u,v
      open(newunit=iu, file=output_file('S_alpha_proj_'//shock//'.dat'))
      read(iu,*) uvS
      close(iu)
      if (i_phi == 1 .or. i_phi == mphi) then
         write(iu_map,'(3es17.8e3)') uvS
      else
         write(iu_map,'(3es17.8e3)') (uvS(:,il), uvS(1,il), -uvS(2,il), uvS(3,il), il=1,2*m_theta)
      end if
      phi = phi + delta_phi
   end do
   close(iu_map)
   deallocate( uvS )

   if (convolution_alpha) call rad_H_alpha_run(shock)

End Subroutine H_alpha_map

!============================================================================
!============================================================================

Subroutine abs_maps(shock, radiative)
   use global
   use controls
   use map_freqs, only: maps
   use CD_abssyn, only: CD_abssyn_run
   use CD_rot, only: CD_rot_run
   use SW_conv, only: SW_conv_run
   implicit none

   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative            ! radiative FS (isothermal layer)
   interface
      pure real(kind(1.d0)) function phi_weight(i_phi, mphi)
         integer, intent(in) :: i_phi, mphi
      end function phi_weight
   end interface
   integer, parameter :: m_l = 2*m_theta     ! CD position, both sides
   real(dp) :: phi, w
   integer  :: i_phi, i_E, im, n_map, iu
   integer  :: iu_syn_avg, iu_th_avg, iu_hot_avg
   integer, allocatable :: iu_syn(:), iu_th(:), iu_hot(:)   ! one output file per map frequency
   character(len=:), allocatable :: tag_cd, tag_out
   logical :: hot
   real(dp), allocatable :: E_sy(:), nu_E(:), Lsy_phi(:), Lsy_phi_abs(:), Lsy_sh_unabs(:), Lsy_sh_abs(:)

   allocate( E_sy(m_nu), nu_E(m_nu), Lsy_phi(m_nu), Lsy_phi_abs(m_nu), Lsy_sh_unabs(m_nu), Lsy_sh_abs(m_nu) )
   Lsy_sh_unabs = 0.d0
   Lsy_sh_abs   = 0.d0

   n_map = size(maps)
   allocate( iu_syn(n_map), iu_th(n_map), iu_hot(n_map) )
   iu_hot = -1
   iu_hot_avg = -1
   hot = radiative .and. ((shock == 'BOTH') .or. (shock == 'FS'))   ! the hot layer only exists in the FS

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
      open(newunit=iu_syn(im), file=output_file('sync_radiomap_'//tag_out//trim(maps(im)%label)//'.dat'))
      open(newunit=iu_th(im),  file=output_file('th_ff_map_'//tag_out//trim(maps(im)%label)//'.dat'))
      if (hot) open(newunit=iu_hot(im), file=output_file('hot_layer_'//trim(maps(im)%label)//'.dat'))
   end do
   open(newunit=iu_syn_avg, file=output_file('sync_radiomap_'//tag_out//'avg.dat'))
   open(newunit=iu_th_avg,  file=output_file('th_ff_map_'//tag_out//'avg.dat'))
   if (hot) open(newunit=iu_hot_avg, file=output_file('hot_layer_avg.dat'))

   phi = 0.d0
   if (rotate) then
      DO i_phi = 1, mphi
         print *, 'i(phi)=',i_phi,'/',mphi

         call CD_rot_run(phi)

         if (abssyn_calc) then
            print *, 'Calculate the free-free synchrotron absorption'
            call CD_abssyn_run(shock, radiative=radiative, E_sy=E_sy, Lsy_sum=Lsy_phi, Lsy_sum_abs=Lsy_phi_abs)
            w = phi_weight(i_phi, mphi)
            Lsy_sh_unabs = Lsy_sh_unabs + w*Lsy_phi        ! SED --> sum emission from different phi
            Lsy_sh_abs   = Lsy_sh_abs   + w*Lsy_phi_abs
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
            call merge_map('CD_sync_radiomap_avg'//tag_cd//'.dat', &
                           'CD_th_radiomap_avg_RS.dat',            &
                           'CD_th_radiomap_avg_FS.dat',            &
                           'CD_hot_radiomap_avg.dat',              &
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

   if (abssyn_calc .and. rotate) then       ! This has the full 3D geometry in the abs
      nu_E = E_sy*eV/h
      open(newunit=iu, file=output_file('radsynabs_'//shock//'_tot.dat'))      ! E [eV], nu*L_nu unabs, abs [erg/s]
      write(iu,'(3es17.8e3)') (E_sy(i_E), Lsy_sh_unabs(i_E), Lsy_sh_abs(i_E), i_E=1,m_nu)
      close(iu)
      open(newunit=iu, file=output_file('S_syn_abs_'//system//'_'//shock//'.dat'))   ! nu [Hz], S_unabs, S_abs [mJy]
      write(iu,'(3es17.8e3)') (nu_E(i_E), Lsy_sh_unabs(i_E)/(4.d0*pi*distance**2*nu_E(i_E))*1.d26, &
                                          Lsy_sh_abs(i_E)  /(4.d0*pi*distance**2*nu_E(i_E))*1.d26, i_E=1,m_nu)
      close(iu)
   endif

   if (convolution) then
      print *, 'Convolution of the emission map with the instrument response'
      do im = 1, n_map
         call SW_conv_run(shock, maps(im), radiative)
      end do
      if (map_avg%on) call SW_conv_run(shock, map_avg, radiative)
   endif

   deallocate( iu_syn, iu_th, iu_hot )
   deallocate( E_sy, nu_E, Lsy_phi, Lsy_phi_abs, Lsy_sh_unabs, Lsy_sh_abs )

contains

   !-------------------------------------------------------------------------
   ! Appends the maps of the current phi (and their mirror image, except
   ! for phi = 0 and pi/2)
   !-------------------------------------------------------------------------
   subroutine merge_map(f_syn, f_th_RS, f_th_FS, f_hot, iu_s, iu_t, iu_h)
      character(len=*), intent(in) :: f_syn, f_th_RS, f_th_FS, f_hot
      integer, intent(in) :: iu_s, iu_t, iu_h
      real(dp) :: uve(4,m_l)     ! u_proj, v_proj, emi, emi_abs
      real(dp) :: S_RS(m_l), S_FS(m_l), S_th(m_l), S_hot(m_l)
      integer :: iu1, il

      open(newunit=iu1, file=output_file(f_syn));   read(iu1,*) uve;  close(iu1)
      open(newunit=iu1, file=output_file(f_th_RS)); read(iu1,*) S_RS; close(iu1)
      open(newunit=iu1, file=output_file(f_th_FS)); read(iu1,*) S_FS; close(iu1)
      S_th = S_RS + S_FS
      S_hot = 0.d0
      if (hot) then
         open(newunit=iu1, file=output_file(f_hot)); read(iu1,*) S_hot; close(iu1)
      end if

      if (i_phi == 1 .or. i_phi == mphi) then
         write(iu_s,'(4es17.8e3)') uve
         write(iu_t,'(es17.8e3)') S_th
         if (hot) write(iu_h,'(es17.8e3)') S_hot
      else
         write(iu_s,'(4es17.8e3)') (uve(:,il), uve(1,il), -uve(2,il), uve(3:4,il), il=1,m_l)
         write(iu_t,'(es17.8e3)') (S_th(il), S_th(il), il=1,m_l)
         if (hot) write(iu_h,'(es17.8e3)') (S_hot(il), S_hot(il), il=1,m_l)
      end if
   end subroutine merge_map

END subroutine abs_maps