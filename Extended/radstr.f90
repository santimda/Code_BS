module radstr
   use global
   use controls, only: verbose
   use diffusion, only: Diff
   implicit none
   private
   public :: radstr_run

   real(dp), parameter :: B_str = 10.d-6                ! magnetic field in the HII region [G]
   character(len=*), parameter :: diff_regime = 'G'     ! Gabici+ 2006 (diffusion.f90)

contains

   !========================================================================
   !        radstr_run
   !========================================================================
   ! Order-of-magnitude estimate of the p-p and relativistic bremsstrahlung
   ! luminosities of the particles that escape from the BS into the HII
   ! region (Stromgren sphere). Only particles with E = E_break are considered.
   !  - Escape time from the HII region: 1/t_esc = 1/t_adv + 1/t_diff, with
   !    t_adv = R_str/v_st and t_diff = R_str**2/(6*D(E_break, B_str)).
   !  - Power injected into the HII region = power injected in the BS minus
   !    the power radiated in the BS:
   !      protons:   L_inj_p_HII = L_inj_p - L_pp
   !      electrons: L_inj_e_HII = L_inj_e - (L_br + L_IC,st + L_sy)
   !  - Fraction radiated before escaping: f = t_esc/(t_esc + t_loss)
   !    (-> t_esc/t_loss for t_esc << t_loss, and <= 1).
   ! Input files (output/): Stromgren.dat, L_inj_p_/L_inj_e_<shock>.dat (one
   ! linear emitter), radpp_tot_<shock>.dat, radbr_<shock>.dat (one linear
   ! emitter), rad_IC_st_tot_<shock>.dat and radsyn_tot_<shock>.dat (totals).
   ! Missing radiation files are taken as zero losses (with a warning).
   ! Output (output/): rad_pp_str_<shock>.dat and rad_Br_str_<shock>.dat,
   ! E [eV] and E*L_E [erg/s].
   !========================================================================

   subroutine radstr_run(shock)
   character(LEN=*), intent(in) :: shock
   integer  :: iu, j
   real(dp) :: R_str, D_diff, t_diff, t_adv, t_esc, t_pp, t_br, n_i
   real(dp) :: L_inj_p_BS, L_inj_e_BS, L_inj_p, L_inj_e
   real(dp) :: Lpp_sh, Lbr_sh, LIC_sh, Lsy_sh, L_pp_Etran, L_br_Etran
   real(dp), allocatable :: E_br(:), EL_br(:)

   80 format(A,X,Es10.2)
   81 format(A,X,Es10.2,A,Es10.2,A)

   call validate_shock(shock, 'radstr_run')
   if (verbose) print *, 'Estimate the emission of the escaping particles in the HII region'

   !------------------------------------------------------------
   ! Escape time from the HII region

   open(newunit=iu, file=output_file('Stromgren.dat'))
   read(iu,*) R_str
   close(iu)

   D_diff = Diff(E_break, B_str, diff_regime)
   t_diff = R_str**2/(6.d0*D_diff)
   t_adv  = R_str/v_st
   t_esc  = 1.d0/(1.d0/t_adv + 1.d0/t_diff)

   write(*,80) 'D_diff =', D_diff
   write(*,80) 't_diff =', t_diff
   write(*,80) 't_adv  =', t_adv
   write(*,80) 't_esc  =', t_esc

   n_i = rho_ISM/(mu_i*mp)            ! ion density of the HII region (n_ISM: total particle density)
   write(*,80) 'n_i    =', n_i

   !============================================================
   ! Protons

   open(newunit=iu, file=output_file('L_inj_p_'//shock//'.dat'))
   read(iu,*) L_inj_p_BS                         ! one linear emitter
   close(iu)

   Lpp_sh = read_total('radpp_tot_'//shock//'.dat')   ! already x n_emit

   L_inj_p = n_emit*L_inj_p_BS - Lpp_sh
   write(*,81) 'L_inj_p_HII =', L_inj_p, '  (L_inj_p =', n_emit*L_inj_p_BS, ')'
   if (L_inj_p < 0.d0) then
      print *, 'WARNING (radstr): L_pp in the BS > L_inj_p; L_inj_p_HII set to 0'
      L_inj_p = 0.d0
   end if

   t_pp = 1.d15/n_i                              ! p-p cooling time (order of magnitude)

   ! We estimate that approx 0.1*L_pp is radiated by E_break; protons with
   ! E = E_break radiate mostly at ~0.1*E_break
   L_pp_Etran = 0.1d0*t_esc/(t_esc + t_pp)*L_inj_p
   write(*,80) 'L_pp(0.1 E_break) =', L_pp_Etran

   open(newunit=iu, file=output_file('rad_pp_str_'//shock//'.dat'))
   write(iu,'(2es17.8e3)') 0.1d0*E_break/eV, L_pp_Etran
   close(iu)

   !============================================================
   ! Electrons

   open(newunit=iu, file=output_file('L_inj_e_'//shock//'.dat'))
   read(iu,*) L_inj_e_BS                         ! one linear emitter
   close(iu)

   ! Bremsstrahlung of one linear emitter: E [eV], E*L_E [erg/s] in a log grid
   Lbr_sh = 0.d0
   if (file_ok('radbr_'//shock//'.dat')) then
      allocate( E_br(mE_f), EL_br(mE_f) )
      open(newunit=iu, file=output_file('radbr_'//shock//'.dat'))
      read(iu,*) (E_br(j), EL_br(j), j=1,mE_f)
      close(iu)
      Lbr_sh = n_emit* sum(EL_br)*log(E_br(2)/E_br(1))   ! int L_E dE = int E*L_E dlnE
      deallocate( E_br, EL_br )
   end if

   LIC_sh = read_total('rad_IC_st_tot_'//shock//'.dat')   ! total (both halves, all phi)
   Lsy_sh = read_total('radsyn_tot_'//shock//'.dat')      ! already x n_emit

   L_inj_e = n_emit*L_inj_e_BS - (Lbr_sh + LIC_sh + Lsy_sh)
   write(*,81) 'L_inj_e_HII =', L_inj_e, '  (L_inj_e =', n_emit*L_inj_e_BS, ')'
   if (L_inj_e < 0.d0) then
      print *, 'WARNING (radstr): radiated power in the BS > L_inj_e; L_inj_e_HII set to 0'
      L_inj_e = 0.d0
   end if

   t_br = 1.d15/n_i                              ! relativistic bremsstrahlung cooling time (order of magnitude)
   write(*,80) 't_br   =', t_br

   ! Electrons with E = E_break radiate mostly at ~E_break
   L_br_Etran = t_esc/(t_esc + t_br)*L_inj_e
   write(*,80) 'L_br(E_break) =', L_br_Etran

   open(newunit=iu, file=output_file('rad_Br_str_'//shock//'.dat'))
   write(iu,'(2es17.8e3)') E_break/eV, L_br_Etran
   close(iu)

   contains

      ! The file exists in output/ (warning otherwise)
      logical function file_ok(name)
      character(len=*), intent(in) :: name
      inquire(file=output_file(name), exist=file_ok)
      if (.not. file_ok) print *, 'WARNING (radstr): '//name//' not found; that loss is taken as 0'
      end function file_ok

      ! First value of a file (0 if it does not exist)
      real(dp) function read_total(name)
      character(len=*), intent(in) :: name
      integer :: iu2
      read_total = 0.d0
      if (file_ok(name)) then
         open(newunit=iu2, file=output_file(name))
         read(iu2,*) read_total
         close(iu2)
      end if
      end function read_total


   end subroutine radstr_run

END module radstr