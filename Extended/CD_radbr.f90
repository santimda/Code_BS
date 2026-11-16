module CD_radbr
   use global
   use controls, only: verbose
   use numerical_utils, only: vector_log
   implicit none
   private

   public :: CD_radbr_run

contains

   !============================================================================
   ! CD_radbr_run
   !============================================================================
   ! Relativistic bremsstrahlung of one linear emitter (Blumenthal & Gould 1970),
   ! screened or unscreened nuclei; z_br from Padovani et al. (2018)
   !============================================================================

   subroutine CD_radbr_run(screened, shock)
      logical, intent(in) :: screened
      character(LEN=*), intent(in) :: shock
      real(dp) :: Ef_int, Eint, coef_sigma, sum_br
      integer  :: l, i, j, iu
      real(dp), allocatable :: rho(:), n_a(:)
      real(dp), allocatable :: E_e(:), gamma_e(:), dE_e(:), N_etot(:,:)
      real(dp), allocatable :: Ef(:), Phi(:,:), Lbr(:,:), Lbr_tot(:)

      allocate( rho(m_theta), n_a(m_theta) )
      allocate( E_e(mE_e), gamma_e(mE_e), dE_e(mE_e), N_etot(m_theta,mE_e) )
      allocate( Ef(mE_f), Phi(mE_f,mE_e), Lbr(m_theta,mE_f), Lbr_tot(mE_f) )

      call validate_shock(shock, 'CD_radbr_run')
      if (verbose) print*, '\n Calculate the (unabsorbed) rel. brem. emission in the '//shock

      !------------------------------------------------------------------------
      !  Read information (thermodynamical quantities and relativistic electrons)

      open(newunit=iu, file=output_file('thermo_'//shock//'.dat'))
      read(iu,*) (dummy, rho(l), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
                  l=1,m_theta)    ! all the 12 columns
      close(iu)
      if (shock == 'RS') then
         n_a = rho/mp/mu_sw       ! ambient number density
         coef_sigma = 4.d0*alpha_fs*re**2*Zq_sw**2*z_br
      else
         n_a = rho/mp/mu_ISM
         coef_sigma = 4.d0*alpha_fs*re**2*Zq_ISM**2*z_br
      end if

      ! Electron distribution N_etot(l,i)
      open(newunit=iu, file=output_file('diste_'//system//'_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu,*) (dummy, E_e(i), N_etot(l,i), i=1,mE_e)
         read(iu,*)
      end do
      close(iu)
      E_e     = E_e*eV            ! energy from eV to erg
      Eint    = E_e(2)/E_e(1)
      gamma_e = E_e/mec2
      dE_e    = E_e*(Eint-1.d0)

      ! Photon energy grid, from 1 keV
      call vector_log(1.d3*eV, E_e(mE_e), Ef_int, Ef)

      !===========================================================
      ! Relativistic bremsstrahlung spectrum, L_E [erg/s/erg]

      Phi = Phi_Br( spread(Ef,2,mE_e)/spread(E_e,1,mE_f), spread(gamma_e,1,mE_f), screened )
      Lbr = spread(n_a*c*coef_sigma, 2, mE_f) * matmul( N_etot*spread(dE_e,1,m_theta), transpose(Phi) )
      Lbr_tot = sum(Lbr, dim=1)
      sum_br  = sum(Lbr_tot*Ef)*(Ef_int-1.d0)

      write(*,'(A7,X,Es10.2)') 'L_Brem=', sum_br*n_emit   ! mphi emitters, 4 quadrants

      !---------------------------------------
      !   Save data (one linear emitter)

      open(newunit=iu, file=output_file('radbr_'//shock//'.dat'))
      write(iu,'(2es17.8e3)') (Ef(j)/eV, Ef(j)*Lbr_tot(j), j=1,mE_f)
      close(iu)

      deallocate( rho, n_a, E_e, gamma_e, dE_e, N_etot, Ef, Phi, Lbr, Lbr_tot )

   END subroutine CD_radbr_run

   !========================================

   elemental function Phi_Br(u,g_e,screened)
      ! Phi function of the rel. Br. cross-section (Blumenthal & Gould 1970); 0 for u >= 1
      real(dp) :: Phi_Br
      real(dp), intent(in) :: u,g_e
      logical, intent(in) :: screened

      if (u >= 1.d0) then
         Phi_Br = 0.d0
      else if (screened) then      ! screened nuclei
         Phi_Br = ( 1.d0 + (1.d0-u)**2 - (2.d0/3.d0)*(1.d0-u) )*log(191.d0/z_br**(1.d0/3.d0)) + (1.d0-u)/9.d0
      else                         ! un-screened nuclei
         Phi_Br = ( 1.d0 + (1.d0-u)**2 - (2.d0/3.d0)*(1.d0-u) )*( log(2.d0*g_e*(1.d0/u-1.d0)) - 0.5d0 )
      endif
      Phi_Br = max(Phi_Br, 0.d0)
   End function Phi_Br

END module CD_radbr