module CD_radpp
  use global
  use controls, only: verbose
  use numerical_utils, only: vector_log
  implicit none
  private

  public :: CD_radpp_run

  ! delta-function approximation (Kelner et al. 2006, Sect. 4.3)
  real(dp), parameter :: K_pi = 0.17d0   ! fraction of the proton kinetic energy given to the pion
  real(dp), parameter :: n_pp = 0.86d0   ! n~ (depends on the proton spectrum; see Kelner+2006)

contains

  !============================================================================
  ! CD_radpp_run
  !============================================================================
  ! Intrinsic gamma-ray emission from pi0 decay (p-p) of one linear emitter,
  ! using the proton distribution distp_<system>_<shock>.dat (Kelner et al. 2006):
  !  - E_f > 10 GeV: Eq. 71, Phi(E_f) = c*n*int sigma(E_p) N(E_p) F(E_f/E_p,E_p) dE_p/E_p
  !  - E_f < 10 GeV: delta-function approximation, Eqs. 77-78, written as an
  !    integral over the proton energy (E_pi = K_pi*(E_p - mpc2), dE_pi = K_pi*dE_p):
  !    Q(E_f) = 2*n~*c*n * int sigma(E_p) N(E_p) dE_p / sqrt(E_pi**2 - m_pi**2),
  !    for E_pi > E_f + m_pi**2/(4*E_f)
  ! In both ranges Q(l,E_f) = c*n_a(l)*z_pp * sum_j N(l,j)*dE(j)*Kern(E_f,j), with a
  ! kernel that does not depend on the position: it is built once.
  ! Output (output/): radpp_<shock>.dat with E_f [eV] and E_f*L_E [erg/s] (one
  ! linear emitter) and radpp_tot_<shock>.dat with the total luminosity (x n_emit).
  !============================================================================

  subroutine CD_radpp_run(shock)
    character(LEN=*), intent(in) :: shock
    integer  :: l, i, j, iu
    real(dp) :: Ef_int, Eint, sum_pp, Epi_min, L1
    real(dp), allocatable :: rho(:), n_a(:)
    real(dp), allocatable :: E_p(:), dE(:), Epi(:), sig(:), N_ptot(:,:)
    real(dp), allocatable :: Ef(:), Kern(:,:), Lpp(:,:), Lpp_tot(:)

    allocate( rho(m_theta), n_a(m_theta) )
    allocate( E_p(mE_p), dE(mE_p), Epi(mE_p), sig(mE_p), N_ptot(m_theta,mE_p) )
    allocate( Ef(mE_f), Kern(mE_f,mE_p), Lpp(m_theta,mE_f), Lpp_tot(mE_f) )

    call validate_shock(shock, 'CD_radpp_run')
    if (verbose) print*, 'Calculate the p-p emission in the '//shock

    !------------------------------------------------------------------
    ! Read the density and the proton distribution

    open(newunit=iu, file=output_file('termo_'//shock//'.dat'))
    read(iu,*) (dummy, rho(l), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
                l=1,m_theta)    ! all the 12 columns
    close(iu)
    if (shock == 'RS') then
      n_a = rho/mp/mu_sw
    else
      n_a = rho/mp/mu_ISM
    end if

    ! Number of protons at r(l) with energy E_p(i): N_ptot(l,i)
    ! (one block of mE_p lines per position, separated by a blank line)
    open(newunit=iu, file=output_file('distp_'//system//'_'//shock//'.dat'))
    do l = 1, m_theta
      read(iu,*) (dummy, E_p(i), N_ptot(l,i), i=1,mE_p)
      read(iu,*)
    end do
    close(iu)
    E_p  = E_p*eV                     ! energy from eV to erg
    Eint = E_p(2)/E_p(1)
    dE   = E_p*(Eint-1.d0)

    ! Photon energy grid, from 10 MeV
    call vector_log(10.d6*eV, E_p(mE_p), Ef_int, Ef)

    !=======================================================================
    ! Kernel Kern(E_f,E_p) [cm2/erg]

    sig = [( sigma_pp(log(E_p(j)/TeV), E_p(j)), j=1,mE_p )]   ! zero below E_th
    Epi = K_pi*(E_p - mpc2)                                    ! pion energy (delta approx.)

    Kern = 0.d0
    do i = 1, mE_f
      if (Ef(i) < 0.01d0*TeV) then            ! delta-function approximation
        Epi_min = 1.005d0*( Ef(i) + m_pic2**2/(4.d0*Ef(i)) )   ! 1.005: avoid the integrable singularity
        do j = 1, mE_p
          if (Epi(j) >= Epi_min) Kern(i,j) = 2.d0*n_pp*sig(j)/sqrt(Epi(j)**2 - m_pic2**2)
        end do
      else                                    ! Kelner+2006, Eq. 71
        do j = 1, mE_p
          if (Ef(i) < E_p(j)) then
            L1 = log(E_p(j)/TeV)
            Kern(i,j) = sig(j)*F_gamma(Ef(i)/E_p(j), L1)/E_p(j)
          end if
        end do
      end if
    end do

    !=======================================================================
    ! L_E(l,E_f) = E_f*Q(l,E_f) [erg/s/erg]

    Lpp = spread(c*n_a*z_pp, 2, mE_f) * matmul(N_ptot*spread(dE,1,m_theta), transpose(Kern)) &
          * spread(Ef, 1, m_theta)
    Lpp_tot = sum(Lpp, dim=1)
    sum_pp  = sum(Lpp_tot*Ef)*(Ef_int-1.d0)

    write(*,'(A5,X,Es10.2)') 'L_pp=', sum_pp*n_emit     ! whole ring

    !---------------------------------------
    !   Save data

    open(newunit=iu, file=output_file('radpp_tot_'//shock//'.dat'))
    write(iu,'(es17.8e3)') sum_pp*n_emit
    close(iu)

    open(newunit=iu, file=output_file('radpp_'//shock//'.dat'))     ! one linear emitter
    write(iu,'(2es17.8e3)') (Ef(i)/eV, Ef(i)*Lpp_tot(i), i=1,mE_f)
    close(iu)

    deallocate( rho, n_a, E_p, dE, Epi, sig, N_ptot, Ef, Kern, Lpp, Lpp_tot )

  END subroutine CD_radpp_run

  !========================================

  elemental function F_gamma(x, L)
    ! Gamma-ray spectrum of a single p-p interaction, Kelner et al. (2006), Eq. 58;
    ! x = E_f/E_p, L = ln(E_p/1 TeV)
    real(dp), intent(in) :: x, L
    real(dp) :: F_gamma, Bg, beta, kg, xb
    Bg   = 1.30d0 + 0.14d0*L + 0.011d0*L**2
    beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L**2)
    kg   = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L**2)
    xb   = x**beta
    F_gamma = Bg*log(x)/x * ( (1.d0-xb)/(1.d0+kg*xb*(1.d0-xb)) )**4 &
            * ( 1.d0/log(x) - 4.d0*beta*xb/(1.d0-xb) - 4.d0*kg*beta*xb*(1.d0-2.d0*xb)/(1.d0+kg*xb*(1.d0-xb)) )
  end function F_gamma

END module CD_radpp