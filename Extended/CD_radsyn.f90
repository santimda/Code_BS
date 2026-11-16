module CD_radsyn
  ! This module calculates the intrinsic synchrotron emission
  use global
  use controls, only: verbose
  use map_freqs, only: maps, nu_index

  implicit none
  private
  public :: CD_radsyn_run

  ! We save the total emission (i.e., integrated in theta) in radsyn.dat.
  ! The position dependent emission, Eps x L_sy(l,Eps_i), is stored in emisyn.dat.
  ! The latter is useful for correcting for absorption and producing emission maps.

contains

  !============================================================================
  ! CD_radsyn_run
  !============================================================================
  ! Synchrotron emission of one linear emitter, with Razin-Tsytovich
  ! suppression. SSA is negligible (T_b << kinetic T of the electrons).
  !============================================================================

  subroutine CD_radsyn_run(shock)
    character(LEN=*), intent(in) :: shock
    integer  :: l, i, im, iu, j_max
    real(dp) :: Eps_int, Eint, sum_sy, nuc_max, ne, nu_RT
    real(dp) :: x_sy(mE_e)        ! private copy per OpenMP thread (small)
    integer,  allocatable :: i_map(:)
    real(dp), allocatable :: nu(:), Eps(:), dEps(:), Lsy(:,:), Lsy_tot(:), S_nu(:)
    real(dp), allocatable :: rho(:), B(:), B_par(:)
    real(dp), allocatable :: E_e(:), gamma_e(:), dE_e(:), N_etot(:,:)

    allocate( nu(m_nu), Eps(m_nu), dEps(m_nu), Lsy(m_theta,m_nu), Lsy_tot(m_nu), S_nu(m_nu) )
    allocate( rho(m_theta), B(m_theta), B_par(m_theta) )
    allocate( E_e(mE_e), gamma_e(mE_e), dE_e(mE_e), N_etot(m_theta,mE_e) )

    call validate_shock(shock, 'CD_radsyn_run')
    if (verbose) print *, 'Calculate the (unabsorbed) synchrotron emission in the '//shock


    !////////////////////////////////////////////////////////////////////////
    !
    ! Read thermodynamical quantites at each position

    open(newunit=iu, file=output_file('thermo_'//shock//'.dat'))
    read(iu,*) (dummy, rho(l), dummy, dummy, dummy, B(l), dummy, dummy, dummy, dummy, dummy, dummy, &
                l=1,m_theta)    ! CGS, all the 12 columns
    close(iu)
    B_par = B*par

    ! Electron distribution N_etot(l,i)
    open(newunit=iu, file=output_file('diste_'//system//'_'//shock//'.dat'))
    do l = 1, m_theta
       read(iu,*) (dummy, E_e(i), N_etot(l,i), i=1,mE_e)
       read(iu,*)
    enddo
    close(iu)
    E_e     = E_e*eV          ! energy from eV to erg
    Eint    = E_e(2)/E_e(1)
    gamma_e = E_e/mec2
    dE_e    = E_e*(Eint - 1.d0)

    ! Frequency grid
    open(newunit=iu, file=output_file('nu.dat'))
    read(iu,*) nu
    close(iu)
    Eps     = nu*h/mec2       ! dimensionless photon energy
    Eps_int = Eps(2)/Eps(1)
    dEps    = Eps*(Eps_int - 1.d0)

    ! Check that the grid reaches 10*nu_c of the most energetic electrons
    j_max = findloc(any(N_etot > 0.d0, dim=1), .true., dim=1, back=.true.)
    if (j_max > 0) then
       nuc_max = 4.22d6*maxval(B_par)*gamma_e(j_max)**2
       if (10.d0*nuc_max > nu(m_nu)) print *, 'WARNING (CD_radsyn): nu_max_syn < 10*nu_c,max =', 10.d0*nuc_max
    end if

    !-----------------------------------------------------------------------!
    !                                                                       !
    !              Calculate Sycnhrotron emission                           !
    !                                                                       !
    !-----------------------------------------------------------------------!

    ! P = 4.39d-22*B_par*x**(1/3)*exp(-x)/h, x = nu/nu_c (neglected for x > 10)
    ! R-T suppression: exp(-nu_RT/nu), nu_RT = 20*n_e/B

    !$omp parallel do default(shared) schedule(dynamic) private(l, i, ne, nu_RT, x_sy)
    DO l = 1, m_theta
      ne    = rho(l)/(mu_e*mp)              ! electron density of the shocked gas
      nu_RT = 20.d0*ne/B(l)                 ! characteristic frequency for R-T effect
      do i = 1,m_nu                         ! Eps(i), nu(i)
        x_sy = nu(i)/(4.22d6*B_par(l)*gamma_e**2)
        Lsy(l,i) = mec2*4.39d-22*B_par(l)/h &
                   * sum( N_etot(l,:)*dE_e*merge(x_sy**(1.d0/3.d0)*exp(-x_sy), 0.d0, x_sy <= 10.d0) ) &
                   * exp(-nu_RT/nu(i))
      enddo
    ENDDO
    !$omp end parallel do

    Lsy_tot = Eps*sum(Lsy, dim=1)                                  ! [erg/s]
    sum_sy  = sum(matmul(Lsy, dEps))
    S_nu    = Lsy_tot*n_emit/(4.*pi*distance**2*nu)*1.d26      ! mphi emitters, 4 quadrants [mJy]

    write(*,'(A19,X,Es10.3)') 'L_syn(unabsorbed)=',sum_sy*n_emit ! mphi emitters, 4 quadrants

    ! Fluxes at the radio-map frequencies (controls.f90)
    i_map = nu_index(maps%nu, nu(1), log(nu(2)/nu(1)), m_nu)
    do im = 1, size(maps)
       write(*,'(A,A,A,X,Es10.3)') ' S_', trim(maps(im)%label), '[mJy]=', S_nu(i_map(im))
    end do

    !---------------------------------------
    !
    !   Save data

    open(newunit=iu, file=output_file('radsyn_tot_'//shock//'.dat'))
    write(iu,'(es17.8e3)') sum_sy*n_emit
    close(iu)

    ! Total SED (one linear emitter): E [eV], nu*L_nu [erg/s], L_E [1/s]
    open(newunit=iu, file=output_file('radsyn_'//shock//'.dat'))
    write(iu,'(3es25.16e3)') (Eps(i)*mec2/eV, Lsy_tot(i), Lsy_tot(i)/(Eps(i)*mec2), i=1,m_nu)
    close(iu)

    ! Fluxes of the whole shock: nu [Hz], S [mJy]
    open(newunit=iu, file=output_file('S_syn_'//system//'_'//shock//'.dat'))
    write(iu,'(2es25.16e3)') (nu(i), S_nu(i), i=1,m_nu)
    close(iu)

    ! Emission of each position: l, E [eV], nu*L_nu [erg/s]
    open(newunit=iu, file=output_file('emisyn_'//shock//'.dat'))
    do l = 1, m_theta
       write(iu,'(i6,2es25.16e3)') (l, Eps(i)*mec2/eV, Eps(i)*Lsy(l,i), i=1,m_nu)
       write(iu,*)
    enddo
    close(iu)

    deallocate( nu, Eps, dEps, Lsy, Lsy_tot, S_nu )
    deallocate( rho, B, B_par )
    deallocate( E_e, gamma_e, dE_e, N_etot )

  END subroutine CD_radsyn_run

END module CD_radsyn