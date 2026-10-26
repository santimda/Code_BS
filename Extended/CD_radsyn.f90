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
  ! Intrinsic (unabsorbed) synchrotron emission of one linear emitter, using
  ! the electron distribution diste_<system>_<shock>.dat and the frequency
  ! grid in nu.dat (vec_nu). Includes the Razin-Tsytovich (R-T) suppression.
  ! Lsy(l,i) is the luminosity per unit Eps = h*nu/mec2 at position l [erg/s].
  ! Output: radsyn_<shock>.dat, S_syn_<system>_<shock>.dat, radsyn_tot_<shock>.dat
  ! and emisyn_<shock>.dat, which has Eps*Lsy(l,i) = nu*L_nu [erg/s] and is
  ! used by CD_abssyn for the absorption and the maps.
  !============================================================================

  subroutine CD_radsyn_run(shock)
    !$ use omp_lib
    character(LEN=*), intent(in) :: shock
    integer  :: l, i, il, im, iu, j_max
    real(dp) :: Eps_int, Eint, sum_sy, nuc_max, ne, nu_RT
    real(dp) :: x_sy(mE_e)
    integer,  allocatable :: i_map(:)
    real(dp), allocatable :: nu(:), Eps(:), dEps(:), Lsy(:,:), Lsy_tot(:), S_nu(:)
    real(dp), allocatable :: rho(:), B(:), B_par(:)
    real(dp), allocatable :: E_e(:), gamma_e(:), dE_e(:), N_etot(:,:)

    allocate( nu(m_nu), Eps(m_nu), dEps(m_nu), Lsy(m_theta,m_nu), Lsy_tot(m_nu), S_nu(m_nu) )
    allocate( rho(m_theta), B(m_theta), B_par(m_theta) )
    allocate( E_e(mE_e), gamma_e(mE_e), dE_e(mE_e), N_etot(m_theta,mE_e) )

    if (verbose) print *, 'Calculate the (unabsorbed) synchrotron emission in the '//shock

    ! set omp threads number
    !$ call omp_set_num_threads(nth)

    !////////////////////////////////////////////////////////////////////////
    !
    ! Read thermodynamical quantites at each position

    open(11,file=output_file('termo_'//shock//'.dat'))
    do l = 1, m_theta
       read(11,*) il,rho(l),dummy,dummy,dummy,B(l) ! CGS
    enddo
    close(11)
    B_par = B*par

    ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i)
    open (14,file=output_file('diste_'//system//'_'//shock//'.dat'))
    do l = 1, m_theta
       do i = 1, mE_e
          read(14,*) dummy, E_e(i),N_etot(l,i)
       enddo
       read(14,*)
    enddo
    close(14)
    E_e     = E_e*eV          ! energy from eV to erg
    Eint    = E_e(2)/E_e(1)
    gamma_e = E_e/mec2
    dE_e    = E_e*(Eint - 1.d0)

    ! Frequency grid (vec_nu)
    open(12,file=output_file('nu.dat'))
      read(12,*) nu
    close(12)
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

    ! Lsy(Eps) = int_Emin^Emax N(E_e) P(E_e,Eps) dE, with
    ! P = 4.39d-22*B_par*x**(1/3)*exp(-x)/h [erg/s/erg], x = nu/nu_c, nu_c = 4.22d6*B_par*gamma**2,
    ! neglected for x > 10. At each position we calculate the synchrotron emission.

    !$omp parallel do default(shared) schedule(dynamic) private(l, i, ne, nu_RT, x_sy)
    DO l = 1, m_theta
      ne    = 1.395d0*rho(l)/mp/2.085d0     ! 2.085 = x + 4*y + 14*z, 1.395 = x + 2*y + 7*z  ?????????
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
    S_nu    = Lsy_tot*(mphi*4.)/(4.*pi*distance**2*nu)*1.d26      ! mphi emitters, 4 quadrants [mJy]

    write(*,'(A19,X,Es10.3)') 'L_syn(unabsorbed)=',sum_sy*(mphi*4.) ! mphi emitters, 4 quadrants

    ! Fluxes at the radio-map frequencies (controls.f90)
    i_map = nu_index(maps%nu, nu(1), log(nu(2)/nu(1)), m_nu)
    do im = 1, size(maps)
       write(*,'(A,A,A,X,Es10.3)') ' S_', trim(maps(im)%label), '[mJy]=', S_nu(i_map(im))
    end do

    !---------------------------------------
    !
    !   Save data

    open(100, file='radsyn_tot_'//shock//'.dat')
      write(100,*) sum_sy*(mphi*4.)
    close(100)

    ! Total SED
    open(17,file='radsyn_'//shock//'.dat')
    write(17,'(3es25.16e3)') (Eps(i)*mec2/eV, Lsy_tot(i), Lsy_tot(i)/(Eps(i)*mec2), i=1,m_nu)
    close(17)

    ! Fluxes
    open(18,file='S_syn_'//system//'_'//shock//'.dat')
    write(18,'(2es25.16e3)') (nu(i), S_nu(i), i=1,m_nu)
    close(18)

    ! 3-D emission map
    open(newunit=iu,file='emisyn_'//shock//'.dat')
    do l = 1, m_theta
       write(iu,'(i6,2es25.16e3)') (l, Eps(i)*mec2/eV, Eps(i)*Lsy(l,i), i=1,m_nu)
       write(iu,*)
    enddo
    close(iu)

    deallocate( nu, Eps, dEps, Lsy, Lsy_tot, S_nu )
    deallocate( rho, B, B_par )
    deallocate( E_e, gamma_e, dE_e, N_etot )

  END subroutine CD_radsyn_run

End module CD_radsyn