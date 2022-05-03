module shock_thermal
  use global
  implicit none

  public :: shock_thermal_run_total, shock_thermal_run_spectrum

  contains

  !============================================================================
  ! shock_thermal_run_total (START)
  !============================================================================
  ! This module calculates the shock's total thermal bremsstrahlung luminosity
  !============================================================================

  subroutine shock_thermal_run_total(shock)
    character(LEN=*), intent(in) :: shock
    real(dp) :: ne,ni,n_a(m_theta),emi_ff
    real(dp) :: f_int(m_theta),dtheta,integral, L_tot
    real(dp) :: dV_ring(m_theta)
    integer :: nl
    real(dp), parameter :: g_average = 1.2 ! Rybicki, Eq. 5.15 b

    allocate( H_theta(m_theta), rho(m_theta),T(m_theta) )

    80  format(A12,X,Es10.2)

    print*, 'Calculate the total thermal luminosity from each shock'

    !------------------------------------------------------------
    !    Read shock position and thermodynamical quantities   

    open(13,file='Hnorm_'//shock//'.dat') 
    do l=1,m_theta
      read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l)
    enddo
    close(13)

    dtheta = theta(2) - theta(1) ! This is equally spaced

    open(11,file='termo_'//shock//'.dat')
    do l=1,m_theta
      read(11,*) nl,rho(l),dummy,dummy,T(l),dummy,dummy,dummy,dummy,dummy,dummy! CGS
      n_a(l) = rho(l)/mp/mu_i_sw   
    enddo
    close(11)

    !------------------------------------------------------------
    !    Calculate the shocked wind (sw) bremsstrahlung total power

    ! Calculate the luminosity using an adaptation of Eq. 24 from Christie+2016:
    ! L = integral(f_int(theta), with f_int(theta) = emi_ff * dV_ring(theta),
    ! emi_ff = dW/dt/dV the total emissivity given by Eq. 5.15b from Rybicki,
    ! and dV_ring(l) the volume of the emitting region

    do l=1,m_theta
      ni = n_a(l) 
      ne = gamma_sw*n_a(l)   
      emi_ff = 1.4d-27*Zq_sw**2*ne*ni*sqrt(T(l))*g_average  
      f_int(l) = emi_ff * dV_ring(l)
    enddo

    call integrate(m_theta,dtheta,f_int,integral)
    L_tot = integral
    write(*,80) 'L_ff_tot =',L_tot

    deallocate ( H_theta, rho, T )

  END subroutine shock_thermal_run_total

  !============================================================================
  ! shock_thermal_run_total (END)
  !============================================================================


  !============================================================================
  ! shock_thermal_run_spectrum (START)
  !============================================================================
  ! This module calculates the shock's thermal bremsstrahlung SED
  !============================================================================

  subroutine shock_thermal_run_spectrum(shock)
    character(LEN=*), intent(in) :: shock
    integer, parameter :: mE_ph = 256
    real(dp) :: n_a(m_theta),emi_ff_nu, cte_l(m_theta)
    real(dp) :: v_sh
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph), dEph
    real(dp) :: dV_ring(m_theta)
    real(dp) :: f_int(m_theta),dtheta,integral
    real(dp) :: LX_sh(mE_ph), L_tot
    integer :: nl, i_ph
    integer, parameter :: nu_mu=146,nu_mga=81 ! # data in table gauntff.dat
    real(dp) :: gamma2,dg,dmu,mu_ff,gff(nu_mu,nu_mga),gffinter
    integer :: j_mu,j_g
    interface
    subroutine interpol2(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
      real(dp), intent(out) :: finterpol
    end subroutine interpol2
    end interface
    interface
    subroutine extrapol2(dx,dy,fx1y1,fx2y1,fx1y2,fextrapol)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
      real(dp), intent(out) :: fextrapol
    end subroutine extrapol2
    end interface  

    allocate( H_theta(m_theta), rho(m_theta), T(m_theta) )

    print*, 'Calculate the spectrum of the thermal emission from each shock'

    LX_sh = 0.d0
    L_tot = 0.d0

    80  format(A12,X,Es10.2)

    if (shock == 'RS') then
      v_sh = v_w
    else if (shock == 'FS') then
      v_sh = v_ISM
    end if

    !------------------------------------------------------------
    ! Extract the Gaunt factor from the table given by van Hoof (2014):
    ! <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, with step 1

    open(13,file='gauntff.dat')
    do i=1,16
      read(13,*)
    enddo

    do i=1,nu_mu
      read(13,*) (gff(i,j), j=1,nu_mga)
    enddo
    close(13)

    !------------------------------------------------------------     
    !    Read shock position and thermodynamical quantities   

    open(13,file='Hnorm_'//shock//'.dat') 
    do l=1,m_theta
      read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l)
    enddo
    close(13)

    dtheta = theta(2) - theta(1) ! This is equally spaced

    open(11,file='termo_'//shock//'.dat')
    do l=1,m_theta
      read(11,*) nl,rho(l),dummy,dummy,T(l),dummy,dummy,dummy,dummy,dummy,dummy ! CGS
      n_a(l) = rho(l)/mp/mu_i_sw   
      !ni = n_a(l) 
      !ne = gamma_sw*n_a(l)    
      cte_l(l) = 6.8d-38*Zq_sw**2*gamma_sw*n_a(l)**2/sqrt(T(l))
    enddo
    close(11)

    !------------------------------------------------------------
    !    
    !    Calculate the shocked wind (sw) bremsstrahlung spectrum
    !
    !------------------------------------------------------------

    Eph_0 = mp*v_sh**2/32.d0   ! Use the k*T_sw(apex) as a reference
    Eph_min = 10.d6*h          ! 10 MHz
    Eph_max = Eph_0*4.d1

    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    ! Calculate the luminosity using an adaptation of Eq. 24 from Christie+2016:
    ! L = integral(f_int(theta), with f_int(theta) = emi_ff_nu * dV_ring(theta),
    ! emi_ff_nu= dW/dt/dV/dnu the specific emissivity given by Eq. 5.14b from Rybicki,
    ! and dV_ring(l) the volume of the shocked ring element

    open(19,file=''//shock//'_ff.dat')
    do i_ph=1,mE_ph

      do l=1,m_theta     
        ! Eq. 5.14b from Rybicki gives the spectrum in erg/s/cm3/Hz, so I divide by h
        ! emi_ff_nu = (6.8d-38*Zq_sw**2*ne*ni*T(l)**(-0.5d0)*exp(-Eph(i_ph)/(k*T(l)))*g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg

        ! Interpolate Gaunt factor
        gamma2 = Zq_sw**2.*Ry/k/T(l)
        j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
        dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

        mu_ff = Eph(i_ph)/k/T(l)
        j_mu = int( (log10(mu_ff)+16.d0)/0.2d0 ) + 1
        if (j_mu < nu_mga) then
          dmu = ((log10(mu_ff)+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
          call interpol2(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
        else
          dmu = ((log10(mu_ff)+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
          call extrapol2(dg,dmu,gff(j_g,nu_mga-1),gff(j_g+1,nu_mga-1),gff(j_g,nu_mga),gffinter)
        endif
        ! write(*,*) i_ph, l, gffinter, 9.77d0*(1.d0+0.13d0*log10(T(l)**1.5d0/Zq_sw/(Eph(i_ph)/h))) ! aprox 
        
        emi_ff_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * gffinter )/h ! erg/s/cm3/erg

        ! emi_ff_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg
        f_int(l) = emi_ff_nu*dV_ring(l)
      enddo

      call integrate(m_theta,dtheta,f_int,integral)
      LX_sh(i_ph) = integral
      write(19,*) Eph(i_ph)/eV, Eph(i_ph)*LX_sh(i_ph), LX_sh(i_ph)

      dEph = Eph(i_ph)*(Eph_int - 1.d0)
      L_tot = L_tot + LX_sh(i_ph)*dEph
    enddo
    close(19)

    write(*,80) 'L_tot =',L_tot

    deallocate ( H_theta, rho, T )

  END subroutine shock_thermal_run_spectrum

  !============================================================================
  ! shock_thermal_run_spectrum (END)
  !============================================================================

END module shock_thermal


!////////////////////////////////////////////////////////////////////////////////

Subroutine interpol2(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
  real(dp), intent(out) :: finterpol
  real(dp) :: fx,fy
  fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
  fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
  finterpol = fx1y1 + fx*dx + fy*dy
END Subroutine interpol2

!////////////////////////////////////////////////////////////////////////////////

Subroutine extrapol2(dx,dy,fx1y1,fx2y1,fx1y2,fextrapol)
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
  real(dp), intent(out) :: fextrapol
  real(dp) :: fx,fy
  fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
  fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
  fextrapol = fx1y2 + fx*dx + fy*dy
END Subroutine extrapol2

