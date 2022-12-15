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

  subroutine shock_thermal_run_total(shock,radiative)
    character(LEN=*), intent(in) :: shock
    logical, intent(in) :: radiative
    real(dp) :: ne,ni,n_a(m_theta),emi_ff
    real(dp) :: f_int(m_theta),dtheta,integral, L_tot
    real(dp) :: dV_ring(m_theta)
    integer :: nl
    real(dp), parameter :: g_average = 1.2 ! Rybicki, Eq. 5.15 b

    allocate( H_theta(m_theta), rho(m_theta),T(m_theta), H_1(m_theta) )

    80  format(A12,X,Es10.2)

    print*, 'Calculate the total thermal luminosity from each shock'

    !------------------------------------------------------------
    !    Read shock position and thermodynamical quantities   

    open(13,file='Hnorm_'//shock//'.dat') 
    do l=1,m_theta
      if ( .not. radiative ) then
        read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l)
      else if ( radiative ) then
        read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l), H_1(l)
      end if
    enddo
    close(13)

    dtheta = theta(2) - theta(1) ! This is equally spaced

    open(11,file='termo_'//shock//'.dat')
    do l=1,m_theta
      read(11,*) nl,rho(l),dummy,dummy,T(l),dummy,dummy,dummy,dummy,dummy,dummy,dummy! CGS
      if (shock == 'RS') then
        n_a(l) = rho(l)/mp/mu_i_sw
      else if (shock == 'FS') then
        n_a(l) = rho(l)/mp/mu_ISM
      end if
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

    deallocate ( H_theta, rho, T, H_1 )

  END subroutine shock_thermal_run_total

  !============================================================================
  ! shock_thermal_run_total (END)
  !============================================================================


  !============================================================================
  ! shock_thermal_run_spectrum (START)
  !============================================================================
  ! This module calculates the shock's thermal bremsstrahlung SED
  !============================================================================

  subroutine shock_thermal_run_spectrum(shock, radiative)
    character(LEN=*), intent(in) :: shock
    logical, intent(in) :: radiative
    integer, parameter :: mE_ph = 750
    integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
    real(dp) :: n_a(m_theta), n_1(m_theta), cte_l(m_theta), cte_l_1(m_theta)
    real(dp) :: v_sh
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph), dEph
    real(dp) :: nu(mE_ph)
    real(dp) :: log_nu_int,flux,flux0,flux1,flux2,flux3
    real(dp) :: dV_ring(m_theta)
    real(dp) :: integral,integ,integ_ion,integ_HL
    real(dp) :: A(m_theta)
    real(dp) :: emi_ff_nu(m_theta,mE_ph), emi_ion_nu
    real(dp) :: alpha_ff_nu1(m_theta,mE_ph), alpha_ff_nu2(m_theta,mE_ph)
    real(dp) :: Tau_1(m_theta,mE_ph), Tau_2(m_theta,mE_ph)
    real(dp) :: L_inj(m_theta)
    real(dp) :: f_int(m_theta),f_vol_ion(m_theta),dtheta
    real(dp) :: dV_ring_ff, H_ff, H_RH, dR(m_theta)
    real(dp) :: r(m_theta), v_perp(m_theta), t_th(m_theta)
    real(dp) :: F_nu_iph(m_l,mE_ph), F_nu_tot(m_theta), F_nu_ion_tot(m_theta)
    real(dp) :: F_nu_avg(m_theta), F_nu_ion_avg(m_theta)
    real(dp) :: dnu_avg, h_avg
    real(dp), parameter :: T_f = T_ISM
    !real(dp) :: F_nu(m_l,mE_ph),F_nu_phi(m_l,mE_ph)
    real(dp) :: n_BB(m_theta,mE_ph), T_min(m_theta)
    real(dp) :: L_tot,L_th(mE_ph),L_ion(mE_ph),L_ion_tot
    real(dp) :: S(mE_ph),S1,S2,S3,S_ion(mE_ph)
    integer :: nl, i_ph
    Integer :: i_890M,i_1G,i_3G,i_5G
    integer :: N_avg
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

    allocate( H_theta(m_theta), H_1(m_theta), rho(m_theta), T(m_theta), rho_f(m_theta) )
    allocate( F_nu(m_l,mE_ph),F_nu_phi(m_l,mE_ph),F_nu_ion(m_l,mE_ph),F_nu_phi_ion(m_l,mE_ph))

  90  format(A19,X,Es10.3)
  95  format(A18,X,Es10.3)
  96  format(A36,X,F4.0,X,A13,X,F4.0,X,A3)

    print*, 'Calculate the spectrum of the thermal emission from each shock'

    L_ion_tot = 0.d0
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

    open(12,file='gauntff.dat')
    do i=1,16
      read(12,*)
    enddo

    do i=1,nu_mu
      read(12,*) (gff(i,j), j=1,nu_mga)
    enddo
    close(12)

    !------------------------------------------------------------     
    ! Read shock position and thermodynamical quantities  

    open(10,file='CD.dat')
      do l = 1, m_theta
         read(10,*) theta(l), r_norm(l), dR(l), dummy, dummy
      end do
    close(10) 

    open(13,file='Hnorm_'//shock//'.dat') 
    do l=1,m_theta
      if ( .not. radiative ) then
        read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l)
      else if (radiative) then
        read(13,*) dummy, H_theta(l), dummy, dummy, dV_ring(l), H_1(l)
      end if
    enddo
    close(13)

    dtheta = theta(2) - theta(1) ! This is equally spaced

    open(11,file='termo_'//shock//'.dat')
    do l=1,m_theta
      read(11,*) nl,rho(l),dummy,dummy,T(l),dummy,v_perp(l),dummy,dummy,r(l),dummy,rho_f(l) ! CGS
      if (shock == 'FS') then
        n_1(l) = rho(l)/mp/mu_ISM
        n_a(l) = rho_f(l)/mp/mu_ISM  
        !ni = n_a(l) 
        !ne = gamma_sw*n_a(l)
        cte_l_1(l) = 6.8d-38*Zq_sw**2*gamma_sw*n_1(l)**2/sqrt(T(l)) 
        cte_l(l) = 6.8d-38*Zq_sw**2*gamma_sw*n_a(l)**2/sqrt(T_f)
      else if (shock == 'RS') then
        n_a(l) = rho(l)/mp/mu_sw
        cte_l(l) = 6.8d-38*Zq_sw**2*gamma_sw*n_a(l)**2/sqrt(T(l))
      end if
    enddo
    close(11)

    open(14,file='Linj_'//shock//'.dat')
    do l=1,m_theta
      read(14,*) dummy, dummy, L_inj(l), dummy, dummy
    end do
    close(14)

    open(15,file='T_min_'//shock//'.dat')
    do l=1,m_theta
      read(15,*) dummy, T_min(l)
    end do
    close(15)

    open(16,file='t_th_'//shock//'.dat')
    do l=1,m_theta
      read(16,*) dummy, t_th(l)
    end do
    close(16)
     

    !------------------------------------------------------------
    !    
    !    Calculate the shocked wind (sw) bremsstrahlung spectrum
    !
    !------------------------------------------------------------

    Eph_0 = mp*v_sh**2/32.d0   ! Use the k*T_sw(apex) as a reference
    Eph_min = 10.d6*h          ! 10 MHz
    Eph_max = Eph_0*1.d2


    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    !------------------------------------------------------------
    do i_ph=1,mE_ph
      nu(i_ph) = Eph(i_ph)/h
    end do
    
    log_nu_int = log(nu(2)/nu(1))

    flux = 1.d26/(4.d0*pi*(distance)**2) ! luminosity to flux in mJy

    i_890M = int(log(890.d6/nu(1))/log_nu_int) + 1        ! 1.4 GHz 
    i_1G = int(log(1.42d9/nu(1))/log_nu_int) + 1        ! 1.4 GHz
    i_3G = int(log(3.00d9/nu(1))/log_nu_int) + 1        ! 3.0 GHz
    i_5G = int(log(4.86d9/nu(1))/log_nu_int) + 1        ! 4.8 GHz
    flux0 = flux/nu(i_890M)
    flux1 = flux/nu(i_1G)   
    flux2 = flux/nu(i_3G) 
    flux3 = flux/nu(i_5G)

    !------------------------------------------------------------
    ! Calculate the luminosity from the hot layer in the FS
    ! L_ff = (1 - T_2/T_1) * Linj_perp

    if (radiative) then 

    open(19,file=''//shock//'_ff.dat')
    open(20,file='S_ff_'//shock//'.dat')

    do l=1,m_theta  
      integ = 0.d0
      do i_ph=1,mE_ph   
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
        ! Thermal emission from the post shock layer in a radiative shock
        ! L_ff = ( 1-(8000K/T(l)) ) * L_inj_perp
        emi_ff_nu(l,i_ph) = ( cte_l_1(l) * exp(-Eph(i_ph)/(k*T(l))) * gffinter )/h  ! erg/s/cm3/erg

        alpha_ff_nu1(l,i_ph) = 3.7d8 * (Zq_sw**2*(rho(l)/mp/mu_ISM)**2/sqrt(T(l)) ) * (1.d0-exp(-h*nu(i_ph)/(k*T(l))) )/nu(i_ph)**3 * gffinter
        Tau_1(l,i_ph) = alpha_ff_nu1(l,i_ph)*H_1(l)
        if ( Tau_1(l,i_ph)>= 1.d0 ) then
          write(*,96) 'Region 1 is optically thick at angle', theta(l)*180./pi, 'and frequency', nu(i_ph)/1.d6, 'MHz'
        end if 
      
        dEph = Eph(i_ph)*(Eph_int - 1.d0)
        integ = integ + emi_ff_nu(l,i_ph)*dEph
      end do

      if (T(l) >= T_ISM) then
        A(l) = L_inj(l)*(1-T_min(l)/T(l))/integ   ! cm3
      else if (T(l) <= T_ISM) then
        A(l) = 0.d0
      end if
    end do

    do i_ph=1,mE_ph
      L_th(i_ph) = 0.d0
      dEph = Eph(i_ph)*(Eph_int - 1.d0)
      do l=1,m_theta
        f_int(l) = A(l)*emi_ff_nu(l,i_ph)*4.d0*mphi    ! 1/s
        L_th(i_ph) = L_th(i_ph) + f_int(l)
        F_nu(l,i_ph) = f_int(l)*Eph(i_ph)              ! erg/s
        F_nu(l+m_theta,i_ph) = F_nu(l,i_ph)
      end do

      !call integrate(m_theta,dtheta,f_int(l),integ_HL)
      !L_th(i_ph) = integ_HL
      S(i_ph) = L_th(i_ph)*Eph(i_ph)*flux/nu(i_ph)  ! Total flux at a certain freq. in mJy
      L_tot = L_tot + L_th(i_ph)*dEph

      write(19,*) Eph(i_ph)/eV, L_th(i_ph)*Eph(i_ph)
      write(20,*) nu(i_ph), S(i_ph)
    end do

    close(19)
    close(20)

    print*, '\n f-f emission from post shock layer:'

    write(*,90) 'S_890MHz[mJy]=',L_th(i_890M)*Eph(i_890M)*flux0
    write(*,90) 'S_1.4GHz[mJy]=',L_th(i_1G)*Eph(i_1G)*flux1
    write(*,90) 'S_3.0GHz[mJy]=',L_th(i_3G)*Eph(i_3G)*flux2
    write(*,90) 'S_4.86GHz[mJy]=',L_th(i_5G)*Eph(i_5G)*flux3

    write(*,80) 'L_tot =',L_tot

    end if

    !------------------------------------------------------------
    !------------------------------------------------------------
    ! Calculate the luminosity using an adaptation of Eq. 24 from Christie+2016:
    ! L = integral(f_int(theta), with f_int(theta) = emi_ff_nu * dV_ring(theta),
    ! emi_ff_nu= dW/dt/dV/dnu the specific emissivity given by Eq. 5.14b from Rybicki,
    ! and dV_ring(l) the volume of the shocked ring element

    open(39,file='ionized_ff_'//shock//'.dat')
    open(40,file='S_ionized'//shock//'.dat')

    do i_ph=1,mE_ph
      L_ion(i_ph) = 0.d0
      do l = 1,m_theta
        ! Eq. 5.14b from Rybicki gives the spectrum in erg/s/cm3/Hz, so I divide by h
        ! emi_ff_nu = (6.8d-38*Zq_sw**2*ne*ni*T(l)**(-0.5d0)*exp(-Eph(i_ph)/(k*T(l)))*g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg

        ! Interpolate Gaunt factor
        if (shock == 'FS') then
            gamma2 = Zq_ISM**2.*Ry/k/T_f
        else if (shock == 'RS') then
            gamma2 = Zq_sw**2*Ry/k/T(l)
        end if

        j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
        dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

        if (shock == 'FS') then
            mu_ff = Eph(i_ph)/k/T_f
        else if (shock == 'RS') then
            mu_ff = Eph(i_ph)/k/T(l)
        end if

        j_mu = int( (log10(mu_ff)+16.d0)/0.2d0 ) + 1
        if (j_mu < nu_mga) then
          dmu = ((log10(mu_ff)+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
          call interpol2(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
        else
          dmu = ((log10(mu_ff)+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
          call extrapol2(dg,dmu,gff(j_g,nu_mga-1),gff(j_g+1,nu_mga-1),gff(j_g,nu_mga),gffinter)
        endif
        ! write(*,*) i_ph, l, gffinter, 9.77d0*(1.d0+0.13d0*log10(T(l)**1.5d0/Zq_sw/(Eph(i_ph)/h))) ! aprox 
        ! Emission from the ionized layer by the star's radiation field. T = 8.000 K (FS)
        ! Thermal emission from the RS
        if (shock == 'FS') then
            emi_ion_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T_f)) * gffinter )/h ! erg/s/cm3/erg
        else if (shock == 'RS') then
            emi_ion_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * gffinter )/h ! erg/s/cm3/erg
        end if

        if (shock == 'FS') then
          alpha_ff_nu2(l,i_ph) = 3.7d8 * (Zq_sw**2*n_a(l)**2/sqrt(T_f)) * (1.d0-exp(-h*nu(i_ph)/(k*T_f)) )/nu(i_ph)**3 * gffinter
        else if (shock == 'RS') then
          alpha_ff_nu2(l,i_ph) = 3.7d8 * (Zq_sw**2*n_a(l)**2/sqrt(T(l))) * (1.d0-exp(-h*nu(i_ph)/(k*T(l))) )/nu(i_ph)**3 * gffinter
        end if
        Tau_2(l,i_ph) = alpha_ff_nu2(l,i_ph)*H_theta(l)
        !if ( Tau_2(l,i_ph)>= 1.d0 ) then
        !  write(*,96) 'Region 2 is optically thick at angle', theta(l)*180./pi, 'and frequency', nu(i_ph)/1.d6, 'MHz'
        !end if

        f_vol_ion(l) = emi_ion_nu*dV_ring(l)
        F_nu_ion(l,i_ph) = f_vol_ion(l)*Eph(i_ph)    ! erg/s
        F_nu_ion(l+m_theta,i_ph) = F_nu_ion(l,i_ph)
        L_ion(i_ph) = L_ion(i_ph) + f_vol_ion(l)
      end do
      !call integrate(m_theta,dtheta,f_vol_ion,integ_ion)

      !L_ion(i_ph) = integ_ion
      S_ion(i_ph) = L_ion(i_ph)*Eph(i_ph)*flux/nu(i_ph)  ! Total flux at a certain freq. in mJy

      write(39,*) Eph(i_ph)/eV, Eph(i_ph)*L_ion(i_ph), L_ion(i_ph)
      write(40,*) nu(i_ph), S_ion(i_ph)

      dEph = Eph(i_ph)*(Eph_int - 1.d0)
      L_ion_tot = L_ion_tot + L_ion(i_ph)*dEph
      !print*, L_ion(i_ph), dEph, Eph_int, Eph(i_ph)
    end do

    close(39)
    close(40)

    if (shock == 'FS') then
      print*, '\n f-f emission from the ionized layer by the star:'
    else if (shock == 'RS') then
      print*, '\n f-f emission from the RS:'
    end if

    write(*,90) 'S_890MHz[mJy]=',L_ion(i_890M)*Eph(i_890M)*flux0
    write(*,90) 'S_1.4GHz[mJy]=',L_ion(i_1G)*Eph(i_1G)*flux1
    write(*,90) 'S_3.0GHz[mJy]=',L_ion(i_3G)*Eph(i_3G)*flux2
    write(*,90) 'S_4.86GHz[mJy]=',L_ion(i_5G)*Eph(i_5G)*flux3

    write(*,80) 'L_ion_tot =',L_ion_tot

    !------------------------------------------------------------
    !------------------------------------------------------------
    !------------------------------------------------------------
    ! Emission maps

    ! Hot layer
    if (radiative) then
      open(25, file='CD_hot_layer_1GHz.dat')
      open(26, file='CD_hot_layer_3GHz.dat')
      open(27, file='CD_hot_layer_5GHz.dat')

      do l = 1,m_l
        do i_ph = 1, mE_ph
          F_nu_phi(l,i_ph) = F_nu(l,i_ph)/mphi/2.d0   ! 4 quadrants, but already mirrored
          !if (l==2 .and. (nu(i_ph) >= nu_av_min)) print*, 'F_HOT =', F_nu(l,i_ph), nu(i_ph)/1.d9
        end do
        S1 = F_nu_phi(l,i_1G)*flux1
        S2 = F_nu_phi(l,i_3G)*flux2
        S3 = F_nu_phi(l,i_5G)*flux3
        write(25,*) S1
        write(26,*) S2
        write(27,*) S3
      end do

      close(25)
      close(26)
      close(27)
    end if

    !------------------------------------------------------------
    ! Emission from the RS and the isothermal region in the FS
    if (shock /= 'BOTH') then
      open(30,file='CD_th_'//shock//'_1GHz.dat') ! 1.4 GHz
      open(31,file='CD_th_'//shock//'_3GHz.dat') ! 3.0 GHz
      open(32,file='CD_th_'//shock//'_5GHz.dat') ! 4.8 GHz
    end if

    do l = 1,m_l
      do i_ph=1,mE_ph
        F_nu_phi_ion(l,i_ph) = F_nu_ion(l,i_ph)/mphi/2.d0   ! 4 quadrants
        !if (l==2 .and. (nu(i_ph) >= nu_av_min)) print*, 'F_ion =', F_nu_ion(l,i_ph), nu(i_ph)/1.d9
      end do
      S1 = F_nu_phi_ion(l,i_1G)*flux1
      S2 = F_nu_phi_ion(l,i_3G)*flux2
      S3 = F_nu_phi_ion(l,i_5G)*flux3
      write(30,*) S1
      write(31,*) S2
      write(32,*) S3
    end do

    close(30)
    close(31)
    close(32)

    !------------------------------------------------------------
    ! Avg flux btw nu_av_min and nu_av_max

    if (radiative) then
      open(50,file='CD_hot_radiomap_avg.dat')
    end if
    open(51,file='CD_th_radiomap_avg_'//shock//'.dat')

    ! Count the number of steps
    N_avg = 0
    do i_ph=1,mE_ph
      if ((nu(i_ph) >= nu_av_min) .and. (nu(i_ph) <= nu_av_max)) then
        N_avg = N_avg + 1
      end if
    end do

    h_avg = (nu_av_max/nu_av_min)**(1.d0/(N_avg-1))

    do l = 1,m_l
      if (radiative) then
        F_nu_tot(l) = 0.d0      ! hot layer
      end if
      F_nu_ion_tot(l) = 0.d0  ! isothermal layer and RS
      do i_ph=1,mE_ph
        if ((nu(i_ph) >= nu_av_min) .and. (nu(i_ph) <= nu_av_max)) then
          dnu_avg = nu(i_ph)*(h_avg-1.d0)
          if (radiative) then
            F_nu_tot(l) = F_nu_tot(l) + F_nu_phi(l,i_ph)*flux/nu(i_ph)*dnu_avg
          end if
          F_nu_ion_tot(l) = F_nu_ion_tot(l) + F_nu_phi_ion(l,i_ph)*flux/nu(i_ph)*dnu_avg
        end if  
      end do

      if (radiative) then
        F_nu_avg(l) = F_nu_tot(l)/(nu_av_max-nu_av_min)
        write(50,*) F_nu_avg(l)
      end if
      F_nu_ion_avg(l) = F_nu_ion_tot(l)/(nu_av_max-nu_av_min)
      write(51,*) F_nu_ion_avg(l)
    end do

    if (radiative) then
      close(50)
    end if
    close(51)

    !------------------------------------------------------------
    !------------------------------------------------------------
    ! Save the opacity at 1 GHz
    open(50, file='Tau.dat')
      do l = 1, m_theta
        write(50,*) theta(l), alpha_ff_nu1(l,i_1G), alpha_ff_nu2(l,i_1G)
      end do
    close(50)

    !------------------------------------------------------------  

    deallocate ( H_theta, H_1, rho, T, rho_f )
    deallocate ( F_nu,F_nu_phi,F_nu_ion,F_nu_phi_ion )
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

