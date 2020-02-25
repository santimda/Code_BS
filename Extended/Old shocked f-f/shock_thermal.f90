module shock_thermal
  ! This module calculates the shock's thermal bremsstrahlung emission
  use global
  
  implicit none
  private
  public :: shock_thermal_run_total, shock_thermal_run_spectrum
  
contains
  
  !============================================================================
  ! shock_thermal_run_total (START)
  !============================================================================
  
  subroutine shock_thermal_run_total  
    use global
    implicit none
    real(dp) :: dummy,ne,ni,n_a(m_theta),emi_ff
    real(dp) :: f_R,f_int(m_theta),R_max,dtheta,integral, L_tot
    integer :: nl
    real(dp), parameter :: g_average = 1.2 ! Rybicki, Eq. 5.15 b
    allocate ( theta(m_theta),r_norm(m_theta),R_theta(m_theta),H_norm(m_theta),H_theta(m_theta) )
    allocate ( rho(m_theta),T(m_theta) )

80  format(A12,X,Es10.2)

    !------------------------------------------------------------
    !     
    !    Read shock position and thermodynamical quantities   
    !
    !------------------------------------------------------------

    open(12,file='CD.dat')
    open(13,file='Hnorm_sw.dat') 
    do l=1,m_theta
       read(12,*) theta(l),r_norm(l)
       read(13,*) dummy,H_norm(l)
       ! R(l) = distance to the star, H(l) = width of the shocked region
       R_theta(l) = r_norm(l)*R0
       H_theta(l) = H_norm(l)*R_theta(l)
    enddo
    close(12)
    close(13)
    
    dtheta = (theta_max-theta(1))/float(m_theta)

    open(11,file='termo_sw.dat')
    do l=1,m_theta
       read(11,*) nl,rho(l),dummy,dummy,T(l),dummy ! CGS
       n_a(l) = rho(l)/mp/mu_i_sw   ! da del orden de 0.2 (lejos) - 1 (apex)
    enddo
    close(11)

    !------------------------------------------------------------
    !    
    !    Calculate the shocked wind (sw) bremsstrahlung total power
    !
    !------------------------------------------------------------
    
    ! Calculate the luminosity using Eq. 24 from Christie+2016:
    ! L = 2*pi*integral(f1_theta), with: f1_theta = sin(theta)*f_R, with:
    ! f_R = integral(f2_theta), with: f2_theta = R'(theta)**2.*emi_ff(theta), with:
    ! emi_ff= dW/dt/dV the total emissivity given by Eq. 5.15b from Rybicki
    ! Then: f_R = emi_ff*(Rmax**3.-R**3.)/3.
    
    do l=1,m_theta
       ni = n_a(l) 
       ne = gamma_sw*n_a(l)   
       emi_ff = 1.4d-27*Zq_sw**2.*ne*ni*sqrt(T(l))*g_average  
       R_max = R_theta(l) + H_theta(l)
       f_R = emi_ff*(R_max**3.-R_theta(l)**3.)/3.
       f_int(l) = sin(theta(l))*f_R       
    enddo
    
    call integrate(m_theta,dtheta,f_int,integral)
    L_tot = 2.*pi*integral
    write(*,80) 'L_tot =',L_tot
      
    deallocate ( theta,r_norm,R_theta,H_norm,H_theta )
    deallocate ( rho,T )

 
  END subroutine shock_thermal_run_total

  !============================================================================
  ! shock_thermal_run_total (END)
  !============================================================================
  

  !============================================================================
  ! shock_thermal_run_spectrum (START)
  !============================================================================
  
  subroutine shock_thermal_run_spectrum  
    use global
    implicit none
    integer, parameter :: mE_ph = 256
    real(dp) :: dummy,n_a(m_theta),emi_ff_nu, cte_l(m_theta)
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph), dEph
    real(dp) :: f_R,f_int(m_theta),R_max,dtheta,integral
    real(dp) :: LX_sw(mE_ph), L_tot
    integer :: nl, i_ph

    interface
       function g_ff(E_photon,T_emitter)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: g_ff
         real(dp), intent(in) :: E_photon,T_emitter
       end function g_ff
    end interface
    
    allocate ( theta(m_theta),r_norm(m_theta),R_theta(m_theta),H_norm(m_theta),H_theta(m_theta) )
    allocate ( rho(m_theta),T(m_theta) )
   
    LX_sw = 0.d0
    L_tot = 0.d0

80  format(A12,X,Es10.2)

    !------------------------------------------------------------
    !     
    !    Read shock position and thermodynamical quantities   
    !
    !------------------------------------------------------------

    open(12,file='CD.dat')
    open(13,file='Hnorm_sw.dat') 
    do l=1,m_theta
       read(12,*) theta(l),r_norm(l)
       read(13,*) dummy,H_norm(l)
       ! R(l) = distance to the star, H(l) = width of the shocked region
       R_theta(l) = r_norm(l)*R0
       H_theta(l) = H_norm(l)*R_theta(l)
    enddo
    close(12)
    close(13)

    dtheta = (theta_max-theta(1))/float(m_theta)

    open(11,file='termo_sw.dat')
    do l=1,m_theta
       read(11,*) nl,rho(l),dummy,dummy,T(l),dummy ! CGS
       n_a(l) = rho(l)/mp/mu_i_sw   ! da del orden de 0.2 (lejos) - 1 (apex)
!          ni = n_a(l) 
!          ne = gamma_sw*n_a(l)    
       cte_l(l) = 6.8d-38*Zq_sw**2.*gamma_sw*n_a(l)**2.*T(l)**(-0.5d0)
    enddo
    close(11)

    !------------------------------------------------------------
    !    
    !    Calculate the shocked wind (sw) bremsstrahlung spectrum
    !
    !------------------------------------------------------------

    Eph_0 = mp*v_w**2./32.    ! Use the k*T_sw(apex) as a reference
    Eph_min = Eph_0/8.d2
    Eph_max = Eph_0*4.d1
    
    call vector_log(Eph_min,Eph_max,Eph_int,Eph)
    
    ! Calculate the luminosity using Eq. 24 from Christie+2016:
    ! L = 2*pi*integral(f1_theta), with: f1_theta = sin(theta)*f_R, with:
    ! f_R = integral(f2_theta), with: f2_theta = R'(theta)**2.*emi_ff(theta), with:
    ! emi_ff_nu= dW/dt/dV/dnu the specific emissivity given by Eq. 5.14b from Rybicki
    ! Then: f_R = emi_ff_nu*(Rmax**3.-R**3.)/3.

    open(19,file='sw_ff.dat')
    do i_ph=1,mE_ph
 
       do l=1,m_theta     
          ! Eq. 5.14b from Rybicki gives the spectrum in erg/s/cm3/Hz, so I divide by h
!          emi_ff_nu = (6.8d-38*Zq_sw**2*ne*ni*T(l)**(-0.5d0)*exp(-Eph(i_ph)/(k*T(l)))*g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg
          emi_ff_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg
          R_max = R_theta(l) + H_theta(l)
          f_R = emi_ff_nu*(R_max**3.-R_theta(l)**3.)/3.
          f_int(l) = sin(theta(l))*f_R 
       enddo
    
       call integrate(m_theta,dtheta,f_int,integral)
       LX_sw(i_ph) = 2.*pi*integral
       write(19,*) Eph(i_ph)/eV,Eph(i_ph)*LX_sw(i_ph)

       dEph = Eph(i_ph)*(Eph_int - 1.d0)
       L_tot = L_tot + LX_sw(i_ph)*dEph
    enddo
    close(19)

    write(*,80) 'L_tot =',L_tot

    !--------------------------------
    !
    !  Save the results 
    !
    !--------------------------------

!    open(19,file='WCR_ff.dat')
!    do i=1, mE_ph2
!       write(19,*) Eph(i)/eV,L1_X(i),L2_X(i),L_X(i)
!    enddo
!    close(19)
    
    !---------------------------------
      
    deallocate ( theta,r_norm,R_theta,H_norm,H_theta )
    deallocate ( rho,T )

 
  END subroutine shock_thermal_run_spectrum

  !============================================================================
  ! shock_thermal_run_spectrum (END)
  !============================================================================
  
END module shock_thermal


 
   !-----------------------------------------------------------------
  
  Function g_ff(Eph_X,T_X)
    use global
    real(dp) :: g_ff,Eph_X,T_X,zeda

    ! See Fig. 5.2 from Rybicki, these ranges are not correct. 

    if (Eph_X/(k*T_X) >= 1.d0 .and. Eph_X/(k*T_X) < 1.d2) then
       g_ff = sqrt( (3.d0/pi) * ((k*T_X)/Eph_X) )
    else if (Eph_X/(k*T_X) < 1.d0 .and. Eph_X/(k*T_X) >= 1.d-2 ) then
       zeda = 1.781d0  ! Es el gamma en el astrophysical formulae
       g_ff = (sqrt(3.d0)/pi) * log(4.d0*k*T_X/(zeda*Eph_X))
    else if (Eph_X/(k*T_X) > 1.d2 .or. Eph_X/(k*T_X) < 1.d-2) then
       g_ff = 0.54*log(5.d7*(T_X**1.5 / (Eph_X/h)))
!       g_ff = 0.d0  ! FRUTAAAAAAAAAA     
    endif
    
!    write(*,*) g_ff, 0.54*log(5.d7*(T_X**1.5 / (Eph_X/h)))
!    g_ff = 0.54*log(5.d7*(T_X**1.5 / (Eph_X/h)))

    RETURN
  end Function g_ff
  
  !-----------------------------------------------------------------
  
