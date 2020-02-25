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

    allocate( theta(m_theta),r_norm(m_theta),R_theta(m_theta),H_norm(m_theta),H_theta(m_theta) )
    allocate( rho(m_theta),T(m_theta) )

80  format(A12,X,Es10.2)


     write(*,*) 'Calculate the total thermal emission from the shocked wind'
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
    integer, parameter :: nu_mu=146,nu_mga=81 ! # data in table gauntff.dat
    real(dp) :: Ry,gamma2,dg,dmu,mu_ff,gff(nu_mu,nu_mga),gffinter
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
    
    Ry = 2.17987d-11          ! Rydberg unit = 2.18e-18 J    

    allocate( theta(m_theta),r_norm(m_theta),R_theta(m_theta),H_norm(m_theta),H_theta(m_theta) )
    allocate( rho(m_theta),T(m_theta) )
   
    write(*,*) 'Calculate the spectrum of the thermal emission from the shocked wind'

    LX_sw = 0.d0
    L_tot = 0.d0

80  format(A12,X,Es10.2)


!///////////////////////////////////////////////////////////////////////////
!
!     We extract the Gaunt factor from the table given by van Hoof (2014):
!     <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, with step 1
!   
      open(13,file='gauntff.dat')
      do i=1,16
         read(13,*)
      enddo
      
      do i=1,nu_mu
         read(13,*) (gff(i,j), j=1,nu_mga)
      enddo
      close(13)
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
    Eph_min = 1.d-2*eV !Eph_0/1.d3
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
          
!          write(*,*) i_ph,l,j_g,dg,j_mu,dmu,gffinter
!          write(*,*) i_ph, l, gffinter, 9.77d0*(1.d0+0.13d0*log10(T(l)**1.5d0/Zq_sw/(Eph(i_ph)/h))) ! aprox 
          emi_ff_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * gffinter )/h ! erg/s/cm3/erg

!          emi_ff_nu = ( cte_l(l) * exp(-Eph(i_ph)/(k*T(l))) * g_ff(Eph(i_ph),T(l)) )/h ! erg/s/cm3/erg
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

  !////////////////////////////////////////////////////////////////////////////////
  !

  Subroutine interpol2(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
    implicit none 
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
    real(dp), intent(out) :: finterpol
    real(dp) :: fx,fy
    fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
    fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
    finterpol = fx1y1 + fx*dx + fy*dy
    return
  END Subroutine interpol2
  

  !////////////////////////////////////////////////////////////////////////////////
  !

  Subroutine extrapol2(dx,dy,fx1y1,fx2y1,fx1y2,fextrapol)
    implicit none 
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
    real(dp), intent(out) :: fextrapol
    real(dp) :: fx,fy
    fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
    fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
    fextrapol = fx1y2 + fx*dx + fy*dy
    return
  END Subroutine extrapol2
  
