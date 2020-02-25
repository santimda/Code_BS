module global

implicit none
public

!===================================================

integer, parameter :: dp = kind(1.d0)
integer, parameter :: nth = 2   ! number of threads (paralelizar)

real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: par = sqrt(2.d0/3.d0)

!Physical constants in CGS units:
real(dp), parameter :: eV = 1.602176d-12  ! electron-volt
real(dp), parameter :: c = 2.99792458d10  ! speed of light in vacuum
real(dp), parameter :: k = 1.380658d-16   ! boltzmann's constant 
real(dp), parameter :: qe = 4.8032068d-10 ! electron charge
real(dp), parameter :: mec2 = 0.511d6*eV  ! electron rest energy
real(dp), parameter :: re = 2.8d-13       ! electron classical radius
real(dp), parameter :: sigma = 5.6705d-5  ! stephan-boltzmann's constant
real(dp), parameter :: h = 6.6260755d-27  ! planck's constant
real(dp), parameter :: h_bar = h/(2.*pi)  ! 
real(dp), parameter :: alfa = 7.29735308d-3 ! fine structure constant
real(dp), parameter :: M_sun_yr = 1.989d33/31556926.d0 ! M_sun/yr
real(dp), parameter :: R_sun = 6.96d10    ! solar radius 
real(dp), parameter :: L_sun = 4.d33      ! solar luminosity (CGS)
real(dp), parameter :: AU = 1.496d13      ! Astronomical Unit
real(dp), parameter :: pc = 3.0856d18     ! parsec
real(dp), parameter :: kpc = 3.0856d21    ! kiloparsec
real(dp), parameter :: mp = 1.673d-24     ! proton mass
real(dp), parameter :: g_ad = 5.d0/3.d0   ! adiabatic coefficient 
real(dp), parameter :: m_pic2 = 134.98d6*eV ! pion rest energy
real(dp), parameter :: mpc2 = 1836.d0*mec2  ! proton rest energy

!===================================================


! System Parameters:
! ******************
!real(dp), parameter :: distance = 1.45d0*kpc       ! source distance Peri (?)
real(dp), parameter :: distance = 1.32d0*kpc        ! source distance K2017
real(dp), parameter :: inclination = 75.d0*pi/180.d0  ! angle between the star motion and the LOS in rad   
! i = 0° is BS in front of the star, i = 90° is edge-on, , i = 180° is BS behind the star
real(dp), parameter :: incli = pi/2.d0 - inclination   ! angle that I use
! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star
!real(dp), parameter :: R0_proj = 3.5d0*60.d3      ! projected BS distance to the star (3.5') in mas  
!real(dp), parameter :: R0 = R0_proj*(distance/kpc)*UA/abs(cos(incli)) ! linear distance to the BS stagnation point
real(dp), parameter :: theta_max = (137.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

!real(dp), parameter :: ab_X=0.705d0, ab_Y=0.275d0, ab_Z=0.02 ! stellar abundances 

! Star:
! *****

!real(dp), parameter :: v_st = 68.d5           ! Peri 2014
!real(dp), parameter :: v_st = 30.d5           ! K2018
real(dp), parameter :: v_st = 40.d5           ! Brookes 2016 da v_t=38.4 
!real(dp), parameter :: T_st = 4.25d4           ! Benaglia 2010  
real(dp), parameter :: T_st = 4.07d4           ! K2017  
real(dp), parameter :: T_st_ad = k*T_st/mec2   ! Adimensional stellar temperature
!real(dp), parameter :: L_st = 7.d5*L_sun       ! del Valle 2012
!real(dp), parameter :: L_st = 8.5d5*L_sun       ! Wikipedia
!real(dp), parameter :: R_st = sqrt( L_st/(4.*pi*sigma*T_st**4.) ) 
real(dp), parameter :: R_st = 19.*R_sun         ! K2017  
!real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4. ! Stellar luminosity 
real(dp), parameter :: L_st = 8.7d5*L_sun      ! K2018
!real(dp), parameter :: Mdot = 1.d-5*M_sun_yr  ! del Valle 2012
!real(dp), parameter :: Mdot = 6.5d-6*M_sun_yr  ! Peri 2014
real(dp), parameter :: Mdot = 9.d-6*M_sun_yr  ! 
!real(dp), parameter :: Mdot = 5.5d-6*M_sun_yr  ! Brookes 2016
!real(dp), parameter :: Mdot = 1.3d-6*M_sun_yr  ! K2018, 
!scales with: *(v_st/30.d5)**2. * (n_ISM/19)
!real(dp), parameter :: v_w = 2200.d5          ! del Valle 2012
real(dp), parameter :: v_w = 2300.d5          ! Brookes2016,Ben2010
!real(dp), parameter :: v_w = 2325.d5          ! Peri 2014
!real(dp), parameter :: v_w = 3000.d5          ! K2018

real(dp), parameter :: v_st_rot = 0.1d0*v_w   ! Rule-of-thumb
real(dp), parameter :: mu_w = 1.3_dp          ! Leitherer 1995 
real(dp), parameter :: Zq_w = 1.0_dp          ! Z = rms ionic charge
real(dp), parameter :: gamma_w = 1.0_dp       ! Leitherer 1995
real(dp), parameter :: filling = 1.d-1        ! ..... 10% ?

! ISM:
! ****
real(dp), parameter :: v_ISM = v_st           ! 
!real(dp), parameter :: n_ISM = 1.d0           ! del Valle 2012
!real(dp), parameter :: n_ISM = 0.2d0           ! Peri 2014
!real(dp), parameter :: n_ISM = 3.5d0           ! Peri 2015
!real(dp), parameter :: n_ISM = 13.9d0           ! Brookes2016
!real(dp), parameter :: n_ISM = 19.d0           ! K2018
real(dp), parameter :: n_ISM = 15.d0           ! 
real(dp), parameter :: T_ISM = 8000.d0          ! warm Brookes p.136
real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2! Adimensional ISM temperature 
!real(dp), parameter :: U_ISM = 4.d0*sigma*T_ISM**4/c ! average energy density
real(dp), parameter :: alpha_ISM = 4.d-2        ! alpha = P_thermal_ISM/P_kin_ISM
real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)

real(dp), parameter :: mu_ISM = 1.37_dp        ! K2018 
real(dp), parameter :: Zq_ISM = 1.0_dp         ! Z = rms ionic charge
real(dp), parameter :: gamma_ISM = 1.0_dp      ! ?

real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2.) ) 
real(dp), parameter :: R0_proj = (R0/AU)/(distance/pc)*abs(cos(incli)) ! angular distance ["] to the BS stagnation point
! Peri 2014: R0_proj ~ 3.5 arcmin
! K2018: R0_proj ~ 193" ~ 3.2 arcmin

! Shocked stellar-wind (sw):
! **************************
real(dp), parameter :: mu_sw = 0.615_dp    ! mean atomic weight = 4/(6x+y+2) for complete ionization 
real(dp), parameter :: mu_e_sw = 1.173_dp  ! mu_e = 2/(1+x)
real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
real(dp), parameter :: Zq_sw = 1.23_dp        
real(dp), parameter :: gamma_sw = 1.34_dp

! Shocked ISM (sism):
! *******************    -----------_> ????????????????
real(dp), parameter :: mu_sism = 0.615_dp  
real(dp), parameter :: mu_e_sism = 1.173_dp  ! mu_e = 2/(1+x)
real(dp), parameter :: mu_i_sism = 1.d0/(1/mu_sism-1/mu_e_sism)   
real(dp), parameter :: Zq_sism = 1.23_dp        
real(dp), parameter :: gamma_sism = 1.34_dp

! Heated dust (IR):
! *****************    -----------_> ????????????????
!real(dp), parameter :: T_IR = 27.*(0.2)**(-1./6.)*(L_st/1.d38)**(1./6.)*(R0/pc)**(-1./3.) ! Eq.10 del Valle 2012
!real(dp), parameter :: T_IR_ad = k*T_IR/mec2  ! Adimensional dust temperature
!real(dp), parameter :: L_IR = 4.d0*pi*(0.2*R0)**2*sigma*T_IR**4. ! Dust luminosity ?????
real(dp), parameter :: L_IR = L_st/350.d0 ! Dust luminosity, K2017
real(dp), dimension(:), allocatable :: T_IR,T_IR_ad


! Non-thermal particles:
! **********************
real(dp), dimension(:), allocatable :: efi_ac_sw         ! particle accel. efficiency at the sw  (l)
!real(dp), parameter :: efi_ac_sism = 2.d0*pi*(c/v_ism)**2 ! particle accel. efficiency at the sism

real(dp), parameter :: frac_sw = sqrt(0.999_dp) ! fraction sqrt(UB_sw/UB_equipartition_sw)
!real(dp), parameter :: frac_sism = 0.5_dp   ! fraction B_sism/B_equipartition_sism
real(dp), parameter :: alpha_sw = 2.2_dp    ! injection spectral index in the sw 
!real(dp), parameter :: alpha_sism = 2.0_dp  ! injection spectral index in the sism

!real(dp), parameter :: iny_e = frac_NT*0.015d0       ! Fraction L_NT_e/L_NT
real(dp), parameter :: iny_e = 0.0042d0       ! Fraction L_NT_e/L_inj_perp
!real(dp), parameter :: iny_p = frac_NT*(1.d0-iny_e) ! Fraction L_NT_p/L_iny_perp 
real(dp), parameter :: iny_p = 0.5d0                 ! Fraction L_NT_p/L_iny_perp
real(dp), parameter :: frac_NT = iny_e+iny_p !0.5_dp              ! Fraction L_NT/L_inj_perp 

!===================================================


! grid/vectors size
integer, parameter :: m_theta = 400 ! angular resolution in the CD/shock
integer, parameter :: mE_e = 256    ! energy resolution for electrons
integer, parameter :: mE_p = 200    ! energy resolution for protons < mE_f 
integer, parameter :: mphi = 4     ! amount of linear-emitters
real(dp) :: delta_phi = 0.5*pi/(mphi-1.d0) ! angular separation between linear-emitters
real(dp) :: d_phi = 0.5*pi/mphi ! angular width of each linear-emitter or wedge

real(dp), parameter :: f_Ny = 0.8d0                 ! Nyquist frequency for convolved maps

!===================================================


integer :: i,j,l

! CD
real(dp), dimension(:), allocatable :: x,y,z
real(dp), dimension(:), allocatable :: theta,r_norm,R_theta,H_norm,H_theta

! termo
real(dp), dimension(:), allocatable :: P,rho,T,v_t,v_perp,B,B_par

! dist
real(dp), parameter :: Emin_e = 2.d0*mec2   ! ~1 MeV
real(dp), parameter :: Emax_e = 1.d14*eV    ! 100 TeV
real(dp) :: Eint_e,dE_e
real(dp), dimension(:,:), allocatable :: N_etot ! Ne(theta,E_i)
real(dp), dimension(:), allocatable :: E_e,N_e  ! N_e(E_i) is for a fixed "l"

real(dp), parameter :: Emin_p = 1.d9*eV ! 1 GeV 
!real(dp), parameter :: Emin_p = 2.d0*mpc2   ! ~ 2 GeV
real(dp), parameter :: Emax_p = 5.d15*eV    ! 5 PeV
real(dp) :: Eint_p,dE_p
real(dp), dimension(:,:), allocatable :: N_ptot ! Np(theta,E_i)
real(dp), dimension(:), allocatable :: E_p,N_p  ! N_p(E_i) is for a fixed "l" 

!===================================================
! emission

real(dp) :: P_Eps, dP_Eps, P_Ef, dP_Ef

! syn
integer, parameter :: m_nu = 600  ! frequency resolution for synchrotron photons
real(dp), dimension(:,:), allocatable :: Lsy      ! L_sy(l,E_i)
real(dp), dimension(:), allocatable :: nu,Lsy_tot,Lsy_tot_abs ! Lsy_tot(i) = sum_l L_sy(l,i)

! IC
real(dp), dimension(:), allocatable :: LIC,LIC_tot,LIC_tot_abs
real(dp) :: Kappa_st, Kappa_IR(m_theta) ! = (R/2D)**2. 
!real(dp), parameter :: Kappa_IR = 1.d0 ! Is it ok ????????
real(dp) :: X_st(3),X_emi(3),X_eE(3),X_obs(3) ! Position vectors for absorption integrals
real(dp), parameter :: cte_Ndot = 2.*re**2.*mec2**3./(pi*h_bar**3.*c**2.)

! br
integer, parameter :: mE_f = 300 ! it is the same for p-p and IC
real(dp), dimension(:), allocatable :: Lbr,Lbr_tot,Lbr_tot_abs

! pp
real(dp), dimension(:), allocatable :: Lpp,Lpp_tot,Lpp_tot_abs

!===================================================
! absorption

! g-g
integer, parameter :: mE_gg = 600, m_ang_gg = 400
real(dp), dimension(:,:), allocatable :: tau
real(dp), parameter :: Ef_min_gg = 5.d7*eV  ! 50 MeV 
real(dp), parameter :: Ef_max_gg = 5.d13*eV ! 50 TeV 
real(dp) :: Ef_gg(mE_gg),Ef_int_gg

! f-f
integer, parameter :: mE_ff = m_nu, mr_ff = 1640, mr_map = 1000
real(dp), dimension(:,:), allocatable :: tauff
real(dp) :: alphaff_w(mr_ff,m_nu), r_w(mr_ff)
!real(dp) :: Lsy_sw(m_theta,m_nu),Lsy_sism(m_theta,m_nu),Lsy_tot_abs(m_nu) 

!===================================================

contains

!============================================================================
! vectorE
!============================================================================
  
  SUBROUTINE vector_log(Emin,Emax,Eint,E_vec)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: Emin,Emax
    integer :: N_puntos, i_vec
    real(dp), intent(out) :: Eint
    real(dp), intent(out) :: E_vec(:)

    N_puntos = SIZE(E_vec)
    Eint = (Emax/Emin)**(1.d0/dfloat(N_puntos-1))

 !   write(*,*) Emin, Emax, Eint,N_puntos

    E_vec(1) = Emin
    do i_vec=2,N_puntos
       E_vec(i_vec)=E_vec(i_vec-1)*Eint    
    enddo

    RETURN

  END SUBROUTINE vector_log


!============================================================================
! integrate
!============================================================================
  
  SUBROUTINE integrate(N,H,FI,S)
    ! Subroutine for integration over f(x) with the Simpson rule.
    ! FI: integrand f(x); H: interval; S: integral.
    integer, intent(in) :: N
    real*8, intent(in) :: H,FI(N)
    real*8, intent(out) :: S
    real*8 :: S0,S1,S2

    S  = 0.d0
    S0 = 0.d0
    S1 = 0.d0
    S2 = 0.d0
!    DO 100 I = 2, N-1, 2
    DO I = 2, N-1, 2
       S1 = S1 + FI(I-1)
       S0 = S0 + FI(I)
       S2 = S2 + FI(I+1)
!100    CONTINUE
    END DO
    S = H*(S1 + 4.*S0 + S2)/3.
    ! If N is even, add the last slice separately
    IF(MOD(N,2).eq.0) S = S + H*(5.d0*FI(N) + 8.d0*FI(N-1) - FI(N-2))/12.d0
    RETURN
    
  End SUBROUTINE integrate

!============================================================================
!  Cooling_IC
    
  Function F_ani(u)  ! Used for anisotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_ani
    real(dp), intent(in):: u
    real(dp), parameter:: c_ani = 6.13

    F_ani = (c_ani*u*log(1.+2.16*u/c_ani)) / (1.+c_ani*u/0.822)

    RETURN
  end Function F_ani

  Function F_iso(u)  ! Used for isotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_iso,g_u
    real(dp), intent(in):: u
    real(dp), parameter:: c_iso = 5.68
    
    g_u = 1./(1. + (-0.362*u**0.682)/(1. + 0.826*u**1.281))
    F_iso = (c_iso*u*log(1.+0.722*u/c_iso)) / (1.+c_iso*u/0.822)*g_u ! Corregido
    
    RETURN
  end Function F_iso

  !   N^dot_IC
    
  Function F_1(x0)  ! Used for anisotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_1,g
    real(dp), intent(in):: x0
    g = 1./(1. + (0.153*x0**0.857)/(1. + 0.254*x0**1.84))
    F_1 = ( (pi**2./6. + x0)*exp(-x0) )*g   
    RETURN
  end Function F_1

  Function F_2(x0)  ! Used for anisotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_2,g
    real(dp), intent(in):: x0
    g = 1./(1. + (1.33*x0**0.691)/(1. + 0.534*x0**1.668))
    F_2 = ( (pi**2./6. + x0)*exp(-x0) )*g   
    RETURN
  end Function F_2

  Function F_3(x0)  ! Used for isotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_3,g
    real(dp), intent(in):: x0
    real(dp), parameter:: c_3 = 0.319
    g = 1./(1. + (0.443*x0**0.606)/(1. + 0.54*x0**1.481))
    F_3 = (pi**2./6.*(1.+c_3*x0)/(1.+pi**2.*c_3*x0/6.)*exp(-x0) )*g   
    RETURN
  end Function F_3

  Function F_4(x0)  ! Used for isotropic IC (Khangulyan 2014)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: F_4,g
    real(dp), intent(in):: x0
    real(dp), parameter:: c_4 = 6.62
    g = 1./(1. + (0.726*x0**0.461)/(1. + 0.382*x0**1.457))
    F_4 = (pi**2./6.*(1.+c_4*x0)/(1.+pi**2.*c_4*x0/6.)*exp(-x0) )*g   
    RETURN
  end Function F_4
  
 !============================================================================
     
  Function nBB(Epho,T)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: nBB,Epho,T,pi,k,c,h
    pi = acos(-1.d0)
    k = 1.380658d-16
    c = 2.99792458d10
    h = 6.6260755d-27
    nBB = 2.d0*pi/((h*c)**3.)*(Epho**2./(exp(Epho/(k*T))-1.))
    RETURN
  end Function nBB

  !============================================================================
  
  
  SUBROUTINE versorresta(V1,V2,N,versor)
    ! subroutine to calculate the versor: (V1-V2)/|V1-V2|
    integer, intent(in) :: N
    real(dp), intent(in) :: V1,V2
    real(dp), intent(out) :: versor
    real(dp) :: resta,suma
    DIMENSION V1(N),V2(N),versor(N),resta(N)
    
    suma=0.d0
    do i=1,N
       resta(i)=V1(i)-V2(i)
       suma=suma+resta(i)**2.
    enddo
    suma=sqrt(suma)          ! suma = || versor_resta ||
    do j=1,N
       versor(j)=resta(j)/suma
    enddo
    
    RETURN
  END SUBROUTINE versorresta
  
  
  !============================================================================
  


END module global
