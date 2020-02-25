module global

implicit none
public

!===================================================

integer, parameter :: dp = kind(1.d0)
integer, parameter :: nth = 4   ! number of threads (paralelizar)

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
real(dp), parameter :: T_CMB = 2.725d0             
real(dp), parameter :: T_CMB_ad = k*T_CMB/mec2
real(dp), parameter :: u_CMB = 0.26*eV
real(dp), parameter :: Ry = 2.17987d-11  ! Rydberg unit = 2.18e-18 J

!===================================================

! System Parameters:

real(dp), parameter :: distance = 0.155d0*kpc   ! source distance 

! Star:

real(dp), parameter :: v_st = 72.d5           ! System spatial velocity 
real(dp), parameter :: T_st = 26.5d4          ! WD temp
real(dp), parameter :: T_st_ad = k*T_st/mec2  ! Adimensional stellar temperature
real(dp), parameter :: R_st = 0.011d0*R_sun   ! WD radius
real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4. ! Stellar luminosity

!Wind:

real(dp), parameter :: Mdot = 1.0d-10*M_sun_yr  !
real(dp), parameter :: v_w = 1800.d5         	! 
real(dp), parameter :: L_wind = 0.5*Mdot*v_w**2 ! Wind kinetic power
real(dp), parameter :: mu_w = 1.3_dp          ! Leitherer 1995 , CHECK
real(dp), parameter :: Zq_w = 1.0_dp          ! Z = rms ionic charge , CHECK
real(dp), parameter :: gamma_w = 1.0_dp       ! Leitherer 1995 , CHECK
real(dp), parameter :: filling = 1.d-1        ! ..... 10% ?
real(dp), parameter :: R0_proj = 4.5d0/60.         ! Bubble radius in arcmin
real(dp), parameter :: R0 = (R0_proj*(2.*pi/(60.*360.)))*distance ! R0 in cm
real(dp), parameter :: H0 = R0*0.3            ! shock width
real(dp), parameter :: Vol = R0**3 ! Bow shock volume

! Shocked stellar-wind (sw):
! **************************
real(dp), parameter :: mu_sw = 0.615_dp    ! mean atomic weight = 4/(6x+y+2) for complete ionization, CHECK
real(dp), parameter :: mu_e_sw = 1.173_dp  ! mu_e = 2/(1+x), CHECK
real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
real(dp), parameter :: Zq_sw = 1.23_dp        
real(dp), parameter :: gamma_sw = 1.34_dp
real(dp), parameter :: rho_w = Mdot/(4.d0*pi*R0**2*v_w)
real(dp), parameter :: n_a = 4.*rho_w/(mp*mu_sw) ! ambient nuclei number density
real(dp), parameter :: T_6 = 10.d0         ! temperature in MK ?????

! For IC (both diste.f90 and radic.f90)
real(dp), parameter :: Kappa_st = (R_st/(2.*R0))**2
real(dp), parameter :: Kappa_CMB = 1.d0             


! Non-thermal particles:

real(dp), parameter :: R_acel = H0  
real(dp), parameter :: efi_ac_sw = 2.d0*pi*(c/v_w)**2 ! particle accel. efficiency at the sw 
real(dp), parameter :: alpha_sw = 2.01_dp    ! injection spectral index in the sw 

real(dp), parameter :: Kep = 1.d-2          ! Kep = U_e/U_p ~ 0.01 - 1
real(dp), parameter :: eta_e = Kep/(1.+Kep) ! U_e = eta_e*U_NT
real(dp), parameter :: eta_p = 1.d0-eta_e   ! U_p = eta_p*U_NT
real(dp), parameter :: eta_mag = 0.75d0     ! U_NT = eta_mag*U_mag

real(dp), parameter :: S_obs = 0.03d-26 	! MeerKAT upper-limit in mJy
real(dp), parameter :: nu_obs = 1.284d9     ! Observing frequency

real(dp), parameter :: a = 0.103d0           ! a(p=2)=0.103,a(2.2)=0.0933,a(2.5)=0.085,a(3.0)=0.074
real(dp), parameter :: div = (a*qe**3/h/mec2)*(3.*qe*c/(4.*pi*mec2**3*nu_obs))**((alpha_sw-1)/2.)
real(dp), parameter :: Emin_e = 2.d0*mec2   
real(dp), parameter :: Ecooled_e = 0.5d12*eV 
real(dp), parameter :: Emax_e_aprox = 1.d12*eV
real(dp), parameter :: alpha2 = alpha_sw - 2.d0
real(dp), parameter :: int_aprox1 = ( Emin_e**(-alpha2)-Ecooled_e**(-alpha2) )/alpha2
real(dp), parameter :: int_aprox2 = ( Ecooled_e**(-alpha2-1.)-Emax_e_aprox**(-alpha2-1.) )/(alpha2-1.)
real(dp), parameter :: int_aprox = int_aprox1+int_aprox2
!There is an inconsistency in using the expression below (for constant spectral index) as well as the 
!broken power-law distribution 
real(dp), parameter :: B = ((S_obs*8.*pi*distance**2*int_aprox)/(h*Vol*div*eta_e*eta_mag))**(2./(alpha_sw+5.))
real(dp), parameter :: B_par = B*sqrt(2./3.)

!===================================================

! grid/vectors size
integer, parameter :: mE_e = 300    ! energy resolution for electrons
integer, parameter :: mE_f = 512    ! energy resolution for br and IC > mE_p
integer, parameter :: mE_p = 300    ! energy resolution for electrons
integer, parameter :: mE_ff = 300    ! energy resolution for electrons

!===================================================

integer :: i,j,l

! dist
real(dp) :: Eint_e,dE_e
real(dp), dimension(:), allocatable :: E_e,N_e 

real(dp), parameter :: Emin_p = 1.d9*eV     ! 1 GeV
real(dp) :: Eint_p,dE_p
real(dp), dimension(:), allocatable :: E_p,N_p 


!===================================================
! emission

real(dp) :: P_Eps, dP_Eps

! syn
integer, parameter :: m_nu = 800   ! frequency resolution
real(dp), dimension(:), allocatable :: nu, Lsy, Lsy_abs

! br
real(dp), dimension(:), allocatable :: Lbr

! IC
real(dp) :: Lic(mE_f)
!real(dp), dimension(:), allocatable :: Lic

! p-p
real(dp), dimension(:), allocatable :: Lpp

real(dp), parameter :: nu_min_ff = 1.d6    ! 1 MHz 
real(dp), parameter :: nu_max_ff = 1.d13   ! 10 THz
real(dp) :: nu_int_ff

!===================================================

contains

!============================================================================
! vector_log
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
! integrar
!============================================================================
  
  SUBROUTINE integrar(N,H,FI,S)
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
       
  End SUBROUTINE integrar

 !============================================================================
 !  Black-Body spectrum

  Function nBB(Epho,T)
    integer, parameter :: dp = kind(1.d0)
    real(dp) :: nBB,Epho,T,pi,k,c,h
    pi = acos(-1.d0)
    k = 1.380658d-16
    c = 2.99792458d10
    h = 6.6260755d-27
    nBB = 2.d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(k*T))-1))
    RETURN
  end Function nBB
    
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
     


END module global
