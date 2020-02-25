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
real(dp), parameter :: T_CMB = 2.73d0             
real(dp), parameter :: T_CMB_ad = k*T_CMB/mec2 ! h*(10.)**(34./3.)/mec2 (Ani)
real(dp), parameter :: u_CMB = 10.d0**(-236.d0/19.d0) ! Según Ani

!===================================================

! System Parameters:
! ******************
real(dp), parameter :: distance = 1.d0*kpc       ! source distance K2017
real(dp), parameter :: inclination = 90.d0*pi/180.d0     ! angle between the star motion and the LOS in rad   
! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
real(dp), parameter :: incli = pi/2.d0 - inclination   ! angle that I use
! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star

!real(dp), parameter :: ab_X=0.705d0, ab_Y=0.275d0, ab_Z=0.02 ! stellar abundances 

! Star:
! *****
real(dp), parameter :: v_st = 30.d5           ! Generic velocity for runaway star 
real(dp), parameter :: T_st = 4.d4            ! Generic T_star 
real(dp), parameter :: T_st_ad = k*T_st/mec2  ! Adimensional stellar temperature
real(dp), parameter :: R_st = 15.d0*R_sun     ! Generic R_star 
real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4. ! Stellar luminosity
real(dp), parameter :: Mdot = 1.d-6*M_sun_yr  ! Generic M^dot
real(dp), parameter :: v_w =  2000.d5         ! Generic wind terminal velocity

real(dp), parameter :: v_st_rot = 0.1d0*v_w   ! Rule-of-thumb
real(dp), parameter :: mu_w = 1.3_dp          ! Leitherer 1995 
real(dp), parameter :: Zq_w = 1.0_dp          ! Z = rms ionic charge
real(dp), parameter :: gamma_w = 1.0_dp       ! Leitherer 1995
real(dp), parameter :: filling = 1.d-1        ! ..... 10% ?

! ISM:
! ****
real(dp), parameter :: v_ISM = v_st           ! 
real(dp), parameter :: n_ISM = 10.d0           ! 
real(dp), parameter :: T_ISM = 40.d0          ! ?????
real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2! Adimensional ISM temperature 
real(dp), parameter :: alpha_ISM = 1.d-2      ! P_thermal_ISM = alpha*P_kin_ISM
real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)

real(dp), parameter :: mu_ISM = 1.3_dp         ! ? 
real(dp), parameter :: Zq_ISM = 1.0_dp         ! Z = rms ionic charge
real(dp), parameter :: gamma_ISM = 1.0_dp      ! ?

real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2.) ) 
real(dp), parameter :: H0 = R0*0.2             ! shock width


! Shocked stellar-wind (sw):
! **************************
real(dp), parameter :: mu_sw = 0.615_dp    ! mean atomic weight = 4/(6x+y+2) for complete ionization 
real(dp), parameter :: mu_e_sw = 1.173_dp  ! mu_e = 2/(1+x)
real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
real(dp), parameter :: Zq_sw = 1.23_dp        
real(dp), parameter :: gamma_sw = 1.34_dp

real(dp), parameter :: rho_w = Mdot/(4.d0*pi*R0**2.*v_w)
real(dp), parameter :: n_a = 4.*rho_w/(mp*mu_sw) ! ambient nuclei number density


! Heated dust (IR):
! *****************    
real(dp), parameter :: T_IR = 100.d0           ! K2017
real(dp), parameter :: T_IR_ad = k*T_IR/mec2  ! Adimensional dust temperature
real(dp), parameter :: L_IR = L_st/1.d3     ! Dust luminosity, K2017
real(dp), parameter :: R_IR = sqrt( L_IR/(4.*pi*sigma*T_IR**4.) ) ! cavity -> 1/4*pi, not 16*pi

! For IC (both diste.f90 and radic.f90)
real(dp), parameter :: Kappa_st = (R_st/(2.*R0))**2.
real(dp), parameter :: Kappa_IR = ( R_IR / R0 )**2. ! Is it ok ????????
real(dp), parameter :: Kappa_CMB = 1.d0             ! ????????


! Non-thermal particles:
! **********************
real(dp), parameter :: efi_ac_sw = 2.d0*pi*(c/v_w)**2   ! particle accel. efficiency at the sw 
real(dp), parameter :: frac_sw = sqrt(0.1_dp)     ! fraction UB_sw/UB_equipartition_sw
real(dp), parameter :: alpha_sw = 2.0_dp    ! injection spectral index in the sw 

real(dp), parameter :: frac_NT = 0.1_dp            ! Fraction of L_inj_perp converted to L_NT
real(dp), parameter :: iny_e = frac_NT*0.5d0       ! Fraction L_NT_e/L_NT
real(dp), parameter :: iny_p = frac_NT*0.5d0       ! Fraction L_NT_p/L_NT 

!real(dp), parameter :: B = frac_sw*sqrt(4.d0*pi*rho_w*4.*v_w**2.d0) 
       ! If U_mag = frac*U_th -> B^2/4pi = frac*(g_ad/g_ad - 1)*P 
       ! If P_mag = frac*P_th -> B^2/8pi = frac*P 
real(dp), parameter :: B = frac_sw*sqrt(Mdot*v_w/R0**2.) ! del Valle 2012
real(dp), parameter :: B_par = B*par
real(dp), parameter :: R_acel = R0 !H0  ! NOT SURE

!===================================================

! grid/vectors size
integer, parameter :: mE_e = 256    ! energy resolution for electrons
integer, parameter :: mE_f = 512    ! energy resolution for br and IC > mE_p

integer, parameter :: mE_p = 256    ! energy resolution for electrons

!===================================================

integer :: i,j,l

! dist
real(dp), parameter :: Emin_e = 2.d0*mec2   ! ~1 MeV
real(dp) :: Eint_e,dE_e
real(dp), dimension(:), allocatable :: E_e,N_e 

real(dp), parameter :: Emin_p = 1.d9*eV     ! 1 GeV
real(dp) :: Eint_p,dE_p
real(dp), dimension(:), allocatable :: E_p,N_p 


!===================================================
! emission

real(dp) :: P_Eps, dP_Eps

! syn
integer, parameter :: m_nu = 600   ! frequency resolution
real(dp), dimension(:), allocatable :: nu, Lsy, Lsy_abs

! br
real(dp), dimension(:), allocatable :: Lbr

! IC
real(dp) :: Lic(mE_f)
!real(dp), dimension(:), allocatable :: Lic

! p-p
real(dp), dimension(:), allocatable :: Lpp


!===================================================
! absorption 

!f-f
integer, parameter :: mE_ff = m_nu, mr_ff = 1700
real(dp), dimension(:), allocatable :: tauff
real(dp) :: alphaff_w(mr_ff,m_nu), r_w(mr_ff)

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
