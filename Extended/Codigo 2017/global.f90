module global

implicit none
public

!===================================================

integer, parameter :: dp = kind(1.d0)
integer, parameter :: nth = 2   ! number of threads (paralelizar)

!Constantes físicas:
real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: par = sqrt(2.d0/3.d0)
real(dp), parameter :: eV = 1.602176d-12
real(dp), parameter :: c = 2.99792458d10
real(dp), parameter :: k = 1.380658d-16    ! erg/K
real(dp), parameter :: qe = 4.8032068d-10
real(dp), parameter :: mec2 = 0.511d6*eV
real(dp), parameter :: re = 2.8d-13      ! radio clásico del electrón
real(dp), parameter :: sigma = 5.6705d-5
real(dp), parameter :: h = 6.6260755d-27
real(dp), parameter :: h_bar = h/(2.*pi)
real(dp), parameter :: alfa = 7.29735308d-3 
real(dp), parameter :: Msolyr = 1.989d33/31556926.d0 ! (m_sol/yr) en cgs
real(dp), parameter :: Rsol = 6.96d10    ! radio solar en cm
real(dp), parameter :: UA = 1.496d13     ! 1 UA en cm
real(dp), parameter :: kpc = 3.0856d21   ! 1 kpc en cm
real(dp), parameter :: mp = 1.673d-24    ! masa del protón en CGS
real(dp), parameter :: ga = 5.d0/3.d0    ! el coeficiente adiabático
real(dp), parameter :: m_pic2 = 134.98d6*eV
real(dp), parameter :: mpc2 = 1836.d0*mec2

!===================================================


! Parámetros del sistema:
! ***********************
real(dp), parameter :: distance = 2.3d0   ! 1 kpc = 3.0856d21 (la dejo en kpc)
real(dp), parameter :: incli = -15.d0*pi/180.d0 ! °  
real(dp), parameter :: closer = 1.d0 ! 36.d0/55.d0
real(dp), parameter :: D_proj = 36.d0/closer ! Sep orbital de la bin en mas en Aug2008 
! D_proj = 26.5 mas (epoch 2013.09 from Sana et al. 2014)
real(dp), parameter :: D = D_proj*distance*UA/abs(cos(incli)) ! en cm
real(dp), parameter :: x_max = 0.95d0*(D/UA)   ! Extensión de la CD
 real(dp), parameter :: eta = 0.5_dp            ! wind momentum rates ratio
! Es el valor que dan en Cohen 2011
! luminosidad del sistema: L_bol ~ 6.34e39 erg/s

!real(dp), parameter :: ab_X=0.705d0, ab_Y=0.275d0, ab_Z=0.02 ! abundancias asumidas
real(dp), parameter :: filling = 0.5_dp ! filling factor: (M_radio/M_X)^2 ~ 0.02
! LO RE MENTÍ
real(dp), parameter :: frac1 = 0.014_dp  ! fracción del B_equipartición que tomo
real(dp), parameter :: frac2 = 0.014_dp  ! fracción del B_equipartición que tomo
real(dp), parameter :: alpha = 3.2_dp  ! spectral index de iny (de obs en radio)
real(dp), parameter :: alpha2 = 2.0_dp ! spectral index de iny a HE

real(dp), parameter :: frac_NT = 0.114_dp!6d-3 ! Fracción de L_iny normal que va a L_NT
real(dp), parameter :: iny_e = frac_NT*0.5d0        ! Fracción de la L_NT en e-
real(dp), parameter :: iny_p = frac_NT*(1.d0-iny_e) ! Fracción de L_NT en p+

! Estrella 1 (O2I):
! *****************
real(dp), parameter :: T1 = 42500.d0          ! Benaglia 2015, Cohen 2011
!real(dp), parameter :: T1 = 44850.d0         ! Muijres 2012 extrapolado
real(dp), parameter :: Rst1 = 18.3d0*Rsol     ! Muijres 2012 extrapolado
!real(dp), parameter :: Mdot1 = 7.d-6*Msolyr  ! valor derivado de R-X
!real(dp), parameter :: Mdot1 = 3.2d-5*Msolyr ! valor derivado de radio
real(dp), parameter :: Mdot1 = 1.d-5*Msolyr   ! promedio
real(dp), parameter :: vinf1 =  3200.d5       ! valor que usan en Cohen 2011
! Puedo extrapolar vinf1 de Muijres 2012 y me daría 3365.d5
real(dp), parameter :: T1mcc = k*T1/mec2
real(dp), parameter :: Lst1 = 4.d0*pi*Rst1**2*sigma*T1**4
real(dp), parameter :: v_rot1 = 0.1d0*vinf1   ! 380.d5 ! Ekstrom para M > 85 Mo
! White & Chen 1995, genérico usan 250 km/s
real(dp), parameter :: mu_w1 = 1.3_dp         ! Leitherer 1995 
real(dp), parameter :: Zq_w1 = 1.0_dp         ! Z = rms ionic charge
real(dp), parameter :: gamma_w1 = 1.0_dp      ! Leitherer 1995

! Estrella 2 (O3.5V):
! *******************
real(dp), parameter :: T2 = 44000.d0          ! Muijres 2012 interpolado
real(dp), parameter :: Rst2 = 14.d0*Rsol    ! Muijres 2012 interpolado
real(dp), parameter :: vinf2 = 2800.d5        ! Prinja 1990 interpolado
!real(dp), parameter :: vinf2 = 3600.d5        ! Muijres 2012 interpolado
!real(dp), parameter :: Mdot2 = 3.2d-6*Msolyr  ! Muijres 2012 interpolado(Vink)
!real(dp), parameter :: Mdot2 = 1.7d-5*Msolyr ! Radio (Benaglia2015) 
!real(dp), parameter :: Mdot2 = 5.2d-6*Msolyr   ! medio usando el eta 
real(dp), parameter :: Mdot2 = Mdot1*(vinf1/vinf2)*eta ! usando el eta 
real(dp), parameter :: T2mcc = k*T2/mec2
real(dp), parameter :: Lst2 = 4.d0*pi*Rst2**2*sigma*T2**4
real(dp), parameter :: v_rot2 = 0.1d0*vinf2   ! 350.d5  ! Ekstrom para M > 40 Mo
! White & Chen 1995, genérico usan 250 km/s
real(dp), parameter :: mu_w2 = 1.3_dp         ! Leitherer 1995 
real(dp), parameter :: Zq_w2 = 1.0_dp         ! Z = rms ionic charge
real(dp), parameter :: gamma_w2 = 1.0_dp      ! Leitherer 1995

! WCR:
!*****
real(dp), parameter :: mu_1 = 0.615_dp   ! mean atomic weight = 4/(6x+y+2) for  
! complete ionization (se usa en la CD, no en el wind)
real(dp), parameter :: mu_e1 = 1.173_dp  ! mu_e = 2/(1+x)
real(dp), parameter :: mu_i1 = 1.d0/(1/mu_1-1/mu_e1)  
real(dp), parameter :: Zq_1 = 1.23_dp        
real(dp), parameter :: gamma_1 = 1.34_dp
real(dp), parameter :: mu_2 = 0.615_dp  
real(dp), parameter :: mu_e2 = 1.173_dp  ! mu_e = 2/(1+x)
real(dp), parameter :: mu_i2 = 1.d0/(1/mu_2-1/mu_e2)   
real(dp), parameter :: Zq_2 = 1.23_dp        
real(dp), parameter :: gamma_2 = 1.34_dp

real(dp), parameter :: efi_ac1 = 2.d0*pi*(c/vinf1)**2 ! eficiencia de aceleración 1
real(dp), parameter :: efi_ac2 = 2.d0*pi*(c/vinf2)**2 ! eficiencia de aceleración 2

!===================================================


! tamaño de las grillas
integer, parameter :: ml = 1300   ! resolución en la CD
integer, parameter :: mE = 160   ! resolución en E de los electrones
integer, parameter :: mE_p = 450 ! resolución en E de los protones > mE_f 
integer, parameter :: mphi = 2  ! cantidad de emisores línea
real(dp) :: delta_phi = (pi/2.d0)/(mphi-1.d0) ! ancho angular de c/emisor línea


!===================================================


integer :: i,j,l

! CD
real(dp), dimension(:), allocatable :: x,y,z

! termo
real(dp), dimension(:), allocatable :: P,rho,T,v,B,B_par

! dist
real(dp), parameter :: Emin_e = 100.d0*mec2  
real(dp), parameter :: Estiff_e = 7.d14*mec2 !
real(dp), parameter :: Emax_e = 1.d14*eV    ! 100 TeV
real(dp) :: Eint_e,dE_e
real(dp), dimension(:,:), allocatable :: N_etot ! Este es Ne(l,E_i)
real(dp), dimension(:), allocatable :: E_e,N_e  ! N_e(E_i) lo uso a un l fijo

real(dp), parameter :: Emin_p = 1.d0*mpc2   !1.d9*eV  ! 1 GeV
real(dp), parameter :: Estiff_p = 1.1d0*mpc2!7.d3*mec2 !
real(dp), parameter :: Emax_p = 1.d16*eV    ! 10 PeV
real(dp) :: Eint_p,dE_p
real(dp), dimension(:,:), allocatable :: N_ptot ! Este es Np(l,E_i)
real(dp), dimension(:), allocatable :: E_p,N_p  ! N_p(E_i) lo uso a un l fijo

!===================================================
! emisión

real(dp) :: P_Eps, dP_Eps

! syn
integer, parameter :: mnu = 600  ! resolución en nu de los fotones
real(dp), dimension(:,:), allocatable :: Lsy      ! L_sy(l,E_i)
real(dp), dimension(:), allocatable :: nu,Lsy_tot ! Lsy_tot(i) = sum_l L_sy(l,i)

! IC
integer, parameter :: mE_low = 20
integer, parameter :: mE_tot = mE_low + mE
integer, parameter :: mE_ph = 16    ! Para la planckeana de cada estrella
real(dp), dimension(:), allocatable :: LIC,LIC_tot,LIC_tot_abs
real(dp) :: Xst1(3),Xst2(3),Xemi(3),XeE(3),Xobs(3) ! También sirve para Br
real(dp) :: Kappa1,Kappa2 ! = (R/2D)**2.
real(dp), parameter :: cte_Ndot = 2.*re**2.*mec2**3./(pi*h_bar**3.*c**2.)

! br
integer, parameter :: mE_f = 440 ! Uso el mismo para p-p
real(dp), dimension(:), allocatable :: Lbr,Lbr_tot,Lbr_tot_abs

! pp
real(dp), dimension(:), allocatable :: Lpp,Lpp_tot,Lpp_tot_abs

!===================================================
! absorción

! g-g
integer, parameter :: mE_gg = 600, m_ang_gg = 400
real(dp), dimension(:,:), allocatable :: tau1,tau2
real(dp), parameter :: Ef_min_gg = 5.d7*eV ! 50 MeV 
real(dp), parameter :: Ef_max_gg = 5.d13*eV ! 50 TeV 
real(dp) :: Ef_gg(mE_gg),Ef_int_gg

! f-f
integer, parameter :: mE_ff = mnu, mr_ff = 800
real(dp), dimension(:,:), allocatable :: tauff
real(dp) :: alphaffw1(mr_ff,mnu),alphaffw2(mr_ff,mnu)
real(dp) :: r_1(mr_ff),r_2(mr_ff)
real(dp) :: Lsy1(ml,mnu),Lsy2(ml,mnu),Lsy_tot_abs(mnu) 

!===================================================

contains

!============================================================================
! vectorE
!============================================================================
  
  SUBROUTINE vectorE(Emin,Emax,Eint,E_vec)
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

  END SUBROUTINE vectorE


!============================================================================
! integrar
!============================================================================
  
  SUBROUTINE integrar(N,H,FI,S)
    ! Subroutine for integration over f(x) with the Simpson rule.
    ! FI: integrand f(x); H: interval; S: integral.
    real*8 :: S,S0,S1,S2,H
    integer :: N
    real*8 :: FI(N)
    S  = 0.d0
    S0 = 0.d0
    S1 = 0.d0
    S2 = 0.d0
    DO 100 I = 2, N-1, 2
       S1 = S1 + FI(I-1)
       S0 = S0 + FI(I)
       S2 = S2 + FI(I+1)
100    CONTINUE
       S = H*(S1 + 4.*S0 + S2)/3.
       ! If N is even, add the last slice separately
       IF(MOD(N,2).eq.0) S = S + H*(5.d0*FI(N) + 8.d0*FI(N-1) - FI(N-2))/12.d0
       RETURN

     End SUBROUTINE integrar

  !============================================================================
  

!============================================================================
!  Cooling_IC
    
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
    nBB = 2.d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(k*T))-1))
    RETURN
  end Function nBB

  
  !============================================================================
  
  
  SUBROUTINE versorresta(V1,V2,N,versor)
    !     subrutina para obtener el versor: (V1-V2)/|V1-V2|
    real*8 :: V1,V2,versor,resta,suma
    integer :: N
    DIMENSION V1(N),V2(N),versor(N),resta(N)
    
    suma=0.d0
    do i=1,N
       resta(i)=V1(i)-V2(i)
       suma=suma+resta(i)**2
    enddo
    suma=dsqrt(suma)          ! en suma guardo el módulo del vector resta
    do j=1,N
       versor(j)=resta(j)/suma
    enddo
    
    RETURN
  END SUBROUTINE versorresta
  
  
  !============================================================================
  


END module global
