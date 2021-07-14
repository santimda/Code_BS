module global

   use constants
   use system_parameters
   implicit none
   public

   !====================================================================
   ! This module contains constants (from constants.f90), variables and 
   ! parameters that all other modules use. 
   ! It also has defined some subroutines that many modules call.
   !====================================================================

   integer, parameter :: nth = 4   ! number of threads (parallelization)

   !===================================================
   ! simulation grid/vectors size

   integer, parameter :: mE_e = 500    ! energy resolution for electrons
   integer, parameter :: mE_f = 512    ! energy resolution for br and IC > mE_p
   integer, parameter :: mE_p = 400    ! energy resolution for electrons

   !===================================================

   integer :: i,j,l

   ! dist
   real(dp), parameter :: Emin_e = 1.d6*ev   ! ~1 MeV
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
   real(dp), dimension(:), allocatable :: nu,L_sy,L_sy_RT,L_sy_abs

   ! br
   real(dp), dimension(:), allocatable :: L_br, Ef 
   !logical, parameter :: screened = .True.   ! screened nuclei
   logical, parameter :: screened = .False. ! unscreened nuclei

   ! IC
   real(dp), dimension(:), allocatable :: L_ic
   logical, parameter :: ani = .True., iso = .False. 

   ! p-p
   real(dp), dimension(:), allocatable :: Lpp, Epi


   !===================================================
   ! absorption 

   !f-f
   integer, parameter :: mE_ff = m_nu, mr_ff = 1000
   real(dp), dimension(:), allocatable :: tauff
   real(dp) :: alphaff_w(mr_ff,m_nu), r_w(mr_ff)

   !===================================================

contains

   !============================================================================
   ! vector_log
   !============================================================================

   SUBROUTINE vector_log(Emin,Emax,Eint,E_vec)
      implicit none
      real(dp), intent(in) :: Emin,Emax
      integer :: N_puntos, i_vec
      real(dp), intent(out) :: Eint
      real(dp), intent(out) :: E_vec(:)

      N_puntos = SIZE(E_vec)
      Eint = (Emax/Emin)**(1.d0/dfloat(N_puntos-1))

      E_vec(1) = Emin
      do i_vec=2,N_puntos
         E_vec(i_vec)=E_vec(i_vec-1)*Eint
      end do

      RETURN

   END SUBROUTINE vector_log


   !============================================================================
   ! integrar
   !============================================================================

   SUBROUTINE integrar(N,H,FI,S)
      ! Subroutine for integration over f(x) with the Simpson rule.
      ! FI: integrand f(x); H: interval; S: integral.
      integer, intent(in) :: N
      real(dp), intent(in) :: H,FI(N)
      real(dp), intent(out) :: S
      real(dp) :: S0,S1,S2

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
      real(dp) :: nBB,Epho,T
      nBB = 2.d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(k*T))-1.))
      RETURN
   end Function nBB

   !============================================================================
   !  Cooling_IC

   Function F_ani(u)  ! Used for anisotropic IC (Khangulyan 2014)
      real(dp) :: F_ani
      real(dp), intent(in):: u
      real(dp), parameter:: c_ani = 6.13
      F_ani = (c_ani*u*log(1.+2.16*u/c_ani)) / (1.+c_ani*u/0.822)
      RETURN
   end Function F_ani

   Function F_iso(u)  ! Used for isotropic IC (Khangulyan 2014)
      real(dp) :: F_iso,g_u
      real(dp), intent(in):: u
      real(dp), parameter:: c_iso = 5.68
      g_u = 1./(1. + (-0.362*u**0.682)/(1. + 0.826*u**1.281))
      F_iso = (c_iso*u*log(1.+0.722*u/c_iso)) / (1.+c_iso*u/0.822)*g_u ! Corregido
      RETURN
   end Function F_iso

   !   N^dot_IC

   Function F_1(x0)  ! Used for anisotropic IC (Khangulyan 2014)
      real(dp) :: F_1,g
      real(dp), intent(in):: x0
      g = 1./(1. + (0.153*x0**0.857)/(1. + 0.254*x0**1.84))
      F_1 = ( (pi**2/6. + x0)*exp(-x0) )*g   
      RETURN
   end Function F_1

   Function F_2(x0)  ! Used for anisotropic IC (Khangulyan 2014)
      real(dp) :: F_2,g
      real(dp), intent(in):: x0
      g = 1./(1. + (1.33*x0**0.691)/(1. + 0.534*x0**1.668))
      F_2 = ( (pi**2/6. + x0)*exp(-x0) )*g   
      RETURN
   end Function F_2

   Function F_3(x0)  ! Used for isotropic IC (Khangulyan 2014)
      real(dp) :: F_3,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_3 = 0.319
      g = 1./(1. + (0.443*x0**0.606)/(1. + 0.54*x0**1.481))
      F_3 = (pi**2/6.*(1.+c_3*x0)/(1.+pi**2*c_3*x0/6.)*exp(-x0) )*g   
      RETURN
   end Function F_3

   Function F_4(x0)  ! Used for isotropic IC (Khangulyan 2014)
      real(dp) :: F_4,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_4 = 6.62
      g = 1./(1. + (0.726*x0**0.461)/(1. + 0.382*x0**1.457))
      F_4 = (pi**2/6.*(1.+c_4*x0)/(1.+pi**2*c_4*x0/6.)*exp(-x0) )*g   
      RETURN
   end Function F_4

   !============================================================================
   !  Diffusion coefficients
   !============================================================================

   Function D_diff(E,B)
      real(dp) :: D_diff
      real(dp), intent(in) :: E, B 
      real(dp), parameter :: delta=1./3.d0, E0=1.d10*eV !Kolm MW normalization   
      !character :: diff_regime --> Defined in system_parameters
      if (diff_regime == 'K') then           !Kolmogorov   
         D_diff = 3.86d28*(E/E0)**delta      !cm^2 s^-1
      else if (diff_regime == 'B') then      ! Bohm diff
         D_diff = E*c/(3.d0*qe*B)            !cm^2 s^-1
      else if (diff_regime == 'S') then      ! del Valle+2018
         D_diff = 1.d25*(E/E0)**0.5d0      !cm^2 s^-1
      else if (diff_regime == 'F') then      ! del Valle+2018
         D_diff = 1.d27*(E/E0)**0.5d0      !cm^2 s^-1
      else 
         write(*,*) 'ERROR, UNKNOWN DIFFUSION'
         stop
      endif

      RETURN
   end Function D_diff

   !============================================================================



END module global
