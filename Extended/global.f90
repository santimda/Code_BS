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

   !===========================================================
   ! Simulation grid/vectors size

   integer, parameter :: m_theta = 200 ! angular resolution in the CD/shock
   integer, parameter :: mE_e = 200    ! energy resolution for electrons
   integer, parameter :: mE_p = 200    ! energy resolution for protons < mE_f 
   integer, parameter :: mE_f = 300    ! it is the same for br, p-p and IC
   integer, parameter :: mphi = 4      ! amount of linear-emitters
   real(dp) :: delta_phi = 0.5*pi/(mphi-1.d0) ! angular separation between linear-emitters
   real(dp) :: d_phi = 0.5*pi/mphi     ! angular width of each linear-emitter or wedge

   real(dp), parameter :: f_Ny = 0.6d0 ! Nyquist frequency for convolved maps

   !===========================================================


   integer :: i,j,l

   ! CD
   real(dp), dimension(:), allocatable :: x,y,z
   real(dp), dimension(:), allocatable :: theta,r_norm,R_theta,H_norm,H_theta

   ! termo
   real(dp), dimension(:), allocatable :: P,rho,T,v_t,v_perp,B,B_par,alpha_sh,P_th_w,P_med,C_s_med,Mach,C_s_sh,r

   ! dist
   real(dp), parameter :: Emin_e = 2.d0*mec2   ! ~1 MeV
   real(dp), parameter :: Emin_e_inj = Emin_e  ! Emin_e_inj >= Emin_e
   real(dp), parameter :: Emax_e = 1.d14*eV    ! 100 TeV
   real(dp) :: Eint_e,dE_e
   real(dp), dimension(:,:), allocatable :: N_etot ! Ne(theta,E_i)
   real(dp), dimension(:), allocatable :: E_e,N_e  ! N_e(E_i) is for a fixed "l"

   real(dp), parameter :: Emin_p = 1.d9*eV     ! 1 GeV 
   real(dp), parameter :: Emax_p = 5.d15*eV    ! 5 PeV
   real(dp) :: Eint_p,dE_p
   real(dp), dimension(:,:), allocatable :: N_ptot ! Np(theta,E_i)
   real(dp), dimension(:), allocatable :: E_p,N_p  ! N_p(E_i) is for a fixed "l" 

   !=======================================================
   ! emission

   real(dp) :: P_Eps, dP_Eps, P_Ef, dP_Ef

   ! syn
   integer, parameter :: m_nu = 600  ! frequency resolution for synchrotron photons
   real(dp), dimension(:,:), allocatable :: Lsy      ! L_sy(l,E_i)
   real(dp), dimension(:), allocatable :: nu,Lsy_tot,Lsy_tot_abs ! Lsy_tot(i) = sum_l L_sy(l,i)

   ! IC
   real(dp), dimension(:), allocatable :: LIC,LIC_tot,LIC_tot_abs
   real(dp) :: Kappa_st, Kappa_IR(m_theta) ! = (R/2D)**2 
   real(dp) :: X_st(3),X_emi(3),X_eE(3),X_obs(3) ! Position vectors for absorption integrals
   real(dp), parameter :: cte_Ndot = 2.*re**2*mec2**3/(pi*h_bar**3*c**2)
   logical, parameter :: ani = .True., iso = .False. 

   ! br
   real(dp), dimension(:), allocatable :: Lbr,Lbr_tot,Lbr_tot_abs
   logical, parameter :: scr = .True., unscr = .False.

   ! pp
   real(dp), dimension(:), allocatable :: Lpp,Lpp_tot,Lpp_tot_abs

   !=======================================================

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

   !=======================================================

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
      E_vec(1) = Emin
      do i_vec = 2, N_puntos
         E_vec(i_vec) = E_vec(i_vec-1)*Eint    
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
      DO I = 2, N-1, 2
         S1 = S1 + FI(I-1)
         S0 = S0 + FI(I)
         S2 = S2 + FI(I+1)
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
      real(dp), parameter:: c_ani = 6.13d0

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
      F_1 = ( (pi**2/6. + x0)*exp(-x0) )*g   
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
      F_3 = (pi**2/6.*(1.+c_3*x0)/(1.+pi**2*c_3*x0/6.)*exp(-x0) )*g   
      RETURN
   end Function F_3

   Function F_4(x0)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_4,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_4 = 6.62
      g = 1./(1. + (0.726*x0**0.461)/(1. + 0.382*x0**1.457))
      F_4 = (pi**2/6.*(1.+c_4*x0)/(1.+pi**2*c_4*x0/6.)*exp(-x0) )*g   
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
      nBB = 2.d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(k*T))-1.))
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

      suma = 0.d0
      do i = 1, N
         resta(i) = V1(i) - V2(i)
         suma = suma + resta(i)**2
      enddo
      suma = sqrt(suma)          ! suma = || versor_resta ||
      do j = 1,N
         versor(j) = resta(j)/suma
      enddo

      RETURN
   END SUBROUTINE versorresta

   !============================================================================
   ! p-p cross section

   Function sigma_pp(L1,Ep)
      real(dp), intent(in) :: L1, Ep
      real(dp) :: sigma_pp
      if (Ep < E_th) then !E_th= 1.22 GeV
         sigma_pp = 0.d0
      else           
         sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
         sigma_pp = sigma_pp*(1.d0-(E_th/Ep)**4)**2 ! Add factor for better fit at low energies 
      endif
      return
   end Function sigma_pp


   !============================================================================


END module global
