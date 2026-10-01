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

   integer, parameter :: m_theta = 300 ! angular resolution in the CD/shock
   integer, parameter :: mE_e = 300    ! energy resolution for electrons
   integer, parameter :: mE_p = 300    ! energy resolution for protons < mE_f 
   integer, parameter :: mE_f = 500    ! it is the same for br, p-p and IC
   integer, parameter :: mE_ph = 500   ! energy resolution for thermal f-f emission in the shocks 
   integer, parameter :: mphi = 30      ! amount of linear-emitters
   real(dp) :: delta_phi = 0.5*pi/(mphi-1.d0) ! angular separation between linear-emitters
   real(dp) :: d_phi = 0.5*pi/mphi     ! angular width of each linear-emitter or wedge

   !real(dp), parameter :: f_Ny = 0.7d0 ! Nyquist frequency for convolved maps
   real(dp), parameter :: f_Ny = 0.7d0

   !===========================================================

   integer :: i,j,l
   real(dp) :: dummy

   ! CD
   real(dp), dimension(:), allocatable :: x,y,z
   real(dp), dimension(:), allocatable :: theta,r_norm,R_theta,H_norm

   ! termo
   real(dp), dimension(:), allocatable :: P,rho,T,v_t,v_perp,B,B_par,alpha_sh,P_th_w
   real(dp), dimension(:), allocatable :: P_med, rho_med, C_s_med,Mach,C_s_sh,C_s_sh_new,r
   real(dp), dimension(:), allocatable :: v_t_RH, P_RH, rho_RH, B_RH, rho_min, rho_IL
   real(dp), dimension(:), allocatable :: v_t_new, P_new, rho_new, B_new
   real(dp), dimension(:), allocatable :: H_theta, H_new, y_eff, H_1

   ! dist
   real(dp), parameter :: Emin_e = 2.d0*mec2   ! ~1 MeV
   real(dp), parameter :: Emin_e_inj = Emin_e  ! Emin_e_inj >= Emin_e
   real(dp), parameter :: E_break = 1.9d9*eV    ! Divides both injection regimes
   !real(dp), parameter :: E_break = 2.d9*eV     ! May be this value?
   real(dp), parameter :: Emax_e = 1.d15*eV    ! 1 PeV
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
   integer, parameter :: m_nu = mE_ph  ! frequency resolution for synchrotron photons
   real(dp), dimension(:,:), allocatable :: Lsy,Lsy_RS,Lsy_FS      ! L_sy(l,E_i)
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

   ! average flux maps
   real(dp), parameter :: nu_av_min = 4.d9
   real(dp), parameter :: nu_av_max = 12.d9

   !=======================================================

   ! g-g
   integer, parameter :: mE_gg = 600, m_ang_gg = 400
   real(dp), dimension(:,:), allocatable :: tau
   real(dp), parameter :: Ef_min_gg = 5.d7*eV  ! 50 MeV 
   real(dp), parameter :: Ef_max_gg = 5.d13*eV ! 50 TeV 
   real(dp) :: Ef_gg(mE_gg),Ef_int_gg

   ! f-f
   integer, parameter :: mE_ff = m_nu, mr_ff = 1640, mr_map = 1000
   real(dp), dimension(:,:), allocatable :: tauff, tau_IL, tau_tot
   real(dp), dimension(:), allocatable :: tau_str
   real(dp), parameter :: nu_min_ff = 10.d6    ! 10 MHz 
   real(dp), parameter :: nu_max_ff = 1.d12   ! 1 THz
   real(dp) :: alphaff_w(mr_ff,m_nu), r_w(mr_ff)
   real(dp), dimension(:), allocatable :: nu_ISM
   real(dp), dimension(:,:), allocatable :: F_nu,F_nu_phi,F_nu_ion,F_nu_phi_ion
   !real(dp) :: Lsy_sw(m_theta,m_nu),Lsy_sism(m_theta,m_nu),Lsy_tot_abs(m_nu) 

   !=======================================================

   contains

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
   end Function F_1

   Function F_2(x0)  ! Used for anisotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_2,g
      real(dp), intent(in):: x0
      g = 1./(1. + (1.33*x0**0.691)/(1. + 0.534*x0**1.668))
      F_2 = ( (pi**2./6. + x0)*exp(-x0) )*g   
   end Function F_2

   Function F_3(x0)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_3,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_3 = 0.319
      g = 1./(1. + (0.443*x0**0.606)/(1. + 0.54*x0**1.481))
      F_3 = (pi**2/6.*(1.+c_3*x0)/(1.+pi**2*c_3*x0/6.)*exp(-x0) )*g   
   end Function F_3

   Function F_4(x0)  ! Used for isotropic IC (Khangulyan 2014)
      integer, parameter :: dp = kind(1.d0)
      real(dp) :: F_4,g
      real(dp), intent(in):: x0
      real(dp), parameter:: c_4 = 6.62
      g = 1./(1. + (0.726*x0**0.461)/(1. + 0.382*x0**1.457))
      F_4 = (pi**2/6.*(1.+c_4*x0)/(1.+pi**2*c_4*x0/6.)*exp(-x0) )*g   
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
   end Function nBB


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
   end Function sigma_pp


   !============================================================================
   ! lambda_cool

   Function lambda_cool(T)
      real(dp), intent(in) :: T
      real(dp) :: lambda_cool
      
      if (T < 1.d4) then
            lambda_cool = 4.d-29 * T**0.8     ! Muller+ 2020
         else if ((1.d4 <= T) .AND. (T <= 1.d5)) then    ! Myasnikov+ 1998
            lambda_cool = 7.d-27 * T
         else if ((1.d5 < T) .AND. (T <= 4.d7)) then
            lambda_cool = 7.d-19 * T**(-0.6) 
         else if (T > 4.d7) then
            lambda_cool = 3.d-27 * sqrt(T)
      end if

   end Function lambda_cool


   !============================================================================
   ! lambda_ff

   Function lambda_ff(T)
      real(dp), intent(in) :: T
      real(dp) :: lambda_ff

      lambda_ff = 3.d-27 * sqrt(T)

   end Function lambda_ff

   !============================================================================
   !============================================================================

END module global
