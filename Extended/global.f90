module global

   use constants
   use common_parameters
   implicit none
   public

   !====================================================================
   ! This module contains constants (from constants.f90), variables and 
   ! parameters that all other modules use. 
   ! It also has defined some subroutines that many modules call.
   !====================================================================

   integer, parameter :: nth = 4   ! number of threads (parallelization)
   character(len=*), parameter :: output_dir = 'output/'   ! directory for output files

   !===========================================================
   ! Simulation grid/vectors size

   integer, parameter :: m_theta = 300 ! angular resolution in the CD/shock
   integer, parameter :: mE_e = 300    ! energy resolution for electrons
   integer, parameter :: mE_p = 300    ! energy resolution for protons < mE_f 
   integer, parameter :: mE_f = 500    ! it is the same for br, p-p and IC
   integer, parameter :: mE_ph = 500   ! energy resolution for thermal f-f emission in the shocks 
   !integer, parameter :: mphi = 30      ! amount of linear-emitters
   !real(dp), parameter :: delta_phi = 0.5d0*pi/(mphi-1.d0) ! angular separation between linear-emitters
   !real(dp), parameter :: d_phi = 0.5d0*pi/mphi            ! angular width of each linear-emitter or wedge
   integer,  parameter :: mphi = 50                    ! rotation angles phi = 0, ..., pi/2 (linear emitters)
   integer,  parameter :: n_emit = 4*(mphi-1)          ! distinct linear emitters in a ring (phi, phi+pi and mirror images)
   real(dp), parameter :: delta_phi = 2.d0*pi/n_emit   ! angular separation between linear emitters = 0.5*pi/(mphi-1)
   real(dp), parameter :: d_phi = delta_phi            ! angular width of each linear emitter (wedge)

   real(dp), parameter :: f_Ny = 0.7d0 ! Nyquist frequency for convolved maps

   !===========================================================

   integer :: i,j,l
   real(dp) :: dummy

   ! CD
   real(dp), dimension(:), allocatable :: x,y,z
   real(dp), dimension(:), allocatable :: theta,r_norm,R_theta

   ! termo
   real(dp), dimension(:), allocatable :: P,rho,T,v_t,v_perp,B,B_par,alpha_sh
   real(dp), dimension(:), allocatable :: rho_med, C_s_med,Mach,C_s_sh,C_s_sh_new,r
   real(dp), dimension(:), allocatable :: v_t_RH, P_RH, rho_RH, B_RH, rho_min, rho_IL
   real(dp), dimension(:), allocatable :: v_t_new, P_new, rho_new, B_new
   real(dp), dimension(:), allocatable :: H_theta, y_eff, H_1
   real(dp), parameter :: T_ion = 8.d3   ! photoionization equilibrium T of the shocked ISM (minimum T of the IL)

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
   real(dp), dimension(:), allocatable :: nu

   ! IC
   real(dp), dimension(:), allocatable :: LIC,LIC_tot
   real(dp) :: Kappa_st, Kappa_IR(m_theta) ! = (R/2D)**2 
   !real(dp) :: X_st(3),X_emi(3),X_eE(3),X_obs(3) ! Position vectors for absorption integrals
   real(dp), parameter :: e_obs(3) = [ -sin(incli), cos(incli), 0.d0 ]   ! unit vector towards the observer (apex at -x; incli = 90 deg: BS in front of the star)
   real(dp), parameter :: cte_Ndot = 2.*re**2*mec2**3/(pi*h_bar**3*c**2)
   logical, parameter :: ani = .True., iso = .False. 

   ! br
   real(dp), dimension(:), allocatable :: Lbr,Lbr_tot
   logical, parameter :: scr = .True., unscr = .False.

   ! pp
   real(dp), dimension(:), allocatable :: Lpp,Lpp_tot

   ! average flux maps
   real(dp), parameter :: nu_av_min = 4.d9
   real(dp), parameter :: nu_av_max = 12.d9

   ! frequency grid nu(m_nu) for synchrotron and f-f absorption (written to nu.dat in initialize.f90)
   real(dp), parameter :: nu_min_syn = 10.d6
   real(dp), parameter :: nu_max_syn = 1.d24

   !=======================================================

   ! f-f
   integer, parameter :: mE_ff = m_nu, mr_ff = 1640, mr_map = 1000
   real(dp), parameter :: nu_min_ff = 10.d6    ! 10 MHz 
   real(dp), parameter :: nu_max_ff = 1.d12   ! 1 THz
   real(dp), dimension(:), allocatable :: nu_ISM

   !=======================================================

   contains

   !============================================================================
   ! Output files and input checks

   ! Path of an output file: output_dir//name
   function output_file(name) result(path)
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: path
      path = output_dir // name
   end function output_file

   ! Creates output_dir if it does not exist (call once, at the start of main)
   subroutine create_output_dir()
      call execute_command_line('mkdir -p '//output_dir)
   end subroutine create_output_dir

   ! Stops the program if shock is not 'RS' or 'FS' ('BOTH' is also
   ! accepted when allow_both = .true.)
   subroutine validate_shock(shock, routine, allow_both)
      character(len=*), intent(in) :: shock, routine
      logical, intent(in), optional :: allow_both
      logical :: both_ok

      both_ok = .false.
      if (present(allow_both)) both_ok = allow_both

      if (shock == 'RS' .or. shock == 'FS') return
      if (both_ok .and. shock == 'BOTH') return

      if (both_ok) then
         write(*,*) 'ERROR: invalid shock in ', trim(routine), &
            ': expected RS, FS or BOTH, received ', trim(shock)
      else
         write(*,*) 'ERROR: invalid shock in ', trim(routine), &
            ': expected RS or FS, received ', trim(shock)
      end if
      error stop
   end subroutine validate_shock

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
   !  Black-Body spectrum

   Function nBB(Epho,T)
      ! Black-body photon number density per unit energy divided by 4 [1/cm3/erg]:
      ! c*nBB is the photon flux per unit energy at the surface, and
      ! L_N = 4*pi*R**2*c*nBB the photon luminosity [photons/s/erg]
      real(dp) :: nBB
      real(dp), intent(in) :: Epho,T
      nBB = 2.d0*pi/((h*c)**3) * (Epho**2/(exp(Epho/(k*T))-1.d0))
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

   elemental function lambda_cool(T)
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

   elemental function lambda_ff(T)
      real(dp), intent(in) :: T
      real(dp) :: lambda_ff

      lambda_ff = 3.d-27 * sqrt(T)

   end Function lambda_ff

   !============================================================================
   !============================================================================

END module global