module common_parameters

   use constants
   use system_parameters
   implicit none
   public

   !====================================================================
   ! This module contains common parameters defined for every system
   !====================================================================

   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination ! angle that I use

   ! Relative abundances (for solar composition)
   real(dp), parameter :: X_ab = 0.73d0, Y_ab = 0.25d0
   real(dp), parameter :: Z_ab = 1.d0 - X_ab - Y_ab ! X + Y + Z = 1

   !--------------------------------------------------------------------------------------
   ! Star
   real(dp), parameter :: T_st_ad = k*T_st/mec2 ! Adimensional stellar temperature
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: T_wind = 0.3*T_st ! Drew1990
   real(dp), parameter :: v_st_rot = 0.1d0*v_w ! Rule-of-thumb
   real(dp), parameter :: mu = 4.d0/(3.d0 + 5.d0*X_ab - Z_ab) ! mean atomic weight = 4/(3+5X-Z) for complete ionization
   real(dp), parameter :: mu_e = 2.d0/(1.d0 + X_ab) ! mu_e = 2/(1+X) only for Z=0; mu_e > mu!
   real(dp), parameter :: mu_i = 1.d0/(1.d0/mu - 1.d0/mu_e)
   real(dp), parameter :: mu_w = 1.5_dp ! Leitherer 1995
   real(dp), parameter :: Zq_w = 1.0_dp ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp ! Leitherer 1995
   real(dp), parameter :: filling = 1.d-1 ! ..... 10% ?
   real(dp), parameter :: r_A_st = 1.0_dp ! r_A/R_st

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st
   real(dp), parameter :: alpha_ISM = k*T_ISM/(v_ISM**2*mp*mu_ISM) ! = P_thermal_ISM/P_kin_ISM
   real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)
   real(dp), parameter :: zq_ISM = 1.d0 ! Z = rms ionic charge

   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM
   real(dp), parameter :: P_kin_ISM = rho_ISM*v_ISM**2 ! ISM ram pressure
   real(dp), parameter :: P_th_ISM = alpha_ISM*P_kin_ISM ! ISM thermal pressure
   real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2) )
   real(dp), parameter :: R0_proj = (R0/AU)/(distance/pc)*abs(cos(incli)) ! angular distance ["] to the BS stagnation point
   real(dp), parameter :: L_ISM = 0.5*pi*rho_ISM*R0**2*v_ISM**3

   !--------------------------------------------------------------------------------------
   ! Shocked stellar-wind (sw)
   real(dp), parameter :: mu_sw = mu ! mean atomic weight for complete ionization
   real(dp), parameter :: Zq_sw = 1.11_dp    ! later squared 1.23_dp

   real(dp), parameter :: z_br = 2.24 ! Padovani+2018 for Wilms+2000 abundances
   real(dp), parameter :: z_pp = 2.17 ! Padovani+2018 for Wilms+2000 abundances

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), dimension(:), allocatable :: T_IR,T_IR_ad

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), dimension(:), allocatable :: efi_ac_sh
   real(dp), parameter :: lambda_diff = 1.d0
   real(dp), parameter :: inj_p = frac_NT-inj_e ! Fraction L_NT_p/L_inj_perp

   !============================================================================


END module common_parameters
