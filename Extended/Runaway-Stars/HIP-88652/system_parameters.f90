module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'HIP-88652'
   real(dp), parameter :: distance = 1.46d0*kpc ! +- 0.04 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: theta_max = (110.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O9.7Ib, Sota+2011)
   real(dp), parameter :: v_st = 38.d5 ! Gaia3, Gontcharov 2006 (radial v)
   real(dp), parameter :: T_st = 3.19d4 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 7.5*R_sun
   real(dp), parameter :: L_st = 5.2d4*L_sun
   real(dp), parameter :: Mdot = 5.d-7*M_sun_yr ! Peri+2012
   real(dp), parameter :: v_w = 1535.d5 ! Peri+2012

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 4.d0 ! Peri+2012: 2
   real(dp), parameter :: T_ISM = 8.d3
   real(dp), parameter :: mu_ISM = 0.61d0
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/206.d0 ! Dust luminosity: no measurement for this system, L_st/L_IR of BD+43 (Kobulnicky+2017) scaled to 1.72 kpc

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.1) ! Magnetic field value
   real(dp), parameter :: alpha_1 = -1.d0 ! not in the previous file
   real(dp), parameter :: frac_NT = 1.d-1
   real(dp), parameter :: inj_e = frac_NT*0.1d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.06d0

   !============================================================================


END module system_parameters
