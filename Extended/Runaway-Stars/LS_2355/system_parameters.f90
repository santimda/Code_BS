module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'LS_2355'
   real(dp), parameter :: distance = 2.22d0*kpc ! +- 0.02 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: dec = -62.65272506566931d0*pi/180.
   real(dp), parameter :: theta_max = (120.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O6.5IV)
   real(dp), parameter :: v_st = 25.d5
   real(dp), parameter :: T_st = 38.3d3 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 10.*R_sun
   real(dp), parameter :: L_st = 5.2d5*L_sun
   real(dp), parameter :: Mdot = 1.d-6*M_sun_yr
   real(dp), parameter :: v_w = 2000.d5

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 40.d0
   real(dp), parameter :: T_ISM = 1.d2
   real(dp), parameter :: mu_ISM = 0.6033d0 ! = mu (complete ionization, solar X, Z)
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/206.d0 ! Dust luminosity: no measurement for this system, L_st/L_IR of BD+43 (Kobulnicky+2017) scaled to 1.72 kpc

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   real(dp), parameter :: alpha_1 = -1.d0
   real(dp), parameter :: frac_NT = 0.13_dp
   real(dp), parameter :: inj_e = frac_NT*1.d-1

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.01d0

   !============================================================================


END module system_parameters
