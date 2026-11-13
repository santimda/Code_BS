module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'BD+60'
   real(dp), parameter :: distance = 2.99d0*kpc ! +- 0.05 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 56.d0*pi/180.d0 ! Green+2019
   real(dp), parameter :: dec = 61.19459355028689d0*pi/180.
   real(dp), parameter :: theta_max = (125.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O6.5V, Sota+2011)
   real(dp), parameter :: v_st = 42.d5 ! Gaia3 + Moutzouri+2022 + Green+2013
   real(dp), parameter :: T_st = 3.83d4 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 9.82*R_sun
   real(dp), parameter :: L_st = 1.9d5*L_sun
   real(dp), parameter :: Mdot = 1.3d-6*M_sun_yr ! Toala+2020
   real(dp), parameter :: v_w = 2000.d5 ! Moutzouri+2022

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 27.d0
   real(dp), parameter :: T_ISM = 8.d3
   real(dp), parameter :: mu_ISM = 0.6033d0 ! = mu (complete ionization, solar X, Z)
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/200.d0 ! Dust luminosity: no measurement for this system, L_st/L_IR similar to BD+43 (Kobulnicky+2017)

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   real(dp), parameter :: alpha_1 = -1.d0
   real(dp), parameter :: frac_NT = 0.2d0
   real(dp), parameter :: inj_e = frac_NT*0.2d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.1d0

   !============================================================================


END module system_parameters
