module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'BD+43'
   real(dp), parameter :: distance = 1.72d0*kpc ! +- 0.02 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0 ! angle between the star motion and the LOS
   real(dp), parameter :: dec = 43.98538212854551d0*pi/180.
   real(dp), parameter :: theta_max = (120.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star
   real(dp), parameter :: v_st = 50.d5 ! Benaglia2021
   real(dp), parameter :: T_st = 4.07d4 ! K2017
   real(dp), parameter :: R_st = 19.*R_sun ! K2017
   real(dp), parameter :: L_st = 8.7d5*L_sun ! K2018
   real(dp), parameter :: Mdot = 4.d-6*M_sun_yr
   real(dp), parameter :: v_w = 3000.d5 ! K2018

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 6.d0
   real(dp), parameter :: T_ISM = 8.d3 ! warm ISM
   real(dp), parameter :: mu_ISM = 0.6033d0 ! = mu (complete ionization, solar X, Z)
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/350.!4.d0*pi*distance**2*5.0d-8 ! Dust luminosity, F_IR = 5e-8 erg/s/cm2 (Kobulnicky+2017, Table 5)

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   real(dp), parameter :: alpha_1 = -1.d0 ! Injection index below E_break
   real(dp), parameter :: frac_NT = 0.13_dp ! NT fraction
   real(dp), parameter :: inj_e = frac_NT*1.d-1 ! Fraction L_NT_e/L_inj_perp

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.25d0

   !============================================================================


END module system_parameters
