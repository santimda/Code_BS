module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'RCW-49-S1' ! CHECK: identified by elimination
   real(dp), parameter :: distance = 8.d0*kpc ! +- 3.4 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: theta_max = (145.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O5III, Peri+2015)
   real(dp), parameter :: v_st = 95.d5 ! Gaia3: 130 +- 60 km/s
   real(dp), parameter :: T_st = 4.1d4 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 9.82*R_sun
   real(dp), parameter :: L_st = 3.5d5*L_sun
   real(dp), parameter :: Mdot = 3.5d-6*M_sun_yr ! Peri+2015
   real(dp), parameter :: v_w = 2800.d5 ! Peri+2015

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 2.d0
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
   real(dp), parameter :: frac_H = 0.1d0

   !============================================================================


END module system_parameters
