module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   ! G3 (HD 319881) according to Peri+2015
   !====================================================================

   character(len=*), parameter :: system = 'G3'
   real(dp), parameter :: distance = 1.79d0*kpc ! +- 0.06 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: dec = -34.541799065885996d0*pi/180.
   real(dp), parameter :: theta_max = (90.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O5V, Peri+2015)
   real(dp), parameter :: v_st = 22.d5 ! Gaia3
   real(dp), parameter :: T_st = 41.1d3 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 11.45*R_sun ! Pecaut & Mamajek
   real(dp), parameter :: L_st = 3.5d5*L_sun ! Pecaut & Mamajek
   real(dp), parameter :: Mdot = 4.d-7*M_sun_yr ! Peri+2015
   real(dp), parameter :: v_w = 2000.d5 ! Peri+2015

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 30.d0
   real(dp), parameter :: T_ISM = 8.d3
   real(dp), parameter :: mu_ISM = 0.6033d0 ! = mu (complete ionization, solar X, Z)
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/206.d0 ! Dust luminosity: no measurement for this system, L_st/L_IR of BD+43 (Kobulnicky+2017) scaled to 1.72 kpc

   !--------------------------------------------------------------------------------------
   ! FS: .true. = radiative (isothermal layer, no NT particles), .false. = adiabatic
   logical, parameter :: radiative_FS = .true.

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   ! Injection spectral index: Q_inj \propto E^-alpha
   real(dp), parameter :: alpha_1 = -1.d0 ! Injection index below E_break
   real(dp), parameter :: alpha_2 = 3.d0 ! Injection index above E_break
   real(dp), parameter :: frac_NT = 0.1d0
   real(dp), parameter :: inj_e = frac_NT*0.1d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.45d0

   !============================================================================


END module system_parameters
