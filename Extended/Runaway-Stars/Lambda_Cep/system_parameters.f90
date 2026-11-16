module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'Lambda_Cep'
   real(dp), parameter :: distance = 0.85d0*kpc ! +- 0.02 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: dec = 59.41443924203745d0*pi/180.
   real(dp), parameter :: theta_max = (120.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O6.5I(n)fp)
   real(dp), parameter :: v_st = 70.d5 ! S-A
   real(dp), parameter :: T_st = 38.3d3 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 10.*R_sun ! K2017
   real(dp), parameter :: L_st = 5.2d5*L_sun ! S-A+2018
   real(dp), parameter :: Mdot = 1.d-6*M_sun_yr ! S-A+2018
   real(dp), parameter :: v_w = 2000.d5 ! S-A+2018

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
   ! FS: .true. = radiative (isothermal layer, no NT particles), .false. = adiabatic
   logical, parameter :: radiative_FS = .true.

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   ! Injection spectral index: Q_inj \propto E^-alpha
   real(dp), parameter :: alpha_1 = -1.d0 ! Injection index below E_break
   real(dp), parameter :: alpha_2 = 3.d0 ! Injection index above E_break
   real(dp), parameter :: frac_NT = 0.25_dp
   real(dp), parameter :: inj_e = frac_NT*1.d-1

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.01d0

   !============================================================================


END module system_parameters
