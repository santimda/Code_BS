module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'Vela_X-1'
   real(dp), parameter :: distance = 2.02d0*kpc ! +- 0.06 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 45.d0*pi/180.d0
   real(dp), parameter :: dec = -40.5546529611754d0*pi/180.
   real(dp), parameter :: theta_max = (135.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (B0.5Ia, Houk 1978)
   real(dp), parameter :: v_st = 90.d5 ! Gaia3: v_t = 59.7 +- 2.2 km/s; Kretschmar+2021
   real(dp), parameter :: T_st = 2.9d4 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 6.48*R_sun
   real(dp), parameter :: L_st = 2.6d4*L_sun
   real(dp), parameter :: Mdot = 1.d-6*M_sun_yr ! Grinberg+2017
   real(dp), parameter :: v_w = 700.d5 ! Grinberg+2017

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 2.d0
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
   real(dp), parameter :: frac_NT = 0.13_dp
   real(dp), parameter :: inj_e = frac_NT*0.1d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.1d0

   !============================================================================


END module system_parameters
