module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   ! G1 according to Peri+2015
   !====================================================================

   character(len=*), parameter :: system = 'G1'
   real(dp), parameter :: distance = 1.79d0*kpc ! +- 0.03 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: dec = -34.24376115929199d0*pi/180.
   real(dp), parameter :: theta_max = (80.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O7V, Peri+2015)
   real(dp), parameter :: v_st = 20.d5 ! Gaia3: 19 km/s
   real(dp), parameter :: T_st = 37.1d3 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 9.42*R_sun ! Pecaut & Mamajek
   real(dp), parameter :: L_st = 1.5d5*L_sun ! Pecaut & Mamajek
   real(dp), parameter :: Mdot = 2.2d-7*M_sun_yr ! Peri+2015: 2e-7
   real(dp), parameter :: v_w = 2100.d5 ! Peri+2015

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 30.d0
   real(dp), parameter :: T_ISM = 8.d3
   real(dp), parameter :: mu_ISM = 0.6033d0 ! = mu (complete ionization, solar X, Z)
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/102.d0*(distance/(1.717d0*kpc))**2 ! Dust luminosity, L_st/L_IR = 102 at 1.717 kpc (Kobulnicky+2019, #692)

   !--------------------------------------------------------------------------------------
   ! FS: .true. = radiative (isothermal layer, no NT particles), .false. = adiabatic
   logical, parameter :: radiative_FS = .true.
   
   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   ! Injection spectral index: Q_inj \propto E^-alpha
   real(dp), parameter :: alpha_1 = -1.d0 ! Injection index below E_break
   real(dp), parameter :: alpha_2 = 3.d0 ! Injection index above E_break
   real(dp), parameter :: frac_NT = 0.35_dp
   real(dp), parameter :: inj_e = frac_NT*0.1d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 1.d0

   !============================================================================


END module system_parameters
