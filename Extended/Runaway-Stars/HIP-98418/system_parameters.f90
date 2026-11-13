module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'HIP-98418' ! CHECK: identified by elimination
   real(dp), parameter :: distance = 2.21d0*kpc ! +- 0.03 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: theta_max = (130.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star (O6.5V, Maiz Apellaniz+2016)
   real(dp), parameter :: v_st = 60.d5 ! Gaia3, Peri+2015
   real(dp), parameter :: T_st = 3.83d4 ! Pecaut & Mamajek
   real(dp), parameter :: R_st = 9.82*R_sun
   real(dp), parameter :: L_st = 1.9d5*L_sun
   real(dp), parameter :: Mdot = 2.4d-7*M_sun_yr ! Peri+2015
   real(dp), parameter :: v_w = 2545.d5 ! Peri+2015

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 10.d0 ! VdE2021
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
