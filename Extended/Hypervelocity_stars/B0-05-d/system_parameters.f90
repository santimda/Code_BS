module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character(len=*), parameter :: system = 'B0-05-d'
   real(dp), parameter :: distance = 5.d0*kpc
   real(dp), parameter :: inclination = 60.d0*pi/180.d0
   real(dp), parameter :: theta_max = (170.d0/180.d0)*pi ! extension of the BS modelled

   !--------------------------------------------------------------------------------------
   ! Star
   real(dp), parameter :: v_st = 500.d5
   real(dp), parameter :: T_st = 2.9d4
   real(dp), parameter :: R_st = 5.5*R_sun
   real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4
   real(dp), parameter :: Mdot = 1.d-8*M_sun_yr
   real(dp), parameter :: v_w = 1500.d5

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: n_ISM = 10.d0
   real(dp), parameter :: T_ISM = 1.d2
   real(dp), parameter :: mu_ISM = 1.28d0
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)
   real(dp), parameter :: L_IR = L_st/206.d0 ! Dust luminosity: no measurement for this system, L_st/L_IR of BD+43 (Kobulnicky+2017) scaled to 1.72 kpc

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.1_dp) ! Magnetic field value
   real(dp), parameter :: alpha_1 = 2.d0
   real(dp), parameter :: frac_NT = 1.d-1
   real(dp), parameter :: inj_e = frac_NT*0.05d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.15d0

   !============================================================================


END module system_parameters
