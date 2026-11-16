module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   real(dp), parameter :: distance = 1.72d0*kpc          ! source distance (Gaia)
   real(dp), parameter :: inclination = 70.d0*pi/180.d0  ! angle between the star motion and the LOS in rad   
   real(dp), parameter :: theta_max = (120.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

   !--------------------------------------------------------------------------------------
   ! Star
   ! Catalogued as O9V by Comeron & Pasquali 2012
   real(dp), parameter :: v_st = 17.7d5            ! Gaia3 - Peri+2015
   real(dp), parameter :: T_st = 3.33d4            ! http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
   real(dp), parameter :: L_st = 6.6d4*L_sun      ! 
   real(dp), parameter :: R_st = 7.72*R_sun        ! 
   real(dp), parameter :: Mdot = 5.d-8*M_sun_yr    ! Peri+2015
   real(dp), parameter :: v_w = 1500.d5                ! Peri+2015

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st             ! 
   real(dp), parameter :: n_ISM = 2.d0             ! Peri+ 2015
   real(dp), parameter :: T_ISM = 1.d4             ! 
   real(dp), parameter :: mu_ISM = 1.28             ! MC?
   real(dp), parameter :: mu_SISM = 0.615d0

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)   
   real(dp), parameter :: L_IR = L_st/200.d0     ! Dust luminosity, K2017

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: frac_sh = sqrt(0.2) ! Magnetic field value
   r! Injection spectral index: Q_inj \propto E^-alpha
   real(dp), parameter :: alpha_1 = -1.d0 ! Injection index below E_break
   real(dp), parameter :: alpha_2 = 3.d0 ! Injection index above E_break
   real(dp), parameter :: frac_NT = 0.13_dp ! NT fraction
   real(dp), parameter :: inj_e = frac_NT*1.d-1 ! Fraction L_NT_e/L_inj_perp

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 0.25d0


   !============================================================================


END module system_parameters