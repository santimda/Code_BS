module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   real(dp), parameter :: distance = 2.21d0*kpc          ! source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0  ! angle between the star motion and the LOS in rad   
   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination  ! angle that I use
   ! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star
   !real(dp), parameter :: R0_proj = 5.d0*60.d3          ! projected BS distance to the star (5') in mas  
   !real(dp), parameter :: R0 = R0_proj*(distance/kpc)*UA/abs(cos(incli)) ! linear distance to the BS stagnation point
   real(dp), parameter :: theta_max = (170.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

   !--------------------------------------------------------------------------------------
   ! Star
   ! Catalogued as O6.5V by Maíz Apellániz+ 2016
   real(dp), parameter :: v_st = 138.d5            ! Gaia3 - Peri+2015
   real(dp), parameter :: T_st = 3.83d4            ! http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
   real(dp), parameter :: T_st_ad = k*T_st/mec2    ! Adimensional stellar temperature
   real(dp), parameter :: L_st = 1.9d5*L_sun      ! 
   real(dp), parameter :: R_st = 9.82*R_sun        ! 
   real(dp), parameter :: Mdot = 2.4d-8*M_sun_yr    ! Van den Eijnden 2021
   real(dp), parameter :: v_w = 2545.d5                ! Van den Eijnden 2021
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: T_wind = 0.3*T_st        ! Drew1990
   real(dp), parameter :: v_st_rot = 0.1d0*v_w     ! Rule-of-thumb
   real(dp), parameter :: mu_w = 1.5_dp            ! Leitherer 1995 
   real(dp), parameter :: Zq_w = 1.0_dp            ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp         ! Leitherer 1995
   real(dp), parameter :: filling = 1.d-1          ! ..... 10% ?
   real(dp), parameter :: r_A_st = 1.0_dp          ! r_A/R_st

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st             ! 
   real(dp), parameter :: n_ISM = 380.d0           ! VdE 2021
   real(dp), parameter :: T_ISM = 1.d2             ! MC?
   real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2  ! Adimensional ISM temperature 
   real(dp), parameter :: alpha_ISM = 5.d-2        ! P_thermal_ISM = alpha*P_kin_ISM ~ kT_ISM/(mu*m_p*v_ISM**2)
   real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)
   real(dp), parameter :: mu_ISM = 1.38_dp          ! K2018; 0.5 for pure H ionized
   real(dp), parameter :: Zq_ISM = 1.0_dp           ! Z = rms ionic charge
   real(dp), parameter :: gamma_ISM = 1.0_dp        ! Pure H

   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
   real(dp), parameter :: P_kin_ISM = rho_ISM*v_ISM**2 ! ISM ram pressure
   real(dp), parameter :: P_th_ISM = alpha_ISM*P_kin_ISM ! ISM thermal pressure 
   real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2) ) 
   real(dp), parameter :: R0_proj = (R0/AU)/(distance/pc)*abs(cos(incli)) ! angular distance ["] to the BS stagnation point
   ! Peri 2014: R0_proj ~ 3.5 arcmin; K2018: R0_proj ~ 193" ~ 3.2 arcmin; Ben2020: ~3.9 arcmin
   real(dp), parameter :: L_ISM = 0.5*pi*rho_ISM*R0**2*v_ISM**3 ! 

   !--------------------------------------------------------------------------------------
   ! Shocked stellar-wind (sw)
   real(dp), parameter :: mu_sw = 0.615_dp         ! mean atomic weight = 4/(6x+y+2) for complete ionization 
   real(dp), parameter :: mu_e_sw = 1.173_dp       ! mu_e = 2/(1+x)
   real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
   real(dp), parameter :: Zq_sw = 1.23_dp        
   real(dp), parameter :: gamma_sw = 1.34_dp
   real(dp), parameter :: T_6 = 1.d0               ! temperature in MK, should scale with v_w?

   real(dp), parameter :: z_br = 2.24               ! Padovani+2018 for Wiilms+2000 abundances
   real(dp), parameter :: z_pp = 2.17               ! Padovani+2018 for Wiilms+2000 abundances

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)   
   !real(dp), parameter :: L_IR = 5.7d36 ! Dust luminosity, observed ?????
   real(dp), parameter :: L_IR = L_st/350.d0     ! Dust luminosity, K2017
   real(dp), parameter :: mu_SISM = 0.615d0        ! Assume solar abundance (?)
   real(dp), dimension(:), allocatable :: T_IR,T_IR_ad

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), dimension(:), allocatable :: efi_ac_sh! particle accel. efficiency at the shocks  (l)
   !real(dp), parameter :: frac_sw = sqrt(0.999_dp) ! fraction sqrt(UB_sw/UB_equipartition_sw)
   real(dp), parameter :: frac_sh = sqrt(0.1)      ! fraction sqrt(UB_sw/UB_equipartition_sw), f_NT=1
   !real(dp), parameter :: alpha_sw = 2.2_dp        ! injection spectral index in the sw 
   !real(dp), parameter :: frac_NT = 1.53d-1       ! Fraction L_NT/L_inj_perp, alpha=0.05, eta_B=1
   real(dp), parameter :: frac_NT = 1.14d-1          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=0.1
   !real(dp), parameter :: frac_NT = 7.8d-1          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=0.01
   !real(dp), parameter :: frac_NT = 4.6d-2          ! Fraction L_NT/L_inj_perp, alpha=0.01, eta_B=1
   !real(dp), parameter :: frac_NT = 4.28d-2          ! Fraction L_NT/L_inj_perp, alpha=0.05, eta_B=1
   !real(dp), parameter :: frac_NT = 2.86d-2          ! Fraction L_NT/L_inj_perp, alpha=0.2, eta_B=1
   !real(dp), parameter :: frac_NT = 1.33d-2          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=1
   real(dp), parameter :: inj_e = frac_NT*0.1d0   ! Fraction L_NT_e/L_iny_perp 
   real(dp), parameter :: inj_p = frac_NT-inj_e    ! Fraction L_NT_p/L_iny_perp 


   !============================================================================


END module system_parameters
