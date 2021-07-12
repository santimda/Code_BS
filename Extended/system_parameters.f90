module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   real(dp), parameter :: distance = 1.d0*kpc           ! source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0  ! angle between the star motion and the LOS in rad   
   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination  ! angle that I use
   ! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star
   !real(dp), parameter :: R0_proj = 5.d0*60.d3          ! projected BS distance to the star (5') in mas  
   !real(dp), parameter :: R0 = R0_proj*(distance/kpc)*UA/abs(cos(incli)) ! linear distance to the BS stagnation point
   real(dp), parameter :: theta_max = (170.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

   !--------------------------------------------------------------------------------------
   ! Star
   !real(dp), parameter :: v_st = 40.d5             ! del Palacio+2018
   !real(dp), parameter :: v_st = 30.d5             ! del Palacio+2018 (generic)
   !real(dp), parameter :: v_st = 500.d5             ! B0-05-d
   !real(dp), parameter :: v_st = 1000.d5             ! B1-1-d/B1-1-h/B1-1-nm
   real(dp), parameter :: v_st = 60000.d5             ! B2-60-d
   !real(dp), parameter :: T_st = 2.9d4              ! B0-05-d
   !real(dp), parameter :: T_st = 2.5d4              ! B1-1-d/B1-1-h/B1-1-nm
   real(dp), parameter :: T_st = 2.d4               ! B2-60-d
   !real(dp), parameter :: T_st = 4.07d4            ! K2017  
   real(dp), parameter :: T_st_ad = k*T_st/mec2    ! Adimensional stellar temperature
   !real(dp), parameter :: L_st = 7.d5*L_sun       ! del Valle 2012
   !real(dp), parameter :: R_st = sqrt( L_st/(4.*pi*sigma*T_st**4) )
   !real(dp), parameter :: R_st = 5.5*R_sun         ! B0-05-d
   !real(dp), parameter :: R_st = 4.8*R_sun         ! B1-1-d/B1-1-h/B1-1-nm
   real(dp), parameter :: R_st = 4.*R_sun          ! B2-60-d
   !real(dp), parameter :: R_st = 5.5*R_sun         ! K2017  
   real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4 
   !real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4 ! Stellar luminosity 
   !real(dp), parameter :: L_st = 5.d5*L_sun       ! del Palacio+ 2018 (generic)
   !real(dp), parameter :: L_st = 8.7d5*L_sun       ! K2018
   !real(dp), parameter :: Mdot = 1.d-8*M_sun_yr     ! B0-05-d
   !real(dp), parameter :: Mdot = 1.d-9*M_sun_yr      ! B1-1-d/B1-1-h/B1-1-nm
   real(dp), parameter :: Mdot = 1.d-10*M_sun_yr     ! B2-60-d
   !real(dp), parameter :: Mdot = 9.d-6*M_sun_yr    ! Peri 2014
   !real(dp), parameter :: Mdot = 1.d-8*M_sun_yr    ! Benaglia2020 (2.d-5 for alpha>0.3)
   !real(dp), parameter :: Mdot = 5.5d-6*M_sun_yr  ! Brookes 2016
   !real(dp), parameter :: v_w =  2200.d5          ! del Valle 2012
   !real(dp), parameter :: v_w =  2325.d5          ! Peri 2014
   !real(dp), parameter :: v_w = 2300.d5            ! Brookes2016,Ben2010
   !real(dp), parameter :: v_w = 1500.d5            ! K2018
   !real(dp), parameter :: v_w = 2000.d5            ! del Palacio+2018 (generic)
   !real(dp), parameter :: v_w = 1500.d5             ! B0-05-d
   !real(dp), parameter :: v_w = 1200.d5             ! B1-1-d/B1-1-h/B1-1-nm
   real(dp), parameter :: v_w = 1000.d5             ! B2-60-d
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: T_wind = 0.3*T_st        ! Drew1990
   real(dp), parameter :: v_st_rot = 0.1d0*v_w     ! Rule-of-thumb
   real(dp), parameter :: mu_w = 1.5_dp
   !real(dp), parameter :: mu_w = 1.5_dp            ! Leitherer 1995 
   real(dp), parameter :: Zq_w = 1.0_dp            ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp         ! Leitherer 1995
   real(dp), parameter :: filling = 1.d-1          ! ..... 10% ?
   real(dp), parameter :: r_A_st = 1.0_dp          ! r_A/R_st

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st             ! 
   !real(dp), parameter :: n_ISM = 1.d0            ! del Valle 2012
   !real(dp), parameter :: n_ISM = 3.5d0           ! Peri 2015
   !real(dp), parameter :: n_ISM = 13.9d0          ! Brookes2016
   !real(dp), parameter :: n_ISM = 19.d0           ! K2018
   !real(dp), parameter :: n_ISM = 15.d0           ! dP2018 = K2018*(1.32/1.7)
   real(dp), parameter :: n_ISM = 10.d0            ! Ben2020
   !real(dp), parameter :: n_ISM = 0.1d0            ! Galactic halo 
   !real(dp), parameter :: n_ISM = 100.d0           ! Molecular cloud
   real(dp), parameter :: mu_ISM = 1.28_dp         ! K2018; 0.5 for pure H ionized
   !real(dp), parameter :: mu_ISM = 2.35_dp         ! Molecular cloud 
   !real(dp), parameter :: T_ISM = 1.d4             ! 10000 K = warm ISM
   real(dp), parameter :: T_ISM = 1.d2             ! Cold ISM
   !real(dp), parameter :: T_ISM = 1.d0             ! del Palacio+2018 (generic)
   !real(dp), parameter :: T_ISM = 1.d4             ! warm ISM (Brookes p.136 uses T=8e3)
   real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2  ! Adimensional ISM temperature 
   !real(dp), parameter :: U_ISM = 4.d0*sigma*T_ISM**4/c ! average energy density
   !real(dp), parameter :: alpha_ISM = 1.d-10
   real(dp), parameter :: alpha_ISM = k*T_ISM/(v_ISM**2*mp*mu_ISM)! = P_thermal_ISM/P_kin_ISM
   !real(dp), parameter :: alpha_ISM =  5.d-2      ! P_thermal_ISM = alpha*P_kin_ISM ~ kT_ISM/(mu*m_p*v_ISM**2)
   !real(dp), parameter :: alpha_ISM = 0.1
   real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)
   real(dp), parameter :: Zq_ISM = 1.0_dp           ! Z = rms ionic charge
   real(dp), parameter :: gamma_ISM = 1.0_dp        ! Pure H

   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
   real(dp), parameter :: P_kin_ISM = rho_ISM*v_ISM**2 ! ISM ram pressure
   real(dp), parameter :: P_th_ISM = alpha_ISM*P_kin_ISM ! ISM thermal pressure 
   !real(dp), parameter :: C_s_ISM = sqrt(g_ad*P_th_ISM/rho_ISM) ! ISM sound speed
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
   real(dp), dimension(:), allocatable :: T_IR,T_IR_ad

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), dimension(:), allocatable :: efi_ac_sh! particle accel. efficiency at the shocks  (l)
   !real(dp), parameter :: frac_sw = sqrt(0.999_dp) ! fraction sqrt(UB_sw/UB_equipartition_sw)
   real(dp), parameter :: frac_sh = sqrt(0.1_dp) ! fraction sqrt(UB_sh/UB_equipartition_sh)
   !real(dp), parameter :: frac_sw = sqrt(0.128_dp) ! fraction sqrt(UB_sw/UB_equipartition_sw), f_NT=1
   !real(dp), parameter :: alpha_sh = 2._dp        ! injection spectral index in the shocks
   real(dp), parameter :: frac_NT = 1.d-1
   !real(dp), parameter :: frac_NT = 1.d0          ! Fraction L_NT/L_inj_perp, alpha=0.05, eta_B=0.15
   !real(dp), parameter :: frac_NT = 1.53d-1       ! Fraction L_NT/L_inj_perp, alpha=0.05, eta_B=1
   !real(dp), parameter :: frac_NT = 1.14d-1          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=0.1
   !real(dp), parameter :: frac_NT = 7.8d-1          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=0.01
   !real(dp), parameter :: frac_NT = 4.6d-2          ! Fraction L_NT/L_inj_perp, alpha=0.01, eta_B=1
   !real(dp), parameter :: frac_NT = 4.28d-2          ! Fraction L_NT/L_inj_perp, alpha=0.05, eta_B=1
   !real(dp), parameter :: frac_NT = 2.86d-2          ! Fraction L_NT/L_inj_perp, alpha=0.2, eta_B=1
   !real(dp), parameter :: frac_NT = 1.33d-2          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=1
   real(dp), parameter :: inj_e = frac_NT*0.05d0   ! Fraction L_NT_e/L_iny_perp
   !real(dp), parameter :: inj_e = frac_NT*0.5d0   ! Fraction L_NT_e/L_iny_perp, del Palacio+2018 (generic)
   real(dp), parameter :: inj_p = frac_NT-inj_e    ! Fraction L_NT_p/L_iny_perp 


   !============================================================================


END module system_parameters
