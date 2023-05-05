
module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   ! G3 according to Peri+ 2015
   !====================================================================

   character(len=*), parameter :: system = 'G3'
   real(dp), parameter :: distance = 1.79d0*kpc  ! +- 0.06 kpc source distance (Gaia)
   real(dp), parameter :: inclination = 60.d0*pi/180.d0  ! angle between the star motion and the LOS in rad   
   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination  ! angle that I use
   ! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star
   real(dp), parameter :: theta_max = (90.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

   ! Relative abundances (for solar composition)
   real(dp), parameter :: X_ab = 0.73d0, Y_ab = 0.25d0
   real(dp), parameter :: Z_ab = 1.d0 - X_ab - Y_ab ! X + Y + Z = 1

   !--------------------------------------------------------------------------------------
   ! Star
   ! Classified as O5V (Peri+ 2015)
   real(dp), parameter :: v_st = 22.d5            ! Gaia3
   real(dp), parameter :: T_st = 41.1d4            ! http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
   real(dp), parameter :: T_st_ad = k*T_st/mec2    ! Adimensional stellar temperature
   real(dp), parameter :: L_st = 3.5d5*L_sun       ! http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
   real(dp), parameter :: R_st = 11.45*R_sun       ! http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
   real(dp), parameter :: Mdot = 4.d-7*M_sun_yr    ! Peri+2015
   real(dp), parameter :: v_w = 2000.d5                ! Peri+ 2015
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: T_wind = 0.3*T_st        ! Drew1990
   real(dp), parameter :: v_st_rot = 0.1d0*v_w     ! Rule-of-thumb
   real(dp), parameter :: mu = 4.d0/(3.d0 + 5.d0*X_ab - Z_ab) ! mean atomic weight = 4/(3+5X-Z) for complete ionization
   real(dp), parameter :: mu_e = 2.d0/(1.d0 + X_ab)    ! mu_e = 2/(1+X) only for Z=0; mu_e > mu! 
   real(dp), parameter :: mu_i = 1.d0/(1.d0/mu - 1.d0/mu_e)
   real(dp), parameter :: mu_w = 1.5_dp            ! Leitherer 1995 
   real(dp), parameter :: Zq_w = 1.0_dp            ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp         ! Leitherer 1995
   real(dp), parameter :: filling = 1.d-1          ! ..... 10% ?
   real(dp), parameter :: r_A_st = 1.0_dp          ! r_A/R_st

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st             ! 
   real(dp), parameter :: n_ISM = 30.d0            ! Derived from R0_proj (Peri+ 2012)
   real(dp), parameter :: T_ISM = 8.d3             ! Warm ISM
   real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2  ! Adimensional ISM temperature
   real(dp), parameter :: mu_ISM = mu             ! 
   real(dp), parameter :: alpha_ISM = k*T_ISM/(v_ISM**2*mp*mu_ISM)! = P_thermal_ISM/P_kin_ISM
   real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)
   real(dp), parameter :: Zq_ISM = 0.79_dp           ! Z = rms ionic charge
   real(dp), parameter :: gamma_ISM = 1.0_dp        ! Pure H

   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
   real(dp), parameter :: P_kin_ISM = rho_ISM*v_ISM**2 ! ISM ram pressure
   real(dp), parameter :: P_th_ISM = alpha_ISM*P_kin_ISM ! ISM thermal pressure 
   real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2) ) 
   real(dp), parameter :: R0_proj = (R0/AU)/(distance/pc)*abs(cos(incli)) ! angular distance ["] to the BS stagnation point
   real(dp), parameter :: L_ISM = 0.5*pi*rho_ISM*R0**2*v_ISM**3 ! 

   !--------------------------------------------------------------------------------------
   ! Shocked stellar-wind (sw)
   real(dp), parameter :: mu_sw = mu! 0.615_dp         ! mean atomic weight = 4/(6x+y+2) for complete ionization 
   real(dp), parameter :: mu_e_sw = mu_e! 1.173_dp       ! mu_e = 2/(1+x)
   real(dp), parameter :: mu_i_sw = mu_i !1.d0/(1/mu_sw-1/mu_e_sw)  
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
   real(dp), parameter :: frac_sh = sqrt(0.1)      ! fraction sqrt(UB_sw/UB_equipartition_sw), f_NT=1
   real(dp), parameter :: frac_NT = 1.d-1          ! Fraction L_NT/L_inj_perp, alpha=0.35, eta_B=0.1
   real(dp), parameter :: alpha_1 = -3.d0
   real(dp), parameter :: inj_e = frac_NT*0.1d0   ! Fraction L_NT_e/L_iny_perp 
   real(dp), parameter :: inj_p = frac_NT-inj_e    ! Fraction L_NT_p/L_iny_perp
   character(len=*), parameter :: diff_regime = 'Bohm'      ! del Valle+2018
   real(dp), parameter :: lambda_diff = 1.d0

   !--------------------------------------------------------------------------------------
   ! Isothermal shock width parameter
   real(dp), parameter :: frac_H = 1.d0


   !============================================================================


END module system_parameters