module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   real(dp), parameter :: distance = 1.32d0*kpc       ! source distance K2017
   real(dp), parameter :: inclination = 30.d0*pi/180.d0     ! angle between the star motion and the LOS in rad   
   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination   ! angle that I use
   ! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star
   !real(dp), parameter :: R0_proj = 5.d0*60.d3      ! projected BS distance to the star (5') in mas  
   !real(dp), parameter :: R0 = R0_proj*(distance/kpc)*UA/abs(cos(incli)) ! linear distance to the BS stagnation point
   real(dp), parameter :: theta_max = (110.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2

   !--------------------------------------------------------------------------------------
   ! Star
   real(dp), parameter :: v_st = 70.d5           ! Peri 2014 
   real(dp), parameter :: T_st = 4.07d4           ! K2017  
   real(dp), parameter :: T_st_ad = k*T_st/mec2   ! Adimensional stellar temperature
   !real(dp), parameter :: L_st = 7.d5*L_sun       ! del Valle 2012
   !real(dp), parameter :: R_st = sqrt( L_st/(4.*pi*sigma*T_st**4.) ) 
   real(dp), parameter :: R_st = 19.*R_sun         ! K2017  
   real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4 ! Stellar luminosity 
   real(dp), parameter :: Mdot = 6.5d-6*M_sun_yr  ! Peri 2014
   !real(dp), parameter :: Mdot = 5.5d-6*M_sun_yr  ! Brookes 2016
   !real(dp), parameter :: v_w =  2200.d5         ! del Valle 2012
   real(dp), parameter :: v_w =  2325.d5          ! Peri 2014
   real(dp), parameter :: frac_BS = 0.5d0        ! Fraction of wind kinetic power injeted in the BS
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: L_inj = frac_BS*L_wind ! Wind kinetic power injeted in the BS
   real(dp), parameter :: v_st_rot = 0.1d0*v_w   ! Rule-of-thumb
   real(dp), parameter :: mu_w = 1.3_dp          ! Leitherer 1995 
   real(dp), parameter :: Zq_w = 1.0_dp          ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp       ! Leitherer 1995
   !real(dp), parameter :: ab_X=0.705d0, ab_Y=0.275d0, ab_Z=0.02 ! stellar abundances 
   real(dp), parameter :: filling = 1.d-1        ! ..... 10% ?

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st           ! 
   !real(dp), parameter :: n_ISM = 1.d0           ! del Valle 2012 
   !real(dp), parameter :: n_ISM = 3.5d0           ! Peri 2015
   real(dp), parameter :: n_ISM = 2.d0           ! ?????
   real(dp), parameter :: T_ISM = 40.d0          ! ?????
   real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2! Adimensional ISM temperature 
   !real(dp), parameter :: U_ISM = 4.d0*sigma*T_ISM**4/c ! average energy density
   real(dp), parameter :: alpha_ISM = 1.d-2      ! P_thermal_ISM = alpha*P_kin_ISM
   real(dp), parameter :: lambda = alpha_ISM/(alpha_ISM + 1.d0)

   real(dp), parameter :: mu_ISM = 1.3_dp         ! ? 
   real(dp), parameter :: Zq_ISM = 1.0_dp         ! Z = rms ionic charge
   real(dp), parameter :: gamma_ISM = 1.0_dp      ! ?

   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM ! CHECK
   real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2.) ) 
   real(dp), parameter :: H0 = R0*0.3             ! shock width

   !--------------------------------------------------------------------------------------
   ! Shocked stellar-wind (sw)
   real(dp), parameter :: mu_sw = 0.615_dp      ! mean atomic weight = 4/(6x+y+2) for complete ionization 
   real(dp), parameter :: mu_e_sw = 1.173_dp    ! mu_e = 2/(1+x)
   real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
   real(dp), parameter :: Zq_sw = 1.23_dp        
   real(dp), parameter :: gamma_sw = 1.34_dp
   real(dp), parameter :: T_6 = 1.d0            ! temperature in MK, should scale with v_w?

   real(dp), parameter :: rho_w = Mdot/(4.d0*pi*R0**2*v_w)
   real(dp), parameter :: n_a = 4.*rho_w/(mp*mu_sw) ! ambient nuclei number density
   real(dp), parameter :: z_br = 2.24               ! Padovani+2018 for Wiilms+2000 abundances
   real(dp), parameter :: z_pp = 2.17               ! Padovani+2018 for Wiilms+2000 abundances

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)   
   real(dp), parameter :: T_IR = 98.d0           ! K2017   
   !real(dp), parameter :: T_IR = 27.*(0.2)**(-1./6.)*(L_st/1.d38)**(1./6.)*(R0/pc)**(-1./3.) ! Eq.10 del Valle 2012
   real(dp), parameter :: T_IR_ad = k*T_IR/mec2  ! Adimensional dust temperature
   !real(dp), parameter :: L_IR = 4.d0*pi*H0**2.*sigma*T_IR**4. ! Dust luminosity ?????
   !real(dp), parameter :: L_IR = 5.7d36 ! Dust luminosity, observed ?????
   real(dp), parameter :: L_IR = L_st/350.d0     ! Dust luminosity, K2017
   real(dp), parameter :: R_IR = sqrt( L_IR/(4.*pi*sigma*T_IR**4) ) ! cavity -> 1/4*pi, not 16*pi

   ! For IC (both diste.f90 and radic.f90)
   real(dp), parameter :: Kappa_st = (R_st/(2.*R0))**2
   real(dp), parameter :: Kappa_IR = ( R_IR / R0 )**2  
   real(dp), parameter :: Kappa_CMB = 1.d0             

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: efi_ac_sw = 2.d0*pi*(c/v_w)**2   ! particle accel. efficiency at the sw 
   real(dp), parameter :: frac_sw = sqrt(0.1_dp)       ! fraction UB_sw/UB_equipartition_sw
   real(dp), parameter :: alpha_sw = 2.0_dp            ! injection spectral index in the sw 

   real(dp), parameter :: frac_NT = 0.1_dp             ! Fraction of L_inj_perp converted to L_NT
   real(dp), parameter :: inj_e = frac_NT*0.05d0       ! Fraction L_NT_e/L_NT
   real(dp), parameter :: inj_p = frac_NT*0.5d0        ! Fraction L_NT_p/L_NT 

   !real(dp), parameter :: B = frac_sw*sqrt(4.d0*pi*rho_w*4.*v_w**2) 
   ! If U_mag = frac*U_th -> B^2/4pi = frac*(g_ad/(g_ad - 1))*P 
   ! If P_mag = frac*P_th -> B^2/8pi = frac*P 
   real(dp), parameter :: B = frac_sw*sqrt(Mdot*v_w/R0**2) ! del Valle 2012
   real(dp), parameter :: B_par = B*par
   real(dp), parameter :: R_acel = R0  !


   !============================================================================

END module system_parameters
