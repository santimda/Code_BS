module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   character, parameter :: diff_regime = 'B' ! diffusion regime (F, S, K, B)

   real(dp), parameter :: distance = 1.d0*kpc         ! source distance K2017
   real(dp), parameter :: inclination = 30.d0*pi/180.d0 ! angle between v_star and the LOS [rad]
   ! i = 0° is face-on, i = 90° corresponds to BS in front of the star, i = -90° is BS behind the star
   real(dp), parameter :: incli = pi/2.d0 - inclination ! angle that I use
   ! incli = 90° is BS in front of the star, incli = 0° is edge-on,  incli = -90° is BS behind the star real(dp), parameter :: theta_max = (110.d0/180.d0)*pi ! ~ extension of the BS modelled ~ pi/2
   real(dp), parameter :: cte_cs = sqrt(g_ad*3./16.d0)  ! c_s ~ sqrt(g_ad*P/rho_sh) = cte_cs*v (P~3/4*rho*v^2, rho_sh=4*rho)
   real(dp), parameter :: frac_BS = 0.5d0   ! Fraction of kinetic power injeted in the BS (geometric factor)

   !--------------------------------------------------------------------------------------
   ! Star
   real(dp), parameter :: v_st = 60000.d5             ! Peri 2014 
   real(dp), parameter :: T_st = 2.d4            ! K2017  
   real(dp), parameter :: T_st_ad = k*T_st/mec2     ! Adimensional stellar temperature
   !real(dp), parameter :: L_st = 7.d5*L_sun        ! del Valle 2012
   !real(dp), parameter :: R_st = sqrt( L_st/(4.*pi*sigma*T_st**4.) ) 
   real(dp), parameter :: R_st = 4.*R_sun            ! K2017  
   real(dp), parameter :: L_st = 4.d0*pi*R_st**2*sigma*T_st**4 ! Stellar luminosity 
   real(dp), parameter :: Mdot = 1.d-10*M_sun_yr   ! Peri 2014
   real(dp), parameter :: v_w =  1000.d5           ! Peri 2014
   real(dp), parameter :: L_wind = 0.5d0*(Mdot*v_w**2) ! Wind kinetic power
   real(dp), parameter :: L_inj_RS = frac_BS*L_wind    ! Wind kinetic power injeted in the RS
   real(dp), parameter :: v_st_rot = 0.1d0*v_w      ! Rule-of-thumb
   real(dp), parameter :: mu_w = 1.3_dp             ! Leitherer 1995 
   real(dp), parameter :: Zq_w = 1.0_dp             ! Z = rms ionic charge
   real(dp), parameter :: gamma_w = 1.0_dp          ! Leitherer 1995
   real(dp), parameter :: filling = 1.d-1           ! ..... 10% ?

   !--------------------------------------------------------------------------------------
   ! ISM
   real(dp), parameter :: v_ISM = v_st             ! ISM speed in the star's ref. frame
   real(dp), parameter :: n_ISM = 10.d0            ! 
   real(dp), parameter :: T_ISM = 1.d4             ! 10000 K = warm ISM
   real(dp), parameter :: T_ISM_ad = k*T_ISM/mec2  ! Adimensional ISM temperature 
   real(dp), parameter :: mu_ISM = 1.28_dp         ! mean atomic weight for an atomic cloud
   real(dp), parameter :: gamma_ISM = 1.0_dp       ! abundance parameter ???
   real(dp), parameter :: rho_ISM = n_ISM*mp*mu_ISM
   real(dp), parameter :: alpha_ISM = k*T_ISM/(v_ISM**2*mp*mu_ISM)! = P_thermal_ISM/P_kin_ISM
   real(dp), parameter :: R0 = sqrt( (Mdot*v_w)/(4.*pi*(1.+alpha_ISM)*rho_ISM*v_ISM**2) ) 
   real(dp), parameter :: H0 = R0*0.3             ! shock width
   real(dp), parameter :: L_ISM = 0.5*pi*rho_ISM*R0**2*v_ISM**3 ! Kinetic power injected in the ISM shock
   real(dp), parameter :: L_inj_FS = frac_BS*L_ISM          ! Wind kinetic power injeted in the BS

   !--------------------------------------------------------------------------------------
   ! Shocked stellar-wind (Reverse Shock)
   real(dp), parameter :: mu_sw = 0.615_dp      ! mean atomic weight = 4/(6x+y+2) for complete ionization 
   real(dp), parameter :: mu_e_sw = 1.173_dp    ! mu_e = 2/(1+x)
   real(dp), parameter :: mu_i_sw = 1.d0/(1/mu_sw-1/mu_e_sw)  
   real(dp), parameter :: Zq_sw = 1.23_dp        
   real(dp), parameter :: gamma_sw = 1.34_dp
   real(dp), parameter :: rho_w = Mdot/(4.d0*pi*R0**2*v_w)    
   real(dp), parameter :: n_sh_RS = 4.*rho_w/(mp*mu_sw)  ! number density (reverse shock)
   real(dp), parameter :: T_sh_RS = 1.36d7*(v_w/1.d8)**2  ! temperature (RS)
   real(dp), parameter :: c_s_RS = v_w*cte_cs            ! velocity of the shocked fluid (RS)
   real(dp), parameter :: z_br = 2.24               ! Padovani+2018 for Wilms+2000 abundances
   real(dp), parameter :: z_pp = 2.17               ! Padovani+2018 for Wilms+2000 abundances
   real(dp), parameter :: z_ion = 2.01              ! Padovani+2018 for Wilms+2000 abundances

   !--------------------------------------------------------------------------------------
   ! Shocked ISM (Forward Shock)
   real(dp), parameter :: n_sh_FS = 4.*n_ISM                ! number density (forward shock)
   real(dp), parameter :: T_sh_FS = 1.36d7*(v_ISM/1.d8)**2  ! temperature (FS)
   real(dp), parameter :: c_s_FS = v_ISM*cte_cs             ! velocity of the shocked fluid (FS)

   !--------------------------------------------------------------------------------------
   ! Heated dust (IR)   
   real(dp), parameter :: T_IR = 98.d0           ! K2017   
   !real(dp), parameter :: T_IR = 27.*(0.2)**(-1./6.)*(L_st/1.d38)**(1./6.)*(R0/pc)**(-1./3.) ! Eq.10 del Valle 2012
   real(dp), parameter :: T_IR_ad = k*T_IR/mec2  ! Adimensional dust temperature
   real(dp), parameter :: L_IR = L_st/350.d0     ! Dust luminosity (K2017, Struck 2020)
   real(dp), parameter :: R_IR = sqrt( L_IR/(4.*pi*sigma*T_IR**4) ) ! cavity -> 1/4*pi, not 16*pi

   ! For IC (both diste.f90 and radic.f90)
   real(dp), parameter :: Kappa_st = (R_st/(2.*R0))**2
   real(dp), parameter :: Kappa_IR = ( R_IR / R0 )**2  
   real(dp), parameter :: Kappa_CMB = 1.d0             

   !--------------------------------------------------------------------------------------
   ! Non-thermal particles
   real(dp), parameter :: efi_ac_RS = 2.d0*pi*(c/v_w)**2    ! Particle accel. efficiency at the RS
   real(dp), parameter :: efi_ac_FS = 2.d0*pi*(c/v_ISM)**2  ! Particle accel. efficiency at the FS
   real(dp), parameter :: alpha_NT = 2.d0            ! Injection spectral index
   real(dp), parameter :: R_acel = R0  !

   real(dp), parameter :: frac_NT = 0.1_dp             ! Fraction of L_inj_perp converted to L_NT
   real(dp), parameter :: inj_e = frac_NT*0.05d0       ! Fraction L_NT_e/L_NT
   real(dp), parameter :: inj_p = frac_NT*0.95d0       ! Fraction L_NT_p/L_NT 

   !real(dp), parameter :: B = frac_B*sqrt(4.d0*pi*rho_w*4.*v_w**2) 
   ! If U_mag = frac*U_th -> B^2/4pi = frac*(g_ad/(g_ad - 1))*P 
   ! If P_mag = frac*P_th -> B^2/8pi = frac*P 
   real(dp), parameter :: frac_B = sqrt(0.1_dp)        ! Fraction B/B_equipartition
   real(dp), parameter :: B_RS = frac_B*sqrt(Mdot*v_w/R0**2) ! del Valle 2012
   real(dp), parameter :: B_FS = frac_B*sqrt(6.d0*pi*rho_ISM*v_ISM**2) ! P ~ 3/4*rho*v^2


   !============================================================================

END module system_parameters
