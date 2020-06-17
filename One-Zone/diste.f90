module diste

   private 
   public :: diste_run

contains

   !============================================================================
   ! diste_run
   !============================================================================
   ! Calculate the steady-state electron energy distribution as: 
   ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE', with |dE/dt|= E/t_loss.
   ! Including escape as an energy loss timescale is suitable por alpha_inj ~ 2, 
   ! so that t_loss = min(t_rad,t_esc)
   !============================================================================


   subroutine diste_run(shock)   
      use global
      implicit none
      character(LEN=*), intent(in) :: shock

      integer :: mE2
      real(dp) :: gamm,const_ic,C_ic_st,C_ic_IR,C_ic_CMB
      real(dp) :: t_ic_st,t_ic_IR,t_ic_CMB,t_ic,t_sy,t_ac,t_adv,t_br,t_diff,t_adi
      real(dp) :: log_coul,gdot_coul,t_coul,cte_sy,cte_br,p_br
      real(dp) :: t_coole,t_esc(mE_e),t_loss
      real(dp) :: log10Emax,D_Bohm
      real(dp) :: E2_e,dE2_e,Emin2_e,Emax2_e,Eint2_e,integral,S,Q
      real(dp) :: E_punto,Emax_e,Ka,E_NT,L_inj_e,c_s
      real(dp) :: T_6,n_a,L_inj

      interface
         subroutine Biseccion(logEmax,c1,c2,c3,t_br,t_adv)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp), intent(in) :: c1,c2,c3,t_br,t_adv
            real(dp), intent(out) :: logEmax
         end subroutine Biseccion
      end interface

      allocate( N_e(mE_e),E_e(mE_e) )
      810  format(A15,X,Es10.2)

      write(*,*) '\n Calculate the non-thermal electron distribution in the SW'

      !=======================================================================
      ! Define and check some values first variables

      if (shock=='RS') then        
         n_a = n_sh_RS
         T_6 = T_sh_RS/1.d6
         L_inj = L_inj_RS
         c_s = c_s_RS
         else if (shock=='FS') then     
         n_a = n_sh_FS
         T_6 = T_sh_FS/1.d6
         L_inj = L_inj_FS
         c_s = c_s_FS
      end if

      L_inj_e = inj_e*L_inj 
      write(*,810) 'L_inje =',L_inj_e
      write(*,810) 'U_IR [eV/cm3] =',L_IR/(pi*R0**2*c)/eV
      write(*,810) 'U_st [eV/cm3] =',L_st/(4.*pi*R0**2*c)/eV
      !write(*,810) 'U_IR/U_st =',16.*pi*R0**2*T_IR**4/L_st
      !write(*,810) 'R_IR/R0 =',R_IR/R0
      !write(*,810) 'Kappa_IR =',Kappa_IR
      !write(*,810) 'Kappa_IR =',L_IR/(4.*sigma*T_IR**4*R0**2 ) 

      !=======================================================================
      ! Define some relevant constants and timescales

      cte_sy = 6.*pi*mec2**2/(c*sigma_T*B**2)   ! Valid for isotropic e- distribution
      log_coul = 29.7d0 + log( T_6/sqrt(n_a) ) 

      const_ic = pi*h_bar**3*c**2/(2.*re**2*mec2**3)
      C_ic_st = const_ic/(T_st_ad**2*Kappa_st)
      C_ic_IR = const_ic/(T_IR_ad**2*Kappa_IR)
      C_ic_CMB = const_ic/(T_CMB_ad**2*Kappa_CMB)

      t_adv = R0/c_s 
      t_adi = 3.d0*R0/c_s*5.d0! adiabatic losses: dE/dt = 1/3*d(log(rho))/dl *E*v, delta_rho/rho ~ 0.2 

      cte_br = 4.*n_a*z_br*re**2*alfa*c   
      if (screened) p_br = log(183./z_br**(1./6)) - 1./18. ! for screened nuclei
      t_br = 1.d15/n_a/z_br   ! added z_be according to Padovani+2018; very approximated expression

      ! Search for maximum electron energy
      call Biseccion(log10Emax,C_ic_st,C_ic_IR,C_ic_CMB,t_br,t_adv)
      Emax_e = MIN(10.d0**log10Emax,qe*B_par*R_acel)! Hillas criterion

      ! Build a logarithmic-spaced energy vector
      call vector_log(Emin_e,10.d0*Emax_e,Eint_e,E_e)

      ! Q = Ka*E**-alpha*e^(-E/Emax), with Ka given by: int Q(E')*E'dE' = L_inj
      ! [Linj]=erg/s, [Q]=1/erg/s
      integral=0.d0
      do i=1,mE_e
         dE_e = E_e(i)*(Eint_e-1.d0)
         integral = integral + E_e(i)**(-alpha_NT+1.d0)*exp(-E_e(i)/Emax_e)*dE_e
      end do
      Ka = L_inj_e/integral

      !=======================================================================
      ! Calculate the electron steady-state energy distribution

      open(70,file='tcoole_'//shock//'.dat')
      write(70,*) '#E_e  t_ic  t_sy  t_br  t_adv  t_ac  t_coole  t_diff  t_ic_st  t_ic_IR  t_ic_CMB  t_esc  t_loss  t_coul  t_adi'

      ! E_e is the energy of the particle energy distribution and E2_e is an integration variable
      mE2 = mE_e                 ! Initial number of steps in E2
      do i =1,mE_e 
         dE_e = E_e(i)*(Eint_e - 1.d0)  
         S = 0.d0 ! S(E) is the result of the integral: int_E^Emax Q(E')dE'/Ka
         Emin2_e = E_e(i)
         Emax2_e = E_e(mE_e)
         E2_e = Emin2_e
         Eint2_e = (10.d0*Emax2_e/Emin2_e)**(1.d0/dfloat(mE2))
         do j=1,mE2
            dE2_e = E2_e * (Eint2_e - 1.d0)
            Q = E2_e**(-alpha_NT)*exp(-E2_e/Emax_e) !normalize later
            S = S + Q*dE2_e
            E2_e = E2_e*Eint2_e
         end do
         mE2 = MAX(mE2-1,32) ! make it faster (less steps when E2 closer to Emax_e)

         t_sy = cte_sy/E_e(i)

         if (.NOT. screened) p_br = log(E_e(i)/mec2 + 0.36) ! for unscreened nuclei
         t_br = 1./(cte_br*p_br)

         gamm = E_e(i)/mec2
         t_ic_st = C_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
         t_ic_IR = C_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad)
         t_ic_CMB = C_ic_CMB*gamm/F_iso(4.*gamm*T_CMB_ad)       
         t_ic = 1./(1.d0/t_ic_st + 1.d0/t_ic_IR + 1.d0/t_ic_CMB)      
       
         gdot_coul = 3.d-14*n_a*log_coul/sqrt(gamm**2-1.d0) ! Pittard2006
         t_coul = gamm/gdot_coul*z_ion    ! Not sure, Coulomb losses = ionization losses ?

         t_coole = 1./(1.d0/t_sy + 1.d0/t_ic + 1.d0/t_br + 1.d0/t_adi + 1.d0/t_coul) 

         D_Bohm = E_e(i)*c/(3.*qe*B)
         t_diff = R_acel**2/(2.*D_Bohm) ! del Valle 2012

         t_esc(i) = 1./(1./t_adv + 1.d0/t_diff) ! escape

         t_loss = min(t_esc(i),t_coole)    ! general losses      

         t_ac = efi_ac_sw*E_e(i)/(B_par*c*qe)

         ! N(E)=1/|dE/dt|*S, with |dE/dt|= E/t_loss
         E_punto = abs(E_e(i)/t_loss)
         N_e(i) = Ka*S/E_punto

         ! Save the relevant timescales for plotting
         write(70,*) E_e(i)/eV,t_ic,t_sy,t_br,t_adv,t_ac,t_coole,t_diff,t_ic_st,t_ic_IR,t_ic_CMB,t_esc(i),t_loss,t_coul,t_adi

      end do

      !===================================================
      !   Save the distribution 

      E_NT = 0.d0
      open(79,file='dist_e_'//shock//'.dat') ! Electron energy distribution
      open(69,file='inj_e_'//shock//'.dat')  ! Approx. electron E-distribution neglecting radiative losses
      do i=1,mE_e
         dE_e = E_e(i)*(Eint_e - 1.d0)  
         E_NT = E_NT + E_e(i)*N_e(i)*dE_e  ! Control
         write(79,*) E_e(i)/eV,N_e(i)
         write(69,*) E_e(i)/eV,Ka*E_e(i)**(-alpha_NT)*exp(-E_e(i)/Emax_e)*t_esc(i)
      end do
      !write(*,810) 'E_NT =',E_NT
      close(79)
      close(69)

      deallocate( N_e,E_e ) 

   END subroutine diste_run


END module diste



!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                           SUBROUTINES                                   //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine Biseccion(x3,c1,c2,c3,t_br,t_adv)
   use global
   implicit none
   integer, parameter :: t_err = 4
   real(dp) :: a_,b_,c_,err
   real(dp), intent(in) :: c1,c2,c3,t_br,t_adv
   real(dp), intent(out) :: x3
   interface
      function resta(a,ic1,ic2,ic3,t_br,t_adv)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: resta
         real(dp), intent(in) :: a,ic1,ic2,ic3,t_br,t_adv
      end function resta
   end interface

   ! Error in nth step < |a-b|/2**n
   a_ = log10(1.d9*eV)     ! log10(Emin [CGS]) ~ 1 GeV
   b_ = log10(1.d13*eV)    ! log10(Emax [CGS]) ~ 10 TeV
   err = 10.d0**(-t_err)   ! error in Emax < 10**-t
   i = 0
   do while ( abs(b_ - a_) > err )
      i = i+1
      c_ = (a_+b_)/2.
      if ( abs(resta(c_,c1,c2,c3,t_br,t_adv)) <= 1.d-8 ) then
         b_ = c_
         elseif (resta(a_,c1,c2,c3,t_br,t_adv)*resta(c_,c1,c2,c3,t_br,t_adv) < 0) then
         b_ = c_
         elseif (resta(a_,c1,c2,c3,t_br,t_adv)*resta(c_,c1,c2,c3,t_br,t_adv) > 0) then
         a_ = c_
      endif
   end do
   x3=b_
   return
end Subroutine Biseccion

!/////////////////////////////////////////////////////////////////////////

Function resta(log10E,c1,c2,c3,t_br,t_adv)
   use global
   implicit none
   real(dp), intent(in) :: log10E,c1,c2,c3,t_br,t_adv
   real(dp) :: resta,E_l,gamm
   real(dp) :: t_sy,t_ic_st,t_ic_IR,t_ic_CMB,t_ac,t_coole
   real(dp) :: D_Bohm,t_diff,t_loss,t_esc

   E_l = 10.d0**log10E

   D_Bohm = E_l*c/(3.d0*qe*B)
   t_diff = R_acel**2/(2.*D_Bohm) ! del Valle 2012
   t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

   gamm = E_l/mec2
   t_ic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
   t_ic_IR = c2*gamm/F_iso(4.*gamm*T_IR_ad)    
   t_ic_CMB = c3*gamm/F_iso(4.*gamm*T_CMB_ad)
   t_sy = 1.d0/1.6d-3/B**2/E_l 

   t_coole = 1./(1.d0/t_sy + 1.d0/t_ic_st + 1.d0/t_ic_IR + 1./t_br + 1.d0/t_ic_CMB)

   t_loss = min(t_coole,t_esc) ! general losses

   t_ac = efi_ac_sw*E_l/(B_par*c*qe)
   resta = t_loss - t_ac

   return

end Function resta

!/////////////////////////////////////////////////////////////////////////

