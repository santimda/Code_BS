module distp

   private 
   public :: distp_run

contains

   !============================================================================
   ! distp_run
   !============================================================================
   ! Calculate the steady-state proton energy distribution as: 
   ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE', with |dE/dt|= E/t_loss.
   ! Including escape as an energy loss timescale is suitable por alpha_inj ~ 2, 
   ! so that t_loss = min(t_rad,t_esc)
   !============================================================================

   subroutine distp_run(shock)   
      use global
      implicit none

      character(LEN=*), intent(in) :: shock
      real(dp) :: L_inj_p,integral,Ka_p,Emax_p,log10Emax
      real(dp) :: y, y1, y2
      real(dp) :: Q,S,t_adv,t_pp,t_diff,t_ac,t_loss,t_esc,t_adi,t_cool_p
      real(dp) :: b_pp,sigma_pp,L1,E_punto,c_s,E_NT,D_Bohm,b_ion,t_ion
      real(dp) :: dE2_p,E2_p,Eint2_p,Emax2_p,Emin2_p,efi_ac,cte_ac
      real(dp) :: n_a,L_inj,B,B_par,cte_Bohm
      real(dp), parameter :: E_th=1.22d9*eV, TeV=1.d12*eV, kappa_pp=0.5d0
      integer :: mE2

      interface
         subroutine Biseccion_p(logEmax,t_pp,t_adv,cte_ac,cB)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp), intent(in) :: t_pp,t_adv,cte_ac,cB
            real(dp), intent(out) :: logEmax
         end subroutine Biseccion_p
      end interface

      allocate( N_p(mE_p),E_p(mE_p) )
      810  format(A11,X,Es10.2)

      write(*,*) '\n Calculate the non-thermal proton distribution in the SW'

      !=======================================================================
      ! Define and check some values first

      if (shock=='RS') then        
         n_a = n_sh_RS
         L_inj = L_inj_RS
         c_s = c_s_RS
         efi_ac = efi_ac_RS
         B = B_RS
         else if (shock=='FS') then     
         n_a = n_sh_FS
         L_inj = L_inj_FS
         c_s = c_s_FS
         efi_ac = efi_ac_FS
         B = B_FS
      end if
      B_par = B*par 

      L_inj_p = inj_p*L_inj  
      write(*,810) 'L_injp =',L_inj_p

      !=======================================================================
      ! Define some relevant constants and timescales
      
      cte_ac = efi_ac/(B_par*c*qe)
      cte_Bohm = c/(3.d0*qe*B)

      t_adv = R0/c_s 
      t_adi = 3.d0*R0/c_s*5.d0! adiabatic losses: dE/dt = 1/3*d(log(rho))/dl *E*v, delta_rho/rho ~ 0.2 
      t_pp = 1.d15/n_a/z_pp   ! approx expression for p-p losses, used only to calculate Emax_p

      if (screened) then   ! interactions with neutral matter
         b_ion = z_ion*1.82d-7*n_a*eV  !Schlickeiser 2002 (Eq. 5.3.51), valid for ~2 GeV << E_p < 1 TeV
         else              ! interactions with fully ionized plasma
         b_ion = z_ion*3.1d-7*n_a*eV   !Schlickeiser 2002 (Eq. 5.3.40)
      endif

      ! Search for maximum proton energy
      call Biseccion_p(log10Emax,t_pp,t_adv,cte_ac,cte_Bohm)
      Emax_p = MIN(10.d0**log10Emax,qe*B_par*R_acel)! Hillas criterion

      ! Build a logarithmic-spaced energy vector
      call vector_log(Emin_p,10.d0*Emax_p,Eint_p,E_p)

      ! Q = Ka*E**-alpha*e^(-E/Emax), with Ka given by: int Q(E')*E'dE' = L_inj. 
      ! [Linj]=erg/s, [Q]=1/erg/s. The integral is solved with the trapezoid method
      integral = 0.d0
      y1 = E_p(1)**(-alpha_NT+1.d0)*exp(-E_p(1)/Emax_p)
      do i = 1, mE_p-1
         y2 = E_p(i+1)**(-alpha_NT+1.d0)*exp(-E_p(i+1)/Emax_p)
         y = (y1 + y2)/2.d0
         dE_p = E_p(i)*(Eint_p-1.d0)
         integral = integral + y*dE_p
         y1 = y2
      end do
      Ka_p = L_inj_p/integral

      !=======================================================================
      ! Calculate the proton steady-state energy distribution

      open(70,file='tcoolp_'//shock//'.dat')
      write(70,*) '#E_p  t_ac  t_pp  t_adv  t_diff  t_cool_p  t_esc  t_loss  t_adi  t_ion'

      ! E_p is the energy of the particle energy distribution and E2_e is an integration variable
      mE2 = mE_p                 ! Initial number of steps in E2
      do i =1,mE_p 

         dE_p = E_p(i)*(Eint_p - 1.d0)  
         S = 0.d0     ! S(E) is the result of the integral: int_E^Emax Q(E')dE'/Ka
         Emin2_p = E_p(i)
         Emax2_p = E_p(mE_p)
         E2_p = Emin2_p
         Eint2_p = (10.d0*Emax2_p/Emin2_p)**(1.d0/dfloat(mE2-1))
         do j=1,mE2
            dE2_p = E2_p * (Eint2_p - 1.d0)
            Q = E2_p**(-alpha_NT)*exp(-E2_p/Emax_p) !normalize later
            S = S + Q*dE2_p
            E2_p = E2_p*Eint2_p
         end do
         mE2 = MAX(mE2-1,32) ! make it faster (less steps when E2 closer to Emax_p)

         ! calculate t_pp       
         if (E_p(i) < E_th) then !E_th= 1.22 GeV
            sigma_pp = 0.d0
            b_pp = 1.d-30        ! small number to avoid dividing by zero
         else            
            L1 = log(E_p(i)/TeV)    ! L = ln(Ep/1TeV)
            sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 
            sigma_pp = sigma_pp*(1.d0-(E_th/E_p(i))**4)**2  ! Factor to improve at low energies
            b_pp = c*kappa_pp*n_a*sigma_pp*E_p(i)*z_pp      ! z_pp from Padovani+2018
         endif
         t_pp = E_p(i)/b_pp
         t_ion = E_p(i)/b_ion
         t_cool_p = 1./(1./t_pp + 1./t_adi)

         D_Bohm = E_p(i)*cte_Bohm
         t_diff = R_acel**2/(2.*D_Bohm)   ! del Valle 2012
         t_esc = 1./( 1./t_adv + 1.d0/t_diff ) ! escape

         t_loss = min(t_cool_p,t_esc)     ! general losses

         t_ac = cte_ac*E_p(i)

         E_punto = abs(E_p(i)/t_loss) 
         N_p(i) = Ka_p*S/E_punto          ! N(E)=1/|dE/dt|*S, with |dE/dt|= E/t_loss

         ! Save the relevant timescales for plotting
         write(70,*) E_p(i)/eV,t_ac,t_pp,t_adv,t_diff,t_cool_p,t_esc,t_loss,t_adi,t_ion

      end do

      !===================================================
      !   Save the distribution 

      E_NT = 0.d0
      open(79,file='dist_p_'//shock//'.dat')
      do i=1,mE_p
         dE_p = E_p(i)*(Eint_p-1.d0)
         E_NT = E_NT + E_p(i)*N_p(i)*dE_p
         write(79,*) E_p(i)/eV,N_p(i) 
      end do
      close(79)
      !write(*,810) 'E_NT(p)=', E_NT

      deallocate( N_p,E_p ) 

   END subroutine distp_run



END module distp



!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                            SUBROUTINES                                  //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine Biseccion_p(x3,t_pp,t_adv,c_ac,cB)
   use global
   implicit none
   integer, parameter :: t_err = 4
   real(dp), intent(in) :: t_pp,t_adv,c_ac,cB
   real(dp), intent(out) :: x3
   real(dp) :: a_,b_,c_,err

   interface
      function resta_p(a,t1,t2,c_ac,cB)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: resta_p
         real(dp), intent(in) :: a,t1,t2,c_ac,cB
      end function resta_p
   end interface

   !     El error en el n-esimo será menor a |a-b|/2**n
   a_ = log10(1.d10*eV) ! log10(Emin [CGS]) ~ 10 GeV
   b_ = log10(1.d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
   err = 10.d0**(-t_err) ! error menor que 10**-t
   i = 0
   do while ( abs(b_ - a_) > err )
      i = i+1
      c_ = (a_+b_)/2.
      if ( abs(resta_p(c_,t_pp,t_adv,c_ac,cB)) <= 1.d-8 ) then
         b_ = c_
         else if (resta_p(a_,t_pp,t_adv,c_ac,cB)*resta_p(c_,t_pp,t_adv,c_ac,cB) < 0) then
         b_ = c_
         else if (resta_p(a_,t_pp,t_adv,c_ac,cB)*resta_p(c_,t_pp,t_adv,c_ac,cB) > 0) then
         a_ = c_
      endif
   end do
   x3=b_
   return
end Subroutine Biseccion_p

!/////////////////////////////////////////////////////////////////////////

Function resta_p(log10E,t_pp,t_adv,c_ac,cB)
   use global
   implicit none
   real(dp), intent(in) :: log10E,t_pp,t_adv,c_ac,cB
   real(dp) :: resta_p,E_l,t_diff,t_ac,t_coolp
   real(dp) :: D_Bohm,t_loss,t_esc

   E_l = 10.d0**log10E

   D_Bohm = E_l*cB
   t_diff = R_acel**2/(2.*D_Bohm) ! del Valle 2012
   t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

   t_coolp = t_pp 

   t_loss = min(t_coolp,t_esc) ! general losses

   t_ac = c_ac*E_l
   resta_p = t_loss - t_ac

   return

end Function resta_p

!/////////////////////////////////////////////////////////////////////////

