module CD_diste

   private
   public :: CD_diste_run 

   contains

   !============================================================================
   ! CD_diste_run (START)
   !============================================================================
   !  This program calculates the electron energy distribution at each position
   !  of the CD using CD.dat and termo_sw.dat. The program generates diste_sw.dat, 
   !  in which Ne(l,E(i)) is saved, with m_theta*mE_e data; diste_sw_total.dat, 
   !  which stores Ne(E(i)), and ENTe.dat, with E_NT(l). Emin_e and Emax_e are 
   !  fixed in global.
   !============================================================================

   subroutine CD_diste_run  

      use global
      implicit none
      integer :: l_min,l2
      real(dp) :: t_cell(m_theta),t_adv(m_theta),t_coole(m_theta,mE_e),t_loss(m_theta,mE_e)
      real(dp) :: t_ic_st,t_ic_IR,t_ic,t_sy,t_ac,t_br,D_diff,t_diff,t_esc,t_adi(m_theta)
      real(dp) :: cte_sy, gamm, cte_IC, C_ic_st,kC_ic_st, C_ic_IR(m_theta),kC_ic_IR
      real(dp) :: log_coul,gdot_coul,cte_coul,t_coul
      real(dp) :: log10Emaxl,R_acel,Q,integral
      real(dp) :: E_dot(m_theta,mE_e),Emax_e0(m_theta),Ka,E_NT
      real(dp) :: N_e0(mE_e),n_a(m_theta),L_inj(m_theta),dR(m_theta)
      real(dp) :: L_inj_e(m_theta), sum,dummy
      real(dp) :: dE_aprox, cool
      integer :: index1,index2,count,steps
      logical :: new_count

      interface
      subroutine Biseccion(c_sy,Bl,logEmax,c1,c2,T1,tbr,tad,xi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: c_sy,Bl,c1,c2,T1,tbr,tad,xi
         real(dp), intent(out) :: logEmax
      end subroutine Biseccion
      end interface

      allocate( theta(m_theta),R_theta(m_theta),H_theta(m_theta) )
      allocate( P(m_theta),rho(m_theta),T(m_theta),B(m_theta),v_t(m_theta) )
      allocate( N_e(mE_e),E_e(mE_e),B_par(m_theta),N_etot(m_theta,mE_e) )
      allocate( T_IR(m_theta),T_IR_ad(m_theta) )
      allocate( efi_ac_sw(m_theta) )

      write(*,*) 'Calculate the non-thermal electron distribution in the SW'

      !---------------------------------------------------------------------------------
      !  Read the position of the shock and cell size: theta, r_norm(theta), dR(theta); 
      !  the thermodynamical quantities and injected luminosity at each cell

      open(12,file='Linj_sw.dat')
      do l=1,m_theta
         read(12,*) theta(l),R_theta(l),L_inj(l),dR(l)
      end do
      close(12)
      R_theta = R0*R_theta     ! The file contained r_norm, not R_theta
      dR = R0*dR               ! The file contained dR/R0, not dR [cm]
      L_inj_e = inj_e*L_inj    ! The cell size is already taken into account in L_inj

      open(11,file='termo_sw.dat') 
      do l=1,m_theta
         read(11,*) l2,rho(l),P(l),v_t(l),T(l),B(l) ! all in CGS units
      end do
      close (11)
      B_par = B*par        ! For synchrotron, only the parallel component of B is relevant 
      n_a = rho/mp/mu_sw   ! Ambient number density

      open(103,file='Hnorm_sw.dat') 
      do l=1,m_theta
         read(103,*) dummy,H_theta(l) ! normalized
         H_theta(l) = H_theta(l)*R_theta(l)
      end do
      close(103)

      open(113,file='T_IR.dat') 
      do l=1,m_theta
         !       T_IR(l) = 27.*(0.2)**(-1./6.)*(L_st/1.d38)**(1./6.)*(R_theta(l)/pc)**(-1./3.) ! Eq.10 del Valle 2012
         T_IR(l) =  100.d0            ! K2017
         T_IR_ad(l) = k*T_IR(l)/mec2  ! Adimensional dust temperature
         Kappa_IR(l) = L_IR/(4.*pi*sigma*T_IR(l)**4*R_theta(l)**2) ! 4*pi, not 16*pi, for a cavity
         write(113,*) theta(l),T_IR(l),T_IR_ad(l),Kappa_IR(l) ! normalized
      end do
      close(113)

      open(114,file='efi_ac_sw.dat')
      do l=1,m_theta
         read(114,*) dummy,efi_ac_sw(l)
      end do
      close(114)

      !--------------------------------------------------------------------
      ! Build an energy vector that is common for all linear-emitters: 

      call vector_log(Emin_e,Emax_e,Eint_e,E_e)

      open(70,file='CD_tcoole.dat')
      cte_IC = pi*h_bar**3*c**2/(2.*re**2*mec2**3)
      C_ic_st = (cte_IC)/(T_st_ad**2)
      do l=1,m_theta
         C_ic_IR(l) = (cte_IC)/(T_IR_ad(l)**2)
      end do

      !========================================================================
      ! First calculate the cooling time at all energies and at each position 
      ! (and also the maximum energy)

      do l = 1, m_theta
         t_cell(l) = dR(l)/v_t(l)      ! Cell advection time
         t_adv(l) = R_theta(l)/v_t(l)  ! Characteristic advection time
         if (l < m_theta) then
            t_adi(l) = 3.*dR(l)/(v_t(l)*log(rho(l)/rho(l+1))) ! Characteristic adiabatic losses time
         else if (l == m_theta) then
            t_adi(l) = t_adi(l-1)
         end if
      end do

      do l=1,m_theta
         Kappa_st = ( R_st/(2.*R_theta(l)) )**2
         kC_ic_st = C_ic_st/Kappa_st
         kC_ic_IR = C_ic_IR(l)/Kappa_IR(l)
         R_acel = H_theta(l)           ! Accelerator size  ????????????
         t_br = 1.d15/n_a(l)/z_br      ! added z_be according to Padovani+2018
         cte_sy = 6.*pi*mec2**2/(c*sigma_T*B(l)**2)   ! Valid for isotropic e- distribution

         call Biseccion(cte_sy,B(l),log10Emaxl,kC_ic_st,kC_ic_IR,T_IR_ad(l),t_br,t_adv(l),efi_ac_sw(l))
         Emax_e0(l) = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efi_ac_sw(l))) 
         if (Emax_e0(l) >= 0.5*Emax_e) write(*,*) 'Emax_e in global.f90 is too small!'

         log_coul = 29.7d0 + log( (T(l)/1.d6)/sqrt(n_a(l)) ) 
         cte_coul = 3.d-14*n_a(l)*log_coul

         do i=1,mE_e

            D_diff = E_e(i)*c/(3.d0*qe*B(l)) ! Bohm
            t_diff = H_theta(l)**2/(2.d0*D_diff) ! Diffusion losses
            t_esc = 1./( 1./t_adv(l) + 1.d0/t_diff ) ! escape

            t_sy = cte_sy/E_e(i) 
            gamm = E_e(i)/mec2
            t_ic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
            t_ic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad(l))
            t_ic = 1./(1./t_ic_st + 1./t_ic_IR)

            gdot_coul = cte_coul/sqrt(gamm**2-1.d0) ! Pittard2006
            t_coul = gamm/gdot_coul

            t_coole(l,i) = 1./(1./t_sy + 1./t_ic + 1./t_br + 1./t_adi(l) + 1./t_coul)

            t_loss(l,i) = 1.d0/(1.d0/t_coole(l,i) + 1.d0/t_esc) ! general losses

            t_ac = efi_ac_sw(l)*E_e(i)/(B(l)*c*qe)

            write(70,*) l,E_e(i)/eV,t_ic,t_sy,t_br,t_diff,t_ac,t_coole(l,i),t_adv(l),t_cell(l),t_ic_st,t_ic_IR,t_esc,t_loss(l,i),t_adi(l),t_coul

            E_dot(l,i) = E_e(i)/t_coole(l,i)

         end do

      end do
      close(70)

      !========================================================================

      ! Consider different starting points for each linear emitter

      N_etot = 0.d0

      ! At the starting point of the line we calculate the distribution 
      ! Ne_0(E_0), the cooling rate Edot0, and we sum it to Ne_tot(l,i)

      Do l_min=1,m_theta

         ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
         ! by the condition: int Q(E')*E'dE' = Linj. [Lrel]=erg/seg, [Ka]=erg/seg, [Q]=erg/erg**2/seg

         integral=0.d0
         i = 1
         do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0(l_min))
            dE_e = E_e(i)*(Eint_e-1.d0)
            integral = integral + E_e(i)**(-alpha_sw+1.d0)*exp(-E_e(i)/Emax_e0(l_min))*dE_e
            i=i+1
         end do
         Ka = L_inj_e(l_min)/integral ! Injected at the starting point of the line

         !------------------------------------------------------------------------

         N_e0 = 0.d0            ! Initialize to zero at each loop in l_min

         i=1

         do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0(l_min)) 
            ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
            ! in [erg/seg], dt = min(t_cell,t_cool) so that [N(E)]=erg/erg**2

            Q = Ka*E_e(i)**(-alpha_sw)*exp(-E_e(i)/Emax_e0(l_min))
            N_e0(i) = Q*min(t_cell(l_min),t_coole(l_min,i))
            N_etot(l_min,i) = N_etot(l_min,i)+N_e0(i)  

            ! We calculate the evolved version of the injected electron distribution at each
            ! position "l" as: N_l(l) = N_l-1* (Edot_l-1/Edot_l), where the energy also evolves
            ! as E0(l) given by: t_cell(l) = int_E0^El Edot * dE

            index1 = i            ! For the cell "l"
            index2 = i            ! For the previous cell, "l-1"
            N_e = 0.d0            ! Initialize to zero at each loop in l_min
            N_e(index1) = N_e0(i) ! This will be evolved
            l = l_min + 1  
            count = 0
            cool = 0.d0
            new_count = .true.     

            do while( l <= m_theta .and. index1 > 1) ! if index1 = 0, the e- are completely cooled
               ! At the position l, the particle has E_e(index1) instead of E_e(i)
               ! We find index1 <= i according to the cooling time:

               if (l == l_min+1) then
                  do l2=l,m_theta
                     cool = cool + t_cell(l2)/t_coole(l2,index1)
                  end do
               else
                  cool = cool - t_cell(l-1)/t_coole(l-1,index1)
               end if
               dE_aprox = (1.-t_cell(l)/t_coole(l,index1))

               if ( 1.-cool > 0.99/Eint_e ) then ! negligible cooling
                  continue

               elseif ( 1.-cool <= 0.99/Eint_e .and. dE_aprox > 1./Eint_e ) then
                  if (new_count) then
                     sum = 1.d0
                     steps = 0 ! number of steps until it changes E(index1) -> E(index-1)
                     l2 = l
                     do while (sum > 1./Eint_e .and. l2 <= m_theta)
                        steps = steps + 1
                        sum = sum*(1.-t_cell(l2)/t_coole(l2,index1))
                        l2 = l2 + 1
                        if (sum > 1./Eint_e .and. l2 == m_theta) steps = m_theta ! NOT COOLING
                     end do
                     new_count = .false.
                  end if
                  count = count + 1

                  if (count >= steps) then
                     count = 0
                     index1 = index1 - 1
                     new_count = .true. 
                  end if

               elseif ( dE_aprox <= 1./Eint_e ) then
                  sum = 0.d0
                  do while (sum < t_cell(l) .and. index1 > 1)
                     index1 = index1 - 1    
                     dE_e = E_e(index1)*(Eint_e-1.d0)
                     sum = sum + dE_e/E_dot(l,index1)        
                  end do
               end if

               N_e(index1) = N_e(index2) * ( E_dot(l-1,index2)/E_dot(l,index1) ) 

               N_etot(l,index1) = N_etot(l,index1)+N_e(index1)  ! The energies match by construction

               index2 = index1  ! index2 keeps the value of the previous step

               l = l + 1 
            end do              ! End loop in shock position (l) 
            i = i + 1
         end do                 ! End loop in electron starting energies (i)

      end do                    ! End loop in starting position (l_min = # linear-emitters)

      !----------------------------------------------------------------------
      ! Save the electron distribution at each position 
      open(79,file='diste_sw.dat')
      do l=1,m_theta
         do i=1,mE_e     
            write(79,*) theta(l),E_e(i)/eV,N_etot(l,i)
         end do
         write(79,*) 
      end do
      close(79)

      sum = 0.d0
      open(80,file='ENTe.dat')
      do l=1,m_theta
         E_NT = 0.d0
         do i=1,mE_e
            dE_e = E_e(i)*(Eint_e - 1.d0)  
            E_NT = E_NT + E_e(i)*N_etot(l,i)*dE_e
         end do
         E_NT = E_NT * 4.d0*mphi ! Total NT energy in each BS "ring"
         write(80,*) l, E_NT
         sum = sum + E_NT
      end do
      close(80)
      !177  format(A10,X,Es10.2)
      !write(*,177) 'E_tot_e=', sum ! Total NT energy in the BS

      !------------------------------------------------------------
      ! Group in bins according to position angle
      call Ne_bin(E_e,N_etot,theta)


      deallocate( theta,R_theta,H_theta,B,B_par,rho,P,v_t,T )
      deallocate( N_e,E_e,N_etot ) 
      deallocate( T_IR,T_IR_ad, efi_ac_sw )

   END subroutine CD_diste_run

   !============================================================================
   ! CD_diste_run (END)
   !============================================================================

END module CD_diste

!--------------------------------------------------------------------------!
!                                                                          !
!                            SUBROUTINES                                   !
!                                                                          !
!--------------------------------------------------------------------------!

subroutine Ne_bin(Ee,Ne_tot,ang)

   use global
   implicit none    
   real(dp), intent(in) :: Ee(mE_e),Ne_tot(m_theta,mE_e),ang(m_theta)
   integer, parameter :: N_bins = 5 ! For the plot diste_sw_bins
   real(dp) :: N_bin_i(N_bins,mE_e), N_total_i

   N_bin_i = 0.d0
   open(81,file='diste_sw_total.dat')
   do i=1,mE_e
      N_total_i = 0.d0
      do l=1,m_theta
         ! Make some distributions binned in intervals of 0.1*D
         if (ang(l) < 0.15*pi) then
            N_bin_i(1,i) = N_bin_i(1,i) + Ne_tot(l,i)
         elseif (ang(l) >= 0.15*pi .and. ang(l) < 0.30*pi) then
            N_bin_i(2,i) = N_bin_i(2,i) + Ne_tot(l,i)
         elseif (ang(l) >= 0.30*pi .and. ang(l) < 0.45*pi) then
            N_bin_i(3,i) = N_bin_i(3,i) + Ne_tot(l,i)
         elseif (ang(l) >= 0.45*pi .and. ang(l) < 0.60*pi) then
            N_bin_i(4,i) = N_bin_i(4,i) + Ne_tot(l,i)
         elseif (ang(l) >= 0.60*pi .and. ang(l) < 0.75*pi) then
            N_bin_i(5,i) = N_bin_i(5,i) + Ne_tot(l,i)
         end if  
         N_total_i = N_total_i + Ne_tot(l,i)
      end do
      write(81,*) Ee(i)/eV, N_total_i * mphi*1.d0 ! Total distribution at each energy for 1 of the 4 quadrants
   end do
   close(81)

   open(82,file='diste_sw_bin.dat')
   do j = 1,N_bins
      do i = 1, mE_e
         write(82,*) j,Ee(i)/eV, N_bin_i(j,i)* mphi*1.d0
      end do
      write(82,*) 
   end do
   close(82)

end subroutine Ne_bin


!=========================================================================


Subroutine Biseccion(csy,Bl,x3,c1,c2,T1,t_br,t_adv,xi)
   use global
   implicit none
   real(dp) :: a_,b_,c_
   real(dp), intent(in) :: csy,Bl,c1,c2,T1,t_br,t_adv,xi
   real(dp), intent(out) :: x3
   real(dp), parameter :: err = 1.d-3  ! error in Emax < err
   interface
   function resta(a,csy,Bl,ic1,ic2,T1,t_br,t_adv,efi)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp) :: resta
      real(dp), intent(in) :: a,csy,Bl,ic1,ic2,T1,t_br,t_adv,efi
   end function resta
   end interface

   ! Error in nth step < |a-b|/2**n
   a_ = log10(2.d10*eV)    ! log10(Emin [CGS]) ~ 20 GeV
   b_ = log10(Emax_e)   
   i = 0
   do while ( abs(b_ - a_) > err )
      i = i+1
      c_ = (a_+b_)/2.
      if (resta(c_,csy,Bl,c1,c2,T1,t_br,t_adv,xi) < 1.d-8 ) then 
         b_ = c_
      elseif (resta(a_,csy,Bl,c1,c2,T1,t_br,t_adv,xi)*resta(c_,csy,Bl,c1,c2,T1,t_br,t_adv,xi) < 0) then 
         b_ = c_
      elseif (resta(a_,csy,Bl,c1,c2,T1,t_br,t_adv,xi)*resta(c_,csy,Bl,c1,c2,T1,t_br,t_adv,xi) > 0) then 
         a_ = c_          
      end if
   end do
   x3=b_
   return

end Subroutine Biseccion

!=========================================================================

Function resta(log10E,csy,Bl,c1,c2,T1,t_br,t_adv,xi)
   use global
   implicit none
   real(dp), intent(in) :: log10E,csy,Bl,c1,c2,T1,t_br,t_adv,xi
   real(dp) :: resta
   real(dp) :: E_l,gamm,t_sy,t_ic_st,t_ic_IR,t_ac,t_coole
   real(dp) :: D_diff,t_diff,t_loss,t_esc

   E_l = 10.d0**log10E

   D_diff = E_l*c/(3.d0*qe*Bl)         ! Bohm
   t_diff = (0.2*R0)**2/(2.d0*D_diff)  ! Diffusion losses APPROX
   t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

   gamm = E_l/mec2
   t_sy = csy/E_l
   t_ic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
   t_ic_IR = c2*gamm/F_iso(4.*gamm*T1)  
   t_coole = 1./(1.d0/t_sy + 1.d0/t_ic_st + 1.d0/t_ic_IR + 1.d0/t_br) ! radiative cooling

   t_loss = 1.d0/(1.d0/t_coole + 1.d0/t_esc) ! general losses

   t_ac = xi*E_l/(Bl*c*qe)  !      t_ac = (efiac*B_l*c*qe/E_l)**(-1)
   resta = t_loss - t_ac

   return
end Function resta

!=========================================================================

