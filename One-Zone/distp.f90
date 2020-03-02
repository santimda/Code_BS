module distp
  
  use global
  
  implicit none
  private 
  public :: distp_run
  
contains
  
  !============================================================================
  ! distp_run
  !============================================================================
  
  subroutine distp_run  
    use global
    implicit none

    real(dp) :: L_inj_p,L_wind,integral,Ka_p,Emax_p,log10Emax
    real(dp) :: Q,S,t_adv,t_pp,t_diff,t_ac,t_loss,t_esc,E_NT,D_Bohm
    real(dp) :: b_pp,sigma_pp,L1,E_punto,t_cool_p, c_s
    real(dp) :: dE2_p,E2_p,Eint2_p,Emax2_p,Emin2_p
    real(dp), parameter :: E_th=1.22d9*eV, TeV=1.d12*eV
    integer :: mE2
    
    interface
       subroutine Biseccion_p(logEmax,t_pp,t_adv)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: t_pp,t_adv
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion_p
    end interface
    
    allocate ( N_p(mE_p),E_p(mE_p) )
80  format(A8,X,Es10.2)
    
    write(*,*) '\n Calculate the non-thermal proton distribution in the SW'

    L_wind = (Mdot*v_w**2)/2.d0
    L_inj_p = inj_p*L_wind*0.1d0 ! L_inj_wind ~ 10% L_wind

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                    C
    !                            PROTONS                                 C
    !                                                                    C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
    c_s = v_w/sqrt(8.d0)  ! sound speed ~ sqrt(g_ad*P_sw/rho_sw), P_sw ~ 3/4 rho_w*v_w^2
    t_adv = R0/c_s 
    t_pp = 1.d15/n_a      ! approx expression for p-p losses, used only to calculate Emax_p

    call Biseccion_p(log10Emax,t_pp,t_adv)
    Emax_p = MIN(10.d0**log10Emax,qe*B_par*R_acel)! Hillas criterion

    ! Build a logarithmic-spaced energy vector
    call vector_log(Emin_p,10.d0*Emax_p,Eint_p,E_p)

    ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculate Ka from: 
    ! int Q(E')*E'dE' = Linj. [Lrel]=erg/s, [Ka]=erg/s, [Q]=erg/erg**2/s
    
    integral=0.d0
    do i=1,mE_p
       dE_p = E_p(i)*(Eint_p-1.d0)
       integral = integral + E_p(i)**(-alpha_sw+1.d0)*exp(-E_p(i)/Emax_p)*dE_p
    enddo
    
    Ka_p = L_inj_p/integral
    
    ! Calculate the proton steady-state energy distribution
    
    open(70,file='tcoolp.dat')

    !     E_e es la de los electrones; E2_e es la variable de integración.
    mE2 = mE_p                 ! Initial number of steps in E2
    do i =1,mE_p 

       dE_p = E_p(i)*(Eint_p - 1.d0)  
       S = 0.d0
       Emin2_p = E_p(i)
       Emax2_p = E_p(mE_p)
       E2_p = Emin2_p
       Eint2_p = (10.d0*Emax2_p/Emin2_p)**(1.d0/dfloat(mE2-1))
       do j=1,mE2
          dE2_p = E2_p * (Eint2_p - 1.d0)
          Q = E2_p**(-alpha_sw)*exp(-E2_p/Emax_p) !normalize later
          S = S + Q*dE2_p
          E2_p = E2_p*Eint2_p
       enddo

       if (mE2 >= 32) then  ! so it is faster
          mE2 = mE2-1
       endif
       ! S is the result of the integral in [erg/erg/s]. We calculate: 
       ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE', with |dE/dt|= E/t_cool 
       ! in [erg/s], so that [N(E)]=erg/erg**2
       
       if (E_p(i) < E_th) then !E_th= 1.22 GeV
          sigma_pp = 0.d0
          b_pp = 1.d-30         ! invento para no dividir por cero 
       else            
          L1 = log(E_p(i)/TeV)    ! L = ln(Ep/1TeV)
          sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
          ! Agregamos un factor para que ajuste mejor a bajas energías
          sigma_pp = sigma_pp*(1.d0-(E_th/E_p(i))**4)**2
          b_pp = c*0.5*n_a*sigma_pp*E_p(i)
       endif
       t_pp = E_p(i)/b_pp
       t_cool_p = t_pp !1./(1.d0/t_pp + 1.d0/t_adv + 1.d0/t_diff) 

       D_Bohm = E_p(i)*c/(3.d0*qe*B)
       t_diff = H0**2./(2.*D_Bohm) ! del Valle 2012
       t_esc = 1./( 1./t_adv + 1.d0/t_diff ) ! escape
       
       t_loss = min(t_cool_p,t_esc) ! general losses
       
       t_ac = efi_ac_sw*E_p(i)/(B_par*c*qe)

       E_punto = abs(E_p(i)/t_loss) ! NOT SURE
       N_p(i) = Ka_p*S/E_punto

       write(70,*) E_p(i)/eV,t_ac,t_pp,t_adv,t_diff,t_cool_p,t_esc,t_loss

    enddo

    !===================================================
    !   Save the distribution 

 
    E_NT = 0.d0
    open(79,file='dist_p.dat')
    do i=1,mE_p
       dE_p = E_p(i)*(Eint_p-1.d0)
       E_NT = E_NT + E_p(i)*N_p(i)*dE_p
       write(79,*) E_p(i)/eV,N_p(i) 
    enddo
    close(79)
    write(*,80) 'E_NT(p)=', E_NT

    deallocate ( N_p,E_p ) 
 
  END subroutine distp_run
  


END module distp




    !///////////////////////////////////////////////////////////////////////////
    !                                                                         //
    !                            SUBRUTINAS                                   //
    !                                                                         //
    !///////////////////////////////////////////////////////////////////////////

    Subroutine Biseccion_p(x3,t_pp,t_adv)
      use global
      implicit none
      integer t_err
      real(dp), intent(in) :: t_pp,t_adv
      real(dp), intent(out) :: x3
      real(dp) :: a_,b_,c_,err

      interface
         function resta_p(a,t1,t2)
           use iso_fortran_env, only: dp => real64
           implicit none
           real(dp) :: resta_p
           real(dp), intent(in) :: a,t1,t2
         end function resta_p
      end interface
      
      !     El error en el n-esimo será menor a |a-b|/2**n
      a_ = log10(1.d10*eV) ! log10(Emin [CGS]) ~ 10 GeV
      b_ = log10(1.d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
      t_err = 4
      err = 10.d0**(-t_err) ! error menor que 10**-t
      i = 0
      do while ( abs(b_ - a_) > err )
         i = i+1
         c_ = (a_+b_)/2.
         if (resta_p(c_,t_pp,t_adv)==0) then
            b_ = c_
         else if (resta_p(a_,t_pp,t_adv)*resta_p(c_,t_pp,t_adv) < 0) then
            b_ = c_
         else if (resta_p(a_,t_pp,t_adv)*resta_p(c_,t_pp,t_adv) > 0) then
            a_ = c_
         endif
      enddo
      x3=b_
      return
    end Subroutine Biseccion_p
    
    !/////////////////////////////////////////////////////////////////////////
    
    Function resta_p(log10E,t_pp,t_adv)
      use global
      implicit none
      real(dp), intent(in) :: log10E,t_pp,t_adv
      real(dp) :: resta_p,E_l,t_diff,t_ac,t_coolp
      real(dp) :: D_Bohm,t_loss,t_esc

      E_l = 10.d0**log10E

      D_Bohm = E_l*c/(3.d0*qe*B)
      t_diff = H0**2./(2.*D_Bohm) ! del Valle 2012
      t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

      t_coolp = t_pp !1./(1.d0/tloss + 1.d0/t_diff)
      
      t_loss = min(t_coolp,t_esc) ! general losses

      t_ac = efi_ac_sw*E_l/(B_par*c*qe)
      resta_p = t_loss - t_ac

      return

    end Function resta_p
    
    !/////////////////////////////////////////////////////////////////////////

