module diste
  
  use global
  
  implicit none
  private
  public :: diste_run
  
contains
  
  !============================================================================
  ! diste_run
  !============================================================================
  
  subroutine diste_run  
    use global
    implicit none

    integer :: mE2
    real(dp) :: gamm, C_ic_st,kC_ic_st, C_ic_IR,kC_ic_IR, C_ic_CMB,kC_ic_CMB
    real(dp) :: t_ic_st,t_ic_IR,t_ic_CMB,t_ic,t_sy,t_ac,t_adv,t_br,t_diff,t_adi
    real(dp) :: t_coole,t_esc(mE_e),t_loss
    real(dp) :: log10Emax,D_Bohm
    real(dp) :: E2_e,dE2_e,Emin2_e,Emax2_e,Eint2_e,integral,S,Q
    real(dp) :: E_punto,Emax_e,Ka,E_NT,L_inj_e,L_wind,c_s

    interface
       subroutine Biseccion(logEmax,c1,c2,c3,t_br,t_adv)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: c1,c2,c3,t_br,t_adv
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion
    end interface
    
    allocate( N_e(mE_e),E_e(mE_e) )
80  format(A7,X,Es10.2)
810  format(A15,X,Es10.2)

    L_wind = (Mdot*v_w**2.)/2.d0
    L_inj_e = iny_e*L_wind*0.5 ! L_wind_perp ~ 50% L_wind

    write(*,810) 'L_inj =',L_wind*0.5
    write(*,810) 'L_inje =',L_inj_e
    write(*,810) 'L_injp =',L_wind*0.5*iny_p
    write(*,810) 'L_st =',L_st
!    write(*,810) 'U_IR [eV/cm3] =',(4.*sigma/c)*T_IR**4/eV
    write(*,810) 'U_IR [eV/cm3] =',L_IR/(pi*R0**2.*c)/eV
    write(*,810) 'U_st [eV/cm3] =',L_st/(4.*pi*R0**2.*c)/eV
    write(*,810) 'U_IR/U_st =',( L_IR / L_st ) !* ( R0 / H0 ) **2.
!    write(*,810) 'U_IR/U_st =',16.*pi*R0**2.*T_IR**4./L_st
!    write(*,810) 'L_IR/L_st =',L_IR/L_st
    write(*,810) 'R_IR/R0 =',R_IR/R0
    write(*,810) 'Kappa_IR =',Kappa_IR
!    write(*,810) 'Kappa_IR =',L_IR/(4.*sigma*T_IR**4.*R0**2. ) 


    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                     C
    !                              ELECTRONS                              C
    !                                                                     C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Search for maximum electron energy
    C_ic_st = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T_st_ad**2.)
    C_ic_IR = ( pi*h_bar**3*c**2/(2.d0*re**2*mec2**3) )/(T_IR_ad**2)
    C_ic_CMB = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T_CMB_ad**2.)
    kC_ic_st = C_ic_st/Kappa_st
    kC_ic_IR = C_ic_IR/Kappa_IR  
    kC_ic_CMB = C_ic_CMB/Kappa_CMB
    
    c_s = v_w/sqrt(8.d0)  ! sound speed ~ sqrt(g_ad*P_sw/rho_sw), P_sw ~ 3/4 rho_w*v_w^2
    t_adv = R0/c_s 
    t_adi = R0/c_s*1.d2   ! adiabatic losses: dE/dt = 1/3 * d(log(rho))/dl * E*v
    t_br = 1.d15/n_a  

    call Biseccion(log10Emax,kC_ic_st,kC_ic_IR,kC_ic_CMB,t_br,t_adv)
    Emax_e = MIN(10.d0**log10Emax,qe*B_par*R_acel)! Hillas criterion

    ! Build a logarithmic-spaced energy vector
    call vector_log(Emin_e,10.d0*Emax_e,Eint_e,E_e)

    ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculate Ka from: 
    ! int Q(E')*E'dE' = Linj. [Lrel]=erg/s, [Ka]=erg/s, [Q]=erg/erg**2/s
    
    integral=0.d0
    do i=1,mE_e
       dE_e = E_e(i)*(Eint_e-1.d0)
       integral = integral + E_e(i)**(-alpha_sw+1.d0)*exp(-E_e(i)/Emax_e)*dE_e
    enddo
    
    Ka = L_inj_e/integral
    
    ! Calculate the electron steady-state energy distribution

    open(70,file='tcoole.dat')

    !     E_e es la de los electrones; E2_e es la variable de integración.
    mE2 = mE_e                 ! Initial number of steps in E2
    do i =1,mE_e 
       dE_e = E_e(i)*(Eint_e - 1.d0)  
       S = 0.d0
       Emin2_e = E_e(i)
       Emax2_e = E_e(mE_e)
       E2_e = Emin2_e
       Eint2_e = (10.d0*Emax2_e/Emin2_e)**(1.d0/dfloat(mE2))
       do j=1,mE2
          dE2_e = E2_e * (Eint2_e - 1.d0)
          Q = E2_e**(-alpha_sw)*exp(-E2_e/Emax_e) !normalize later
          S = S + Q*dE2_e
          E2_e = E2_e*Eint2_e
       enddo
       
       mE2 = MAX(mE2-1,32) ! make it faster

       ! S is the result of the integral in [erg/erg/s]. We calculate: 
       ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE', with |dE/dt|= E/t_cool 
       ! in [erg/s], so that [N(E)]=erg/erg**2
       
       t_sy = 1.d0/1.6d-3/B_par**2.d0/E_e(i) 

       gamm = E_e(i)/mec2
       t_ic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
       t_ic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad)
       t_ic_CMB = kC_ic_CMB*gamm/F_iso(4.*gamm*T_CMB_ad)       
       t_ic = 1./(1.d0/t_ic_st + 1.d0/t_ic_IR + 1.d0/t_ic_CMB)       
       t_coole = 1./(1.d0/t_sy + 1.d0/t_ic + 1.d0/t_br + 1.d0/t_adi) 

       D_Bohm = E_e(i)*c/(3.d0*qe*B)
       t_diff = H0**2./(2.*D_Bohm) ! del Valle 2012
       
       t_esc(i) = 1./(1./t_adv + 1.d0/t_diff) ! escape
       
       t_loss = min(t_esc(i),t_coole)

       t_ac = efi_ac_sw*E_e(i)/(B_par*c*qe)

       E_punto = abs(E_e(i)/t_loss)
       N_e(i) = Ka*S/E_punto
       write(70,*) E_e(i)/eV,t_ic,t_sy,t_br,t_adv,t_ac,t_coole,t_diff,t_ic_st,t_ic_IR,t_ic_CMB,t_esc(i),t_loss

    enddo

    !===================================================
    !   Save the distribution 

    E_NT = 0.d0
    open(79,file='dist_e.dat')
    open(69,file='iny_e.dat')
    do i=1,mE_e
       dE_e = E_e(i)*(Eint_e - 1.d0)  
       E_NT = E_NT + E_e(i)*N_e(i)*dE_e  ! Control
       write(79,*) E_e(i)/eV,N_e(i)
       write(69,*) E_e(i)/eV,Ka*E_e(i)**(-alpha_sw)*exp(-E_e(i)/Emax_e)*t_esc(i)
    enddo
    write(*,80) 'E_NT=',E_NT
    close(79)
    close(69)
 
    deallocate ( N_e,E_e ) 
 
  END subroutine diste_run
  

END module diste



    !///////////////////////////////////////////////////////////////////////////
    !                                                                         //
    !                            SUBRUTINAS                                   //
    !                                                                         //
    !///////////////////////////////////////////////////////////////////////////

    Subroutine Biseccion(x3,c1,c2,c3,t_br,t_adv)
      use global
      implicit none
      integer t_err
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
      
      
      !     El error en el n-esimo será menor a |a-b|/2**n
      a_ = log10(1.d9*eV) ! log10(Emin [CGS]) ~ 1 GeV
      b_ = log10(1.d13*eV) ! log10(Emax [CGS]) ~ 10 TeV
      t_err = 4
      err = 10.d0**(-t_err) ! error menor que 10**-t
      i = 0
      do while ( abs(b_ - a_) > err )
         i = i+1
         c_ = (a_+b_)/2.
         if (resta(c_,c1,c2,c3,t_br,t_adv)==0) then
            b_ = c_
         else if (resta(a_,c1,c2,c3,t_br,t_adv)*resta(c_,c1,c2,c3,t_br,t_adv) < 0) then
            b_ = c_
         else if (resta(a_,c1,c2,c3,t_br,t_adv)*resta(c_,c1,c2,c3,t_br,t_adv) > 0) then
            a_ = c_
         endif
      enddo
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
      t_diff = H0**2./(2.*D_Bohm) ! del Valle 2012
      t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

      gamm = E_l/mec2
      t_ic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
      t_ic_IR = c2*gamm/F_iso(4.*gamm*T_IR_ad)    
!      t_ic_CMB = c3*gamm/F_iso(4.*gamm*T_CMB_ad)
      t_sy = 1.d0/1.6d-3/B_par**2.d0/E_l 

      t_coole = 1./(1.d0/t_sy + 1.d0/t_ic_st + 1.d0/t_ic_IR + 1./t_br) !+ 1.d0/t_ic_CMB 

      t_loss = min(t_coole,t_esc) ! general losses

      t_ac = efi_ac_sw*E_l/(B_par*c*qe)
      resta = t_loss - t_ac

      return

    end Function resta
    
    !/////////////////////////////////////////////////////////////////////////

