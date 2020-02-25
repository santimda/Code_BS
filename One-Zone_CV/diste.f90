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

    real(dp) :: gamm, C_ic_st,kC_ic_st, C_ic_CMB,kC_ic_CMB
    real(dp) :: t_ic_st,t_ic_CMB,t_ic,t_sy,t_ac,t_br,t_diff,t_adv
    real(dp) :: t_coole,t_esc(mE_e),t_loss,log_coul,gdot_coul,t_coul
    real(dp) :: log10Emax,D_Bohm,integral
    real(dp) :: Emax_e,K_e,E_NT,c_s

    interface
       subroutine Biseccion(logEmax,c1,c2,t_br,t_adi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: c1,c2,t_br,t_adi
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion
    end interface
    
    allocate( N_e(mE_e),E_e(mE_e) )
80  format(A13,X,Es10.2)

    write(*,*) 'Calculate the non-thermal electron distribution in the SW'

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                     C
    !                              ELECTRONS                              C
    !                                                                     C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! Search for maximum electron energy
    C_ic_st = ( pi*h_bar**3*c**2/(2.*re**2*mec2**3) )/(T_st_ad**2)
    C_ic_CMB = ( pi*h_bar**3*c**2/(2.*re**2*mec2**3) )/(T_CMB_ad**2)
    kC_ic_st = C_ic_st/Kappa_st
    kC_ic_CMB = C_ic_CMB/Kappa_CMB
    
    c_s = v_w/sqrt(8.d0)  ! sound speed ~ sqrt(g_ad*P_sw/rho_sw), P_sw ~ 3/4 rho_w*v_w^2
    t_adv = R0/c_s       
    t_br = 1.d15/n_a  

    call Biseccion(log10Emax,kC_ic_st,kC_ic_CMB,t_br,t_adv)
    if (10.d0**log10Emax > qe*B_par*R_acel) then
      write(*,*) 'Hillas dominates Emax!'
      Emax_e = qe*B_par*R_acel! Hillas criterion
    endif
    write(*,80) 'Emax(e)[eV] =',Emax_e/eV
    write(*,80) 'Ecool_e[eV] =',1.d0/(1.6d-3*B_par**2*t_adv*eV)
    write(*,80) 'Ebreak[eV] =',H0**2*3*qe*B/(2.d0*c*2*t_adv*eV)

    ! Build a logarithmic-spaced energy vector
    call vector_log(Emin_e,10.d0*Emax_e,Eint_e,E_e)

    ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculate Ka from: 
    ! int Q(E')*E'dE' = Linj. [Lrel]=erg/s, [Ka]=erg/s, [Q]=erg/erg**2/s
    
    integral=0.d0
    do i=1,mE_e
    	if (E_e(i) <= Ecooled_e) then
       		N_e(i) = E_e(i)**(-alpha_sw)
       	else
       		! E_corr = Ecooled_e**(alpha2-alpha)=Ecooled_e, alpha2=alpha+1
       		N_e(i) = E_e(i)**(-alpha_sw-0.99)*exp(-E_e(i)/Emax_e)*Ecooled_e**0.99
       	endif
       dE_e = E_e(i)*(Eint_e-1.d0)
       integral = integral + E_e(i)*N_e(i)*dE_e!E_e(i)**(-alpha_sw+1.d0)*exp(-E_e(i)/Emax_e)*dE_e
    enddo
    K_e = (B**2./(8.*pi)/integral)*eta_mag*eta_e
!    write(*,80) 'Ke =',K_e
    N_e = K_e*N_e*Vol
!    write(*,80) 'U_NT =',K_e*integral/eta_e
!    write(*,80) 'U_NT_e =',K_e*integral
    write(*,80) 'L_inj,NT =',Vol*K_e*integral/t_adv
    write(*,80) 'L_NT/L_w[%] =',Vol*K_e*integral/t_adv/L_wind*100.d0
    
    !=========================================================================================
    ! Calculate the electron steady-state energy distribution

    log_coul = 29.7d0 + log( T_6/sqrt(n_a) ) 
    
    open(70,file='tcoole.dat')

    do i =1,mE_e 

       t_sy = 1.d0/1.6d-3/B_par**2/E_e(i) 

       gamm = E_e(i)/mec2
       t_ic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
!       t_ic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad)
       t_ic_CMB = kC_ic_CMB*gamm/F_iso(4.*gamm*T_CMB_ad)       
       t_ic = 1./(1.d0/t_ic_st + 1.d0/t_ic_CMB)       
       gdot_coul = 3.d-14*n_a*log_coul/sqrt(gamm**2-1.d0) ! Pittard2006
       t_coul = gamm/gdot_coul

       t_coole = 1./(1.d0/t_sy + 1.d0/t_ic + 1.d0/t_br + 1.d0/t_coul) 

       D_Bohm = E_e(i)*c/(3.d0*qe*B)
       t_diff = H0**2/(2.*D_Bohm) ! del Valle 2012
       
       t_esc(i) = 1./(1./t_adv + 1.d0/t_diff) ! escape
       
       t_loss = min(t_esc(i),t_coole)

       t_ac = efi_ac_sw*E_e(i)/(B_par*c*qe)

       write(70,*) E_e(i)/eV,t_ic,t_sy,t_br,t_adv,t_ac,t_coole,t_diff,t_ic_st,t_ic_CMB,t_esc(i),t_loss,t_coul

    enddo

    !===================================================
    !   Save the distribution 

    open(79,file='dist_e.dat')
    E_NT = 0.d0
    do i=1,mE_e
       dE_e = E_e(i)*(Eint_e-1.d0)
       E_NT = E_NT + E_e(i)*N_e(i)*dE_e
       write(79,*) E_e(i)/eV,N_e(i)
    enddo
    close(79)
    write(*,80) 'E_NT(e)=', E_NT
 
    deallocate( N_e,E_e ) 
 
  END subroutine diste_run
  

END module diste



    !///////////////////////////////////////////////////////////////////////////
    !                                                                         //
    !                            SUBRUTINAS                                   //
    !                                                                         //
    !///////////////////////////////////////////////////////////////////////////

    Subroutine Biseccion(x3,c1,c2,t_br,t_adv)
      use global
      implicit none
      integer t_err
      real(dp) :: a_,b_,c_,err
      real(dp), intent(in) :: c1,c2,t_br,t_adv
      real(dp), intent(out) :: x3
      interface
         function resta(a,ic1,ic2,t_br,t_adv)
           use iso_fortran_env, only: dp => real64
           implicit none
           real(dp) :: resta
           real(dp), intent(in) :: a,ic1,ic2,t_br,t_adv
         end function resta
      end interface
            
      !     El error en el n-esimo será menor a |a-b|/2**n
      a_ = log10(1.d10*eV) ! log10(Emin [CGS]) ~ 1 GeV
      b_ = log10(1.d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
      t_err = 4
      err = 10.d0**(-t_err) ! error menor que 10**-t
      i = 0
      do while ( abs(b_ - a_) > err )
         i = i+1
         c_ = (a_+b_)/2.
         if (resta(c_,c1,c2,t_br,t_adv)==0) then
            b_ = c_
         else if (resta(a_,c1,c2,t_br,t_adv)*resta(c_,c1,c2,t_br,t_adv) < 0) then
            b_ = c_
         else if (resta(a_,c1,c2,t_br,t_adv)*resta(c_,c1,c2,t_br,t_adv) > 0) then
            a_ = c_
         endif
      enddo
      x3=b_
      return
    end Subroutine Biseccion
    
    !/////////////////////////////////////////////////////////////////////////
    
    Function resta(log10E,c1,c2,t_br,t_adv)
      use global
      implicit none
      real(dp), intent(in) :: log10E,c1,c2,t_br,t_adv
      real(dp) :: resta,E_l,gamm
      real(dp) :: t_sy,t_ic_st,t_ic_CMB,t_ac,t_coole
      real(dp) :: D_Bohm,t_diff,t_loss,t_esc

      E_l = 10.d0**log10E

      D_Bohm = E_l*c/(3.d0*qe*B)
      t_diff = H0**2./(2.*D_Bohm) ! del Valle 2012
      t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape

      gamm = E_l/mec2
      t_ic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
      t_ic_CMB = c2*gamm/F_iso(4.*gamm*T_CMB_ad)
      t_sy = 1.d0/1.6d-3/B_par**2.d0/E_l 

      t_coole = 1./(1.d0/t_sy + 1.d0/t_ic_st + 1.d0/t_ic_CMB + 1./t_br)

      t_loss = min(t_coole,t_esc) ! general losses

      t_ac = efi_ac_sw*E_l/(B_par*c*qe)
      resta = t_loss - t_ac

      return

    end Function resta
    
    !/////////////////////////////////////////////////////////////////////////

