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

    integer :: mE1,mE2,i1,i2
    real(dp) :: gamm, C_ic_st,kC_ic_st, C_ic_IR,kC_ic_IR, C_ic_CMB,kC_ic_CMB
    real(dp) :: tcoole,tic_st,tic_IR,tic_CMB,tic,tsy,tac,tad,tbr,tdif,t_loss
    real(dp) :: log10Emax,D_Bohm 
    real(dp) :: integral,S,Q_e,b_E,t_esc,Tau,cte_sy
    real(dp) :: E1_e,dE1_e,Emin1_e,Emax1_e,Eint1_e,sum1
    real(dp) :: E2_e,dE2_e,Emin2_e,Emax2_e,Eint2_e,sum2
    real(dp) :: E_punto,Emax_e,Ka,E_NT,L_inj_e,L_wind

    interface
       subroutine Biseccion(logEmax,c1,c2,c3,tloss)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: c1,c2,c3,tloss
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion
    end interface
    
    allocate ( N_e(mE_e),E_e(mE_e) )
80  format(A7,X,Es10.2)
810  format(A15,X,Es10.2)

    L_wind = (Mdot*v_w**2.)/2.d0
    L_inj_e = iny_e*L_wind*0.1d0 ! L_inj_wind ~ 10% L_wind

!    write(*,810) 'U_IR [eV/cm3] =',(4.*sigma/c)*T_IR**4/eV
    write(*,810) 'U_IR [eV/cm3] =',(5.7d36)/(4.*pi*H0**2.*c)/eV
    write(*,810) 'U_st [eV/cm3] =',L_st/(4.*pi*R0**2.*c)/eV
    write(*,810) 'U_IR/U_st =',( L_IR / L_st ) * ( R0 / H0 ) **2.
!    write(*,810) 'U_IR/U_st =',16.*pi*R0**2.*T_IR**4./L_st
    write(*,810) 'L_IR/L_st =',L_IR/L_st
    write(*,810) 'R_IR/R0 =',R_IR/R0
    write(*,810) 'Kappa_IR =',Kappa_IR


    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                     C
    !                              ELECTRONS                              C
    !                                                                     C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Search for maximum electron energy
    C_ic_st = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T_st_ad**2.)
    C_ic_IR = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T_IR_ad**2.)
    C_ic_CMB = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T_CMB_ad**2.)
    kC_ic_st = C_ic_st/Kappa_st
    kC_ic_IR = C_ic_IR/Kappa_IR  
    kC_ic_CMB = C_ic_CMB/Kappa_CMB
    
    tad = R0/v_w 
    tbr = 1.d15/n_a  
    t_loss = min(tbr,tad)

    call Biseccion(log10Emax,kC_ic_st,kC_ic_IR,kC_ic_CMB,t_loss)
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
    
    ! Calculate the electron steady-state energy distribution:
    ! N(E) = 1/|E^dot|_loss * integral_E^E_max dE1*Q(E1)*exp(-Tau(E,E1)/t_esc) ,
    ! with: Tau(E,E1) = integral_E^E1 dE2/|E2^dot|_loss

    open(70,file='tcoole.dat')

    S = 0.d0
    t_esc = tad                ! Ignore diffusion
    cte_sy = 1.d0/1.6d-3/B_par**2.d0
    mE1 = mE_e*16                 ! Initial number of steps in E1
    do i = 1,mE_e              ! E
       sum1 = 0.d0
       Emin1_e = E_e(i)
       Emax1_e = E_e(mE_e)
       Eint1_e = (Emax1_e/Emin1_e)**(1.d0/dfloat(mE1-1))
       E1_e = Emin1_e
       mE2 = mE_e

       do i1 = 1,mE1                   ! E1 -> N(E) = integral Q(E1)...
          sum2 = 0.d0 
          Emin2_e = E_e(i)
          Emax2_e = E1_e
          Eint2_e = (Emax2_e/Emin2_e)**(1.d0/dfloat(mE2-1))   
          E2_e = Emin2_e
    
          do i2 = 1,mE2                 ! E2 --> Tau(E,E1) = integral(E2)
             dE2_e = E2_e*(Eint2_e-1.d0)

             tsy = cte_sy/E2_e 

             gamm = E2_e/mec2
             tic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
             tic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad)
             tic_CMB = kC_ic_CMB*gamm/F_iso(4.*gamm*T_CMB_ad)       
             tic = 1./(1.d0/tic_st + 1.d0/tic_IR + 1.d0/tic_CMB)
             
!             D_Bohm = E2_e*c/(3.d0*qe*B)
!             tdif = H0**2./(2.*D_Bohm) ! del Valle 2012
             
             tcoole = 1./(1.d0/tsy + 1.d0/tic + 1/tbr) 
             b_E = E2_e/tcoole
             
             sum2 = sum2 + dE2_e/b_E
             E2_e = E2_e*Eint2_e
!             if (i > 80) then
!                write(*,*) i,dE2_e,b_E
!             endif
          enddo
          Tau = sum2
!          write(*,*) i,Tau!, (1./E_e(i)-1./E1_e)/(1.6d-3*B**2.)
!          write(*,*) 

          mE2 = MAX(mE2-1,32) !  Lower down the # of points for speed
          
          dE1_e = E1_e*(Eint1_e-1.d0)
          Q_e = E1_e**(-alpha_sw)*exp(-E1_e/Emax_e)
          sum1 = sum1 + Q_e*dE1_e*exp(-Tau/t_esc)
          E1_e = E1_e*Eint1_e
!          write(*,*) i1,Tau,Q_e,sum1,exp(-Tau/t_esc)
!          write(*,*) 
       enddo
              
       S = sum1
!       write(*,*) i,S

       mE1 = MAX(mE1-1,64) !  Lower down the # of points for speed
       
       tsy = cte_sy/E_e(i) 

       gamm = E_e(i)/mec2
       tic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
       tic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad)
       tic_CMB = kC_ic_CMB*gamm/F_iso(4.*gamm*T_CMB_ad)       
       tic = 1./(1.d0/tic_st + 1.d0/tic_IR + 1.d0/tic_CMB)
       
       D_Bohm = E_e(i)*c/(3.d0*qe*B)
       tdif = H0**2./(2.*D_Bohm) ! del Valle 2012

       tcoole = 1./(1.d0/tsy + 1.d0/tic + 1/tbr)! + 1.d0/tad + 1.d0/tdif) 
       tac = efi_ac_sw*E_e(i)/(B_par*c*qe)

       E_punto = abs(E_e(i)/tcoole)
       N_e(i) = Ka*S/E_punto
       if (tcoole > 10.*t_esc) then
          N_e(i) = Ka*E_e(i)**(-alpha_sw)*t_esc
       endif
       write(70,*) E_e(i)/eV,tic,tsy,tbr,tad,tac,tcoole,tdif,tic_st,tic_IR,tic_CMB

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
       write(69,*) E_e(i)/eV,Ka*E_e(i)**(-alpha_sw)*exp(-E_e(i)/Emax_e)*t_esc
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

    Subroutine Biseccion(x3,c1,c2,c3,tloss)
      use global
      implicit none
      integer t_err
      real(dp) :: a_,b_,c_,err
      real(dp), intent(in) :: c1,c2,c3,tloss
      real(dp), intent(out) :: x3
      interface
         function resta(a,ic1,ic2,ic3,t)
           use iso_fortran_env, only: dp => real64
           implicit none
           real(dp) :: resta
           real(dp), intent(in) :: a,ic1,ic2,ic3,t
         end function resta
      end interface
      
      
      !     El error en el n-esimo será menor a |a-b|/2**n
      a_ = log10(1.d10*eV) ! log10(Emin [CGS]) ~ 10 GeV
      b_ = log10(2.d14*eV) ! log10(Emax [CGS]) ~ 200 TeV
      t_err = 3
      err = 10.d0**(-t_err) ! error menor que 10**-t
      i = 0
      do while ( abs(b_ - a_) > err )
         i = i+1
         c_ = (a_+b_)/2.
         if (resta(c_,c1,c2,c3,tloss)==0) then
            b_ = c_
         else if (resta(a_,c1,c2,c3,tloss)*resta(c_,c1,c2,c3,tloss) < 0) then
            b_ = c_
         else if (resta(a_,c1,c2,c3,tloss)*resta(c_,c1,c2,c3,tloss) > 0) then
            a_ = c_
         endif
      enddo
      x3=b_

      return
    end Subroutine Biseccion
    
    !/////////////////////////////////////////////////////////////////////////
    
    Function resta(log10E,c1,c2,c3,tloss)
      use global
      implicit none
      real(dp), intent(in) :: log10E,c1,c2,c3,tloss
      real(dp) :: resta,E_l,gamm
      real(dp) :: tsy,tic_st,tic_IR,tic_CMB,tac,tcoole
      real(dp) :: D_Bohm,tdif

      E_l = 10.d0**log10E
      gamm = E_l/mec2

      tsy = 1.d0/1.6d-3/B_par**2.d0/E_l 
     
      tic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
      tic_IR = c2*gamm/F_iso(4.*gamm*T_IR_ad)    
      tic_CMB = c3*gamm/F_iso(4.*gamm*T_CMB_ad)

      D_Bohm = E_l*c/(3.d0*qe*B)
      tdif = H0**2./(2.*D_Bohm) ! del Valle 2012

      tcoole = 1./(1.d0/tsy + 1.d0/tic_st + 1.d0/tic_IR + &
           1.d0/tic_CMB + 1.d0/tdif + 1.d0/tloss)
      tac = efi_ac_sw*E_l/(B_par*c*qe)
      resta = tcoole - tac

      return

    end Function resta
    
    !/////////////////////////////////////////////////////////////////////////

