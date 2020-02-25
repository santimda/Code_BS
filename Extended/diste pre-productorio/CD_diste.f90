module CD_diste
  
  use global
  
  implicit none
  private
  public :: CD_diste_run 
  
contains
  
  !============================================================================
  ! CD_diste_run (START)
  !============================================================================
  
  subroutine CD_diste_run  
    use global
    implicit none
    !***************************************************************************
    !  This program calculates the electron energy distribution at each position
    !  of the CD using CD.dat and termo_sw.dat. The program generates diste_sw.dat, 
    !  in which Ne(l,E(i)) is saved, with m_theta*mE_e data; diste_sw_total.dat, 
    !  which stores Ne(E(i)), and ENTe.dat, with E_NT(l). Emin_e and Emax_e are 
    !  fixed in global.
    !***************************************************************************
    integer :: l_min,l2
    real(dp) :: D_diff,tdiff,integral
    real(dp) :: tcoole,tic_st,tic_IR,tic,tsy,tac,tadi,tadv,tbr,efiac
    real(dp) :: gamm, cte_IC, C_ic_st,kC_ic_st, C_ic_IR(m_theta),kC_ic_IR
    real(dp) :: log10Emaxl,R_acel,Q
    real(dp) :: E_dot,Emax_e0,Ka,E_NT,N_total_i
    real(dp) :: E_dot0(mE_e),N_e0(mE_e),n_a(m_theta),L_inj(m_theta),dR(m_theta)
    real(dp) :: L_inj_e(m_theta), sum,dummy
    
    interface
       subroutine Biseccion(Be,logEmax,c1,c2,t1,tad,xi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: Be,c1,c2,t1,tad,xi
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion
    end interface
    
    allocate ( theta(m_theta),R_theta(m_theta),H_theta(m_theta) )
    allocate ( P(m_theta),rho(m_theta),T(m_theta),B(m_theta),v_t(m_theta) )
    allocate ( N_e(mE_e),E_e(mE_e),B_par(m_theta),N_etot(m_theta,mE_e) )
    allocate ( T_IR(m_theta),T_IR_ad(m_theta) )
    allocate ( efi_ac_sw(m_theta) )
    N_etot = 0.d0

70    format(A10,X,Es10.2)

    !---------------------------------------------------------------------------------
    !  Read the position of the shock and cell size: theta, r_norm(theta), dR(theta); 
    !  the thermodynamical quantities and injected luminosity at each cell
    !---------------------------------------------------------------------------------

    open(12,file='Linj_sw.dat')
    do l=1,m_theta
       read(12,*) theta(l),R_theta(l),L_inj(l),dR(l)
    enddo
    close(12)
    R_theta = R0*R_theta     ! The file contained r_norm, not R_theta
    dR = R0*dR               ! The file contained dR/R0, not dR [cm]
    L_inj_e = iny_e*L_inj    ! The cell size is already taken into account in L_inj
    
    open(11,file='termo_sw.dat') 
    do l=1,m_theta
       read(11,*) l2,rho(l),P(l),v_t(l),T(l),B(l) ! all in CGS units
    enddo
    close (11)
    B_par = B*par        ! For synchrotron, only the parallel component of B is relevant 
    n_a = rho/mp/mu_sw   ! Ambient number density

    open(103,file='Hnorm_sw.dat') 
    do l=1,m_theta
       read(103,*) dummy,H_theta(l) ! normalized
       H_theta(l) = H_theta(l)*R_theta(l)
    enddo
    close(103)
    
    open(113,file='T_IR.dat') 
    do l=1,m_theta
!       T_IR(l) = 27.*(0.2)**(-1./6.)*(L_st/1.d38)**(1./6.)*(R_theta(l)/pc)**(-1./3.) ! Eq.10 del Valle 2012
       T_IR(l) =  98.d0             ! K2017
       T_IR_ad(l) = k*T_IR(l)/mec2  ! Adimensional dust temperature
       Kappa_IR(l) = L_IR/(4.*pi*sigma*T_IR(l)**4.*R_theta(l)**2.)
       write(113,*) theta(l),T_IR(l),T_IR_ad(l),Kappa_IR(l) ! normalized
    enddo
    close(113)

    open(114,file='efi_ac_sw.dat')
    do l=1,m_theta
       read(114,*) dummy,efi_ac_sw(l)
    enddo
    close(114)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
    !                                                                       !
    !                 ELECTRON ENERGY DISTRIBUTION                          !
    !                                                                       !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
    
    ! Build an energy vector that is common for all linear-emitters: 
    call vector_log(Emin_e,Emax_e,Eint_e,E_e)

    ! Consider different starting points for each linear emitter
    
    open(70,file='CD_tcoole.dat')
    cte_IC = pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.)
    C_ic_st = (cte_IC)/(T_st_ad**2.)
    do l=1,m_theta
       C_ic_IR(l) = (cte_IC)/(T_IR_ad(l)**2.)
    enddo
    
    Do l_min=1,m_theta
       
       ! In the starting point of the line we calculate E_max_e, the distribution 
       ! Ne_0(E_0), the cooling rate Edot0, and we sum it to Ne_tot(l,i)
       
       l = l_min
       
       Kappa_st = ( R_st/(2.*R_theta(l)) )**2.
       !Kappa_IR = L_IR/(4.*sigma*T_IR**4.*R0**2. ) ! Is it ok ????????
       !Kappa_IR = 1.d0             ! ????????
       kC_ic_st = C_ic_st/Kappa_st
       kC_ic_IR = C_ic_IR(l)/Kappa_IR(l)     

       R_acel = H_theta(l)         ! Accelerator size  ????????????
       efiac = efi_ac_sw(l)
       
       tadi=0.d0
       do l2=l,m_theta
          tadi = tadi + dR(l2)/v_t(l2)
       enddo
!       tadi = R_acel/v_t(l)        ! Adiabatic losses ?????????
       tbr = 1.d15/n_a(l)          ! It does not depend on the electron energy
       tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)
       call Biseccion(B(l),log10Emaxl,kC_ic_st,kC_ic_IR,T_IR_ad(l),tcoole,efiac)
       Emax_e0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac)) 
       if (Emax_e0 >= 0.6*Emax_e) then
          write(*,*) ''
          write(*,*) 'Emax_e in global.f90 is too small!'
          write(*,*) ''
       endif
       
       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
       ! by the condition: int Q(E')*E'dE' = Linj. [Lrel]=erg/seg, [Ka]=erg/seg, [Q]=erg/erg**2/seg
       
       integral=0.d0
       i = 1
       do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0)
          dE_e = E_e(i)*(Eint_e-1.d0)
          integral = integral + E_e(i)**(-alpha_sw+1.d0)*exp(-E_e(i)/Emax_e0)*dE_e
          i=i+1
       enddo
       Ka = L_inj_e(l)/integral ! Injected at the starting point of the line
       
       !------------------------------------------------------------------------

       N_e0 = 0.d0            ! Initialize to zero at each loop in l_min
       tadv = dR(l)/v_t(l)    ! Advection time (used only for the injection)
                
       i=1
       do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0) 
          
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
          ! in [erg/seg], dt = t_adv, so that [N(E)]=erg/erg**2
            
          Q = Ka*E_e(i)**(-alpha_sw)*exp(-E_e(i)/Emax_e0)
          N_e0(i) = Q*tadv
            
          D_diff = E_e(i)*c/(3.d0*qe*B(l)) ! Bohm
!          tdiff = R_acel**2./(2.d0*D_diff) ! Diffusion losses
          tdiff = H_theta(l)**2./(2.d0*D_diff) ! Diffusion losses

          tsy = 1.d0/1.6d-3/B_par(l)**2./E_e(i) 

          gamm = E_e(i)/mec2
          tic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
          tic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad(l))
          tic = 1./(1./tic_st + 1./tic_IR)
          tcoole = 1./(1./tsy + 1./tic + 1./tbr)! + 1./tadi)
          tac = efiac*E_e(i)/(B(l)*c*qe)
          
          E_dot0(i) = abs(E_e(i)/tcoole)
          N_etot(l_min,i) = N_etot(l_min,i) + N_e0(i)
          write(70,*) l_min,E_e(i)/eV,tic,tsy,tbr,tdiff,tac,tcoole,tadi,tadv,tic_st,tic_IR
          i = i + 1
!          if (tadv/tcoole > 0.2d0) then
!             write(*,*) l_min,i,tadv/tcoole
!          endif
       enddo         ! End loop in E_e(i)
         

!/////////////////////////////////////////////////////////////////////////


!     We now calculate the evolved version of the injected electron distribution at each
!     position "l" as: N(E,l) = N0(E0)*(Edot0/Edot)*theta(E-Emax) , for which we only need 
!     the value of Edot(l) and Emax(l)
         
       do l=l_min+1,m_theta
          
          Kappa_st = (R_st/(2.*R_theta(l)))**2.
          !Kappa_IR = 1.d0             ! ????????
          kC_ic_st = C_ic_st/Kappa_st
          kC_ic_IR = C_ic_IR(l)/Kappa_IR(l) 
          
!          R_acel = H_theta(l)
!          efiac = efi_ac_sw(l)

          ! Obtain Emax_e solving for: t_cool(E) - tac(E) = 0. 
          tadi=0.d0
          do l2=l,m_theta
             tadi = tadi + dR(l2)/v_t(l2)
          enddo
          !          tadi = R_acel/v_t(l)        ! Adiabatic losses ??????
          tbr = 1.d15/n_a(l)       ! Independent of E_e
          tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)

!          call Biseccion(B(l),log10Emaxl,kC_ic_st,kC_ic_IR,T_IR_ad(l),tcoole,efiac)
!          Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5/efiac))
          
          i = 1
          do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0) ! NOT SURE
             tsy = 1.d0/1.6d-3/B_par(l)**2.d0/E_e(i) 
             gamm = E_e(i)/mec2
             tic_st = kC_ic_st*gamm/F_iso(4.*gamm*T_st_ad)
             tic_IR = kC_ic_IR*gamm/F_iso(4.*gamm*T_IR_ad(l))
             tic = (1.d0/tic_st + 1.d0/tic_IR)**(-1.)

             D_diff = E_e(i)*c/(3.d0*qe*B(l)) ! Bohm
             !          tdiff = R_acel**2./(2.d0*D_diff) ! Diffusion losses
             tdiff = H_theta(l)**2./(2.d0*D_diff) ! Diffusion losses
             
             tcoole = (1.d0/tsy + 1.d0/tic + 1.d0/tbr + 1./tdiff)**(-1.)
!             tcoole = (1.d0/tsy + 1.d0/tic + 1.d0/tbr + 1.d0/tadi + 1./tdiff)**(-1.)
!             tac = efiac*E_e(i)/(B(l)*c*qe)
               
             E_dot = abs(E_e(i)/tcoole)
             N_e(i) = N_e0(i)*E_dot0(i)/E_dot !*exp(-E_e(i)/Emaxl) ! NOT SURE
             
             N_etot(l,i) = N_etot(l,i)+N_e(i)  ! The energies match by construction
             i = i + 1
          enddo               ! End loop in electron energies (i)

       enddo                  ! End loop in shock position (l)
  !     write(*,*) N_etot(m_theta-2,1)
    enddo                     ! End loop in starting position (l_min = # linear-emitters)
    close(70)

    ! Save the electron distribution at each position 
    open(79,file='diste_sw.dat')
    do l=1,m_theta
       do i=1,mE_e     
          write(79,*) theta(l),E_e(i)/eV,N_etot(l,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    sum = 0.d0
    open(80,file='ENTe.dat')
    do l=1,m_theta
       E_NT = 0.d0
       do i=1,mE_e
          dE_e = E_e(i)*(Eint_e - 1.d0)  
          E_NT = E_NT + E_e(i)*N_etot(l,i)*dE_e
       enddo
       E_NT = E_NT * 4.d0*mphi ! Total NT energy in each BS "ring"
       write(80,*) l, E_NT
       sum = sum + E_NT
    enddo
    close(80)
    write(*,70) 'E_tot_e=', sum ! Total NT energy in the BS
    
    open(81,file='diste_sw_total.dat')
    do i=1,mE_e
       N_total_i = 0.d0
       do l=1,m_theta
          N_total_i = N_total_i + N_etot(l,i)
       enddo
       N_total_i = N_total_i * mphi ! Total distribution at each energy for one quadrant
       write(81,*) E_e(i)/eV, N_total_i   
    enddo
    close(81)
    
    !________________________________________________
    
    deallocate ( theta,R_theta,H_theta,B,B_par,rho,P,v_t,T )
    deallocate ( N_e,E_e,N_etot ) 
    deallocate ( T_IR,T_IR_ad, efi_ac_sw )
    
  END subroutine CD_diste_run
  
  
  !============================================================================
  ! CD_diste_run (END)
  !============================================================================
    
  
END module CD_diste



!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                            SUBROUTINE                                   //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine Biseccion(Be,x3,c1,c2,T1,tad,xi)
  use global
  implicit none
  integer t_err
  real(dp) :: a_,b_,c_,err
  real(dp), intent(in) :: Be,c1,c2,T1,tad,xi
  real(dp), intent(out) :: x3
  interface
     function resta(a,B,ic1,ic2,T1,t,efi)
       use iso_fortran_env, only: dp => real64
       implicit none
       real(dp) :: resta
       real(dp), intent(in) :: a,B,ic1,ic2,T1,t,efi
     end function resta
  end interface
  
  !     El error en el n-esimo será menor a |a-b|/2**n
  a_ = log10(1.d11*eV) ! log10(Emin [CGS]) ~ 100 GeV
  b_ = log10(Emax_e)   !
  t_err = 2
  err = 10.d0**(-t_err) ! error < 10**-t
  i = 0
  do while ( abs(b_ - a_) > err )
     i = i+1
     c_ = (a_+b_)/2.
     if (resta(c_,Be,c1,c2,T1,tad,xi)==0) then
        b_ = c_
     else if (resta(a_,Be,c1,c2,T1,tad,xi)*resta(c_,Be,c1,c2,T1,tad,xi) < 0) then
        b_ = c_
     else if (resta(a_,Be,c1,c2,T1,tad,xi)*resta(c_,Be,c1,c2,T1,tad,xi) > 0) then
        a_ = c_
     endif
  enddo
  x3=b_
  return
end Subroutine Biseccion

!/////////////////////////////////////////////////////////////////////////

Function resta(log10E,B_l,c1,c2,T1,tad,xi)
  use global
  implicit none
  real(dp), intent(in) :: log10E,B_l,c1,c2,T1,tad,xi
  real(dp) :: resta
  real(dp) :: E_l,gamm,tsy,tic_st,tic_IR,tac,tcoole
  real(dp) :: D_diff,tdiff

  E_l = 10.d0**log10E
  gamm = E_l/mec2
  tsy = 1.d0/1.6d-3/(B_l*par)**2.d0/E_l   
  tic_st = c1*gamm/F_iso(4.*gamm*T_st_ad)
  tic_IR = c2*gamm/F_iso(4.*gamm*T1)  
  
  D_diff = E_l*c/(3.d0*qe*B_l) ! Bohm
  tdiff = (0.2*R0)**2./(2.d0*D_diff) ! Diffusion losses APPROX
  
  tcoole = (1.d0/tsy + 1.d0/tic_st + 1.d0/tic_IR &
       + 1.d0/tad + 1.d0/tdiff)**(-1.) ! in tad is also br

  tac = xi*E_l/(B_l*c*qe)  !      tac = (efiac*B_l*c*qe/E_l)**(-1)
  resta = tcoole - tac
!  write(*,*) tsy,tic_st,tic_IR
  return
  
end Function resta

!/////////////////////////////////////////////////////////////////////////

