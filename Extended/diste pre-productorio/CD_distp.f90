module CD_distp
  
  use global
  
  implicit none
  private
  public :: CD_distp_run

contains
  
  !============================================================================
  ! CD_distp_run
  !============================================================================
  
  subroutine CD_distp_run  
    use global
    implicit none
    !***************************************************************************
    !  This program calculates the proton energy distribution at each position
    !  of the CD using CD.dat and termo_sw.dat. The program generates distp_sw.dat, 
    !  in which Np(l,E(i)) is saved, with m_theta*mE_e data; distp_sw_total.dat, 
    !  which stores Np(E(i)), and ENTp.dat, with E_NT(l). Emin_p and Emax_p are 
    !  fixed in global.
    !***************************************************************************
    integer :: l_min,l2
    real(dp) :: t_coolp,t_ac,t_adi,t_adv,t_pp,efiac,t_diff,D_diff,t_adpp
    real(dp) :: log10Emaxl,R_acel,integral,Q
    real(dp) :: E_dot,Emax_p0,Ka,E_NT,N_totalp_i
    real(dp) :: E_dot0(mE_p),N_p0(mE_p),n_a(m_theta),L_inj(m_theta),dR(m_theta)
    real(dp) :: L_inj_p(m_theta), sum,dummy
    
    interface
       subroutine Biseccion_p(Be,logEmax,tloss,xi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: Be,tloss,xi
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion_p
    end interface

    allocate ( theta(m_theta),R_theta(m_theta),H_theta(m_theta) )
    allocate ( P(m_theta),rho(m_theta),T(m_theta),B(m_theta),v_t(m_theta) )
    allocate ( N_p(mE_p),E_p(mE_p),B_par(m_theta),N_ptot(m_theta,mE_p) )
    allocate ( efi_ac_sw(m_theta) )
    N_ptot = 0.d0

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
    L_inj_p = iny_p*L_inj    ! The cell size is already taken into account in L_inj
    
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
    
    open(114,file='efi_ac_sw.dat')
    do l=1,m_theta
       read(114,*) dummy,efi_ac_sw(l)
    enddo
    close(114)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                     C
    !                              PROTONS                                C
    !                                                                     C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Build an energy vector that is common for all linear-emitters: 
    call vector_log(Emin_p,Emax_p,Eint_p,E_p)

    ! Consider different starting points for each linear emitter
    
    open(70,file='CD_tcoolp.dat')
        
    Do l_min=1,m_theta

       ! In the starting point of the line we calculate E_max_p, the distribution 
       ! Np_0(E_0), the cooling rate Edot0, and we sum it to Np_tot(l,i)
       
       l = l_min
       
       R_acel = H_theta(l)         ! Accelerator size  ????????????
       efiac = efi_ac_sw(l)

       t_adi=0.d0
       do l2=l,m_theta
          t_adi = t_adi + dR(l2)/v_t(l2)
       enddo
!       t_adi = R_acel/v(l)       ! Adiabatic losses ?????????
       t_pp = 1.d15/n_a(l)       ! p-p losses (approx)
       t_adpp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      
       
       call Biseccion_p(B(l),log10Emaxl,t_adpp,efiac)
       Emax_p0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))! Pérdidas difusivas
       if (Emax_p0 >= 0.5*Emax_p) then
          write(*,*) ''
          write(*,*) 'Emax_p in global.f90 is too small!'
          write(*,*) ''
       endif

       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
       ! by the condition: int Q(E')*E'dE' = Linj. [Lrel]=erg/seg, [Ka]=erg/seg, [Q]=erg/erg**2/seg
              
       integral=0.d0
       i = 1
       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0)
          dE_p = E_p(i)*(Eint_p - 1.d0)
          integral = integral + E_p(i)**(-alpha_sw+1.d0)*exp(-E_p(i)/Emax_p0)*dE_p
          i=i+1
       enddo
       
       Ka = L_inj_p(l)/integral ! En cada celda inyecto esa fracción de L_w

       !------------------------------------------------------------------------

       N_p0 = 0.d0           ! Initialize to zero at each loop in l_min 
       t_adv = dR(l)/v_t(l)  ! Advection time (used only for the injection)
       
       i=1
       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0) 
          
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
          ! in [erg/seg], dt = t_adv, so that [N(E)]=erg/erg**2
          
          Q = Ka*E_p(i)**(-alpha_sw)*exp(-E_p(i)/Emax_p0) 
          N_p0(i) = Q*t_adv
          
          D_diff = E_p(i)*c/(3.d0*qe*B(l))
          t_diff = H_theta(l)**2./(2.d0*D_diff) ! Diffusion losses
          t_coolp = 1.d0/(1.d0/t_diff + 1.d0/t_pp) !t_adpp

          t_ac = efiac*E_p(i)/(B(l)*c*qe)
            
          E_dot0(i) = abs(E_p(i)/t_coolp)
          N_ptot(l_min,i) = N_ptot(l_min,i) + N_p0(i)
          write(70,*) l_min,E_p(i)/eV,t_ac,t_pp,t_adv,t_adi,t_diff,t_coolp
          i = i + 1
          
       enddo
       

!/////////////////////////////////////////////////////////////////////////


!     We now calculate the evolved version of the injected electron distribution at each
!     position "l" as: N(E,l) = N0(E0)*(Edot0/Edot)*theta(E-Emax) , for which we only need 
!     the value of Edot(l) and Emax(l)
       
       do l=l_min+1,m_theta
          
          ! Obtain Emax_e solving for: t_cool(E) - tac(E) = 0. 
          
!          R_acel = H_theta(l)
!          efiac = efi_ac_sw(l)

          t_adi=0.d0
          do l2=l,m_theta
             t_adi = t_adi + dR(l2)/v_t(l2)
          enddo
          !          t_adi = R_acel/v_t(l)        ! Adiabatic losses ??????
          t_pp = 1.d15/n_a(l)       ! pérdidas protón-protón
          t_adpp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      
!          call Biseccion_p(B(l),log10Emaxl,t_coolp,efiac)       
!          Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))
          
          i = 1
          do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0) ! NOT SURE

             !               t_ac = efiac*E_p(i)/(B(l)*c*qe)  
             D_diff = E_p(i)*c/(3.d0*qe*B(l)) ! Bohm
             !          tdiff = R_acel**2./(2.d0*D_diff) ! Diffusion losses
             t_diff = H_theta(l)**2./(2.d0*D_diff) ! Diffusion losses
!             t_coolp = 1.d0/(1.d0/t_adpp + 1.d0/t_diff)   
             t_coolp = 1.d0/(1.d0/t_diff + 1.d0/t_pp) !t_adpp   
             
             E_dot = abs(E_p(i)/t_coolp)
             N_p(i) = N_p0(i)*E_dot0(i)/E_dot !*exp(-E_p(i)/Emaxl)! NOT SURE
             
             N_ptot(l,i) = N_ptot(l,i)+N_p(i)  ! The energies match by construction
             i = i + 1

          enddo                ! End loop in proton energies (i)
          
       enddo                   ! End loop in shock position (l)
       
    Enddo                      ! End loop in starting position (l_min = # linear-emitters)
    
    ! Save the proton distribution at each position 
    open(79,file='distp_sw.dat')
    do l=1,m_theta
       do i=1,mE_p
          write(79,*) theta(l),E_p(i)/eV,N_ptot(l,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    sum = 0.d0
    open(80,file='ENTp.dat')
    do l=1,m_theta
       E_NT = 0.d0
       do i=1,mE_p
          dE_p = E_p(i)*(Eint_p - 1.d0)  
          E_NT = E_NT + E_p(i)*N_ptot(l,i)*dE_p
       enddo
       E_NT = E_NT * 4.d0*mphi ! Total NT energy in each BS "ring"
       write(80,*) l, E_NT
       sum = sum + E_NT
    enddo
    close(80)
    write(*,70) 'E_tot_p=', sum ! Total NT energy in the BS
    
    open(81,file='distp_sw_total.dat')
    do i=1,mE_p
       N_totalp_i = 0.d0
       do l=1,m_theta
          N_totalp_i = N_totalp_i + N_ptot(l,i)
       enddo
       N_totalp_i = N_totalp_i * mphi ! Total distribution at each energy for one quadrant
       write(81,*) E_p(i)/eV, N_totalp_i
    enddo
    close(81)
    
    !________________________________________________
    
    deallocate ( theta,R_theta,H_theta,B,B_par,rho,P,v_t,T )
    deallocate ( N_p,E_p,N_ptot, efi_ac_sw )
    
  END subroutine CD_distp_run
  
  
  !============================================================================
  ! CD_distp_run (END)
  !============================================================================
       
    
END module CD_distp

  
  
  !///////////////////////////////////////////////////////////////////////////
  !                                                                         //
  !                            SUBROUTINE                                   //
  !                                                                         //
  !///////////////////////////////////////////////////////////////////////////
  
  Subroutine Biseccion_p(Be,x3,tloss,xi)
    use global
    implicit none
    integer :: t_err
    real(dp) :: a_,b_,c_,err
    real(dp), intent(in) :: Be,tloss,xi
    real(dp), intent(out) :: x3
    interface
       function resta_p(a,B,tloss,efi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: resta_p
         real(dp), intent(in) :: a,B,tloss,efi
       end function resta_p
    end interface
    
    !     El error en el n-esimo será menor a |a-b|/2**n
    a_ = log10(1.d11*eV) ! log10(Emin [CGS]) ~ 100 GeV
    b_ = log10(1.d15*eV) ! log10(Emax [CGS]) ~ 1 PeV
    t_err = 3
    err = 10.d0**(-t_err) ! error < 10**-t
    i = 0
    do while ( abs(b_ - a_) > err )
       i = i+1
       c_ = (a_+b_)/2.
       if (resta_p(c_,Be,tloss,xi)==0) then
          b_ = c_
       else if (resta_p(a_,Be,tloss,xi)*resta_p(c_,Be,tloss,xi) < 0) then
          b_ = c_
       else if (resta_p(a_,Be,tloss,xi)*resta_p(c_,Be,tloss,xi) > 0) then
          a_ = c_
       endif
    enddo
    x3=b_
    return
  end Subroutine Biseccion_p
  
  !/////////////////////////////////////////////////////////////////////////
  
  Function resta_p(log10E,B_l,t_loss,xi)
    use global
    implicit none
    real(dp), intent(in) :: log10E,B_l,t_loss,xi
    real(dp) :: resta_p,E_l,t_diff,t_ac,t_coolp
    real(dp) :: D_diff 
    
    E_l = 10.d0**log10E

    D_diff = E_l*c/(3.d0*qe*B_l) ! Bohm
    t_diff = (0.2*R0)**2./(2.d0*D_diff) ! Diffusion losses APPROX
    
    t_coolp = 1./(1.d0/t_loss + 1.d0/t_diff) ! tloss --> adv,p-p

    t_ac = xi*E_l/(B_l*c*qe)    
    
    resta_p = t_coolp - t_ac
    
    return
    
  end Function resta_p
  
  !/////////////////////////////////////////////////////////////////////////

