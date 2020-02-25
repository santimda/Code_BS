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
    real(dp) :: t_cell(m_theta),t_adv(m_theta),t_coolp(m_theta,mE_p)
    real(dp) :: t_ac,t_pp,t_diff,D_diff,t_esc,t_loss(m_theta,mE_p)
    real(dp) :: log10Emaxl,R_acel,integral,Q
    real(dp) :: E_dot(m_theta,mE_p),Emax_p0(m_theta),Ka,E_NT,N_totalp_i
    real(dp) :: N_p0(mE_p),n_a(m_theta),L_inj(m_theta),dR(m_theta)
    real(dp) :: L_inj_p(m_theta), sum,dummy
    real(dp) :: dE_aprox ,cool,t_adi(m_theta)
    integer :: index1,index2,count,steps 
    logical :: new_count
    integer,parameter :: N_bins = 5 ! For the plot diste_sw_bins
    real(dp) :: N_bin_i(N_bins,mE_p)
    
    interface
       subroutine Biseccion_p(Be,logEmax,tpp,tadv,xi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: Be,tpp,tadv,xi
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion_p
    end interface
    
    allocate ( theta(m_theta),R_theta(m_theta),H_theta(m_theta) )
    allocate ( P(m_theta),rho(m_theta),T(m_theta),B(m_theta),v_t(m_theta) )
    allocate ( N_p(mE_p),E_p(mE_p),B_par(m_theta),N_ptot(m_theta,mE_p) )
    allocate ( efi_ac_sw(m_theta) )
    N_ptot = 0.d0
    
70  format(A10,X,Es10.2)
    
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
        
    !========================================================================
    ! First calculate the cooling time at all energies and at each position 
    ! (and also de maximum energy)
    
    do l = 1, m_theta
       t_cell(l) = dR(l)/v_t(l) ! Cell advection time
       t_adv(l) = R_theta(l)/v_t(l) ! Characteristic advection time
       if (l < m_theta) then
          t_adi(l) = 3.*dR(l)/(v_t(l)*log(rho(l)/rho(l+1))) ! Characteristic adiabatic losses time
       else
          t_adi(l) = t_adi(l-1)
       endif
    enddo
    
!    t_adv = 0.d0    ! advection time
    do l=1,m_theta
!       if (l == 1) then
!          do l2=l,m_theta
!             t_adv(l) = t_adv(l) + t_cell(l2)
!          enddo
!       else
!          t_adv(l) = t_adv(l-1) - t_cell(l-1) ! PREGUNTAR V!!!
!       endif
       
       R_acel = H_theta(l)         ! Accelerator size  ????????????
       
       t_pp = 1.d15/n_a(l)       ! p-p losses (approx)
       
       call Biseccion_p(B(l),log10Emaxl,t_pp,t_adv(l),efi_ac_sw(l))
       Emax_p0(l) = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efi_ac_sw(l)))! Pérdidas difusivas
       if (Emax_p0(l) >= 0.5*Emax_p) then
          write(*,*) 'Emax_p in global.f90 is too small!'
       endif
       
       do i=1,mE_p
          
          D_diff = E_p(i)*c/(3.d0*qe*B(l))
          t_diff = H_theta(l)**2./(2.d0*D_diff) ! Diffusion losses
          t_esc = 1./( 1./t_adv(l) + 1.d0/t_diff ) ! escape

          t_coolp(l,i) = 1./(1.d0/t_pp + 1./t_adi(l))

          t_loss(l,i) = min(t_coolp(l,i),t_esc) ! general losses

          t_ac = efi_ac_sw(l)*E_p(i)/(B(l)*c*qe)
          
          write(70,*) l,E_p(i)/eV,t_ac,t_pp,t_adv(l),t_diff,t_coolp(l,i),t_cell(l),t_esc,t_loss(l,i),t_adi(l)
         
          E_dot(l,i) = E_p(i)/t_coolp(l,i)
          
       enddo
       
    enddo
    close(70)
    
    !========================================================================
    
    N_ptot = 0.d0
    
    ! At the starting point of the line we calculate the distribution 
    ! Np_0(E_0), the cooling rate Edot0, and we sum it to Np_tot(l,i)
       
    Do l_min=1,m_theta

       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. The normalization constant Ka is obtained 
       ! by the condition: int Q(E')*E'dE' = Linj. [Lrel]=erg/seg, [Ka]=erg/seg, [Q]=erg/erg**2/seg
       
       integral=0.d0
       i = 1
       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0(l_min))
          dE_p = E_p(i)*(Eint_p - 1.d0)
          integral = integral + E_p(i)**(-alpha_sw+1.d0)*exp(-E_p(i)/Emax_p0(l_min))*dE_p
          i=i+1
       enddo
       
       Ka = L_inj_p(l_min)/integral ! Injected at the starting point of the line

       !------------------------------------------------------------------------

       N_p0 = 0.d0            ! Initialize to zero at each loop in l_min
       
       i=1

       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0(l_min)) 
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, with |dE/dt|= E/t_cool 
          ! in [erg/seg], dt = min(t_cell,t_cool) so that [N(E)]=erg/erg**2
                    
          Q = Ka*E_p(i)**(-alpha_sw)*exp(-E_p(i)/Emax_p0(l_min)) 
          N_p0(i) = Q*min(t_cell(l_min),t_coolp(l_min,i))
          N_ptot(l_min,i) = N_ptot(l_min,i) + N_p0(i)
            
          ! We calculate the evolved version of the injected electron distribution at each
          ! position "l" as: N_l(l) = N_l-1* (Edot_l-1/Edot_l), where the energy also evolves
          ! as E0(l) given by: t_cell(l) = int_E0^El Edot * dE
          
          index1 = i            ! For the cell "l"
          index2 = i            ! For the previous cell, "l-1"
          N_p = 0.d0            ! Initialize to zero at each loop in l_min
          N_p(index1) = N_p0(i) ! This will be evolved
          l = l_min + 1  
          count = 0
          cool = 0.d0
          new_count = .true.     
     
          do while( l <= m_theta .and. index1 > 1) ! if index1 = 0, the e- are completely cooled
             ! At the position l, the particle has E_e(index1) instead of E_e(i)
             ! We find index1 <= i according to the cooling time:

             if (l == l_min+1) then
                do l2=l,m_theta
                   cool = cool + t_cell(l2)/t_coolp(l2,index1)
                enddo
             else
                cool = cool - t_cell(l-1)/t_coolp(l-1,index1)
             endif
             dE_aprox = (1.-t_cell(l)/t_coolp(l,index1))

             if ( 1.-cool > 0.99/Eint_p ) then
                continue
!                write(*,*) 'hola1',l_min,i,l,index1

             elseif ( 1.-cool <= 0.99/Eint_p .and. dE_aprox > 1./Eint_p ) then
                if (new_count .eqv. .true.) then
                   sum = 1.d0
                   steps = 0
                   l2 = l
                   do while (sum > 1./Eint_p .and. l2 <= m_theta)
                      steps = steps + 1
                      sum = sum*(1.-t_cell(l2)/t_coolp(l2,index1))
                      l2 = l2 + 1
                      if (sum > 1./Eint_p .and. l2 == m_theta) then
                         steps = m_theta ! NOT COOLING
                      endif
                   enddo
                   new_count = .false.
                endif
                count = count + 1
!                write(*,*) 'hola2',i,l,int(1./(Eint_p*dE_aprox)),steps
                if (count >= steps) then
                   count = 0
                   index1 = index1 - 1
                   new_count = .true. 
                endif

             elseif ( dE_aprox <= 1./Eint_p ) then
                sum = 0.d0
                do while (sum < t_cell(l) .and. index1 > 1)
                   index1 = index1 - 1    
                   dE_p = E_p(index1)*(Eint_p-1.d0)
                   sum = sum + dE_p/E_dot(l,index1)        
                enddo
 !               write(*,*) 'hola3',i,l,index2,index1
                
             endif

             N_p(index1) = N_p(index2) * ( E_dot(l,index2)/E_dot(l,index1) )

             N_ptot(l,index1) = N_ptot(l,index1)+N_p(index1)  ! The energies match by construction
             
             index2 = index1

             l = l + 1 
          enddo               ! End loop in shock position (l) 
          i = i + 1
       enddo                  ! End loop in electron starting energies (i)

    enddo                     ! End loop in starting position (l_min = # linear-emitters)
    
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
    N_bin_i = 0.d0
    do i=1,mE_p
       N_totalp_i = 0.d0
       do l=1,m_theta
          ! Make some distributions binned in intervals of (0.75/N_bins)*pi
          if (theta(l) < 0.15*pi) then
             N_bin_i(1,i) = N_bin_i(1,i) + N_ptot(l,i)
          elseif (theta(l) >= 0.15*pi .and. theta(l) < 0.30*pi) then
             N_bin_i(2,i) = N_bin_i(2,i) + N_ptot(l,i)
          elseif (theta(l) >= 0.30*pi .and. theta(l) < 0.45*pi) then
             N_bin_i(3,i) = N_bin_i(3,i) + N_ptot(l,i)
          elseif (theta(l) >= 0.45*pi .and. theta(l) < 0.6*pi) then
             N_bin_i(4,i) = N_bin_i(4,i) + N_ptot(l,i)
          elseif (theta(l) >= 0.60*pi .and. theta(l) < 0.75*pi) then
             N_bin_i(5,i) = N_bin_i(5,i) + N_ptot(l,i)
          endif             
          N_totalp_i = N_totalp_i + N_ptot(l,i)
       enddo
       N_totalp_i = N_totalp_i * mphi*1.d0 ! Total distribution at each energy for 1 of the 4 quadrants
       write(81,*) E_p(i)/eV, N_totalp_i
    enddo
    close(81)
    
    open(82,file='distp_sw_bin.dat')
    do j = 1,N_bins
       do i = 1, mE_p
          write(82,*) j,E_p(i)/eV, N_bin_i(j,i)* mphi*1.d0
       enddo
       write(82,*) 
    enddo
    close(82)
    
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
  
  Subroutine Biseccion_p(Be,x3,t_pp,t_adv,xi)
    use global
    implicit none
    integer :: t_err
    real(dp) :: a_,b_,c_,err
    real(dp), intent(in) :: Be,t_pp,t_adv,xi
    real(dp), intent(out) :: x3
    interface
       function resta_p(a,B,t_pp,t_adv,efi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: resta_p
         real(dp), intent(in) :: a,B,t_pp,t_adv,efi
       end function resta_p
    end interface
    
    !     El error en el n-esimo será menor a |a-b|/2**n
    a_ = log10(1.d12*eV) ! log10(Emin [CGS]) ~ 100 GeV
    b_ = log10(1.d15*eV) ! log10(Emax [CGS]) ~ 1 PeV
    t_err = 2
    err = 10.d0**(-t_err) ! error < 10**-t
    i = 0
    do while ( abs(b_ - a_) > err )
       i = i+1
       c_ = (a_+b_)/2.
       if (resta_p(c_,Be,t_pp,t_adv,xi)==0) then
          b_ = c_
       else if (resta_p(a_,Be,t_pp,t_adv,xi)*resta_p(c_,Be,t_pp,t_adv,xi) < 0) then
          b_ = c_
       else if (resta_p(a_,Be,t_pp,t_adv,xi)*resta_p(c_,Be,t_pp,t_adv,xi) > 0) then
          a_ = c_
       endif
    enddo
    x3=b_
    return
  end Subroutine Biseccion_p
  
  !/////////////////////////////////////////////////////////////////////////
  
  Function resta_p(log10E,B_l,t_pp,t_adv,xi)
    use global
    implicit none
    real(dp), intent(in) :: log10E,B_l,t_pp,t_adv,xi
    real(dp) :: resta_p,E_l,t_diff,t_ac,t_coolp
    real(dp) :: D_diff,t_loss,t_esc
    
    E_l = 10.d0**log10E

    D_diff = E_l*c/(3.d0*qe*B_l) ! Bohm
    t_diff = (0.2*R0)**2./(2.d0*D_diff) ! Diffusion losses APPROX
    t_esc = 1./(1./t_adv + 1.d0/t_diff) ! escape
    
    t_coolp = t_pp

    t_loss = min(t_coolp,t_esc) ! general losses

    t_ac = xi*E_l/(B_l*c*qe)    
    resta_p = t_loss - t_ac
    
    return
    
  end Function resta_p
  
  !/////////////////////////////////////////////////////////////////////////

