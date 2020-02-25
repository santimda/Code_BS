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
    real(dp) :: t_cell(m_theta),t_adv(m_theta),t_coole(m_theta,mE_e)
    real(dp) :: tsy,tac
    real(dp) :: log10Emaxl,R_acel,Q,integral
    real(dp) :: E_dot(m_theta,mE_e),Emax_e0(m_theta),Ka,E_NT,N_total_i
    real(dp) :: N_e0(mE_e),L_inj(m_theta)
    real(dp) :: L_inj_e(m_theta), sum !, E_aux(m_theta)!E_dot0(mE_e),
    real(dp) :: dE_aprox ,cool,N_NT !,E_aux2
    integer :: index1,index2,count,steps !,index1_aux
    logical :: new_count
    
    interface
       subroutine Biseccion(Be,logEmax,tad,xi)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: Be,tad,xi
         real(dp), intent(out) :: logEmax
       end subroutine Biseccion
    end interface
    
    allocate ( N_e(mE_e),E_e(mE_e),N_etot(m_theta,mE_e) )

70    format(A10,X,Es10.2)

    !---------------------------------------------------------------------------------
    !  Read the position of the shock and cell size: theta, r_norm(theta), dR(theta); 
    !  the thermodynamical quantities and injected luminosity at each cell
    !---------------------------------------------------------------------------------

    open(12,file='Linj_sw.dat')
    do l=1,m_theta
       read(12,*) i,z(l),L_inj(l),dR(l)
    enddo
    close(12)
    dR = R0*dR               ! The file contained dR/R0, not dR [cm]
    L_inj_e = iny_e*L_inj    ! The cell size is already taken into account in L_inj
    
    open(11,file='termo_sw.dat') 
    do l=1,m_theta
       read(11,*) l2,B(l) ! all in CGS units
    enddo
    close (11)
    B_par = B*par        ! For synchrotron, only the parallel component of B is relevant 

    !!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!!
    !!                                                                       !!
    !!                 ELECTRON ENERGY DISTRIBUTION                          !!
    !!                                                                       !!
    !!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!!
    
    ! Build an energy vector that is common for all linear-emitters: 
    call vector_log(Emin_e,Emax_e,Eint_e,E_e)
    ! Consider different starting points for each linear emitter
    
    open(70,file='CD_tcoole.dat')

    !========================================================================
    ! First calculate the cooling time at all energies and at each position 
    ! (and also de maximum energy)

    do l = 1, m_theta
       t_cell(l) = dR(l)/v_w ! Cell advection time
    enddo

    t_adv = 0.d0    ! advection time
    do l=1,m_theta
     !  if (l == 1) then
     !     do l2=l,m_theta
     !        t_adv(l) = t_adv(l) + t_cell(l2)
     !     enddo
     !  else
     !     t_adv(l) = t_adv(l-1) !- t_cell(l-1)
     !  endif
       t_adv(l) = z(l)/v_w
       
       R_acel = R0*(z(l)/z(1))         ! Accelerator size  ????????????

       call Biseccion(B(l),log10Emaxl,t_adv(l),efi_ac)
       Emax_e0(l) = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efi_ac)) 
       if (Emax_e0(l) >= 0.5*Emax_e) then
          write(*,*) 'Emax_e in global.f90 is too small!'
       endif
!       write(*,70) 'Emax=',Emax_e0(l)/eV
       do i=1,mE_e
            
          tsy = 1.d0/1.6d-3/B_par(l)**2./E_e(i) 

          t_coole(l,i) = 1./(1./tsy)
          tac = efi_ac*E_e(i)/(B(l)*c*qe)
          
          write(70,*) l,E_e(i)/eV,tsy,tac,t_coole(l,i),t_adv(l),t_cell(l)

          E_dot(l,i) = E_e(i)/t_coole(l,i)
          
       enddo
       
    enddo
    close(70)
    
    !========================================================================
    
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
       enddo
       Ka = L_inj_e(l_min)/integral ! Injected at the starting point of the line
       
       !------------------------------------------------------------------------
       
       N_e0 = 0.d0            ! Initialize to zero at each loop in l_min
       i=1

       do while (i < mE_e+1 .and. E_e(i) < 10.d0*Emax_e0(l_min)) 
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' ~ Q(E)*dt, with |dE/dt|= E/t_cool 
          ! in [erg/seg], so that [N(E)]=erg/erg**2; dt ~ min(t_cell,t_cool)
            
          Q = Ka*E_e(i)**(-alpha_sw)*exp(-E_e(i)/Emax_e0(l_min))
          N_e0(i) = Q*min(t_cell(l_min),t_coole(l_min,i))
!         write(*,*) l_min,i,Q,N_e0(i),N_etot(l_min,i)
          N_etot(l_min,i) = N_etot(l_min,i)+N_e0(i)  

          ! We calculate the evolved version of the injected electron distribution at each
          ! position "l" as: N_l(l) = N_l-1* (Edot_l-1/Edot_l), where the energy also evolves
          ! as E0(l) given by: t_cell(l) = int_E0^El Edot * dE
          
!          E_aux(l_min) = E_e(i)
          index1 = i            ! For the cell "l"
          index2 = i            ! For the previous cell, "l-1"
          N_e = 0.d0            ! Initialize to zero at each loop in "i"
          N_e(index1) = N_e0(i) ! This will be evolved
          l = l_min + 1       
          count = 0
          cool = 0.d0
          new_count = .true.

          do while( l <= m_theta .and. index1 > 0) ! if index1 = 0, the e- are completely cooled
             ! At the position l, the particle has E_e(index1) instead of E_e(i)
             ! We find index1 <= i according to the cooling time:
             
             if (l == l_min+1) then
                do l2=l,m_theta
                   cool = cool + t_cell(l2)/t_coole(l2,index1)
                enddo
             else
                cool = cool - t_cell(l-1)/t_coole(l-1,index1)
             endif
             dE_aprox = (1.-t_cell(l)/t_coole(l,index1))
             
             if ( 1.-cool > 0.99/Eint_e ) then
                continue
!                write(*,*) 'hola1',l_min,i,l,index1

             elseif ( 1.-cool <= 0.99/Eint_e .and. dE_aprox > 1./Eint_e ) then
                if (new_count .eqv. .true.) then
!                   steps = int( -log(Eint_e)/log(dE_aprox))
!                   write(*,*) 'hola2',i,l,steps
                   sum = 1.d0
                   steps = 0
                   l2 = l
                   do while (sum > 1./Eint_e .and. l2 <= m_theta)
                      steps = steps + 1
                      sum = sum*(1.-t_cell(l2)/t_coole(l2,index1))
                      l2 = l2 + 1
                      if (sum > 1./Eint_e .and. l2 == m_theta) then
                         steps = m_theta ! NOT COOLING
                      endif
                   enddo
!                   steps = steps - 1
!                   write(*,*) 'hola2',i,l,steps
                   new_count = .false.
                endif
                count = count + 1
!                steps = int(1./(Eint_e*dE_aprox))
!                write(*,*) 'hola2',i,l,int(1./(Eint_e*dE_aprox)),steps
                if (count >= steps) then
                   count = 0
                   index1 = index1 - 1
                   new_count = .true. 
                endif

             elseif ( dE_aprox <= 1./Eint_e ) then
!                E_aux2 = E_e(index1) / ( 1. + t_cell(l)/t_coole(l,index1) )
!                index1_aux = max( int( log(E_aux2/E_e(1))/log(Eint_e) ) + 1 , 1)
!                write(*,*) 'hola3',i,l,index2,index1_aux
                sum = 0.d0
                do while (sum < t_cell(l) .and. index1 > 1)
                   index1 = index1 - 1    
                   dE_e = E_e(index1)*(Eint_e-1.d0)
                   sum = sum + dE_e/E_dot(l,index1)
!                   E_aux(l) = E_e(index1)
!                   write(*,*) l_min,i,l,index2,index1          
                enddo
 !               write(*,*) 'hola3',i,l,index2,index1
                
             endif
!             write(*,*) l_min,i,l,index2,index1     
!             write(*,*) i,l,index2,index1,E_dot(l,index1)/E_dot(l,index2) 
!             write(*,*) l_min,i,l,E_e(index1)/E_e(index2),N_e(index1),N_e(index1)*( E_dot(l,index2)/E_dot(l,index1) )

             N_e(index1) = N_e(index2) * ( E_dot(l,index2)/E_dot(l,index1) )

             N_etot(l,index1) = N_etot(l,index1)+N_e(index1)  ! The energies match by construction
             
             index2 = index1

             l = l + 1 
          enddo               ! End loop in shock position (l) 
          i = i + 1
       enddo                  ! End loop in electron starting energies (i)

    enddo                     ! End loop in starting position (l_min = # linear-emitters)

    ! Save the electron distribution at each position 
    open(79,file='diste_sw.dat')
    do l=1,m_theta
       do i=1,mE_e     
          write(79,*) l,E_e(i)/eV,N_etot(l,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    sum = 0.d0
    open(80,file='ENTe.dat')
    do l=1,m_theta
       E_NT = 0.d0
       N_NT = 0.d0
       do i=1,mE_e
          dE_e = E_e(i)*(Eint_e - 1.d0)  
          E_NT = E_NT + E_e(i)*N_etot(l,i)*dE_e
          N_NT = N_NT + N_etot(l,i)*dE_e
       enddo
!       write(*,*) l, N_NT
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
       write(81,*) E_e(i)/eV, N_total_i   
    enddo
    close(81)
    
    !________________________________________________
    
    deallocate ( N_e,E_e,N_etot ) 
    
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

Subroutine Biseccion(Be,x3,tad,xi)
  use global
  implicit none
  integer t_err
  real(dp) :: a_,b_,c_,err
  real(dp), intent(in) :: Be,tad,xi
  real(dp), intent(out) :: x3
  interface
     function resta(a,B,t,efi)
       use iso_fortran_env, only: dp => real64
       implicit none
       real(dp) :: resta
       real(dp), intent(in) :: a,B,t,efi
     end function resta
  end interface
  
  !     El error en el n-esimo será menor a |a-b|/2**n
  a_ = log10(5.d10*eV) ! log10(Emin [CGS]) ~ 50 GeV
  b_ = log10(Emax_e)   !
  t_err = 2
  err = 10.d0**(-t_err) ! error < 10**-t
  i = 0
  do while ( abs(b_ - a_) > err )
     i = i+1
     c_ = (a_+b_)/2.
     if (resta(c_,Be,tad,xi)==0) then
        b_ = c_
     else if (resta(a_,Be,tad,xi)*resta(c_,Be,tad,xi) < 0) then
        b_ = c_
     else if (resta(a_,Be,tad,xi)*resta(c_,Be,tad,xi) > 0) then
        a_ = c_
     endif
  enddo
  x3=b_
  return
end Subroutine Biseccion

!/////////////////////////////////////////////////////////////////////////

Function resta(log10E,B_l,tad,xi)
  use global
  implicit none
  real(dp), intent(in) :: log10E,B_l,tad,xi
  real(dp) :: resta
  real(dp) :: E_l,tsy,tac,tcoole

  E_l = 10.d0**log10E
  tsy = 1.d0/1.6d-3/(B_l*par)**2.d0/E_l   
  
  tcoole = min(tsy,tad) 
  tac = xi*E_l/(B_l*c*qe) 
  resta = tcoole - tac
  
  return
  
end Function resta

!/////////////////////////////////////////////////////////////////////////

