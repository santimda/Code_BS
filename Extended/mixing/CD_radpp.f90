module CD_radpp
  !Este módulo calcula la emisión pp ya corregida por abs gg
  use global
  
  implicit none
  private
  public :: CD_radpp_run

contains
  
  !============================================================================
  ! CD_radpp_run
  !============================================================================
  
  subroutine CD_radpp_run  
    use global
    Implicit none
    Real(dp) :: dummy,n_a(m_theta),Q_Ef,Q_phi,dx_pp,x_pp,qpi
    Real(dp) :: Ef_int,dEf,Ef(mE_f),sum_pp,sigma_pp,Epp,L1
    Real(dp) :: E_th= 1.22d9*eV, TeV = 1.d12*eV
    Real(dp) :: Epi_min,Epi_max,Epi_int,dE_pi,integral
    Integer :: nl,mE_pi,in
    ! m_theta - CD, mE_p - protons, mE_f - photons, mE_pi(variable) - pions
    Real(dp), dimension(:), allocatable :: Epi                      
    interface
       function F(a,b)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: F
         real(dp), intent(in) :: a,b
       end function F
    end interface
    
    allocate ( Lpp(mE_f),Lpp_tot(mE_f) )
    allocate ( B(m_theta),rho(m_theta),P(m_theta),v_t(m_theta),T(m_theta) )
    allocate ( B_par(m_theta),E_p(mE_p),N_ptot(m_theta,mE_p),N_p(mE_p) )
    
80    format(A7,X,Es10.2)

    !////////////////////////////////////////////////////////////////////////
    !
    ! Read thermodynamical quantites at each position
    
    open(11,file='termo_sw.dat')
    do l=1,m_theta
       read(11,*) nl,rho(l),P(l),v_t(l),T(l),B(l) 
!       n_a(l) = rho(l)/mp/mu_sw  
       n_a(l) = 4.*rho_ISM/mp/mu_sw 
    enddo
    close(11)
    B_par = B*par

    ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
    open (14,file='distp_sw.dat')
    do l=1,m_theta
       do i=1,mE_p
          read(14,*) dummy, E_p(i),N_ptot(l,i)
          E_p(i) = E_p(i)*eV    !   energy from eV to erg
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_p = E_p(2)/E_p(1)
        
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                  Calculate p-p emission                               C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
       
    Lpp_tot = 0.d0

    ! Build logarithmic spaced photon energy vector, Ef(i), with Ef_min = 10 MeV
    
    call vector_log(10.d6*eV,E_p(mE_p),Ef_int,Ef)

    ! At each position we calculate the p-p emission; then we sum all.
    sum_pp=0.d0
    Epi_max = E_p(mE_p)

    DO l=1,m_theta          ! l
       do i=1,mE_p
          N_p(i)=N_ptot(l,i)
       enddo

       mE_pi = mE_f    ! Initial value
       do i = 1,mE_f           ! Ef(i) 
          dEf = Ef(i)*(Ef_int-1.d0)
          ! We use a different formalism in 2 energy ranges: 

          if (Ef(i) < 0.1d0*TeV) then    ! Ef < 100 GeV
             allocate(Epi(mE_pi))
             Epi_min = Ef(i)+(m_pic2)**2.d0/4.d0/Ef(i)
             Epi_min = Epi_min*(1.d0+1.d0/mE_pi)
             call vector_log(Epi_min,Epi_max,Epi_int,Epi)
             integral = 0.d0
             do j=1,mE_pi       ! E_pi
                dE_pi = Epi(j)*(Epi_int-1.d0)
                Epp = Epi(j)/0.17d0 + mpc2
                in = int(log(Epp/Emin_p)/log(Eint_p))+1
                if (in > mE_pi .or. Epp > E_p(mE_p)) then
                   goto 0979
                endif

                if (Epp < E_th) then ! E_th= 1.22 GeV
                   sigma_PP = 0.d0
                else            
                   L1 = log(Epp/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
                   ! Add factor for better fit at low energies 
                   sigma_PP = sigma_PP*(1.d0-(E_th/Epp)**4.)**2.  
                endif
                qpi = N_p(in)*sigma_pp/0.17d0 ! 1/erg/seg
                integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-(m_pic2)**2)
             enddo              ! E_pi
0979         continue
             
             Q_Ef = 2.3d0*c*n_a(l)*integral/1.3d0! ñ=1.1 corresponds to alfa=2 (Kelner2006?)
             ! SED = Ef * L_Ef =  Ef^2 * Q_Ef 
             Lpp(i) = Ef(i)*Q_Ef
             sum_PP = sum_PP + Lpp(i)*dEf
             
             deallocate(Epi)
             mE_pi = MAX(mE_pi-1,32) !  Lower down the # of points for speed
             
          elseif (Ef(i) > 0.1d0*TeV) then  ! Ef > 100 GeV > Emin_p
             integral=0
             do j=1,mE_p
                if (Ef(i) < E_p(j)) then
                   dE_p = E_p(j)*(Eint_p-1.d0)
                   x_pp = Ef(i)/E_p(j)
                   dx_pp = Ef(i)*dE_p/(E_p(j)**2)
                   L1 = log(E_p(j)/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2.)*1.d-27
                   ! Add factor for better fitting at low E:
                   sigma_PP = sigma_PP*(1.d0-(E_th/E_p(j))**4)**2.
                   integral = integral + sigma_PP*N_p(j)*F(x_pp,L1)*dx_pp/x_pp
                endif
             enddo
             Q_phi = c*n_a(l)*integral
             Lpp(i) = Ef(i)*Q_phi          ! [erg/seg]

             sum_PP = sum_PP + Lpp(i)*dEf              
          endif          
                     
          Lpp_tot(i) = Lpp_tot(i) + Lpp(i)
                             
       enddo    ! i
         
    ENDDO       ! l

    write(*,80) 'L_pp=',sum_pp*mphi*2.

    !---------------------------------------
    !
    !   Save data

    open(16,file='radpp_sw.dat')
    do i=1,mE_f
       write(16,*) Ef(i)/eV,Ef(i)*Lpp_tot(i)
    enddo
    close(16)

    deallocate ( Lpp,Lpp_tot)
    deallocate ( rho,P,v_t,T,B,B_par )
    deallocate ( E_p,N_ptot,N_p )


  END subroutine CD_radpp_run
  
  
  !========================================
  
  
END module CD_radpp


!========================================

Function F(x,L)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x,L
  real(dp) :: B,k,beta,F,L2

  L2 = L**2.d0

  B = 1.30d0 + 0.14d0*L + 0.011d0*L2
  beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L2)
  k = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L2)
  F = B*log(x)/x*((1.d0-x**beta)/(1.d0+k*x**beta*(1.d0-x**beta)))**4.d0* &
       (1.d0/log(x)-(4.d0*beta*x**beta)/(1.d0-x**beta)- &
       (4.d0*k*beta*x**beta*(1.d0-2.d0*x**beta)/(1.d0+k*x**beta*(1.d0-x**beta))))

  return
end Function F

