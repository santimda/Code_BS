module CD_radpp

  !============================================================================
  ! This module calculates the intrinsic and absorbed p-p emission using as an 
  ! input the proton distribution in distp(number).dat. The integrated 
  ! luminosity (unabs and abs) for the linear emitter is stored in radpp.dat.
  !============================================================================

  use global
  implicit none
  private
  public :: CD_radpp_run

  contains

  !============================================================================
  ! CD_radpp_run (START)
  !============================================================================

  subroutine CD_radpp_run  
    use global
    implicit none
    real(dp) :: dummy,n_a(m_theta),Q_Ef,Q_phi,dx_pp,x_pp,qpi
    real(dp) :: Ef_int,dEf,Ef(mE_f),sum_pp,Epp,L1
    real(dp) :: Epi_min,Epi_max,Epi_int,dE_pi,integral
    real(dp), parameter :: K_pi=0.17d0, n_pp=0.86d0, kappa_pp=K_pi*n_pp! n_pp=1.1 for alfa=2, n_pp=0.86 for alpha=2.5
    integer :: nl,mE_pi,i_p
    ! m_theta - CD, mE_p - protons, mE_f - photons, mE_pi(variable) - pions
    real(dp), dimension(:), allocatable :: Epi                      
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

    write(*,*) 'Calculate the p-p emission in the SW'

    !------------------------------------------------------------------
    ! Read thermodynamical quantites at each position

    open(11,file='termo_sw.dat')
    do l=1,m_theta
      read(11,*) nl,rho(l),P(l),v_t(l),T(l),B(l) 
      n_a(l) = rho(l)/mp/mu_sw  
    end do
    close(11)
    B_par = B*par

    ! Number of electrons at r(l) with energy E_p(i): N_ptot(l,i) 
    open (14,file='distp_sw.dat')
    do l=1,m_theta
      do i=1,mE_p
        read(14,*) dummy, E_p(i),N_ptot(l,i)
        E_p(i) = E_p(i)*eV    !   energy from eV to erg
      end do
      read(14,*)
    end do
    close(14) 
    Eint_p = E_p(2)/E_p(1)


    !=======================================================================
    !  Begin p-p calculations 
    !=======================================================================

    ! Build logarithmic spaced photon energy vector, Ef(i), with Ef_min = 10 MeV

    call vector_log(10.d6*eV,E_p(mE_p),Ef_int,Ef)

    ! At each position we calculate the p-p emission; then we sum all.
    Epi_max = E_p(mE_p)

    Lpp_tot = 0.d0
    sum_pp=0.d0

    do l=1,m_theta          ! l
      do i=1,mE_p
        N_p(i)=N_ptot(l,i)
      end do

      mE_pi = mE_f    ! Initial value
      do i = 1,mE_f           ! Ef(i) 
        dEf = Ef(i)*(Ef_int-1.d0)
        ! We use a different formalism in 2 energy ranges: 

        if (Ef(i) < 0.01d0*TeV) then     ! Ep < 100 GeV -> Ef < 0.01 TeV
          allocate(Epi(mE_pi))
          Epi_min = Ef(i)+(m_pic2)**2/4.d0/Ef(i)
          Epi_min = Epi_min*1.005d0     ! Factor to avoid numerical errors very close to Epi_min
          call vector_log(Epi_min,Epi_max,Epi_int,Epi)
          integral = 0.d0
          do j=1,mE_pi       ! E_pi
            dE_pi = Epi(j)*(Epi_int-1.d0)
            Epp = n_pp*Epi(j)/kappa_pp + mpc2       ! Eq.75 Kelner+2006
            i_p = int(log(Epp/Emin_p)/log(Eint_p))+1
            if (i_p <= mE_p .and. Epp <= E_p(mE_p) .and. Epp > E_th) then 
              L1 = log(Epp/TeV)    
              qpi = n_pp*N_p(i_p)*sigma_pp(L1,Epp)/K_pi ! Eq.77 Kelner+2006 [1/erg/s]
              integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-(m_pic2)**2) 
            end if
          end do              ! E_pi
          Q_Ef = 2.d0*c*n_a(l)*integral*z_pp   !Eq.78 Kelner+2006
          Lpp(i) = Ef(i)*Q_Ef
          deallocate(Epi)
          mE_pi = MAX(mE_pi-1,32) ! Diminish the number of points 


        elseif (Ef(i) > 0.01d0*TeV) then  ! Ef > 100 GeV > Emin_p
          integral=0
          do j=1,mE_p
            if (Ef(i) < E_p(j)) then
              dE_p = E_p(j)*(Eint_p-1.d0)
              x_pp = Ef(i)/E_p(j)
              dx_pp = Ef(i)*dE_p/(E_p(j)**2)
              L1 = log(E_p(j)/TeV)      ! L = ln(Ep/1TeV)
              integral = integral + sigma_pp(L1,E_p(j))*N_p(j)*F(x_pp,L1)*dx_pp/x_pp
            end if
          end do
          Q_phi = c*n_a(l)*integral*z_pp
          Lpp(i) = Ef(i)*Q_phi          ! [erg/seg]

        end if          

        sum_pp = sum_pp + Lpp(i)*dEf

        Lpp_tot(i) = Lpp_tot(i) + Lpp(i)

      end do      ! i

    end do        ! l

    write(*,80) 'L_pp=',sum_pp*mphi*2.d0

    !---------------------------------------
    !
    !   Save data

    open(16,file='radpp_sw.dat')
    do i=1,mE_f
      write(16,*) Ef(i)/eV,Ef(i)*Lpp_tot(i)
    end do
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

  L2 = L**2
  B = 1.30d0 + 0.14d0*L + 0.011d0*L2
  beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L2)
  k = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L2)
  F = B*log(x)/x*((1.d0-x**beta)/(1.d0+k*x**beta*(1.d0-x**beta)))**4.d0* &
    (1.d0/log(x)-(4.d0*beta*x**beta)/(1.d0-x**beta)- &
    (4.d0*k*beta*x**beta*(1.d0-2.d0*x**beta)/(1.d0+k*x**beta*(1.d0-x**beta))))

  return
end Function F

