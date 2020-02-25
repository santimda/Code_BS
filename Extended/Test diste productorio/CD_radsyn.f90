module CD_radsyn
  ! This module calculates the intrinsic synchrotron emission 
  use global
  
  implicit none
  private
  public :: CD_radsyn_run
  
  ! We save the total emission (i.e., integrated in theta) in radsyn.dat. 
  ! The position dependent emission, Eps x L_sy(l,Eps_i), is stored in emisyn.dat.
  ! The latter is useful for correcting for absorption and producing emission maps.

contains
  
  !============================================================================
  ! CD_radsyn_run
  !============================================================================
  
  subroutine CD_radsyn_run  
    use global
    !$ use omp_lib
    implicit none
    ! m_theta - CD, mE_e - electrons, m_nu - emitted photons 
    integer :: il
    real(dp) :: Eps(m_nu),dEps,Eps_int, gamma(mE_e), dummy
    real(dp) :: factor,nuc,nuc_max,nuc_min,sum,sum_sy,x_sy
    real(dp) :: nu_RT,nuc_l
    allocate ( N_e(mE_e),E_e(mE_e),N_etot(m_theta,mE_e) )
    allocate ( Lsy_tot(m_nu),nu(m_nu),Lsy(m_theta,m_nu) )
    Lsy_tot = 0.d0
  
80  format(A19,X,Es10.2)

    ! set omp threads number
    !$ call omp_set_num_threads(nth) 

    !////////////////////////////////////////////////////////////////////////
    !
    ! Read shock position (x,y) and thermodynamical quantites at each position
    
    open(11,file='termo_sw.dat')
    do l=1,m_theta
       read(11,*) il,B(l) ! todo en CGS
    enddo
    close(11)
    B_par = B*par

    open(12,file='CD.dat')
    do l=1,m_theta
       read(12,*) il,z(l) ! la trayectoria en cm
    enddo
    close(12)

    ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
    open (14,file='diste_sw.dat')
    do l=1,m_theta
       do i=1,mE_e
          read(14,*) dummy, E_e(i),N_etot(l,i)
          E_e(i) = E_e(i)*eV    !   energy from eV to erg
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma = E_e/mec2
    
    !-----------------------------------------------------------------------!
    !                                                                       !
    !              Calculate Sycnhrotron emission                           !
    !                                                                       !
    !-----------------------------------------------------------------------!

    factor = h/mec2
    ! Use the same frequency range for all position
    nuc_max = 4.22d6*B_par(1)*gamma(mE_e)**2. *2.d0  ! A bit arbitrary
    nuc_min = 10.d6  !4.22d6*B_par(m_theta)*gamma(1)**2 / 1.d2  ! 10 MHz
    
    ! Build logarithmic-spaced energy vector for emitted photons, Eps(i):
    call vector_log(factor*nuc_min,10.d0*factor*nuc_max,Eps_int,Eps)
    nu = Eps/factor ! [nu] = Hz

    open(66,file='nu.dat')
    do i=1,m_nu
       write(66,*) nu(i)
    enddo
    close(66)

    ! At each position we calculate the synchrotron emission; then we sum all.
    sum_sy=0.d0
    DO l=1,m_theta
       do i=1,mE_e
          N_e(i)=N_etot(l,i)
       enddo
       nu_RT = 1.d3/B(l)   ! characteristic frequency for R-T effect 
       nuc_l = 4.22d6*B_par(l) ! to speed things up in the j loop
       ! Loop in the photon energy (Eps), given by nu 
       do i = 1,m_nu           ! Eps(i), nu(i)
          dEps = Eps(i)*(Eps_int-1.d0)
          sum = 0.d0
          ! Lsy(Eps) = int_Emin^Emax N(E_e) P(E_e,eps) dE [erg/s/eps]

          !$omp parallel do private(nuc,dP_Eps,P_Eps,dE_e) &
          !$omp reduction(+ : sum)
          do j = 1,mE_e         ! Ne(j), nuc(j), gamma(j)
             nuc = nuc_l*gamma(j)**2.
             x_sy = nu(i)/nuc
             if (x_sy > 10.d0) then
                dP_Eps = 0.d0
             else
                ! P_eps [erg/s/Hz]; divide by h to transform units
                P_Eps = ( 4.39d-22*B_par(l)*x_sy**(1.d0/3.d0)*exp(-x_sy) )/h
                dE_e = E_e(j)*(Eint_e - 1.d0)
                dP_Eps = N_e(j)*P_eps*dE_e
             endif
             sum = dP_Eps + sum
          enddo    
          !$omp end parallel do 
    
          Lsy(l,i) = sum*mec2*exp(-nu_RT/nu(i))     !
          Lsy_tot(i) = Lsy_tot(i) + Eps(i)*Lsy(l,i) ! [erg/s]
          
          sum_sy = sum_sy + Lsy(l,i)*dEps

       enddo  ! Eps(i),nu(i)
       
    ENDDO     ! l
    
    write(*,80) 'L_syn(unabsorbed)=',sum_sy

    !---------------------------------------
    !
    !   Save data

    ! Total SED
    open(17,file='radsyn_sw.dat')
    do i=1,m_nu
       write(17,*) Eps(i)*mec2/eV,Lsy_tot(i)
    enddo
    close(17)
    
    deallocate ( Lsy_tot,nu,Lsy )
    deallocate ( N_e,E_e,N_etot )

  END subroutine CD_radsyn_run
  
  
!==============================

  
END module CD_radsyn



