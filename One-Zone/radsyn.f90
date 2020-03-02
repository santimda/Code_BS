module radsyn
  ! This module calculates the synchrotron emission
  use global
  
  implicit none
  private
  public :: radsyn_run

contains
  
  !============================================================================
  ! radsyn_run
  !============================================================================
  
  subroutine radsyn_run  
    use global
    !$ use omp_lib
    implicit none
    ! mE_e - electrons, m_nu - emitted photons 
    real(dp) :: Eps(m_nu),dEps,Eps_int, gamma(mE_e)
    real(dp) :: factor,nuc,nuc_max,nuc_min,sum,sum_sy
    real(dp) :: ne,nu_RT,x_sy
    real(dp) :: log_nu_int
    allocate ( N_e(mE_e),E_e(mE_e) )
    allocate ( nu(m_nu),Lsy(m_nu) )

80  format(A14,X,Es10.2)

    write(*,*) '\n Calculate the (unabsorbed) sync. emission in the SW'

    ! set omp threads number
    !$ call omp_set_num_threads(nth) 

    !  Number of electrons with energy between E_e(i) and E_e(i)+dE_e: N_e(i) 
    open (14,file='dist_e.dat')
    do i=1,mE_e
       read(14,*) E_e(i),N_e(i)
       E_e(i) = E_e(i)*eV   !   energy from eV to erg
    end do
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma = E_e/mec2
    
    !-----------------------------------------------------------------------!
    !                                                                       !
    !              Calculate Sycnhrotron emission                           !
    !                                                                       !
    !-----------------------------------------------------------------------!

    factor = h/mec2

    nuc_max = 4.22d6*B_par*gamma(mE_e)**2. *9.d0    ! Gives an extra
    nuc_min = 10.d6 ! 4.22d6*B_par*gamma(1)**2*10.d0 ! 10 MHz

    ! Build logarithmic-spaced energy vector for emitted photons, Eps(i):
    call vector_log(factor*nuc_min,factor*nuc_max,Eps_int,Eps)
    nu = Eps/factor ! [nu] = Hz
    log_nu_int = log(nu(2)/nu(1))
    open (88,file='nu.dat')
    do i=1,m_nu
       write(88,*) nu(i)
    end do
    close(88)

    ! Para cada posición voy a calcular la emisión sinc y luego sumo todas
    sum_sy=0.d0
    ne = n_a*mu_sw/mu_e_sw
    nu_RT = 20.d0*ne/B  ! frecuencia característica para efecto R-T
    write(*,80) 'nu_RT[GHz] =',nu_RT/1.d9  ! frecuencia característica para efecto R-T
    
    ! Barremos todos los rangos de energía de los fotones, Epsilon (Eps), 
    ! entre Eps_min y Eps_max, dados por nu_c(gamma_min) y nu_c(gamma_max)
    do i = 1,m_nu           ! Eps(i), nu(i)
       dEps = Eps(i)*(Eps_int-1.d0)
       sum = 0.d0
       ! Usamos E como variable de integración, P_eps = P(E,epsilon),
       ! Lsy = P(epsilon) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps],
       ! x = nu/nu_c va entre 0 y 10
       
       !$omp parallel do private(nuc,dP_Eps,P_Eps,dE_e) &
       !$omp reduction(+ : sum)
       do j = 1,mE_e         ! Ne(j), nuc(j), gamma(j)
          nuc = 4.22d6*B_par*gamma(j)**2.
             x_sy = nu(i)/nuc
             if (x_sy > 10.d0) then
             dP_Eps = 0.d0
          else
             ! Calculamos P_eps [erg/s/Hz] y pasamos a [1/s] dividiendo por h
             P_Eps = ( 4.39d-22*B_par*x_sy**(1.d0/3.d0)*exp(-x_sy) )/h
             dE_e = E_e(j)*(Eint_e - 1.d0)
             dP_Eps = N_e(j)*P_eps*dE_e
          endif
          sum = dP_Eps + sum
       end do
       !$omp end parallel do 
       
       Lsy(i) = sum*mec2*exp(-nu_RT/nu(i)) ! Lsy corrected for R-T effect

       sum_sy = sum_sy + Lsy(i)*dEps
       
    end do  ! Eps(i),nu(i)
    
    
    write(*,80) 'L_syn(unabs)=',sum_sy

    !   Ya tengo la emisión total no absorbida:
    
    open(17,file='radsyn.dat')
    do i=1,m_nu
       write(17,*) Eps(i)*mec2/eV,Eps(i)*Lsy(i) ! [erg/s]
    end do
    close(17)

    deallocate ( nu,Lsy )
    deallocate ( N_e,E_e )
    
  END subroutine radsyn_run
   
  
!==============================

  
END module radsyn



