module radbr
  !Este módulo calcula la emisión bremsstrahlung 
  use global
  
  implicit none
  private
  public :: radbr_run
  
contains
  
  !============================================================================
  ! radbr_run
  !============================================================================
  
  subroutine radbr_run  
    use global
    implicit none
    real(dp) :: Ef_int,dEf,sum,sum_br,Phi_Br,dsigma_Br,u_br
    real(dp) :: gamma(mE_e),Ef(mE_f),cte_sigma,dq_Ef,q_Ef,dP_Ef
    !     mE - electrones, mE_f - fotones
    
    allocate( E_e(mE_e),N_e(mE_e),Lbr(mE_f) )
    
    write(*,*) '\n Calculate the (unabsorbed) rel. brem. emission in the SW'

    cte_sigma = 4.d0*alfa*re**2*Zq_sw**2 
    
    !     Ahora la distribución de e- con energía Ee(i): N_e(i)
    open (14,file='dist_e.dat')
    do i=1,mE_e
       read(14,*) E_e(i),N_e(i)
       E_e(i) = E_e(i)*eV    !   Pasamos la energía a ergios
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                        CÁLCULOS                                       C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías de fotones, Eps(i), con Eps_min del orden 
    ! del keV (si gamma 2 ~ 1MeV, el Brems es NR) 
    
    call vector_log(E_e(1)/1.d2,E_e(mE_e),Ef_int,Ef)
    gamma = E_e/mec2 ! gamma es adimensional
    sum_br = 0.d0
    do i = 1,mE_f             ! Ef(i) --> fotones 
       dEf = Ef(i)*(Ef_int-1.d0)
       sum = 0.d0
       do j = 1,mE_e            ! N_e(j),gamma(j) --> electrones
          !Pedimos que Ef <= Ee por conservación de energía
          u_br = Ef(i)/E_e(j) ! energy ratio
          if (u_br > 1.d0) then
             dsigma_Br = 0.d0  ![cm**2/m_e c2]
          else
             ! Apantallado:
             !Phi_Br=(1.d0+(1.d0-u_br)**2-(2.d0/3.d0)*(1.d0-u_br))*&
             !     log(191.d0/Z**(1.d0/3.d0)) + 1.d0/9.d0*(1.d0-u_br)
             ! No apantallado (desnudo):
             Phi_Br=(1.d0+(1.d0-u_br)**2-(2.d0/3.d0)*(1.d0-u_br))*&
                  ( log(2.*gamma(j)*(1./u_br-1.)) -0.5d0)
             dsigma_Br = cte_sigma/Ef(i)*Phi_Br
             q_Ef = n_a*c*dsigma_Br
             dE_e = E_e(j)*(Eint_e-1.d0)
             dq_Ef = N_e(j)*q_ef*dE_e
             sum = dP_Ef + sum ! q = int dq_Eps
          endif
       enddo ! j
       ! To transform from photon emissivity to specific luminosity:
       Lbr(i) = Ef(i)*sum            ! [erg/seg/erg]
       sum_br = sum_br + Lbr(i)*dEf
    enddo ! i
    
    write(*,81) 'L_Br(tot)=',sum_br
81  format(A11,X,Es10.2)
    
    open(16,file='radbr.dat')
    do i=1,mE_f
       write(16,*) Ef(i)/eV,Ef(i)*Lbr(i)
    enddo
    close(16)
    
    deallocate( E_e,N_e,Lbr )
    
  END subroutine Radbr_run
  
  
  !========================================
  
  
END module Radbr


!========================================

