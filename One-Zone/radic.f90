module radic
  ! This module calculates the intrinsic IC emission 
  use global
  
  implicit none
  private
  public :: radic_run

contains
  
  !============================================================================
  ! radic_run
  !============================================================================
  
  subroutine radic_run(T,Kappa,Ef,L_ic) 
    use global
    implicit none
    real(dp), intent(in) :: T,Kappa
    real(dp), intent(out) :: Ef(mE_f),L_ic(mE_f)
    real(dp) :: Eps(mE_f),dEps,Eps_int,Ef_int,gamma(mE_e),x0,z,t_K
    real(dp) :: Ndot,sum,sum_ic,cte_Ndot,T_ad
    real(dp) :: T_ad_ang
    real(dp), parameter :: aniso = 1.d1 ! if = 1.d0, then anisotropic IC-star
    allocate( N_e(mE_e),E_e(mE_e) )

80  format(A10,X,Es10.2)

    T_ad = k*T/mec2  ! adimensional temperature

    !------------------------------------------------------------------------
    !
    !  Reading information here  
    !

    ! N_e(i) = number of e- with energy E_e(i)
    open (14,file='dist_e.dat')
    do i=1,mE_e
       read(14,*) E_e(i),N_e(i)
       E_e(i) = E_e(i)*eV    !   [erg]
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma = E_e/mec2

    !///////////////////////////////////////////////////////////////////////
    !
    !                     Inverse Compton 
    !
    !///////////////////////////////////////////////////////////////////////
    
    sum_ic = 0.d0
 
    cte_Ndot = 2.d0*re**2*mec2**3*Kappa*T_ad**2/(pi*h_bar**3*c**2)

    if (T == ANISO*T_st) then
       T_ad_ang = T_ad*(1.-cos(inclination))  !anisotropic IC
    endif

    call vector_log(k*T,E_e(mE_e),Ef_int,Ef)
    Eps = Ef/mec2
    Eps_int = Ef_int
    do i = 1,mE_f
       dEps=Eps(i)*(Eps_int-1.d0)
       sum = 0.d0
       do j = 1,mE_e
          z = Eps(i)/gamma(j)
          if (z < 1.d0) then 
             if (T == ANISO*T_st) then
                t_K = 4.*gamma(j)*T_ad_ang  !anisotropic
             else
                t_K = 4.*gamma(j)*T_ad      !isotropic
             endif
             x0 = z/( (1.-z)*t_K )
             if (T == ANISO*T_st) then
                Ndot = (cte_Ndot/gamma(j)**2) * ( (z**2/(2.*(1.-z)))*F_1(x0) + F_2(x0) ) !anisotropic
             else
                Ndot = (cte_Ndot/gamma(j)**2) * ( (z**2/(2.*(1.-z)))*F_3(x0) + F_4(x0) ) !isotropic
             endif
             P_eps = Eps(i)*Ndot !N^dot = nf*c*dsigmaIC
             dE_e = E_e(j)*(Eint_e-1.d0)
             dP_eps = N_e(j)*P_eps*dE_e
          else
             dP_eps = 0.d0
          endif
          sum = dP_eps + sum
       enddo
       L_ic(i) = sum 
       !   P=[m_e c**2/s], lo pasamos a [erg/s]. Eps lo pasamos a eV.
       sum_ic = sum_ic + L_ic(i)*dEps*mec2
    enddo
    
    write(*,80) 'L_IC =',sum_ic
    
    deallocate( N_e,E_e )
    
  END subroutine radic_run
   
  
!==============================

  
END module radic



