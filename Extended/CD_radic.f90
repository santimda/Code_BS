module CD_radic
  use global
  implicit none

  public :: CD_radic_st_run,CD_radic_IR_run

  ! radic_st uses the stellar photon field, while radic_IR uses the IR 
  ! (theta-dependent) one.

  contains

  !============================================================================
  ! CD_radic_st_run (start)
  !============================================================================
  ! This module calculates the intrinsic IC emission using as an 
  ! input the electron distribution in diste_<shock>.dat and the  
  ! (anisotropic) stellar photon field
  !============================================================================

  subroutine CD_radic_st_run(LIC_int,shock)  
    character(LEN=*), intent(in) :: shock
    real(dp),intent(out) :: LIC_int
    real(dp) :: ang,L_unabs,L_unabs_tot,Ef(mE_f),Ef_int
    integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
    integer :: l2

    allocate ( LIC(mE_f),LIC_tot(mE_f) )
    allocate ( E_e(mE_e),N_etot(m_theta,mE_e),N_e(mE_e) )
    allocate ( x(m_l),y(m_l),z(m_l) )

    80  format(A9,X,Es10.2)

    LIC = 0.d0 
    LIC_tot = 0.d0
    L_unabs_tot = 0.d0

    print *, 'Calculate the IC(st) emission'

    !------------------------------------------------------------------------
    !  Reading information here  
  
    open(11,file='CDrot.dat')
    do l=1,m_l     ! "up" and "down"
      read(11,*) x(l),y(l),z(l)
    enddo
    close(11)
    x = R0*x
    y = R0*y
    z = R0*z

    ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
    open (14,file='diste_'//shock//'.dat')
    do l=1,m_theta
      do i=1,mE_e
        read(14,*) dummy, E_e(i),N_etot(l,i)
        E_e(i) = E_e(i)*eV    !   energy from eV to erg
      enddo
      read(14,*)
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)

    !------------------------------------------------------------------------
    !  Calculate IC radiation next

    ! Define a logarithmically-spaced photon energy vector
    call vector_log( k*T_st,E_e(mE_e),Ef_int,Ef)

    X_obs = (/ sin(incli),cos(incli),0.d0 /)
    X_st = (/ 0.d0, 0.d0, 0.d0 /)

    ! We calculate the IC emission for a linear emitter. The "l" loop is for 
    ! the different positions of the emitter along the CD. The stellar photon
    ! field is anisotropic, and its temperature is constant.

    ! "l" is only for y > 0, the variable "l2" is used to take into account
    ! the "mirror" side with y < 0
    DO l2 = 1,m_l

      X_emi = (/ x(l2),y(l2),z(l2) /) 
      if (l2 <= m_theta) then
        l = l2
      else
        l = l2 - m_theta
      endif
      Kappa_st = ( R_st/(2.d0*R_theta(l)) )**2

      ! Calculate the versor: XEe/|XEe|
      call versorresta(X_emi,X_st,3,X_Ee)  ! from the 2nd to the 1st 
      ang = acos(DOT_PRODUCT(X_obs,X_Ee)) ! angle between Xobs and XEe

      ! Obtain the SED, L_IC(Ef), and the integrated luminosity, L_unabs, at "l":
      N_e = N_etot(l,:)      ! N_e = electron distribution at position "l"
      call IC_ani(E_e,N_e,Ef,Kappa_st,T_st_ad,ang,LIC,L_unabs) 

      L_unabs_tot = L_unabs_tot + L_unabs 

      ! Add the emission from the different positions to make one SED
      do i = 1,mE_f
        LIC_tot(i) = LIC_tot(i)+LIC(i)
      enddo

    ENDDO

    write(*,80) 'L_IC_st = ', L_unabs_tot
    LIC_int = L_unabs_tot ! To calculate the total LIC in main.f90

    ! Save the SED      
    open(16,file='radic_st_'//shock//'.dat')
    do i = 1,mE_f     
      write(16,*) Ef(i)/eV,Ef(i)*LIC_tot(i)
    enddo
    close(16)

    deallocate ( LIC,LIC_tot )
    deallocate ( E_e,N_etot,N_e )
    deallocate ( x,y,z )


  END subroutine CD_radic_st_run

  !============================================================================
  ! CD_radic_st_run (end)
  !============================================================================



  !============================================================================
  ! CD_radic_IR_run (start)
  !============================================================================
  ! This module calculates the intrinsic IC emission using as an 
  ! input the electron distribution in diste_<shock>.dat and the  
  ! (isotropic) IR photon field
  !============================================================================

  subroutine CD_radic_IR_run(shock)  
    character(LEN=*), intent(in) :: shock
    real(dp) :: L_unabs,L_unabs_tot,Ef(mE_f),Ef_int

    allocate ( LIC(mE_f),LIC_tot(mE_f) )
    allocate ( E_e(mE_e),N_etot(m_theta,mE_e),N_e(mE_e) )
    allocate ( T_IR_ad(m_theta) )

    80  format(A9,X,Es10.2)

    LIC = 0.d0 
    LIC_tot = 0.d0
    L_unabs_tot = 0.d0

    print *, 'Calculate the IC(IR) emission'
    !------------------------------------------------------------------------
    !
    !  Reading information here  
    !

    open(13,file='T_IR.dat')
    do l=1,m_theta
      read(13,*) dummy,dummy,T_IR_ad(l),Kappa_IR(l)
    enddo
    close(13)

    ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
    open (14,file='diste_'//shock//'.dat')
    do l=1,m_theta
      do i=1,mE_e
        read(14,*) dummy, E_e(i),N_etot(l,i)
        E_e(i) = E_e(i)*eV    !   energy from eV to erg
      enddo
      read(14,*)
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)

    !------------------------------------------------------------------------
    !
    !  Calculate IC radiation next
    !

    ! Define a logarithmically-spaced photon energy vector
    call vector_log( T_IR_ad(1)*mec2,E_e(mE_e),Ef_int,Ef)

    X_obs = (/ sin(incli),cos(incli),0.d0 /)
    X_st = (/ 0.d0, 0.d0, 0.d0 /)

    ! We calculate the IC emission for a linear emitter. The "l" loop is for 
    ! the different positions of the emitter along the CD. The IR photon
    ! field is isotropic, and its temperature depends on "l".

    DO l = 1,m_theta
      do i = 1,mE_e
        N_e(i)=N_etot(l,i)
      enddo

      call IC_iso(E_e,N_e,Ef,Kappa_IR(l),T_IR_ad(l),LIC,L_unabs) 

      L_unabs_tot = L_unabs_tot + L_unabs 

      ! Add the emission from the different positions to make one SED
      do i = 1,mE_f
        LIC_tot(i) = LIC_tot(i)+LIC(i)
      enddo

    ENDDO

    write(*,80) 'L_IC_IR = ', L_unabs_tot *mphi*4.d0 ! isotropic, 4 quadrants

    ! Save the SED      
    open(16,file='radic_ir_'//shock//'.dat')
    do i = 1,mE_f     
      write(16,*) Ef(i)/eV,Ef(i)*LIC_tot(i)
    enddo
    close(16)

    deallocate ( LIC,LIC_tot )
    deallocate ( E_e,N_etot,N_e )
    deallocate ( T_IR_ad )


  END subroutine CD_radic_IR_run

  !============================================================================
  ! CD_radic_IR_run (end)
  !============================================================================


END module CD_radic


!---------------------------------------------------------------------------!
!                                                                           !
!                           SUBROUTINES                                     !
!                                                                           !
!---------------------------------------------------------------------------!

SUBROUTINE IC_ani(Ee,Ne,Ef,Kappa,T_ad,ang,L_IC,L_unabs)
  ! Calculate anisotropic IC emission from interactions with BB photon field,
  ! using the formalism of Khangulyan+ 2014
  use global
  implicit none
  Real(dp), intent(in) :: Ee(mE_e),Ne(mE_e),Ef(mE_f),Kappa,T_ad,ang
  Real(dp), intent(out) :: L_IC(mE_f),L_unabs
  Real(dp) :: Eps(mE_f),dEps,Eps_int
  Real(dp) :: gamm(mE_e),Eeint,dEe,NedEe(mE_e)
  Real(dp) :: T_ad_ang,T_ad_ang2,x0_K,z_ic,t_ic
  Real(dp) :: Ndot,sum,cte_Ndot_ani
  !      Dimension :: Ee(mE_e),Ne(mE_e),Ef(mE_f),L_IC(mE_f)

  T_ad_ang = T_ad*(1.-cos(ang))
  T_ad_ang2 = T_ad_ang*2.d0

  Eeint = Ee(2)/Ee(1)
  gamm = Ee/mec2
  Eps = Ef/mec2
  Eps_int = Eps(2)/Eps(1) 

  cte_Ndot_ani = cte_Ndot*T_ad**2*Kappa

  L_unabs = 0.d0

  ! This helps to speed up future calculations
  do j = 1, mE_e 
    dEe = Ee(j)*(Eint_e-1.d0)
    NedEe(j) = Ne(j)*dEe
  end do 

  DO i = 1,mE_f              ! Ef
    dEps=Eps(i)*(Eps_int-1.d0)
    sum=0.d0

    do j = 1,mE_e          ! Ee
      z_ic = Eps(i)/gamm(j)
      if (z_ic < 1.d0) then
        t_ic = gamm(j)*T_ad_ang2
        x0_K = z_ic/( (1.-z_ic)*t_ic )
        Ndot = (cte_Ndot_ani/gamm(j)**2) * ( (z_ic**2./(2.*(1.-z_ic)))*F_1(x0_K) + F_2(x0_K) )
        P_eps = Eps(i)*Ndot !N^dot = nf*c*dsigmaIC
        dP_eps = P_eps*NedEe(j) ! P_eps*Ne(Ee)*dEe
        sum = dP_eps + sum
      endif
    enddo

    L_IC(i) = sum  ! Especific luminosity (plot E(j)*LIC(Ej))
    L_unabs = L_unabs + sum*dEps*mec2
  ENDDO

END SUBROUTINE IC_ani


!/////////////////////////////////////////////////////////////////////////


SUBROUTINE IC_iso(Ee,Ne,Ef,Kappa,T_ad,L_IC,L_unabs)
  ! Calculate isotropic IC emission from interactions with BB photon field,
  ! using the formalism of Khangulyan+ 2014
  use global
  implicit none
  Real(dp), intent(in) :: Ee,Ne,Ef,Kappa,T_ad
  Real(dp), intent(out) :: L_IC(mE_f),L_unabs
  Real(dp) :: Eps(mE_f),dEps,Eps_int
  Real(dp) :: Eeint,gamm(mE_e),dEe,NedEe(mE_e)
  Real(dp) :: x0_K,z_ic,t_ic,T_ad4
  Real(dp) :: Ndot,sum,cte_Ndot_iso
  Dimension :: Ee(mE_e),Ne(mE_e),Ef(mE_f)

  T_ad4 = T_ad * 4.d0
  Eeint = Ee(2)/Ee(1)
  gamm = E_e/mec2
  Eps = Ef/mec2
  Eps_int = Eps(2)/Eps(1) 

  cte_Ndot_iso = cte_Ndot*T_ad**2*Kappa

  L_unabs = 0.d0

  ! This helps to speed up future calculations
  do j = 1, mE_e
    dEe = Ee(j)*(Eint_e-1.d0)
    NedEe(j) = Ne(j)*dEe
  end do 

  DO i=1,mE_f              ! Ef
    dEps=Eps(i)*(Eps_int-1.d0)
    sum=0.d0

    do j = 1,mE_e          ! Ee
      z_ic = Eps(i)/gamm(j)
      if (z_ic < 1.d0) then
        t_ic = gamm(j)*T_ad4
        x0_K = z_ic/( (1.-z_ic)*t_ic )
        Ndot = (cte_Ndot_iso/gamm(j)**2) * ( (z_ic**2/(2.*(1.-z_ic)))*F_3(x0_K) + F_4(x0_K) )
        P_eps = Eps(i)*Ndot !N^dot = nf*c*dsigmaIC
        dEe = Ee(j)*(Eint_e-1.d0)
        dP_eps = P_eps*NedEe(j) ! P_eps*Ne(Ee)*dEe
        sum = dP_eps + sum
      endif
    enddo

    L_IC(i) = sum  ! Especific luminosity (plot E(j)*LIC(Ej))
    L_unabs = L_unabs + sum*dEps*mec2
  ENDDO

END SUBROUTINE IC_iso
