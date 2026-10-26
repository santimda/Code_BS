module CD_radic
  use global
  use controls, only: verbose
  use numerical_utils, only: vector_log, unit_vector_diff, safe_acos_dot
  implicit none
  private

  public :: CD_radic_st_run, CD_radic_IR_run

  ! radic_st uses the stellar photon field (anisotropic, constant T), while
  ! radic_IR uses the IR one (isotropic, theta-dependent T).
  ! IC emission: formalism of Khangulyan et al. (2014), with the functions
  ! F_1, F_2 (anisotropic) and F_3, F_4 (isotropic) defined in global.

  ! Tables of F_1..F_4, built once and interpolated linearly in log(x);
  ! outside [x_lo, x_hi] the functions are evaluated directly
  real(dp), parameter :: x_lo = 1.d-6, x_hi = 1.d3
  integer,  parameter :: n_tab = 20001         ! relative interpolation error ~1e-7
  real(dp), allocatable, save :: F_tab(:,:)    ! (n_tab, 4)
  real(dp), save :: dlx_tab

contains

  !============================================================================
  ! CD_radic_st_run
  !============================================================================
  ! Intrinsic IC emission with the (anisotropic) stellar photon field, using
  ! the electron distribution diste_<system>_<shock>.dat and the positions of
  ! the rotated CD (CDrot.dat, both halves: "up" and "down").
  ! Output (output/): radic_st_<shock>.dat with E_f [eV] and E_f*L_E [erg/s]
  !============================================================================

  subroutine CD_radic_st_run(LIC_int, shock)
    !$ use omp_lib
    character(LEN=*), intent(in) :: shock
    real(dp), intent(out) :: LIC_int
    integer,  parameter :: m_l = 2*m_theta  ! CD positions, both sides
    integer  :: l, l2, i, iu
    real(dp) :: ang, Kappa_st, L_unabs, L_unabs_tot, Ef_int
    real(dp) :: X_emi(3), X_Ee(3), X_st(3)
    real(dp) :: LIC(mE_f), N_e(mE_e)          ! private copies per OpenMP thread (small)
    real(dp), allocatable :: Ef(:), LIC_tot(:), E_e(:), N_etot(:,:), x(:), y(:), z(:)

    allocate( Ef(mE_f), LIC_tot(mE_f), E_e(mE_e), N_etot(m_theta,mE_e), x(m_l), y(m_l), z(m_l) )

    if (verbose) print *, 'Calculate the IC(st) emission'

    ! set omp threads number
    !$ call omp_set_num_threads(nth)

    !------------------------------------------------------------------------
    !  Read the rotated CD positions and the electron distribution

    open(newunit=iu, file=output_file('CDrot.dat'))
    read(iu,*) (x(l), y(l), z(l), l=1,m_l)     ! "up" and "down"
    close(iu)
    x = R0*x
    y = R0*y
    z = R0*z

    call read_diste(shock, E_e, N_etot)

    call init_F_table()

    !------------------------------------------------------------------------
    !  IC emission of one linear emitter. The stellar photon field is
    !  anisotropic and its temperature is constant. l2 runs over both
    !  halves; l is the position in the upper half (same N_e in the mirror).

    call vector_log(k*T_st, E_e(mE_e), Ef_int, Ef)

    X_st  = (/ 0.d0, 0.d0, 0.d0 /)

    LIC_tot = 0.d0
    L_unabs_tot = 0.d0

    !$omp parallel do default(shared) schedule(dynamic) reduction(+:LIC_tot,L_unabs_tot) &
    !$omp private(l2, l, X_emi, X_Ee, ang, Kappa_st, N_e, LIC, L_unabs)
    DO l2 = 1, m_l
      l = merge(l2, l2-m_theta, l2 <= m_theta)
      X_emi = (/ x(l2), y(l2), z(l2) /)
      Kappa_st = ( R_st/(2.d0*norm2(X_emi-X_st)) )**2   ! dilution at the emitter

      call unit_vector_diff(X_emi, X_st, 3, X_Ee)  ! (X_emi - X_st)/|X_emi - X_st|
      ang = safe_acos_dot(e_obs, X_Ee)             ! scattering angle: LOS (global.f90) vs X_Ee

      N_e = N_etot(l,:)
      call IC_ani(E_e, N_e, Ef, Kappa_st, T_st_ad, ang, LIC, L_unabs)

      LIC_tot = LIC_tot + LIC
      L_unabs_tot = L_unabs_tot + L_unabs
    ENDDO
    !$omp end parallel do

    write(*,'(A9,X,Es10.2)') 'L_IC_st =', L_unabs_tot
    LIC_int = L_unabs_tot     ! To calculate the total LIC in main.f90

    open(newunit=iu, file=output_file('radic_st_'//shock//'.dat'))
    write(iu,'(2es17.8e3)') (Ef(i)/eV, Ef(i)*LIC_tot(i), i=1,mE_f)
    close(iu)

    deallocate( Ef, LIC_tot, E_e, N_etot, x, y, z )

  END subroutine CD_radic_st_run

  !============================================================================
  ! CD_radic_IR_run
  !============================================================================
  ! Intrinsic IC emission with the (isotropic) IR photon field, whose
  ! temperature and dilution depend on theta (T_IR.dat, from CD_diste).
  ! Output (output/): radic_ir_<shock>.dat with E_f [eV] and E_f*L_E [erg/s]
  ! (one linear emitter)
  !============================================================================

  subroutine CD_radic_IR_run(shock)
    !$ use omp_lib
    character(LEN=*), intent(in) :: shock
    integer  :: l, i, iu
    real(dp) :: L_unabs, L_unabs_tot, Ef_int
    real(dp) :: LIC(mE_f), N_e(mE_e)          ! private copies per OpenMP thread (small)
    real(dp), allocatable :: Ef(:), LIC_tot(:), E_e(:), N_etot(:,:), T_IR_ad(:), Kappa_IR(:)

    allocate( Ef(mE_f), LIC_tot(mE_f), E_e(mE_e), N_etot(m_theta,mE_e), T_IR_ad(m_theta), Kappa_IR(m_theta) )

    if (verbose) print *, 'Calculate the IC(IR) emission'

    ! set omp threads number
    !$ call omp_set_num_threads(nth)

    !------------------------------------------------------------------------
    !  Read the IR field and the electron distribution

    open(newunit=iu, file=output_file('T_IR.dat'))
    read(iu,*) (dummy, dummy, T_IR_ad(l), Kappa_IR(l), l=1,m_theta)
    close(iu)

    call read_diste(shock, E_e, N_etot)

    call init_F_table()

    !------------------------------------------------------------------------
    !  IC emission of one linear emitter

    call vector_log(T_IR_ad(1)*mec2, E_e(mE_e), Ef_int, Ef)

    LIC_tot = 0.d0
    L_unabs_tot = 0.d0

    !$omp parallel do default(shared) schedule(dynamic) reduction(+:LIC_tot,L_unabs_tot) &
    !$omp private(l, N_e, LIC, L_unabs)
    DO l = 1, m_theta
      N_e = N_etot(l,:)
      call IC_iso(E_e, N_e, Ef, Kappa_IR(l), T_IR_ad(l), LIC, L_unabs)
      LIC_tot = LIC_tot + LIC
      L_unabs_tot = L_unabs_tot + L_unabs
    ENDDO
    !$omp end parallel do

    write(*,'(A9,X,Es10.2)') 'L_IC_IR =', L_unabs_tot*mphi*4.d0   ! isotropic: mphi emitters, 4 quadrants

    open(newunit=iu, file=output_file('radic_ir_'//shock//'.dat'))
    write(iu,'(2es17.8e3)') (Ef(i)/eV, Ef(i)*LIC_tot(i), i=1,mE_f)
    close(iu)

    deallocate( Ef, LIC_tot, E_e, N_etot, T_IR_ad, Kappa_IR )

  END subroutine CD_radic_IR_run


  !---------------------------------------------------------------------------!
  !                                                                           !
  !                           SUBROUTINES                                     !
  !                                                                           !
  !---------------------------------------------------------------------------!

  subroutine read_diste(shock, E_e, N_etot)
    ! Electron distribution N_etot(l,i) at position l and energy E_e(i) [erg]
    ! (one block of mE_e lines per position, separated by a blank line)
    character(LEN=*), intent(in) :: shock
    real(dp), intent(out) :: E_e(mE_e), N_etot(m_theta,mE_e)
    integer :: l, i, iu

    open(newunit=iu, file=output_file('diste_'//system//'_'//shock//'.dat'))
    do l = 1, m_theta
      read(iu,*) (dummy, E_e(i), N_etot(l,i), i=1,mE_e)
      read(iu,*)
    enddo
    close(iu)
    E_e = E_e*eV        ! energy from eV to erg
  end subroutine read_diste

  !============================================================================

  SUBROUTINE IC_ani(Ee,Ne,Ef,Kappa,T_ad,ang,L_IC,L_unabs)
    ! Anisotropic IC emission from interactions with a BB photon field
    ! (Khangulyan et al. 2014): t = 2*gamma*T*(1-cos(ang))
    real(dp), intent(in)  :: Ee(mE_e),Ne(mE_e),Ef(mE_f),Kappa,T_ad,ang
    real(dp), intent(out) :: L_IC(mE_f),L_unabs
    real(dp) :: gamm(mE_e),NedEe(mE_e),Eps(mE_f),Eps_int,T_ad_ang2,cte_Ndot_ani
    real(dp) :: z_ic,x0_K,Ndot,suma
    integer  :: i, j

    T_ad_ang2 = 2.d0*T_ad*(1.d0-cos(ang))
    gamm  = Ee/mec2
    NedEe = Ne*Ee*(Ee(2)/Ee(1)-1.d0)      ! N(E)*dE, to speed up the calculation
    Eps   = Ef/mec2
    Eps_int = Eps(2)/Eps(1)
    cte_Ndot_ani = cte_Ndot*T_ad**2*Kappa

    do i = 1, mE_f              ! Ef
      suma = 0.d0
      do j = 1, mE_e            ! Ee
        z_ic = Eps(i)/gamm(j)
        if (z_ic < 1.d0) then
          x0_K = z_ic/( (1.d0-z_ic)*gamm(j)*T_ad_ang2 )
          Ndot = (cte_Ndot_ani/gamm(j)**2) * ( (z_ic**2/(2.d0*(1.d0-z_ic)))*F_tab_eval(1,x0_K) + F_tab_eval(2,x0_K) )
          suma = suma + Eps(i)*Ndot*NedEe(j)   ! N^dot = nf*c*dsigma_IC
        endif
      enddo
      L_IC(i) = suma            ! specific luminosity (plot E*L_IC(E))
    enddo
    L_unabs = sum(L_IC*Eps)*(Eps_int-1.d0)*mec2

  END SUBROUTINE IC_ani

  !============================================================================

  SUBROUTINE IC_iso(Ee,Ne,Ef,Kappa,T_ad,L_IC,L_unabs)
    ! Isotropic IC emission from interactions with a BB photon field
    ! (Khangulyan et al. 2014): t = 4*gamma*T
    real(dp), intent(in)  :: Ee(mE_e),Ne(mE_e),Ef(mE_f),Kappa,T_ad
    real(dp), intent(out) :: L_IC(mE_f),L_unabs
    real(dp) :: gamm(mE_e),NedEe(mE_e),Eps(mE_f),Eps_int,T_ad4,cte_Ndot_iso
    real(dp) :: z_ic,x0_K,Ndot,suma
    integer  :: i, j

    T_ad4 = 4.d0*T_ad
    gamm  = Ee/mec2
    NedEe = Ne*Ee*(Ee(2)/Ee(1)-1.d0)      ! N(E)*dE, to speed up the calculation
    Eps   = Ef/mec2
    Eps_int = Eps(2)/Eps(1)
    cte_Ndot_iso = cte_Ndot*T_ad**2*Kappa

    do i = 1, mE_f              ! Ef
      suma = 0.d0
      do j = 1, mE_e            ! Ee
        z_ic = Eps(i)/gamm(j)
        if (z_ic < 1.d0) then
          x0_K = z_ic/( (1.d0-z_ic)*gamm(j)*T_ad4 )
          Ndot = (cte_Ndot_iso/gamm(j)**2) * ( (z_ic**2/(2.d0*(1.d0-z_ic)))*F_tab_eval(3,x0_K) + F_tab_eval(4,x0_K) )
          suma = suma + Eps(i)*Ndot*NedEe(j)   ! N^dot = nf*c*dsigma_IC
        endif
      enddo
      L_IC(i) = suma            ! specific luminosity (plot E*L_IC(E))
    enddo
    L_unabs = sum(L_IC*Eps)*(Eps_int-1.d0)*mec2

  END SUBROUTINE IC_iso

  !============================================================================

  subroutine init_F_table()
    ! Tabulates F_1..F_4 on a logarithmic grid in [x_lo, x_hi] (only once)
    integer  :: i
    real(dp) :: x

    if (allocated(F_tab)) return
    allocate( F_tab(n_tab,4) )
    dlx_tab = log(x_hi/x_lo)/real(n_tab-1,dp)
    do i = 1, n_tab
      x = x_lo*exp(dlx_tab*real(i-1,dp))
      F_tab(i,:) = (/ F_1(x), F_2(x), F_3(x), F_4(x) /)
    enddo
  end subroutine init_F_table

  !============================================================================

  function F_tab_eval(n, x)
    ! F_n(x) interpolated linearly in log(x); direct evaluation outside the table
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: F_tab_eval, u, w
    integer  :: i

    if (x <= x_lo .or. x >= x_hi) then
      select case (n)
      case (1)
        F_tab_eval = F_1(x)
      case (2)
        F_tab_eval = F_2(x)
      case (3)
        F_tab_eval = F_3(x)
      case default
        F_tab_eval = F_4(x)
      end select
    else
      u = log(x/x_lo)/dlx_tab
      i = min(int(u) + 1, n_tab - 1)
      w = u - real(i-1,dp)
      F_tab_eval = (1.d0-w)*F_tab(i,n) + w*F_tab(i+1,n)
    endif
  end function F_tab_eval

END module CD_radic