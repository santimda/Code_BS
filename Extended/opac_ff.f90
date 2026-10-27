module opac_ff
  !=========================================================================
  ! f-f absorption and emission of the unshocked stellar wind and of the
  ! ionized ISM (Stromgren sphere).
  ! Wind: n_e = gamma_w*n_i (Leitherer et al. 1995). ISM: n_ISM = total particle density.
  !=========================================================================
  use global
  use numerical_utils, only: vector_log
  use gaunt_factor, only: gaunt_ff_calc
  use map_freqs, only: maps, nu_index
  implicit none
  private
  public :: opac_ff_run

contains

  subroutine opac_ff_run
    real(dp) :: dr, log_nu_int, Mdot_eff, A_ff
    real(dp) :: nu_S, B_ff_D2, cte_Tau, convert, Delta_r
    real(dp) :: E_st_0, E_st_int, Q_ion, alpha_B, R_str, vol, vol_str
    real(dp) :: ne_ISM, ni_ISM, nH_ISM, nu_int_ff
    real(dp) :: P_trans, M_tot, V_trans
    ! stellar wind
    real(dp), allocatable :: dif_R(:), r_ff(:), g_ff_w(:), cte_alpha(:), kappa_nu(:), S(:)
    real(dp), allocatable :: alphaff(:,:)
    ! wind radio maps
    real(dp), allocatable :: r_map(:), S_ff_cm2(:), S_map(:), phi_map(:)
    ! Stromgren sphere and ionized ISM
    real(dp), allocatable :: E_st(:), g_ff_ISM(:), p_ff(:), factor_opac(:), Lnu(:)
    integer  :: im, i_nu, iu
    character(len=32) :: fmt

    allocate( nu(mE_ff), nu_ISM(mE_ff) )
    allocate( dif_R(m_theta), r_ff(mr_ff), g_ff_w(mE_ff), cte_alpha(mr_ff), kappa_nu(mE_ff), S(mE_ff) )
    allocate( alphaff(mr_ff,mE_ff) )
    allocate( r_map(mr_map), S_ff_cm2(mr_map), S_map(mr_map), phi_map(mphi) )
    allocate( E_st(mE_ph), g_ff_ISM(mE_ff), p_ff(mE_ff), factor_opac(mE_ff), Lnu(mE_ff) )

    !-----------------------------------------------------------------
    ! Shock position and cell size, and frequency grid

    open(newunit=iu, file=output_file('CD.dat'))
    read(iu,*) (theta(l), r_norm(l), dif_R(l), dummy, dummy, l=1,m_theta)
    close(iu)

    open(newunit=iu, file=output_file('nu.dat'))
    read(iu,*) nu
    close(iu)
    log_nu_int = log(nu(2)/nu(1))

    !==================================================================
    !
    !   UNSHOCKED STELLAR WIND
    !
    !==================================================================

    print*, 'Calculate the f-f abs. and emi. coefficients in the stellar wind'

    ! Radial grid (spherically symmetric wind, r >> R_st)
    dr = (5.d0*R0 - 10.d0*R_st)/real(mr_ff-1,dp)
    r_ff = 10.d0*R_st + dr*[( real(l,dp), l=0,mr_ff-1 )]

    open(newunit=iu, file=output_file('r_ff.dat'))
    write(iu,'(es25.16e3)') r_ff
    close(iu)

    call gaunt_ff_calc(nu, mE_ff, T_wind, Zq_w, g_ff_w)    ! T_wind: constant T in the wind

    ! Opacity alpha_ff(r,nu) [1/cm] (Rybicki & Lightman, Eq. 5.18b)
    Mdot_eff = Mdot/sqrt(filling)              ! effective mass-loss rate (clumping)
    A_ff = Mdot_eff/(4.d0*pi*v_w*mu_w*mp)      ! r >> R_st: v_w = v_infty
    cte_alpha = 3.7d8*Zq_w**2*gamma_w*(A_ff/r_ff**2)**2/sqrt(T_wind)
    kappa_nu  = nu**(-3)*(1.d0-exp(-h*nu/(k*T_wind)))*g_ff_w
    alphaff   = spread(cte_alpha,2,mE_ff)*spread(kappa_nu,1,mr_ff)

    write(fmt,'(a,i0,a)') '(', mE_ff, 'es17.8e3)'
    open(newunit=iu, file=output_file('windalphaff.dat'))
    write(iu,fmt) ((alphaff(l,i), i=1,mE_ff), l=1,mr_ff)
    close(iu)
    deallocate( alphaff )                      ! not needed anymore (largest array)

    ! Radio flux density of the wind [Jy] (Wright & Barlow 1975, Eq. 8)
    S = 23.2d0/(distance/kpc)**2 * ( ((Mdot_eff/M_sun_yr)/(v_w/1.d5)/mu_w)**2 &
        * nu*gamma_w*g_ff_w*Zq_w**2 )**(2.d0/3.d0)
    S = S*merge(exp(-h*nu/(k*T_wind)), 1.d0, nu > 1.d-3*k*T_wind/h)   ! ad hoc high-frequency cutoff

    open(newunit=iu, file=output_file('windemiff.dat'))
    write(iu,'(2es17.8e3)') (nu(i), S(i)*1.d3, i=1,mE_ff)               ! [mJy]
    close(iu)

    !------------------------------------------------------------------
    ! Wind radio maps at the frequencies of maps(:)

    convert = 206265.d0/distance             ! linear to angular size ["]
    phi_map = delta_phi*[( real(j,dp), j=0,mphi-1 )]

    do im = 1, size(maps)
       i_nu = nu_index(maps(im)%nu, nu(1), log_nu_int, mE_ff)
       nu_S = nu(i_nu)
       B_ff_D2 = 2.d0*nu_S**2*k*T_wind/c**2/distance**2*1.d26
       cte_Tau = 3.7d8*h*Zq_w**2*g_ff_w(i_nu)/(k*T_wind**1.5d0*nu_S**2) * gamma_w*A_ff**2*pi/2.d0

       call S_ff_map(B_ff_D2, cte_Tau, R_st, R0, S_ff_cm2, r_map)
       Delta_r = r_map(2)-r_map(1)
       S_map = S_ff_cm2*r_map*Delta_r*0.5d0*pi/mphi          ! flux of each cell [mJy]

       open(newunit=iu, file=output_file('S_ff_map_'//trim(maps(im)%label)//'.dat'))
       write(iu,'(3es17.8e3)') ((r_map(i)*cos(phi_map(j))*convert, r_map(i)*sin(phi_map(j))*convert, &
                                 S_map(i), j=1,mphi), i=1,mr_map)
       close(iu)
    end do

    !==================================================================
    !
    !   IONIZED UNSHOCKED ISM (STROMGREN SPHERE)
    !
    !==================================================================

    print*, 'Calculate the emission from the ionized unshocked ISM'

    ! Ionizing photon rate of the star (E > 13.6 eV)
    E_st_0 = 2.7d0*k*T_st
    call vector_log(13.6d0*eV, 20.d0*E_st_0, E_st_int, E_st)
    Q_ion = 4.d0*pi*R_st**2*c*(E_st_int-1.d0) * sum([( nBB(E_st(i),T_st)*E_st(i), i=1,mE_ph )])

    ! Densities of the ionized ISM
    ne_ISM = rho_ISM/(mu_e*mp)
    ni_ISM = rho_ISM/(mu_i*mp)
    nH_ISM = X_ab*rho_ISM/mp

    ! Stromgren radius (case B)
    alpha_B = 2.56d-13*(T_ISM/1.d4)**(-0.83d0)
    R_str = ( 3.d0*Q_ion/(4.d0*pi*ne_ISM*nH_ISM*alpha_B) )**(1.d0/3.d0)
    write(*,'(a, Es10.2)') ' Q_ion(ph/s) = ', Q_ion
    write(*,'(a, F7.3)')   ' R_str(pc) = ', R_str/pc
    write(*,'(a, F7.2)')   ' R_str/R0 = ', R_str/R0

    !----------------------------------------------------------
    ! f-f emission of the sphere (Olnon 1975)

    call vector_log(nu_min_ff, nu_max_ff, nu_int_ff, nu_ISM)
    call gaunt_ff_calc(nu_ISM, mE_ff, T_ISM, Zq_ISM, g_ff_ISM)

    ! Optical depth through the center of the sphere (Olnon 1975, Eqs. 4 and 10)
    p_ff = 2.d0*ne_ISM*ni_ISM*(R_str/pc) * 8.235d-2*T_ISM**(-1.35d0)*(nu_ISM/1.d9)**(-2.1d0)

    ! Olnon (1975), Eq. 11
    where (p_ff > 1.d-4)
       factor_opac = 1.d0 - (2.d0/p_ff**2)*(1.d0-(p_ff+1.d0)*exp(-p_ff))
    elsewhere
       factor_opac = 2.d0*p_ff/3.d0
    end where
    Lnu = 8.d0*pi**2*k*T_ISM*R_str**2/c**2 * nu_ISM**2 * factor_opac      ! [erg/s/Hz]

    open(newunit=iu, file=output_file('gaunt_ISM.dat'))
    write(iu,'(2es17.8e3)') (h*nu_ISM(i)/(k*T_ISM), g_ff_ISM(i), i=1,mE_ff)
    close(iu)

    open(newunit=iu, file=output_file('radff.dat'))
    write(iu,'(2es17.8e3)') (h*nu_ISM(i)/eV, Lnu(i)/h, i=1,mE_ff)            ! E [eV], L_E [1/s]
    close(iu)

    open(newunit=iu, file=output_file('flux_ISM.dat'))
    write(iu,'(2es17.8e3)') (nu_ISM(i), Lnu(i)*1.d26/(4.d0*pi*distance**2), i=1,mE_ff)   ! [mJy]
    close(iu)

    !----------------------------------------------------------
    ! Volume of the HII region outside the BS

    vol = sum( pi*(R_theta*sin(theta))**2*dif_R )
    vol_str = (4.d0/3.d0)*pi*R_str**3 - vol

    open(newunit=iu, file=output_file('Stromgren.dat'))
    write(iu,'(2es25.16e3)') R_str, vol_str
    close(iu)

    !----------------------------------------------------------
    ! Velocity transferred to the ISM by the stellar radiation
    ! V = F*t/m, F = L_st/c, t = R_str/v_st (extreme case)

    P_trans = (L_st/c)*(R_str/v_st)
    M_tot = rho_ISM*(4.d0/3.d0)*pi*R_str**3
    V_trans = P_trans/M_tot
    write(*,'(a, Es10.2)') ' V_trans[cm/s] =', V_trans

    deallocate( nu, nu_ISM )
    deallocate( dif_R, r_ff, g_ff_w, cte_alpha, kappa_nu, S )
    deallocate( r_map, S_ff_cm2, S_map, phi_map )
    deallocate( E_st, g_ff_ISM, p_ff, factor_opac, Lnu )

  END subroutine opac_ff_run

  !=========================================================================
  ! Brightness of the wind at impact parameter r_map (Wright & Barlow 1975, Eq. 5)
  !=========================================================================

  subroutine S_ff_map(B_ff_D2, cte_Tau, r_i, r_f, S_ff_cm2, r_map)
    real(dp), intent(in)  :: B_ff_D2, cte_Tau, r_i, r_f
    real(dp), intent(out) :: S_ff_cm2(mr_map), r_map(mr_map)
    real(dp) :: R_c, x_c(mr_map), Tau_max(mr_map)
    integer  :: ir

    R_c = 2.d0*r_i ! NOT SURE!
    r_map = r_i + (r_f-r_i)/real(mr_map-1,dp)*[( real(ir,dp), ir=0,mr_map-1 )]
    x_c = r_map/R_c

    where (r_map >= R_c)
       Tau_max = cte_Tau/r_map**3
    elsewhere
       Tau_max = cte_Tau/r_map**3*( 0.5d0 - acos(x_c)/pi - sqrt(1.d0-x_c**2)*x_c/pi )
    end where
    S_ff_cm2 = B_ff_D2*(1.d0-exp(-Tau_max))

  END subroutine S_ff_map

END module opac_ff