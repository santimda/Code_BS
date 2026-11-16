module shock_thermal
  !==========================================================================
  ! Thermal (f-f) emission of the shocks (Martinez et al. 2023, Sect. 3.2.2):
  ! RS, and hot (Eq. 14) and isothermal (Eq. 15) layers of a radiative FS.
  ! Emissivity: Rybicki & Lightman, Eq. 5.14b
  !==========================================================================
  use global
  use gaunt_factor, only: gaunt_ff_calc
  use map_freqs, only: maps, nu_index
  implicit none
  private
  public :: shock_thermal_run

contains

  subroutine shock_thermal_run(shock, radiative)
    character(LEN=*), intent(in) :: shock
    logical, intent(in) :: radiative
    integer, parameter :: m_l = 2*m_theta        ! CD positions, both halves
    real(dp) :: nu_int, log_nu_int, flux, zq, L_v_tot, L_HL_tot
    real(dp), allocatable :: nu_th(:), dnu(:), nuL_v_tot(:), nuL_HL_tot(:), w_avg(:)
    real(dp), allocatable :: dV_ring(:), L_inj(:), T_min(:), rho_sh(:), T_sh(:), rho_IL_s(:)
    ! RS (shocked wind) or FS isothermal layer: "volume" component
    real(dp), allocatable :: ne_v(:), ni_v(:), T_v(:), V_v(:)
    ! FS hot layer
    real(dp), allocatable :: integ(:), P_rad(:)
    ! maps and diagnostics
    real(dp), allocatable :: S_avg(:), S_cum(:), alpha_HL(:), alpha_v(:)
    ! (theta, nu) arrays
    real(dp), allocatable :: g_v(:,:), nuL_v(:,:), g_HL(:,:), nuL_HL(:,:)
    real(dp), allocatable :: S_map(:,:)
    integer,  allocatable :: i_map(:)
    integer  :: im, n_map, iu, i_1G, i_3G

    call validate_shock(shock, 'shock_thermal_run')
    print*, 'Calculate the thermal emission from the '//shock

    allocate( nu_th(m_nu), dnu(m_nu), nuL_v_tot(m_nu), nuL_HL_tot(m_nu), w_avg(m_nu) )
    allocate( dV_ring(m_theta), L_inj(m_theta), T_min(m_theta), rho_sh(m_theta), T_sh(m_theta), &
              rho_IL_s(m_theta), ne_v(m_theta), ni_v(m_theta), T_v(m_theta), V_v(m_theta), &
              integ(m_theta), P_rad(m_theta), S_avg(m_theta), S_cum(m_theta), &
              alpha_HL(m_theta), alpha_v(m_theta) )
    allocate( g_v(m_theta,m_nu), nuL_v(m_theta,m_nu), g_HL(m_theta,m_nu), nuL_HL(m_theta,m_nu) )

    if (shock == 'RS') then
      zq = Zq_sw
    else
      zq = Zq_ISM
    end if

    !------------------------------------------------------------
    ! Frequency grid

    open(newunit=iu, file=output_file('nu.dat'))
    read(iu,*) nu_th
    close(iu)
    nu_int = nu_th(2)/nu_th(1)
    log_nu_int = log(nu_int)
    dnu = nu_th*(nu_int-1.d0)
    flux = 1.d26/(4.d0*pi*distance**2)          ! luminosity to flux in mJy

    !------------------------------------------------------------
    ! Shock position, volume, injected power and thermodynamics

    open(newunit=iu, file=output_file('CD.dat'))
    read(iu,*) (theta(l), dummy, dummy, dummy, dummy, l=1,m_theta)
    close(iu)

    open(newunit=iu, file=output_file('Hnorm_'//shock//'.dat'))
    if (radiative) then
      read(iu,*) (dummy, dummy, dummy, dummy, dV_ring(l), dummy, l=1,m_theta)
    else
      read(iu,*) (dummy, dummy, dummy, dummy, dV_ring(l), l=1,m_theta)
    end if
    close(iu)

    open(newunit=iu, file=output_file('Linj_'//shock//'.dat'))
    read(iu,*) (dummy, dummy, L_inj(l), dummy, dummy, l=1,m_theta)    ! one linear emitter
    close(iu)

    open(newunit=iu, file=output_file('thermo_'//shock//'.dat'))
    read(iu,*) (dummy, rho_sh(l), dummy, dummy, T_sh(l), dummy, dummy, dummy, dummy, dummy, dummy, &
                rho_IL_s(l), l=1,m_theta)
    close(iu)

    if (radiative) then
      open(newunit=iu, file=output_file('T_IL_'//shock//'.dat'))
      read(iu,*) (dummy, T_min(l), l=1,m_theta)
      close(iu)
    end if

    !------------------------------------------------------------
    ! Emission of the shocked wind (RS) or of the isothermal layer (FS)

    if (radiative) then
      ne_v = rho_IL_s/(mu_e*mp)
      ni_v = rho_IL_s/(mu_i*mp)
      T_v  = T_min                             ! T of the IL (CD_thermo)
      V_v  = frac_H*dV_ring                    ! IL width reduced by eta_H
    else
      ne_v = rho_sh/(mu_e*mp)
      ni_v = rho_sh/(mu_i*mp)
      T_v  = T_sh
      V_v  = dV_ring
    end if

    do l = 1, m_theta
      call gaunt_ff_calc(nu_th, m_nu, T_v(l), zq, g_v(l,:))
    end do

    ! nu*L_nu of each ring [erg/s]
    nuL_v = spread(6.8d-38*zq**2*ne_v*ni_v/sqrt(T_v)*V_v, 2, m_nu) &
          * exp(-h*spread(nu_th,1,m_theta)/(k*spread(T_v,2,m_nu))) * g_v * spread(nu_th,1,m_theta)
    nuL_v_tot = sum(nuL_v, dim=1)
    L_v_tot = sum(nuL_v_tot)*(nu_int-1.d0)

    !------------------------------------------------------------
    ! Hot layer of the FS: L_ff,HL = (Lambda_ff/Lambda_tot)*(1-T_IL/T_HL)*L_perp

    nuL_HL = 0.d0
    if (radiative) then
      do l = 1, m_theta
        call gaunt_ff_calc(nu_th, m_nu, T_sh(l), zq, g_HL(l,:))
      end do
      nuL_HL = exp(-h*spread(nu_th,1,m_theta)/(k*spread(T_sh,2,m_nu))) * g_HL * spread(nu_th,1,m_theta)
      integ = matmul(nuL_HL, dnu/nu_th)                     ! int j_nu dnu (up to a constant)
      P_rad = merge( lambda_ff(T_sh)/lambda_cool(T_sh)*(1.d0-T_min/T_sh)*L_inj*n_emit, &
                     0.d0, T_sh > T_min )                   ! whole ring: n_emit linear emitters
      nuL_HL = nuL_HL*spread(P_rad/integ, 2, m_nu)
    end if
    nuL_HL_tot = sum(nuL_HL, dim=1)
    L_HL_tot = sum(nuL_HL_tot)*(nu_int-1.d0)

    !------------------------------------------------------------
    ! SEDs and fluxes

    if (radiative) then
      open(newunit=iu, file=output_file(shock//'_ff_'//system//'.dat'))
      write(iu,'(2es17.8e3)') (h*nu_th(i)/eV, nuL_HL_tot(i), i=1,m_nu)             ! E [eV], nu*L_nu [erg/s]
      close(iu)
      open(newunit=iu, file=output_file('S_ff_'//system//'_'//shock//'.dat'))
      write(iu,'(2es17.8e3)') (nu_th(i), nuL_HL_tot(i)/nu_th(i)*flux, i=1,m_nu)     ! nu [Hz], S [mJy]
      close(iu)
    end if

    open(newunit=iu, file=output_file('ionized_ff_'//system//'_'//shock//'.dat'))
    write(iu,'(3es17.8e3)') (h*nu_th(i)/eV, nuL_v_tot(i), nuL_v_tot(i)/(h*nu_th(i)), i=1,m_nu)  ! E [eV], nu*L_nu, L_E
    close(iu)
    open(newunit=iu, file=output_file('S_ionized_'//system//'_'//shock//'.dat'))
    write(iu,'(2es17.8e3)') (nu_th(i), nuL_v_tot(i)/nu_th(i)*flux, i=1,m_nu)         ! nu [Hz], S [mJy]
    close(iu)

    n_map = size(maps)
    allocate( i_map(n_map), S_map(m_l,n_map) )
    i_map = nu_index(maps%nu, nu_th(1), log_nu_int, m_nu)

    if (radiative) then
      print*, '\n f-f emission from the hot layer:'
      do im = 1, n_map
        write(*,'(A,A,A,X,Es10.3)') ' S_', trim(maps(im)%label), '[mJy]=', nuL_HL_tot(i_map(im))/nu_th(i_map(im))*flux
      end do
      write(*,'(A,X,Es10.2)') ' L_HL_tot =', L_HL_tot
      write(*,'(A,X,Es10.2)') ' L_inj*(1-T_IL/T_HL)/L_HL_tot =', sum(L_inj*(1.d0-T_min/T_sh))*n_emit/L_HL_tot
      print*, '\n f-f emission from the isothermal layer:'
    else
      print*, '\n f-f emission from the '//shock//':'
    end if
    do im = 1, n_map
      write(*,'(A,A,A,X,Es10.3)') ' S_', trim(maps(im)%label), '[mJy]=', nuL_v_tot(i_map(im))/nu_th(i_map(im))*flux
    end do
    write(*,'(A,X,Es10.2)') ' L_ion_tot =', L_v_tot

    ! In an adiabatic shock the thermal emission must be a small fraction of the injected power
    if (.not. radiative .and. L_v_tot > sum(L_inj)*n_emit) &
      write(*,*) 'WARNING (shock_thermal): thermal luminosity of the '//shock// &
                 ' larger than the injected power; the shock is not adiabatic'

    !------------------------------------------------------------
    ! Emission maps at the frequencies of maps(:) [mJy]

    if (radiative) then
      S_map(1:m_theta,:) = nuL_HL(:,i_map)/n_emit * spread(flux/nu_th(i_map), 1, m_theta)
      S_map(m_theta+1:,:) = S_map(1:m_theta,:)
      do im = 1, n_map
        open(newunit=iu, file=output_file('CD_hot_layer_'//trim(maps(im)%label)//'.dat'))
        write(iu,'(es17.8e3)') S_map(:,im)
        close(iu)
      end do
    end if

    S_map(1:m_theta,:) = nuL_v(:,i_map)/n_emit * spread(flux/nu_th(i_map), 1, m_theta)
    S_map(m_theta+1:,:) = S_map(1:m_theta,:)
    do im = 1, n_map
      open(newunit=iu, file=output_file('CD_th_'//shock//'_'//trim(maps(im)%label)//'.dat'))
      write(iu,'(es17.8e3)') S_map(:,im)
      close(iu)
    end do

    ! Average flux between nu_av_min and nu_av_max
    w_avg = merge(flux*(nu_int-1.d0)/(nu_av_max-nu_av_min), 0.d0, &
                  nu_th >= nu_av_min .and. nu_th <= nu_av_max) / n_emit

    if (radiative) then
      S_avg = matmul(nuL_HL, w_avg)
      open(newunit=iu, file=output_file('CD_hot_radiomap_avg.dat'))
      write(iu,'(es17.8e3)') S_avg, S_avg                    ! both halves
      close(iu)
    end if

    S_avg = matmul(nuL_v, w_avg)
    open(newunit=iu, file=output_file('CD_th_radiomap_avg_'//shock//'.dat'))
    write(iu,'(es17.8e3)') S_avg, S_avg                      ! both halves
    close(iu)

    !------------------------------------------------------------
    ! Diagnostics

    ! Cumulative flux at 3 GHz along the shock (FS)
    if (radiative) then
      i_3G = nu_index(3.d9, nu_th(1), log_nu_int, m_nu)
      S_cum = [( sum(nuL_v(1:l,i_3G)), l=1,m_theta )]/nu_th(i_3G)*flux
      open(newunit=iu, file=output_file('L_3GHz.dat'))
      write(iu,'(2es17.8e3)') (theta(l), S_cum(l), l=1,m_theta)
      close(iu)
    end if

    ! f-f absorption coefficient at 1.4 GHz [1/cm]: HL and RS/IL
    i_1G = nu_index(1.42d9, nu_th(1), log_nu_int, m_nu)
    alpha_v = 3.7d8*zq**2*ne_v*ni_v/sqrt(T_v) * (1.d0-exp(-h*nu_th(i_1G)/(k*T_v)))/nu_th(i_1G)**3 * g_v(:,i_1G)
    alpha_HL = 0.d0
    if (radiative) alpha_HL = 3.7d8*zq**2*(rho_sh/(mu_e*mp))*(rho_sh/(mu_i*mp))/sqrt(T_sh) &
                              * (1.d0-exp(-h*nu_th(i_1G)/(k*T_sh)))/nu_th(i_1G)**3 * g_HL(:,i_1G)
    open(newunit=iu, file=output_file('Tau_'//shock//'.dat'))
    write(iu,'(3es17.8e3)') (theta(l), alpha_HL(l), alpha_v(l), l=1,m_theta)
    close(iu)

    deallocate( i_map, S_map )
    deallocate( nu_th, dnu, nuL_v_tot, nuL_HL_tot, w_avg, dV_ring, L_inj, T_min, rho_sh, T_sh, &
                rho_IL_s, ne_v, ni_v, T_v, V_v, integ, P_rad, S_avg, S_cum, alpha_HL, alpha_v )
    deallocate( g_v, nuL_v, g_HL, nuL_HL )

  END subroutine shock_thermal_run

END module shock_thermal