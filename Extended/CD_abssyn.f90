module CD_abssyn
   use global
   use controls, only: verbose
   use gaunt_factor, only: gaunt_ff_calc
   use map_freqs, only: maps, nu_index
   implicit none
   private
   public :: CD_abssyn_run

   !========================================================================
   ! Regions crossed by a ray (see region below)
   integer,  parameter :: REG_OUT  = 0     ! shocked wind, shocked ISM or unperturbed ISM: no absorption here
   integer,  parameter :: REG_WIND = 1     ! unshocked stellar wind (inside the RS, or the open tail)
   integer,  parameter :: REG_IL   = 2     ! isothermal (ionized) layer of the radiative FS
   integer,  parameter :: n_bis = 40       ! bisection iterations to locate a boundary
   integer,  parameter :: n_sub = 16       ! midpoint sub-steps for each crossing of the IL
   integer,  parameter :: n_step_max = 10000000
   real(dp), parameter :: tau_occ = 1.d3   ! optical depth of a photon occulted by the star

   ! Geometry of the shell, set in CD_abssyn_run and only read inside the parallel loop
   logical  :: IL                         ! the FS has an isothermal layer
   integer  :: n_abs                      ! wind and IL absorption computed for nu(1:n_abs)
   real(dp) :: th_min, th_max, dtheta     ! polar angle grid (measured from the apex, -x)
   real(dp) :: ds_in                      ! integration step for r < R_far
   real(dp) :: R_far                      ! largest radius of the shell [cm]
   real(dp) :: r_stop                     ! the ray marching stops at r > r_stop [cm]
   real(dp), allocatable :: r_CD(:)       ! CD radius [cm]
   real(dp), allocatable :: r_RS(:)       ! inner edge of the shocked wind (= RS) [cm]
   real(dp), allocatable :: r_ILo(:)      ! outer edge of the IL [cm]
   real(dp), allocatable :: nnk_IL(:,:)   ! alpha_ff of the IL, (n_abs, m_theta) [1/cm]

contains

   !========================================================================
   !        CD_abssyn_run
   !========================================================================
   ! Free-free absorption of the synchrotron emission of the shock(s) for the
   ! linear emitters at the current rotation angle phi (CDrot.dat).
   ! Absorbers:
   !  - unshocked stellar wind (tau_w). alpha_w(r,nu) = K_w(nu)/r**4, so the
   !    integral along each straight segment of the ray is analytic (G_w).
   !    The wind fills r < r_RS(theta) = r_CD - (radial width of the shocked
   !    wind) and the open tail theta > theta(m_theta).
   !  - isothermal layer of the radiative FS (tau_il), r_CD <= r < r_CD + eta_H*H
   !    (radial width), with n_e*n_i and alpha_ff evaluated at T_IL(theta).
   !    Used whenever the FS is radiative (T_IL_FS.dat exists), for any shock.
   !  - Stromgren sphere (tau_str), radius from Stromgren.dat.
   ! Wind and IL absorption are computed for nu < nu_max_ff (global.f90).
   ! Photons whose ray crosses the star are occulted (tau = tau_occ).
   ! The f-f absorption in the hot shocked gas is neglected given the large T.
   !
   ! Output (arguments): E_sy, photon energy [eV], and the SED summed over the
   ! linear emitters, unabsorbed (Lsy_sum) and absorbed (Lsy_sum_abs), as
   ! nu*L_nu [erg/s]. Files (output/): radio maps, average map, radsyn_abs*.dat
   ! and, if verbose, CD_tauff*.dat.
   !========================================================================

   subroutine CD_abssyn_run(shock, radiative, E_sy, Lsy_sum, Lsy_sum_abs)
   !$ use omp_lib
   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative
   real(dp), intent(out) :: E_sy(m_nu), Lsy_sum(m_nu), Lsy_sum_abs(m_nu)
   integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
   real(dp) :: log_nu_int, flux, convert, R_str, r_ff1, ne_ISM, ni_ISM, e_u(3)
   integer :: l, i, im, n_map, iu
   character(len=:), allocatable :: tag
   integer,  allocatable :: i_nu(:)                       ! index of each map frequency in freq(:)
   real(dp), allocatable :: pos(:,:)                      ! (x,y,z) of each linear emitter [cm]
   real(dp), allocatable :: freq(:)                       ! frequency grid (nu.dat) [Hz]
   real(dp), allocatable :: L_em(:,:)                     ! nu*L_nu of each linear emitter (position, nu)
   real(dp), allocatable :: L_FS(:,:)
   real(dp), allocatable :: tau_w(:,:), tau_il(:,:), tau_str(:), tau_tot(:,:)
   real(dp), allocatable :: gaunt(:), kappa_ISM(:), K_w(:)
   real(dp), allocatable :: u_map(:), v_map(:), emis(:,:), emis_abs(:,:), w_avg(:)

   allocate( pos(3,m_l), freq(m_nu), L_em(m_l,m_nu) )
   allocate( gaunt(m_nu), kappa_ISM(m_nu), K_w(m_nu) )
   allocate( tau_w(m_l,m_nu), tau_il(m_l,m_nu), tau_str(m_nu), tau_tot(m_l,m_nu) )

   call validate_shock(shock, 'CD_abssyn_run')
   if (verbose) print *, 'Calculate the f-f absorption of the synchrotron emission'

   convert = 206265/distance             ! image to angular size ["]
   flux = 1.d26/(4.d0*pi*(distance)**2)  ! luminosity to flux in mJy
   ! set omp threads number
   !$ call omp_set_num_threads(nth)

   !======================================================
   ! Geometry and frequency grid

   open(newunit=iu, file=output_file('CDrot.dat'))
   read(iu,*) pos                         ! one (x,y,z) per line, "up" and "down"
   close(iu)
   pos = R0*pos

   open(newunit=iu, file=output_file('nu.dat'))
   read(iu,*) freq
   close(iu)
   log_nu_int = log(freq(2)/freq(1))
   n_abs = count(freq < nu_max_ff)

   !======================================================
   ! Wind: alpha_w(r,nu) = K_w(nu)/r**4 (opac_ff.f90). K_w from the first row
   ! of the table (r_ff(1)), which is exact because T_wind is constant.

   open(newunit=iu, file=output_file('r_ff.dat'))
   read(iu,*) r_ff1
   close(iu)
   open(newunit=iu, file=output_file('windalphaff.dat'))
   read(iu,*) K_w                         ! first row: alpha_w(r_ff(1),nu)
   close(iu)
   K_w = K_w*r_ff1**4

   !======================================================
   ! Shell: CD, RS and IL (radial widths), alpha_ff of the IL

   call shell_geometry(shock, radiative, freq)

   !======================================================
   ! Opacities: stellar wind and IL (nu < nu_max_ff), Stromgren sphere

   call tau_calc(shock, pos, K_w, tau_w, tau_il)

   open(newunit=iu, file=output_file('Stromgren.dat'))
   read(iu,*) R_str
   close(iu)
   ne_ISM = rho_ISM/(mu_e*mp)            ! n_ISM: total particle density (rho_ISM = n_ISM*mu_ISM*mp)
   ni_ISM = rho_ISM/(mu_i*mp)
   call gaunt_ff_calc(freq, m_nu, T_ISM, Zq_ISM, gaunt)
   kappa_ISM = 3.7d8 * Zq_ISM**2/sqrt(T_ISM) * (1.d0-exp(-h*freq/(k*T_ISM)))/freq**3 * gaunt
   tau_str = kappa_ISM * ne_ISM*ni_ISM * R_str

   tau_tot = tau_w + tau_il + spread(tau_str, 1, m_l)

   if (shock == 'BOTH') then
      tag = ''
   else
      tag = shock
   end if

   if (verbose) then
      open(newunit=iu, file=output_file('CD_tauff'//tag//'.dat'))
      write(iu,'(i6,4es17.8e3)') ((l, freq(i), tau_w(l,i), tau_il(l,i), tau_tot(l,i), i=1,m_nu), l=1,m_l)
      close(iu)
   end if

   !/////////////////////////////////////////////////////////////////////
   ! Synchrotron emission of each linear emitter: emisyn_*.dat has
   ! E [eV] and nu*L_nu [erg/s] (already multiplied by E)

   select case (shock)
   case ('RS', 'FS')
      call read_emisyn(shock, E_sy, L_em)
   case ('BOTH')
      call read_emisyn('RS', E_sy, L_em)
      allocate( L_FS(m_l,m_nu) )
      call read_emisyn('FS', E_sy, L_FS)
      !L_em = L_em + L_FS                         ! FS synchrotron emission (disabled)
      deallocate( L_FS )
   end select
   L_em(m_theta+1:m_l,:) = L_em(1:m_theta,:)       ! symmetry ("down" side)

   !-----------------------------------------------------------------------
   ! Total unabsorbed and absorbed emission (sum over linear emitters)

   Lsy_sum     = sum(L_em, dim=1)
   Lsy_sum_abs = sum(L_em*exp(-tau_tot), dim=1)

   if (shock == 'BOTH') then
      open(newunit=iu, file=output_file('radsyn_abs.dat'))
   else
      open(newunit=iu, file=output_file('radsyn_abs_'//shock//'.dat'))
   end if
   write(iu,'(3es17.8e3)') (E_sy(i), Lsy_sum(i), Lsy_sum_abs(i), i=1,m_nu)
   close(iu)

   !-----------------------------------------------------------------------
   ! Emission maps at the frequencies selected in controls (maps(:)),
   ! projected in the plane of the sky: u along e_u = (cos i, sin i, 0) and
   ! v along z, both perpendicular to the LOS e_obs = (-sin i, cos i, 0)

   n_map = size(maps)
   allocate( i_nu(n_map), u_map(m_l), v_map(m_l), emis(m_l,max(n_map,1)), emis_abs(m_l,max(n_map,1)) )

   e_u   = [ cos(incli), sin(incli), 0.d0 ]
   u_map = matmul(e_u, pos)*convert
   v_map = pos(3,:)*convert

   i_nu = nu_index(maps%nu, freq(1), log_nu_int, m_nu)
   emis(:,1:n_map)     = L_em(:,i_nu) * spread(flux/freq(i_nu), 1, m_l)
   emis_abs(:,1:n_map) = emis(:,1:n_map) * exp(-tau_tot(:,i_nu))

   do im = 1, n_map
      open(newunit=iu, file=output_file('CD_sync_radiomap_'//tag//trim(maps(im)%label)//'.dat'))
      write(iu,'(4es17.8e3)') (u_map(l), v_map(l), emis(l,im), emis_abs(l,im), l=1,m_l)
      close(iu)
   end do

   !------------------------------------------------------------
   ! Average flux between nu_av_min and nu_av_max. In the log grid
   ! dnu = nu*(nu(2)/nu(1)-1), so (L_em/nu)*flux*dnu = L_em*flux*(nu(2)/nu(1)-1)

   allocate( w_avg(m_nu) )
   w_avg = merge(flux*(freq(2)/freq(1)-1.d0)/(nu_av_max-nu_av_min), 0.d0, &
                 (freq >= nu_av_min) .and. (freq <= nu_av_max))

   emis(:,1)     = matmul(L_em, w_avg)                  ! reuse emis(:,1) for the average map
   emis_abs(:,1) = matmul(L_em*exp(-tau_tot), w_avg)

   open(newunit=iu, file=output_file('CD_sync_radiomap_avg'//tag//'.dat'))
   write(iu,'(4es17.8e3)') (u_map(l), v_map(l), emis(l,1), emis_abs(l,1), l=1,m_l)
   close(iu)

   deallocate( i_nu, u_map, v_map, emis, emis_abs, w_avg )
   deallocate( pos, freq, L_em, gaunt, kappa_ISM, K_w )
   deallocate( tau_w, tau_il, tau_str, tau_tot )
   deallocate( r_CD, r_RS, r_ILo, nnk_IL )

   END subroutine CD_abssyn_run

   !========================================================================
   !        read_emisyn
   !========================================================================
   ! nu*L_nu [erg/s] of each position of the upper half (CD_radsyn.f90):
   ! one block of m_nu lines (l, E [eV], nu*L_nu) per position + a blank line

   subroutine read_emisyn(sh, E, Lv)
   character(LEN=*), intent(in) :: sh
   real(dp), intent(out) :: E(:), Lv(:,:)
   integer :: l, i, iu

   open(newunit=iu, file=output_file('emisyn_'//sh//'.dat'))
   do l = 1, m_theta
      read(iu,*) (dummy, E(i), Lv(l,i), i=1,m_nu)
      read(iu,*)
   end do
   close(iu)

   END subroutine read_emisyn

   !========================================================================
   !        shell_geometry
   !========================================================================
   ! Radial positions of the boundaries of the shell as a function of theta.
   ! H (H_<shock>.dat) is the width normal to the CD; the radial width is
   ! H*sqrt(1 + (dln r/dtheta)**2), since tan(psi) = (dr/dtheta)/r is the
   ! angle between the normal and the radial direction.
   !  - r_RS = r_CD - (radial width of the shocked wind), H_RS.dat
   !  - r_ILo = r_CD + eta_H*(radial width of the IL), H_FS.dat, radiative FS
   ! alpha_ff of the IL: n_e*n_i = (rho_IL/(mu_e*mp))*(rho_IL/(mu_i*mp)),
   ! at T_IL(theta) (T_IL_FS.dat) and Zq_ISM.
   !========================================================================

   subroutine shell_geometry(shock, radiative, freq)
   character(LEN=*), intent(in) :: shock
   logical,  intent(in) :: radiative
   real(dp), intent(in) :: freq(:)
   integer  :: l, j, iu
   logical  :: ex_RS, ex_IL
   real(dp), allocatable :: g_rad(:), dlnr(:), Hw(:), rhoIL(:), T_IL(:), gaunt(:)

   allocate( r_CD(m_theta), r_RS(m_theta), r_ILo(m_theta), nnk_IL(n_abs,m_theta) )
   allocate( g_rad(m_theta), dlnr(m_theta), Hw(m_theta), gaunt(size(freq)) )

   th_min = theta(1)
   th_max = theta(m_theta)
   dtheta = theta(2)-theta(1)                 ! equiespaciated

   r_CD = r_norm*R0

   ! radial width / normal width
   dlnr(2:m_theta-1) = (log(r_norm(3:m_theta)) - log(r_norm(1:m_theta-2)))/(2.d0*dtheta)
   dlnr(1)       = (log(r_norm(2)) - log(r_norm(1)))/dtheta
   dlnr(m_theta) = (log(r_norm(m_theta)) - log(r_norm(m_theta-1)))/dtheta
   g_rad = sqrt(1.d0 + dlnr**2)

   !------------------------------------------------------------
   ! Shocked wind (RS side)

   inquire(file=output_file('H_RS.dat'), exist=ex_RS)
   if (ex_RS) then
      open(newunit=iu, file=output_file('H_RS.dat'))
      do l = 1, m_theta
         read(iu,*) dummy, Hw(l)
      end do
      close(iu)
      r_RS = max(r_CD - Hw*g_rad, 0.1d0*r_CD)
   else
      print *, 'WARNING (CD_abssyn): H_RS.dat not found; the wind fills r < r_CD'
      r_RS = r_CD
   end if

   !------------------------------------------------------------
   ! Isothermal layer (FS side)

   inquire(file=output_file('T_IL_FS.dat'), exist=ex_IL)
   IL = radiative .and. ex_IL
   if (radiative .and. .not. ex_IL) print *, 'WARNING (CD_abssyn): T_IL_FS.dat not found; no IL absorption'

   r_ILo  = r_CD
   nnk_IL = 0.d0
   ds_in  = 0.01d0*R0
   if (IL) then
      allocate( rhoIL(m_theta), T_IL(m_theta) )
      open(newunit=iu, file=output_file('H_FS.dat'))
      do l = 1, m_theta
         read(iu,*) dummy, Hw(l)
      end do
      close(iu)
      open(newunit=iu, file=output_file('T_IL_FS.dat'))
      do l = 1, m_theta
         read(iu,*) dummy, T_IL(l)
      end do
      close(iu)
      open(newunit=iu, file=output_file('termo_FS.dat'))
      do l = 1, m_theta
         read(iu,*) (dummy, j=1,11), rhoIL(l)        ! column 12
      end do
      close(iu)

      Hw     = frac_H*Hw                                ! IL width reduced by eta_H
      r_ILo = r_CD + Hw*g_rad
      ds_in = min(ds_in, 0.25d0*minval(Hw))            ! resolve the thinnest IL

      do l = 1, m_theta
         call gaunt_ff_calc(freq, size(freq), T_IL(l), Zq_ISM, gaunt)
         nnk_IL(:,l) = (rhoIL(l)/(mu_e*mp))*(rhoIL(l)/(mu_i*mp)) &
                     * 3.7d8*Zq_ISM**2/sqrt(T_IL(l)) * (1.d0-exp(-h*freq(1:n_abs)/(k*T_IL(l)))) &
                     / freq(1:n_abs)**3 * gaunt(1:n_abs)
      end do
      deallocate( rhoIL, T_IL )
   end if
   ds_in = max(ds_in, 1.d-4*R0)

   R_far  = 1.01d0*maxval(r_ILo)
   r_stop = max(5.d0*R0, 1.5d0*R_far)

   if (verbose) then
      write(*,'(a,l2,a,es10.2,a)') ' IL absorption:', IL, ';  step inside the shell = ', ds_in/R0, ' R0'
   end if
   if (verbose .and. shock == 'FS' .and. IL) print *, 'FS emitters placed at the outer edge of the IL'

   deallocate( g_rad, dlnr, Hw, gaunt )

   END subroutine shell_geometry

   !========================================================================
   !        tau_calc
   !========================================================================
   ! Optical depth of the wind (tau_w) and of the IL (tau_il) for a photon
   ! emitted at pos(:,l) towards the observer (e_obs, global.f90), for
   ! nu(1:n_abs). The ray geometry does not depend on nu: each ray is followed
   ! once (trace_ray) and the wind gives a single geometric integral I_w.
   ! RS emitters start at the CD; FS emitters at the outer edge of the IL
   ! (the FS emission comes from outside the IL).
   !========================================================================

   subroutine tau_calc(shock, pos, K_w, tau_w, tau_il)
   character(LEN=*), intent(in) :: shock
   real(dp), intent(in)  :: pos(:,:), K_w(:)
   real(dp), intent(out) :: tau_w(:,:), tau_il(:,:)
   integer  :: l, it0, n_occ
   real(dp) :: p0(3), I_w
   logical  :: occ

   tau_w  = 0.d0
   tau_il = 0.d0
   n_occ  = 0

   !$omp parallel do default(shared) schedule(dynamic) private(l, it0, p0, I_w, occ) reduction(+:n_occ)
   do l = 1, size(pos,2)
      it0 = merge(l, l-m_theta, l <= m_theta)
      p0  = pos(:,l)
      if (shock == 'FS' .and. IL) p0 = p0*(r_ILo(it0)/r_CD(it0))

      call trace_ray(p0, I_w, tau_il(l,1:n_abs), occ)

      if (occ) then
         tau_w(l,:)  = tau_occ
         tau_il(l,:) = 0.d0
         n_occ = n_occ + 1
      else
         tau_w(l,1:n_abs) = K_w(1:n_abs)*I_w
      end if
   end do
   !$omp end parallel do

   if (verbose .and. n_occ > 0) write(*,'(a,i0,a)') ' CD_abssyn: ', n_occ, ' emitters occulted by the star'

   END subroutine tau_calc

   !========================================================================
   !        trace_ray
   !========================================================================
   ! Follows the ray p0 + s*e_obs (s >= 0). The ray is divided into segments
   ! of constant region; each boundary is located by bisection between two
   ! consecutive steps (ds_in for r < R_far, 2% of r beyond). For each segment:
   !  - wind: I_w += int ds/r**4 = G_w(u_a) - G_w(u_b), u = s - s_c
   !  - IL:   tl  += int alpha_IL ds (midpoint rule, n_sub sub-steps)
   ! At r > r_stop (moving outwards) the march stops; if the photon is still
   ! in the wind (open tail) the wind integral is extended to infinity.
   !========================================================================

   subroutine trace_ray(p0, I_w, tl, occ)
   real(dp), intent(in)  :: p0(3)
   real(dp), intent(out) :: I_w          ! int ds/r**4 over the wind [cm^-3]
   real(dp), intent(out) :: tl(:)        ! tau of the IL, nu(1:n_abs)
   logical,  intent(out) :: occ          ! the ray crosses the star
   real(dp) :: s_c, b_imp, s, s_a, s_new, lo, hi, mid, q(3), rq
   integer  :: reg_a, reg_new, n, j

   I_w = 0.d0
   tl  = 0.d0

   ! closest approach to the star
   s_c = -dot_product(p0, e_obs)
   b_imp   = sqrt(max(dot_product(p0,p0) - s_c**2, 0.d0))
   occ = (s_c > 0.d0) .and. (b_imp < R_st)
   if (occ) return

   s     = 0.d0
   s_a   = 0.d0
   reg_a = region(p0)

   do n = 1, n_step_max
      q = p0 + s*e_obs
      rq = norm2(q)
      if (rq > r_stop .and. dot_product(q, e_obs) > 0.d0) exit

      s_new   = s + merge(ds_in, max(ds_in, 0.02d0*rq), rq < R_far)
      reg_new = region(p0 + s_new*e_obs)

      if (reg_new == reg_a) then
         s = s_new
      else
         ! boundary between s (reg_a) and s_new (not reg_a)
         lo = s
         hi = s_new
         do j = 1, n_bis
            mid = 0.5d0*(lo + hi)
            if (region(p0 + mid*e_obs) == reg_a) then
               lo = mid
            else
               hi = mid
            end if
         end do
         call add_segment(reg_a, s_a, hi)
         s_a   = hi
         s     = hi                       ! continue from the boundary (another one may follow)
         reg_a = region(p0 + hi*e_obs)
      end if
   end do

   ! last segment
   if (reg_a == REG_WIND) then
      I_w = I_w + G_w(s_a - s_c, b_imp)       ! up to infinity
   else
      call add_segment(reg_a, s_a, s)
   end if

   contains

      subroutine add_segment(reg, sa, sb)
      integer,  intent(in) :: reg
      real(dp), intent(in) :: sa, sb
      real(dp) :: ds_sub, w
      integer  :: jj, i1
      select case (reg)
      case (REG_WIND)
         I_w = I_w + G_w(sa - s_c, b_imp) - G_w(sb - s_c, b_imp)
      case (REG_IL)
         ds_sub = (sb - sa)/n_sub
         do jj = 1, n_sub
            call theta_interp(theta_of(p0 + (sa + (jj-0.5d0)*ds_sub)*e_obs), i1, w)
            tl = tl + ds_sub*( (1.d0-w)*nnk_IL(:,i1) + w*nnk_IL(:,i1+1) )
         end do
      end select
      end subroutine add_segment

   END subroutine trace_ray

   !========================================================================
   ! Region of the point q (see REG_* above)

   pure integer function region(q)
   real(dp), intent(in) :: q(3)
   real(dp) :: rq, w, th
   integer  :: i1

   th = theta_of(q)
   if (th > th_max) then
      region = REG_WIND                   ! open tail: the shell is not modelled beyond theta(m_theta)
      return
   end if

   rq = norm2(q)
   call theta_interp(th, i1, w)
   if (rq < (1.d0-w)*r_RS(i1) + w*r_RS(i1+1)) then
      region = REG_WIND
   else if (rq < (1.d0-w)*r_CD(i1) + w*r_CD(i1+1)) then
      region = REG_OUT                    ! shocked wind
   else if (IL .and. rq < (1.d0-w)*r_ILo(i1) + w*r_ILo(i1+1)) then
      region = REG_IL
   else
      region = REG_OUT
   end if

   end function region

   !========================================================================
   ! Polar angle of q, measured from the apex direction (-x)

   pure real(dp) function theta_of(q)
   real(dp), intent(in) :: q(3)
   theta_of = atan2(sqrt(q(2)**2 + q(3)**2), -q(1))
   end function theta_of

   !========================================================================
   ! Linear interpolation in theta: f(th) = (1-w)*f(i1) + w*f(i1+1)

   pure subroutine theta_interp(th, i1, w)
   real(dp), intent(in)  :: th
   integer,  intent(out) :: i1
   real(dp), intent(out) :: w
   real(dp) :: xi

   xi = (th - th_min)/dtheta
   i1 = min( max( int(xi) + 1, 1 ), m_theta-1 )
   w  = min( max( xi - real(i1-1,dp), 0.d0 ), 1.d0 )

   end subroutine theta_interp

   !========================================================================
   ! G_w(u,b) = int_u^infinity du'/(u'**2 + b**2)**2, i.e. int ds/r**4 from
   ! the point at distance u (along the ray) from the closest approach to the
   ! star, with impact parameter b:
   !   G_w(u) = [atan(b/u) - u*b/(u**2+b**2)]/(2*b**3)    for u > 0
   !   G_w(u) = pi/(2*b**3) - G_w(|u|)                    for u < 0
   ! For b << u the series (1/3 - 2t**2/5 + 3t**4/7)/u**3, t = b/u, avoids the
   ! cancellation (b -> 0: rays aligned with the star, e.g. the apex at i = 90 deg).

   pure real(dp) function G_w(ub, bb)
   real(dp), intent(in) :: ub, bb
   real(dp) :: ua, tt

   ua = abs(ub)
   if (bb <= 1.d-2*ua) then
      tt   = bb/ua
      G_w = (1.d0/3.d0 - 0.4d0*tt**2 + 3.d0/7.d0*tt**4)/ua**3
   else
      G_w = (atan2(bb, ua) - ua*bb/(ua**2 + bb**2))/(2.d0*bb**3)
   end if
   if (ub < 0.d0) G_w = pi/(2.d0*bb**3) - G_w

   end function G_w

   END module CD_abssyn