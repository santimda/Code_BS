module CD_abssyn
   use global
   use controls, only: verbose
   use gaunt_factor
   use map_freqs, only: maps, nu_index
   implicit none
   private
   public :: CD_abssyn_run

   contains

   !========================================================================
   !        CD_abssyn_run
   !========================================================================
   ! Free-free absorption of the synchrotron emission of the shock(s) for the
   ! linear emitters at the current rotation angle phi (CDrot.dat).
   ! Absorbers:
   !  - unshocked stellar wind inside the BS (tau_w), using alpha(r,nu) from
   !    windalphaff.dat (opac_ff.f90)
   !  - ionized layer of the FS (tau_il), only for shock = 'BOTH' and a
   !    radiative FS (termo_FS.dat, H_FS.dat)
   !  - Stromgren sphere (tau_str), radius from Stromgren.dat
   ! Wind and IL absorption are computed for nu < nu_max_ff (global.f90).
   ! The f-f absorption in the shocks is neglected given the large T.
   !
   ! Output (arguments): E_sy, photon energy [eV], and the SED summed over the
   ! linear emitters, unabsorbed (Lsy_sum) and absorbed (Lsy_sum_abs), as
   ! nu*L_nu [erg/s]. Files: radio maps, average map, CD_tauff.dat and
   ! radsyn_abs*.dat (diagnostic).
   !========================================================================

   subroutine CD_abssyn_run(shock, radiative, E_sy, Lsy_sum, Lsy_sum_abs)
   !$ use omp_lib
   character(LEN=*), intent(in) :: shock
   logical, intent(in) :: radiative
   real(dp), intent(out) :: E_sy(m_nu), Lsy_sum(m_nu), Lsy_sum_abs(m_nu)
   integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
   real(dp) :: log_nu_int, flux, convert, tginc, sinc, cinc, R_str, h_avg, L_FS
   integer :: l, i, nl, n_abs, N_avg, im, n_map, iu
   character(len=:), allocatable :: tag
   integer,  allocatable :: i_nu(:)                       ! index of each map frequency in nu(:)
   real(dp), allocatable :: L_em(:,:)                     ! nu*L_nu of each linear emitter (position, nu)
   real(dp), allocatable :: tau_w(:,:), tau_il(:,:), tau_str(:), tau_tot(:,:)
   real(dp), allocatable :: gaunt(:), kappa_ff(:)
   real(dp), allocatable :: u_map(:), v_map(:), emis(:,:), emis_abs(:,:), w_avg(:)

   allocate( x(m_l), y(m_l), z(m_l), nu(m_nu) )
   allocate( L_em(m_l,m_nu), gaunt(m_nu), kappa_ff(m_nu) )
   allocate( tau_w(m_l,m_nu), tau_il(m_l,m_nu), tau_str(m_nu), tau_tot(m_l,m_nu) )

   if (verbose) print *, 'Calculate the f-f absorption of the synchrotron emission'

   convert = 206265/distance             ! image to angular size ["]
   flux = 1.d26/(4.d0*pi*(distance)**2)  ! luminosity to flux in mJy
   tginc = tan(incli)
   sinc = sin(incli)
   cinc = cos(incli)
   ! set omp threads number
   !$ call omp_set_num_threads(nth)

   !======================================================
   ! Geometry and frequency grid

   open(11,file=output_file('CDrot.dat'))
   do l=1,m_l     ! "up" and "down"
      read(11,*) x(l),y(l),z(l)
   end do
   close(11)
   x = R0*x
   y = R0*y
   z = R0*z

   open(12,file=output_file('nu.dat'))
      read(12,*) nu
   close(12)
   log_nu_int = log(nu(2)/nu(1))

   !======================================================
   ! f-f absorption coefficient per unit n_e*n_i at T_ISM: alpha = kappa_ff*n**2

   call gaunt_ff_calc(nu, m_nu, T_ISM, Zq_ISM, gaunt)
   kappa_ff = 3.7d8 * Zq_ISM**2/sqrt(T_ISM) * (1.d0-exp(-h*nu/(k*T_ISM)))/nu**3 * gaunt

   !======================================================
   ! Opacities: stellar wind and ionized layer (nu < nu_max_ff), Stromgren sphere

   n_abs = count(nu < nu_max_ff)
   call tau_calc(shock, radiative, n_abs, kappa_ff, tau_w, tau_il)

   open(30,file=output_file('Stromgren.dat'))
      read(30,*) R_str
   close(30)
   tau_str = kappa_ff * n_ISM**2 * R_str

   tau_tot = tau_w + tau_il + spread(tau_str, 1, m_l)

   open(25,file='CD_tauff.dat')
   write(25,'(i6,2es17.8e3)') ((l, nu(i), tau_tot(l,i), i=1,m_nu), l=1,m_l)
   close(25)

   !/////////////////////////////////////////////////////////////////////
   ! Synchrotron emission of each linear emitter: emisyn_*.dat has
   ! E [eV] and nu*L_nu [erg/s] (already multiplied by E)

   if (shock == 'RS') then
      open (32,file=output_file('emisyn_RS.dat'))
      do l=1,m_theta
         do i=1,m_nu
            read(32,*) nl,E_sy(i),L_em(l,i)
         end do
      end do
      close(32)

   else if (shock == 'FS') then
      open (33,file=output_file('emisyn_FS.dat'))
      do l=1,m_theta
         do i=1,m_nu
            read(33,*) nl,E_sy(i),L_em(l,i)
         end do
      end do
      close(33)

   else if (shock == 'BOTH') then
      open (32,file=output_file('emisyn_RS.dat'))
      open (33,file=output_file('emisyn_FS.dat'))
      do l=1,m_theta
         do i=1,m_nu
            read(32,*) nl,E_sy(i),L_em(l,i)
            read(33,*) nl,dummy,L_FS
            !L_em(l,i) = L_em(l,i) + L_FS          ! FS synchrotron emission (disabled)
         end do
      end do
      close(32)
      close(33)
   end if
   L_em(m_theta+1:m_l,:) = L_em(1:m_theta,:)       ! symmetry ("down" side)

   !-----------------------------------------------------------------------
   ! Total unabsorbed and absorbed emission (sum over linear emitters)

   Lsy_sum     = sum(L_em, dim=1)
   Lsy_sum_abs = sum(L_em*exp(-tau_tot), dim=1)

   if (shock == 'BOTH') then
      tag = ''
   else
      tag = shock
   end if

   if (shock == 'BOTH') then
      open(34,file='radsyn_abs.dat')
   else
      open(34,file='radsyn_abs_'//shock//'.dat')
   end if
   write(34,'(3es17.8e3)') (E_sy(i), Lsy_sum(i), Lsy_sum_abs(i), i=1,m_nu)
   close(34)

   !-----------------------------------------------------------------------
   ! Emission maps at the frequencies selected in controls (maps(:)),
   ! projected in the plane of the sky (u_map,v_map)

   n_map = size(maps)
   allocate( i_nu(n_map), u_map(m_l), v_map(m_l), emis(m_l,max(n_map,1)), emis_abs(m_l,max(n_map,1)) )

   i_nu     = nu_index(maps%nu, nu(1), log_nu_int, m_nu)
   u_map    = (y + x/abs(tginc))*abs(sinc)*convert
   v_map    = z*convert
   emis(:,1:n_map)     = L_em(:,i_nu) * spread(flux/nu(i_nu), 1, m_l)
   emis_abs(:,1:n_map) = emis(:,1:n_map) * exp(-tau_tot(:,i_nu))

   do im = 1, n_map
      open(newunit=iu, file='CD_sync_radiomap_'//tag//trim(maps(im)%label)//'.dat')
      write(iu,'(4es17.8e3)') (u_map(l), v_map(l), emis(l,im), emis_abs(l,im), l=1,m_l)
      close(iu)
   end do

   !------------------------------------------------------------
   ! Average flux between nu_av_min and nu_av_max. In the log grid
   ! dnu = nu*(h_avg-1), so (L_em/nu)*flux*dnu = L_em*flux*(h_avg-1)

   N_avg = count((nu >= nu_av_min) .and. (nu <= nu_av_max))
   h_avg = (nu_av_max/nu_av_min)**(1.d0/(N_avg-1))
   allocate( w_avg(m_nu) )
   w_avg = merge(flux*(h_avg-1.d0)/(nu_av_max-nu_av_min), 0.d0, (nu >= nu_av_min) .and. (nu <= nu_av_max))

   emis(:,1)     = matmul(L_em, w_avg)                  ! reuse emis(:,1:2) for the average map
   emis_abs(:,1) = matmul(L_em*exp(-tau_tot), w_avg)

   open(76,file='CD_sync_radiomap_avg'//tag//'.dat')
   write(76,'(4es17.8e3)') (u_map(l), v_map(l), emis(l,1), emis_abs(l,1), l=1,m_l)
   close(76)

   deallocate( i_nu, u_map, v_map, emis, emis_abs, w_avg )
   deallocate( x, y, z, nu, L_em, gaunt, kappa_ff )
   deallocate( tau_w, tau_il, tau_str, tau_tot )

   END subroutine CD_abssyn_run

   !========================================================================

   !========================================================================
   !        tau_calc
   !========================================================================
   ! Optical depth of the stellar wind (tau_w) and of the ionized layer of
   ! the FS (tau_il) for a photon emitted at (x(l),y(l),z(l)) towards the
   ! observer, for nu(1:n_abs). The ray geometry does not depend on nu, so
   ! each ray is followed once and all frequencies are accumulated together.
   !  - Wind: integrated along the LOS while the photon is inside the BS,
   !    starting dr0 (estimate of the shock width) away from the emitter.
   !  - IL: crossed at most once, where the photon leaves the BS through the
   !    shell; alpha_IL is integrated for r_norm*R0 <= r < r_norm*R0 + H,
   !    from the last point inside the BS and with a finer step.
   !========================================================================

   subroutine tau_calc(shock, radiative, n_abs, kappa_ff, tau_w, tau_il)
   character(LEN=*), intent(in) :: shock
   logical,  intent(in)  :: radiative
   integer,  intent(in)  :: n_abs                ! absorption computed for nu(1:n_abs)
   real(dp), intent(in)  :: kappa_ff(:)          ! alpha_ff/n**2 at T_ISM, Zq_ISM
   real(dp), intent(out) :: tau_w(:,:), tau_il(:,:)
   integer, parameter :: m_l = 2*m_theta
   real(dp) :: r_ff(mr_ff), dr_tab, dtheta, k_los(3), p(3), p_in(3)
   real(dp) :: dr0, ds, ds_IL, r1, rho_l
   integer  :: l, it, ir
   logical  :: IL, crosses_IL
   real(dp), allocatable :: alphaffw(:,:), n_IL(:), H_IL(:)

   allocate( alphaffw(mr_ff,m_nu), n_IL(m_theta), H_IL(m_theta) )

   tau_w  = 0.d0
   tau_il = 0.d0
   IL = (shock == 'BOTH') .and. radiative

   !------------------------------------------------------------
   ! Wind opacity table alpha(r_ff,nu) (opac_ff.f90)

   open(15,file=output_file('r_ff.dat'))
   read(15,*) r_ff
   close(15)
   dr_tab = r_ff(2)-r_ff(1)          ! equiespaciated table

   open(18,file=output_file('windalphaff.dat'))
   read(18,*) ((alphaffw(ir,it), it=1,m_nu), ir=1,mr_ff)
   close(18)

   !------------------------------------------------------------
   ! Ionized layer: density and width (CD_termo.f90)

   if (IL) then
      open(20,file=output_file('termo_FS.dat'))
      open(21,file=output_file('H_FS.dat'))
      do l=1, m_theta
         read(20,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,rho_l
         read(21,*) dummy, H_IL(l), dummy, dummy
         n_IL(l) = rho_l/mu_SISM/mp
      end do
      close(20)
      close(21)
   end if

   !------------------------------------------------------------

   dtheta = theta(2)-theta(1)                 ! equiespaciated
   k_los  = [-sin(incli), cos(incli), 0.d0]   ! LOS direction; i = 0 is face-on

   !$omp parallel do default(shared) schedule(dynamic) &
   !$omp   private(l, dr0, ds, ds_IL, p, p_in, r1, it, ir, crosses_IL)
   do l = 1, m_l        ! Position in the CD

      dr0 = 0.2d0*r_norm(merge(l, l-m_theta, l <= m_theta))*R0   ! Estimate of the shock width, USE H0 for consistency!!!!
      ds  = dr0/16.d0                                            ! integration step along the LOS

      ! Stellar wind: integrate along the LOS while the photon is inside the BS
      p    = [x(l), y(l), z(l)] + dr0*k_los
      p_in = [x(l), y(l), z(l)]       ! last point inside the BS (starts at the emitter)
      crosses_IL = .false.
      do
         r1 = norm2(p)
         if (r1 > r_ff(mr_ff)) exit                                ! too far (escapes through the open tail)
         it = theta_index(p)
         if (r1 < r_norm(it)*R0 .or. it == m_theta) then           ! "inside" the BS
            ir = min( max( int((r1-r_ff(1))/dr_tab) + 1, 1 ), mr_ff )
            tau_w(l,1:n_abs) = tau_w(l,1:n_abs) + alphaffw(ir,1:n_abs)*ds
            p_in = p
         else                                                      ! leaves the BS through the shell at theta(it)
            crosses_IL = .true.
            exit
         end if
         p = p + ds*k_los
      end do

      ! Ionized layer (FS side)
      if (IL .and. crosses_IL) then
         if (H_IL(it) > 0.d0) then
            ds_IL = min(ds, H_IL(it)/10.d0)
            p = p_in
            do
               r1 = norm2(p)
               if (r1 > r_ff(mr_ff)) exit
               it = theta_index(p)
               if (r1 >= r_norm(it)*R0 + H_IL(it) .and. it < m_theta) exit   ! beyond the IL
               if (r1 >= r_norm(it)*R0 .and. it < m_theta) then              ! inside the IL
                  tau_il(l,1:n_abs) = tau_il(l,1:n_abs) + kappa_ff(1:n_abs)*n_IL(it)**2*ds_IL
               end if
               p = p + ds_IL*k_los
            end do
         end if
      end if

   end do
   !$omp end parallel do

   deallocate( alphaffw, n_IL, H_IL )

   contains

   ! Index of the polar angle (measured from the apex direction, -x) of q
   pure integer function theta_index(q)
      real(dp), intent(in) :: q(3)
      theta_index = min( max( int((atan2(sqrt(q(2)**2 + q(3)**2), -q(1)) - theta(1))/dtheta) + 1, 1 ), m_theta )
   end function theta_index

   END subroutine tau_calc

   END module CD_abssyn