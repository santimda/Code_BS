module rad_H_alpha
   !============================================================================
   ! H-alpha recombination emission of the isothermal (ionized) layer of a
   ! radiative FS: map projected in the plane of the sky and its convolution
   ! with a Gaussian beam.
   !============================================================================
   use global
   use controls, only: verbose
   implicit none
   private
   public :: rad_H_alpha_run, project_H_alpha_map

   real(dp), parameter :: sigma_x = 2.d0, sigma_y = 2.d0   ! beam widths ["] (arbitrary)
   real(dp), parameter :: R_factor = 5.66d-18              ! 1 Rayleigh at H-alpha [erg/s/cm2/arcsec2]
   real(dp), parameter :: n_sig_margin = 3.d0              ! margin of the convolution grid around the data [sigma]

contains

   !============================================================================
   !        project_H_alpha_map
   !============================================================================
   ! H-alpha flux of each linear emitter at the rotation angle phi (CDrot.dat),
   ! for an optically thin IL:
   !   S = kappa_alpha * n_e*n_p * V_cell * convert**2   [erg/s/cm2]
   ! kappa_alpha = 2.85e-33*T**(-0.9) [erg/s/cm3/arcsec2 per n_e*n_p], with
   ! T = T_IL(theta) (T_IL_FS.dat), n_e = rho_IL/(mu_e*mp), n_p = X*rho_IL/mp,
   ! and V_cell = eta_H*dV_ring/n_emit: the n_emit = 4*(mphi-1) distinct
   ! azimuths of the map (phi = 0, ..., pi/2, + up/down and the mirror in
   ! main.f90) cover the ring once. S does not depend on phi; only the
   ! projected position does.
   ! Output (output/): S_alpha_proj_<shock>.dat with u ["], v ["], S.
   !============================================================================

   subroutine project_H_alpha_map(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_l = 2*m_theta    ! CD position, both sides
      real(dp), parameter :: convert = 206265/distance
      real(dp) :: e_u(3)
      integer  :: l, j, iu
      logical  :: ex
      real(dp), allocatable :: pos(:,:), rhoIL(:), T_IL(:), dV_ring(:), S_alpha(:)

      inquire(file=output_file('T_IL_'//shock//'.dat'), exist=ex)
      if (.not. ex) then
         print *, 'WARNING (rad_H_alpha): T_IL_'//shock//'.dat not found (non-radiative shock?); no H-alpha map'
         return
      end if

      allocate( pos(3,m_l), rhoIL(m_theta), T_IL(m_theta), dV_ring(m_theta), S_alpha(m_l) )

      !-----------------------------------------
      ! IL: density, temperature and volume of each ring (CD_termo.f90)

      open(newunit=iu, file=output_file('termo_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu,*) (dummy, j=1,11), rhoIL(l)        ! column 12
      end do
      close(iu)

      open(newunit=iu, file=output_file('T_IL_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu,*) dummy, T_IL(l)
      end do
      close(iu)

      open(newunit=iu, file=output_file('Hnorm_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu,*) dummy, dummy, dummy, dummy, dV_ring(l)   ! whole ring (2*pi)
      end do
      close(iu)

      !-----------------------------------------
      ! Flux of each linear emitter (same for the "up" and "down" sides)

      S_alpha(1:m_theta) = 2.85d-33*T_IL**(-0.9d0) * (rhoIL/(mu_e*mp))*(X_ab*rhoIL/mp) &
                         * frac_H*dV_ring/n_emit * convert**2
      S_alpha(m_theta+1:m_l) = S_alpha(1:m_theta)

      !-----------------------------------------
      ! Projection in the plane of the sky (as in CD_abssyn): u along
      ! e_u = (cos i, sin i, 0), v along z; both perpendicular to e_obs

      open(newunit=iu, file=output_file('CDrot.dat'))
      read(iu,*) pos
      close(iu)
      pos = R0*pos

      e_u = [ cos(incli), sin(incli), 0.d0 ]
      open(newunit=iu, file=output_file('S_alpha_proj_'//shock//'.dat'))
      write(iu,'(3es17.8e3)') (dot_product(e_u, pos(:,l))*convert, pos(3,l)*convert, S_alpha(l), l=1,m_l)
      close(iu)

      deallocate( pos, rhoIL, T_IL, dV_ring, S_alpha )

   END subroutine project_H_alpha_map

   !============================================================================
   !        rad_H_alpha_run
   !============================================================================
   ! Convolution of the H-alpha map (S_alpha_map_<shock>.dat, written in
   ! main.f90 for all phi) with a Gaussian beam of widths sigma_x, sigma_y.
   ! The map is read as (u, v, S); the convolution grid uses x = v (symmetric
   ! with respect to the axis) and y = u (along the projected axis).
   ! Surface brightness [R] = sum_l S_l/R_factor * exp(-d2) / (2*pi*sigma_x*sigma_y),
   ! d2 = (x_l-x_c)**2/(2*sigma_x**2) + (y_l-y_c)**2/(2*sigma_y**2), computed
   ! row by row in conv_rows. The grid covers the data plus n_sig_margin*sigma
   ! on each side (x symmetric with respect to the axis), spacing f_Ny*sigma.
   ! Output (output/): R0_proj_<system>.dat and S_alpha_convmap_<system>.dat
   ! (x ["], y ["], brightness [R]; one block per x).
   !============================================================================

   subroutine rad_H_alpha_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = n_emit*m_theta            ! number of points of the map
      real(dp) :: A_beam, D_x, D_y, w_x, w_y, half_x, y_min, y_max
      integer  :: i, j, l, m_x, m_y, iu
      real(dp), allocatable :: x_map(:), y_map(:), F_alpha(:), x_c(:), y_c(:), F_alpha_conv(:,:), EY(:,:)

      if (verbose) print *, 'Convolve the H-alpha emission map'

      A_beam = 2.d0*pi*sigma_x*sigma_y     ! Gaussian normalization factor

      open(newunit=iu, file=output_file('R0_proj_'//system//'.dat'))
      write(iu,'(2es17.8e3)') 0.d0, R0_proj/abs(cos(incli))
      close(iu)

      !------------------------------------------------------
      ! Map: u -> y_map, v -> x_map; flux from cgs to R*arcsec2

      allocate( x_map(m_tot), y_map(m_tot), F_alpha(m_tot) )
      open(newunit=iu, file=output_file('S_alpha_map_'//shock//'.dat'))
      read(iu,*) (y_map(l), x_map(l), F_alpha(l), l=1,m_tot)
      close(iu)
      F_alpha = F_alpha/R_factor

      !-------------------------------------------------------
      ! Convolution grid ["]: the data plus n_sig_margin*sigma on each side
      ! (x symmetric with respect to the axis), spacing f_Ny*sigma

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      half_x = maxval(abs(x_map)) + n_sig_margin*sigma_x
      m_x    = max( ceiling(2.d0*half_x/D_x), 1 )
      x_c    = ([(real(i,dp), i=1,m_x)] - 0.5d0*(m_x+1))*D_x

      y_min = minval(y_map) - n_sig_margin*sigma_y
      y_max = maxval(y_map) + n_sig_margin*sigma_y
      m_y   = max( ceiling((y_max - y_min)/D_y), 1 )
      y_c   = y_min + ([(real(j,dp), j=1,m_y)] - 0.5d0)*D_y

      !------------------------------------------
      ! Convolve (see conv_rows)

      w_x = 2.d0*sigma_x**2
      w_y = 2.d0*sigma_y**2
      allocate( EY(m_tot,m_y), F_alpha_conv(m_x,m_y) )
      EY = exp( -(spread(y_map,2,m_y) - spread(y_c,1,m_tot))**2/w_y )
      call conv_rows(x_map, F_alpha, x_c, w_x, EY, F_alpha_conv)
      F_alpha_conv = F_alpha_conv/A_beam
      deallocate( EY )

      !------------------------------------------
      ! Save

      open(newunit=iu, file=output_file('S_alpha_convmap_'//system//'.dat'))
      do i = 1, m_x
         write(iu,'(F11.4,X,F11.4,X,ES12.4E3)') (x_c(i), y_c(j), F_alpha_conv(i,j), j=1,m_y)
         write(iu,*)
      end do
      close(iu)

      deallocate( x_map, y_map, F_alpha, x_c, y_c, F_alpha_conv )

   END subroutine rad_H_alpha_run


   !============================================================================
   ! F(i,j) = sum_l S_l*exp(-(x_l-x_c(i))**2/w_x)*EY(l,j), EY(l,j) = exp(-dy2).
   ! Row by row (OpenMP), using only the points with dx2 < 40 (exp(-40) ~ 4e-18).

   subroutine conv_rows(x_map, S, x_c, w_x, EY, F)
      real(dp), intent(in)  :: x_map(:), S(:), x_c(:), w_x, EY(:,:)
      real(dp), intent(out) :: F(:,:)
      real(dp), allocatable :: dx2(:)
      integer,  allocatable :: idx(:)
      integer :: i, l

      !$omp parallel do default(shared) schedule(dynamic) private(i, dx2, idx)
      do i = 1, size(x_c)
         dx2 = (x_map - x_c(i))**2/w_x
         idx = pack([(l, l=1,size(x_map))], dx2 < 40.d0)
         if (size(idx) == 0) then
            F(i,:) = 0.d0
         else
            F(i,:) = matmul(S(idx)*exp(-dx2(idx)), EY(idx,:))
         end if
      end do
      !$omp end parallel do

   end subroutine conv_rows

END module rad_H_alpha