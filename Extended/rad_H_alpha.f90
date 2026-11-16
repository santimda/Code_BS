module rad_H_alpha
   !============================================================================
   ! H-alpha emission map of the isothermal layer of the FS
   !============================================================================
   use global
   use controls, only: verbose
   implicit none
   private
   public :: rad_H_alpha_run, project_H_alpha_map

   real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0   ! beam widths ["] (arbitrary)
   real(dp), parameter :: R_factor = 5.66d-18              ! 1 Rayleigh at H-alpha [erg/s/cm2/arcsec2]
   real(dp), parameter :: n_sig_margin = 3.d0              ! margin of the convolution grid around the data [sigma]

contains

   !============================================================================
   !        project_H_alpha_map
   !============================================================================
   ! H-alpha flux [erg/s/cm2] of each linear emitter at the current phi
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
      ! IL: density, temperature and volume of each ring

      open(newunit=iu, file=output_file('thermo_'//shock//'.dat'))
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
      ! Flux of each linear emitter

      S_alpha(1:m_theta) = 2.85d-33*T_IL**(-0.9d0) * (rhoIL/(mu_e*mp))*(X_ab*rhoIL/mp) &
                         * frac_H*dV_ring/n_emit * convert**2
      S_alpha(m_theta+1:m_l) = S_alpha(1:m_theta)

      !-----------------------------------------
      ! Projection in the plane of the sky

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
   ! Convolution of the H-alpha map with a Gaussian beam; brightness in R
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
      ! Map (u -> y_map, v -> x_map), in R*arcsec2

      allocate( x_map(m_tot), y_map(m_tot), F_alpha(m_tot) )
      open(newunit=iu, file=output_file('S_alpha_map_'//shock//'.dat'))
      read(iu,*) (y_map(l), x_map(l), F_alpha(l), l=1,m_tot)
      close(iu)
      F_alpha = F_alpha/R_factor

      !-------------------------------------------------------
      ! Convolution grid ["]

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
      ! Convolve

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
   ! Gaussian convolution, row by row

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