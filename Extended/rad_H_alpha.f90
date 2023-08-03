module rad_H_alpha
   ! This module calculates the H_alpha line emission map
   use global
   implicit none

   private
   public :: rad_H_alpha_run

   contains

   subroutine rad_H_alpha_run(shock)

      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      real(dp), parameter :: sigma_x = 20.d0, sigma_y = 20.d0  ! Arbitrary
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_alpha(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_alpha_conv
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam
      real(dp) :: S_map(m_tot)
      real(dp) :: w_x, w_y
      real(dp), parameter :: R_factor = 5.66d-18

      print*, '\n Calculate the H_alpha emission map'

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 

      ci = abs(cos(incli))

      open(188,file='R0_proj.dat')
         write(188,*) 0.d0,R0_proj/abs(cos(incli))
      close(188)

      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 5.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 5.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_alpha_conv(m_x,m_y)  )

      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      open(10,file='S_alpha_map_'//shock//'.dat')
      do l = 1, m_tot
         read(10,*) y_map(l),x_map(l),F_alpha(l)
         F_alpha(l) = F_alpha(l)/R_factor    ! Convert from cgs to R units
      enddo
      close(10)

      !------------------------------------------
      ! Convolve

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then  ! flux drops < 1e-3 farther from there, no need to sum    
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_alpha(l)*atenua
                  suma = suma + D_F
               endif
            end do

            F_alpha_conv(i,j) = suma/A_beam
         end do

      end do

      !------------------------------------------
      ! Data storage + consistency check:

      sum = 0.d0
      open (20,file='S_alpha_map.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(20,60) x_c(i), y_c(j), F_alpha_conv(i,j)
         enddo
         write(20,*) 
      enddo
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A12,X,Es10.2,X,A1)

      deallocate (x_c,y_c,F_alpha_conv)

   END subroutine rad_H_alpha_run

   !============================================================================
   ! rad_H_alpha_run (END)
   !============================================================================

END module rad_H_alpha