module SW_conv
   ! This program makes the convolution of the emission in the plane of the sky 
   ! with the instrument response 
   use global
   implicit none

   private
   public :: SW_conv1G_run,SW_conv110M_run,SW_conv190M_run, SW_conv3G_run, SW_conv5G_run, SW_conv6G_run, SW_conv10G_run, SW_convavg_run,SW_spex_110_190_run, SW_spex_run, SW_spex_6_10_run
   public :: SW_conv150M_run,SW_conv325M_run

   contains

   !============================================================================
   ! SW_conv110M_run
   !============================================================================

   subroutine SW_conv110M_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 16.d0, sigma_y = 4.d0  ! LOFAR https://science.astron.nl/telescopes/lofar/lofar-system-overview/observing-modes/lofar-imaging-capabilities-and-sensitivity/ 
      !real(dp), parameter :: sigma_x = 77.d0, sigma_y = 77.d0
      real(dp), parameter :: sigma_x = 4.d0, sigma_y = 4.d0
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_110Mabs(m_tot),F_110Mun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map110Mabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) ::  suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      open(188,file='R0_proj.dat')
         write(188,*) 0.d0,R0_proj/abs(cos(incli))
      close(188)

      !///////////////////////////////////////////////////////
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
      allocate( x_c(m_x), y_c(m_y), F_map110Mabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_110MHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_110MHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_110Mun(l),F_110Mabs(l)
         sum = sum + F_110Mabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(110MHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_110M.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'110M.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(110M)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_110MHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(110MHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_110Mabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map110Mabs(i,j) = suma + suma2 + suma3
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_110Mabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map110Mabs(i,j)
            sum = sum + F_map110Mabs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(110MHz)=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map110Mabs)

   END subroutine SW_conv110M_run

   !============================================================================
   ! SW_conv110M_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv150M_run
   !============================================================================

   subroutine SW_conv150M_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      real(dp), parameter :: theta_x = 25.d0, theta_y = 25.d0  ! The FWHM
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))  ! Conversion factor between FWHM and sigma
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc  ! The sigma
      !real(dp), parameter :: sigma_x = 25.d0, sigma_y = 25.d0
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_150Mabs(m_tot),F_150Mun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map150Mabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) ::  suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      open(188,file='R0_proj.dat')
         write(188,*) 0.d0,R0_proj/abs(cos(incli))
      close(188)

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map150Mabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_150MHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_150MHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_150Mun(l),F_150Mabs(l)
         sum = sum + F_150Mabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(150MHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_150M.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'150M.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(150M)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_150MHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(150MHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_150Mabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map150Mabs(i,j) = suma + suma2 + suma3
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_150Mabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map150Mabs(i,j)
            sum = sum + F_map150Mabs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(150MHz)=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map150Mabs)

   END subroutine SW_conv150M_run

   !============================================================================
   ! SW_conv150M_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv190M_run
   !============================================================================

   subroutine SW_conv190M_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 16.d0, sigma_y = 4.d0  ! LOFAR https://science.astron.nl/telescopes/lofar/lofar-system-overview/observing-modes/lofar-imaging-capabilities-and-sensitivity/ 
      real(dp), parameter :: sigma_x = 4.d0, sigma_y = 4.d0
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_190Mabs(m_tot),F_190Mun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map190Mabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) ::  suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map190Mabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_190MHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_190MHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_190Mun(l),F_190Mabs(l)
         sum = sum + F_190Mabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(190MHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_190M.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'190M.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(190M)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_190MHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(190MHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_190Mabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map190Mabs(i,j) = suma + suma2 + suma3
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_190Mabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map190Mabs(i,j)
            sum = sum + F_map190Mabs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(190MHz)=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map190Mabs)

   END subroutine SW_conv190M_run

   !============================================================================
   ! SW_conv190M_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv325M_run
   !============================================================================

   subroutine SW_conv325M_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 77.d0, sigma_y = 77.d0 ! Brookes
      real(dp), parameter :: theta_x = 54., theta_y = 54.d0/sin(dec)  ! The FWHM (full axis, not semi-axis)
      !real(dp), parameter :: theta_x = 77., theta_y = 54.  ! The FWHM (full axis, not semi-axis)
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))  ! Conversion factor between FWHM and sigma
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc  ! The sigma
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_325Mabs(m_tot),F_325Mun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map325Mabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) ::  suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      open(188,file='R0_proj.dat')
         write(188,*) 0.d0,R0_proj/abs(cos(incli))
      close(188)

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 2.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 2.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map325Mabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_325MHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_325MHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_325Mun(l),F_325Mabs(l)
         sum = sum + F_325Mabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(325MHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_325M.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'325M.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(325M)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_325MHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(325MHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_325Mabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map325Mabs(i,j) = suma + suma2 + suma3
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_325Mabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map325Mabs(i,j)
            sum = sum + F_map325Mabs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(325MHz)=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map325Mabs)

   END subroutine SW_conv325M_run

   !============================================================================
   ! SW_conv325M_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv1G_run
   !============================================================================

   subroutine SW_conv1G_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0  ! Benaglia2010
      real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 18.3d0, sigma_y = 9.8d0    ! Moutzouri2022
      !real(dp), parameter :: sigma_x = 9.2d0, sigma_y = 7.8d0    ! Moutzouri2022
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0
      !    real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !    real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)  
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_1Gabs(m_tot),F_1Gun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map1Gabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) ::  suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map1Gabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_1GHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_1GHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_1Gun(l),F_1Gabs(l)
         sum = sum + F_1Gabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(1.4GHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_1G.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'1G.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(1.4G)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_1GHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(1.4GHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0  !1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_1Gabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map1Gabs(i,j) = suma + suma2 + suma3
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_1Gabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map1Gabs(i,j)
            sum = sum + F_map1Gabs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(1.4GHz)=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map1Gabs)

   END subroutine SW_conv1G_run

   !============================================================================
   ! SW_conv1G_run (END)
   !============================================================================


   !============================================================================
   ! SW_conv3G_run
   !============================================================================

   subroutine SW_conv3G_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)    ! Total number of data
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))         ! theta to sigma
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: theta_x = 20.2d0, theta_y = 12.5d0   ! HPBW = Resolution 
      !real(dp), parameter :: theta_x = 20.d0, theta_y = 10.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc   ! Benaglia2010
      !real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0   ! Benaglia2020
      !real(dp), parameter :: sigma_x = 20.2d0/coc, sigma_y = 12.5d0/coc  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 13.d0, sigma_y = 10.d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 5.d0
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      !real(dp), parameter :: sigma_x = 18.3d0, sigma_y = 9.8d0    ! Moutzouri2022
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_3Gabs(m_tot),F_3Gun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map3Gabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) :: S_map3(m_tot),suma3
      real(dp) :: w_x, w_y

      A_beam = 2.d0*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      open(188,file='R0_proj.dat')
         write(188,*) 0.d0,R0_proj/abs(cos(incli))
      close(188)

      print*, '\n 3 GHz'

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      !ancho_x = 4.d0*R0_proj/ci * 1.3d1 ! [mas](extra margin)
      !ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d1

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y
      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map3Gabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_3GHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_3GHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_3Gun(l),F_3Gabs(l)
         sum = sum + F_3Gabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(3.0GHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_3G.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'3G.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(3.0G)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_3GHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(3.0GHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y ! y-coord of the center of the (x_i,y_j) ellipse
            !if ((i==m_x/2 .and. j==m_y/2) .or. (i==m_x/2+1 .and. j==m_y/2) &
            !         .or. (i==m_x/2 .and. j==m_y/2+1) .or. (i==m_x/2+1 .and. j==m_y/2+1)) then
            suma = 0.d0!1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then 
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_3Gabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
               D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_F
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_F
            enddo

            F_map3Gabs(i,j) = suma + suma2 + suma3! >>>>> CHECK! FF FLUX TOO HIGH!
         enddo

      enddo


      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (18,file='SW_conv_3Gabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(18,61) x_c(i),y_c(j),F_map3Gabs(i,j)
            sum = sum + F_map3Gabs(i,j)*mJy_beam
         enddo
         write(18,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(3.0GHz)=',sum,'mJy'
      write(*,*) ''
      close(18)

      61  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map3Gabs)

   END subroutine SW_conv3G_run

   !============================================================================
   ! SW_conv3G_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv5G_run
   !============================================================================

   subroutine SW_conv5G_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      !    real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !    real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)  
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_5Gabs(m_tot),F_5Gun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map5Gabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) :: suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y
      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map5Gabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_5GHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_5GHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_5Gun(l),F_5Gabs(l)
         sum = sum + F_5Gabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(4.8GHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_5G.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'5G.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(4.8G)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_5GHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(4.8GHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0!1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then 
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_5Gabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map5Gabs(i,j) = suma !+ suma2 + suma3
         enddo

      enddo


      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open(18,file='SW_conv_5Gabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(18,61) x_c(i),y_c(j),F_map5Gabs(i,j)
            sum = sum + F_map5Gabs(i,j)*mJy_beam
         enddo
         write(18,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(4.8GHz)=',sum,'mJy'
      write(*,*) ''
      close(18)

      61  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map5Gabs)

   END subroutine SW_conv5G_run

   !============================================================================
   ! SW_conv5G_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv6G_run
   !============================================================================

   subroutine SW_conv6G_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      !real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      real(dp), parameter :: theta_x = 16.d0, theta_y = 10.d0     ! BD+ 60?
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)  
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_6Gabs(m_tot),F_6Gun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map6Gabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) :: suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y
      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map6Gabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_6GHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_6GHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_6Gun(l),F_6Gabs(l)
         sum = sum + F_6Gabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(6.0GHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_6G.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'6G.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(6.0G)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_6GHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(6.0GHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0!1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then 
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_6Gabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map6Gabs(i,j) = suma + suma2 + suma3
         enddo

      enddo


      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open(18,file='SW_conv_6Gabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(18,61) x_c(i),y_c(j),F_map6Gabs(i,j)
            sum = sum + F_map6Gabs(i,j)*mJy_beam
         enddo
         write(18,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(6.0GHz)=',sum,'mJy'
      write(*,*) ''
      close(18)

      61  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map6Gabs)

   END subroutine SW_conv6G_run

   !============================================================================
   ! SW_conv6G_run (END)
   !============================================================================

   !============================================================================
   ! SW_conv10G_run
   !============================================================================

   subroutine SW_conv10G_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      !real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      real(dp), parameter :: theta_x = 16.d0, theta_y = 10.d0     ! BD+ 60?
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)  
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_10Gabs(m_tot),F_10Gun(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map10Gabs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) :: suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y
      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map10Gabs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_10GHz.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_10GHz.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_10Gun(l),F_10Gabs(l)
         sum = sum + F_10Gabs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S(10.GHz)=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_10G.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'10G.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff(10.G)=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_10GHz.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S(10.GHz)=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0!1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then 
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_10Gabs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map10Gabs(i,j) = suma + suma2 + suma3
         enddo

      enddo


      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open(18,file='SW_conv_10Gabs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(18,61) x_c(i),y_c(j),F_map10Gabs(i,j)
            sum = sum + F_map10Gabs(i,j)*mJy_beam
         enddo
         write(18,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(10.GHz)=',sum,'mJy'
      write(*,*) ''
      close(18)

      61  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map10Gabs)

   END subroutine SW_conv10G_run

   !============================================================================
   ! SW_conv10G_run (END)
   !============================================================================


   !============================================================================
   ! SW_convavg_run
   !============================================================================

   subroutine SW_convavg_run(shock)
      character(LEN=*), intent(in) :: shock
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      !real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 18.3d0, sigma_y = 9.8d0    ! Moutzouri2022
      !real(dp), parameter :: sigma_x = 9.2d0, sigma_y = 7.8d0    ! Moutzouri2022
      !real(dp), parameter :: theta_x = 16.d0, theta_y = 10.d0     ! BD+ 60?
      real(dp), parameter :: theta_x = 19.7d0, theta_y = 8.25d0    ! BD+43?
      !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc ! HPBW = Resolution (Ben2010)  
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_avg_abs(m_tot),F_avg_un(m_tot)
      real(dp) :: ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map_avg_abs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, suma2,x_map2(m_map),y_map2(m_map),S_map2(m_tot)
      real(dp) :: suma3,S_map3(m_tot)
      real(dp) :: w_x, w_y

      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))

      print*, '\n Average'

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 2.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      !ancho_x = 4.d0*R0_proj/ci * 1.3d1 ! [mas](extra margin)
      !ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d1

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y
      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map_avg_abs(m_x,m_y)  )

      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_avg.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_avg.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_avg_un(l),F_avg_abs(l)
         sum = sum + F_avg_abs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S_sync_av=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_avg.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'avg.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_f--f_av=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_avg.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S_hot_lay=',sum,'mJy'
      end if

      !------------------------------------------

      w_x = 2.d0*sigma_x**2! Weight in x distance
      w_y = 2.d0*sigma_y**2! Weight in y distance

      do i = 1, m_x
         x_c(i) = x_c0 + (i-0.5)*D_x ! x-coord of the center of the i-th ellipse

         do j = 1, m_y
            y_c(j) = y_c0 + (j-0.5)*D_y ! y-coord of the center of the (x_i,y_j) ellipse

            suma = 0.d0!1.d-9  ! simulate noise (?)
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then 
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = F_avg_abs(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_f
            enddo

            suma2 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map2(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma2 = suma2 + D_f
            enddo

            suma3 = 0.d0
            do l = 1, m_tot
               d_pesada2 = ((x_map(l)-x_c(i))**2)/w_x + ((y_map(l)-y_c(j))**2)/w_y
               if (d_pesada2 < 6.d0) then      
                  ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
                  atenua = exp( -d_pesada2 ) !/ norm
                  D_F = S_map3(l)*atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma3 = suma3 + D_f
            enddo

            F_map_avg_abs(i,j) = suma + suma2 + suma3
         enddo

      enddo


      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open(18,file='SW_conv_avg_abs.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(18,61) x_c(i),y_c(j),F_map_avg_abs(i,j)
            sum = sum + F_map_avg_abs(i,j)*mJy_beam
         enddo
         write(18,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S(4-12av)=',sum,'mJy'
      write(*,*) ''
      close(18)

      61  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map_avg_abs)

   END subroutine SW_convavg_run

   !============================================================================
   ! SW_convavg_run (END)
   !============================================================================

   !============================================================================
   ! SW_spex_110_190_run
   !============================================================================

   subroutine SW_spex_110_190_run  
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      real(dp), parameter :: sigma_x = 4.d0, sigma_y = 4.d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      !    real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !    real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)
      real(dp) :: ancho_x, ancho_y, D_x, D_y, lognu2nu1, ci
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map110Mabs,F_map190Mabs,spex
      integer :: m_x,m_y

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 10.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 10.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction

      allocate( x_c(m_x), y_c(m_y), F_map110Mabs(m_x,m_y), F_map190Mabs(m_x,m_y), spex(m_x,m_y)  )


      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      open (17,file='SW_conv_110Mabs.dat')
      open (18,file='SW_conv_190Mabs.dat')
      open (19,file='SW_spex_110_190.dat')

      lognu2nu1 = log(190./110.)

      do i = 1, m_x
         do j = 1, m_y
            read(17,61) x_c(i),y_c(j),F_map110Mabs(i,j)
            read(18,61) x_c(i),y_c(j),F_map190Mabs(i,j)
            if ( F_map110Mabs(i,j) > 1.d-3 .and. F_map190Mabs(i,j) > 2.d-2 ) then  ! at least 20 uJy
               spex(i,j) = log(F_map190Mabs(i,j)/F_map110Mabs(i,j)) / lognu2nu1
               !             write(*,*) F_map1Gabs(i,j),F_map5Gabs(i,j),spex(i,j)
               write(19,61) x_c(i),y_c(j),spex(i,j)
            endif
         enddo
         read(17,*) 
         read(18,*) 
         write(19,*) 
      enddo
      close(17)
      close(18)
      close(19)

      61  format(F11.4,X,F11.4,X,E11.4)

      deallocate (x_c,y_c,F_map110Mabs,F_map190Mabs,spex)

   END subroutine SW_spex_110_190_run

   !============================================================================
   ! SW_spex_run (END)
   !============================================================================


   !============================================================================
   ! SW_spex_run
   !============================================================================

   subroutine SW_spex_run  
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      !    real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !    real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)
      real(dp) :: ancho_x, ancho_y, D_x, D_y, lognu2nu1, ci
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map1Gabs,F_map5Gabs,spex
      integer :: m_x,m_y

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction

      allocate( x_c(m_x), y_c(m_y), F_map1Gabs(m_x,m_y), F_map5Gabs(m_x,m_y), spex(m_x,m_y)  )


      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      open (17,file='SW_conv_1Gabs.dat')
      open (18,file='SW_conv_5Gabs.dat')
      open (19,file='SW_spex.dat')

      lognu2nu1 = log(4.86/1.42)

      do i = 1, m_x
         do j = 1, m_y
            read(17,61) x_c(i),y_c(j),F_map1Gabs(i,j)
            read(18,61) x_c(i),y_c(j),F_map5Gabs(i,j)
            if ( F_map1Gabs(i,j) > 1.d-3 .and. F_map5Gabs(i,j) > 2.d-2 ) then  ! at least 20 uJy
               spex(i,j) = log(F_map5Gabs(i,j)/F_map1Gabs(i,j)) / lognu2nu1
               !             write(*,*) F_map1Gabs(i,j),F_map5Gabs(i,j),spex(i,j)
               write(19,61) x_c(i),y_c(j),spex(i,j)
            endif
         enddo
         read(17,*) 
         read(18,*) 
         write(19,*) 
      enddo
      close(17)
      close(18)
      close(19)

      61  format(F11.4,X,F11.4,X,E11.4)

      deallocate (x_c,y_c,F_map1Gabs,F_map5Gabs,spex)

   END subroutine SW_spex_run

   !============================================================================
   ! SW_spex_run (END)
   !============================================================================

   !============================================================================
   ! SW_spex_6_10_run
   !============================================================================

   subroutine SW_spex_6_10_run  
      !real(dp), parameter :: sigma_x = 12.d0, sigma_y = 12.d0   ! Benaglia2010
      !real(dp), parameter :: sigma_x = 20.2d0, sigma_y = 12.5d0  ! Benaglia2021
      !real(dp), parameter :: sigma_x = 10.d0, sigma_y = 10.d0  ! Moutzouri2022
      real(dp), parameter :: theta_x = 16.d0, theta_y = 10.d0     ! BD+ 60?
      real(dp), parameter :: coc = 1.d0
      real(dp), parameter :: sigma_x = theta_x/coc, sigma_y = theta_y/coc
      !    real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))        ! Benaglia gives sigma, not theta!
      !    real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution (Ben2010)
      real(dp) :: ancho_x, ancho_y, D_x, D_y, lognu2nu1, ci
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map6Gabs,F_map10Gabs,spex_6_10
      integer :: m_x,m_y

      ci = abs(cos(incli))

      !///////////////////////////////////////////////////////
      !-------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 4.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      ancho_y = 2.d0*(1.+ci)*R0_proj/ci * 1.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction

      allocate( x_c(m_x), y_c(m_y), F_map6Gabs(m_x,m_y), F_map10Gabs(m_x,m_y), spex_6_10(m_x,m_y)  )


      !//////////////////////////////////////////////////////
      !------------------------------------------------------
      ! Data read

      open (17,file='SW_conv_6Gabs.dat')
      open (18,file='SW_conv_10Gabs.dat')
      open (19,file='SW_spex_6_10.dat')

      lognu2nu1 = log(10./6.)

      do i = 1, m_x
         do j = 1, m_y
            read(17,61) x_c(i),y_c(j),F_map6Gabs(i,j)
            read(18,61) x_c(i),y_c(j),F_map10Gabs(i,j)
            if ( F_map6Gabs(i,j) > 1.d-3 .and. F_map10Gabs(i,j) > 2.d-2 ) then  ! at least 20 uJy
               spex_6_10(i,j) = log(F_map10Gabs(i,j)/F_map6Gabs(i,j)) / lognu2nu1
               write(19,61) x_c(i),y_c(j),spex_6_10(i,j)
            endif
         enddo
         read(17,*) 
         read(18,*) 
         write(19,*) 
      enddo
      close(17)
      close(18)
      close(19)

      61  format(F11.4,X,F11.4,X,E11.4)

      deallocate (x_c,y_c,F_map6Gabs,F_map10Gabs,spex_6_10)

   END subroutine SW_spex_6_10_run

   !============================================================================
   ! SW_spex_6_10_run (END)
   !============================================================================


END module SW_conv
