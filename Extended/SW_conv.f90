module SW_conv
   ! This program makes the convolution of the emission in the plane of the sky 
   ! with the instrument response 
   use global
   use map_freq_type, only: map_freq_t
   implicit none

   public :: SW_conv_run, SW_spex_110_190_run, SW_spex_run, SW_spex_6_10_run

   contains

   !======================================================================================
   ! SW_conv_run
   !======================================================================================
   ! This routine makes a convolved map taking into account the reverse shock(s) emission,
   ! and the FS (HL and IL) emission. 
   ! mf = one entry of maps(:) (map_freqs.f90) or map_avg (controls.f90). It
   ! gives the file-name label (mf%label) and the resolution (mf%sigma).
   !--------------------------------------------------------------------------------------
   subroutine SW_conv_run(shock, mf)
      character(LEN=*), intent(in) :: shock
      type(map_freq_t), intent(in) :: mf
      character(len=:), allocatable :: freq
      integer, parameter :: m_tot = (2*(mphi-1))*(m_theta*2)   ! Total number of data
      real(dp) :: x_map(m_tot),y_map(m_tot)
      real(dp) :: F_abs(m_tot),F_un(m_tot)
      real(dp) :: sigma_x, sigma_y, ancho_x, ancho_y, D_x, D_y,ci
      real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
      real(dp), dimension(:), allocatable :: x_c,y_c
      real(dp), dimension(:,:), allocatable :: F_map_abs
      integer :: m_x,m_y
      integer, parameter:: m_map = mphi*mr_map
      real(dp) :: A_beam,A_cell, x_map2(m_map),y_map2(m_map),S_map2(m_tot), S_map3(m_tot)
      real(dp) :: w_x, w_y

      !----------------------------------------------------------------
      ! Obtain the resolution at the given frequency

      freq    = trim(mf%label)
      sigma_x = mf%sigma(1)
      sigma_y = mf%sigma(2)
      A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
      A_cell = f_Ny**2*sigma_x*sigma_y ! cell sample size

      ci = abs(cos(incli))
      open(188,file='R0_proj_'//system//'.dat')
         write(188,*) 0.d0,R0_proj/ci
      close(188)

      !---------------------------------------------------------------
      ! Initialize the data for the convolution region

      ancho_x = 8.d0*R0_proj/ci * 1.3d0 ! [mas](extra margin)
      if (freq == '110MHz') ancho_x = ancho_x * 4.08d0 ! [mas](extra margin)
      ancho_y = 4.d0*(1.+ci)*R0_proj/ci * 5.3d0 ! [mas](extra margin)

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
      m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
      x_c0 = -ancho_x/2.d0   ! Center of the first ellipse 
      y_c0 = -ancho_y*(3.-ci)/4.d0
      allocate( x_c(m_x), y_c(m_y), F_map_abs(m_x,m_y)  )

      !------------------------------------------------------
      ! Data read

      sum = 0.d0
      if (shock == 'BOTH') then
         open(13,file='sync_radiomap_'//freq//'.dat')
      else if (shock /= 'BOTH') then
         open(13,file='sync_radiomap_'//shock//'_'//freq//'.dat')
      end if
      do l = 1, m_tot
         read(13,*) y_map(l),x_map(l),F_un(l),F_abs(l)
         sum = sum + F_abs(l)
      enddo
      write(*,*) 'Pre-conv: sync'
      write(*,80) 'S('//freq//')=',sum,'mJy'
      close(13)

      sum = 0.d0
      if (shock == 'BOTH') then
         open(133,file='th_ff_map_'//freq//'.dat')
      else if (shock /= 'BOTH') then
         open(133,file='th_ff_map_'//shock//'_'//freq//'.dat')
      end if
      do l = 1, m_tot
         read(133,*) S_map2(l)
         sum = sum + S_map2(l)
      enddo
      write(*,*) 'Pre-conv: ff'
      write(*,80) 'S_ff('//freq//')=',sum,'mJy'
      close(133)

      if ((shock == 'BOTH') .or. (shock == 'FS')) then 
         sum = 0.d0
         open(233,file='hot_layer_'//freq//'.dat')
         do l = 1, m_tot
            read(233,*) S_map3(l)
            sum = sum + S_map3(l)
         end do
         write(*,*) 'Pre-conv: hot layer'
         write(*,80) 'S('//freq//')=',sum,'mJy'
      end if

      !------------------------------------------
      ! Make the convolution

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
                  D_F = ( F_abs(l) + S_map2(l) + S_map3(l) ) * atenua
               else
                  D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
               endif
               suma = suma + D_F
            enddo

            F_map_abs(i,j) = suma 
         enddo

      enddo

      !///////////////////////////////////////////////////////////////
      !
      ! Data storage + consistency check:

      mJy_beam = A_cell/A_beam ! Transform mJy/beam --> mJy (for consistency checks)

      sum = 0.d0
      open (17,file='SW_conv_'//freq//'abs_'//system//'.dat')
      do i = 1, m_x
         do j = 1, m_y
            write(17,60) x_c(i), y_c(j), F_map_abs(i,j)
            sum = sum + F_map_abs(i,j)*mJy_beam
         enddo
         write(17,*) 
      enddo
      write(*,*) 'Post-conv:'
      write(*,80) 'S('//freq//')=',sum,'mJy'
      write(*,*) ''
      close(17)

      60  format(F11.4,X,F11.4,X,E11.4)
      80  format(A11,X,Es10.2,X,A3)

      deallocate (x_c,y_c,F_map_abs)

   END subroutine SW_conv_run

   !============================================================================
   ! SW_conv_run (END)
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

      open (17,file='SW_conv_110Mabs_'//system//'.dat')
      open (18,file='SW_conv_190Mabs_'//system//'.dat')
      open (19,file='SW_spex_110_190_'//system//'.dat')

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

      open (17,file='SW_conv_1Gabs_'//system//'.dat')
      open (18,file='SW_conv_5Gabs_'//system//'.dat')
      open (19,file='SW_spex_'//system//'.dat')

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

      open (17,file='SW_conv_6Gabs_'//system//'.dat')
      open (18,file='SW_conv_10Gabs_'//system//'.dat')
      open (19,file='SW_spex_6_10_'//system//'.dat')

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