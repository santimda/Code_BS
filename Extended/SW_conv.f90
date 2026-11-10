module SW_conv
   !============================================================================
   ! Convolution of the radio maps with the instrument response
   !============================================================================
   use global
   use controls, only: verbose
   use map_freq_type, only: map_freq_t
   implicit none
   private
   public :: SW_conv_run, SW_spex_run

   real(dp), parameter :: n_sig_margin = 3.d0      ! margin of the convolution grid around the data [sigma]

contains

   !============================================================================
   !        SW_conv_run
   !============================================================================
   ! Convolved map [mJy/beam] of the absorbed synchrotron + thermal (RS, FS
   ! and hot layer) emission. mf: map frequency (label and beam ["]).
   !============================================================================

   subroutine SW_conv_run(shock, mf, radiative)
      character(LEN=*), intent(in) :: shock
      logical, intent(in) :: radiative            ! radiative FS (hot layer)
      type(map_freq_t), intent(in) :: mf
      character(len=:), allocatable :: freq, tag
      real(dp) :: sigma_x, sigma_y, A_beam, A_cell, D_x, D_y, w_x, w_y, half_x, y_min, y_max
      integer  :: i, j, l, m_tot, m_x, m_y, iu
      logical  :: hot
      real(dp), allocatable :: x_map(:), y_map(:), F_un(:), F_abs(:), S_th(:), S_hot(:)
      real(dp), allocatable :: x_c(:), y_c(:), F_map_abs(:,:), EY(:,:)

      80 format(A,X,Es10.2,X,A)

      freq = trim(mf%label)
      sigma_x = mf%sigma(1)
      sigma_y = mf%sigma(2)
      A_beam = 2.d0*pi*sigma_x*sigma_y     ! Gaussian normalization factor
      A_cell = f_Ny**2*sigma_x*sigma_y     ! cell sample size
      hot  = radiative .and. ((shock == 'BOTH') .or. (shock == 'FS'))   ! the hot layer only exists in the FS

      if (shock == 'BOTH') then
         tag = ''
      else
         tag = shock//'_'
      end if

      open(newunit=iu, file=output_file('R0_proj_'//system//'.dat'))
      write(iu,'(2es17.8e3)') 0.d0, R0_proj/abs(cos(incli))
      close(iu)

      !------------------------------------------------------
      ! Maps (u -> y_map, v -> x_map)

      m_tot = count_records('sync_radiomap_'//tag//freq//'.dat')
      allocate( x_map(m_tot), y_map(m_tot), F_un(m_tot), F_abs(m_tot), S_th(m_tot), S_hot(m_tot) )

      open(newunit=iu, file=output_file('sync_radiomap_'//tag//freq//'.dat'))
      read(iu,*) (y_map(l), x_map(l), F_un(l), F_abs(l), l=1,m_tot)
      close(iu)

      open(newunit=iu, file=output_file('th_ff_map_'//tag//freq//'.dat'))
      read(iu,*) S_th
      close(iu)

      S_hot = 0.d0
      if (hot) then
         open(newunit=iu, file=output_file('hot_layer_'//freq//'.dat'))
         read(iu,*) S_hot
         close(iu)
      end if

      write(*,80) 'Pre-conv:  S_sync('//freq//') =', sum(F_abs), 'mJy'
      write(*,80) '           S_ff('//freq//')   =', sum(S_th),  'mJy'
      if (hot) write(*,80) '           S_hot('//freq//')  =', sum(S_hot), 'mJy'

      !---------------------------------------------------------------
      ! Convolution grid

      D_x = f_Ny*sigma_x
      D_y = f_Ny*sigma_y

      half_x = maxval(abs(x_map)) + n_sig_margin*sigma_x
      m_x = max( ceiling(2.d0*half_x/D_x), 1 )
      x_c = ([(real(i,dp), i=1,m_x)] - 0.5d0*(m_x+1))*D_x          ! symmetric with respect to x = 0

      y_min = minval(y_map) - n_sig_margin*sigma_y
      y_max = maxval(y_map) + n_sig_margin*sigma_y
      m_y = max( ceiling((y_max - y_min)/D_y), 1 )
      y_c = y_min + ([(real(j,dp), j=1,m_y)] - 0.5d0)*D_y

      if (verbose) write(*,'(a,i0,a,i0)') ' Convolution grid: m_x = ', m_x, ', m_y = ', m_y

      !------------------------------------------
      ! Convolution

      w_x = 2.d0*sigma_x**2
      w_y = 2.d0*sigma_y**2

      allocate( EY(m_tot,m_y), F_map_abs(m_x,m_y) )
      EY = exp( -(spread(y_map,2,m_y) - spread(y_c,1,m_tot))**2/w_y )
      call conv_rows(x_map, F_abs + S_th + S_hot, x_c, w_x, EY, F_map_abs)
      deallocate( EY )

      !------------------------------------------
      ! Save

      open(newunit=iu, file=output_file('SW_conv_'//freq//'abs_'//system//'.dat'))
      do i = 1, m_x
         write(iu,'(F11.4,X,F11.4,X,ES12.4E3)') (x_c(i), y_c(j), F_map_abs(i,j), j=1,m_y)
         write(iu,*)
      end do
      close(iu)

      write(*,80) 'Post-conv: S('//freq//')      =', sum(F_map_abs)*A_cell/A_beam, 'mJy'
      write(*,*) ''

      deallocate( x_map, y_map, F_un, F_abs, S_th, S_hot, x_c, y_c, F_map_abs )

   END subroutine SW_conv_run

   !============================================================================
   !        SW_spex_run
   !============================================================================
   ! Spectral index map between two convolved maps (same beam)
   !============================================================================

   subroutine SW_spex_run(mf1, mf2)
      type(map_freq_t), intent(in) :: mf1, mf2
      real(dp), parameter :: F1_min = 1.d-3, F2_min = 2.d-2    ! [mJy/beam]
      character(len=:), allocatable :: file1, file2, line
      real(dp) :: x1, y1, F1, x2, y2, F2, lognu2nu1
      integer  :: iu1, iu2, iu3, ios1, ios2
      character(len=256) :: buf1, buf2

      file1 = 'SW_conv_'//trim(mf1%label)//'abs_'//system//'.dat'
      file2 = 'SW_conv_'//trim(mf2%label)//'abs_'//system//'.dat'
      if (any(abs(mf1%sigma - mf2%sigma) > 1.d-6*mf1%sigma)) then
         print *, 'WARNING (SW_spex_run): different beams for '//trim(mf1%label)//' and '//trim(mf2%label)//'; no map'
         return
      end if

      lognu2nu1 = log(mf2%nu/mf1%nu)

      open(newunit=iu1, file=output_file(file1))
      open(newunit=iu2, file=output_file(file2))
      open(newunit=iu3, file=output_file('SW_spex_'//trim(mf1%label)//'_'//trim(mf2%label)//'_'//system//'.dat'))
      do
         read(iu1,'(a)',iostat=ios1) buf1
         read(iu2,'(a)',iostat=ios2) buf2
         if (ios1 /= 0 .or. ios2 /= 0) exit
         if (len_trim(buf1) == 0) then         ! end of a block
            write(iu3,*)
            cycle
         end if
         line = buf1; read(line,*) x1, y1, F1
         line = buf2; read(line,*) x2, y2, F2
         if (abs(x1-x2) > 1.d-3 .or. abs(y1-y2) > 1.d-3) then
            print *, 'WARNING (SW_spex_run): the grids of '//file1//' and '//file2//' differ'
            exit
         end if
         if (F1 > F1_min .and. F2 > F2_min) write(iu3,'(F11.4,X,F11.4,X,ES12.4E3)') x1, y1, log(F2/F1)/lognu2nu1
      end do
      close(iu1)
      close(iu2)
      close(iu3)

   END subroutine SW_spex_run


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

   !============================================================================
   ! Number of non-empty lines of a file

   integer function count_records(name)
      character(len=*), intent(in) :: name
      character(len=512) :: buf
      integer :: iu, ios
      count_records = 0
      open(newunit=iu, file=output_file(name), status='old')
      do
         read(iu,'(a)',iostat=ios) buf
         if (ios /= 0) exit
         if (len_trim(buf) > 0) count_records = count_records + 1
      end do
      close(iu)
   end function count_records

END module SW_conv