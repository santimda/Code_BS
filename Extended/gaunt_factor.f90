module gaunt_factor
  
  !==========================================================================
  ! Free-free Gaunt factor averaged over a Maxwellian, <g_ff>(u, gamma2),
  ! interpolated from the table of van Hoof et al. 2014 (gauntff.dat):
  !    u      = h*nu/(k*T)
  !    gamma2 = Zq**2*Ry/(k*T)
  ! Table: log10(u)      = -16 ... 13 (n_u  = 146 values, rows)
  !        log10(gamma2) =  -6 ... 10 (n_g2 =  81 values, columns)
  ! both with a step of 0.2 dex. The table is read once, on the first call.
  !==========================================================================
  use constants, only: dp, h, k, Ry
  implicit none
  private
  public :: gaunt_ff_calc

  integer,  parameter :: n_header = 16         ! header lines in gauntff.dat
  integer,  parameter :: n_u = 146, n_g2 = 81  ! table size
  real(dp), parameter :: log_u_min = -16.d0    ! first value of log10(u)
  real(dp), parameter :: log_g2_min = -6.d0    ! first value of log10(gamma2)
  real(dp), parameter :: dlog = 0.2d0          ! table step [dex]

  real(dp) :: gff(n_u, n_g2)                   ! gff(i_u, i_g2)
  logical  :: table_loaded = .false.

contains

  !============================================================================
  ! gaunt_ff_calc: gaunt_ff(i) = <g_ff>(h*nu_ff(i)/(k*T_l), Zq_l**2*Ry/(k*T_l))
  ! Bilinear interpolation in (log10 u, log10 gamma2). Stops if a value falls
  ! outside the table.
  !============================================================================

  subroutine gaunt_ff_calc(nu_ff, m_nu_ff, T_l, Zq_l, gaunt_ff)
    integer,  intent(in)  :: m_nu_ff
    real(dp), intent(in)  :: nu_ff(m_nu_ff), T_l, Zq_l
    real(dp), intent(out) :: gaunt_ff(m_nu_ff)
    integer  :: i, i_u, i_g
    real(dp) :: w_u, w_g

    if (.not. table_loaded) call load_table()

    ! gamma2 does not depend on nu: locate it once
    call locate(log10(Zq_l**2*Ry/(k*T_l)), log_g2_min, n_g2, 'gamma2', i_g, w_g)

    do i = 1, m_nu_ff
       call locate(log10(h*nu_ff(i)/(k*T_l)), log_u_min, n_u, 'u', i_u, w_u)
       gaunt_ff(i) = (1.d0-w_u)*(1.d0-w_g)*gff(i_u  ,i_g  ) &
                   +       w_u *(1.d0-w_g)*gff(i_u+1,i_g  ) &
                   + (1.d0-w_u)*      w_g *gff(i_u  ,i_g+1) &
                   +       w_u *      w_g *gff(i_u+1,i_g+1)
    end do

  end subroutine gaunt_ff_calc

  !============================================================================
  ! locate: for x = log10 of a table variable, returns the cell idx, such that
  ! x_min + (idx-1)*dlog <= x <= x_min + idx*dlog, and the weight w in [0,1]
  ! of the upper node (x = x_min + (idx-1+w)*dlog).
  !============================================================================

  subroutine locate(x, x_min, n, name, idx, w)
    real(dp), intent(in) :: x, x_min
    integer,  intent(in) :: n
    character(len=*), intent(in) :: name
    integer,  intent(out) :: idx
    real(dp), intent(out) :: w
    real(dp) :: s

    s = (x - x_min)/dlog          ! position in units of the step (0 at the first node)
    if (.not. (s >= 0.d0 .and. s <= real(n-1,dp))) then   ! also catches NaN
       write(*,*) 'ERROR in gaunt_ff_calc: log10(', name, ') = ', x, &
          ' is outside the gauntff.dat table'
       error stop
    end if
    idx = min(int(s) + 1, n-1)    ! idx+1 <= n also at the last node
    w = s - real(idx-1,dp)

  end subroutine locate

  !============================================================================
  ! load_table: reads gauntff.dat (input file, in the working directory)
  !============================================================================

  subroutine load_table()
    integer :: iu, i

    open(newunit=iu, file='gauntff.dat', status='old', action='read')
    do i = 1, n_header
       read(iu,*)
    end do
    do i = 1, n_u
       read(iu,*) gff(i,:)
    end do
    close(iu)
    table_loaded = .true.

  end subroutine load_table

END module gaunt_factor