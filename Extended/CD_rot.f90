module CD_rot
  !============================================================================
  ! Rotates the CD by an angle phi around the x-axis (the symmetry axis of
  ! the BS). Input: CD_xy.dat (upper half, z = 0, in units of R0, from CD.f90).
  ! Output: CDrot.dat with both halves, upper (y >= 0) and mirror (y <= 0),
  ! in units of R0. It is read by CD_radic, CD_abssyn and rad_H_alpha.
  ! Rotation of (x, y, 0) around the x-axis: (x, y*cos(phi), y*sin(phi)).
  !============================================================================
  use global
  use controls, only: verbose
  implicit none
  private
  public :: CD_rot_run

contains

  subroutine CD_rot_run(phi)
    real(dp), intent(in) :: phi
    real(dp) :: x_cd(m_theta), y_cd(m_theta)
    integer  :: iu

    if (verbose) print *, 'Rotate the coords'

    open(newunit=iu, file=output_file('CD_xy.dat'))
    read(iu,*) (x_cd(l), y_cd(l), l=1,m_theta)
    close(iu)

    open(newunit=iu, file=output_file('CDrot.dat'))
    ! upper half
    write(iu,'(3es25.16e3)') (x_cd(l),  y_cd(l)*cos(phi),  y_cd(l)*sin(phi), l=1,m_theta)
    ! lower half (mirror, y -> -y)
    write(iu,'(3es25.16e3)') (x_cd(l), -y_cd(l)*cos(phi), -y_cd(l)*sin(phi), l=1,m_theta)
    close(iu)

  END subroutine CD_rot_run

END module CD_rot