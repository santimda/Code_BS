module CD_rot
  !============================================================================
  ! Rotates the CD by an angle phi around the x-axis (the symmetry axis of
  ! the BS). Input: CD_xy.dat (upper half, z = 0, in units of R0, from CD.f90).
  ! Output: CDrot.dat with both halves, upper (y >= 0) and mirror (y <= 0),
  ! in units of R0. It is read by CD_radic, CD_abssyn and rad_H_alpha.
  !============================================================================
  use global
  use controls, only: verbose
  implicit none
  private
  public :: CD_rot_run

contains

  !============================================================================
  ! CD_rot_run
  !============================================================================

  subroutine CD_rot_run(phi)
    real(dp), intent(in) :: phi
    real(dp) :: x_cd(m_theta), y_cd(m_theta)
    real(dp) :: xrot, yrot, zrot
    integer  :: iu

    if (verbose) print *, 'Rotate the coords'

    open(newunit=iu, file=output_file('CD_xy.dat'))
    do l = 1, m_theta
       read(iu,*) x_cd(l), y_cd(l)
    end do
    close(iu)

    open(newunit=iu, file=output_file('CDrot.dat'))
    do l = 1, m_theta                ! upper half
       call rotate_vector(x_cd(l), y_cd(l), phi, xrot, yrot, zrot)
       write(iu,*) xrot, yrot, zrot
    end do
    do l = 1, m_theta                ! lower half (mirror, y -> -y)
       call rotate_vector(x_cd(l), -y_cd(l), phi, xrot, yrot, zrot)
       write(iu,*) xrot, yrot, zrot
    end do
    close(iu)

  END subroutine CD_rot_run

  !============================================================================
  ! rotate_vector: rotates the vector (Vx, Vy, 0) by an angle ang around the
  ! x-axis
  !============================================================================

  subroutine rotate_vector(Vx, Vy, ang, VRx, VRy, VRz)
    real(dp), intent(in)  :: Vx, Vy, ang
    real(dp), intent(out) :: VRx, VRy, VRz

    VRx = Vx                ! unaltered
    VRy = cos(ang)*Vy       ! Vz = 0
    VRz = sin(ang)*Vy       ! Vz = 0

  END subroutine rotate_vector

END module CD_rot