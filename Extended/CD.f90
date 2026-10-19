module CD
   use global
   implicit none

   public :: CD_run

contains

   !============================================================================
   ! CD_run
   !============================================================================
   ! This module calulates the 2D curve that represents the position of the shock
   ! according to Christie+2016. Only the "upper" half is calculated, as there is
   ! azimuthal symmetry. The "rotated" positions are later calculated by CD_rot.f90
   !============================================================================

   subroutine CD_run
      real(dp), parameter :: theta_0 = 1.d-4   ! theta(1) = 0 yields numerical errors
      real(dp) :: dtheta, dR(m_theta), Length(m_theta), R_acel(m_theta), r_aprox(m_theta)
      integer  :: iu

      interface
         Subroutine RK(x1,y1,dx,y2)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp), intent(out) :: y2
            real(dp), intent(in) :: x1,y1,dx
         end Subroutine RK
      end interface

      allocate( x(m_theta), y(m_theta) )

      print*, 'Calculate the CD position'

      !-----------------------------------------------------------------------------
      !  Scheme: We calculate the (theta,r(theta)) positions of the CD. We build
      !  arrays theta(n),r(theta(n)), such that the star is located at r=0
      !  and r(0)=R0/R0=1 (theta(0)=0) is the stagnation point obtained from pressure balance.
      !  We assume that the stellar wind has reached its terminal velocity.
      !  The outputs are: a table with theta, r_norm, dR (CD.dat for the Christie+ solution,
      !  CD_aprox.dat for the W96 sol), and a table with cartesian coords x, y (CD_xy.dat)
      !-----------------------------------------------------------------------------

      80  format(A12,X,F5.2)
      90  format(A8,X,Es10.2)
      100 format(A17,X,F5.2,X,A12)
      110 format(A21,X,F5.2)

      ! Angular grid, from theta_0 to theta_max
      dtheta = (theta_max-theta_0)/real(m_theta-1,dp)
      theta = theta_0 + dtheta*[( real(i,dp), i=0,m_theta-1 )]

      ! CD shape: dr/dtheta from Eq. 11 of Christie+16 (function g), integrated
      ! with RK; each step depends on the previous one
      r_norm(1) = 1.d0
      do i = 2, m_theta
         call RK(theta(i-1),r_norm(i-1),dtheta,r_norm(i))
      end do
      R_theta = r_norm*R0

      ! Cell length (cosine theorem) and length along the CD from the apex
      dR(1) = R0*r_norm(1)*theta(1)
      dR(2:) = R0*sqrt( r_norm(1:m_theta-1)**2 + r_norm(2:)**2 - 2.d0*r_norm(2:)*r_norm(1:m_theta-1)*cos(dtheta) )
      Length = [( sum(dR(2:i)), i=1,m_theta )]

      write(*,100) 'R0_proj[arcmin] =', R0_proj/60./abs(cos(incli)), ' (no cos(i))'
      write(*,90) 'R0[pc] =', R0/pc
      write(*,80) 'Length[R0] =', Length(m_theta)/R0

      ! Cartesian coordinates (z = 0, rotated later), Wilkin (1996) solution for
      ! alpha = 0 (consistency check) and accelerator size for diffusive escape
      x = -r_norm*cos(theta)
      y =  r_norm*sin(theta)
      r_aprox = sqrt(3.d0*(1.d0-theta/tan(theta)))/sin(theta)
      R_acel = max(y, 1.d0)          ! = R0 near the apex, y further away

      open(newunit=iu, file=output_file('CD.dat'))
      write(iu,'(5es25.16e3)') (theta(l), r_norm(l), dR(l), Length(l)/R0, R_acel(l), l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('CD_aprox.dat'))
      write(iu,'(2es25.16e3)') (theta(l), r_aprox(l), l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('CD_xy.dat'))
      write(iu,'(2es25.16e3)') (x(l), y(l), l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('CD_xy_aprox.dat'))
      write(iu,'(2es25.16e3)') (-r_aprox(l)*cos(theta(l)), r_aprox(l)*sin(theta(l)), l=1,m_theta)
      close(iu)

      write(*,110) 'Length_proj[arcsec] =', (sqrt((x(m_theta)*R0+R0)**2+y(m_theta)**2)/AU)/(distance/pc)/abs(cos(incli))

      deallocate(x,y)

   End subroutine CD_run


End module CD


!============================================================================

Subroutine RK(xi,yi,dx,yf)
   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), intent(in) :: xi,yi,dx
   real(dp), intent(out) :: yf
   real(dp) :: y,k1,k2,k3,k4,h
   integer :: i
   integer, parameter :: n=1024

   dimension y(n)
   interface
      function g(a,b)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp) :: g
         real(dp), intent(in) :: a,b
      end function g
   end interface

   y(1) = yi
   h = dx/real(n-1,dp)   ! n points span dx in n-1 substeps
   do i = 1, n-1
      k1 = h*g(xi+(i-1)*h,y(i))
      k2 = h*g(xi+(i-1)*h+h/2.d0,y(i)+k1/2.d0)
      k3 = h*g(xi+(i-1)*h+h/2.d0,y(i)+k2/2.d0)
      k4 = h*g(xi+(i)*h,y(i)+k3)
      y(i+1) = y(i) + (k1 + 2.*k2 + 2.*k3 + k4)/6.d0
   end do
   yf = y(n)
End subroutine RK

!============================================================================

Function g(theta_i,r_i)
   use global
   implicit none
   real(dp) :: g
   real(dp),intent(in) :: theta_i,r_i
   real(dp) :: num,den,st,ct

   st = sin(theta_i)
   ct = cos(theta_i)

   num = theta_i*(1.-lambda*r_i**2) + r_i**2*(lambda-1.)*st*ct
   num = r_i*st*num

   den = (1.-lambda*r_i**2) * (theta_i - st*ct)*ct - st**3*(1.-r_i**2)

   g = num/den ! dr/dtheta
End function g

!============================================================================