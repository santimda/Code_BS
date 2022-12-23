  module gaunt_factor
    use global  
    implicit none
    public

  contains

  !============================================================================
  !  Gaunt factor calculation
  !============================================================================

  Subroutine gaunt_ff_calc(nu_ff,m_nu_ff,T_l,Zq_l,gaunt_ff)
    integer, intent(in) :: m_nu_ff
    real(dp), intent(in) :: nu_ff,T_l,Zq_l
    real(dp), intent(out) :: gaunt_ff
    integer, parameter :: nu_mu = 146, nu_mga=81 ! # data in gauntff.dat
    real(dp) :: gff(nu_mu,nu_mga),gamma2,dg,dmu,gff_inter
    integer :: j_mu,j_g
    real(dp), allocatable :: mu_ff(:)
    DIMENSION nu_ff(m_nu_ff), gaunt_ff(m_nu_ff)

    interface
    subroutine interpol2(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
      real(dp), intent(out) :: finterpol
    end subroutine interpol2
    end interface
    interface
    subroutine extrapol2(dx,dy,fx1y1,fx2y1,fx1y2,fextrapol)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
      real(dp), intent(out) :: fextrapol
    end subroutine extrapol2
    end interface  

    !--------------------------------------------------------------------------
    ! gauntff.dat (van Hoof 2014) allows to compute the Gaunt factor:
    ! <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, step 0.2
      
    open(13,file='gauntff.dat')
    do i=1,16
      read(13,*)
    end do
    do i=1,nu_mu
      read(13,*) (gff(i,j), j=1,nu_mga)
    end do
    close(13)

    gamma2 = Zq_l**2*Ry/k/T_l
    j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
    dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

    allocate( mu_ff(m_nu_ff) )
    mu_ff = h*nu_ff/(k*T_l)

    do i=1,m_nu_ff 
      !if (h*nu_ff(i) > 3*k*T_l) print*, 'WARNING: h nu > kT in opac_ff (gaunt_ff_calc)'
      ! index = (value - table_value(1))/(table_step) + 1, then do linear interpolation

      j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
      dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp

      if (j_mu < nu_mga) then
        call interpol2(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gff_inter)
      else
        call extrapol2(dg,dmu,gff(j_g,nu_mga-1),gff(j_g+1,nu_mga-1),gff(j_g,nu_mga),gff_inter)
      endif
      gaunt_ff(i) = gff_inter
      !gaunt_ff(i) = 9.77d0*(1.d0+0.13d0*log10(T_l**1.5d0/Zq_l/nu_ff(i))) ! Approx.
      !print*, i, h*nu_ff(i)/k/T_l, gff_inter, gaunt_ff(i), gff_inter/gaunt_ff(i)
    end do  

    deallocate( mu_ff )

  End subroutine gaunt_ff_calc

!============================================================================

END module gaunt_factor


!////////////////////////////////////////////////////////////////////////////////

Subroutine interpol2(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
  real(dp), intent(out) :: finterpol
  real(dp) :: fx,fy
  fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
  fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
  finterpol = fx1y1 + fx*dx + fy*dy
END Subroutine interpol2

!////////////////////////////////////////////////////////////////////////////////

Subroutine extrapol2(dx,dy,fx1y1,fx2y1,fx1y2,fextrapol)
  implicit none 
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
  real(dp), intent(out) :: fextrapol
  real(dp) :: fx,fy
  fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
  fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
  fextrapol = fx1y2 + fx*dx + fy*dy
END Subroutine extrapol2
  