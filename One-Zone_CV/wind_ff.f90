module wind_ff
 
  use global
  
  implicit none
  private
  public :: wind_ff_run

  contains

    !============================================================================
    ! wind_ff_run
    !============================================================================
    ! Calculate f-f emission for ionised stellar winds
    ! Assume spherically symmetric stellar winds --> tau_ff(r), emi_ff(r)
    ! Output files:
  
    subroutine wind_ff_run

      use global   
      !$ use omp_lib
      implicit none
      real(dp) :: nu_ff(mE_ff)

      interface
        subroutine gaunt_ff_calc(nu_ff,Tw,Zw,gaunt)
          use global
          implicit none
          real(dp), intent(in) :: nu_ff(mE_ff),Tw,Zw
          real(dp), intent(out) :: gaunt(mE_ff)
        end subroutine gaunt_ff_calc
      end interface 

      interface
        subroutine wind_ff_calc(nu_ff)
          use global
          implicit none
          real(dp), intent(in) :: nu_ff(mE_ff)
        end subroutine wind_ff_calc
      end interface 

      interface
        subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
          use iso_fortran_env, only: dp => real64
          implicit none
          real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
          real(dp), intent(out) :: finterpol
        end subroutine interpol
      end interface


      !============================================================

      write(*,*) 'Calculate f-f emission of the stellar wind'
      
      ! Frequencies at which tau_ff is calculated
      call vector_log(nu_min_ff,nu_max_ff,nu_int_ff,nu_ff)

      call wind_ff_calc(nu_ff)     
      
    END subroutine wind_ff_run

  END module wind_ff


  !============================================================================
  ! gaunt_ff_calc 
  !============================================================================

  Subroutine gaunt_ff_calc(nu_ff,T_w,gaunt_ff)
    use global  
    implicit none
    real(dp), intent(in) :: nu_ff(mE_ff),T_w
    real(dp), intent(out) :: gaunt_ff(mE_ff)
    integer :: j_mu,j_g
    integer, parameter :: nu_mu = 146,nu_mga=81 ! # data in gauntff.dat
    real(dp) :: gff(nu_mu,nu_mga),gamma2,dg,dmu,mu_ff(mE_ff),gff_inter
    
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

    gamma2 = Zq_w**2*Ry/k/T_w
    j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
    dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0
 
    mu_ff = h*nu_ff/k/T_w   

    do i=1,mE_ff 
      ! index = (value - table_value(1))/(table_step) + 1, then do linear interpolation
      j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
      dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
      call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gff_inter)      
      gaunt_ff(i) = gff_inter
      gaunt_ff(i) = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w/nu_ff(i))) ! Approx.
    end do  

    RETURN

  End subroutine gaunt_ff_calc

  !============================================================================
  ! wind_ff_calc
  !============================================================================
  ! The output of the program is:
  ! windemiff.dat, containing: nu [Hz], S_nu [mJy]

  subroutine wind_ff_calc(nu_ff)
    use global
    implicit none
    real(dp),intent(in) :: nu_ff(mE_ff)
    real(dp) :: Mdot_eff,vinf,T_w
    real(dp) :: g_ff(mE_ff),S(mE_ff)

    Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)
    T_w = 0.4d0*T_st    ! assume T_w = cte for simplicity (the distances are r >> R_*)
    vinf=v_w

    !------------------------------------------------------------------------------   
    ! The total radio flux from the wind is given by Eq. 8 in Wright 1975:

    open(19,file='windemiff.dat')
    S = 0.d0

    call gaunt_ff_calc(nu_ff,T_w,g_ff)

    do i = 1, mE_ff     
      S(i) = (23.2d0/(distance/kpc)**2) * ( ((Mdot_eff/M_sun_yr)/&
        (vinf/1.d5)/mu_w)**2 *nu_ff(i)*gamma_w*g_ff(i)*Zq_w**2 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
      if (nu_ff(i) > 1.d-3*k*T_w/h) then
        S(i) = S(i)*exp(-h*nu_ff(i)/k/T_w)  ! Invented, formula valid for h*nu << k*T
      endif
      write(19,*) nu_ff(i),S(i)*1.d3 ! Stored in mJy
      !write(*,*) g_ff(i)
    end do
    close(19)
    
  RETURN

  END subroutine wind_ff_calc

  !===============================================================================

  Subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
    implicit none 
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
    real(dp), intent(out) :: finterpol
    real(dp) :: fx,fy

    fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
    fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
    finterpol = fx1y1 + fx*dx + fy*dy
    return
  end Subroutine interpol
    
  !===============================================================================
