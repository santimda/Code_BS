module radadi

  private
  public :: radadi_run

contains
  
  !=====================================================================
  ! radadi_run
  !=====================================================================
  ! This module determines if the shock is either adiabatic or radiative
  ! based on Stevens et al, 1991.
  !=====================================================================

  subroutine radadi_run
  use constants
  use system_parameters
  implicit none

  real(dp) :: t_cool, t_adv
  real(dp) :: c_s   !sound speed
  real(dp) :: n_sh  !n_sh=4*n_pre_sh  number density of the wind at the shock
  real(dp) :: n_sh_f  !number density (forward shock)
  real(dp) :: n_sh_r  !number density (reverse shock)
  real(dp) :: quot
!  real :: tsw(vel) = 1.36d7*(vel/1.d8)**2 !temperature of the shock-heated wind
  real(dp) :: tsw_f = 1.36d7*(v_ISM/1.d8)**2
  real(dp) :: tsw_r = 1.36d7*(v_w/1.d8)**2
  real(dp) :: c_s_r = v_w/sqrt(8.)
  real(dp) :: c_s_f = v_ISM
  interface
     subroutine quotient(n_sh, tsw, c_s, quot)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: n_sh, tsw, c_s
        real(dp), intent(out) :: quot
      end subroutine quotient
   end interface

  !=====================================================================
  

  !=====================================================================
   !Forward shock
   
   n_sh_f = 4*n_ISM
  write(*,*) "Forward shock is:"
  call quotient(n_sh_f, tsw_f, c_s_f, quot)

  
  !====================================================================
  !   Reverse shock:
  
  n_sh_r = n_a! 4.*rho_w/(mp*mu_sw), ambient nuclei number density
 
  write(*,*) "Reverse shock is:"
  call quotient(n_sh_r, tsw_r, c_s_r, quot)

  
  !====================================================================

END subroutine radadi_run

  !====================================================================

END module radadi


!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                           SUBROUTINES                                   //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine quotient(n_sh, tsw, c_s, quot)   ! Relation between t_cool and t_adv
  use constants
  use system_parameters
    real(dp) :: n_sh
    real(dp) :: quot
    real(dp) :: lambda
    real(dp) :: tsw
    real(dp) :: c_s
    interface
       
     function lambda_t(tw)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp) :: lambda_t
        real(dp), intent(in) :: tw
      end function lambda_t
   end interface

   lambda = lambda_t(tsw)
 !  write(*,*) "lambda=", lambda
    
    t_cool = k*tsw/(4.*n_sh*lambda) !radiative losses
 !   write(*,*) "k=", k, "tsw=", tsw, "n_sh=", n_sh
 !   write(*,*) "t_cool=", t_cool
    t_adv = R0*pc/c_s
  !  write(*,*) "t_adv=", t_adv

    quot = t_cool/t_adv  !comparison between time scales

    if (quot .GT. 1) then
     write(*,*) 'adiabatic.'
    else if (quot .LT. 1) then
     write(*,*) 'radiative.'
    end if
  

   Return  
   END subroutine quotient
   


   !/////////////////////////////////////////////////////////////////////
   
Function lambda_t(tw)   ! Myasnikov et al, 1998
   use constants
   use system_parameters
   implicit none
   real(dp), intent(in) :: tw  !variable to calculate lambda
   real(dp) :: lambda_t
   

   write(*,*)
    if ((1.d4 .LE. tw) .AND. (tw .LE. 1.d5)) then
       lambda_t = 7.d-27 * tw

    else if ((1.d5 .LT. tw) .AND. (tw .LE. 4.d7)) then
       lambda_t = 7.d-19 *(tw)**(-0.6) 

    else if (tw .GT. 4.d7) then
       lambda_t = 3.d-27 * (tw)**(0.5)
    end if
    
    return

    End function lambda_t

    !/////////////////////////////////////////////////////////////////////
  
