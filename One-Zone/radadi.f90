module radadi

  private
  public :: radadi_run

contains

  !=======================================================================
  ! radadi_run
  !=======================================================================
  ! This module determines if the shock is either adiabatic or radiative
  ! based on Stevens et al. 1991.
  !=======================================================================

  subroutine radadi_run
    use global
    implicit none

    real(dp), parameter :: n_sh_f = 4.*n_ISM             !number density (forward shock)
    real(dp), parameter :: n_sh_r = 4.*rho_w/(mp*mu_sw)  !number density (reverse shock)
    real(dp), parameter :: c_s_f = v_ISM                 ! velocity of the shocked fluid (FS)
    real(dp), parameter :: c_s_r = v_w/sqrt(8.)          ! velocity of the shocked fluid (RS)
    real(dp), parameter :: T_sh_f = 1.36d7*(v_ISM/1.d8)**2
    real(dp), parameter :: T_sh_r = 1.36d7*(v_w/1.d8)**2
    real(dp) :: quot
    character*3 :: shock_type
    interface
      subroutine quotient(n_sh, T_sh, v_sh, quot, shock_type)
        use iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: n_sh, T_sh, v_sh
        real(dp), intent(out) :: quot
        character*3, intent(out) :: shock_type
      end subroutine quotient
    end interface

    !=====================================================================

    800 FORMAT (A17,X,A3,X,A15,Es10.2,A1)

    !=====================================================================
    ! Forward shock

    call quotient(n_sh_f, T_sh_f, c_s_f, quot, shock_type)
    write(*,800) "Forward shock is:",shock_type,'(t_cool/t_adv = ',quot,')'

    !====================================================================
    ! Reverse shock:

    call quotient(n_sh_r, T_sh_r, c_s_r, quot, shock_type)
    write(*,800) "Reverse shock is:",shock_type,'(t_cool/t_adv = ',quot,')'

  !====================================================================

  END subroutine radadi_run

!======================================================================

END module radadi


!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                           SUBROUTINES                                   //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine quotient(n_sh, T_sh, v_sh, quot, shock_type)   
  use global
  real(dp), intent(in) :: n_sh, T_sh, v_sh
  real(dp), intent(out) :: quot 
  character*3, intent(out) :: shock_type
  real(dp) :: t_cool, t_adv
  interface
    function lambda_T(T)
      use iso_fortran_env, only: dp => real64
      implicit none
      real(dp) :: lambda_T
      real(dp), intent(in) :: T
    end function lambda_T
  end interface

  t_cool = k*T_sh/(4.*n_sh*lambda_T(T_sh)) !radiative losses

  t_adv = R0/v_sh

  quot = t_cool/t_adv   !comparison between time scales

  if (quot > 1.d0) then
    shock_type = 'adi'
  else 
    shock_type = 'rad'
  end if

  Return  
END subroutine quotient


!/////////////////////////////////////////////////////////////////////

Function lambda_T(T)   ! Myasnikov et al, 1998
  use global
  implicit none
  real(dp), intent(in) :: T
  real(dp) :: lambda_T

  if ((1.d4 <= T) .AND. (T <= 1.d5)) then
    lambda_T = 7.d-27 * T

  else if ((1.d5 < T) .AND. (T <= 4.d7)) then
    lambda_T = 7.d-19 * T**(-0.6) 

  else if (T > 4.d7) then
    lambda_T = 3.d-27 * sqrt(T)

  else
    write(*,*) '\n ERROR IN CALCULATION OF lambda_T, TEMPERATURE OUT OF RANGE \n'
    lambda_T = 1.d0 ! Dummy value

  end if

  return

End function lambda_T

!/////////////////////////////////////////////////////////////////////

