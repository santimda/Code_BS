  module diffusion

  use constants
  use system_parameters
  implicit none
  public

  contains

  !============================================================================
  !  Diffusion coefficients
  !============================================================================
  
    function Diff(E,B)

      real(dp) :: Diff, D_10_GeV, delta
      real(dp), intent(in) :: E, B

      !character :: diff_regime --> Defined in system_parameters
      if (diff_regime == 'Valle') then       ! del Valle+2018
        delta=0.5d0
        D_10_GeV = 1.d28
        Diff = lambda_diff*D_10_GeV*(E/(1.d9*eV))**delta  !cm^2 s^-1
      else if (diff_regime == 'Bohm') then  ! Bohm diff
        Diff = lambda_diff*E*c/(3.d0*qe*B) !cm^2 s^-1
      else 
        write(*,*) 'ERROR, UNKNOWN DIFFUSION'
        stop
      endif

     end Function Diff


END module diffusion
  !============================================================================