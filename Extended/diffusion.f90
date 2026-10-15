  module diffusion

  use constants
  use system_parameters
  implicit none
  public

  contains

  !============================================================================
  !  Diffusion coefficients
  !============================================================================
  
    function Diff(E,B,diff_regime)

      real(dp) :: Diff, D_1_GeV, delta, xhi
      real(dp), intent(in) :: E, B
      character(len=*), intent(in) :: diff_regime

      !character :: diff_regime --> Defined in system_parameters
      if (diff_regime == 'G') then       ! Gabici+ 2007
        delta=0.5d0
        xhi = 0.1d0
        D_1_GeV = 3.d27
        Diff = xhi*D_1_GeV*( E/(1.d9*eV) / (B/(3.d-6)) )**delta  !cm^2 s^-1
      else if (diff_regime == 'B') then  ! Bohm diff
        Diff = lambda_diff*E*c/(3.d0*qe*B) !cm^2 s^-1
      else 
        write(*,*) 'ERROR in Diff: unknown diffusion regime ', trim(diff_regime), ' (expected G or B)'
        error stop
      endif

     end Function Diff


END module diffusion
  !============================================================================