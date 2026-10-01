module map_freq_type
   !====================================================================
   ! Defines the data of one radio-map frequency. The list of frequencies (map_catalogue) is set in controls.f90
   !====================================================================
   use constants, only: dp
   implicit none
   private
   public :: map_freq_t

   type :: map_freq_t
      real(dp)         :: nu          ! frequency [Hz]
      character(len=8) :: label       ! file-name suffix: '110MHz', '1GHz', ...
      real(dp)         :: sigma(2)    ! instrument resolution (x,y) for SW_conv
      logical          :: on          ! yes = make maps at this frequency
   end type map_freq_t

end module map_freq_type