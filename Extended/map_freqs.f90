module map_freqs
   !====================================================================
   ! maps(:) = frequencies of map_catalogue (controls.f90)
   ! All modules that make or read radio maps loop over maps(:).
   !====================================================================
   use constants, only: dp
   use map_freq_type, only: map_freq_t
   use controls, only: map_catalogue
   implicit none
   private
   public :: maps, maps_init, nu_index

   type(map_freq_t), allocatable :: maps(:)

contains

   subroutine maps_init()
      ! Builds maps from the frequency catalogue
      maps = pack(map_catalogue, map_catalogue%on)
   end subroutine maps_init

   ! Index of frequency f in a log-spaced grid with first value nu1,
   ! log-step dlnu = log(nu(2)/nu(1)) and n points
   elemental integer function nu_index(f, nu1, dlnu, n)
      real(dp), intent(in) :: f, nu1, dlnu
      integer,  intent(in) :: n
      nu_index = min(max(int(log(f/nu1)/dlnu) + 1, 1), n)
   end function nu_index

end module map_freqs