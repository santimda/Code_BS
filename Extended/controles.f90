module controls

use map_freq_type, only: map_freq_t
implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: verbose = yes

logical, parameter :: RS = yes
logical, parameter :: FS = yes
logical, parameter :: BOTH = yes    

logical, parameter :: BBrad_calc = no!yes
logical, parameter :: opacff_calc = yes  !also calculates emi_ff(wind), needs also electrons = yes

logical, parameter :: CD_calc = yes
logical, parameter :: termo_calc = yes
logical, parameter :: Xrad_calc = yes ! shock thermal emission
logical, parameter :: rotate = yes

!----------

logical, parameter :: electrons = yes
logical, parameter :: calc_N_bin = yes

logical, parameter :: diste_calc = yes
logical, parameter :: radsyn_sw_calc = yes
logical, parameter :: radbr_sw_calc = no
logical, parameter :: radic_st_calc = no
logical, parameter :: radic_IR_calc = no
logical, parameter :: abssyn_calc = yes !needs maps = yes
logical, parameter :: maps_calc = yes
logical, parameter :: convolution = yes! no
! Radio-map frequencies
! nu: freq in Hz
! label: str
! sigma: instrument resolution (x,y) ["]
type(map_freq_t), parameter :: map_catalogue(*) = [ &
   !           nu [Hz]   label      sigma (x,y)      make maps?
   map_freq_t( 1.10d8,  '110MHz', [ 4.0d0,  4.0d0], no ), &
   map_freq_t( 1.50d8,  '150MHz', [25.0d0, 25.0d0], yes ), &
   map_freq_t( 1.90d8,  '190MHz', [ 4.0d0,  4.0d0], no ), &
   map_freq_t( 3.25d8,  '325MHz', [54.0d0, 77.0d0], yes ), &
   map_freq_t( 1.42d9,  '1GHz',   [20.2d0, 12.5d0], no  ), &
   map_freq_t( 3.00d9,  '3GHz',   [20.2d0, 12.5d0], yes ), &
   map_freq_t( 4.86d9,  '5GHz',   [20.2d0, 12.5d0], yes ), &
   map_freq_t( 6.00d9,  '6GHz',   [16.0d0, 10.0d0], no ), &
   map_freq_t( 1.00d10, '10GHz',  [16.0d0, 10.0d0], yes ) ]

type(map_freq_t), parameter :: map_avg = &
   map_freq_t( 0.d0, 'avg', [19.7d0, 8.25d0], yes )

!----------

logical, parameter :: escribir_e = yes

!----------

logical, parameter :: protons = no

logical, parameter :: distp_calc = yes
logical, parameter :: radpp_calc = yes

logical, parameter :: escribir_p = yes

!----------

logical, parameter :: stromgren = no

!----------
logical, parameter :: rad_alpha = no
logical, parameter :: convolution_alpha = no

!===================================================

END module controls