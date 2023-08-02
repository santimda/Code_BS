module controls

implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: RS = yes
logical, parameter :: FS = yes
logical, parameter :: BOTH = yes    

logical, parameter :: BBrad_calc = no!yes
logical, parameter :: opacff_calc = yes  !also calculates emi_ff(wind), needs
!also electrons = yes
integer, parameter :: t_err = 4

logical, parameter :: CD_calc = yes
logical, parameter :: termo_calc = yes
logical, parameter :: Xrad_calc = yes ! shock thermal emission
logical, parameter :: rotate = yes

!----------

logical, parameter :: electrons = yes

logical, parameter :: diste_calc = yes
logical, parameter :: radsyn_sw_calc = yes
logical, parameter :: radbr_sw_calc = no
logical, parameter :: radic_st_calc = no
logical, parameter :: radic_IR_calc = no
logical, parameter :: abssyn_calc = no !needs maps = yes
logical, parameter :: maps_calc = no
logical, parameter :: convolution = no

logical, parameter :: escribir_e = yes

!----------

logical, parameter :: protons = yes

logical, parameter :: distp_calc = yes
logical, parameter :: radpp_calc = yes

logical, parameter :: escribir_p = yes

!----------

logical, parameter :: stromgren = no

!----------
logical, parameter :: rad_alpha = yes
logical, parameter :: convolution_alpha = yes

!===================================================

END module controls
