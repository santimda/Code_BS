module controls

implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: BBrad_calc = no!yes
logical, parameter :: opacff_calc = no!yes ! also calculates emi_ff(wind), needs
!also electrons = yes

logical, parameter :: CD_calc = yes
logical, parameter :: termo_calc = yes
logical, parameter :: Xrad_calc = no!yes ! shock thermal emission
logical, parameter :: rotate = yes

!----------

logical, parameter :: electrons = yes

logical, parameter :: diste_sw_calc = yes
logical, parameter :: radsyn_sw_calc = yes
logical, parameter :: radbr_sw_calc = yes
logical, parameter :: radic_st_calc = yes
logical, parameter :: radic_IR_calc = yes
logical, parameter :: abssyn_calc = yes
logical, parameter :: maps_calc = no!yes
logical, parameter :: convolution = no!yes

logical, parameter :: escribir_e = yes

!----------

logical, parameter :: protons = no!yes

logical, parameter :: distp_calc = yes
logical, parameter :: radpp_calc = yes

logical, parameter :: escribir_p = yes

!===================================================

END module controls
