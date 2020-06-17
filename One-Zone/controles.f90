module controls

implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: check_shocks = yes
logical, parameter :: RS = yes
logical, parameter :: FS = yes

logical, parameter :: electrons = yes
logical, parameter :: protons = yes

logical, parameter :: diste_calc = yes
logical, parameter :: distp_calc = yes

logical, parameter :: radsyn_calc = yes
logical, parameter :: radbr_calc = yes
logical, parameter :: radic_calc = yes
logical, parameter :: radpp_calc = yes

logical, parameter :: opac_ff_calc = yes
logical, parameter :: abssyn_calc = yes

!===================================================

END module controls
