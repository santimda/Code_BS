module controls

implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: ff_calc = no!yes

logical, parameter :: electrons = yes
logical, parameter :: protons = no!yes

logical, parameter :: diste_calc = yes
logical, parameter :: distp_calc = yes

logical, parameter :: radsyn_calc = yes
logical, parameter :: radbr_calc = no!yes
logical, parameter :: radic_calc = no!yes
logical, parameter :: radpp_calc = no!yes

!===================================================

END module controls
