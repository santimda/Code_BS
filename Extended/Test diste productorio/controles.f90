module controls

implicit none
public

!===================================================
! Programs to run

logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: CD_calc = yes
logical, parameter :: termo_calc = yes

!----------

logical, parameter :: electrons = yes

logical, parameter :: diste_sw_calc = yes
logical, parameter :: radsyn_sw_calc = no!yes

!===================================================

END module controls
