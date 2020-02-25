module controles

implicit none
public

!===================================================
! Qué programas voy a correr
logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: CD_ = yes
logical, parameter :: termo = yes
logical, parameter :: absgg = no
logical, parameter :: absff = no!yes ! también calcula emi_ff(wind)
logical, parameter :: BBrad = no
logical, parameter :: Xrad = no ! de la WCR
logical, parameter :: rotar = yes

logical, parameter :: electrons = yes
logical, parameter :: dist1 = yes
logical, parameter :: dist2 = yes
logical, parameter :: radsyn1 = no!yes
logical, parameter :: radsyn2 = no!yes
logical, parameter :: radic1 = yes
logical, parameter :: radic2 = yes
logical, parameter :: radbr1 = yes
logical, parameter :: radbr2 = yes
logical, parameter :: abssyn = no!yes
logical, parameter :: mapas = no!yes

logical, parameter :: protons = yes
logical, parameter :: dist1p = yes
logical, parameter :: dist2p = yes
logical, parameter :: radpp1 = yes
logical, parameter :: radpp2 = yes

logical, parameter :: escribir_e = yes ! F de procesos leptónicos en la WCR
logical, parameter :: escribir_p = yes ! F de procesos hadrónicos en la WCR

logical, parameter :: convolucionar = no!yes
logical, parameter :: test = no!yes

!===================================================

END module controles
