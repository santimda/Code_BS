PROGRAM main

  use global
  use controls
  use CD, only: CD_run
  use CD_termo, only: CD_termo_run
  use CD_diste, only: CD_diste_run
  use CD_radsyn, only: CD_radsyn_run

  implicit none
  
  write(*,*) ''
  write(*,*) '============================'
  write(*,*) ' START'
  write(*,*) '============================'
  write(*,*) ''

  !---------------------------------------------------------------!
  !                                                               !
  !                  Shocks                                       !
  !                                                               !
  !---------------------------------------------------------------!

  if (CD_calc.eqv.yes) then
     write(*,*) 'Calculate the geometry of the shocks (CD position)'
     call CD_run  
     write(*,*) ''      
  endif

  if (termo_calc.eqv.yes) then
     write(*,*) 'Calculate the thermodynamical quantities at the shock'
     call CD_termo_run
     write(*,*) ''
  endif
  
  !---------------------------------------------------------------!
  !                                                               !
  !                  Electrons                                    !
  !                                                               !
  !---------------------------------------------------------------!

  IF (electrons.eqv.yes) then
     write(*,*) 'ELECTRONS:'
     
     if (diste_sw_calc.eqv.yes) then
        write(*,*) 'Calculate the non-thermal electron distribution in the SW'
        call CD_diste_run    
        write(*,*) ''
     endif

     if (radsyn_sw_calc.eqv.yes) then
        write(*,*) 'Calculate the (unabsorbed) synchrotron emission in the SW'
        call CD_radsyn_run           ! Just one linear emitter      
        write(*,*) ''
     endif

  ENDIF


  !================================================================

  write(*,*) '============================'
  write(*,*) ' END'
  write(*,*) '============================'
  
END PROGRAM main

   
    
