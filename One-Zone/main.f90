PROGRAM main

   use global
   use controls
   use diste, only: diste_run
   use distp, only: distp_run
   use radbr, only: radbr_run
   use radic, only: radic_run
   use radsyn, only: radsyn_run
   use radpp, only: radpp_run
   use opac_ff, only: opac_ff_run
   use abssyn, only: abssyn_run
   use radadi, only: radadi_run
   
   implicit none  

   write(*,*) '\n ============================'
   write(*,*) ' START'
   write(*,*) '============================ \n'


   88    format(A11,X,Es10.2)
   !write(*,88) 'R_st[Ro]=',R_st/R_sun
   !write(*,88) 'L_st[Lo]=',L_st/L_sun
   !write(*,88) 'R0_pr[s]=',(R0/AU)/(distance/pc)*abs(cos(incli))
   write(*,88) 'R0[pc]=',R0/pc
   write(*,88) 'B[G]=',B
   write(*,88) 'alpha_ISM=',alpha_ISM
   !write(*,88) 'T_IR[K]=',T_IR
   !write(*,88) 'T_IR[K]=',T_IR
   !write(*,88) 'T_IR[K]=',16.4* ( (L_st/(4.*pi*R0**2.))/0.0217 )**(1./6.) ! Kobulnicky 2016
   !write(*,88) 'T_IR[K]=',T_st*(R_st/(2*R0))**(1./3.) ! Meyer 2016
   !write(*,810) 'L_inj =',L_inj

   ! ADD HERE A CHECK OF SHOCK NATURE (ADIABATIC OR RADIATIVE)
   !--------------------------------------------------------------
   
   if (check_shocks) call radadi_run

   !--------------------------------------------------------------

   IF (electrons) then
      write(*,*) '\n ELECTRONS:'

      if (diste_calc) call diste_run('RS')
      if (radbr_calc) call radbr_run('RS')

      if (radic_calc) then
         call radic_run(T_IR,Kappa_IR,iso,'IR','RS') ! IC with photons from the dust
         call radic_run(T_st,Kappa_st,iso,'st','RS') ! IC with stellar photons (the BS is extended, aniso in one-zone model makes no sense)
         call radic_run(T_CMB,Kappa_CMB,iso,'CMB','RS') ! IC with CMB photons, if relevant
      endif

      if (radsyn_calc) call radsyn_run('RS')
      if (opac_ff_calc) call opac_ff_run
      if (abssyn_calc) call abssyn_run('RS')

   ENDIF

   !--------------------------------

   IF (protons) then
      write(*,*) '\n PROTONS:'

      if (distp_calc) call distp_run('RS')      

      if (radpp_calc) call radpp_run('RS')

   ENDIF

   !================================================================


   write(*,*) '\n ============================'
   write(*,*) ' END'
   write(*,*) '============================'


END PROGRAM main



