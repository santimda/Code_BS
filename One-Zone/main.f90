PROGRAM main

   use global
   use controls
   use radadi, only: radadi_run
   
   implicit none  

   write(*,*) '\n ======================================'
   write(*,*) '         START'
   write(*,*) '====================================== \n'


   88    format(A11,X,Es10.2)
   !write(*,88) 'R_st[Ro]=',R_st/R_sun
   !write(*,88) 'L_st[Lo]=',L_st/L_sun
   !write(*,88) 'R0_pr[s]=',(R0/AU)/(distance/pc)*abs(cos(incli))
   write(*,88) 'R0[pc]=',R0/pc
   write(*,88) 'B[G]=',B
   write(*,88) 'alpha_ISM=',alpha_ISM
   !write(*,88) 'T_IR[K]=',T_IR
   !write(*,88) 'T_IR[K]=',16.4* ( (L_st/(4.*pi*R0**2.))/0.0217 )**(1./6.) ! Kobulnicky 2016
   !write(*,88) 'T_IR[K]=',T_st*(R_st/(2*R0))**(1./3.) ! Meyer 2016

   !--------------------------------------------------------------
   
   if (check_shocks) call radadi_run

   !--------------------------------------------------------------

   write(*,*) '\n ----------------------------'
   write(*,*) ' Reverse Shock'
   write(*,*) '----------------------------'

   if (RS) call calculate('RS')

   write(*,*) '\n ----------------------------'
   write(*,*) ' Forward shock'
   write(*,*) '----------------------------'

   if (FS) call calculate('FS')

   write(*,*) '\n ======================================'
   write(*,*) '         END'
   write(*,*) '======================================'



END PROGRAM main


   !============================================================================


   subroutine calculate(shock) 
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

      character(LEN=*), intent(in) :: shock

      !--------------------------------

      IF (electrons) then
         write(*,*) '\n ELECTRONS:'

         if (diste_calc) call diste_run(shock)
         if (radbr_calc) call radbr_run(shock)

         if (radic_calc) then
            call radic_run(T_IR,Kappa_IR,iso,'IR',shock) ! IC with photons from the dust
            call radic_run(T_st,Kappa_st,iso,'st',shock) ! IC with stellar photons (the BS is extended, aniso in one-zone model makes no sense)
            call radic_run(T_CMB,Kappa_CMB,iso,'CMB',shock) ! IC with CMB photons, if relevant
         endif

         if (radsyn_calc) call radsyn_run(shock)
         if (opac_ff_calc) call opac_ff_run
         if (abssyn_calc) call abssyn_run(shock)

      ENDIF

      !--------------------------------

      IF (protons) then
         write(*,*) '\n PROTONS:'

         if (distp_calc) call distp_run(shock)      

         if (radpp_calc) call radpp_run(shock)

      ENDIF

      !--------------------------------

   end subroutine calculate

   !================================================================