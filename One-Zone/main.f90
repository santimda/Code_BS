PROGRAM main

   use global
   use controls
   use radadi, only: radadi_run
   !$ use omp_lib

   implicit none
   integer :: st
   interface
      subroutine calculate(shock) 
         character(LEN=*), intent(in) :: shock
      end subroutine calculate
   end interface


   write(*,*) '\n ======================================'
   write(*,*) '         START'
   write(*,*) '====================================== \n'

   ! OpenMP threads (set once): OMP_NUM_THREADS or nth
   call get_environment_variable('OMP_NUM_THREADS', status=st)
   !$ if (st /= 0) call omp_set_num_threads(nth)


   88    format(A11,X,Es10.2)
   !write(*,88) 'R_st[Ro]=',R_st/R_sun
   !write(*,88) 'L_st[Lo]=',L_st/L_sun
   !write(*,88) 'R0_pr[s]=',(R0/AU)/(distance/pc)*abs(cos(incli))
   write(*,88) 'R0[pc]=',R0/pc
   write(*,88) 'B_RS[G]=',B_RS
   write(*,88) 'B_FS[G]=',B_FS
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
      real(dp) :: E_NT_p, E_NT_e, E_NT
      real(dp) :: U_NT, U_B, eta_mag

825   format(A13,X,Es7.1)

      !--------------------------------

      IF (electrons) then
         write(*,*) '\n ELECTRONS:'

         if (diste_calc) call diste_run(shock, E_NT_e)
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

         if (distp_calc) call distp_run(shock, E_NT_p)      

         if (radpp_calc) call radpp_run(shock)

      ENDIF

      !--------------------------------

      write(*,*) '\n Calculate eta_mag = U_B/U_NT'
      
      E_NT = E_NT_e + E_NT_p

      U_NT = E_NT/(pi*R0**2*H0)      ! Define the NT energy density

      if (shock=='RS') then     
         U_B = B_RS**2/(8.*pi)      ! Define the magnetic energy density
      else if (shock=='FS') then    
         U_B = B_FS**2/(8.*pi)      ! Define the magnetic energy density
      end if
      
      eta_mag = U_B/U_NT            ! De Becker 2017  
      write(*,825) 'eta_mag =', eta_mag

      !--------------------------------

   end subroutine calculate

   !================================================================
