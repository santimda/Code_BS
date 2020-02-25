PROGRAM main

  use global
  use controls
  use diste, only: diste_run
  use distp, only: distp_run
  use radbr, only: radbr_run
  use radic, only: radic_run
  use radsyn, only: radsyn_run
  use radpp, only: radpp_run
  use wind_ff, only: wind_ff_run
  
  implicit none    
  integer :: i_IC
  real(dp) :: Ef_IC(mE_f)

!  interface
!     subroutine radic_run(T,Kappa,LIC_vec)!,n_vec)
!       use iso_fortran_env, only: dp => real64
!       implicit none
!       integer, intent(in) :: n_vec
!       real(dp), intent(in) :: T,Kappa
!       real(dp), intent(out) :: LIC_vec(:) !(n_vec)
!     end subroutine radic_run
!  end interface
  
  write(*,*) ''
  write(*,*) '============================'
  write(*,*) ' START'
  write(*,*) '============================'
  write(*,*) ''


88 format(A12,X,Es10.2)
!  write(*,88) 'R0[pc]=',R0/pc
  write(*,88) 'R0[cm]=',R0
  write(*,88) 'B[G]=',B
!  write(*,88) 'n_a[cm-3]=',n_a
  write(*,88) 'L_w =',L_wind
  write(*,88) 'U_mag[eV/cm3]=',B**2/(8.*pi)/eV
  write(*,88) 'U_st[eV/cm3]=',L_st/(4.*pi*R0**2.*c)/eV
!  write(*,88) 'U_CMB[eV/cm3]=',u_CMB/eV
!  write(*,88) 'F_X(WR)=',(L_st/1.d7)/(4.*pi*distance**2)

  write(*,*) ''
     
  if (ff_calc) call wind_ff_run 

  IF (electrons.eqv.yes) then
     write(*,*) 'ELECTRONS:'
     
     if (diste_calc) call diste_run 
     if (radsyn_calc) call radsyn_run
     if (radbr_calc) call radbr_run
     
     if (radic_calc.eqv.yes) then
        write(*,*) 'Calculate the (unabsorbed) IC-star emission'
        call radic_run(T_st,Kappa_st,Ef_IC,Lic)
        open(18,file='radic_st.dat')
        do i_IC = 1,mE_f
           write(18,*) Ef_IC(i_IC)/eV,Ef_IC(i_IC)*Lic(i_IC)
        enddo
        close(18)
        write(*,*) 'Calculate the (unabsorbed) IC-CMB emission'
        call radic_run(T_CMB,Kappa_CMB,Ef_IC,Lic)
        open(18,file='radic_cmb.dat')
        do i_IC = 1,mE_f
           write(18,*) Ef_IC(i_IC)/eV,Ef_IC(i_IC)*Lic(i_IC)
        enddo
        close(18)
     endif
     
  ENDIF

  !--------------------------------

  IF (protons.eqv.yes) then
     write(*,*) 'PROTONS:'
     
     if (distp_calc) call distp_run      
     if (radpp_calc.eqv.yes) call radpp_run
     
  ENDIF

  !================================================================

!  if (radsyn_calc .and. radic_calc) call fluxes_calc(EphSy,Lsy,EphIC,LIC,Flux)		

  write(*,*) '============================'
  write(*,*) ' END'
  write(*,*) '============================'
  
  
END PROGRAM main

   
    
