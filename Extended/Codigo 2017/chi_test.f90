module CHI_TEST
  !Este módulo junta todo en una sola SED
  use global
  
  implicit none
  private
  public :: chi_test_run
  
contains
  
  !============================================================================
  ! chi_test_run
  !============================================================================
  
  subroutine chi_test_run  
    use global
    implicit none
    Integer, parameter :: N_obs=3 ! 1.0D
!    Integer, parameter :: N_obs=4 ! 1.5D
    Real(dp) :: E_sy(mnu),flujo_mod(N_obs),F_mod(mnu),dummy,flux,chi2
    Real(dp) :: log_nu_int, nu_obs(N_obs),flujo_obs(N_obs),sigma_obs(N_obs)
    Integer :: i_obs(N_obs)
    allocate ( nu(mnu) )
    flux = 1.d26/(4.d0*pi*(distance*3.0856d21)**2)! luminosidad a flujo en mJy

    !_____________________________________________
    !/////
    !/////   Primero leo los datos observacionales

    open(15,file='data_radio2008y2009.dat')
!    open(15,file='data_radio2003y2004.dat')
    do j = 1,N_obs
       read(15,*) nu_obs(j),flujo_obs(j),dummy,sigma_obs(j)
    enddo
    close(15)
    nu_obs = nu_obs*1.d9 ! Lo paso a Hz (estaba todo en GHz)

    !_____________________________________________
    !/////
    !/////    Ahora leo los datos de mi modelo

    open(12,file='nu.dat')
    do i=1,mnu
       read(12,*) nu(i)
    enddo
    close(12)
    log_nu_int = log(nu(2)/nu(1))
    
    do j=1,N_obs
       i_obs(j) = int(log(nu_obs(j)/nu(1))/log_nu_int) + 1
!       write(*,*) i_obs(j),nu(i_obs(j))/1.d9
    enddo
    
    open (99,file='radsyn_WCR.dat')
    do j=1,mnu
       read(99,*) E_sy(j),dummy,F_mod(j) ! E_sy(j)[eV] = h*nu(j)*eV
       F_mod(j) = F_mod(j)*flux/nu(j)
    enddo
    close(99)
    
    !_____________________________________________
    !/////
    !/////    Ahora interpolo los flujos en las frecuencias de obs

    chi2=0.d0
    do j=1,N_obs
       flujo_mod(j) = F_mod(i_obs(j)) + ( F_mod(i_obs(j)+1) - F_mod(i_obs(j)) ) &
            /(nu(i_obs(j)+1)-nu(i_obs(j))) * ( nu_obs(j) - nu(i_obs(j)) )
       chi2 = chi2 + ( (flujo_mod(j)-flujo_obs(j))/sigma_obs(j) )**2.d0 
       write(*,*) flujo_mod(j),flujo_obs(j)
    enddo
    chi2 = chi2/N_obs 
    write(*,*) 'chi2=',chi2,'chi=',sqrt(chi2)
    

  END subroutine chi_test_run

  !============================================================================
  ! chi_test_run (END)
  !============================================================================
  
 
END module CHI_TEST
 
