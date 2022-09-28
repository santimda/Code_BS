module BB
  use global
  implicit none

  public :: BB_run
  
contains
  
  !============================================================================
  ! BB_run
  !============================================================================
  ! This module calculates the BB thermal emission from the star.
  !============================================================================
  
  subroutine BB_run  
    integer, parameter :: mE_ph = 200
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph)
    real(dp) :: L_BB,Rst2

    print*, 'Calculate the BB thermal emission from the star'
    
    Eph_0 = 2.7d0*k*T_st
    Eph_min = 10.d9*h    ! 10 GHz ~ Eph_0/1.d6
    Eph_max = Eph_0*5.d1
    
    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    Rst2 = R_st**2
    open(20,file='LBB.dat')
    do i=1,mE_ph
       L_BB = Rst2*nBB(Eph(i),T_st)
       write(20,*) Eph(i)/eV,Eph(i)*L_BB
    enddo
    close(20)

  END subroutine BB_run

  !============================================================================
  ! BB_run (END)
  !============================================================================
  
END module BB
