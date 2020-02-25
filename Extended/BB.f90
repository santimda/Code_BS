module BB
  !Este módulo calcula la emisión de BB de las estrellas 1 y 2
  use global
  
  implicit none
  private
  public :: BB_run
  
contains
  
  !============================================================================
  ! BB_run
  !============================================================================
  
  subroutine BB_run  
    use global
    implicit none
    integer, parameter :: mE_ph = 200
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph)
    real(dp) :: L_BB,Rst2

    write(*,*) 'Calculate the BB thermal emission from the star'
    
    Eph_0 = 3.d0*k*T_st
    Eph_min = Eph_0/2.d6
    Eph_max = Eph_0*5.d1
    
    call vector_log(Eph_min,Eph_max,Eph_int,Eph)

    Rst2 = R_st**2.
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
