module WCR_thermal
  ! Este programa calcula la emisión de bremsstrahlung térmico originada en la WCR
  use global
  
  implicit none
  private
  public :: WCR_thermal_run
  
contains
  
  !============================================================================
  ! WCR_thermal_run
  !============================================================================
  
  subroutine WCR_thermal_run  
    use global
    implicit none
    integer, parameter :: mE_ph2 = 100*mE_ph
    real(dp) :: g_ff,ne,ni,nmedio(ml),dr,r,emi_ff,Dl(ml),dV
    real(dp) :: Eph_0, Eph_min, Eph_max, Eph_int, Eph(mE_ph2)
    real(dp) :: L1_X(mE_ph2), L2_X(mE_ph2), L_X(mE_ph2)
    integer :: nl
    allocate ( x(ml),y(ml),rho(ml),P(ml),v(ml),T(ml),B(ml) )
    L1_X = 0.d0
    L2_X = 0.d0
    L_X = 0.d0

    open(12,file='CD.dat')
    do l=1,ml/2
       read(12,*) x(l),y(l)   ! la rama de arriba de la trayectoria en cm
       if (l == 1) then
          Dl(l) = y(1)
       else
          Dl(l) = sqrt((x(l)-x(l-1))**2 + (y(l)-y(l-1))**2)
       endif
    enddo
    close(12)

    Eph_0 = 3.d0*k*(T1*1.d3)  ! lo uso de referencia
    Eph_min = Eph_0/5.d4
    Eph_max = Eph_0*1.d1
    
    call vectorE(Eph_min,Eph_max,Eph_int,Eph)

    !///////////////////////////////////////////////////////////////////////////
    !     
    !    Calculamos primero para S1
    !
    !///////////////////////////////////////////////////////////////////////////
    
    open(11,file='termo1.dat')
    do l=1,ml
       read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
       nmedio(l) = rho(l)/mp/mu_i1 ! da del orden de 10^4 (lejos) - 10^6 (apex)
    enddo
    close(11)
    
    !     Calculamos la emisividad a partir de la ec. 5.15b del Rybicki:
    
    do l=1,ml/2

       ni = nmedio(l) 
       ne = gamma_1*nmedio(l) !
       
       r = sqrt(x(l)**2+y(l)**2) ! Dist a la estrella 1
       dr = 0.1d0*r               ! Tamaño lineal de la región de colisión 
       
       do i=1,mE_ph2
          emi_ff = 6.8d-38*Zq_1**2*ne*ni*T(l)**(-0.5d0)*exp(-Eph(i)/(k*T(l)))*g_ff(Eph(i),T(l)) 
          ! erg/s/cm3/Hz
          dV = ( 2.d0*pi*sqrt( y(l)**2 + (x(l)-x(1))**2 ) ) *Dl(l)*dr
          L1_X(i) = L1_X(i) + emi_ff*(Eph(i)/h)*dV
       enddo
       
    enddo
!    L1_X = 2.d0*L1_X! Agrego la mitad de abajo # Ya está contemplada con el 2pi


    !///////////////////////////////////////////////////////////////////////////
    !     
    !    Calculamos primero para S2
    !
    !///////////////////////////////////////////////////////////////////////////
    
    open(11,file='termo2.dat')
    do l=1,ml
       read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
       nmedio(l) = rho(l)/mp/mu_i2 ! da del orden de 10^4 (lejos) - 10^6 (apex)
    enddo
    close(11)
    
    !     Calculamos la emisividad a partir de la ec. 5.15b del Rybicki:
    
    do l=1,ml/2

       ni = nmedio(l) 
       ne = gamma_2*nmedio(l)
       
       r = sqrt((D-x(l))**2+y(l)**2) ! Dist a la estrella 1
       dr = 0.1d0*r               ! Tamaño lineal de la región de colisión 
 
       do i=1,mE_ph2
          emi_ff = 6.8d-38*Zq_2**2*ne*ni*T(l)**(-0.5d0)*exp(-Eph(i)/(k*T(l)))*g_ff(Eph(i),T(l)) 
          ! erg/s/cm3/Hz
          dV = ( 2.d0*pi*sqrt( y(l)**2 + (x(l)-x(1))**2 ) ) *Dl(l)*dr
          L2_X(i) = L2_X(i) + emi_ff*(Eph(i)/h)*dV
       enddo
       
    enddo
!    L2_X = 2.d0*L2_X ! Agrego la mitad de abajo# Ya está contemplada con el 2pi



    !////////////////// Grabo todo:


    open(19,file='WCR_ff.dat')
    do i=1, mE_ph2
       L_X(i) = L1_X(i) + L2_X(i)
       write(19,*) Eph(i)/eV,L1_X(i),L2_X(i),L_X(i)
    enddo
    close(19)
      
    
      
    deallocate ( x,y,rho,P,v,T,B )

 
  END subroutine WCR_thermal_run

  !============================================================================
  ! WCR_thermal_run (END)
  !============================================================================
  

END module WCR_THERMAL



 
   !-----------------------------------------------------------------
  
  Function g_ff(Eph_X,T_X)
    use global
    real(dp) :: g_ff,Eph_X,T_X,zeda
    g_ff = 0.d0
    if (Eph_X/(k*T_X) >= 1.d0 .and. Eph_X/(k*T_X) < 1.d3) then
       g_ff = sqrt( (3.d0/pi) * ((k*T_X)/Eph_X) )
    else if (Eph_X/(k*T_X) < 1.d0 .and. Eph_X/(k*T_X) >= 1.d-3 ) then
       zeda = 1.781d0  ! Es el gamma en el astrophysical formulae
       g_ff = (sqrt(3.d0)/pi) * log(4.d0*k*T_X/(zeda*Eph_X))
    else if (Eph_X/(k*T_X) > 1.d3 .or. Eph_X/(k*T_X) < 1.d-3) then
       g_ff = 0.d0  ! FRUTAAAAAAAAAA     
    endif

    RETURN
  end Function g_ff
  
  !-----------------------------------------------------------------
  
