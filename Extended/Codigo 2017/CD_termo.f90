module CD_termo
  
  use global
  
  implicit none
  private
  public :: CD_termo_run
  
contains
  
  !============================================================================
  ! CD_termo_run
  !============================================================================
  subroutine CD_termo_run
    
    use global
    !$ use omp_lib
    implicit none
    
    !¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    ! Este programa se encarga de calcular las variables termodinámicas
    ! sobre una curva superpuesta a la CD (generada por CD.f90), a un lado
    ! y otro de la misma, generando termo1.dat y termo2.dat. También se encarga
    ! de calcular la luminosidad de inyección que corresponde a cada punto en cada
    ! choque.
    !________________________________________________________________________
    
    real(dp) :: R1(ml),R2(ml),Dy(ml),Dx(ml),Dl(ml)
    real(dp) :: dtheta(ml),dtheta_prom,normal(3),dOmega,perp(ml)
    real(dp) :: Beq0,rhow(ml),L_iny(ml),L_wind1,L_wind2,L_iny_tot
    allocate ( P(ml),rho(ml),T(ml),B(ml),v(ml),x(ml),y(ml) )
    
    !///////////////////////////////////////////////////////////////////////////
    !
    ! Leemos los puntos (x(l),y(l)) de la curva, siendo (x(1),y(1)) las coords 
    ! del apex de la CD
    
    open(11,file='CD.dat')
    do l=1,ml/2
       read(11,*) x(l),y(l) ! la rama de arriba de la trayectoria en cm
    enddo
    close(11)
    
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !/                                                                   C
    !/                          S1                                       C
    !/                                                                   C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !  Luminosidad de inyección:
    !
    !///////////////////////////////////////////////////////////////////

    L_wind1 =(Mdot1*vinf1**2)/2.d0
    Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
    
    do l=1,ml/2
       R1(l) = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
       if (l == 1) then
          Dx(1)=0.d0 ! arranca centrado
          Dy(1)=y(1)
       else
          Dx(l) = x(l)-x(l-1)
          Dy(l) = y(l)-y(l-1)
       endif
       Dl(l) = sqrt(Dx(l)**2+Dy(l)**2)
       dtheta(l) = atan(Dl(l)/R1(l))
    enddo
  
    L_iny_tot = 0.d0
    open(42,file='Liny1.dat')
    open(43,file='n_perp1.dat')
    do l=1,ml/2-1
       
       call vecnormal(Dx(l),Dx(l+1),Dy(l),Dy(l+1),Dl(l),Dl(l+1),normal)
       Xemi = (/ x(l),y(l),0.d0 /)
       call versorresta(Xemi,Xst1,3,XEe)
       dtheta_prom = 0.5d0*(dtheta(l)+dtheta(l+1))
       dOmega = (y(l)/R1(l))*dtheta_prom*delta_phi ! sen(o)*do*dphi
       perp(l) = abs(DOT_PRODUCT(XEe,normal))
       L_iny(l) = L_wind1*(dOmega/4.d0/pi)*perp(l)
       ! Solo me interesa el flujo de energía en la dirección normal al S1
       L_iny_tot = L_iny_tot + L_iny(l)
       
       write(42,*) x(l),y(l),L_iny(l),Dl(l)
       write(43,*) normal(1),normal(2),normal(3),perp(l)
    enddo
    perp(ml/2) = perp(ml/2-1)
    
    write(42,*) x(ml/2),y(ml/2),L_iny(ml/2-1),Dl(ml/2) ! no me complico en el borde
    close(42)    
    write(43,*) normal(1),normal(2),normal(3),perp(ml/2) ! no me complico en el borde
    close(43)
    write(*,99) 'L_iny(S1)=',L_iny_tot*mphi
    write(*,99) 'L_NT/L_W(%)=', frac_NT*L_iny_tot*mphi/L_wind1*100.d0


    !//////////////////////////////////////////////////////////////////////
    !
    !     Variables termodinámicas (S1)
    !
    !//////////////////////////////////////////////////////////////////////

    ! set omp threads number
    ! $ call omp_set_num_threads(nth)   

    !  $omp parallel do private(l) -> Necesito rho(1) para paralelizar =S
    do l=1,ml/2                             ! Calculo para la mitad de arriba
       rhow(l) = Mdot1/(4.d0*pi*R1(l)**2*vinf1)
       rho(l) = 4.d0*rhow(l)
       P(l) = 0.75d0*rhow(l)*vinf1**2*perp(l)     ! Valentí
       v(l) = vinf1*sqrt(1.d0-perp(l)**2) !sqrt(ga*P(1)/rho(1))! ga es el coeficiente adiabático
!       write(*,*) v(l)/sqrt(ga*P(1)/rho(1))
       T(l) = P(l)/(rho(l)*k/mp/mu_1)
       Beq0 = sqrt(8.d0*pi*ga*P(l)) !sqrt(4.d0*pi*rho(l)*v(l)**2.d0)
       B(l) = frac1*Beq0
    enddo
    !  $omp end parallel do
    do l=ml/2+1,ml                       ! Espejo para obtener la mitad de abajo
       rho(l)=rho(l-ml/2)
       P(l)=P(l-ml/2)
       v(l)=v(l-ml/2)
       T(l)=T(l-ml/2)
       B(l)=B(l-ml/2)
    enddo

    open(16,file='termo1.dat')
    open(17,file='termonorm1.dat')
    do l=1,ml
       write(16,*) l,rho(l),P(l),v(l),T(l),B(l) 
       ! Hago uno normalizado para poder hacer chequeos de consistencia:
       write(17,*) l,rho(l)/rho(1),P(l)/P(1),v(l)/vinf1,T(l)/T(1),B(l)/B(1) 
    enddo
    close(16)
    close(17)

    write(*,*) 'Para el apex de S1:'
    write(*,41) 'rho0=',rho(1)
    write(*,41) 'P0=',P(1)
    write(*,41) 'v0=',v(1)
    write(*,41) 'T0=',T(1)
    write(*,41) 'B0=',B(1)
    write(*,41) 'Beq0=',sqrt(8.d0*pi*ga*P(1))
    write(*,41) 'B_st1=',B(1)/4.d0*(x(1)/Rst1)*(vinf1/v_rot1)*1.d0 ! Asumo r_A = 1.0*R_st
41  format (A7,X,ES9.2)



    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !/                                                                   C
    !/                          S2                                       C
    !/                                                                   C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !  Luminosidad de inyección:
    !
    !///////////////////////////////////////////////////////////////////

    L_wind2 =(Mdot2*vinf2**2)/2.d0
    Xst2 = (/ D, 0.d0, 0.d0 /)

    do l=1,ml/2 
       R2(l) = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
!       if (l == 1) then -------> Ya está de antes
!          Dx(1)=x(1)
!          Dy(1)=y(1)
!       else
!          Dx(l) = x(l)-x(l-1)
!          Dy(l) = y(l)-y(l-1)
!       endif
!       Dl(l) = sqrt(Dx(l)**2+Dy(l)**2)
       dtheta(l) = atan(Dl(l)/R2(l))
    enddo
  
    L_iny_tot = 0.d0
    open(44,file='Liny2.dat')
    open(45,file='n_perp2.dat')
    do l=1,ml/2-1
       
       call vecnormal(Dx(l),Dx(l+1),Dy(l),Dy(l+1),Dl(l),Dl(l+1),normal)
       Xemi = (/ x(l),y(l),0.d0 /)
       call versorresta(Xemi,Xst2,3,XEe)
       dtheta_prom = 0.5d0*(dtheta(l)+dtheta(l+1))
       dOmega = (y(l)/R2(l))*dtheta_prom*delta_phi ! sen(o)*do*dphi
       perp(l) = abs(DOT_PRODUCT(XEe,normal))
       L_iny(l) = L_wind2*(dOmega/2.d0/pi)*perp(l)
       ! Solo me interesa el flujo de energía en la dirección normal al S2
       L_iny_tot = L_iny_tot + L_iny(l)
       
       write(44,*) x(l),y(l),L_iny(l),Dl(l)
       write(45,*) normal(1),normal(2),normal(3),perp(l)
    enddo
    perp(ml/2) = perp(ml/2-1)   
    
    write(44,*) x(ml/2),y(ml/2),L_iny(ml/2-1),Dl(ml/2) ! no me complico en el borde
    close(44)
    write(45,*) normal(1),normal(2),normal(3),perp(ml/2) ! no me complico en el borde
    close(45)
    write(*,99) 'L_iny(S2)=',L_iny_tot*mphi
    write(*,99) 'L_NT/L_W(%)=', frac_NT*L_iny_tot*mphi/ L_wind2*100.d0
99  format (A11,X,ES12.4)



    !//////////////////////////////////////////////////////////////////////
    !
    !     Variables termodinámicas (S2)
    !
    !//////////////////////////////////////////////////////////////////////
 
    ! Si estoy a la der de la CD la 2 es la más cercana
    !  $omp parallel do private(l) -> Necesito rho(1) para paralelizar =S
    do l=1,ml/2                             ! Calculo para la mitad de arriba
       rhow(l) = Mdot2/(4.d0*R2(l)**2*vinf2)
       rho(l) = 4.d0*rhow(l)
       P(l) = 0.75d0*rhow(l)*vinf2**2*perp(l)  ! Valentí
       v(l) = vinf2*sqrt(1.d0-perp(l)**2) !sqrt(ga*P(1)/rho(1))! ga es el coeficiente adiabático
!       write(*,*) v(l)/sqrt(ga*P(1)/rho(1))
       T(l)=P(l)/(rho(l)*k/mp/mu_1)
       Beq0 = sqrt(8.d0*pi*ga*P(l)) !sqrt(4.d0*pi*rho(l)*v(l)**2.d0)
       B(l)=frac1*Beq0
    enddo
    !  $omp end parallel do
    do l=ml/2+1,ml
       rho(l)=rho(l-ml/2)
       P(l)=P(l-ml/2)
       v(l)=v(l-ml/2)
       T(l)=T(l-ml/2)
       B(l)=B(l-ml/2)
    enddo
   
    open(18,file='termo2.dat')
    open(19,file='termonorm2.dat')
    do l=1,ml
       write(18,*) l,rho(l),P(l),v(l),T(l),B(l) 
       ! Hago uno normalizado para poder hacer chequeos de consistencia:
       write(19,*) l,rho(l)/rho(1),P(l)/P(1),v(l)/vinf2,T(l)/T(1),B(l)/B(1) 
    enddo
    close(18)
    close(19)

    write(*,*) 'Para el lado de la estrella 2:'
    write(*,41) 'rho0=',rho(1)
    write(*,41) 'P0=',P(1)
    write(*,41) 'v0=',v(1)
    write(*,41) 'T0=',T(1)
    write(*,41) 'B0=',B(1)
    write(*,41) 'Beq0=',sqrt(8.d0*pi*ga*P(1))
    write(*,41) 'B_st2=',B(1)/4.d0*((D-x(1))/Rst2)*(vinf2/v_rot2)*1.d0 ! Asumo r_A = 1*R_st

    deallocate ( P,rho,T,B,v,x,y )
  END subroutine CD_termo_run
  
END module CD_termo



!============================================================================

SUBROUTINE vecnormal(Dx1,Dx2,Dy1,Dy2,Dl1,Dl2,n)
  integer, parameter :: dp = kind(1.d0)
  real(dp) :: Dx1,Dx2,Dy1,Dy2,Dl1,Dl2
  real(dp) :: n1(3),n2(3),n(3)
  
  n1 = (/ -Dy1/Dl1,Dx1/Dl1,0.d0 /)
  
  n2 = (/ -Dy2/Dl2,Dx2/Dl2,0.d0 /)
  
  n = n1/Dl1 + n2/Dl2 ! Peso más al que está más cerca
  n = n/sqrt(n(1)**2+n(2)**2)
  
  RETURN
END SUBROUTINE vecnormal


!============================================================================
  
