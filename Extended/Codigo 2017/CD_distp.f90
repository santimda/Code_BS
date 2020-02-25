module CD_distp
  
  use global
  
  implicit none
  private
  public :: CD_distp1_run, CD_distp2_run
  
  ! La única diferencia entre dist1 y dist2 es que uno lee los datos de termo1
  ! y el otro de termo2, además de los nombres de los archivos de salida

contains
  
  !============================================================================
  ! CD_distp1_run
  !============================================================================
  
  subroutine CD_distp1_run  
    use global
    implicit none
    !*********************************************************************
    !     Este programa calcula, a partir de CD.dat y termo(SHOCK).dat, la 
    !     distribución de electrones en energías para cada punto de la 
    !     curva de un lado de la CD. La salida es Np(SHOCK).dat, en la que
    !     se guarda Np(l,E(i)), con ml y mE datos respectivamente. Los valores 
    !     de las energías posibles los acoto de entrada.
    !*********************************************************************
    integer :: l_min,npaso,l2
    real(dp) :: t_coolp,t_ac,t_adi,t_adv,t_pp,efiac,t_diff,D_diff
    real(dp) :: log10Emaxl,Emaxl,R_acel,R1,R2, E_corr,integral,Q
    real(dp) :: E_punto,Emax_p0,Ka,E_NT,N_totalp_i
    real(dp) :: E_punto0(mE_p),N_p0(mE),nmedio(ml),L_iny(ml),Dl(ml)
    allocate ( x(ml),y(ml),z(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
    allocate ( N_p(mE_p),E_p(mE_p),B_par(ml),N_ptot(ml,mE_p) )
    N_ptot = 0.d0
    E_corr = Estiff_p**(alpha2-alpha)
    
    !////////////////////////////////////////////////////////////////////////
    !     Lectura de datos previos
    
    !     Los puntos (x,y) de la curva, con (x(1),y(1)) las coords del apex,
    !     junto con las variables termodinámicas en esos  puntos (x(l),y(l))
    
    open(11,file='termo1.dat') ! De acá saco los datos
    do l=1,ml
       read(11,*) l2,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
    enddo
    close (11)
    B_par = B*par ! Sólo la componente paralela del B importa para el sy
    nmedio = rho/mp/mu_1

    open(12,file='Liny1.dat')
    do l=1,ml/2
       read(12,*) x(l),y(l),L_iny(l),Dl(l) ! la rama de arriba de la trayectoria en cm
    enddo
    close(12)
    L_iny = iny_p*L_iny ! Uso la misma fracción en cada emisor línea
    
    efiac = efi_ac1


    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                              PROTONES                                 C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías (común para todas las posibles líneas):
    
    call vectorE(Emin_p,Emax_p,Eint_p,E_p)
    
    ! Simulo distintas trayectorias arrancando desde distintos puntos 
    ! (distinto valor de l_min). Si tomo un paso más chico gano más detalle.
    
    
    !////////////////////////////////////////////////////////////////////////
    ! Voy a calcular sólo para la mitad (ml/2) de los puntos, para el resto
    ! serán las mismas distribuciones, tan sólo debo "espejarlas"
    
    open(70,file='CD_tcoolp1.dat')
    
    npaso = 1   !      npaso = (ml/2)/128
    
    Do l_min=1,ml/2,npaso  
       
       !  Calculo para el primer punto de la línea la energía máxima, 
       !  la distribución Ne_0(E_0) y el Epunto0, y ya sumo a Netot(l,i)
       
       l = l_min
       
       R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
       R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
       ! Tomo tamaño del acelerador como ~ la distancia a la estrella más 
       ! cercana, ya que las prop de la region postshock varían poco a escalas del
       ! tamaño del choque, es más allá que las densidades de diluirán considerably
       if (R1 <= R2) then
          R_acel = R1
       else
          R_acel = R2
       endif

       t_adv = Dl(l)/v(l)           ! Tiempo de difusión (uso para la inyección)
       t_adi = R_acel/v(l)          ! Pérdidas adiabáticas
       t_pp = 1.d15/nmedio(l)       ! pérdidas protón-protón
       t_coolp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      
       ! SI NO TENGO EN CUENTA p-gamma, ESTE YA ES EL t_coolp
       
       call Biseccion_p(B(l),log10Emaxl,t_coolp,efiac)
       Emax_p0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))! Pérdidas difusivas
       if (Emax_p0 > Emax_p) then
          write(*,*) 'TE QUEDASTE CORTO CON EMAX_P !!!'
       endif
       
       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculamos Ka normalizando: 
       ! int Q(E')*E'dE' = Linj y despejando. [Lrel]=erg/seg, [Ka]=erg/seg, 
       ! [Q]=erg/erg**2/seg
       
       integral=0.d0
       i = 1
       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0)
          dE_p = E_p(i)*(Eint_p - 1.d0)
               if (E_p(i) <= Estiff_p) then
                  integral = integral + E_p(i)**(-alpha+1.d0)*exp(-E_p(i)/Emax_p0)*dE_p
               elseif (E_p(i) > Estiff_p) then
                  integral = integral + E_p(i)**(-alpha2+1.d0)*exp(-E_p(i)/Emax_p0)*&
                       dE_p*E_corr
               endif
          i=i+1
       enddo
       
       Ka = L_iny(l)/integral ! En cada celda inyecto esa fracción de L_w
       
       N_p0 = 0.d0       ! Inicializo a cero para cada loop en l_min
       
       i=1
       do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0) 

          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, con |dE/dt|= E/t_cool 
          ! en [erg/seg], dt = t_adv, con lo cual [N(E)]=erg/erg**2

          if (E_p(i) <= Estiff_p) then
             Q = E_p(i)**(-alpha)*exp(-E_p(i)/Emax_p0) !normalizo después
          elseif  (E_p(i) > Estiff_p) then
             Q = E_p(i)**(-alpha2)*exp(-E_p(i)/Emax_p0)*E_corr
          endif
          N_p0(i) = Ka*Q*t_adv

          D_diff = E_p(i)*c/(3.d0*qe*B(l))
          t_diff = R_acel**2.d0/(2.d0*D_diff) ! Pérdidas por difusión

          t_ac = efiac*E_p(i)/(B(l)*c*qe)
            
          E_punto0(i) = abs(E_p(i)/t_coolp)
          N_ptot(l_min,i) = N_p0(i)
          write(70,*) l_min,E_p(i)/eV,t_pp,t_adv,t_ac,t_coolp,t_adi,t_diff
          i = i + 1
       enddo
         



!/////////////////////////////////////////////////////////////////////////


!     Ahora calculamos la versión evolucionada en cada punto como:
!     N(E,l) = N0(E0)*(Epunto0/Epunto)*theta(E-Emax) , para lo cual sólo 
!     necesito el valor de Epunto y Emax en "l"
         
         do l=l_min+1,ml/2
            R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
            R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
          
!     Busco la Emax como aquella en que t_cool = t_ac. Para ello busco las 
!     raíces de la función resta mediante una subrutina:
            
            open(78,file='Emax.dat') ! Acá guardo la Emax en cada locación

            t_adi = R_acel/v(l)          ! Pérdidas adiabáticas
            t_pp = 1.d15/nmedio(l)       ! pérdidas protón-protón
            t_coolp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      

            call Biseccion_p(B(l),log10Emaxl,t_coolp,efiac)       
            Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))! Pérdidas difusivas
            write(78,*) (R1-x(1))/UA,Emaxl/(1.d9*eV)
            
            i = 1
            do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emaxl)
!               t_ac = efiac*E_p(i)/(B(l)*c*qe)              
               E_punto = abs(E_p(i)/t_coolp)
               N_p(i) = N_p0(i)*E_punto0(i)/E_punto
               N_ptot(l,i) = N_ptot(l,i)+N_p(i) ! Las energías coinciden
               i = i + 1
            enddo               ! Cierro el do en las energías
            
         enddo                  ! Cierro el do en l
         
      Enddo                     ! Cierro el do en l_min (# trayectorias)

!     Felizmente guardo los datos de la distribución del lado 1
      open(79,file='distp1.dat')
      do l=1,ml/2
         do i=1,mE_p
            write(79,*) l,E_p(i)/eV,N_ptot(l,i)
         enddo
         write(79,*) 
      enddo
      do l=ml/2+1,ml
         do i=1,mE_p
            write(79,*) l,E_p(i)/eV,N_ptot(l-ml/2,i)
         enddo
         write(79,*) 
      enddo
      close(79)
      
      open(80,file='ENTp1.dat')
      do l=1,ml/2
         E_NT = 0.d0
         do i=1,mE_p
            E_NT = E_NT + E_p(i)**2*N_ptot(l,i)
         enddo
         write(80,*) l, E_NT
      enddo
      close(80)

      open(81,file='distp1_total.dat')
      do i=1,mE_p
         N_totalp_i = 0.d0
         do l=1,ml/2
            N_totalp_i = N_totalp_i + N_ptot(l,i)
         enddo
         N_totalp_i = N_totalp_i * 2.d0*mphi ! La total de la WCR
         write(81,*) E_p(i)/eV, N_totalp_i
      enddo
      close(81)

      deallocate ( x,y,z,B,rho,P,v,T )
      deallocate ( N_p,E_p,B_par,N_ptot )
      
    END subroutine CD_distp1_run




    !============================================================================
    ! CD_distp2_run
    !============================================================================
    
    subroutine CD_distp2_run  
      use global
      implicit none
      !*********************************************************************
      !     Este programa calcula, a partir de CD.dat y termo(SHOCK).dat, la 
      !     distribución de electrones en energías para cada punto de la 
      !     curva de un lado de la CD. La salida es Ne(SHOCK).dat, en la que
      !     se guarda Ne(l,E(i)), con ml y mE datos respectivamente. Los valores 
      !     de las energías posibles los acoto de entrada.
      !*********************************************************************
      integer :: l_min,npaso,l2
      real(dp) :: t_coolp,t_ac,t_adi,t_adv,t_pp,efiac,t_diff,D_diff
      real(dp) :: log10Emaxl,Emaxl,R_acel,R1,R2,E_corr,integral,Q
      real(dp) :: E_punto,Emax_p0,Ka,E_NT,N_totalp_i
      real(dp) :: E_punto0(mE_p),N_p0(mE),nmedio(ml),L_iny(ml),Dl(ml)
      allocate ( x(ml),y(ml),z(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
      allocate ( N_p(mE_p),E_p(mE_p),B_par(ml),N_ptot(ml,mE_p) )
      N_ptot = 0.d0
      E_corr = Estiff_p**(alpha2-alpha)
      
      
      !////////////////////////////////////////////////////////////////////////
      !     Lectura de datos previos
      
      !     Los puntos (x,y) de la curva, con (x(1),y(1)) las coords del apex,
      !     junto con las variables termodinámicas en esos  puntos (x(l),y(l))
      
      open(11,file='termo2.dat') ! De acá saco los datos
      do l=1,ml
         read(11,*) l2,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
      enddo
      close (11)
      B_par = B*par ! Sólo la componente paralela del B importa para el sy
      nmedio = rho/mp/mu_2
      
      open(12,file='Liny2.dat')
      do l=1,ml/2
         read(12,*) x(l),y(l),L_iny(l),Dl(l) ! la rama de arriba de la trayectoria en cm
      enddo
      close(12)
      L_iny = iny_p*L_iny ! Uso la misma fracción en cada emisor línea
      
      efiac = efi_ac2
      

      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !                                                                       C
      !                              PROTONES                                 C
      !                                                                       C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      ! Armo el vector de energías (común para todas las posibles líneas):
      
      call vectorE(Emin_p,Emax_p,Eint_p,E_p)
      
      ! Simulo distintas trayectorias arrancando desde distintos puntos 
      ! (distinto valor de l_min). Si tomo un paso más chico gano más detalle.
      
      
      !////////////////////////////////////////////////////////////////////////
      ! Voy a calcular sólo para la mitad (ml/2) de los puntos, para el resto
      ! serán las mismas distribuciones, tan sólo debo "espejarlas"
      
      open(70,file='CD_tcoolp2.dat')
      
      npaso = 1   !      npaso = (ml/2)/128
      
      
      Do l_min=1,ml/2,npaso  
         
         !  Calculo para el primer punto de la línea la energía máxima, 
         !  la distribución Ne_0(E_0) y el Epunto0, y ya sumo a Netot(l,i)
         
         l = l_min
         
         R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
         R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
         ! Tomo tamaño del acelerador como ~ la distancia a la estrella más 
         ! cercana, ya que las prop de la region postshock varían poco a escalas del
         ! tamaño del choque, es más allá que las densidades de diluirán considerably
         if (R1 <= R2) then
            R_acel = R1
         else
            R_acel = R2
         endif
         
         t_adv = Dl(l)/v(l)           ! Tiempo de difusión (uso para la inyección)
         t_adi = R_acel/v(l)          ! Pérdidas adiabáticas
         t_pp = 1.d15/nmedio(l)       ! pérdidas protón-protón
         t_coolp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      
         ! SI NO TENGO EN CUENTA p-gamma, ESTE YA ES EL t_coolp
         
         call Biseccion_p(B(l),log10Emaxl,t_coolp,efiac)
         Emax_p0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))! Pérdidas difusivas
         if (Emax_p0 > Emax_p) then
            write(*,*) 'TE QUEDASTE CORTO CON EMAX_P !!!'
         endif
         
         ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculamos Ka normalizando: 
         ! int Q(E')*E'dE' = Linj y despejando. [Lrel]=erg/seg, [Ka]=erg/seg, 
         ! [Q]=erg/erg**2/seg
         
         integral=0.d0
         i = 1
         do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0)
            dE_p = E_p(i)*(Eint_p - 1.d0)
            if (E_p(i) <= Estiff_p) then
               integral = integral + E_p(i)**(-alpha+1.d0)*exp(-E_p(i)/Emax_p0)*dE_p
            elseif (E_p(i) > Estiff_p) then
               integral = integral + E_p(i)**(-alpha2+1.d0)*exp(-E_p(i)/Emax_p0)*&
                    dE_p*E_corr
            endif
            i=i+1
         enddo
         
         Ka = L_iny(l)/integral ! En cada celda inyecto esa fracción de L_w
         
         N_p0 = 0.d0       ! Inicializo a cero para cada loop en l_min
         
         i=1
         do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emax_p0) 
            
            ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, con |dE/dt|= E/t_cool 
            ! en [erg/seg], dt = t_adv, con lo cual [N(E)]=erg/erg**2
            
            if (E_p(i) <= Estiff_p) then
               Q = E_p(i)**(-alpha)*exp(-E_p(i)/Emax_p0) !normalizo después
            elseif  (E_p(i) > Estiff_p) then
               Q = E_p(i)**(-alpha2)*exp(-E_p(i)/Emax_p0)*E_corr
            endif
            N_p0(i) = Ka*Q*t_adv
            
            D_diff = E_p(i)*c/(3.d0*qe*B(l))
            t_diff = R_acel**2.d0/(2.d0*D_diff) ! Pérdidas por difusión
            
            t_ac = efiac*E_p(i)/(B(l)*c*qe)
            
            E_punto0(i) = abs(E_p(i)/t_coolp)
            N_ptot(l_min,i) = N_p0(i)
            write(70,*) l_min,E_p(i)/eV,t_pp,t_adv,t_ac,t_coolp,t_adi,t_diff
            i = i + 1
         enddo
         
         
         
         !/////////////////////////////////////////////////////////////////////////
         
         
         !     Ahora calculamos la versión evolucionada en cada punto como:
         !     N(E,l) = N0(E0)*(Epunto0/Epunto)*theta(E-Emax) , para lo cual sólo 
         !     necesito el valor de Epunto y Emax en "l"
         
         do l=l_min+1,ml/2
            R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
            R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
            
            ! Busco la Emax como aquella en que t_cool = t_ac. Para ello busco las 
            ! raíces de la función resta mediante una subrutina:
            
            open(78,file='Emax.dat') ! Acá guardo la Emax en cada locación
            
            t_adi = R_acel/v(l)          ! Pérdidas adiabáticas
            t_pp = 1.d15/nmedio(l)       ! pérdidas protón-protón
            t_coolp = 1.d0/(1.d0/t_adi + 1.d0/t_pp)      
            
            call Biseccion_p(B(l),log10Emaxl,t_coolp,efiac)       
            Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))! Pérdidas difusivas
            write(78,*) (R1-x(1))/UA,Emaxl/(1.d9*eV)
            
            i = 1
            do while (i < mE_p+1 .and. E_p(i) < 10.d0*Emaxl)
!               t_ac = efiac*E_p(i)/(B(l)*c*qe)              
               E_punto = abs(E_p(i)/t_coolp)
               N_p(i) = N_p0(i)*E_punto0(i)/E_punto
               N_ptot(l,i) = N_ptot(l,i)+N_p(i) ! Las energías coinciden
               i = i + 1
            enddo               ! Cierro el do en las energías
            
         enddo                  ! Cierro el do en l
         
      Enddo                     ! Cierro el do en l_min (# trayectorias)
      
      
      ! Felizmente guardo los datos de la distribución del lado 2
      open(79,file='distp2.dat')
      do l=1,ml/2
         do i=1,mE_p
            write(79,*) l,E_p(i)/eV,N_ptot(l,i)
         enddo
         write(79,*) 
      enddo
      do l=ml/2+1,ml
         do i=1,mE_p
            write(79,*) l,E_p(i)/eV,N_ptot(l-ml/2,i)
         enddo
         write(79,*) 
      enddo
      close(79)
      
      open(80,file='ENTp2.dat')
      do l=1,ml/2
         E_NT = 0.d0
         do i=1,mE_p
            E_NT = E_NT + E_p(i)**2*N_ptot(l,i)
         enddo
         write(80,*) l, E_NT
      enddo
      close(80)

      open(81,file='distp2_total.dat')
      do i=1,mE_p
         N_totalp_i = 0.d0
         do l=1,ml/2
            N_totalp_i = N_totalp_i + N_ptot(l,i)
         enddo
         N_totalp_i = N_totalp_i * 2.d0*mphi ! La total de la WCR
         write(81,*) E_p(i)/eV, N_totalp_i
      enddo
      close(81)

      deallocate ( x,y,z,B,rho,P,v,T )
      deallocate ( N_p,E_p,B_par,N_ptot )

      
      
    END subroutine CD_distp2_run
    
    !==========================================================================
    
    
    
  END module CD_distp

  
  
  !///////////////////////////////////////////////////////////////////////////
  !                                                                         //
  !                            SUBRUTINAS                                   //
  !                                                                         //
  !///////////////////////////////////////////////////////////////////////////
  
  Subroutine Biseccion_p(Be,x3,tcool,xi)
    use global
    implicit none
    integer :: t_err
    real(dp) :: resta_p,a_,b_,c_,err
    real(dp) :: Be,x3,tcool,xi
    external resta_p
    
    !     El error en el n-esimo será menor a |a-b|/2**n
    a_ = log10(1.d11*eV) ! log10(Emin [CGS]) ~ 100 GeV
    b_ = log10(1.d15*eV) ! log10(Emax [CGS]) ~ 1 PeV
    t_err = 3
    err = 10.d0**(-t_err) ! error menor que 10**-t
    i = 0
    do while ( abs(b_ - a_) > err )
       i = i+1
       c_ = (a_+b_)/2.
       if (resta_p(c_,Be,tcool,xi)==0) then
          b_ = c_
       else if (resta_p(a_,Be,tcool,xi)*resta_p(c_,Be,tcool,xi) < 0) then
          b_ = c_
       else if (resta_p(a_,Be,tcool,xi)*resta_p(c_,Be,tcool,xi) > 0) then
          a_ = c_
       endif
    enddo
    x3=b_
    return
  end Subroutine Biseccion_p
  
  !/////////////////////////////////////////////////////////////////////////
  
  Function resta_p(log10E,B_l,tcool,xi)
    use global
    implicit none
    real(dp) :: resta_p,log10E,B_l,E_l
    real(dp) :: tac,tcool,xi
    E_l = 10.d0**log10E
    tac = xi*E_l/(B_l*c*qe)
    resta_p = tcool - tac
    
    return
    
  end Function resta_p
  
  !/////////////////////////////////////////////////////////////////////////

