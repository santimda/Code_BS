module CD_dist
  
  use global
  
  implicit none
  private
  public :: CD_dist1_run, CD_dist2_run
  
  ! La única diferencia entre dist1 y dist2 es que uno lee los datos de termo1
  ! y el otro de termo2, además de los nombres de los archivos de salida

contains
  
  !============================================================================
  ! CD_dist1_run
  !============================================================================
  
  subroutine CD_dist1_run  
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
    real(dp) :: D_diff,tdiff,integral
    real(dp) :: tcoole,tic1,tic2,tic,tsy,tac,tadi,tadv,tbr,efiac
    real(dp) :: gamm,E_corr,C_ic1_k,C_ic1,C_ic2_k,C_ic2
    real(dp) :: log10Emaxl,Emaxl,R_acel,R1,R2,Q
    real(dp) :: E_punto,Emax_e0,Ka,E_NT,N_total_i
    real(dp) :: E_punto0(mE),N_e0(mE),nmedio(ml),L_iny(ml),Dl(ml)
    allocate ( x(ml),y(ml),z(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
    allocate ( N_e(mE),E_e(mE),B_par(ml),N_etot(ml,mE) )
    N_etot = 0.d0
    E_corr = Estiff_e**(alpha2-alpha)

    !////////////////////////////////////////////////////////////////////////
    !     Lectura de datos previos
    
    !     Los puntos (x,y) de la curva, con (x(1),y(1)) las coords del apex,
    !     junto con las variables termodinámicas en esos puntos (x(l),y(l))
    
    open(11,file='termo1.dat') 
    do l=1,ml
       read(11,*) l2,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
    enddo
    close (11)
    B_par = B*par        ! Sólo la componente paralela del B importa para el sy
    nmedio = rho/mp/mu_1

    open(12,file='Liny1.dat')
    do l=1,ml/2
       read(12,*) x(l),y(l),L_iny(l),Dl(l) ! la rama de arriba de la trayectoria en cm
    enddo
    close(12)
    L_iny = iny_e*L_iny  ! En L_iny ya está contemplado el tamaño de cada celda   
    efiac = efi_ac1
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                              ELECTRONES                               C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías (común para todas las posibles líneas):
    call vectorE(Emin_e,Emax_e,Eint_e,E_e)
    
    ! Simulo distintas trayectorias arrancando desde distintos puntos 
    ! (distinto valor de l_min). Calcular sólo para la mitad (ml/2) de los puntos, 
    ! para el resto son las mismas distribuciones, tan sólo debo "espejarlas"
    
    open(70,file='CD_tcoole1.dat')
    
    npaso = 1 !      npaso = (ml/2)/128 ! si quiero menos emisores

    C_ic1 = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T1mcc**2.)
    C_ic2 = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T2mcc**2.)
    
    Do l_min=1,ml/2,npaso  
       
       ! Calculo para el primer punto de la línea la energía máxima, 
       ! la distribución Ne_0(E_0) y el Epunto0, y ya sumo a Netot(l,i)
       
       l = l_min
       
       R1 = sqrt(x(l)**2+y(l)**2)     ! Dist. del punto a la estrella 1
       R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
       Kappa1 = (Rst1/(2.*R1))**2.
       Kappa2 = (Rst2/(2.*R2))**2.
       C_ic1_k = C_ic1/Kappa1
       C_ic2_k = C_ic2/Kappa2       

       ! Tomo el tamaño del acelerador como ~ la distancia a la estrella más cercana
       if (R1 <= R2) then
          R_acel = R1
       else
          R_acel = R2
       endif
       
       tadv = Dl(l)/v(l)           ! Tiempo de difusión (uso para la inyección)
       tadi = R_acel/v(l)          ! Pérdidas adiabáticas
       tbr = 1.d15/nmedio(l)       ! No depende de la energía
       tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)
       call Biseccion(B(l),log10Emaxl,C_ic1_k,C_ic2_k,tcoole,efiac)
       Emax_e0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))  ! Pérdidas difusivas
       if (Emax_e0 > Emax_e) then
          write(*,*) 'TE QUEDASTE CORTO CON EMAX !!!'
       endif
       
       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculamos Ka normalizando: 
       ! int Q(E')*E'dE' = Linj y despejando. [Lrel]=erg/seg, [Ka]=erg/seg, 
       ! [Q]=erg/erg**2/seg
       
       integral=0.d0
       i = 1
       do while (i < mE+1 .and. E_e(i) < 10.d0*Emax_e0)
          dE_e = E_e(i)*(Eint_e-1.d0)
          if (E_e(i) <= Estiff_e) then
             integral = integral + E_e(i)**(-alpha+1.d0)*exp(-E_e(i)/Emax_e0)*dE_e
          elseif  (E_e(i) > Estiff_e) then
             integral = integral + E_e(i)**(-alpha2+1.d0)*exp(-E_e(i)/Emax_e0)*&
                  dE_e*E_corr
          endif
          i=i+1
       enddo
       
       Ka = L_iny(l)/integral ! En cada celda inyecto esa fracción de L_w
       
       N_e0 = 0.d0            ! Inicializo a cero para cada loop en l_min
         
       i=1
       do while (i < mE+1 .and. E_e(i) < 10.d0*Emax_e0) 
          
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, con |dE/dt|= E/t_cool 
          ! en [erg/seg], dt = t_adv, con lo cual [N(E)]=erg/erg**2
            
          if (E_e(i) <= Estiff_e) then
             Q = E_e(i)**(-alpha)*exp(-E_e(i)/Emax_e0) !normalizo después
          elseif  (E_e(i) > Estiff_e) then
             Q = E_e(i)**(-alpha2)*exp(-E_e(i)/Emax_e0)*E_corr
          endif
          N_e0(i) = Ka*Q*tadv
            
          D_diff = E_e(i)*c/(3.d0*qe*B(l))
          tdiff = R_acel**2.d0/(2.d0*D_diff) ! Pérdidas por difusión
          tsy = 1.d0/1.6d-3/B_par(l)**2.d0/E_e(i) 

          gamm = E_e(i)/mec2
          tic1 = C_ic1_k*gamm/F_iso(4.*gamm*T1mcc)
          tic2 = C_ic2_k*gamm/F_iso(4.*gamm*T2mcc)
          tic = (1.d0/tic1 + 1.d0/tic2)**(-1.)
          tcoole = (1.d0/tsy + 1.d0/tic + 1/tbr + 1.d0/tadi)**(-1.) 
          tac = efiac*E_e(i)/(B(l)*c*qe)
          
          E_punto0(i) = abs(E_e(i)/tcoole)
          N_etot(l_min,i) = N_e0(i)
          write(70,*) l_min,E_e(i)/eV,tic,tsy,tbr,tdiff,tac,tcoole,tadi,tadv
          i = i + 1

       enddo         ! Cierro el do en E_e(i)
         

!/////////////////////////////////////////////////////////////////////////


!     Ahora calculamos la versión evolucionada en cada punto como:
!     N(E,l) = N0(E0)*(Epunto0/Epunto)*theta(E-Emax) , para lo cual sólo 
!     necesito el valor de Epunto y Emax en "l"
         
       do l=l_min+1,ml/2
          R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
          R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 

          Kappa1 = (Rst1/(2.*R1))**2.
          Kappa2 = (Rst2/(2.*R2))**2.
          C_ic1_k = C_ic1/Kappa1
          C_ic2_k = C_ic2/Kappa2     

          ! Busco la Emax como aquella en que t_cool = t_ac. Para ello busco las 
          ! raíces de la función resta mediante una subrutina:
          
          tadi = R_acel/v(l)          ! Pérdidas adiabáticas
          tbr = 1.d15/nmedio(l)       ! No depende de la energía
          tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)

          call Biseccion(B(l),log10Emaxl,C_ic1_k,C_ic2_k,tcoole,efiac)
          Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5/efiac))  ! Pérdidas difusivas
          
          i = 1
          do while (i < mE+1 .and. E_e(i) < 10.d0*Emaxl)
             dE_e = E_e(i)*(Eint_e - 1.d0)  
             tsy = 1.d0/1.6d-3/B_par(l)**2.d0/E_e(i) 
             gamm = E_e(i)/mec2
             tic1 = C_ic1_k*gamm/F_iso(4.*gamm*T1mcc)
             tic2 = C_ic2_k*gamm/F_iso(4.*gamm*T2mcc)
             tic = (1.d0/tic1 + 1.d0/tic2)**(-1.)
             tcoole = (1.d0/tsy + 1.d0/tic + 1.d0/tbr + 1.d0/tadi)**(-1.)
!             tac = efiac*E_e(i)/(B(l)*c*qe)
               
             E_punto = abs(E_e(i)/tcoole)
             N_e(i) = N_e0(i)*E_punto0(i)/E_punto
             
             N_etot(l,i) = N_etot(l,i)+N_e(i) ! Las energías coinciden
             i = i + 1
          enddo               ! Cierro el do en las energías
          
       enddo                  ! Cierro el do en l
       
    enddo                     ! Cierro el do en l_min (# trayectorias)
    close(70)

    !     Felizmente guardo los datos de la distribución del lado 1
    open(79,file='dist1.dat')
    do l=1,ml/2
       do i=1,mE
          write(79,*) l,E_e(i)/eV,N_etot(l,i)
       enddo
       write(79,*) 
    enddo
    do l=ml/2+1,ml
       do i=1,mE
          write(79,*) l,E_e(i)/eV,N_etot(l-ml/2,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    open(80,file='ENT1.dat')
    do l=1,ml/2
       E_NT = 0.d0
       do i=1,mE
          E_NT = E_NT + E_e(i)**2*N_etot(l,i)
       enddo
       E_NT = E_NT * 2.d0*mphi ! La total de la WCR
       write(80,*) l, E_NT
    enddo
    close(80)
    
    open(81,file='dist1_total.dat')
    do i=1,mE
       N_total_i = 0.d0
       do l=1,ml/2
          N_total_i = N_total_i + N_etot(l,i)
       enddo
       N_total_i = N_total_i * 2.d0*mphi ! La total de la WCR
       write(81,*) E_e(i)/eV, N_total_i
    enddo
    close(81)
    
    deallocate ( x,y,z,B,B_par,rho,P,v,T )
    deallocate ( N_e,E_e,N_etot ) 
    
  END subroutine CD_dist1_run
  
  
  !==========================================================================

  

  !============================================================================
  ! CD_dist2_run
  !============================================================================
  
  subroutine CD_dist2_run  
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
    real(dp) :: D_diff,tdiff,integral
    real(dp) :: tcoole,tic1,tic2,tic,tsy,tac,tadi,tadv,tbr,efiac
    real(dp) :: gamm,E_corr,C_ic1_k,C_ic1,C_ic2_k,C_ic2
    real(dp) :: log10Emaxl,Emaxl,R_acel,R1,R2,Q
    real(dp) :: E_punto,Emax_e0,Ka,E_NT,N_total_i
    real(dp) :: E_punto0(mE),N_e0(mE),nmedio(ml),L_iny(ml),Dl(ml)
    allocate ( x(ml),y(ml),z(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
    allocate ( N_e(mE),E_e(mE),B_par(ml),N_etot(ml,mE) )
    N_etot = 0.d0
    E_corr = Estiff_e**(alpha2-alpha)

    !////////////////////////////////////////////////////////////////////////
    !     Lectura de datos previos
    
    !     Los puntos (x,y) de la curva, con (x(1),y(1)) las coords del apex,
    !     junto con las variables termodinámicas en esos puntos (x(l),y(l))
    
    open(11,file='termo2.dat') 
    do l=1,ml
       read(11,*) l2,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
    enddo
    close (11)
    B_par = B*par        ! Sólo la componente paralela del B importa para el sy
    nmedio = rho/mp/mu_1

    open(12,file='Liny2.dat')
    do l=1,ml/2
       read(12,*) x(l),y(l),L_iny(l),Dl(l) ! la rama de arriba de la trayectoria en cm
    enddo
    close(12)
    L_iny = iny_e*L_iny  ! En L_iny ya está contemplado el tamaño de cada celda   
    efiac = efi_ac1
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                              ELECTRONES                               C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías (común para todas las posibles líneas):
    call vectorE(Emin_e,Emax_e,Eint_e,E_e)
    
    ! Simulo distintas trayectorias arrancando desde distintos puntos 
    ! (distinto valor de l_min). Calcular sólo para la mitad (ml/2) de los puntos, 
    ! para el resto son las mismas distribuciones, tan sólo debo "espejarlas"
    
    open(70,file='CD_tcoole2.dat')
    
    npaso = 1 !      npaso = (ml/2)/128 ! si quiero menos emisores

    C_ic1 = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T1mcc**2.)
    C_ic2 = ( pi*h_bar**3.*c**2./(2.*re**2.*mec2**3.) )/(T2mcc**2.)
    
    Do l_min=1,ml/2,npaso  
       
       ! Calculo para el primer punto de la línea la energía máxima, 
       ! la distribución Ne_0(E_0) y el Epunto0, y ya sumo a Netot(l,i)
       
       l = l_min
       
       R1 = sqrt(x(l)**2+y(l)**2)     ! Dist. del punto a la estrella 1
       R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 2
       Kappa1 = (Rst1/(2.*R1))**2.
       Kappa2 = (Rst2/(2.*R2))**2.
       C_ic1_k = C_ic1/Kappa1
       C_ic2_k = C_ic2/Kappa2       

       ! Tomo el tamaño del acelerador como ~ la distancia a la estrella más cercana
       if (R1 <= R2) then
          R_acel = R1
       else
          R_acel = R2
       endif
       
       tadv = Dl(l)/v(l)           ! Tiempo de difusión (uso para la inyección)
       tadi = R_acel/v(l)          ! Pérdidas adiabáticas
       tbr = 1.d15/nmedio(l)       ! No depende de la energía
       tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)
       call Biseccion(B(l),log10Emaxl,C_ic1_k,C_ic2_k,tcoole,efiac)
       Emax_e0 = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5d0/efiac))  ! Pérdidas difusivas
       if (Emax_e0 > Emax_e) then
          write(*,*) 'TE QUEDASTE CORTO CON EMAX !!!'
       endif
       
       ! Q = Ka*E**-alpha*e^(-E/Emax), alpha=2. Calculamos Ka normalizando: 
       ! int Q(E')*E'dE' = Linj y despejando. [Lrel]=erg/seg, [Ka]=erg/seg, 
       ! [Q]=erg/erg**2/seg
       
       integral=0.d0
       i = 1
       do while (i < mE+1 .and. E_e(i) < 10.d0*Emax_e0)
          dE_e = E_e(i)*(Eint_e-1.d0)
          if (E_e(i) <= Estiff_e) then
             integral = integral + E_e(i)**(-alpha+1.d0)*exp(-E_e(i)/Emax_e0)*dE_e
          elseif  (E_e(i) > Estiff_e) then
             integral = integral + E_e(i)**(-alpha2+1.d0)*exp(-E_e(i)/Emax_e0)*&
                  dE_e*E_corr
          endif
          i=i+1
       enddo
       
       Ka = L_iny(l)/integral ! En cada celda inyecto esa fracción de L_w
       
       N_e0 = 0.d0            ! Inicializo a cero para cada loop en l_min
         
       i=1
       do while (i < mE+1 .and. E_e(i) < 10.d0*Emax_e0) 
          
          ! N(E)=1/|dE/dt|*int_E^Emax Q(E')dE' = Q(E)*dt, con |dE/dt|= E/t_cool 
          ! en [erg/seg], dt = t_adv, con lo cual [N(E)]=erg/erg**2
            
          if (E_e(i) <= Estiff_e) then
             Q = E_e(i)**(-alpha)*exp(-E_e(i)/Emax_e0) !normalizo después
          elseif  (E_e(i) > Estiff_e) then
             Q = E_e(i)**(-alpha2)*exp(-E_e(i)/Emax_e0)*E_corr
          endif
          N_e0(i) = Ka*Q*tadv
            
          D_diff = E_e(i)*c/(3.d0*qe*B(l))
          tdiff = R_acel**2.d0/(2.d0*D_diff) ! Pérdidas por difusión
          tsy = 1.d0/1.6d-3/B_par(l)**2.d0/E_e(i) 

          gamm = E_e(i)/mec2
          tic1 = C_ic1_k*gamm/F_iso(4.*gamm*T1mcc)
          tic2 = C_ic2_k*gamm/F_iso(4.*gamm*T2mcc)
          tic = (1.d0/tic1 + 1.d0/tic2)**(-1.)
          tcoole = (1.d0/tsy + 1.d0/tic + 1/tbr + 1.d0/tadi)**(-1.) 
          tac = efiac*E_e(i)/(B(l)*c*qe)
          
          E_punto0(i) = abs(E_e(i)/tcoole)
          N_etot(l_min,i) = N_e0(i)
          write(70,*) l_min,E_e(i)/eV,tic,tsy,tbr,tdiff,tac,tcoole,tadi,tadv
          i = i + 1

       enddo         ! Cierro el do en E_e(i)
         

!/////////////////////////////////////////////////////////////////////////


!     Ahora calculamos la versión evolucionada en cada punto como:
!     N(E,l) = N0(E0)*(Epunto0/Epunto)*theta(E-Emax) , para lo cual sólo 
!     necesito el valor de Epunto y Emax en "l"
         
       do l=l_min+1,ml/2
          R1 = sqrt(x(l)**2+y(l)**2) ! Dist. del punto a la estrella 1
          R2 = sqrt((D-x(l))**2+y(l)**2) ! Dist. del punto a la estrella 

          Kappa1 = (Rst1/(2.*R1))**2.
          Kappa2 = (Rst2/(2.*R2))**2.
          C_ic1_k = C_ic1/Kappa1
          C_ic2_k = C_ic2/Kappa2     

          ! Busco la Emax como aquella en que t_cool = t_ac. Para ello busco las 
          ! raíces de la función resta mediante una subrutina:
          
          tadi = R_acel/v(l)          ! Pérdidas adiabáticas
          tbr = 1.d15/nmedio(l)       ! No depende de la energía
          tcoole = 1.d0/(1.d0/tadi + 1.d0/tbr)

          call Biseccion(B(l),log10Emaxl,C_ic1_k,C_ic2_k,tcoole,efiac)
          Emaxl = MIN(10.d0**log10Emaxl,qe*B(l)*R_acel*sqrt(1.5/efiac))  ! Pérdidas difusivas
          
          i = 1
          do while (i < mE+1 .and. E_e(i) < 10.d0*Emaxl)
             dE_e = E_e(i)*(Eint_e - 1.d0)  
             tsy = 1.d0/1.6d-3/B_par(l)**2.d0/E_e(i) 
             gamm = E_e(i)/mec2
             tic1 = C_ic1_k*gamm/F_iso(4.*gamm*T1mcc)
             tic2 = C_ic2_k*gamm/F_iso(4.*gamm*T2mcc)
             tic = (1.d0/tic1 + 1.d0/tic2)**(-1.)
             tcoole = (1.d0/tsy + 1.d0/tic + 1.d0/tbr + 1.d0/tadi)**(-1.)
!             tac = efiac*E_e(i)/(B(l)*c*qe)
               
             E_punto = abs(E_e(i)/tcoole)
             N_e(i) = N_e0(i)*E_punto0(i)/E_punto
             
             N_etot(l,i) = N_etot(l,i)+N_e(i) ! Las energías coinciden
             i = i + 1
          enddo               ! Cierro el do en las energías
          
       enddo                  ! Cierro el do en l
       
    enddo                     ! Cierro el do en l_min (# trayectorias)
    close(70)

    !     Felizmente guardo los datos de la distribución del lado 1
    open(79,file='dist2.dat')
    do l=1,ml/2
       do i=1,mE
          write(79,*) l,E_e(i)/eV,N_etot(l,i)
       enddo
       write(79,*) 
    enddo
    do l=ml/2+1,ml
       do i=1,mE
          write(79,*) l,E_e(i)/eV,N_etot(l-ml/2,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    open(80,file='ENT2.dat')
    do l=1,ml/2
       E_NT = 0.d0
       do i=1,mE
          E_NT = E_NT + E_e(i)**2*N_etot(l,i)
       enddo
       E_NT = E_NT * 2.d0*mphi ! La total de la WCR
       write(80,*) l, E_NT
    enddo
    close(80)
    
    open(81,file='dist2_total.dat')
    do i=1,mE
       N_total_i = 0.d0
       do l=1,ml/2
          N_total_i = N_total_i + N_etot(l,i)
       enddo
       N_total_i = N_total_i * 2.d0*mphi ! La total de la WCR
       write(81,*) E_e(i)/eV, N_total_i
    enddo
    close(81)
    
    deallocate ( x,y,z,B,B_par,rho,P,v,T )
    deallocate ( N_e,E_e,N_etot ) 
    
    
  END subroutine CD_dist2_run

  !==========================================================================
  
  
  
END module CD_dist



!///////////////////////////////////////////////////////////////////////////
!                                                                         //
!                            SUBRUTINAS                                   //
!                                                                         //
!///////////////////////////////////////////////////////////////////////////

Subroutine Biseccion(Be,x3,c1,c2,tad,xi)
  use global
  implicit none
  integer t_err
  real(dp) :: resta,a_,b_,c_,err
  real(dp), intent(in) :: Be,c1,c2,tad,xi
  real(dp), intent(out) :: x3
  external resta
  
  !     El error en el n-esimo será menor a |a-b|/2**n
  a_ = log10(1.d10*eV) ! log10(Emin [CGS]) ~ 10 GeV
  b_ = log10(5.d12*eV) ! log10(Emax [CGS]) ~ 5 TeV
  t_err = 2
  err = 10.d0**(-t_err) ! error menor que 10**-t
  i = 0
  do while ( abs(b_ - a_) > err )
     i = i+1
     c_ = (a_+b_)/2.
     if (resta(c_,Be,c1,c2,tad,xi)==0) then
        b_ = c_
     else if (resta(a_,Be,c1,c2,tad,xi)*resta(c_,Be,c1,c2,tad,xi) < 0) then
        b_ = c_
     else if (resta(a_,Be,c1,c2,tad,xi)*resta(c_,Be,c1,c2,tad,xi) > 0) then
        a_ = c_
     endif
  enddo
  x3=b_
  return
end Subroutine Biseccion

!/////////////////////////////////////////////////////////////////////////

Function resta(log10E,B_l,c1,c2,tad,xi)
  use global
  implicit none
  real(dp), intent(in) :: log10E,B_l,c1,c2,tad,xi
  real(dp) :: resta,E_l,gamm,tsy,tic1,tic2,tac,tcoole
  E_l = 10.d0**log10E
  gamm = E_l/mec2
  tsy = 1.d0/1.6d-3/(B_l*par)**2.d0/E_l   
  tic1 = c1*gamm/F_iso(4.*gamm*T1mcc)
  tic2 = c2*gamm/F_iso(4.*gamm*T2mcc)
  
  tcoole = (1.d0/tsy + 1.d0/tic1 + 1.d0/tic2 + 1.d0/tad)**(-1.) ! en tad ya está br
  tac = xi*E_l/(B_l*c*qe)  !      tac = (efiac*B_l*c*qe/E_l)**(-1)
  resta = tcoole - tac
  
  return
  
end Function resta

!/////////////////////////////////////////////////////////////////////////

