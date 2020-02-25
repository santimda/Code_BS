module CD_radsyn
  !Este módulo calcula la emisión sincrotrón NO CORREGIDA POR ABSORCIÓN
  use global
  
  implicit none
  private
  public :: CD_radsyn1_run, CD_radsyn2_run
  
  ! La única diferencia entre radsyn1 y radsyn2 es que uno lee los datos de 
  ! dist1 y el otro de dist2, además de los nombres de los archivos de salida.
  ! En radsyn.dat queda la emisión integrada en l, mientras que en emisyn.dat
  ! queda la emisión Eps x L_sy(l,Eps_i), que luego corregiré por abs.
  ! radsyn2 aprovecha que radsyn1 arma el vector nu(i) y usa el mismo, así 
  ! me queda calculada la emisión para exactamente las mismas frecuencias

contains
  
  !============================================================================
  ! CD_radsyn1_run
  !============================================================================
  
  subroutine CD_radsyn1_run  
    use global
    !$ use omp_lib
    implicit none
    ! ml - CD, mE - electrones, mnu - fotones emitidos/tabla abs_ff
    integer :: il
    real(dp) :: Eps(mnu),dEps,Eps_int, gamma(mE)
    real(dp) :: factor,nuc,nuc_max,nuc_min,sum,sum_sy
    real(dp) :: ni,ne,nu_RT
    allocate ( x(ml),y(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
    allocate ( N_e(mE),E_e(mE),B_par(ml),N_etot(ml,mE) )
    allocate ( Lsy_tot(mnu),nu(mnu),Lsy(ml,mnu) )
    Lsy_tot = 0.d0

    ! set omp threads number
    !$ call omp_set_num_threads(nth) 

    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y) de la curva y las variables 
    ! termodinámicas en esos puntos (x(l),y(l))
    
    open(11,file='termo1.dat')
    do l=1,ml
       read(11,*) il,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
    enddo
    close(11)
    B_par = B*par

    open(12,file='CD.dat')
    do l=1,ml
       read(12,*) x(l),y(l) ! la trayectoria en cm
    enddo
    close(12)

    ! Distribución de e- en el punto r(l) con energía E_e(i): N_etot(l,i) 
    open (14,file='dist1.dat')
    do l=1,ml
       do i=1,mE
          read(14,*) j, E_e(i),N_etot(l,i)
          E_e(i) = E_e(i)*eV    !   Pasamos la energía a ergios
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma = E_e/mec2
    
    !------------------------------------------------------------------------
    !
    !  Comienzo los cálculos    
    !

    factor = h/mec2
    ! Uso el mismo rango de frecuencias para TODAS las posiciones, para eso
    ! tomo la más máxima y la más mínima de todas las frecuencias, con un margen
    nuc_max = 4.22d6*B_par(1)*gamma(mE)**2 / 1.d2 ! No llega tan lejos
    nuc_min = 4.22d6*B_par(ml)*gamma(1)**2 / 1.d2
    
    ! Armo el vector de energías de fotones, Eps(i):
    call vectorE(factor*nuc_min,10.d0*factor*nuc_max,Eps_int,Eps)
    nu = Eps/factor ! [nu] = Hz ?

    open(66,file='nu.dat')
    do i=1,mnu
       write(66,*) nu(i)
    enddo
    close(66)

    ! Para cada posición voy a calcular la emisión sinc y luego sumo todas
    sum_sy=0.d0
    DO l=1,ml/2
       do i=1,mE
          N_e(i)=N_etot(l,i)
       enddo
       ni = rho(l)/mp/2.085d0 ! 2.085 = x + 4*y + 14*z
       ne = 1.395d0*ni         ! 1.395 = x + 2*y + 7*z
       nu_RT = 20.d0*ne/B(l)  ! frecuencia característica para efecto R-T

       ! Barremos todos los rangos de energía de los fotones, Epsilon (Eps), 
       ! entre Eps_min y Eps_max, dados por nu_c(gamma_min) y nu_c(gamma_max)
       do i = 1,mnu           ! Eps(i), nu(i)
          dEps = Eps(i)*(Eps_int-1.d0)
          sum = 0.d0
          ! Usamos E como variable de integración, P_eps = P(E,epsilon),
          ! Lsy = P(epsilon) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps],
          ! x = nu/nu_c va entre 0 y 10

          !$omp parallel do private(nuc,dP_Eps,P_Eps,dE_e) &
          !$omp reduction(+ : sum)
          do j = 1,mE         ! Ne(j), nuc(j), gamma(j)
             nuc = 4.22d6*B_par(l)*gamma(j)**2
             if (nu(i)/nuc > 10.d0) then
                dP_Eps = 0.d0
             else
                ! Calculamos P_eps [erg/s/Hz] y pasamos a [1/s] dividiendo por h
                P_Eps = (4.39d-22*B_par(l)*(nu(i)/nuc)**(1.d0/3.d0)*exp(-nu(i)/nuc))/h
                dE_e = E_e(j)*(Eint_e - 1.d0)
                dP_Eps = N_e(j)*P_eps*dE_e
             endif
             sum = dP_Eps + sum
          enddo    
          !$omp end parallel do 
    
          Lsy(l,i) = sum*mec2*exp(-nu_RT/nu(i)) ! Lsy en esa posición y a esa 
          !energía, corregida por efecto R-T
          
          Lsy_tot(i) = Lsy_tot(i) + Eps(i)*Lsy(l,i) ! [erg/s]
          
          sum_sy = sum_sy + Lsy(l,i)*dEps

       enddo  ! Eps(i),nu(i)
       
    ENDDO     ! l
    
    write(*,*) 'L_syn(unabsorbed)=',sum_sy

!   Ya tengo la emisión total no absorbida:

    open(17,file='radsyn1.dat')
    do i=1,mnu
       write(17,*) Eps(i)*mec2/eV,Lsy_tot(i)
    enddo
    close(17)

    open(79,file='emisyn1.dat')
    do l=1,ml/2
       do i=1,mnu
          write(79,*) l,Eps(i)*mec2/eV,Eps(i)*Lsy(l,i)
       enddo
       write(79,*) 
    enddo
    do l=ml/2+1,ml
       do i=1,mnu
          write(79,*) l,Eps(i)*mec2/eV,Eps(i)*Lsy(l-ml/2,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    deallocate ( Lsy_tot,nu,Lsy )
    deallocate ( x,y,B_par,B,rho,P,v,T )
    deallocate ( N_e,E_e,N_etot )

  END subroutine CD_radsyn1_run
  
  

  !============================================================================
  ! CD_radsyn2_run
  !============================================================================
  

  subroutine CD_radsyn2_run  
    use global
    !$ use omp_lib
    implicit none
    ! ml - CD, mE - electrones, mnu - fotones emitidos/tabla abs_ff
    integer :: il
    real(dp) :: Eps(mnu),dEps,Eps_int, gamma(mE)
    real(dp) :: factor,nuc,sum,sum_sy,ni,ne,nu_RT
    ! set omp threads number
    !$ call omp_set_num_threads(nth)   

    allocate ( x(ml),y(ml),B(ml),rho(ml),P(ml),v(ml),T(ml) )
    allocate ( N_e(mE),E_e(mE),B_par(ml),N_etot(ml,mE) )
    allocate ( Lsy_tot(mnu),nu(mnu),Lsy(ml,mnu) )
    Lsy_tot = 0.d0
    
    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y) de la curva y las variables 
    ! termodinámicas en esos puntos (x(l),y(l))

    open(11,file='termo2.dat')
    do l=1,ml
       read(11,*) il,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
    enddo
    close(11)
    B_par = B*par
    
    open(12,file='CD.dat')
    do l=1,ml
       read(12,*) x(l),y(l) ! la trayectoria en cm
    enddo
    close(12)
    
    ! Distribución de e- en el punto r(l) con energía E_e(i): N_etot(l,i) 
    open (14,file='dist2.dat')
    do l=1,ml
       do i=1,mE
          read(14,*) j, E_e(i),N_etot(l,i)
          E_e(i) = E_e(i)*eV    !   Pasamos la energía a ergios
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma = E_e/mec2
    
    !------------------------------------------------------------------------
    !
    !  Comienzo los cálculos    
    !

    factor = h/mec2
    ! Uso el mismo rango de frecuencias para TODAS las posiciones, para eso
    ! tomo la más máxima y la más mínima de todas las frecuencias, con un margen
    
    open(66,file='nu.dat') ! este lo crea el radsyn1
    do i=1,mnu
       read(66,*) nu(i)
       Eps(i) = nu(i)*factor
    enddo
    close(66)
    Eps_int = Eps(2)/Eps(1)

    ! Para cada posición voy a calcular la emisión sinc y luego sumo todas
    sum_sy=0.d0
    DO l=1,ml/2
       do i=1,mE
          N_e(i)=N_etot(l,i)
       enddo

       ni = rho(l)/mp/2.085d0 ! 2.085 = x + 4*y + 14*z
       ne = 1.395d0*ni         ! 1.395 = x + 2*y + 7*z
       nu_RT = 20.d0*ne/B(l)  ! frecuencia característica para efecto R-T

       ! Barremos todos los rangos de energía de los fotones, Epsilon (Eps), 
       ! entre Eps_min y Eps_max, dados por nu_c(gamma_min) y nu_c(gamma_max)
       do i = 1,mnu           ! Eps(i), nu(i)
          dEps = Eps(i)*(Eps_int-1.d0)
          sum = 0.d0
          ! Usamos E como variable de integración, P_eps = P(E,epsilon),
          ! Lsy = P(epsilon) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps],
          ! x = nu/nu_c va entre 0 y 10

          !$omp parallel do private(nuc,dP_Eps,P_Eps,dE_e) &
          !$omp reduction(+ : sum)
          do j = 1,mE         ! Ne(j), nuc(j), gamma(j)
             nuc = 4.22d6*B_par(l)*gamma(j)**2
             if (nu(i)/nuc.gt.10) then
                dP_Eps = 0.d0
             else
                ! Calculamos P_eps [erg/s/Hz] y pasamos a [1/s] dividiendo por h
                P_Eps = (4.39d-22*B_par(l)*(nu(i)/nuc)**(1.d0/3.d0)*exp(-nu(i)/nuc))/h
                dE_e = E_e(j)*(Eint_e - 1.d0)
                dP_Eps = N_e(j)*P_eps*dE_e
             endif
             sum = dP_Eps + sum
          enddo  
          !$omp end parallel do

          Lsy(l,i) = sum*mec2*exp(-nu_RT/nu(i)) ! Lsy en esa posición y a esa 
          !energía, corregida por efecto R-T   
          Lsy_tot(i) = Lsy_tot(i) + Eps(i)*Lsy(l,i) ! [erg/s]
          sum_sy = sum_sy + Lsy(l,i)*dEps

       enddo  ! Eps(i),nu(i)
       
    ENDDO     ! l
    
    write(*,*) 'L_syn(unabsorbed)=',sum_sy

!   Ya tengo la emisión total no absorbida:

    open(17,file='radsyn2.dat')
    do i=1,mnu
       write(17,*) Eps(i)*mec2/eV,Lsy_tot(i)
    enddo
    close(17)

    open(79,file='emisyn2.dat')
    do l=1,ml/2
       do i=1,mnu
          write(79,*) l,Eps(i)*mec2/eV,Eps(i)*Lsy(l,i)
       enddo
       write(79,*) 
    enddo
    do l=ml/2+1,ml
       do i=1,mnu
          write(79,*) l,Eps(i)*mec2/eV,Eps(i)*Lsy(l-ml/2,i)
       enddo
       write(79,*) 
    enddo
    close(79)
    
    deallocate ( Lsy_tot,nu,Lsy )
    deallocate ( x,y,B_par,B,rho,P,v,T )
    deallocate ( N_e,E_e,N_etot )

  END subroutine CD_radsyn2_run
  
  
!==============================

  
END module CD_radsyn



