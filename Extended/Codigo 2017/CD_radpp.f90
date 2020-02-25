module CD_radpp
  !Este módulo calcula la emisión pp ya corregida por abs gg
  use global
  
  implicit none
  private
  public :: CD_radpp1_run, CD_radpp2_run
  
contains
  
  !============================================================================
  ! CD_radpp1_run
  !============================================================================
  
  subroutine CD_radpp1_run  
    use global
    Implicit none
    Real(dp) :: d1,d2,nmedio(ml),Q_Ef,Q_phi,dx_pp,x_pp,qpi
    Real(dp) :: Ef_int,dEf,Ef(mE_f),sum_pp,sigma_pp,Epp,L1
    Real(dp) :: E_th= 1.22d9*eV, TeV = 1.d12*eV
    Real(dp) :: Epi_min,Epi_max,Epi_int,dE_pi,integral,F
    Real(dp) :: ang_min,ang_int,ang1,ang2
    Integer :: nl,iang1,iang2,iEf,mE_pi,in
    !     ml - CD, mE_p - protones, mE_f - fotones, mE_pi(variable) - piones
    Real(dp), dimension(:), allocatable :: Epi
    allocate ( Lpp(mE_f),Lpp_tot(mE_f),Lpp_tot_abs(mE_f) )
    allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml) )
    allocate ( E_p(mE_p),N_ptot(ml,mE_p),N_p(mE_p) )
    allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )
    
    Lpp_tot_abs = 0.d0
    Lpp_tot = 0.d0

    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))
    
    open(11,file='termo1.dat')
    do l=1,ml
       read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
       nmedio(l) = rho(l)/mp/mu_1 ! da del orden de 10^4 (lejos) - 10^6 (apex)
    enddo
    close(11)

    open(12,file='CDrot.dat')
    do l=1,ml
       read(12,*) x(l),y(l),z(l)
    enddo
    close(12)

    !     Ahora la distribución de p- en el punto r(l) con energía Ep(i): N(l,i)
    open (14,file='distp1.dat')
    do l=1,ml
       do i=1,mE_p
          read(14,*) j, E_p(i),N_ptot(l,i)
          E_p(i) = E_p(i)*eV    !   Pasamos la energía a ergios
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_p = E_p(2)/E_p(1)
    
    ! Leemos una matriz con los coeficientes de opacidad a la distancia de un 
    ! radio estelar del centro de la st i para diferentes ángulos entre el eje 
    ! que une el emisor con la st y la dirección del observador, y energías de
    ! los fotones gamma: tau(Ef,ángulo)
    
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/dfloat(mE_gg-1))
    
    open (88,file='opac_gg1.dat')
    do i=1,mE_gg
       read(88,*) (tau1(i,j), j=1,m_ang_gg)
    enddo
    close(88)
    open (88,file='opac_gg2.dat')
    do i=1,mE_gg
       read(88,*) (tau2(i,j), j=1,m_ang_gg)
    enddo
    close(88)
    ang_int = pi/float(m_ang_gg)
    ang_min = ang_int/2.d0
    
    Xobs = (/ sin(incli),cos(incli),0.d0 /)
    Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
    Xst2 = (/ D, 0.d0, 0.d0 /)
    
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                        CÁLCULOS                                       C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías de fotones, Ef(i), con Ef_min = 10 MeV
    
    call vectorE(10.d6*eV,E_p(mE_p),Ef_int,Ef)
    
    ! Para cada posición voy a calcular la emisión pp y luego sumo todas
    sum_pp=0.d0
    DO l=1,ml          ! l
       do i=1,mE_f
          N_p(i)=N_ptot(l,i)
       enddo

       Xemi = (/ x(l),y(l),z(l) /) 
       d1 = sqrt(Xemi(1)**2+Xemi(2)**2+Xemi(3)**2)
 
       ! Calculamos el versor: XEe/|XEe|
       call versorresta(Xemi,Xst1,3,XEe) ! va desde el 2do hasta el 1ro

       ang1 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
       
       d2 = sqrt((D-Xemi(1))**2+Xemi(2)**2+Xemi(3)**2)
       ! Calculamos el versor: XEe/|XEe|
       call versorresta(Xemi,Xst2,3,XEe) ! va desde el 2do hasta el 1ro
       ang2 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
 
       iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1
       iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1

       !///////////////////////////////////////////////////////////////////
       !                                                                 //
       !                    p-p                                          //
       !                                                                 //
       !///////////////////////////////////////////////////////////////////
       
       
       mE_pi = mE_f
       do i = 1,mE_f           ! Ef(i) 
          dEf = Ef(i)*(Ef_int-1.d0)

          ! Usamos un formalismo distinto según el rango de energías:
          if (Ef(i) < 0.1d0*TeV) then    ! Ef < 100 GeV
             allocate(Epi(mE_pi))
             Epi_min = Ef(i)+(m_pic2)**2.d0/4.d0/Ef(i)
             Epi_min = Epi_min*(1.d0+1.d0/mE_pi)
             Epi_max = E_p(mE_p)
             call vectorE(Epi_min,Epi_max,Epi_int,Epi)
             integral = 0.d0
             do j=1,mE_pi       ! E_pi
                dE_pi = Epi(j)*(Epi_int-1.d0)
                Epp = Epi(j)/0.17d0 + mpc2
                in = int(log(Epp/Emin_p)/log(Eint_p))+1
                if (in > mE_pi .or. Epp > Emax_p) then
                   goto 0979
                endif
                ! Rango de energías 10 MeV-100TeV; sigmaPP=0 si E < Eminp=1GeV.
                if (Epp < E_th) then !E_th= 1.22 GeV
                   sigma_PP = 0.d0
                else            
                   L1 = log(Epp/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
                   ! Agregamos un factor para que ajuste mejor a bajas energías
                   sigma_PP = sigma_PP*(1.d0-(E_th/Epp)**4)**2  
                endif
                qpi = N_p(in)*sigma_pp/0.17d0 ! 1/erg/seg
                integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-(m_pic2)**2)
             enddo              ! E_pi
0979         continue
             ! write(*,*) Ef(i), nmedio(l), integral
             Q_Ef = 2.3d0*c*nmedio(l)*integral/1.3d0! ñ=1.1 corresponde a alfa=2
             ! SED = Ef * L_Ef =  Ef^2 * Q_Ef 
             Lpp(i) = Ef(i)*Q_Ef
             sum_PP = sum_PP + Lpp(i)*dEf
             
             deallocate(Epi)
!             mE_pi = MAX(mE_pi-1,20) ! Disminuyo la cantidad de puntos 
             
          elseif (Ef(i) > 0.1d0*TeV) then  ! Ef > 100 GeV
             integral=0
             do j=1,mE_p
                if (Ef(i) < E_p(j)) then
                   dE_p = E_p(j)*(Eint_p-1.d0)
                   x_pp = Ef(i)/E_p(j)
                   dx_pp = Ef(i)*dE_p/(E_p(j)**2)
                   ! Como Ef > 100 GeV > E_th, sigma_pp es distinto de cero 
                   L1 = log(E_p(j)/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27
                   ! Agregamos un factor para que ajuste mejor a bajas E:
                   sigma_PP = sigma_PP*(1.d0-(E_th/E_p(j))**4)**2
                   integral = integral + sigma_PP*N_p(j)*F(x_pp,L1)*dx_pp/x_pp
                endif
             enddo
             Q_phi = c*nmedio(l)*integral
             Lpp(i) = Ef(i)*Q_phi
             !             write (19,*) Ef(i)/eV, Ef(i)*Lpp(i)
             sum_PP = sum_PP + Lpp(i)*dEf 
             
          endif
          
                     
          Lpp_tot(i) = Lpp_tot(i) + Lpp(i)
                   
          !Corrijo por absorción:
          if (Ef(i) > Ef_min_gg .and. Ef(i) < Ef_max_gg) then
             iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
             Lpp_tot_abs(i) = Lpp_tot_abs(i) + Lpp(i)*exp( -(Rst1/d1)*&
                  tau1(iEf,iang1) -(Rst2/d2)*tau2(iEf,iang2) )
          else
             Lpp_tot_abs(i) = Lpp_tot_abs(i) + Lpp(i)
          endif
          
       enddo    ! i
       
       
    ENDDO       ! l

    write(*,*) 'L_pp=',sum_pp

    ! Ya tengo la emisión pp y la emisión pp corregida por absorción g-g:
 
    open(16,file='radpp1.dat')
    do i=1,mE_f
       write(16,*) Ef(i)/eV,Ef(i)*Lpp_tot(i),Ef(i)*Lpp_tot_abs(i)
    enddo
    close(16)

    deallocate ( Lpp,Lpp_tot,Lpp_tot_abs)
    deallocate ( x,y,z,rho,P,v,T,B )
    deallocate ( E_p,N_ptot,N_p )
    deallocate ( tau1,tau2 )


  END subroutine CD_radpp1_run
  


  !============================================================================
  ! CD_radpp2_run
  !============================================================================
  
  subroutine CD_radpp2_run  
    use global
    Implicit none
    Real(dp) :: d1,d2,nmedio(ml),Q_Ef,Q_phi,dx_pp,x_pp,qpi
    Real(dp) :: Ef_int,dEf,Ef(mE_f),sum_pp,sigma_pp,Epp,L1
    Real(dp) :: E_th= 1.22d9*eV, TeV = 1.d12*eV
    Real(dp) :: Epi_min,Epi_max,Epi_int,dE_pi,integral,F
    Real(dp) :: ang_min,ang_int,ang1,ang2
    Integer :: nl,iang1,iang2,iEf,mE_pi,in
    !     ml - CD, mE_p - protones, mE_f - fotones, mE_pi(variable) - piones
    Real(dp), dimension(:), allocatable :: Epi
    allocate ( Lpp(mE_f),Lpp_tot(mE_f),Lpp_tot_abs(mE_f) )
    allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml) )
    allocate ( E_p(mE_p),N_ptot(ml,mE_p),N_p(mE_p) )
    allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )
    
    Lpp_tot_abs = 0.d0
    Lpp_tot = 0.d0

    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))
    
    open(11,file='termo2.dat')
    do l=1,ml
       read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
       nmedio(l) = rho(l)/mp/mu_2 ! da del orden de 10^4 (lejos) - 10^6 (apex)
    enddo
    close(11)

    open(12,file='CDrot.dat')
    do l=1,ml
       read(12,*) x(l),y(l),z(l)
    enddo
    close(12)

    !     Ahora la distribución de p- en el punto r(l) con energía Ep(i): N(l,i)
    open (14,file='distp2.dat')
    do l=1,ml
       do i=1,mE_p
          read(14,*) j, E_p(i),N_ptot(l,i)
          E_p(i) = E_p(i)*eV    !   Pasamos la energía a ergios
       enddo
       read(14,*)
    enddo
    close(14) 
    Eint_p = E_p(2)/E_p(1)
    
    ! Leemos una matriz con los coeficientes de opacidad a la distancia de un 
    ! radio estelar del centro de la st i para diferentes ángulos entre el eje 
    ! que une el emisor con la st y la dirección del observador, y energías de
    ! los fotones gamma: tau(Ef,ángulo)
    
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/dfloat(mE_gg-1))
    
    open (88,file='opac_gg1.dat')
    do i=1,mE_gg
       read(88,*) (tau1(i,j), j=1,m_ang_gg)
    enddo
    close(88)
    open (88,file='opac_gg2.dat')
    do i=1,mE_gg
       read(88,*) (tau2(i,j), j=1,m_ang_gg)
    enddo
    close(88)
    ang_int = pi/float(m_ang_gg)
    ang_min = ang_int/2.d0
    
    Xobs = (/ sin(incli),cos(incli),0.d0 /)
    Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
    Xst2 = (/ D, 0.d0, 0.d0 /)
    
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                       C
    !                        CÁLCULOS                                       C
    !                                                                       C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! Armo el vector de energías de fotones, Ef(i), con Ef_min = 10 MeV
    
    call vectorE(10.d6*eV,E_p(mE_p),Ef_int,Ef)
    
    ! Para cada posición voy a calcular la emisión pp y luego sumo todas
    sum_pp=0.d0
    DO l=1,ml          ! l
       do i=1,mE_f
          N_p(i)=N_ptot(l,i)
       enddo

       Xemi = (/ x(l),y(l),z(l) /) 
       d1 = sqrt(Xemi(1)**2+Xemi(2)**2+Xemi(3)**2)
 
       ! Calculamos el versor: XEe/|XEe|
       call versorresta(Xemi,Xst1,3,XEe) ! va desde el 2do hasta el 1ro

       ang1 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
       
       d2 = sqrt((D-Xemi(1))**2+Xemi(2)**2+Xemi(3)**2)
       ! Calculamos el versor: XEe/|XEe|
       call versorresta(Xemi,Xst2,3,XEe) ! va desde el 2do hasta el 1ro
       ang2 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
 
       iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1
       iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1

       !///////////////////////////////////////////////////////////////////
       !                                                                 //
       !                    p-p                                          //
       !                                                                 //
       !///////////////////////////////////////////////////////////////////
       
       
       mE_pi = mE_f
       do i = 1,mE_f           ! Ef(i) 
          dEf = Ef(i)*(Ef_int-1.d0)

          ! Usamos un formalismo distinto según el rango de energías:
          if (Ef(i) < 0.1d0*TeV) then    ! Ef < 100 GeV
             allocate(Epi(mE_pi))
             Epi_min = Ef(i)+(m_pic2)**2.d0/4.d0/Ef(i)
             Epi_min = Epi_min*(1.d0+1.d0/mE_pi)
             Epi_max = E_p(mE_p)
             call vectorE(Epi_min,Epi_max,Epi_int,Epi)      
             integral = 0.d0
             do j=1,mE_pi       ! E_pi
                dE_pi = Epi(j)*(Epi_int-1.d0)
                Epp = Epi(j)/0.17d0 + mpc2
                in = int(log(Epp/Emin_p)/log(Eint_p))+1
                if (in > mE_pi .or. Epp > Emax_p) goto 0979
                ! Rango de energías 10 MeV-100TeV; sigmaPP=0 si E < Eminp=1GeV.
                if (Epp < E_th) then !E_th= 1.22 GeV
                   sigma_PP = 0.d0
                else            
                   L1 = log(Epp/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
                   ! Agregamos un factor para que ajuste mejor a bajas energías
                   sigma_PP = sigma_PP*(1-(E_th/Epp)**4)**2  
                endif
                qpi = N_p(in)*sigma_pp/0.17d0 ! 1/erg/seg
                integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-(m_pic2)**2)
             enddo              ! E_pi
0979         continue
             ! write(*,*) Ef(i), nmedio(l), integral
             Q_Ef = 2.2d0*c*nmedio(l)*integral/1.3d0! ñ=1.1 corresponde a alfa=2
             ! 
             ! SED = Ef * L_Ef =  Ef^2 * Q_Ef 
             Lpp(i) = Ef(i)*Q_Ef
             sum_PP = sum_PP + Lpp(i)*dEf
             
             deallocate(Epi)
!             mE_pi = MAX(mE_pi-1,20) ! Disminuyo la cantidad de puntos 
             
          elseif (Ef(i) > 0.1d0*TeV) then  ! Ef > 100 GeV
             integral=0
             do j=1,mE_p
                if (Ef(i) < E_p(j)) then
                   dE_p = E_p(j)*(Eint_p-1.d0)
                   x_pp = Ef(i)/E_p(j)
                   dx_pp = Ef(i)*dE_p/(E_p(j)**2)
                   ! Como Ef > 100 GeV > E_th, sigma_pp es distinto de cero 
                   L1 = log(E_p(j)/TeV) ! L = ln(Ep/1TeV)
                   sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27
                   ! Agregamos un factor para que ajuste mejor a bajas E:
                   sigma_PP = sigma_PP*(1-(E_th/E_p(j))**4)**2
                   integral = integral + sigma_PP*N_p(j)*F(x_pp,L1)*dx_pp/x_pp
                endif
             enddo
             Q_phi = c*nmedio(l)*integral
             Lpp(i) = Ef(i)*Q_phi
             !             write (19,*) Ef(i)/eV, Ef(i)*Lpp(i)
             sum_PP = sum_PP + Lpp(i)*dEf 
             
          endif
          
                     
          Lpp_tot(i) = Lpp_tot(i) + Lpp(i)
                   
          !Corrijo por absorción:
          if (Ef(i) > Ef_min_gg .and. Ef(i) < Ef_max_gg) then
             iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
             Lpp_tot_abs(i) = Lpp_tot_abs(i) + Lpp(i)*exp( -(Rst1/d1)*&
                  tau1(iEf,iang1) -(Rst2/d2)*tau2(iEf,iang2) )
          else
             Lpp_tot_abs(i) = Lpp_tot_abs(i) + Lpp(i)
          endif
          
       enddo    ! i
       
       
    ENDDO       ! l

    write(*,*) 'L_pp=',sum_pp

    ! Ya tengo la emisión pp y la emisión pp corregida por absorción g-g:
 
    open(16,file='radpp2.dat')
    do i=1,mE_f
       write(16,*) Ef(i)/eV,Ef(i)*Lpp_tot(i),Ef(i)*Lpp_tot_abs(i)
    enddo
    close(16)

    deallocate ( Lpp,Lpp_tot,Lpp_tot_abs)
    deallocate ( x,y,z,rho,P,v,T,B )
    deallocate ( E_p,N_ptot,N_p )
    deallocate ( tau1,tau2 )

  
  END subroutine CD_radpp2_run

  
  
  
  !========================================
  
  
END module CD_radpp


!========================================

Function F(x,L)
  integer, parameter :: dp = kind(1.d0)
  real(dp) :: x,B,k,beta,L,F
  B = 1.30d0 + 0.14d0*L + 0.011d0*L**2
  beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L**2)
  k = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L**2)
  F = B*log(x)/x*((1.d0-x**beta)/(1.d0+k*x**beta*(1.d0-x**beta)))**4.d0* &
       (1.d0/log(x)-(4.d0*beta*x**beta)/(1.d0-x**beta)- &
       (4.d0*k*beta*x**beta*(1.d0-2.d0*x**beta)/(1.d0+k*x**beta*(1.d0-x**beta))))
  return
end Function F

