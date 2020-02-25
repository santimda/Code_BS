module CD_radic
  !Este módulo calcula la emisión IC ya corregida por absorción
  use global
  
  implicit none
  private
  public :: CD_radic1_run, CD_radic2_run
  
  ! La única diferencia entre radic1 y radic2 es que uno lee los datos de 
  ! dist1 y el otro de dist2, además de los nombres de los archivos de salida.

contains
  
  !============================================================================
  ! CD_radic1_run
  !============================================================================
  
  subroutine CD_radic1_run  
    use global
    implicit none
    real(dp) :: Eint_ph1,Eint_ph2,Eph0_1,Eph0_2,E_int_low
    real(dp) :: ang_min,ang_int,ang1,ang2,d1,d2
    real(dp) :: L_unabs,L_unabs_tot,L_abs_tot
    integer :: iang1,iang2,iEf,nl
    real(dp) :: E_ph1(mE_ph),E_ph2(mE_ph),Ef_low(mE_low),Ef(mE)
    allocate ( LIC(mE_tot),LIC_tot(mE_tot),LIC_tot_abs(mE_tot) )
    allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml),B_par(ml) )
    allocate ( E_e(mE),N_etot(ml,mE),N_e(mE) )
    allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

    LIC_tot_abs = 0.d0
    LIC_tot = 0.d0
    LIC = 0.d0 

 
    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))

      open(11,file='termo1.dat')
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
         B_par(l) = sqrt(2.d0/3.d0)*B(l)
      enddo
      close(11)

      open(12,file='CDrot.dat')
      do l=1,ml
         read(12,*) x(l),y(l),z(l)
      enddo
      close(12)

!     Ahora la distribución de e- en el punto r(l) con energía Ee(i): N(l,i)
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

!     Dejamos todo lo necesario para posteriores cálculos de IC y abs g-g

!     Ephmin_i = Eph0_i/9.d0,      Ephmax_i = Eph0_i*5.d0
      Eph0_1 = 3.d0*k*T1       ! Nos interesa solo cerca del pico de la dist.
      call vectorE(Eph0_1/9.d0,Eph0_1*5.d0,Eint_ph1,E_ph1)

      Eph0_2 = 3.d0*k*T2       ! Nos interesa solo cerca del pico de la dist.
      call vectorE(Eph0_2/9.d0,Eph0_2*5.d0,Eint_ph2,E_ph2)

!      Ef_min_low = 5.d0*eV, Ef_max_low = E_e(1) = Ef_min
      call vectorE( 5.d0*eV,0.99d0*E_e(1),E_int_low,Ef_low)
      
      Ef = E_e 

      call vectorE(Ef_min_gg,Ef_max_gg,Ef_int_gg,Ef_gg)

!     Leemos una matriz con los coeficientes de opacidad a la distancia de un 
!     radio estelar del centro de la st i para diferentes ángulos entre el eje 
!     que une el emisor con la st y la dirección del observador, y energías de
!     los fotones gamma: tau(Ef,ángulo)

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
      
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !                                                                      C
      !                       CÁLCULOS                                       C
      !                                                                      C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Tengo que calcular la emisión IC de toda la CD. Para eso hago un do 
      ! en l, en cada punto calculo IC y lo corrijo por abs g-g, y luego sumo.

       DO l=1,ml
         do i=1,mE
            N_e(i)=N_etot(l,i)
         enddo
         
         Xemi = (/ x(l),y(l),z(l) /) 
         
         ! Voy a calcular la emisión por IC de la estrella 1 y luego la emisión
         ! de la estrella 2. Finalmente corrijo por la absorción de las 2 st

         ! E_ph es la de los fotones objetivo (st 1 o 2) y E_f la de los fotones
         ! emitidos por la fuente de rayos gamma. Para disminuir el ruido, 
         ! calculo E*F_E para los mismos valores de E que los valores de E 
         ! para los que calculé N_e(E).


         !!////////////////////////////////////////////////////////////////
         !!                                                              //
         !!               Estrella 1                                     //
         !!                                                              //
         d1 = sqrt(Xemi(1)**2+Xemi(2)**2+Xemi(3)**2)
         ! Calculamos el versor: XEe/|XEe|
         call versorresta(Xemi,Xst1,3,XEe) ! va desde el 2do hasta el 1ro
         ang1 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
         iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1

         !!//////////////////////////////////////////////////////////////////
         !!                                                                //
         !!                 Estrella 2                                     //
         !!                                                                //
         d2 = sqrt((D-Xemi(1))**2+Xemi(2)**2+Xemi(3)**2)
         ! Calculamos el versor: XEe/|XEe|
         call versorresta(Xemi,Xst2,3,XEe) ! va desde el 2do hasta el 1ro
         ang2 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
         iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1

         call IC(E_e,N_e,E_ph1,Ef,Ef_low,ang1,LIC,d1,Rst1,T1,L_unabs) ! St1
         do i=1,mE_tot
            LIC_tot(i) = LIC_tot(i)+LIC(i)
            if (Ef_min_gg < Ef(i-mE_low) .and. Ef(i-mE_low) < Ef_max_gg) then
               iEf = int(log(Ef(i-mE_low)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_tot_abs(i)=LIC_tot_abs(i) +LIC(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) - (Rst2/d2)*tau2(iEf,iang2) )
            else
               LIC_tot_abs(i) = LIC_tot(i)
            endif
         enddo

         call IC(E_e,N_e,E_ph2,Ef,Ef_low,ang2,LIC,d2,Rst2,T2,L_unabs) ! St2
         do i=1,mE_tot
            LIC_tot(i) = LIC_tot(i)+LIC(i)
            if (Ef_min_gg < Ef(i-mE_low) .and. Ef(i-mE_low) < Ef_max_gg) then
               iEf = int(log(Ef(i-mE_low)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_tot_abs(i)=LIC_tot_abs(i) +LIC(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) - (Rst2/d2)*tau2(iEf,iang2) )
            else
               LIC_tot_abs(i) = LIC_tot(i)
            endif
         enddo
         
         
         if (mod(l,200)==0) then
            write(*,*) 'IC, paso:',l
         endif
         
      ENDDO


!     Ya tengo la emisión IC y la emisión IC corregida por absorción g-g:

      L_abs_tot = 0.d0  
      L_unabs_tot = 0.d0
      open(16,file='radic1.dat')
      do i=1,mE_tot
         if (i.le.mE_low) then
            write(16,*) Ef_low(i)/eV,Ef_low(i)*LIC_tot(i),Ef_low(i)*LIC_tot_abs(i)
            L_unabs_tot = L_unabs_tot + LIC_tot(i)*Ef_low(i)*(E_int_low-1.d0)
            L_abs_tot = L_abs_tot + LIC_tot_abs(i)*Ef_low(i)*(E_int_low-1.d0)
         else
            write(16,*) Ef(i-mE_low)/eV,Ef(i-mE_low)*LIC_tot(i),Ef(i-mE_low)*LIC_tot_abs(i)
            L_unabs_tot = L_unabs_tot + LIC_tot(i)*Ef(i-mE_low)*(Eint_e-1.d0)
            L_abs_tot = L_abs_tot + LIC_tot_abs(i)*Ef(i-mE_low)*(Eint_e-1.d0)
         endif
      enddo
      close(16)
      write(*,*) 'L_IC(unabs) = ', L_unabs_tot
      write(*,*) 'L_IC(abs) = ', L_abs_tot

    deallocate ( LIC,LIC_tot,LIC_tot_abs )
    deallocate ( tau1,tau2 )
    deallocate ( x,y,z,rho,P,v,T,B,B_par )
    deallocate ( E_e,N_etot,N_e )

    END subroutine CD_radic1_run





!/////////////////////////////////////////////////////////////////////////////

  
  !============================================================================
  ! CD_radic2_run
  !============================================================================
  
    subroutine CD_radic2_run  
      use global
      implicit none
      real(dp) :: Eint_ph1,Eint_ph2,Eph0_1,Eph0_2,E_int_low
      real(dp) :: ang_min,ang_int,ang1,ang2,d1,d2
      real(dp) :: L_unabs,L_unabs_tot,L_abs_tot
      integer :: iang1,iang2,iEf,nl
      real(dp) :: E_ph1(mE_ph),E_ph2(mE_ph),Ef_low(mE_low),Ef(mE)
      
      allocate ( LIC(mE_tot),LIC_tot(mE_tot),LIC_tot_abs(mE_tot) )
      allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )
      allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml),B_par(ml) )
      allocate ( E_e(mE),N_etot(ml,mE),N_e(mE) )
      
      LIC_tot_abs = 0.d0
      LIC_tot = 0.d0
      LIC = 0.d0 

 
    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))

      open(11,file='termo2.dat')
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
         B_par(l) = sqrt(2.d0/3.d0)*B(l)
      enddo
      close(11)

      open(12,file='CDrot.dat')
      do l=1,ml
         read(12,*) x(l),y(l),z(l)
      enddo
      close(12)

!     Ahora la distribución de e- en el punto r(l) con energía Ee(i): N(l,i)
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

!     Dejamos todo lo necesario para posteriores cálculos de IC y abs g-g

!      Ef_min_low ~ kT, Ef_max_low ~ E_e(1) = Ef_min
      call vectorE( k*T1,0.99d0*E_e(1),E_int_low,Ef_low)
      
      Ef = E_e 

      call vectorE(Ef_min_gg,Ef_max_gg,Ef_int_gg,Ef_gg)

!     Leemos una matriz con los coeficientes de opacidad a la distancia de un 
!     radio estelar del centro de la st i para diferentes ángulos entre el eje 
!     que une el emisor con la st y la dirección del observador, y energías de
!     los fotones gamma: tau(Ef,ángulo)

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

      Xobs = (/ cos(incli),sin(incli),0.d0 /)
      Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
      Xst2 = (/ D, 0.d0, 0.d0 /)
      
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !                                                                      C
      !                       CÁLCULOS                                       C
      !                                                                      C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Tengo que calcular la emisión IC de toda la CD. Para eso hago un do 
      ! en l, en cada punto calculo IC y lo corrijo por abs g-g, y luego sumo.

      DO l=1,ml
         do i=1,mE
            N_e(i)=N_etot(l,i)
         enddo
         
         Xemi = (/ x(l),y(l),z(l) /) 
         
         ! Voy a calcular la emisión por IC de la estrella 1 y luego la emisión
         ! de la estrella 2. Finalmente corrijo por la absorción de las 2 st

         ! E_ph es la de los fotones objetivo (st 1 o 2) y E_f la de los fotones
         ! emitidos por la fuente de rayos gamma. Para disminuir el ruido, 
         ! calculo E*F_E para los mismos valores de E que los valores de E 
         ! para los que calculé N_e(E).


         !!////////////////////////////////////////////////////////////////
         !!                                                              //
         !!               Estrella 1                                     //
         !!                                                              //
         d1 = sqrt(Xemi(1)**2+Xemi(2)**2+Xemi(3)**2)
         ! Calculamos el versor: XEe/|XEe|
         call versorresta(Xemi,Xst1,3,XEe) ! va desde el 2do hasta el 1ro
         ang1 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
         iang1 = int(m_ang_gg*abs(ang1-ang_min)/pi)+1

         !!//////////////////////////////////////////////////////////////////
         !!                                                                //
         !!                 Estrella 2                                     //
         !!                                                                //
         d2 = sqrt((D-Xemi(1))**2+Xemi(2)**2+Xemi(3)**2)
         ! Calculamos el versor: XEe/|XEe|
         call versorresta(Xemi,Xst2,3,XEe) ! va desde el 2do hasta el 1ro
         ang2 = acos(DOT_PRODUCT(Xobs,XEe))! ángulo comprendido entre Xobs y XEe
         iang2 = int(m_ang_gg*abs(ang2-ang_min)/pi)+1

         call IC(E_e,N_e,E_ph1,Ef,Ef_low,ang1,LIC,d1,Rst1,T1,L_unabs) ! St1
         do i=1,mE_tot
            LIC_tot(i) = LIC_tot(i)+LIC(i)
            if (Ef_min_gg < Ef(i-mE_low) .and. Ef(i-mE_low) < Ef_max_gg) then
               iEf = int(log(Ef(i-mE_low)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_tot_abs(i)=LIC_tot_abs(i) +LIC(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) - (Rst2/d2)*tau2(iEf,iang2) )
            else
               LIC_tot_abs(i) = LIC_tot(i)
            endif
         enddo

         call IC(E_e,N_e,E_ph2,Ef,Ef_low,ang2,LIC,d2,Rst2,T2,L_unabs) ! St2
         do i=1,mE_tot
            LIC_tot(i) = LIC_tot(i)+LIC(i)
            if (Ef_min_gg < Ef(i-mE_low) .and. Ef(i-mE_low) < Ef_max_gg) then
               iEf = int(log(Ef(i-mE_low)/Ef_min_gg)/log(Ef_int_gg)) + 1
               LIC_tot_abs(i)=LIC_tot_abs(i) +LIC(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) - (Rst2/d2)*tau2(iEf,iang2) )
            else
               LIC_tot_abs(i) = LIC_tot(i)
            endif
         enddo
         

         if (mod(l,200)==0) then
            write(*,*) 'IC, paso:',l
         endif

      ENDDO


!     Ya tengo la emisión IC y la emisión IC corregida por absorción g-g:


      L_abs_tot = 0.d0  
      L_unabs_tot = 0.d0
      open(17,file='radic2.dat')
      do i=1,mE_tot
         if (i.le.mE_low) then
            write(17,*) Ef_low(i)/eV,Ef_low(i)*LIC_tot(i),Ef_low(i)*LIC_tot_abs(i)
            L_unabs_tot = L_unabs_tot + LIC_tot(i)*Ef_low(i)*(E_int_low-1.d0)
            L_abs_tot = L_abs_tot + LIC_tot_abs(i)*Ef_low(i)*(E_int_low-1.d0)
         else
            write(17,*) Ef(i-mE_low)/eV,Ef(i-mE_low)*LIC_tot(i),Ef(i-mE_low)*LIC_tot_abs(i)
            L_unabs_tot = L_unabs_tot + LIC_tot(i)*Ef(i-mE_low)*(Eint_e-1.d0)
            L_abs_tot = L_abs_tot + LIC_tot_abs(i)*Ef(i-mE_low)*(Eint_e-1.d0)
         endif
      enddo
      close(17)
      write(*,*) 'L_IC(unabs) = ', L_unabs_tot
      write(*,*) 'L_IC(abs) = ', L_abs_tot

      deallocate ( LIC,LIC_tot,LIC_tot_abs )
      deallocate ( tau1,tau2 )
      deallocate ( x,y,z,rho,P,v,T,B,B_par )
      deallocate ( E_e,N_etot,N_e )

    END subroutine CD_radic2_run

    !___________
    
  END module CD_radic




!/////////////////////////////////////////////////////////////////////////////
!                                                                           //
!                              SUBRUTINAS                                   //
!                                                                           //
!/////////////////////////////////////////////////////////////////////////////


    SUBROUTINE IC(Ee,Ne,Eph,Ef,Eflow,ang,L_IC,dist,Rst,Tst,L_unabs)
      use global
      Implicit none
      Real(dp), intent(in) :: Ee,Ne,Eph,Ef,Eflow,ang,dist,Rst,Tst
      Real(dp), intent(out) :: L_IC(mE_tot),L_unabs
      Real(dp) :: Eps,Ephint,Eeint,Epsint
      Real(dp) :: dEps,dEe,dEph,dP_IC,sum,gamma,P_epsIC
      Real(dp) :: Tmcc,Tmcc_ang,Kappa,cte_Ndot_TK,z_K,t_K,x0_K,Ndot
      Dimension :: Ee(mE),Ne(mE),Eph(mE_ph),Ef(mE),Eflow(mE_low)

      Tmcc = k*Tst/mec2
      Tmcc_ang = Tmcc*(1.-cos(ang))
      Kappa = (Rst/(2.*dist))**2.
      cte_Ndot_TK = cte_Ndot*Tmcc**2.*Kappa
      
      Eeint = Ee(2)/Ee(1)
      Ephint = Eph(2)/Eph(1)

!------------------------------------------------------------
!     Calculo primero para energías chicas (desprecio abs):
      Epsint = Eflow(2)/Eflow(1)
!     Eph será la energía del fotón incidente, Eps la del fotón emitido y 
!     gamma la energía del electrón, todas se usan en unidades de mec2
      L_unabs = 0.d0
      ! En IC, la interaccion es con los fotones de distribucion nBB

      do i=1,mE_low              ! Ef
         Eps = Eflow(i)/mec2
         dEps=Eps*(Epsint-1.d0)
         sum=0.d0
         do j = 1,mE            ! Ee
            gamma = Ee(j)/mec2  ! gm
            z_K = Eps/gamma
            t_K = 2.*gamma*Tmcc_ang
            x0_K = z_K/( (1.-z_K)*t_K )
            Ndot = (cte_Ndot_TK/gamma**2.) * ( (z_K**2./(2.*(1.-z_K)))*F_1(x0_K) + F_2(x0_K) )
            P_epsIC = Eps*Ndot
            dEe = Ee(j)*(Eeint-1.d0)
            dP_IC = Ne(j)*P_epsIC*dEe
            sum = dP_IC + sum
         enddo      
         L_IC(i) = sum           ! Es la P , luego grafico E(j)*LIC(Ej)
         L_unabs = L_unabs + sum*dEps*mec2
      enddo


!-----------------------------

!   Ahora calculo para energías > 1MeV
      
      Epsint = Ef(2)/Ef(1)
!     Eph será la energía del fotón incidente, Eps la del fotón emitido y 
!     gamma la energía del electrón, todas se usan en unidades de mec2
      do i=1,mE            ! Ef
         Eps = 0.999*Ef(i)/mec2
         dEps = Eps*(Epsint-1.d0)
         sum = 0.d0
         do j = 1,mE            ! E
            gamma = Ee(j)/mec2  ! gm
            z_K = Eps/gamma
            if (z_K < 1.d0) then
               t_K = 2.*gamma*Tmcc_ang
               x0_K = z_K/( (1.-z_K)*t_K )
               Ndot = (cte_Ndot_TK/gamma**2.) * ( (z_K**2./(2.*(1.-z_K)))*F_1(x0_K) + F_2(x0_K) )
               P_epsIC = Eps*Ndot
               dEe = Ee(j)*(Eeint-1.d0)
               dP_IC = Ne(j)*P_epsIC*dEe
            else
               dP_IC = 0.d0
            endif
            sum = dP_IC + sum
         enddo      
         LIC(mE_low+i) = sum     ! Es la P , luego grafico E(j)*LIC(Ej)
         L_unabs = L_unabs + sum*dEps*mec2
      enddo
      
      RETURN

    end SUBROUTINE IC


    !/////////////////////////////////////////////////////////////////////////

