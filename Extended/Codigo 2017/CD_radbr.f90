module CD_radbr
  !Este módulo calcula la emisión bremsstrahlung ya corregida por abs gg
  use global
  
  implicit none
  private
  public :: CD_radbr1_run, CD_radbr2_run
  
contains
  
  !============================================================================
  ! CD_radbr1_run
  !============================================================================
  
  subroutine CD_radbr1_run  
    use global
    Implicit none
      Real(dp) :: d1,d2,nmedio(ml)
      Real(dp) :: Eps_int,dEps,sum,sum_br,Phi_Br,dsigma_Br
      Real(dp) :: gamma(mE_f), Eps(mE_f), Ef(mE_f)
      Real(dp) :: ang_min,ang_int,ang1,ang2,cte_sigma
      Integer :: nl,iang1,iang2,iEf
!     ml - CD, mE - electrones, mE_f - fotones

      allocate ( Lbr(mE_f),Lbr_tot(mE_f),Lbr_tot_abs(mE_f) )
      allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml),B_par(ml) )
      allocate ( E_e(mE),N_etot(ml,mE),N_e(mE) )
      allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

      Lbr_tot_abs = 0.d0
      Lbr_tot = 0.d0
      
      cte_sigma = 4.d0*alfa*qe**4.d0

    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))

      open(11,file='termo1.dat')
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
         B_par(l) = sqrt(2.d0/3.d0)*B(l)
         nmedio(l) = rho(l)/mp/mu_1 ! da del orden de 10^4 (lejos) - 10^6 (apex)
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
 
!     Leemos una matriz con los coeficientes de opacidad a la distancia de un 
!     radio estelar del centro de la st i para diferentes ángulos entre el eje 
!     que une el emisor con la st y la dirección del observador, y energías de
!     los fotones gamma: tau(Ef,ángulo)

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

      Xobs = (/ cos(incli),sin(incli),0.d0 /)
      Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
      Xst2 = (/ D, 0.d0, 0.d0 /)
 
     
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     !                                                                       C
     !                        CÁLCULOS                                       C
     !                                                                       C
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Armo el vector de energías de fotones, Eps(i), con Eps_min del orden 
      ! del keV (si gamma < 2 ~ 1MeV, el Brems es NR) 

      call vectorE((E_e(1)/1.d3)/mec2,E_e(mE)/mec2,Eps_int,Eps)
      Ef = Eps*mec2
      gamma = E_e/mec2
      
      ! Para cada posición voy a calcular la emisión brem y luego sumo todas
      sum_br=0.d0
      DO l=1,ml
         do i=1,mE
            N_e(i)=N_etot(l,i)
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
         !                    Brem                                         //
         !                                                                 //
         !///////////////////////////////////////////////////////////////////
         
         do i = 1,mE_f           ! Eps(i) 
            dEps = Eps(i)*(Eps_int-1.d0)
            sum = 0.d0
            do j = 1,mE            ! N_e(j),gamma(j)
               !Pedimos que Eps <= gamma por conservación de energía
               if (gamma(j) < Eps(i)) then
                  dsigma_Br = 0.d0  ![cm**2/m_e c2]
               else
                  ! Apantallado:
                  !     Phi=(1+(1-Eps/gamma)**2-2.d0/3.d0*(1-Eps/gamma))*&
                  !     log(191/z**(1.d0/3.d0))+1.d0/9.d0*(1-Eps/gamma)
                  
                  ! No apantallado:
                  Phi_Br=(1.d0+(1.d0-Eps(i)/gamma(j))**2.d0-2.d0/3.d0*&
                       (1.d0-Eps(i)/gamma(j)))*&
                       (log(2.d0*gamma(j)*(gamma(j)-Eps(i))/Eps(i))-0.5d0)
                  dsigma_Br= cte_sigma/(Eps(i)*(mec2**2.d0))*Phi_Br
                  ! [P_eps] = 1/seg
                  ! Factor: Sum ni*Zi**2 ~ 0.705*1**2 + 0.275*2**2 + 0.020*7**2 ~ 2.785
                  P_eps = Eps(i)*2.785d0*nmedio(l)*c*dsigma_Br
                  dE_e = E_e(j)*(Eint_e-1.d0)
                  dP_Eps = N_e(j)*P_eps*dE_e
                  sum = dP_Eps + sum ! P = int dP_Eps
               endif
            enddo ! j
            Lbr(i) = sum!*mec2         ! [erg/seg]
            Lbr_tot(i) = Lbr_tot(i) + Lbr(i)

            if (Ef(i) > Ef_min_gg .and. Ef(i) < Ef_max_gg) then
               iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
               Lbr_tot_abs(i) = Lbr_tot_abs(i) + Lbr(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) -(Rst2/d2)*tau2(iEf,iang2) )
            else
               Lbr_tot_abs(i) = Lbr_tot_abs(i) + Lbr(i)
            endif

            sum_br = sum_br + Lbr(i)*dEps
            
         enddo ! i
               
      ENDDO ! l
      
      write(*,*) 'L_Brem=',sum_br
           
!     Ya tengo la emisión brems y la emisión brems corregida por absorción g-g:
      
      open(16,file='radbr1.dat')
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*Lbr_tot(i),Ef(i)*Lbr_tot_abs(i)
      enddo
      close(16)
      
      deallocate ( x,y,z,rho,P,v,T,B,B_par )
      deallocate ( E_e,N_etot,N_e )
      deallocate ( Lbr,Lbr_tot,Lbr_tot_abs )
      deallocate ( tau1,tau2 )

    END subroutine CD_radbr1_run
    

  
  !============================================================================
  ! CD_radbr2_run
  !============================================================================
  
    subroutine CD_radbr2_run  
      use global
      Implicit none
      Real(dp) :: d1,d2,nmedio(ml)
      Real(dp) :: Eps_int,dEps,sum,sum_br,Phi_Br,dsigma_Br
      Real(dp) :: gamma(mE_f), Eps(mE_f), Ef(mE_f)
      Real(dp) :: ang_min,ang_int,ang1,ang2,cte_sigma
      Integer :: nl,iang1,iang2,iEf
      !     ml - CD, mE - electrones, mE_f - fotones

      allocate ( Lbr(mE_f),Lbr_tot(mE_f),Lbr_tot_abs(mE_f) )
      allocate ( x(ml),y(ml),z(ml),rho(ml),P(ml),v(ml),T(ml),B(ml),B_par(ml) )
      allocate ( E_e(mE),N_etot(ml,mE),N_e(mE) )
      allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

      Lbr_tot_abs = 0.d0
      Lbr_tot = 0.d0
        
      cte_sigma = 4.d0*alfa*qe**4.d0

    !////////////////////////////////////////////////////////////////////////
    ! Lectura de datos previos: Los puntos (x,y,z) de la curva junto con las 
    ! variables termodinámicas en esos puntos (x(l),y(l),z(l))

      open(11,file='termo2.dat')
      do l=1,ml
         read(11,*) nl,rho(l),P(l),v(l),T(l),B(l) ! todo en CGS
         B_par(l) = sqrt(2.d0/3.d0)*B(l)
         nmedio(l) = rho(l)/mp/mu_2 ! da del orden de 10^4 (lejos) - 10^6 (apex)
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
 
!     Leemos una matriz con los coeficientes de opacidad a la distancia de un 
!     radio estelar del centro de la st i para diferentes ángulos entre el eje 
!     que une el emisor con la st y la dirección del observador, y energías de
!     los fotones gamma: tau(Ef,ángulo)

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

      Xobs = (/ cos(incli),sin(incli),0.d0 /)
      Xst1 = (/ 0.d0, 0.d0, 0.d0 /)
      Xst2 = (/ D, 0.d0, 0.d0 /)
 
     
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     !                                                                       C
     !                        CÁLCULOS                                       C
     !                                                                       C
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Armo el vector de energías de fotones, Eps(i), con Eps_min del orden 
      ! del keV (si gamma < 2 ~ 1MeV, el Brems es NR) 

      call vectorE((E_e(1)/1.0d3)/mec2,E_e(mE)/mec2,Eps_int,Eps)
      Ef = Eps*mec2
      gamma = E_e/mec2
      
      ! Para cada posición voy a calcular la emisión brem y luego sumo todas
      sum_br=0.d0
      DO l=1,ml
         do i=1,mE
            N_e(i)=N_etot(l,i)
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
         !                    Brem                                         //
         !                                                                 //
         !///////////////////////////////////////////////////////////////////
         
         do i = 1,mE_f           ! Eps(i) 
            dEps = Eps(i)*(Eps_int-1.d0)
            sum = 0.d0
            do j = 1,mE            ! N_e(j),gamma(j)
               !Pedimos que Eps <= gamma por conservación de energía
               if (gamma(j) < Eps(i)) then
                  dsigma_Br = 0.d0  ![cm**2/m_e c2]
               else
                  ! Apantallado:
                  !     Phi=(1+(1-Eps/gamma)**2-2.d0/3.d0*(1-Eps/gamma))*&
                  !     log(191/z**(1.d0/3.d0))+1.d0/9.d0*(1-Eps/gamma)
                  
                  ! No apantallado:
                  Phi_Br=(1.d0+(1.d0-Eps(i)/gamma(j))**2.d0-2.d0/3.d0*&
                       (1.d0-Eps(i)/gamma(j)))*&
                       (log(2.d0*gamma(j)*(gamma(j)-Eps(i))/Eps(i))-0.5d0)
                  dsigma_Br = cte_sigma/(Eps(i)*(mec2**2.d0))*Phi_Br
                  ! [P_eps] = 1/seg
                  ! Factor: Sum ni*Zi**2 ~ 0.705*1**2 + 0.275*2**2 + 0.020*7**2 ~ 2.785
                  P_eps = Eps(i)*2.785d0*nmedio(l)*c*dsigma_Br
                  dE_e = E_e(j)*(Eint_e-1.d0)
                  dP_Eps = N_e(j)*P_eps*dE_e
                  sum = dP_Eps + sum ! P = int dP_Eps
               endif
            enddo ! j
            Lbr(i) = sum!*mec2         ! [erg/seg]
            Lbr_tot(i) = Lbr_tot(i) + Lbr(i)

            if (Ef(i) > 1.d8*eV .and. Ef(i) < Ef_max_gg) then
               iEf = int(log(Ef(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
               Lbr_tot_abs(i) = Lbr_tot_abs(i) + Lbr(i)*exp( -(Rst1/d1)*&
                    tau1(iEf,iang1) -(Rst2/d2)*tau2(iEf,iang2) )
            else
               Lbr_tot_abs(i) = Lbr_tot_abs(i) + Lbr(i)
            endif

            sum_br = sum_br + Lbr(i)*dEps
            
         enddo ! i
               
      ENDDO ! l
      
      write(*,*) 'L_Brem=',sum_br
           
!     Ya tengo la emisión brems y la emisión brems corregida por absorción g-g:
      
      open(16,file='radbr2.dat')
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*Lbr_tot(i),Ef(i)*Lbr_tot_abs(i)
      enddo
      close(16)
      
      deallocate ( x,y,z,rho,P,v,T,B,B_par )
      deallocate ( E_e,N_etot,N_e )
      deallocate ( Lbr,Lbr_tot,Lbr_tot_abs )
      deallocate ( tau1,tau2 )

    END subroutine CD_radbr2_run




!========================================


  END module CD_radbr


!========================================

