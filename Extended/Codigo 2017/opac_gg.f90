module opacidad_gg
  !Este módulo calcula la absorción gamma-gamma
  use global
  
  implicit none
  private
  public :: opacidad_gg_run

  contains

  !============================================================================
  ! CD_opacidad_run
  !============================================================================
  
    subroutine opacidad_gg_run
      use global

      implicit none
      real(dp) :: Eph0
      allocate ( tau1(mE_gg,m_ang_gg),tau2(mE_gg,m_ang_gg) )

      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !
      !                              ESTRELLA 1
      !
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      !  Vamos a tomar E_ph como los fotones objetivo de la estrella y Ef como
      !  los fotones emitidos por la fuente de rayos gamma
      Eph0 = 3.d0*k*T1    
      
      ! La subrutina calcula una matriz con los coeficientes de opacidad a la 
      ! distancia de un radio estelar del centro de la estrella para diferentes 
      ! ángulos entre el eje que une el emisor con la estrella y la dirección 
      ! del observador, y energías de los fotones gamma: tau(Ef,ángulo)

      write(*,*) 'Calculo la opacidad gamma-gamma de la estrella 1'

      call opacgg(tau1,Eph0,T1,Rst1)
      
      open (88,file='opac_gg1.dat')
      do i=1,mE_gg
         write(88,*) (tau1(i,j), j=1,m_ang_gg)
      enddo
      close(88)
      

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                              ESTRELLA 2
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     Vamos a tomar E_ph como los fotones objetivo de la estrella y Ef como
!     los fotones emitidos por la fuente de rayos gamma
      Eph0 = 3.d0*k*T2  

!     La subrutina calcula una matriz con los coeficientes de opacidad a la
!     distancia de un radio estelar del centro de la estrella para diferentes 
!     ángulos entre el eje que une el emisor con la estrella y la dirección 
!     del observador, y energías de los fotones gamma: tau(Ef,ángulo)
 
      write(*,*) 'Calculo la opacidad gamma-gamma de la estrella 2'

      call opacgg(tau2,Eph0,T2,Rst2)

      open (89,file='opac_gg2.dat')
      do i=1,mE_gg
         write(89,*) (tau2(i,j), j=1,m_ang_gg)
      enddo
      close(89)
            
      deallocate ( tau1,tau2 )
        
    END subroutine opacidad_gg_run
    
  END module opacidad_gg
  
  
  !/////////////////////////////////////////////////////////////////////////
  !
  !     SUBRUTINAS
  !
  !/////////////////////////////////////////////////////////////////////////
  
  ! Efmingg energía inicial fotón gamma
  ! Efintgg paso energía fotón gamma
  ! nfot número de pasos fotón gamma
  ! rstar radio estrella
  ! rorb radio de referencia para decidir donde dejar de integrar en l (del 
  ! orden de varias veces el tamaño del sistema)
  ! taugg opacidad matriz
  ! ee0 3kT
  
  subroutine opacgg(tau_gg,ee0,zT,rstar)
    use global
    implicit none
    integer, parameter :: ir_max=100
    integer :: iang,ie,ir,ie0
    real(dp) :: tau_gg(mE_gg,m_ang_gg),f(ir_max)
    real(dp) :: Cnumc2,tau,ee0,zT,rstar,Efemin,Efe,dang,ang,xi,yi,dr
    real(dp) :: angi, zmu,znu,Bb,Be,deint,e0,e0min,e0max,eint,siggg,sum
    real(dp) :: znbmc2,znst,xe

    Cnumc2 = 1.23d20
    
    ee0 = ee0/mec2
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/dfloat(mE_gg))
    Efemin = Ef_min_gg/mec2
    Efe = Efemin*Ef_int_gg**0.5d0

    dang = pi/m_ang_gg

    do ie=1,mE_gg
       ang = dang/2.d0

       do iang=1,m_ang_gg
          xi=rstar
          yi=0.d0
          dr=rstar/15.d0 ! R_final = Rstar + Rstar*ir_max/15
          tau=0.d0
 
          do ir=1,ir_max
             sum=0.d0
             xi=dr*cos(ang)+xi
             yi=dr*sin(ang)+yi
             zmu=(xi*cos(ang)+yi*sin(ang))/(xi**2.d0+yi**2.d0)**0.5d0
             if (zmu > 1.d0) then 
                zmu=1.d0
             elseif (zmu < -1.d0) then 
                zmu=-1.d0
             endif
             angi=acos(zmu)! ángulo de interacción a lo largo de la integral de camino; i.e. acos e_*e_gamma
             
             ! Integral para la densidad de fotones target             
             e0min=1.d0/Efe/0.5d0/(1.d0-cos(angi))
             e0max=30.d0*ee0
             eint=(e0max/e0min)**(1.d0/dfloat(30))
             if (e0min >= e0max) goto 0111
             e0=e0min*eint**0.5d0  
             
             do ie0=1,30
                deint=e0*(eint-1.d0)
                
                ! Densidad de fotones target
                
                znu=e0*Cnumc2
                Bb = (2.d0*h*znu**3.d0/c**2.d0)/(exp(h*znu/(k*zT))-1.d0)
                znbmc2=Bb*Cnumc2/(znu*h)
                znst=(znbmc2*pi/c)*(rstar**2/(xi**2+yi**2)) ! agregué dilución
                
                ! Sección eficaz gamma-gamma (Gould & Schreck 1967)
                
                xe=Efe*e0*0.5d0*(1.d0-cos(angi))
                Be=(1.d0-1/xe)**0.5d0
                siggg=(0.5d0*pi*(qe**2.d0/mec2)**2.d0)*(1.d0-Be**2.d0) &
                     *((3.d0-Be**4.d0)*log((1+Be)/(1-Be))-2.d0*Be*(2.d0-Be**2.d0))
                if (xe < 1.d0) siggg=0.d0
                
                ! Opacidad: tau(Efe,ang)
                sum=sum+siggg*znst*(1.d0-cos(angi))*deint
                ! tau=tau+siggg*znst*dr*(1.d0-dcos(angi))*deint
                
                e0=e0*eint
             enddo
0111         continue
             f(ir)=sum
          enddo
          call integrar(ir_max,dr,f,tau)
          tau_gg(ie,iang)=tau
          
          ang=ang+dang
       enddo
       
       Efe=Efe*Ef_int_gg
    enddo
    
    return
    
  end subroutine opacgg
  
  
