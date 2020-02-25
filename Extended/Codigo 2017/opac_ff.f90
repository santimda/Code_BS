module opacidad_ff
  !Este módulo calcula la absorción free-free en el viento
  use global
  
  implicit none
  private
  public :: opacidad_ff_run

  contains

  !============================================================================
  ! CD_opacidad_ff_run
  !============================================================================
  
    subroutine opacidad_ff_run
      use global   
      !$ use omp_lib
      implicit none

      !/////////////////////////////////////////////////////////////////////////
      !     Este programa calcula los coeficientes de la emisividad
      !     y la absorción free-free en radio para el viento de las estrellas.
      !     Las salidas de este programa son:
      !     1) windalphaff1.dat y windalphaff2.dat, en los que se da tabla con: 
      !     nu(Hz) y la opacidad alphaff_i(r_i,nu) [1/cm], con r_i la distancia 
      !     a la estrella i.
      !     2) windemiff1.dat y windemiff2.dat, en los que se da tabla con:
      !     nu(Hz) y el coef de emisividad epsff(l,nu)*dl^3
      !     3) r.dat, en el que se da una lista de r_1(n), r_2(n)
      !/////////////////////////////////////////////////////////////////////////

      real(dp) :: Ry,T_w,S(mE_ff),S1(mE_ff),S2(mE_ff),ne,ni,gffinter
      real(dp) :: alphaff(mr_ff,mE_ff),gamma2,dg,dmu
      real(dp) :: x_apex,y_apex,r_1max,r_2max,dr1,dr2,A1,A2
      integer :: j_mu,j_g
      integer, parameter :: nu_mu = 146,nu_mga=81 ! # datos en tabla gauntff.dat
      real(dp) :: gff(nu_mu,nu_mga),mu_ff(mE_ff)
      allocate ( nu(mE_ff) )
      S1 = 0.d0
      S2 = 0.d0
      S = 0.d0
      Ry = 2.17987d-11          ! Rydberg unit = 2.18e-18 J

      ! set omp threads number
      !$ call omp_set_num_threads(nth)
      ! !$ write(*,'(a, i3, a)') '# Using OpenMP with ', nth ,' threads'
      ! !$ write(*,*)
      
!///////////////////////////////////////////////////////////////////////////
!
!     Considero vientos con simetría esférica, así que para calcular tauff
!     y emiff del wind lo hago en la coord r

 
      open(11,file='CD.dat')
      read(11,*) x_apex,y_apex  
      r_1max = 6.d0*x_apex         ! Radio máx desde la * 1 = 6 * dist1 al apex
      r_2max = 2.d0/eta*(D-x_apex) ! Radio máx desde la * 2 ~ 6 * dist2 al apex
      close(11)

      r_1(1) = Rst1
      r_1(mr_ff) = r_1max
      dr1 = (r_1max - Rst1)/float(mr_ff)
      r_2(1) = Rst2
      r_2(mr_ff) = r_2max
      dr2 = (r_2max - Rst2)/float(mr_ff)
      open (12,file='r.dat')
      write(12,*) r_1(1),r_2(1)
      do l = 2,mr_ff
         r_1(l) = r_1(l-1) + dr1
         r_2(l) = r_2(l-1) + dr2
         write(12,*) r_1(l),r_2(l)
      enddo
      close(12)


!///////////////////////////////////////////////////////////////////////////
!
!     Leemos la tabla de van Hoof 2014 para estimar el valor del Gaunt factor:
!     <gff(u,gamma2)>, con -8 < log10u < 8, -4 < log10gamma < 4, todos a paso 1
!
      
      open(13,file='gauntff.dat')
      do i=1,16
         read(13,*)
      enddo
      
      do i=1,nu_mu
         read(13,*) (gff(i,j), j=1,nu_mga)
      enddo
      close(13)


!///////////////////////////////////////////////////////////////////////////
!     
!     Procedemos a calcular los coeficientes tauff(r,E) para la *1
!
!///////////////////////////////////////////////////////////////////////////

      T_w = 0.4d0*T1    ! tomo T = cte para el viento, pues estoy lejos de la * 
      open(14,file='nu.dat')
      do i=1,mE_ff
         read(14,*) nu(i)
         mu_ff(i) = h*nu(i)/k/T_w   ! total tomo T = cte para el viento
      enddo
      close(14)
      
!     Calculamos los coeficientes para el lado 1 (de Falceta 2012):

!     ESTO NO ME CIERRA
      A1 = Mdot1/(4.d0*pi*vinf1*mu_w1*mp)
      gamma2 = Zq_w1**2*Ry/k/T_w
      j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
      dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

      !$omp parallel do &
      !$omp private(l,ni,ne,i,j_mu,dmu,gffinter) 
      do l=1,mr_ff
         ni = A1/r_1(l)**2/(gamma_w1+1.d0)  ! mui = mu*(g_e+1)
         ne = gamma_w1*ni           
         do i=1,mE_ff
       
            ! Para gff: i = (valor - valor inicial tabla) / (paso tabla) + 1
            ! Luego de obtener el punto, interpolo con la aproximación al plano tangente
            j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
            dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
            call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
!            gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w1/nu(i))) ! Una aprox
            alphaff(l,i) = 3.7d8*Zq_w1**2*ne*ni*T_w**(-0.5_dp)*nu(i)**(-3.d0)*&
                 (1.d0-exp(-h*nu(i)/k/T1))*gffinter !1/cm
         enddo
      enddo
      !$omp end parallel do

      open(18,file='windalphaff1.dat')
      do l=1,mr_ff
         write(18,*) (alphaff(l,i), i=1,mE_ff)
      enddo
      close(18)

      ! Para el flujo en radio uso la ec. 8 de Wright 1975:
      open(19,file='windemiff1.dat')

      i=1
      do while (i < mE_ff .and. nu(i) < 1.d-1*k*T_w/h )! la fórmula vale para h*nu << k*T 
         j_mu = int((log10(mu_ff(i))+16.d0)/0.2d0) + 1
         dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2
         call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)    
!         gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w1/nu(i))) ! Una aprox          
         S1(i) = (23.2d0/distance**2)* ( ( (Mdot1/sqrt(filling)/Msolyr)/&
              (vinf1/1.d5)/mu_w1)**2.d0 *nu(i)*gamma_w1*&
              gffinter*Zq_w1**2.d0 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
         if (nu(i) > 1.d-3*k*T_w/h) then
            S1(i) = S1(i)*exp(-h*nu(i)/k/T_w) ! Lo inventé yo 
         endif      

         write(19,*) nu(i),S1(i)*1.d3 ! Lo guardo en mJy
         i=i+1
      enddo
      close(19)


!///////////////////////////////////////////////////////////////////////////
!     
!     Procedemos a calcular los coeficientes tauff(r,E) para la *2
!
!///////////////////////////////////////////////////////////////////////////

      T_w = 0.4d0*T2
      
      do i=1,mE_ff
         mu_ff(i) = h*nu(i)/k/T_w   ! total tomo T = cte para el viento
      enddo
      
!     Calculamos los coeficientes para el lado 2 (de Falceta 2012):

      
!     ESTO NO ME CIERRA
      A2 = Mdot2/(4.d0*pi*vinf2*mu_w2*mp)
      gamma2 = Zq_w2**2*Ry/k/T_w
      j_g = int((log10(gamma2)+6.d0)/0.2d0) + 1
      dg = ((log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g))*0.2

      !$omp parallel do &
      !$omp private(l,ni,ne,i,j_mu,dmu,gffinter) 
      do l=1,mr_ff
         ni = A2/r_2(l)**2/(gamma_w2+1.d0) 
         ne = gamma_w2*ni   
         do i=1,mE_ff
       
            ! Para gff: i = (valor - valor inicial tabla) / (paso tabla) + 1
            ! Luego de obtener el punto, interpolo con la aproximación al plano tangente
            j_mu = int((log10(mu_ff(i))+16.d0)/0.2d0) + 1
            dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
            call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
            !gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w2/nu(i))) ! Una aprox util
            alphaff(l,i) = 3.7d8*Zq_w1**2*ne*ni*T_w**(-0.5_dp)*nu(i)**(-3.d0)*&
                 (1.d0-exp(-h*nu(i)/k/T_w))*gffinter !1/cm

         enddo
      enddo
      !$omp end parallel do

     open(18,file='windalphaff2.dat')
      do l=1,mr_ff
         write(18,*) (alphaff(l,i), i=1,mE_ff)
      enddo
      close(18)

!     Para el flujo en radio uso la ec. 8 de Wright 1975:      
      open(19,file='windemiff2.dat')
      i=1
      do while (i < mE_ff .and. nu(i) < 1.d-1*k*T_w/h )! la fórmula vale para h*nu << k*T 
         j_mu = int((log10(mu_ff(i))+16.d0)/0.2d0) + 1
         dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2d0
         call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
!         gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w2/nu(i)))
         S2(i)= (23.2d0/distance**2)* ( ( (Mdot2/sqrt(filling)/Msolyr)/&
              (vinf2/1.d5)/mu_w2)**2.d0 *nu(i)*gamma_w2*&
              gffinter*Zq_w2**2.d0 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
         if (nu(i) > 1.d-3*k*T_w/h) then
            S2(i) = S2(i)*exp(-h*nu(i)/k/T_w) ! Lo inventé yo 
         endif

         write(19,*) nu(i),S2(i)*1.d3 ! Lo guardo en mJy
         i=i+1
      enddo
      close(19)

      open(19,file='windemiff.dat')
      do i = 1, mE_ff
         S(i) = S1(i) + S2(i)
         write(19,*) nu(i),S(i)*1.d3 ! Lo guardo en mJy
      enddo
      close(19)

      deallocate ( nu )

    END subroutine opacidad_ff_run

  END module opacidad_ff

  !////////////////////////////////////////////////////////////////////////////////
  !

  Subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
      real*8 dx,dy,fx1y1,fx2y1,fx1y2,finterpol,fx,fy
      fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
      fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
      finterpol = fx1y1 + fx*dx + fy*dy
      return
    end Subroutine interpol
    
