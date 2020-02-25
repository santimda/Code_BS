module CD_abssyn
  !Este módulo calcula la absorción en radio de la fuente
  use global
  
 implicit none
  private
  public :: CD_abssyn_run

  contains
    
    !========================================================================
    !        abssyn_run        
    !========================================================================
  
    subroutine CD_abssyn_run
      use global

      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! Este programa calcula la opacidad en radio total para un fotón radio
      ! emitido en un punto (x(l),y(l)) de la CD. Genera como archivo de salida
      ! el archivo radsyn.dat en el que figura E(i), L_sy(E(i), L_sy_abs(E(i))
      ! El programa usa los resultados de la opacidad f-f del viento 1 (tauffw) 
      ! y parámetros de la región chocada. Para el viento, tengo una estimación 
      ! de alpha(r,nu_i), donde "r" está tabulada en r.dat y es la distancia a 
      ! la estrella 1 (pues el observador está a más de 90°). La absorción ff en
      ! la CD es despreciable pues T es muy alta allí.
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      !$ use omp_lib
      Implicit none
      Real(dp) :: dr1,dr2,r,dr,x1,y1,z1,phi,r0,dr0,x1_0,y1_0,z1_0
      Real(dp) :: dr_x,dr_y,tginc,sinc,cinc,a_CD
      Real(dp) :: Ef(mnu),log_nu_int,flux,flux1,flux2
      Real(dp) :: u_proj,v_proj ! coords del plano del cielo 
      Real(dp) :: emision1,emision2,emision1_abs,emision2_abs,sum,convert
      Integer :: ir,nl,inu_max,i_2G,i_8G
      Logical :: integra
      allocate ( tauff(ml,mnu) ) !, tauff1(ml,mnu),tauff2(ml,mnu) )
      allocate ( x(ml),y(ml),z(ml),nu(mnu),Lsy_tot(mnu) )
      tauff = 0.d0
      Lsy_tot = 0.d0
      Lsy_tot_abs = 0.d0
      convert = 206265.d3/(distance*kpc) ! imagen a tamaño angular (mas)
      flux = 1.d26/(4.d0*pi*(distance*kpc)**2) ! luminosidad a flujo en mJy

      tginc = tan(incli)
      sinc = sin(incli)
      cinc = cos(incli)
      ! set omp threads number
      !$ call omp_set_num_threads(nth)  

!     Leemos las posiciones y tablas de opacidad/absorción:

      open(11,file='CDrot.dat')
      do l=1,ml
         read(11,*) x(l),y(l),z(l)
      enddo
      close(11)

      open(12,file='nu.dat')
      do i=1,mnu
         read(12,*) nu(i)
      enddo
      close(12)
      log_nu_int = log(nu(2)/nu(1))

      open(15,file='r.dat')
      do i=1,mr_ff
         read(15,*) r_1(i),r_2(i)
      enddo
      close(15)
      dr1 = r_1(2)-r_1(1)         ! están equiespaciados
      dr2 = r_2(2)-r_2(1)

      open(18,file='windalphaff1.dat')
      open(19,file='windalphaff2.dat')
      do j=1,mr_ff
         read(18,*) (alphaffw1(j,i), i=1,mnu)
         read(19,*) (alphaffw2(j,i), i=1,mnu)
      enddo
      close(18)
      close(19)

      open(11,file='phi.dat')
      read(11,*) phi            ! Está en rad
      close(11)

      !/////////////////////////////////////////////////////////////////////
      !
      ! Ya estoy en condiciones de calcular la opacidad total tauff para 
      ! fotones emitidos desde cada punto de la CD. Para eso tengo en cuenta
      ! la opacidad en el shock por efecto R-T y la f-f en el viento, para lo 
      ! cual integro a d > r_shock siguiendo la línea de la visual.


      ! Calculo sólo la absorción a nu < 60 GHz ( ~ 8.2e-5 eV). Importante: no 
      ! en todos los puntos absorbe el viento de la estrella 2, por lo que me 
      ! fijo si el punto está "adentro" de la superficie (aproximada) de la CD:
      ! x - x_v = a*( (y-y_v)^2 + (z-z_v)**2 ), con un margen

      a_CD = 1.1d0*(x(ml/3-1)-x(1))/( (y(ml/3-1))**2 + z(ml/3-1)**2 ) 
!      write(*,*) 'a_CD [1/UA] =',a_CD*UA, x(1)/UA, y(1)/UA

      inu_max = int( log(500.d9*(closer/1.2d0)**2/nu(1))/log_nu_int ) + 1
!      write(*,*) 'inu_max=', inu_max

      do l=1,ml        ! Posición sobre la CD

         if (incli >= 0.d0) then
            r0 = sqrt(x(l)**2+y(l)**2+z(l)**2)     ! Dist a la estrella 1
         elseif (incli < 0.d0) then
            r0 = sqrt((D-x(l))**2+y(l)**2+z(l)**2) ! Dist a la estrella 2
         endif
         dr0 = 0.05d0*r0            ! estimación del ancho del shock
         x1_0 = x(l) - dr0*sinc    ! i = 0° = face-on, en ese caso dx=0
         y1_0 = y(l) - dr0*cinc
         z1_0 = z(l)               ! no varía, se mueve a z cte 
         dr = dr0/8.d0            ! Achico el paso un poco
         dr_x = dr*sinc
         dr_y = dr*cinc
         ! !$omp parallel do private(i,r,x1,y1,z1,ir,integra)
         do i=1,inu_max            ! Cada frecuencia de emisión       
            x1 = x1_0  
            y1 = y1_0
            z1 = z1_0
            integra = .true.
            do while(integra) 

               if( x1 >= 1.05_dp*x(1)+(y1**2+z1**2)*a_CD ) then  ! Absorbe st2
                  r = sqrt((x1-D)**2 + y1**2 + z1**2)  ! st2
                  ir = int((r-r_2(1))/dr2)             ! st2 adelante
                  if( r < r_2(mr_ff)) then
                     integra = .true.
                     tauff(l,i) = tauff(l,i) + alphaffw2(ir,i)*dr ! st2
                  else
                     integra = .false.
                  endif
               else                                    ! Absorbe st1
                  r = sqrt(x1**2 + y1**2 + z1**2)      ! st1
                  ir = int((r-r_1(1))/dr1)             ! st1 adelante 
                  if( r < r_1(mr_ff)) then
                     integra = .true.  
                     tauff(l,i) = tauff(l,i) + alphaffw1(ir,i)*dr ! st1
                  else
                     integra = .false.
                  endif
               endif
               x1 = x1 - dr_x
               y1 = y1 - dr_y

            enddo
         enddo
         ! !$omp end parallel do

      enddo

      open (25,file='CD_tauff.dat')
      do l=1,ml    
         do i=1,mnu
            write(25,*) l,nu(i),tauff(l,i) 
         enddo
      enddo
      close(25)

      open (31,file='emisyn1.dat')
      open (32,file='emisyn2.dat')
      do l=1,ml
         do i=1,mnu
            read(31,*) nl,Ef(i),Lsy1(l,i) ! Lsy ya está multiplicado por Ef(i)
            read(32,*) nl,Ef(i),Lsy2(l,i) ! Lsy ya está multiplicado por Ef(i)
        enddo
      enddo
      close(31)
      close(32)

      i_2G = int(log(2.3d9/nu(1))/log_nu_int) + 1 ! Bandwidth = 64 o 128 MHz?
      i_8G = int(log(8.6d9/nu(1))/log_nu_int) + 1 ! Bandwidth = 128 MHz
!      flux1 = 128.d6*flux/(nu(i_2G)**2.d0)   ! NO ESTOY SEGURO DEL BANDWIDTH
!      flux2 = 128.d6*flux/(nu(i_8G)**2.d0) 
      flux1 = flux/nu(i_2G)   ! NO ESTOY SEGURO DEL BANDWIDTH
      flux2 = flux/nu(i_8G) 
      ! Chequeo si caen cerca:
      !write(*,*) '2.3GHz aprox', nu(i_2G)/1.d9, '8.6GHz aprox', nu(i_8G)/1.d9
      !write(*,*) i_2G, i_8G

      ! Ahora voy a calcular la emisión total y absorbida para cada emisor 
      ! línea. También voy a guardar, aparte, cuál es la emisión proyectada
      ! en el plano del cielo (u_proj,v_proj) para armar mapas sintéticos
      open(70,file='CD_radiomap_2GHz.dat')
      open(71,file='CD_radiomap_8GHz.dat')
      do l=1,ml
!         u_proj = (y(l)+x(l)/tginc)*sinc *convert ! st1 ?
         u_proj = (y(l)+x(l)/abs(tginc))*abs(sinc) *convert ! st1 ?
         v_proj = z(l)*convert

         do i=1,mnu
            sum = Lsy1(l,i) + Lsy2(l,i)
            Lsy_tot(i) = Lsy_tot(i) + sum

            if (i < inu_max) then ! En estas frecuencias calculé la abs 
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + sum*exp(-tauff(l,i))  
            else
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + sum
            endif         
         enddo
         
         emision1 = (Lsy1(l,i_2G) + Lsy2(l,i_2G))*flux1 ! Para 2.3 GHz
         emision1_abs = emision1*exp(-tauff(l,i_2G))
         emision2 = (Lsy1(l,i_8G) + Lsy2(l,i_8G))*flux2 ! Para 8.6 GHz
         emision2_abs = emision2*exp(-tauff(l,i_8G))
         
         write(70,*) u_proj,v_proj,emision1,emision1_abs
         write(71,*) u_proj,v_proj,emision2,emision2_abs

      enddo
      close(70)
      close(71)

      open(33,file='radsynabs.dat')
      do i=1,mnu
         write(33,*) Ef(i),Lsy_tot(i),Lsy_tot_abs(i)
      enddo
     close(33)

      deallocate ( tauff ) !,tauff1,tauff2 )
      deallocate ( x,y,z,nu,Lsy_tot )
   

    END subroutine CD_abssyn_run
    
  END module CD_abssyn
  
