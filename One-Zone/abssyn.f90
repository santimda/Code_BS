module abssyn
   private
   public :: abssyn_run

contains

   !========================================================================
   !        abssyn_run        
   !========================================================================
   !This module calculates the f-f absorption of radio emission from the shock
   ! This program calculates the total radio opacity for a radio photon 
   ! emitted from a position (x(l),y(l)) in the CD. The output is the file
   ! radsyn.dat which contains E(i), L_sy(E(i), L_sy_abs(E(i)).
   ! The program uses the results of f-f opacity (tauffw); we have an estimate 
   ! if the wind alpha(r,nu_i), with "r" tabulated in r.dat is the distance to 
   ! the star. The f-f absorption in the shock is negligible given the large T
   !========================================================================

   subroutine abssyn_run
      use global
      !$ use omp_lib
      Implicit none
      real(dp) :: alphaffw(mr_ff,mE_ff)
      Real(dp) :: r_ff(mr_ff),dr_ff,R_aprox
      Real(dp) :: theta1,x1,y1,z1,r1,dr0,x1_0,y1_0,z1_0,x,y,z
      Real(dp) :: dr_x,dr_y,sinc,cinc
      Real(dp) :: log_nu_int,flux,flux1,flux2
      Real(dp) :: emision1,emision2,emision1_abs,emision2_abs
      Real(dp) :: dE_nu,Eint_nu,sum
      Integer :: ir,inu_max,i_1G,i_5G,i_ff
      Logical :: integra

      allocate( tauff(m_nu) ) 
      allocate( Ef(m_nu),nu(m_nu),L_sy(m_nu),L_sy_abs(m_nu) )

      write(*,*) '\n Calculate the free-free synchrotron absorption'

      ! set omp threads number
      !$ call omp_set_num_threads(nth)  

      !============================================================
      ! Define some variables first and then read the opacity table

      tauff = 0.d0
      flux = 1.d26/(4.d0*pi*(distance)**2) ! luminosity to flux in mJy
      sinc = sin(incli)
      cinc = cos(incli)

      open(12,file='nu.dat')
      do i=1,m_nu
         read(12,*) nu(i)
      end do
      close(12)
      log_nu_int = log(nu(2)/nu(1))

      open(15,file='r_ff.dat')
      do i_ff=1,mr_ff
         read(15,*) r_ff(i_ff)
      end do
      close(15)
      dr_ff = r_ff(2)-r_ff(1)          ! equiespaciated

      open(18,file='windalphaff.dat')
      do j=1,mr_ff
         read(18,*) (alphaffw(j,i), i=1,m_nu)
      end do
      close(18)

      !============================================================
      ! We calculate the total f-f opacity, tauff, for photons coming from a 
      ! point in the CD integrating at distances d > r_shock along the LOS, and 
      ! for photons with nu < 1 GHz (~4.1e-5 eV). Note that the stellar wind  
      ! absorbs only when the photons are "inside" the CD ( r_xy < r(theta) )
      inu_max = int( log(1.d9/nu(1))/log_nu_int ) + 1

      ! Set the emitter position 
      dr0 = H0
      x = -R0
      y = 0.d0
      z = 0.d0

      x1_0 = x - dr0*sinc                 ! i = 0° = face-on, in that case dx=0
      y1_0 = y + dr0*cinc
      z1_0 = z                            ! it does not vary, it moves at cte z 
      dr_ff = dr0/64.d0                  ! integration step along the LOS
      dr_x = -dr_ff*sinc                  ! If i > 0°, x decreaces
      dr_y = dr_ff*cinc     

      ! $omp parallel do private(i,r,x1,y1,z1,ir,integra)
      do i=1,inu_max           ! Photon frequency
         x1 = x1_0  
         y1 = y1_0
         z1 = z1_0          

         integra = .true.
         do while(integra) 
            r1 = sqrt(x1**2 + y1**2 + z1**2)                     
            theta1 = abs(atan(y1/x1))
            R_aprox = R0 * sqrt(3.*(1.-theta1/tan(theta1)))/sin(theta1) ! R_BS(theta)
            if( r1 > r_ff(mr_ff) .or. r1 > R_aprox ) then   ! Stop integration too far
               integra = .false.
            else     ! "inside" the BS ( r1 <= R_aprox ) 
               ir = int((r1-r_ff(1))/dr_ff) + 1       
               integra = .true.  
               tauff(i) = tauff(i) + alphaffw(ir,i)*dr_ff 
            endif      
            x1 = x1 + dr_x ! CHECK
            y1 = y1 + dr_y          
         end do

      end do
      ! !$omp end parallel do

      open (25,file='tauff.dat')
      do i=1,m_nu
         write(25,*) nu(i),tauff(i) 
      end do
      close(25)

      !=====================================================================
      ! Correct synchrotron emission by the absorption

      open (31,file='radsyn.dat')
      do i=1,m_nu
         read(31,*) Ef(i),L_sy(i)  ! L_sy is E*L(E) [erg/s]; Ef(i) is in eV 
      end do
      close(31)

      i_1G = int(log(1.4d9/nu(1))/log_nu_int) + 1        ! 1.4 GHz
      i_5G = int(log(4.8d9/nu(1))/log_nu_int) + 1        ! 4.8 GHz
      flux1 = flux/nu(i_1G)   
      flux2 = flux/nu(i_5G) 

      !-----------------------------------------------------------------------
      ! Calculate the intrinsic vs absorbed emission 

      do i=1,m_nu
         if (i < inu_max) then !  The abs was calculated at these freqs  
            L_sy_abs(i) = L_sy(i)*exp(-tauff(i))  
         else
            L_sy_abs(i) = L_sy(i)
         endif
      end do

      emision1 = L_sy(i_1G) * flux1 ! At 1.4 GHz
      emision1_abs = emision1 * exp(-tauff(i_1G))
      emision2 = L_sy(i_5G) * flux2 ! At 4.8 GHz
      emision2_abs = emision2 * exp(-tauff(i_5G))
      write(*,*) 'Unabsorbed:'
      write(*,88) 'S(1.4GHz)=',emision1,'mJy'
      write(*,88) 'S(4.8GHz)=',emision2,'mJy'
      write(*,*) 'Absorbed:'
      write(*,88) 'S(1.4GHz)=',emision1_abs,'mJy'
      write(*,88) 'S(4.8GHz)=',emision2_abs,'mJy'

      88  format(A11,X,Es10.2,X,A3)
      89  format(A11,X,Es10.2)

      sum = 0.d0
      Eint_nu = nu(2)/nu(1)
      open(33,file='radsyn_abs.dat')
      do i=1,m_nu
         dE_nu = (h*nu(i))*(Eint_nu - 1.d0)  
         write(33,*) Ef(i),L_sy(i),L_sy_abs(i)
         sum = sum + L_sy_abs(i)/(Ef(i)*eV)*dE_nu
      end do
      close(33)
      write(*,89) 'L_syn(abs)=',sum

      deallocate( tauff ) 
      deallocate( Ef,nu,L_sy,L_sy_abs )


   END subroutine abssyn_run

END module abssyn

