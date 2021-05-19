program main
   implicit none
   real*8 :: pi, c, Mdot, vw, bep, vp, Rorb, rhop, p0, rw0
   real*8 :: x0, l, dl, x, y, psx0, psy0, th10, th20, sum1, sum2, sum3
   real*8 :: rst, rp, rhow, pw, th1, th2, dth1, dth2
   real*8 :: dS1, dS2, psx, psy, ang, ps, L_ISM
   integer :: i

   open(1,file = 'kk.dat')
   open(2,file = 'L_theta.dat')
   pi = acos(-1.d0)
   c = 3.d10      

   Mdot = 1.d-7*2.d33/3.1d7
   vw = 2.d8

   bep = 0.001d0
   vp = bep*c

   rhop = 1.7d-24 ! = m_p [cgs]
   p0 = rhop*vp**2
   rw0 = sqrt(Mdot*vw/4.d0/pi/p0)

   Rorb = 1.d7*rw0 ! Valor exageradamente grande para simular viento plano
   
   ! Aproximación de L_ISM para un one-zone con R = R_0
   L_ISM = (0.5d0*pi*rw0**2*rhop*vp**3)

   x0 = rw0 
   l = 0.d0
   dl = rw0/100.d0 
   x = x0
   y = dl

   psx0 = 0.d0
   psy0 = 1.d0
   th10 = 0.d0
   th20 = 0.d0      
   sum1 = 0.d0
   sum2 = 0.d0  
   sum3 = 0.d0                  
   do i = 1,10000

      rst = sqrt( x**2 + y**2 )
      rhow = Mdot/4.d0/pi/rst**2/vw
      pw = rhow*vw**2 
      
      th1 = atan2(y,x) ! Este es el theta de R(theta)
      th2 = atan2(y,Rorb-x)   
      dth1 = th1-th10
      dth2 = th2-th20
      th10 = th1
      th20 = th2

      dS1 = 2.d0*pi*rst**2*sin(th1)*dth1
      dS2 = 2.d0*pi*Rorb**2*sin(th2)*dth2
      ! DUDA: COMO ESCRIBIR dS2

      !print*, dS2/(2.d0*pi*rst*dl*cos(th1))

      !print*, th2
      !psx = pw*cos(th1)*dS1 - pd*cos(th2)*dS2 + psx0
      !psy = pw*sin(th1)*dS1 + pd*sin(th2)*dS2 + psy0
      psx = pw*cos(th1)*dS1 - p0*dS2 + psx0
      psy = pw*sin(th1)*dS1 + psy0
      ang = atan2(psy,psx)      
      ps = sqrt(psx**2 + psy**2)

      x = x + psx*dl/ps
      y = y + psy*dl/ps
      l = l + dl

      !write(1,*) x,y
      write(1,*) x/rw0, y/rw0

      psx0 = psx
      psy0 = psy

      sum1 = sum1 + 2.d0*pi*y*dl*sin(ang)*rhop*vp**3/2.d0
      sum2 = sum2 + 2.d0*pi*y*dl*sin(ang)**2*rhop*vp**3/2.d0
      sum3 = sum3 + 2.d0*pi*y*dl*sin(ang)**3*rhop*vp**3/2.d0

      write(2,*) th1, sum1/L_ISM, sum2/L_ISM, sum3/L_ISM
   enddo

   print *, 'L:', sum1/L_ISM, sum2/L_ISM, sum3/L_ISM

   close(1)
   close(2)
end

