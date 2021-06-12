program main
   implicit none
   real*8 :: pi, c, Mdot, vw, v_st, Rorb, rhop, p0, rw0
   real*8 :: x0, l, dl, x, y, psx0, psy0, th10, th20, sum1, sum2, sum3
   real*8 :: rst, rp, rhow, pw, th1, th2, dth1, dth2, theta_max
   real*8 :: dS1, dS2, psx, psy, ang, ps, L_ISM, dA
   real*8 :: gamp, h0, gad, fga, ang0, ang1, cdh, cs, dx, dxp, dy, dyp
   real*8 :: Mdotp, n, n0, n1, p, p1, pd, sum, vs, xm, ym, xp, yp
   integer :: i

   open(1,file = 'kk1.dat')
   open(2,file = 'kk2.dat')

   pi = acos(-1.d0)
   c = 2.99792458d10 
   theta_max = (160.d0/180.d0)*pi ! ~ extension of the BS 

   Mdot = 1.d-8*1.989d33/31556926.d0 ! Stellar mass-loss rate
   vw = 2.d8                         ! Stellar wind velocity

   gad = 5.d0/3.d0
   fga = gad/(gad-1.d0)      
   v_st = 500.d5   ! v_star
   gamp = 1d0/sqrt(1.d0-(v_st/c)**2)
   h0 = gamp        ! enthalpy ?
   n0 = fga*p0/(h0-1.d0)/c**2 ! density at the apex?
   rhop = 10.d0 * 1.6726d-24 * 1.28d0 ! = n * m_p [cgs] * mu
   p0 = rhop*v_st**2
   n0 = fga*p0/(h0-1.d0)/c**2
   rw0 = sqrt(Mdot*vw/4.d0/pi/p0) ! = R0

   print*, 'R_0 [pc] =  ', rw0/3.08d18
   Rorb = 1.d7*rw0 ! Valor exageradamente grande para simular viento plano
      
   ! Aproximación de L_ISM para un one-zone con R = R_0
   L_ISM = (0.5d0*pi*rw0**2*rhop*v_st**3)

   x0 = rw0 
   l = 0.d0
   dl = rw0/100.d0 
   x = x0
   y = dl

   psx0 = 0.d0
   psy0 = 1.d0
   th10 = 0.d0
   th20 = 0.d0      
   n = n0
   cdh = 0.d0
   ang = pi/2.d0
   xm = x
   ym = y
   Mdotp = 0.d0

   sum = 0.d0
   sum2 = 0.d0
   vs = 0.d0
   p = 0
   i = 1                
   do while (th10 <= theta_max)

      rst = sqrt( x**2 + y**2 )
      rhow = Mdot/4.d0/pi/rst**2/vw
      pw = rhow*vw**2 

      rp = sqrt((Rorb-x)**2+y**2) ! rp = Rorb? 
      pd = p0*(Rorb/rp)**2     

      th1 = atan2(y,x)  ! Este es el theta de R(theta)
      th2 = atan2(y,Rorb-x)      
      dth1 = th1-th10
      dth2 = th2-th20
      th10 = th1
      th20 = th2

      dS1 = 2.d0*pi*rst**2*sin(th1)*dth1
      dS2 = 2.d0*pi*rp**2*sin(th2)*dth2 

      psx = pw*cos(th1)*dS1 - p0*dS2 + psx0
      psy = pw*sin(th1)*dS1 + psy0
      ang = atan2(psy,psx)
      ps = sqrt(psx**2 + psy**2)

      dx = psx*dl/ps
      dy = psy*dl/ps
      x = x + dx
      y = y + dy
      l = l + dl

      cs = sqrt(gad*p/n)
      !print*, i, p, vs, cs
      if (vs <= cs) then
         p = pd*sin(ang)
         ang0 = ang
         n1 = n
         p1 = p
      else
         p = pd*sin(ang)**2/sin(ang0)
      endif
      n = n0*(p/p0)**(1.d0/gad)
      vs = sqrt( 2.d0 * (0.5d0*v_st**2-fga*p/n) )
      Mdotp = Mdotp + 2.d0*pi*y*dy*rhop*v_st
      cdh = Mdotp/2.d0/pi/y/n/vs

      xp = x+cdh*cos(atan2(-psx,psy))
      yp = y+cdh*sin(atan2(-psx,psy))
      dxp = xp-xm
      dyp = yp-ym
      ang1 = atan2(dyp,dxp)
      xm = xp
      ym = yp

      psx0 = psx
      psy0 = psy

      write(1,*) x/rw0, y/rw0
      write(2,*) yp,xp

      sum = sum + 2.d0*pi*dy*rhop*v_st**3*y/2.d0
      sum2 = sum2 + 2.d0*pi*dy*sin(ang)*rhop*v_st**3*y/2.d0
      i = i + 1
   enddo
   print*, 'i_max = ', i
   sum1 = 2.d0*pi*y*cdh*fga*p*vs

   ! NO SE QUE ES LO QUE IMPRIME ACA...
   print*, sum1/(0.5d0*pi*rw0**2*rhop*v_st**3)
   print*, sum2/(0.5d0*pi*rw0**2*rhop*v_st**3)
   print*, sum/(0.5d0*pi*rw0**2*rhop*v_st**3)

   print *, (0.5d0*pi*2.1d18**2*rhop*v_st**3)/ (0.5d0*pi*rw0**2*rhop*v_st**3)
   print*, y,x
   print*, sqrt(1.d-8*2.d33*2.d8/4.d0/pi/1.7d-24/3.d7**2/3.2d7)

   close(2)      
   close(1)
end

