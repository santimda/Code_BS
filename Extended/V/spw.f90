program main
   implicit none
   real*8 :: pi, c, Mdot, vw, v_st, Rorb, rhop, p0, rw0
   real*8 :: x0, l, dl, x, y, psx0, psy0, th10, th20, sum1, sum2, sum3
   real*8 :: rst, rp, rhow, pw, th1, th2, dth1, dth2, theta_max
   real*8 :: dS1, dS2, psx, psy, ang, ps, L_ISM, dA
   real*8 :: v_t(2), v_p(2), v_vec(2)
   integer :: i

   open(1,file = 'kk.dat')
   open(2,file = 'L_theta.dat')
   open(3,file = 'v.dat')
   open(4,file = 'S.dat')
   pi = acos(-1.d0)
   c = 2.99792458d10 
   theta_max = (170.d0/180.d0)*pi ! ~ extension of the BS 

   Mdot = 1.d-8*1.989d33/31556926.d0
   vw = 1.5d8    

   v_st = 500.d5   ! v_star

   rhop = 10.d0 * 1.6726d-24 * 1.28d0 ! = n * m_p [cgs] * mu
   p0 = rhop*v_st**2
   rw0 = sqrt(Mdot*vw/4.d0/pi/p0) ! = R0

   print*, 'R_0 [pc]= ', rw0/3.08d18
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
   sum1 = 0.d0
   sum2 = 0.d0  
   sum3 = 0.d0  
   i = 1                
   do while (th10 <= theta_max)

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
      
      !psx = pw*cos(th1)*dS1 - pd*cos(th2)*dS2 + psx0
      !psy = pw*sin(th1)*dS1 + pd*sin(th2)*dS2 + psy0
      psx = pw*cos(th1)*dS1 - p0*dS2 + psx0
      psy = pw*sin(th1)*dS1 + psy0
      ang = atan2(psy,psx)      
      ps = sqrt(psx**2 + psy**2)

      x = x + psx*dl/ps
      y = y + psy*dl/ps
      l = l + dl

      psx0 = psx
      psy0 = psy

      v_t = (/ psx/ps, psy/ps /)
      v_p = (/ -psy/ps, psx/ps /)
      v_vec = v_st * (/ cos(th2),sin(th2) /)

      dA = 2.d0*pi*y*dl
      sum1 = sum1 + dA*sin(ang)*rhop*v_st**3/2.d0
      sum2 = sum2 + dA*sin(ang)**2*rhop*v_st**3/2.d0
      !sum2 = sum2 + dA*(DOT_PRODUCT(v_p,v_vec))**2*rhop*v_st/2.d0 --> This is equivalent
      sum3 = sum3 + dA*sin(ang)**3*rhop*v_st**3/2.d0

      write(1,*) x/rw0, y/rw0
      write(2,*) th1, sum1/L_ISM, sum2/L_ISM, sum3/L_ISM
      write(3,*) th1, abs(DOT_PRODUCT(v_t,v_vec)/v_st), abs(DOT_PRODUCT(v_p,v_vec)/v_st)
      write(4,*) th1, dS1/dl, dS2/dl, dA/dl
      i = i + 1
   enddo
   print*, 'i_max = ', i
   print*, 'L:', sum1/L_ISM, sum2/L_ISM, sum3/L_ISM

   close(1)
   close(2)
   close(3)
   close(4)
end

