      program main
      implicit real*8 (a-h,o-z)	
      open(1,file='kk.dat')
      pi=dacos(-1.d0)
      c=3.d10      
             
      zMdot=1.d-7*2.d33/3.1d7
      vw=2.d8
      
      bep=0.001d0
      vp=bep*c
      zLp=1.d50
      
      rhop=1.7d-24
      zp0=rhop*vp**2.d0
      Rorb=(2.d0*zLp/bep/c/4.d0/pi/zp0)**0.5d0
      rw0=(zMdot*vw/4.d0/pi/zp0)**0.5d0  
      print*, 'R_0 [pc] =  ', rw0/3.08d18
      
      zx0=rw0 
      zl=0.d0
      dl=rw0/100.d0 
      zx=zx0
      zy=dl
            
      psx0=0.d0
      psy0=1.d0
      th10=0.d0
      th20=0.d0      
      sum1=0.d0
      sum2=0.d0                  
      do i=1,1000

      zrst=(zx**2.d0+zy**2.d0)**0.5d0
      rp=((Rorb-zx)**2.d0+zy**2.d0)**0.5d0
      rhow=zMdot/4.d0/pi/zrst**2.d0/vw
      pw=rhow*vw**2.d0
      zpd=zp0*(Rorb/rp)**2.d0     

      th1=datan2(zy,zx)
      th2=datan2(zy,Rorb-zx)      
      dth1=th1-th10
      dth2=th2-th20
      th10=th1
      th20=th2

      dS1=2.d0*pi*zrst**2.d0*dsin(th1)*dth1
      dS2=2.d0*pi*rp**2.d0*dsin(th2)*dth2

      psx=pw*dcos(th1)*dS1-zpd*dcos(th2)*dS2+psx0
      psy=pw*dsin(th1)*dS1+zpd*dsin(th2)*dS2+psy0
      ang=datan2(psy,psx)      
      ps=(psx**2.d0+psy**2.d0)**0.5d0

      zx=zx+psx*dl/ps
      zy=zy+psy*dl/ps
      zl=zl+dl

      !write(1,*) zx,zy
      write(1,*) zx/rw0, zy/rw0

      psx0=psx
      psy0=psy

c      the=datan2(zy,zx)
c      zx1=dcos(the)*
c     #rw0*(1.d0/dsin(the))*(3.d0*(1.d0-the/dtan(the)))**0.5d0
c      zy1=dsin(the)*
c     #rw0*(1.d0/dsin(the))*(3.d0*(1.d0-the/dtan(the)))**0.5d0
c      write(3,*) zx1,zy1
c      write(4,*) zx1+cdh*dcos(th1),zy1+cdh*dsin(th1)
c      psx1=zx1-zx10
c      psy1=zy1-zy10
c      zx10=zx1
c      zy10=zy1
      
      sum1=sum1+2.d0*pi*zy*dl*dsin(ang)*rhop*vp**3.d0/2.d0
      sum2=sum2+2.d0*pi*zy*dl*dsin(ang)**2.d0*rhop*vp**3.d0/2.d0
      enddo
      
      print *, 'L:',
     #sum1/(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0),
     #sum2/(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0)
            
      close(1)
      end

