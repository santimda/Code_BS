      program main
      implicit real*8 (a-h,o-z)	
      open(1,file='kk1.dat')
      open(2,file='kk2.dat')
c      open(1,file='kk3.dat')
c      open(2,file='kk4.dat')
      pi=dacos(-1.d0)
      c=2.99792458d10 
             
      zMdot=1.d-8*1.989d33/31556926.d0
      vw=2.d8
      
      gad=5.d0/3.d0
      fga=gad/(gad-1.d0)      
      bep=0.001d0
      vp=bep*c
      gamp=1d0/(1.d0-bep**2.d0)**0.5d0
      zLp=1.d50
      zh0=gamp
      
      rhop=1.7d-24
      zp0=rhop*vp**2.d0
      zn0=fga*zp0/(zh0-1.d0)/c**2.d0
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
c      zx10=0.d0
c      zy10=0.d0
      zn=zn0
      cdh=0.d0
      ang=pi/2.d0
      zxm=zx
      zym=zy
      zMdotp=0.d0
                  
      sum=0.d0
      sum2=0.d0
      vs=0.d0
      do i=1,11250

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
      !print*, psx, psy, ang
      ps=(psx**2.d0+psy**2.d0)**0.5d0
      dx=psx*dl/ps
      dy=psy*dl/ps
      zx=zx+dx
      zy=zy+dy
      zl=zl+dl
      
      cs=(gad*zp/zn)**0.5d0
      
      if (vs.le.cs) then
      zp=zpd*dsin(ang)
      ang0=ang
      zn1=zn
      zp1=zp
      else
      zp=zpd*dsin(ang)**2.d0/dsin(ang0)
      !print*, zpd, ang, ang0
      endif
      zn=zn0*(zp/zp0)**(1.d0/gad)
      vs=2.d0**0.5d0*(0.5d0*vp**2.d0-fga*zp/zn)**0.5d0
      zMdotp=zMdotp+2.d0*pi*zy*dy*rhop*vp
      cdh=zMdotp/2.d0/pi/zy/zn/vs
      
      !print*, i, ang, p, vs, zMdotp, cdh

      zxp=zx+cdh*dcos(datan2(-psx,psy))
      zyp=zy+cdh*dsin(datan2(-psx,psy))
      dxp=zxp-zxm
      dyp=zyp-zym
      ang1=datan2(dyp,dxp)
      zxm=zxp
      zym=zyp

      psx0=psx
      psy0=psy
      
      write(1,*) zy,zx
      write(2,*) zyp,zxp
       
      sum=sum+2.d0*pi*dy*rhop*vp**3.d0*zy/2.d0
      sum2=sum2+2.d0*pi*dy*dsin(ang)*rhop*vp**3.d0*zy/2.d0
      enddo
      sum1=2.d0*pi*zy*cdh*fga*zp*vs
      
      print *, 
     #sum1/(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0),
     #sum2/(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0),
     #sum/(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0)

      print *, 0.5d0*pi*2.1d18**2.d0*rhop*vp**3.d0/
     #(0.5d0*pi*rw0**2.d0*rhop*vp**3.d0)
      print *, zy,zx,
     #(1.d-8*2.d33*2.d8/4.d0/pi/1.7d-24/3.d7**2.d0/3.2d7)**0.5d0


c      print *, cdh/3.08d18,zx/3.08d18,
c     #ang-ang1
c      print *, datan2((gad*zp1/zn1)**0.5d0,vp)
      
c      close(4)      
c      close(3)      
      close(2)      
      close(1)
      end

