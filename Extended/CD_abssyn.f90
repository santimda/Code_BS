module CD_abssyn
  use global      
  implicit none
  public :: CD_abssyn_run

  contains
    
    !========================================================================
    !        abssyn_run        
    !========================================================================
  	 !This module calculates the f-f absorption of radio emission from the shock.
    ! This program calculates the total radio opacity for a radio photon 
    ! emitted from a position (x(l),y(l)) in the CD. The output is the file
    ! radsyn.dat which contains E(i), L_sy(E(i), L_sy_abs(E(i)).
    ! The program uses the results of f-f opacity (tauffw); we have an estimate 
    ! if the wind alpha(r,nu_i), with "r" tabulated in r.dat is the distance to 
    ! the star. The f-f absorption in the shock is negligible given the large T
    !========================================================================
  
    subroutine CD_abssyn_run(shock)
      !$ use omp_lib
!      real(dp), intent(in) :: phi
      character(LEN=*), intent(in) :: shock
      real(dp) :: dtheta, alphaffw(mr_ff,mE_ff)
!      real(dp) :: alpha_FS(mr_ff,mE_ff)
      Real(dp) :: r_ff(mr_ff),dr_ff
      Real(dp) :: theta1,x1,y1,z1,r1,dr0,x1_0,y1_0,z1_0
      Real(dp) :: dr_x,dr_y,tginc,sinc,cinc
      Real(dp) :: Ef(m_nu),log_nu_int,flux,flux1,flux2,flux3,flux4,flux5
      Real(dp) :: u_proj,v_proj  ! coords in the plane of the sky 
      Real(dp) :: emision1,emision2,emision3,emision1_abs,emision2_abs,emision3_abs,convert
      Real(dp) :: emision4,emision4_abs,emision5,emision5_abs
      Integer :: ir,nl,inu_max,i_1G,i_3G,i_5G,i_6G,i_10G,i_ff,i_theta
      integer :: N_avg
      integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
      real(dp) :: dnu_avg,h_avg
      real(dp) :: emi_tot(m_l),emi_avg(m_l),emi_tot_abs(m_l),emi_avg_abs(m_l)
      Logical :: integra
 
      allocate( tauff(m_l,m_nu) ) 
      allocate( x(m_l),y(m_l),z(m_l) )
      allocate( nu(m_nu),Lsy(m_l,m_nu),Lsy_tot(m_nu),Lsy_tot_abs(m_nu) )
      allocate( Lsy_RS(m_l,m_nu),Lsy_FS(m_l,m_nu) )
  
!100  format(A15,X,Es10.2)

      tauff = 0.d0
      Lsy_tot = 0.d0
      Lsy_tot_abs = 0.d0
      convert = 206265/distance             ! image to angular size ["]
      flux = 1.d26/(4.d0*pi*(distance)**2) ! luminosity to flux in mJy

!      write(*,100) 'R0_proj[arcmin] =',R0_proj/60

      tginc = tan(incli)
      sinc = sin(incli)
      cinc = cos(incli)
      dtheta = theta(2)-theta(1)  ! Equiesp.
      ! set omp threads number
      !$ call omp_set_num_threads(nth)  

      !======================================================
      !
      !    Read the opacity table
      
      open(11,file='CDrot.dat')
      do l=1,m_l     ! "up" and "down"
         read(11,*) x(l),y(l),z(l)
      end do
      close(11)
      x = R0*x
      y = R0*y
      z = R0*z

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

      !/////////////////////////////////////////////////////////////////////
      !
      ! We calculate the total f-f opacity, tauff, for photons coming from a 
      ! point in the CD integrating at distances d > r_shock along the LOS, and 
      ! for photons with nu < 2 GHz ( ~ 4.1e-5 eV). Note that the stellar wind  
      ! absorbs only when the photons are "inside" the CD ( r_xy < r(theta) )
  
      inu_max = int( log(0.5d9/nu(1))/log_nu_int ) + 1
!      write(*,*) 'inu_max=', inu_max

      do l=1,m_l        ! Position in the CD 

         if (l <= m_theta) then
            dr0 = 0.2d0*r_norm(l)*R0               ! Estimate of the shock width, USE H0 for consistency!!!!
         else
            dr0 = 0.2d0*r_norm(l-m_theta)*R0       ! Estimate of the shock width, USE H0 for consistency!!!!
         endif

         x1_0 = x(l) - dr0*sinc                 ! i = 0° = face-on, in that case dx=0
         y1_0 = y(l) + dr0*cinc
         z1_0 = z(l)                            ! it does not vary, it moves at cte z 
         dr_ff = dr0/16.d0                       ! integration step along the LOS
         dr_x = -dr_ff*sinc                      ! If i > 0°, x decreaces
         dr_y = dr_ff*cinc     

         ! !$omp parallel do private(i,r,x1,y1,z1,ir,integra)
         do i=1,inu_max           ! Photon frequency
            x1 = x1_0  
            y1 = y1_0
            z1 = z1_0          
            
            integra = .true.
            do while(integra) 
               r1 = sqrt(x1**2. + y1**2. + z1**2.)                     
               
               theta1 = abs(atan(y1/x1))
               i_theta = min( int((theta1 - theta(1))/dtheta) + 1 , m_theta)
               
               if( r1 > r_ff(mr_ff)) then        ! Stop integration too far
                  integra = .false.
               elseif ( r1 < r_norm(i_theta)*R0 .or. i_theta == m_theta) then ! "inside" the BS
                  ir = int((r1-r_ff(1))/dr_ff) + 1       
                  integra = .true.  
                  tauff(l,i) = tauff(l,i) + alphaffw(ir,i)*dr_ff 
               else
                  integra = .false.
               endif
               
               x1 = x1 + dr_x ! CHECK
               y1 = y1 + dr_y
               
            end do
         end do
         ! !$omp end parallel do

      end do

      open (25,file='CD_tauff.dat')
      do l=1,m_l
         do i=1,m_nu
            write(25,*) l,nu(i),tauff(l,i) 
         end do
      end do
      close(25)

      !/////////////////////////////////////////////////////////////////////

      !open (31,file='emisyn_'//shock//'.dat')
      if (shock == 'RS') then
         open (32,file='emisyn_RS.dat')
         do l=1,m_theta
            do i=1,m_nu
               read(32,*) nl,Ef(i),Lsy(l,i)  ! Lsy is already multiplied by Ef(i)
               Lsy(l+m_theta,i) = Lsy(l,i)   ! symmetry 
         end do
         end do
         close(32)
      
      else if (shock == 'FS') then
         open (33,file='emisyn_FS.dat')
         do l=1,m_theta
            do i=1,m_nu
               read(33,*) nl,Ef(i),Lsy(l,i)  ! Lsy is already multiplied by Ef(i)
               Lsy(l+m_theta,i) = Lsy(l,i)   ! symmetry 
         end do
         end do
         close(33)

      else if (shock == 'BOTH') then
         open (32,file='emisyn_RS.dat')
         !open (33,file='emisyn_FS.dat')
         do l=1,m_theta
            do i=1,m_nu
               read(32,*) nl,Ef(i),Lsy_RS(l,i)  ! Lsy is already multiplied by Ef(i)
               !read(33,*) nl,Ef(i),Lsy_FS(l,i)
               Lsy(l,i) = Lsy_RS(l,i)! + Lsy_FS(l,i)
               Lsy(l+m_theta,i) = Lsy(l,i)   ! symmetry 
         end do
         end do
         close(32)
         !close(33)
      end if

      i_1G = int(log(1.42d9/nu(1))/log_nu_int) + 1        ! 1.4 GHz
      i_3G = int(log(3.00d9/nu(1))/log_nu_int) + 1         ! 3.0 GHz
      i_5G = int(log(4.86d9/nu(1))/log_nu_int) + 1        ! 4.8 GHz
      i_6G = int(log(6.d9/nu(1))/log_nu_int) + 1        ! 6.0 GHz
      i_10G = int(log(10.d9/nu(1))/log_nu_int) + 1        ! 10. GHz
      flux1 = flux/nu(i_1G)   
      flux2 = flux/nu(i_3G) 
      flux3 = flux/nu(i_5G)
      flux4 = flux/nu(i_6G)
      flux5 = flux/nu(i_10G)

      !-----------------------------------------------------------------------
      ! Calculate the total and absorbed emission for each linear emitter. 
      ! Also stores the emission projected in the plane of the sky, (u_proj,v_proj)

      if (shock /= 'BOTH') then
         open(70,file='CD_sync_radiomap_'//shock//'1GHz.dat') ! 1.4 GHz
         open(71,file='CD_sync_radiomap_'//shock//'3GHz.dat') ! 3.0 GHz
         open(72,file='CD_sync_radiomap_'//shock//'5GHz.dat') ! 4.8 GHz
         open(74,file='CD_sync_radiomap_'//shock//'6GHz.dat') ! 6.0 GHz
         open(75,file='CD_sync_radiomap_'//shock//'10GHz.dat') ! 10. GHz
      else if (shock == 'BOTH') then
         open(70,file='CD_sync_radiomap_1GHz.dat') ! 1.4 GHz
         open(71,file='CD_sync_radiomap_3GHz.dat') ! 3.0 GHz
         open(72,file='CD_sync_radiomap_5GHz.dat') ! 4.8 GHz
         open(74,file='CD_sync_radiomap_6GHz.dat') ! 6.0 GHz
         open(75,file='CD_sync_radiomap_10GHz.dat') ! 10. GHz
      end if

      do l=1,m_l
!         u_proj = (y(l)+x(l)/tginc)*sinc *convert ! ?
         u_proj = (y(l)+x(l)/abs(tginc))*abs(sinc) *convert 
         v_proj = z(l)*convert

         do i=1,m_nu
            Lsy_tot(i) = Lsy_tot(i) + Lsy(l,i)
            if (i < inu_max) then !  The abs was calculated at these freqs  
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy(l,i)*exp(-tauff(l,i))  
            else
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy(l,i)
            endif

         end do
         
         emision1 = Lsy(l,i_1G) * flux1 ! At 1.4 GHz
         emision1_abs = emision1 * exp(-tauff(l,i_1G))
         emision2 = Lsy(l,i_3G) * flux2 ! At 3.0 GHz
         emision2_abs = emision2 * exp(-tauff(l,i_3G))
         emision3 = Lsy(l,i_5G) * flux3 ! At 4.8 GHz
         emision3_abs = emision3 * exp(-tauff(l,i_5G))
         emision4 = Lsy(l,i_6G) * flux4 ! At 6.0 GHz
         emision4_abs = emision4 * exp(-tauff(l,i_6G))
         emision5 = Lsy(l,i_10G) * flux5 ! At 10. GHz
         emision5_abs = emision5 * exp(-tauff(l,i_10G))
         
         write(70,*) u_proj,v_proj,emision1,emision1_abs
         write(71,*) u_proj,v_proj,emision2,emision2_abs
         write(72,*) u_proj,v_proj,emision3,emision3_abs
         write(74,*) u_proj,v_proj,emision4,emision4_abs
         write(75,*) u_proj,v_proj,emision5,emision5_abs
      end do
      close(70)
      close(71)
      close(72)
      close(74)
      close(75)

      !------------------------------------------------------------
      ! Avg flux btw nu_av_min and nu_av_max
      if (shock /= 'BOTH') then
         open(73,file='CD_sync_radiomap_avg'//shock//'.dat')
      else if (shock == 'BOTH') then
         open(73,file='CD_sync_radiomap_avg.dat')
      end if

      ! Count the number of steps
      N_avg = 0
      do i=1,m_nu
         if ((nu(i) >= nu_av_min) .and. (nu(i) <= nu_av_max)) then
            N_avg = N_avg + 1
         end if
      end do

      h_avg = (nu_av_max/nu_av_min)**(1.d0/(N_avg-1))

      do l=1,m_l
         u_proj = (y(l)+x(l)/abs(tginc))*abs(sinc) *convert 
         v_proj = z(l)*convert
         emi_tot(l) = 0.d0
         emi_tot_abs(l) = 0.d0
         do i=1,m_nu
            if ((nu(i) >= nu_av_min) .and. (nu(i) <= nu_av_max)) then
               dnu_avg = nu(i)*(h_avg-1.d0)
               emi_tot(l) = emi_tot(l) + Lsy(l,i)*flux/nu(i) * dnu_avg
               emi_tot_abs(l) = emi_tot_abs(l) + Lsy(l,i)*exp(-tauff(l,i))*flux/nu(i) * dnu_avg
            end if
         end do
         emi_avg(l) = emi_tot(l)/(nu_av_max-nu_av_min)
         emi_avg_abs(l) = emi_tot_abs(l)/(nu_av_max-nu_av_min)
         write(73,*) u_proj,v_proj,emi_avg(l),emi_avg_abs(l)
      end do
      close(73)

      !------------------------------------------------------------
      if (shock /= 'BOTH') then
         open(34,file='radsyn_abs_'//shock//'.dat')
      else if (shock == 'BOTH') then
         open(34,file='radsyn_abs.dat')
      end if
      do i=1,m_nu
         write(34,*) Ef(i),Lsy_tot(i),Lsy_tot_abs(i)
      end do
     close(34)

      deallocate( tauff ) 
      deallocate( x,y,z,nu,Lsy,Lsy_tot,Lsy_tot_abs,Lsy_RS,Lsy_FS )
   
    END subroutine CD_abssyn_run
    
  END module CD_abssyn
  
