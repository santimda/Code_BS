module CD_abssyn
   use global
   use gaunt_factor  
   use map_freqs, only: maps, nu_index
   implicit none
   public :: CD_abssyn_run

   contains
    
    !========================================================================
    !        abssyn_run        
    !========================================================================
    ! This module calculates the f-f absorption of radio emission from the shock.
    ! This program calculates the total radio opacity for a radio photon 
    ! emitted from a position (x(l),y(l)) in the CD. The output is the file
    ! radsyn.dat which contains E(i), L_sy(E(i), L_sy_abs(E(i)).
    ! The program uses the results of f-f opacity (tauffw); we have an estimate 
    ! if the wind alpha(r,nu_i), with "r" tabulated in r.dat is the distance to 
    ! the star. The f-f absorption in the shock is negligible given the large T
    !========================================================================

    subroutine CD_abssyn_run(shock,radiative)
      !$ use omp_lib
      character(LEN=*), intent(in) :: shock
      logical, intent(in) :: radiative
      real(dp) :: dtheta
      Real(dp) :: r_ff(mr_ff),dr_ff
      Real(dp) :: theta1,x1,y1,z1,r1,dr0,x1_0,y1_0,z1_0
      Real(dp) :: dr_x,dr_y,tginc,sinc,cinc
      Real(dp) :: Ef(m_nu),log_nu_int,flux,convert
      Real(dp) :: u_proj,v_proj  ! coords in the plane of the sky 
      real(dp) :: dnu_avg,h_avg
      real(dp) :: R_str
      Integer :: ir,nl,inu_max,i_ff,i_theta
      ! radio maps at the frequencies of maps(:) (map_freqs.f90)
      integer :: im, n_map, iu
      integer,  allocatable :: i_nu(:)               ! index of each map frequency in nu(:)
      real(dp), allocatable :: u_map(:), v_map(:)    ! coords in the plane of the sky
      real(dp), allocatable :: emis(:,:), emis_abs(:,:)  ! (position, map frequency)
      character(len=:), allocatable :: tag
      integer :: N_avg
      integer, parameter :: m_l = 2*m_theta  ! CD position, both sides
      Logical :: integra, inside

      real(dp), allocatable :: emi_tot(:), emi_avg(:), emi_tot_abs(:), emi_avg_abs(:)
      real(dp), allocatable :: H_theta(:), H_IL(:), rho(:), n(:)
      real(dp), allocatable :: alpha_IL(:,:), alpha_str(:), alphaffw(:,:)
      real(dp), allocatable :: gaunt(:,:), gaunt_str(:)

      allocate( emi_tot(m_l),emi_avg(m_l),emi_tot_abs(m_l),emi_avg_abs(m_l) )
      allocate( H_theta(m_l), H_IL(m_theta), rho(m_theta), n(m_l) )
      allocate( alpha_IL(m_l,m_nu), alpha_str(m_nu), alphaffw(mr_ff,mE_ff) )
      allocate( tauff(m_l,m_nu),tau_IL(m_l,m_nu),tau_str(m_nu),tau_tot(m_l,m_nu) ) 
      allocate( x(m_l),y(m_l),z(m_l) )
      allocate( nu(m_nu),Lsy(m_l,m_nu),Lsy_tot(m_nu),Lsy_tot_abs(m_nu) )
      allocate( Lsy_RS(m_l,m_nu),Lsy_FS(m_l,m_nu) )
      allocate( gaunt(m_l,m_nu),gaunt_str(m_nu) )

      tauff = 0.d0
      Lsy_tot = 0.d0
      Lsy_tot_abs = 0.d0
      convert = 206265/distance             ! image to angular size ["]
      flux = 1.d26/(4.d0*pi*(distance)**2) ! luminosity to flux in mJy

      tginc = tan(incli)
      sinc = sin(incli)
      cinc = cos(incli)
      dtheta = theta(2)-theta(1)  ! Equiesp.
      ! set omp threads number
      !$ call omp_set_num_threads(nth)

      !======================================================
      ! Read data of the IL
      if ((shock == 'BOTH') .and. radiative) then
         open(20,file='termo_FS.dat')
         open(21,file='H_FS.dat')
         do l=1, m_theta
            read(20,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,rho(l)
            read(21,*) dummy, H_IL(l), dummy, dummy
            n(l) = rho(l)/mu_SISM/mp
            n(l+m_theta) = rho(l)/mu_SISM/mp
            H_theta(l) = H_IL(l) 
            H_theta(l+m_theta) = H_IL(l)
         end do
         close(20)
         close(21)
      end if

      !======================================================
      !
      ! Read the opacity table
      
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

      !======================================================
      do l=1,m_l
         do i=1,m_nu
            call gaunt_ff_calc(nu,m_nu,T_ISM,Zq_ISM,gaunt(l,i))
         end do
      end do

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

            inside = .true.
            do while(inside)
               r1 = sqrt(x1**2. + y1**2. + z1**2.)                     
               
               theta1 = abs(atan(y1/x1))
               i_theta = min( int((theta1 - theta(1))/dtheta) + 1 , m_theta)
               
               if ( r1 < r_norm(i_theta)*R0 .or. i_theta == m_theta) then ! "inside" the BS
                  alpha_IL(l,i) = 3.7d8 * (Zq_ISM**2*n(i_theta)**2/sqrt(T_ISM)) * (1.d0-exp(-h*nu(i)/(k*T_ISM)) )/nu(i)**3 * gaunt(i_theta,i)
                  tau_IL(l,i) = alpha_IL(l,i)*H_theta(i_theta)
               else
                  inside = .false.
               endif
               
               x1 = x1 + dr_x ! CHECK
               y1 = y1 + dr_y
            end do
         end do
         ! !$omp end parallel do

      end do

      !==============================================
      open(30,file='Stromgren.dat')
         read(30,*) R_str
      close(30)

      do i=1,m_nu
         call gaunt_ff_calc(nu,m_nu,T_ISM,Zq_ISM,gaunt_str(i))
      end do

      do i=1,m_nu
         alpha_str(i) = 3.7d8 * (Zq_ISM**2*n_ISM**2/sqrt(T_ISM)) * (1.d0-exp(-h*nu(i)/(k*T_ISM)) )/nu(i)**3 * gaunt_str(i)
         tau_str(i) = alpha_str(i)*R_str
      end do

      open (25,file='CD_tauff.dat')
      do l=1,m_l
         do i=1,m_nu
            if (shock == 'BOTH') then
               tau_tot(l,i) = tauff(l,i) + tau_IL(l,i) + tau_str(i)
            else
               tau_tot(l,i) = tauff(l,i)
            end if
            write(25,*) l,nu(i),tau_tot(l,i)
         end do
      end do
      close(25)

      !/////////////////////////////////////////////////////////////////////

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
         open (33,file='emisyn_FS.dat')
         do l=1,m_theta
            do i=1,m_nu
               read(32,*) nl,Ef(i),Lsy_RS(l,i)  ! Lsy is already multiplied by Ef(i)
               Lsy(l,i) = Lsy_RS(l,i) + Lsy_FS(l,i)
               Lsy(l+m_theta,i) = Lsy(l,i)   ! symmetry 
         end do
         end do
         close(32)
         close(33)
      end if

      !-----------------------------------------------------------------------
      ! Calculate the total and absorbed emission (sum over linear emitters)

      do l=1,m_l
         do i=1,m_nu
            Lsy_tot(i) = Lsy_tot(i) + Lsy(l,i)
            if (i < inu_max) then !  The abs was calculated at these freqs  
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy(l,i)*exp(-tau_tot(l,i))  
            else
               Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy(l,i)
            endif
         end do
      end do

      !-----------------------------------------------------------------------
      ! Emission maps at the frequencies selected in controls (maps(:)),
      ! projected in the plane of the sky (u_map,v_map)

      n_map = size(maps)
      allocate( i_nu(n_map), u_map(m_l), v_map(m_l), emis(m_l,n_map), emis_abs(m_l,n_map) )

      i_nu     = nu_index(maps%nu, nu(1), log_nu_int, m_nu)
      u_map    = (y + x/abs(tginc))*abs(sinc)*convert
      v_map    = z*convert
      emis     = Lsy(:,i_nu) * spread(flux/nu(i_nu), 1, m_l)   ! Lsy is already multiplied by Ef
      emis_abs = emis * exp(-tau_tot(:,i_nu))

      if (shock == 'BOTH') then
         tag = ''
      else
         tag = shock
      end if

      do im = 1, n_map
         open(newunit=iu, file='CD_sync_radiomap_'//tag//trim(maps(im)%label)//'.dat')
         write(iu,'(4es16.8)') (u_map(l), v_map(l), emis(l,im), emis_abs(l,im), l=1,m_l)
         close(iu)
      end do

      deallocate( i_nu, u_map, v_map, emis, emis_abs )

      !------------------------------------------------------------
      ! Avg flux btw nu_av_min and nu_av_max
      if (shock /= 'BOTH') then
         open(76,file='CD_sync_radiomap_avg'//shock//'.dat')
      else if (shock == 'BOTH') then
         open(76,file='CD_sync_radiomap_avg.dat')
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
               emi_tot_abs(l) = emi_tot_abs(l) + Lsy(l,i)*exp(-tau_tot(l,i))*flux/nu(i) * dnu_avg
            end if
         end do
         emi_avg(l) = emi_tot(l)/(nu_av_max-nu_av_min)
         emi_avg_abs(l) = emi_tot_abs(l)/(nu_av_max-nu_av_min)
         write(76,*) u_proj,v_proj,emi_avg(l),emi_avg_abs(l)
      end do
      close(76)

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

      deallocate( emi_tot,emi_avg,emi_tot_abs,emi_avg_abs )
      deallocate( H_theta, H_IL, rho, n )
      deallocate( alpha_IL, alpha_str, alphaffw )
      deallocate( tauff,tau_IL,tau_str,tau_tot )
      deallocate( gaunt, gaunt_str )
      deallocate( x,y,z,nu,Lsy,Lsy_tot,Lsy_tot_abs,Lsy_RS,Lsy_FS )

    END subroutine CD_abssyn_run
    
   END module CD_abssyn