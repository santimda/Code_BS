module opac_ff
  !This module calculates the f-f absorption and emission of the unshocked SW 
  use global
  
  implicit none
  private
  public :: opac_ff_run

  contains

  !============================================================================
  ! opac_ff_run
  !============================================================================
  
    subroutine opac_ff_run
      use global   
      !$ use omp_lib
      implicit none

      !/////////////////////////////////////////////////////////////////////////
      !     The output of the program are:
      !     1) windalphaff.dat, which gives a table with nu(Hz) and opacity 
      !     alphaff(r_ff,nu) [1/cm], where r is the distance to the star.
      !     2) windemiff.dat, which gives a table with nu(Hz) and the emissivity 
      !     epsff(l,nu)*dl^3 [mJy]
      !     3) r_ff.dat, a list of r_ff(n)
      !/////////////////////////////////////////////////////////////////////////

      real(dp) :: Ry,T_w,S(mE_ff),ne,ni,gffinter
      real(dp) :: alphaff(mr_ff,mE_ff),gamma2,dg,dmu
      real(dp) :: r_ff(mr_ff),dr,A_ff,cte_alpha_ff
      integer :: j_mu,j_g
      integer, parameter :: nu_mu=146,nu_mga=81 ! # data in table gauntff.dat
      real(dp) :: gff(nu_mu,nu_mga),mu_ff(mE_ff)     
      real(dp) :: Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)
      !Variables for radio maps:
      integer :: i_1G,i_5G
      real(dp) :: log_nu_int,nu_S,B_ff,B_ff_D2,K_ff,cte_Tau
      real(dp) :: r_map(mr_map),S_ff_cm2(mr_map),A_cell(mr_map),S_map(mr_map)
      real(dp) :: convert,Delta_r,x_map(mphi),y_map(mphi),phi_map(mphi)
      interface
         subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
           use iso_fortran_env, only: dp => real64
           implicit none
           real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
           real(dp), intent(out) :: finterpol
         end subroutine interpol
      end interface

      interface
         subroutine S_ff_map(B_ff_D2,cte_Tau,r_i,r_f,S_ff_cm2,r_map)
           use global
           real(dp), intent(in) :: B_ff_D2,cte_Tau,r_i,r_f
           real(dp), intent(out) :: r_map(mr_map),S_ff_cm2(mr_map)
         end subroutine S_ff_map
      end interface

      allocate ( nu(mE_ff) )
      S = 0.d0
      Ry = 2.17987d-11          ! Rydberg unit = 2.18e-18 J

      ! set omp threads number
      !$ call omp_set_num_threads(nth)
      ! !$ write(*,'(a, i3, a)') '# Using OpenMP with ', nth ,' threads'
      ! !$ write(*,*)
      
!///////////////////////////////////////////////////////////////////////////
!
!    Assume spherically symmetric winds
          
      r_ff(1) = 10.*R_st        ! r_apex >> R_st
      r_ff(mr_ff) = 5.d0*R0 
      dr = (r_ff(mr_ff) - r_ff(1))/float(mr_ff)
      open (12,file='r_ff.dat')
      write(12,*) r_ff(1)
      do l = 2,mr_ff
         r_ff(l) = r_ff(l-1) + dr
         write(12,*) r_ff(l)
      enddo
      close(12)


!///////////////////////////////////////////////////////////////////////////
!
!     We extract the Gaunt factor from the table given by van Hoof (2014):
!     <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, with step 1
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
!     Calculate tauff(r,E) 
!
!///////////////////////////////////////////////////////////////////////////

      T_w = 0.4d0*T_st      ! assume constant T in the wind (valid far from the star) 
      open(14,file='nu.dat')
      do i=1,mE_ff
         read(14,*) nu(i)
         mu_ff(i) = h*nu(i)/k/T_w   ! assuming T = cte, no need for mu_ff(r)
      enddo
      close(14)      
      log_nu_int = log(nu(2)/nu(1))
      
      ! Falceta 2012?. Not so sure about this...
      A_ff = Mdot_eff/(4.d0*pi*v_w*mu_w*mp)
      gamma2 = Zq_w**2*Ry/k/T_w
      j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
      dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

      !$omp parallel do &
      !$omp private(l,ni,ne,i,j_mu,dmu,gffinter) 
      do l=1,mr_ff
         ni = A_ff/r_ff(l)**2./(gamma_w+1.d0)  ! mui = mu*(g_e+1)
         ne = gamma_w*ni      
         cte_alpha_ff = 3.7d8*Zq_w**2.*ne*ni*T_w**(-0.5_dp) ! For optimization
         do i=1,mE_ff
       
            ! For gff: i = (value - table initial value) / (table step) + 1
            ! After obtaining the evaluation point, interpolate linearly (2-variables)
            j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
            dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2_dp
            call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
!            gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w1/nu(i))) ! Una aprox
!            write(*,*) gffinter, 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w/nu(i)))
            alphaff(l,i) = cte_alpha_ff*nu(i)**(-3.d0)*(1.d0-exp(-h*nu(i)/k/T_w))*gffinter !1/cm
         enddo
      enddo
      !$omp end parallel do

      open(18,file='windalphaff.dat')
      do l=1,mr_ff
         write(18,*) (alphaff(l,i), i=1,mE_ff)
      enddo
      close(18)

      !///////////////////////////////////////////////////////////////////////////
      !     
      ! The total radio flux is given by Eq. 8 from Wright 1975:

      open(19,file='windemiff.dat')

      i=1
      do while (i < mE_ff .and. nu(i) < 1.d-1*k*T_w/h )! formula valid for h*nu << k*T 
         j_mu = int((log10(mu_ff(i))+16.d0)/0.2d0) + 1
         dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2
         call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)    
!         gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w1/nu(i))) ! aprox          
         S(i) = (23.2d0/distance**2.) * ( ((Mdot_eff/M_sun_yr)/&
              (v_w/1.d5)/mu_w)**2.d0 *nu(i)*gamma_w*&
              gffinter*Zq_w**2.d0 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
         if (nu(i) > 1.d-3*k*T_w/h) then
            S(i) = S(i)*exp(-h*nu(i)/k/T_w)      ! Invented... 
         endif      

         write(19,*) nu(i),S(i)*1.d3 ! Stored in mJy
         i=i+1
      enddo
      close(19)


      !============================================================
      !
      !   Fluxes for radio-maps
      
      phi_map(1)=0.d0
      do j=2,mphi
         phi_map(j) = phi_map(j-1) + delta_phi
      enddo
      convert = 206265.d0/distance             ! image to angular size ["]
                  
      !------------  1.4 GHz

      i_1G = int(log(1.42d9/nu(1))/log_nu_int) + 1        ! 1.4 GHz
      nu_S = nu(i_1G)

      B_ff = 2.*nu_S**2.*k*T_w/c**2.    ! Planck function at the h nu << kT limit
      B_ff_D2 = B_ff/distance**2.*1.d26 ! in mJy

      j_mu = int((log10(mu_ff(i_1G))+16.d0)/0.2d0) + 1
      dmu = ((log10(mu_ff(i_1G))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2
      call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)  
      K_ff = 3.7d8*h*Zq_w**2.*gffinter/(k*T_w**1.5*nu_S**2.)
      cte_Tau = K_ff*gamma_w*A_ff**2.*pi/2.
 !     write(*,*) gffinter,K_ff,cte_Tau

      call S_ff_map(B_ff_D2,cte_Tau,R_st,R0,S_ff_cm2,r_map)
      open(160,file='S_ff_map_1G.dat')
      Delta_r = r_map(2)-r_map(1)
      do i = 1,mr_map
         A_cell(i) = r_map(i)*Delta_r*d_phi
         S_map(i) = S_ff_cm2(i)*A_cell(i)
         do j=1,mphi
            y_map(j) = r_map(i) * sin(phi_map(j)) * convert
            x_map(j) = r_map(i) * cos(phi_map(j)) * convert
            write(160,*) x_map(j),y_map(j),S_map(i)
         enddo
      enddo
      close(160)

      !------------ 4.8 GHz

      i_5G = int(log(4.86d9/nu(1))/log_nu_int) + 1        ! 4.8 GHz
      nu_S = nu(i_5G)

      B_ff = 2.*nu_S**2.*k*T_w/c**2.    ! Planck function at the h nu << kT limit
      B_ff_D2 = B_ff/distance**2.*1.d26 ! in mJy

      j_mu = int((log10(mu_ff(i_5G))+16.d0)/0.2d0) + 1
      dmu = ((log10(mu_ff(i_5G))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2
      call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)  
      K_ff = 3.7d8*h*Zq_w**2.*gffinter/(k*T_w**1.5*nu_S**2.)
      cte_Tau = K_ff*gamma_w*A_ff**2.*pi/2.
!      write(*,*) gffinter,K_ff,cte_Tau
      
      call S_ff_map(B_ff_D2,cte_Tau,R_st,R0,S_ff_cm2,r_map)
      open(161,file='S_ff_map_5G.dat')
      do i = 1,mr_map
         A_cell(i) = r_map(i)*Delta_r*d_phi
         S_map(i) = S_ff_cm2(i)*A_cell(i)
         do j=1,mphi
            y_map(j) = r_map(i) * sin(phi_map(j)) * convert
            x_map(j) = r_map(i) * cos(phi_map(j)) * convert
            write(161,*) x_map(j),y_map(j),S_map(i)
         enddo
      enddo
      close(161)



  !////////////////////////////////////////////////////////////////////////////////

      deallocate ( nu )

    END subroutine opac_ff_run

  END module opac_ff

  !////////////////////////////////////////////////////////////////////////////////
  !

  Subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
    implicit none 
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
    real(dp), intent(out) :: finterpol
    real(dp) :: fx,fy
    fx = (fx2y1-fx1y1)/0.2    ! x2-x1 = 0.2
    fy = (fx1y2-fx1y1)/0.2    ! y2-y1 = 0.2
    finterpol = fx1y1 + fx*dx + fy*dy
    return
  END Subroutine interpol
  

  !////////////////////////////////////////////////////////////////////////////////
  !

  Subroutine S_ff_map(B_ff_D2,cte_Tau,r_i,r_f,S_ff_cm2,r_map)
    use global
    real(dp), intent(in) :: B_ff_D2,cte_Tau,r_i,r_f
    real(dp), intent(out) :: S_ff_cm2(mr_map),r_map(mr_map)
    real(dp) :: delta_r,Tau_max,R_c,coc_r

    R_c = 2.*r_i ! NOT SURE!
    delta_r = (r_f - r_i)/float(mr_map-1)
    r_map(1) = r_i
    do i = 2,mr_map
       r_map(i) = r_map(i-1) + delta_r
    enddo
    
    do i = 1,mr_map
       ! Eq. 5 from Wright1975
       if ( r_map(i) >= R_c ) then 
          Tau_max = cte_Tau/r_map(i)**3.
       else if ( r_map(i) < R_c ) then
          coc_r = r_map(i)/R_c
          Tau_max = cte_Tau/r_map(i)**3.* ( 0.5 - acos(coc_r)/pi -&
               sqrt(1.-coc_r**2.)*coc_r/pi )
       endif
       S_ff_cm2(i) = B_ff_D2*(1.-exp(-Tau_max))
    enddo
    
    return
  END Subroutine S_ff_map
