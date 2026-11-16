module opac_ff

   private
   public :: opac_ff_run

contains

   !============================================================================
   ! opac_ff_run
   !============================================================================
   !This module calculates the f-f absorption and emission of the unshocked SW 
   !The outputs of the program are:
   !1) windalphaff.dat: table with nu(Hz) and opacity alphaff(r_ff,nu) [1/cm] (r is the distance to the star)
   !2) windemiff.dat: table with nu(Hz) and the emissivity epsff(l,nu)*dl^3 [mJy]
   !3) r_ff.dat, a list of r_ff(n)
   !============================================================================

   subroutine opac_ff_run
      use global   
      !$ use omp_lib
      implicit none
      real(dp) :: T_w,S(mE_ff),ne,ni,gffinter
      real(dp) :: alphaff(mr_ff,mE_ff),gamma2,dg,dmu
      real(dp) :: r_ff(mr_ff),dr,A_ff,coef_alpha_ff
      integer :: j_mu,j_g
      integer, parameter :: nu_mu=146,nu_mga=81 ! # data in table gauntff.dat
      real(dp) :: gff(nu_mu,nu_mga),mu_ff(mE_ff)     
      real(dp) :: Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)
      interface
         subroutine interpol(dx,dy,fx1y1,fx2y1,fx1y2,finterpol)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp), intent(in) :: dx,dy,fx1y1,fx2y1,fx1y2
            real(dp), intent(out) :: finterpol
         end subroutine interpol
      end interface

      allocate( nu(mE_ff) )

      write(*,*) '\n Calculate the f-f abs. and emi. coefficients in the stellar wind'

      ! set omp threads number
      !$ call omp_set_num_threads(nth)
      ! !$ write(*,'(a, i3, a)') '# Using OpenMP with ', nth ,' threads'
      ! !$ write(*,*)

      !--------------------------------------------------------------------------
      ! Assume spherically symmetric winds      
      r_ff(1) = 5.*R_st        ! r_min (note that r_apex >> R_st)
      r_ff(mr_ff) = 8.d0*R0    ! r_max
      dr = (r_ff(mr_ff) - r_ff(1))/float(mr_ff)
      open (12,file='r_ff.dat')
      write(12,*) r_ff(1)
      do l = 2,mr_ff
         r_ff(l) = r_ff(l-1) + dr
         write(12,*) r_ff(l)
      end do
      close(12)

      !--------------------------------------------------------------------------
      ! Extract the Gaunt factor from the table given by van Hoof (2014):
      ! <gff(u,gamma2)>, with -8 < log10u < 8, -4 < log10gamma < 4, with step 1  
      open(13,file='gauntff.dat')
      do i=1,16
         read(13,*)
      end do
      do i=1,nu_mu
         read(13,*) (gff(i,j), j=1,nu_mga)
      end do
      close(13)

      !==========================================================================
      ! Calculate tauff(r,E) 

      T_w = 0.4d0*T_st      ! assume constant T in the wind (valid far from the star) 
      open(14,file='nu.dat')
      do i=1,mE_ff
         read(14,*) nu(i)
         mu_ff(i) = h*nu(i)/k/T_w   ! assuming T = cte, no need for mu_ff(r)
      end do
      close(14)

      A_ff = Mdot_eff/(4.d0*pi*v_w*mu_w*mp)
      gamma2 = Zq_w**2*Ry/k/T_w
      j_g = int( (log10(gamma2)+6.d0)/0.2d0 ) + 1
      dg = ( (log10(gamma2)+6.d0)/0.2d0+1.d0 - float(j_g) )*0.2d0

      !$omp parallel do &
      !$omp private(l,ni,ne,i,j_mu,dmu,gffinter) 
      do l=1,mr_ff
         ni = A_ff/r_ff(l)**2/(gamma_w+1.d0)  ! mui = mu*(g_e+1)
         ne = gamma_w*ni      
         coef_alpha_ff = 3.7d8*Zq_w**2*ne*ni/sqrt(T_w) ! For optimization
         do i=1,mE_ff
            ! For gff: i = (value - table initial value) / (table step) + 1
            ! After obtaining the evaluation point, interpolate linearly (2-variables)
            j_mu = int( (log10(mu_ff(i))+16.d0)/0.2d0 ) + 1
            dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2d0
            call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)
            !gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_sw/nu(i))) ! approx
            alphaff(l,i) = coef_alpha_ff*nu(i)**(-3.d0)*(1.d0-exp(-h*nu(i)/k/T_w))*gffinter !1/cm
         end do
      end do
      !$omp end parallel do

      open(18,file='windalphaff.dat')
      do l=1,mr_ff
         write(18,*) (alphaff(l,i), i=1,mE_ff)
      end do
      close(18)

      !==========================================================================
      ! Calculate the total wind radio flux using Eq. 8 from Wright+1975:

      S = 0.d0
      open(19,file='windemiff.dat')
      i=1
      do while (i < mE_ff .and. nu(i) < 1.d-1*k*T_w/h )! formula valid for h*nu << k*T 
         j_mu = int((log10(mu_ff(i))+16.d0)/0.2d0) + 1
         dmu = ((log10(mu_ff(i))+16.d0)/0.2d0+1.d0 - float(j_mu))*0.2
         call interpol(dg,dmu,gff(j_g,j_mu),gff(j_g+1,j_mu),gff(j_g,j_mu+1),gffinter)    
         !gffinter = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w1/nu(i))) ! approx          
         S(i) = (23.2d0/distance**2) * ( ((Mdot_eff/M_sun_yr)/(v_w/1.d5)/mu_w)**2&
            *nu(i)*gamma_w*gffinter*Zq_w**2 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
         if (nu(i) > 1.d-3*k*T_w/h) S(i) = S(i)*exp(-h*nu(i)/k/T_w)  ! artificial cutoff... 
         write(19,*) nu(i),S(i)*1.d3 ! Stored in mJy
         i=i+1
      end do
      close(19)

      deallocate(nu)

   END subroutine opac_ff_run

   !---------------------------------------------------------------------------

END module opac_ff

!=============================================================================


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

