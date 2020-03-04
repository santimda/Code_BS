module radbr

   private
   public :: radbr_run

contains

   !============================================================================
   ! radbr_run
   !============================================================================
   ! This module calculates the intrinsic bremsstrahlung emission. 
   ! If screened = .True. it is screend, if =.False. it is unscreend 
   ! The program reads as an input the electron energy distribution. 
   !============================================================================

   subroutine radbr_run(screened) 
      use global
      implicit none
      logical, intent(in) :: screened
      real(dp) :: Ef_int,dEf,sum,sum_br,dsigma_Br,u_br
      real(dp) :: gamma(mE_e),cte_sigma,dq_Ef,q_Ef
      interface
         function Phi_Br(u_br,gamma_e,screened)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: Phi_Br
            real(dp), intent(in) :: u_br,gamma_e
            logical, intent(in) :: screened
         end function Phi_Br
      end interface

      allocate( E_e(mE_e),N_e(mE_e), Ef(mE_f),L_br(mE_f) )

      81  format(A11,X,Es10.2)

      write(*,*) '\n Calculate the (unabsorbed) rel. brem. emission in the SW'

      cte_sigma = 4.d0*alfa*re**2*z_br

      !------------------------------------------------------------------------
      !  Read information

      ! N_e(i) = number of e- per unit energy between [E_e(i), E_e(i)+dE]; [N] = 1/erg
      open (14,file='dist_e.dat')
      do i=1,mE_e
         read(14,*) E_e(i),N_e(i)
         E_e(i) = E_e(i)*eV    !   Pasamos la energía a ergios
      end do
      close(14) 
      Eint_e = E_e(2)/E_e(1)


      !------------------------------------------------------------
      ! Build Ef vector for the SED with Ef_min = 10 keV
      call vector_log(1.d4*eV,E_e(mE_e),Ef_int,Ef)
      gamma = E_e/mec2  ! electron Lorentz factor

      !===========================================================
      ! Calculate the relativistic Bremsstrahlung spectra

      sum_br = 0.d0
      do i = 1,mE_f                 ! loop in high-energy photons energy (Ef)
         dEf = Ef(i)*(Ef_int-1.d0)
         sum = 0.d0
         do j = 1,mE_e              ! loop in electron energy distribution (Ne(j),gamma(j))
    
            u_br = Ef(i)/E_e(j)     ! energy ratio
            if (u_br > 1.d0) then   ! Ef <= Ee from energy conservation
               dsigma_Br = 0.d0     ![cm**2/m_e c2]
            else
               dsigma_Br = cte_sigma/Ef(i)*Phi_Br(u_br,gamma(j),screened) 
               if (dsigma_Br > 0.d0) then ! Phi_Br can yield negative values in bad E ranges
                  q_Ef = n_a*c*dsigma_Br
                  dE_e = E_e(j)*(Eint_e-1.d0)
                  dq_Ef = N_e(j)*q_ef*dE_e
                  sum = dq_Ef + sum ! q = int dq_Eps
               endif 
            endif
         end do 
         L_br(i) = Ef(i)*sum         ! transform photon emissivity to specific luminosity [erg/s/erg]
         sum_br = sum_br + L_br(i)*dEf
      end do 

      write(*,81) 'L_Br(tot)=',sum_br

      open(16,file='radbr.dat')
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*L_br(i)
      end do
      close(16)

      deallocate( E_e,N_e, Ef,L_br )

   END subroutine Radbr_run


   !========================================


END module radbr


!========================================

Function Phi_Br(u,g_e,screened)
   ! Calculates Phi function in rel.Br. cross-section assuming a factor z_br given in global
   ! It can give negative values for un-screened nuclei in certain E ranges.
   use global  
   implicit none
   real(dp) :: Phi_Br
   real(dp), intent(in) :: u,g_e
   logical, intent(in) :: screened

   if (screened) then      ! screened nuclei
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u))*log(191.d0/z_br**(1./3.d0)) + 1./9.d0*(1.-u)
   else                    ! un-screened nuclei 
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u))*( log(2.*g_e*(1./u-1.)) -0.5d0 )
   endif

   return
End function Phi_Br

!========================================