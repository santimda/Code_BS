module CD_radbr

   private
   public :: CD_radbr_run

contains

   !============================================================================
   ! CD_radbr_run
   !============================================================================
   ! This module calculates the intrinsic bremsstrahlung emission. 
   ! If screened = .True. it is screend, if =.False. it is unscreend 
   ! The program reads as an input the electron energy distribution. 
   !============================================================================

   subroutine CD_radbr_run(screened) 
      use global
      Implicit none
      logical, intent(in) :: screened
      Real(dp) :: dummy,n_a(m_theta),dsigma_Br,cte_sigma,u_br
      Real(dp) :: Ef_int,dEf,sum,sum_br
      Real(dp) :: gamma_e(mE_e), Eps(mE_f), Ef(mE_f)
      Integer :: nl

      interface
         function Phi_Br(u_br,gamma_e,screened)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: Phi_Br
            real(dp), intent(in) :: u_br,gamma_e
            logical, intent(in) :: screened
         end function Phi_Br
      end interface

      allocate( Lbr(mE_f),Lbr_tot(mE_f) )
      allocate( B(m_theta),rho(m_theta),P(m_theta),v_t(m_theta),T(m_theta) )
      allocate( N_e(mE_e),E_e(mE_e),N_etot(m_theta,mE_e) )

      80    format(A7,X,Es10.2)

      write(*,*) '\n Calculate the (unabsorbed) rel. brem. emission in the SW'

      !------------------------------------------------------------------------
      !  Read information (thermodynamical quantities and relativist electrons)

      open(11,file='termo_sw.dat')
      do l=1,m_theta
         read(11,*) nl,rho(l),P(l),v_t(l),T(l),B(l) 
         n_a(l) = rho(l)/mp/mu_sw
      end do
      close(11)

      ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
      open (14,file='diste_sw.dat')
      do l=1,m_theta
         do i=1,mE_e
            read(14,*) dummy, E_e(i),N_etot(l,i)
            E_e(i) = E_e(i)*eV    !   energy from eV to erg
         end do
         read(14,*)
      end do
      close(14) 
      Eint_e = E_e(2)/E_e(1)

      !------------------------------------------------------------
      ! Build Ef vector for the SED with Ef_min = 1 keV
      call vector_log(1.d3*eV,E_e(mE_e),Ef_int,Ef)

      cte_sigma = 4.d0*alfa*re**2*Zq_sw**2*z_br  ! z_br > 1 if not pure H


      !===========================================================
      ! Calculate the relativistic Bremsstrahlung spectra


      Lbr_tot = 0.d0

      Eps = Ef/mec2
      gamma_e = E_e/mec2

      ! At each position we calculate the brem. emission; then we sum all.
      sum_br=0.d0
      DO l=1,m_theta
         do i=1,mE_e
            N_e(i)=N_etot(l,i)
         end do

         do i = 1,mE_f              ! Eps(i) 
            dEf = Ef(i)*(Ef_int-1.d0)
            sum = 0.d0
            do j = 1,mE_e              ! N_e(j),gamma_e(j)
               u_br = Eps(i)/gamma_e(j)     ! energy ratio
               if (u_br > 1.d0) then   ! Ef <= Ee from energy conservation
                  dsigma_Br = 0.d0  ![cm**2/m_e c2]
               else
                  dsigma_Br= cte_sigma*Phi_Br(u_br,gamma_e(j),screened) !/Ef(i), cancels later
                  if (dsigma_Br > 0.d0) then !Phi_Br can yield negative values in bad E ranges
                     P_ef = n_a(l)*c*dsigma_Br ! *Ef(i), cancels in dsigma_Br; [P_eps] = 1/seg
                     dE_e = E_e(j)*(Eint_e-1.d0)
                     dP_Ef = N_e(j)*P_ef*dE_e
                     sum = dP_Ef + sum ! P = int dP_Ef
                  endif
               endif
            end do ! j
            Lbr(i) = sum            ! [erg/seg]
            Lbr_tot(i) = Lbr_tot(i) + Lbr(i)

            sum_br = sum_br + Lbr(i)*dEf

         end do ! i

      end do ! l

      write(*,80) 'L_Brem=',sum_br*mphi*2.

      !---------------------------------------
      !
      !   Save data

      open(16,file='radbr_sw.dat')
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*Lbr_tot(i)
      end do
      close(16)

      deallocate( rho,P,v_t,T,B )
      deallocate( E_e,N_etot,N_e )
      deallocate( Lbr,Lbr_tot )

   END subroutine CD_radbr_run


   !========================================


END module CD_radbr


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
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u))*log(191.d0/z_br**(1./3.d0)) + (1.-u)/9.d0
      else                    ! un-screened nuclei 
      Phi_Br = ( 1. + (1.-u)**2 - (2./3.d0)*(1.-u))*( log(2.*g_e*(1./u-1.)) -0.5d0 )
   endif

   return
End function Phi_Br

!========================================
