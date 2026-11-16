module radsyn
   use global

   implicit none
   private
   public :: radsyn_run

contains

   !============================================================================
   ! radsyn_run
   !============================================================================
   ! This module calculates the intrinsic synchrotron emission corrected at low
   ! frequencis for the Razin-Tsytovich effect. 
   ! Stores the SED in radsyn.dat: nu, E*L_sy(with R-T), E*L_sy(without R-T)
   !============================================================================

   subroutine radsyn_run(shock)  
      use global
      !$ use omp_lib
      implicit none
      character(LEN=*), intent(in) :: shock
      real(dp) :: Eps(m_nu),dEps,Eps_int, gamma(mE_e)
      real(dp) :: factor,nuc,nuc_max,nuc_min,sum,sum_sy
      real(dp) :: n_a,ne,nu_RT,x_sy,coef_Psy
      real(dp) :: log_nu_int,B,B_par
      allocate( N_e(mE_e),E_e(mE_e), nu(m_nu),L_sy(m_nu),L_sy_RT(m_nu) )

      80  format(A14,X,Es10.2)

      write(*,*) '\n Calculate the (unabsorbed) sync. emission in the SW'

      ! set omp threads number
      !$ call omp_set_num_threads(nth) 

      !------------------------------------------------------------------------
      !  Read information and define variables

      if (shock=='RS') then        
         n_a = n_sh_RS
         B = B_RS
         else if (shock=='FS') then     
         n_a = n_sh_FS
         B = B_FS
      end if
      B_par = B*par 

      ! N_e(i) = number of e- per unit energy between [E_e(i), E_e(i)+dE]; [N] = 1/erg
      open (14,file='dist_e_'//shock//'.dat')
      do i=1,mE_e
         read(14,*) E_e(i),N_e(i)
         E_e(i) = E_e(i)*eV   !   energy from eV to erg
      end do
      close(14) 
      Eint_e = E_e(2)/E_e(1)
      gamma = E_e/mec2        ! electron Lorentz factor


      ne = n_a*mu_sw/mu_e_sw
      nu_RT = 20.d0*ne/B  ! characteristic frequency of the Razin-Tsytovich effect
      write(*,80) 'nu_RT[GHz] =',nu_RT/1.d9  ! R-T characteristic frequency

      !------------------------------------------------------------
      ! Build nu vector for the SED
      factor = h/mec2
      nuc_max = 4.22d6*B_par*gamma(mE_e)**2 *9.d0    ! Gives an extra
      nuc_min = 10.d6 ! 4.22d6*B_par*gamma(1)**2*10.d0 ! 10 MHz
      call vector_log(factor*nuc_min,factor*nuc_max,Eps_int,Eps)
      nu = Eps/factor ! [nu] = Hz
      log_nu_int = log(nu(2)/nu(1))

      ! Also save the nu vector for other (f-f) programs
      open (88,file='nu.dat')
      do i=1,m_nu
         write(88,*) nu(i)
      end do
      close(88)


      !===========================================================
      ! Calculate the synchrotron spectra
      ! L_sy(eps) = P(eps) = int_Emin^Emax N(E) P(E,eps) dE [erg/s/eps]

      coef_Psy = 1.85*sqrt(3.)*qe**3*B_par/h/mec2 !4.39d-22*B_par/h
      sum_sy=0.d0
      do i = 1,m_nu           ! loop in emitted photon energies (Eps(i), nu(i))
         dEps = Eps(i)*(Eps_int-1.d0)
         sum = 0.d0

         !$omp parallel do private(nuc,dP_Eps,P_Eps,dE_e) &
         !$omp reduction(+ : sum)
         do j = 1,mE_e         ! loop in electron energy distribution (Ne(j), nuc(j), gamma(j))
            nuc = 4.22d6*B_par*gamma(j)**2
            x_sy = nu(i)/nuc
            if (x_sy <= 10.d0) then  
               P_Eps = coef_Psy*x_sy**(1.d0/3.d0)*exp(-x_sy)  ![erg/s/erg]
               dE_e = E_e(j)*(Eint_e - 1.d0)
               dP_Eps = N_e(j)*P_eps*dE_e
               sum = dP_Eps + sum
            endif
         end do
         !$omp end parallel do 

         L_sy(i) = sum*mec2
         L_sy_RT(i) = L_sy(i)*exp(-nu_RT/nu(i)) ! L_sy corrected for R-T effect
         sum_sy = sum_sy + L_sy_RT(i)*dEps
      end do  ! Eps(i),nu(i)

      write(*,80) 'L_syn(unabs) =',sum_sy

      !---------------------------------------------------------
      ! Write the SED to a file

      open(17,file='radsyn_'//shock//'.dat')
      do i=1,m_nu
         write(17,*) Eps(i)*mec2/eV,Eps(i)*L_sy_RT(i),Eps(i)*L_sy(i)
      end do
      close(17)

      deallocate( N_e,E_e, nu,L_sy,L_sy_RT )

   END subroutine radsyn_run


   !==============================


END module radsyn



