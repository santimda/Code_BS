module radpp

   private
   public :: radpp_run

contains

   !============================================================================
   ! radpp_run
   !============================================================================
   ! This module calculates the p-p emission using Kelner+2006 formulae
   !============================================================================

   subroutine radpp_run(shock)   
      use global
      Implicit none
      character(LEN=*), intent(in) :: shock
      real(dp) :: Q_Ef,Q_phi,dx_pp,x_pp,qpi
      real(dp) :: Ef_int,dEf,sum_pp,Epp,L1,n_a 
      real(dp) :: E_th=1.22d9*eV, TeV=1.d12*eV
      real(dp) :: Epi_min,Epi_max,Epi_int,dE_pi,integral
      real(dp), parameter :: K_pi=0.17d0, n_pp=0.86d0 ! n_pp=1.1 for alfa=2, n_pp=0.86 for alpha=2.5
      integer :: mE_pi,i_p
      interface
         function F(a,b)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: F
            real(dp), intent(in) :: a,b
         end function F
      end interface

      interface
         function sigma_pp(L,Ep)
            use iso_fortran_env, only: dp => real64
            implicit none
            real(dp) :: sigma_pp
            real(dp), intent(in) :: L,Ep
         end function sigma_pp
      end interface

      allocate( E_p(mE_p),N_p(mE_p), Ef(mE_f),Lpp(mE_f) )

      80  format(A11,X,Es10.2)

      write(*,*) '\n Calculate the p-p emission in the SW'

      !------------------------------------------------------------------------
      !  Read information

      if (shock=='RS') then        
         n_a = n_sh_RS
         else if (shock=='FS') then     
         n_a = n_sh_FS
      end if

      ! N_p(i) = number of p+ per unit energy between [E_p(i), E_p(i)+dE]; [N] = 1/erg
      open (14,file='dist_p_'//shock//'.dat')
      do i=1,mE_p
         read(14,*) E_p(i),N_p(i)
         E_p(i) = E_p(i)*eV    !   eV to erg
      end do
      close(14) 
      Eint_p = E_p(2)/E_p(1)

      !------------------------------------------------------------
      ! Build logarithmic spaced photon energy vector, Ef(i), with Ef_min = 10 MeV
      call vector_log(10.d6*eV,E_p(mE_p),Ef_int,Ef)

      Lpp = 0.d0
      sum_pp = 0.d0  
      mE_pi = mE_f    ! Initial value
      do i = 1, mE_f           ! Ef(i) 
         dEf = Ef(i)*(Ef_int-1.d0)

         ! Different formalism according to energy range:

         ! Eq.78 Kelner+2006 for Ep < 100 GeV -> Ef < 0.01 TeV
         if (Ef(i) < 0.01d0*TeV) then   
            allocate(Epi(mE_pi))
            Epi_min = ( Ef(i)+(m_pic2)**2/4.d0/Ef(i) )*1.005 ! Factor to avoid numerical errors very close to Epi_min
            Epi_max = E_p(mE_p)
            call vector_log(Epi_min,Epi_max,Epi_int,Epi)
            integral = 0.d0
            do j=1,mE_pi         ! E_pi
               dE_pi = Epi(j)*(Epi_int-1.d0)
               Epp = Epi(j)/K_pi + mpc2       ! Eq.75 Kelner+2006
               i_p = int(log(Epp/Emin_p)/log(Eint_p))+1
               if (i_p <= mE_p .and. Epp <= E_p(mE_p) .and. Epp > E_th) then 
                  L1 = log(Epp/TeV)    
                  qpi = n_pp*N_p(i_p)*sigma_pp(L1,Epp)/K_pi ! Eq.77 Kelner+2006 [1/erg/s]
                  integral = integral + qpi*dE_pi/sqrt(Epi(j)**2-(m_pic2)**2) 
               endif
            end do               ! E_pi
            Q_Ef = 2.d0*c*n_a*integral*z_pp   !Eq.78 Kelner+2006
            Lpp(i) = Ef(i)*Q_Ef  ! SED = Ef * L_Ef =  Ef^2 * Q_Ef 
            deallocate(Epi)
            mE_pi = MAX(mE_pi-1,32) ! Diminish the number of points 

         ! Eq.72 Kelner+2006 for Ep > 100 GeV -> Ef > 0.01 TeV
         else 
            integral = 0.d0
            do j = 1, mE_p
               if (Ef(i) < E_p(j)) then
                  dE_p = E_p(j)*(Eint_p-1.d0)
                  x_pp = Ef(i)/E_p(j)
                  dx_pp = Ef(i)*dE_p/(E_p(j)**2)
                  L1 = log(E_p(j)/TeV) 
                  integral = integral + sigma_pp(L1,E_p(j))*N_p(j)*F(x_pp,L1)*dx_pp/x_pp 
               endif
            end do
            Q_phi = c*n_a*integral*z_pp
            Lpp(i) = Ef(i)*Q_phi

         endif

         sum_PP = sum_PP + Lpp(i)*dEf 
      end do    ! i

      write(*,80) 'L_pp =',sum_pp

      !---------------------------------------------------------
      ! Write the SED to a file

      open (19,file='radpp_'//shock//'.dat') 
      do i=1, mE_f 
         write(19,*) Ef(i)/eV,Ef(i)*Lpp(i)  
      end do
      close(19)

      deallocate( E_p,N_p, Lpp, Ef )


   END subroutine radpp_run


END module radpp


!========================================

Function F(x,L)
   integer, parameter :: dp = kind(1.d0)
   real(dp), intent(in) :: x,L
   real(dp) :: B,k,beta,F,L2
   L2 = L**2
   B = 1.30d0 + 0.14d0*L + 0.011d0*L2
   beta = 1.d0/(1.79d0 + 0.11d0*L + 0.008d0*L2)
   k = 1.d0/(0.801d0 + 0.049d0*L + 0.014d0*L2)
   F = B*log(x)/x*((1.d0-x**beta)/(1.d0+k*x**beta*(1.d0-x**beta)))**4 &
   *(1.d0/log(x)-(4.d0*beta*x**beta)/(1.d0-x**beta) &
      -(4.d0*k*beta*x**beta*(1.d0-2.d0*x**beta)/(1.d0+k*x**beta*(1.d0-x**beta))))
   return
end Function F

!========================================

Function sigma_pp(L1,Ep)
   use global
   real(dp), intent(in) :: L1,Ep
   real(dp) :: E_th=1.22d9*eV
   real(dp) :: sigma_pp
   sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
   sigma_pp = sigma_pp*(1.d0-(E_th/Ep)**4)**2 ! Add factor for better fit at low energies 
   return
end Function sigma_pp