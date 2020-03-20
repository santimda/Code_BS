module radic

   private
   public :: radic_run

contains

   !============================================================================
   ! radic_run
   !============================================================================
   ! This module calculates the intrinsic IC emission using Khangulyan+2014
   ! formulae for black-body-like photon fields. 
   ! It calculates anisotropic IC-star if aniso=.True., else isotropic. 
   ! The program reads as an input the electron energy distribution. 
   !============================================================================

   subroutine radic_run(T,Kappa,aniso,target_field) 
      use global
      implicit none
      real(dp), intent(in) :: T,Kappa
      logical, intent(in) :: aniso 
      character(LEN=*), intent(in) :: target_field
      real(dp) :: Eps(mE_f),dEps,Eps_int,Ef_int,gamma(mE_e),x0,z,t_K
      real(dp) :: Ndot,sum,sum_ic,cte_Ndot,T_ad
      real(dp) :: T_ad_ang
      allocate( N_e(mE_e),E_e(mE_e), Ef(mE_f),L_ic(mE_f) )

      80  format(A11,X,Es10.2)

      write(*,*) '\n Calculate the (unabsorbed) IC-',target_field,' emission in the SW'

      !------------------------------------------------------------------------
      !  Read information

      ! N_e(i) = number of e- per unit energy between [E_e(i), E_e(i)+dE]; [N] = 1/erg
      open (14,file='dist_e.dat')
      do i=1,mE_e
         read(14,*) E_e(i),N_e(i)
         E_e(i) = E_e(i)*eV    !   [erg]
      end do
      close(14) 
      Eint_e = E_e(2)/E_e(1)
      gamma = E_e/mec2

      !------------------------------------------------------------
      ! Build Ef vector for the SED
      call vector_log(k*T,E_e(mE_e),Ef_int,Ef)
      Eps = Ef/mec2
      Eps_int = Ef_int

      !===========================================================
      ! Calculate the IC spectra

      T_ad = k*T/mec2     ! adimensional temperature
      T_ad_ang = T_ad     ! dummy definition so compiler doesn't warn possible uninitialization
      cte_Ndot = 2.d0*re**2*mec2**3*Kappa*T_ad**2/(pi*h_bar**3*c**2)

      if (aniso) T_ad_ang = T_ad*(1.-cos(incli)) 
      
      sum_ic = 0.d0
      do i = 1,mE_f                 ! loop in high-energy photons energy
         dEps=Eps(i)*(Eps_int-1.d0)
         sum = 0.d0
         do j = 1,mE_e              ! loop in electron energy distribution (Ne(j),gamma(j))
            z = Eps(i)/gamma(j)
            if (z < 1.d0) then      ! Ef <= Ee from energy conservation
               if (aniso) then                !anisotropic IC
                  t_K = 4.*gamma(j)*T_ad_ang  
                  x0 = z/( (1.-z)*t_K )
                  Ndot = (cte_Ndot/gamma(j)**2) * ( (z**2/(2.*(1.-z)))*F_1(x0) + F_2(x0) ) 
               else
                  t_K = 4.*gamma(j)*T_ad      !isotropic IC
                  x0 = z/( (1.-z)*t_K )
                  Ndot = (cte_Ndot/gamma(j)**2) * ( (z**2/(2.*(1.-z)))*F_3(x0) + F_4(x0) ) 
               endif
               P_eps = Eps(i)*Ndot !N^dot = nf*c*dsigmaIC
               dE_e = E_e(j)*(Eint_e-1.d0)
               dP_eps = N_e(j)*P_eps*dE_e
            else
               dP_eps = 0.d0
            endif
            sum = dP_eps + sum
         end do
         L_ic(i) = sum 
         sum_ic = sum_ic + L_ic(i)*dEps*mec2
      end do

      write(*,80) 'L_IC =',sum_ic

      !---------------------------------------------------------
      ! Write the SED to a file

      open(18,file='radic_'//target_field//'.dat')
      do i = 1,mE_f
         write(18,*) Ef(i)/eV,Ef(i)*L_ic(i)
      end do
      close(18)

      deallocate( N_e,E_e, Ef, L_ic )

   END subroutine radic_run


   !==============================


END module radic



