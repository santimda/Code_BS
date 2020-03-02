module CD_radbr

  !  This module calculates the relativistic bremsstrahlung emission

  use global
  
  implicit none
  private
  public :: CD_radbr_run
  
contains
  
  !============================================================================
  ! CD_radbr_run
  !============================================================================
  
  subroutine CD_radbr_run  
    use global
    Implicit none
      Real(dp) :: dummy,n_a(m_theta),dsigma_Br,cte_sigma
      Real(dp) :: Ef_int,dEf,sum,sum_br
      Real(dp) :: gamma(mE_f), Eps(mE_f), Ef(mE_f)
      Integer :: nl
      ! m_theta - CD, mE_e - electrons, mE_f - photons
   
      interface
         function Phi_Br(Ef,Ee,nuc)
           use iso_fortran_env, only: dp => real64
           implicit none
           real(dp) :: Phi_Br
           real(dp), intent(in) :: Ef,Ee
           character, intent(in) :: nuc
         end function Phi_Br
      end interface
      
      allocate ( Lbr(mE_f),Lbr_tot(mE_f) )
      allocate ( B(m_theta),rho(m_theta),P(m_theta),v_t(m_theta),T(m_theta) )
      allocate ( N_e(mE_e),E_e(mE_e),N_etot(m_theta,mE_e) )
      
80    format(A7,X,Es10.2)
      
      
      write(*,*) 'Calculate the (unabsorbed) rel. brem. emission in the SW'

      !////////////////////////////////////////////////////////////////////////
      !
      ! Read thermodynamical quantites at each position
      
      open(11,file='termo_sw.dat')
      do l=1,m_theta
         read(11,*) nl,rho(l),P(l),v_t(l),T(l),B(l) 
         n_a(l) = rho(l)/mp/mu_sw
      enddo
      close(11)
      
      ! Number of electrons at r(l) with energy E_e(i): N_etot(l,i) 
      open (14,file='diste_sw.dat')
      do l=1,m_theta
         do i=1,mE_e
            read(14,*) dummy, E_e(i),N_etot(l,i)
            E_e(i) = E_e(i)*eV    !   energy from eV to erg
         enddo
         read(14,*)
      enddo
      close(14) 
      Eint_e = E_e(2)/E_e(1)
      gamma = E_e/mec2
            
      
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !                                                                       C
      !          Calculate Bremsstrahlung emission                            C
      !                                                                       C
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      Lbr_tot = 0.d0
      
      cte_sigma = 4.d0*alfa*re**2.*Zq_sw**2.*2.785d0 !qe**4.d0 ?????????
      
      ! Build logarithmic-spaced energy vector for emitted photons, Eps(i):
      call vector_log(1.d3*eV,E_e(mE_e),Ef_int,Ef)

      Eps = Ef/mec2
      gamma = E_e/mec2
      
      ! At each position we calculate the brem. emission; then we sum all.
      sum_br=0.d0
      DO l=1,m_theta
         do i=1,mE_e
            N_e(i)=N_etot(l,i)
         enddo
         
         do i = 1,mE_f              ! Eps(i) 
            dEf = Ef(i)*(Ef_int-1.d0)
            sum = 0.d0
            do j = 1,mE_e           ! N_e(j),gamma(j)
               if (gamma(j) < Eps(i)) then
                  dsigma_Br = 0.d0  ![cm**2/m_e c2]
               else
                  dsigma_Br= cte_sigma/Ef(i)*Phi_Br(Eps(i),gamma(j),'U') !U or S
                  ! [P_eps] = 1/seg
                  ! Factor: Sum ni*Zi**2 ~ 0.705*1**2 + 0.275*2**2 + 0.020*7**2 ~ 2.785  ????????
                  if (dsigma_Br > 0.d0) then ! bad energy ranges may yield negative values
                    P_ef = Ef(i)*n_a(l)*c*dsigma_Br
                    dE_e = E_e(j)*(Eint_e-1.d0)
                    dP_Ef = N_e(j)*P_ef*dE_e
                    sum = dP_Ef + sum ! P = int dP_Ef
                  endif
               endif
            enddo ! j
            Lbr(i) = sum            ! [erg/seg]
            Lbr_tot(i) = Lbr_tot(i) + Lbr(i)

            sum_br = sum_br + Lbr(i)*dEf
            
         enddo ! i
               
      ENDDO ! l
      
      write(*,80) 'L_Brem=',sum_br*mphi*2.
      
    !---------------------------------------
    !
    !   Save data

      open(16,file='radbr_sw.dat')
      do i=1,mE_f
         write(16,*) Ef(i)/eV,Ef(i)*Lbr_tot(i)
      enddo
      close(16)
      
      deallocate ( rho,P,v_t,T,B )
      deallocate ( E_e,N_etot,N_e )
      deallocate ( Lbr,Lbr_tot )

    END subroutine CD_radbr_run
    

!========================================


  END module CD_radbr


!========================================



Function Phi_Br(Ef0,Ee,nucl)
  
  use global
  implicit none
  real(dp) :: Phi_Br
  real(dp), intent(in) :: Ef0,Ee  ! Normalized / mec2
  character, intent(in) :: nucl
  real(dp) :: Z_Br

  Z_Br = Zq_sw  ! ????

  if (nucl == 'S') then ! Screened
     
     Phi_Br = ( 1. + (1.-Ef0/Ee)**2.-(2./3.)*(1.-Ef0/Ee))*&
          log(191./Z_Br**(1./3.)) + (1./9.)*(1.-Ef0/Ee)
     
  else if (nucl == 'U') then ! Un-screened:
     
     Phi_Br = ( 1. + (1.-Ef0/Ee)**2.-(2./3.)*(1.-Ef0/Ee))*&
          (log(2.*Ee*(Ee-Ef0)/Ef0)-0.5)

  else
     
     Phi_Br = 0.d0
     write(*,*) 'ERROR IN Phi_BR ARGUMENT'
     
  endif
  
  return
End function Phi_Br
