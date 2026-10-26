module CD_distp
   use global
   use controls, only: verbose, calc_N_bin
   use numerical_utils, only: vector_log
   use diffusion, only: Diff
   implicit none
   private

   public :: CD_distp_run, CD_distp_approx_run

   real(dp), parameter :: k_pp = 0.5d0     ! inelasticity of the p-p interaction

   contains

   !============================================================================
   ! CD_distp_run (START)
   !============================================================================
   !  Proton energy distribution at each position of the shock (same scheme as
   !  CD_diste). Inputs: Linj_<shock>.dat, termo_<shock>.dat, Hnorm_<shock>.dat,
   !  efi_ac_<shock>.dat, alpha_<shock>.dat and t_conv_<shock>.dat.
   !  Outputs: distp_<system>_<shock>.dat, with Np(l,E(i)) (m_theta*mE_p data);
   !  ENTp_<system>_<shock>.dat, with E_NT(l); L_inj_p_<shock>.dat and
   !  several diagnostic files. Emin_p and Emax_p are fixed in global.
   !============================================================================

   subroutine CD_distp_run(shock)
      !$ use omp_lib
      character(LEN=*), intent(in) :: shock
      integer  :: l, i, l2, l_min, l_cut, m_theta_inj_max, n_E, iu
      integer  :: index1, index2, n_count, steps
      logical  :: new_count
      real(dp) :: Eint, L1, b_pp, D_diff, t_pp_l
      real(dp) :: Ka, integral, Np, cool, dE_aprox, suma
      real(dp) :: Q_E(mE_p)        ! private copy per OpenMP thread (small, on each thread's stack)
      ! Shock quantities (one value per position l)
      real(dp), allocatable :: theta(:), rho(:), B(:), v_t(:), n_a(:)
      real(dp), allocatable :: L_inj_p(:), dR(:), H_theta(:), R_acel(:), dV_ring(:)
      real(dp), allocatable :: efi_ac_sh(:), alpha_sh(:), E_corr(:)
      real(dp), allocatable :: t_cell(:), t_adv(:), t_adi(:), t_adi_pp(:)
      real(dp), allocatable :: log10Emax(:), Emax_p0(:), rg_p_max(:)
      real(dp), allocatable :: E_NT(:), U_NT(:), U_B(:)
      ! Energy grid and (l,E) quantities
      real(dp), allocatable :: E_p(:), dE(:), N_ptot(:,:), E_dot(:,:)
      real(dp), allocatable :: t_pp(:,:), t_coolp(:,:), t_diff(:,:), t_esc(:,:), t_loss(:,:), t_ac(:,:)

      allocate( theta(m_theta), rho(m_theta), B(m_theta), v_t(m_theta), n_a(m_theta) )
      allocate( L_inj_p(m_theta), dR(m_theta), H_theta(m_theta), R_acel(m_theta), dV_ring(m_theta) )
      allocate( efi_ac_sh(m_theta), alpha_sh(m_theta), E_corr(m_theta) )
      allocate( t_cell(m_theta), t_adv(m_theta), t_adi(m_theta), t_adi_pp(m_theta) )
      allocate( log10Emax(m_theta), Emax_p0(m_theta), rg_p_max(m_theta) )
      allocate( E_NT(m_theta), U_NT(m_theta), U_B(m_theta) )
      allocate( E_p(mE_p), dE(mE_p), N_ptot(m_theta,mE_p), E_dot(m_theta,mE_p) )
      allocate( t_pp(m_theta,mE_p), t_coolp(m_theta,mE_p), t_diff(m_theta,mE_p) )
      allocate( t_esc(m_theta,mE_p), t_loss(m_theta,mE_p), t_ac(m_theta,mE_p) )

      call validate_shock(shock, 'CD_distp_run')
      if (verbose) print *, 'Calculate the non-thermal proton distribution in the '//shock

      ! set omp threads number
      !$ call omp_set_num_threads(nth)

      !---------------------------------------------------------------------------------
      !  Read the cell size, injected luminosity, thermodynamical quantities, shock
      !  width, acceleration parameters and t_adv at each cell

      call read_shock_p(shock, theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, R_acel, dV_ring, &
                        efi_ac_sh, alpha_sh, t_adv)

      open(newunit=iu, file=output_file('L_inj_p_'//shock//'.dat'))
      write(iu,'(es17.8e3)') sum(L_inj_p)
      close(iu)

      !--------------------------------------------------------------------
      ! Build an energy vector that is common for all linear-emitters:

      call vector_log(Emin_p,Emax_p,Eint,E_p)
      dE = E_p*(Eint - 1.d0)

      !========================================================================
      ! First calculate the cooling time at all energies and at each position
      ! (and also the maximum energy)

      t_cell = dR/v_t                ! Cell advection time

      ! Characteristic adiabatic losses time
      t_adi(1:m_theta-1) = 3.d0*dR(1:m_theta-1)/(v_t(1:m_theta-1)*log(rho(1:m_theta-1)/rho(2:m_theta)))
      t_adi(m_theta) = t_adi(m_theta-1)
      where (t_adi <= 0.d0 .or. t_adi > 1.d20) t_adi = 1.d20   ! To avoid numerical errors (also Inf)

      ! Energy-independent p-p time (very approx, just to calculate E_max)
      t_adi_pp = 1.d0/(1.d0/t_adi + n_a*z_pp/1.d15)

      !$omp parallel do default(shared) schedule(dynamic) private(l, i, L1, b_pp, D_diff, t_pp_l)
      do l = 1, m_theta
         ! Maximum energy (the cut in l is applied below, after the loop)
         call Biseccion_p(B(l),log10Emax(l),t_adi_pp(l),t_adv(l),efi_ac_sh(l),R_acel(l))

         do i = 1, mE_p
            if (E_p(i) <= E_th) then    ! E_th = 1.22 GeV
               t_pp_l = 1.d30           ! no p-p losses below the threshold
            else
               L1 = log(E_p(i)/TeV)     ! L = ln(Ep/1TeV)
               b_pp = c*k_pp*n_a(l)*sigma_pp(L1,E_p(i))*E_p(i)
               t_pp_l = E_p(i)/(b_pp*z_pp)
            endif
            t_pp(l,i) = t_pp_l

            D_diff = Diff(E_p(i),B(l),'B')
            t_diff(l,i) = R_acel(l)**2/(2.d0*D_diff)                   ! Diffusion losses
            t_esc(l,i)  = 1.d0/( 1.d0/t_adv(l) + 1.d0/t_diff(l,i) )    ! escape
         end do
      end do
      !$omp end parallel do

      t_coolp = 1.d0/(1.d0/t_pp + 1.d0/spread(t_adi,2,mE_p))          ! Cooling losses
      t_loss  = 1.d0/(1.d0/t_coolp + 1.d0/t_esc)                       ! general losses
      t_ac    = lambda_diff * spread(efi_ac_sh,2,mE_p)*spread(E_p,1,m_theta)/(spread(B,2,mE_p)*c*qe)
      E_dot   = spread(E_p,1,m_theta)/t_coolp

      ! Biseccion_p returns log10Emax < log10(1.5*Emin_p) when there is no acceleration:
      ! particle acceleration stops at the first such position, l_cut
      l_cut = findloc(log10Emax < log10(1.5d0*Emin_p), .true., dim=1)
      m_theta_inj_max = m_theta              ! Default value is to inject up to m_theta
      if (l_cut > 0) then
         m_theta_inj_max = l_cut
         log10Emax(l_cut) = log10(1.001d0*Emin_p)
      end if
      Emax_p0 = MIN( 10.d0**log10Emax, E_Hillas(B, R_acel, efi_ac_sh) )
      if (l_cut > 0) Emax_p0(l_cut+1:) = 1.0001d0*Emin_p   ! Dummy value, no particle acceleration
      if (any(Emax_p0 >= 0.5d0*Emax_p)) print *, 'Emax_p in global.f90 is too small! max(Emax_p0) [eV] =', maxval(Emax_p0)/eV

      rg_p_max = Emax_p0/qe/B     ! Max gyroradius

      open(newunit=iu, file=output_file('CD_tcoolp_'//system//'_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(i6,10es17.8e3)') (l, E_p(i)/eV, t_ac(l,i), t_pp(l,i), t_adv(l), t_diff(l,i), &
              t_coolp(l,i), t_cell(l), t_esc(l,i), t_loss(l,i), t_adi(l), i=1,mE_p)
      end do
      close(iu)

      open(newunit=iu, file=output_file('Emaxp'//shock//'.dat'))
      write(iu,'(i6,es17.8e3)') (l, Emax_p0(l)/eV, l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('Gyroradius_p_'//shock//'.dat'))
      write(iu,'(i6,3es17.8e3)') (l, theta(l), rg_p_max(l), H_theta(l)/rg_p_max(l), l=1,m_theta)
      close(iu)

      !========================================================================

      ! Consider different starting points for each linear emitter
      E_corr = E_break**(alpha_1-alpha_sh)

      N_ptot = 0.d0
      ! At the starting point of the line we calculate the distribution
      ! Np_0(E_0), the cooling rate Edot0, and we sum it to Np_tot(l,i).
      ! Each starting point l_min is independent: OpenMP over l_min, adding
      ! the contributions to N_ptot with a reduction.

      !$omp parallel do default(shared) schedule(dynamic) reduction(+:N_ptot) &
      !$omp private(l_min, n_E, Q_E, integral, Ka, i, Np, index1, index2, l, l2, &
      !$omp         n_count, steps, new_count, cool, dE_aprox, suma)
      Do l_min = 1, m_theta_inj_max ! After m_theta_inj_max, particles are no longer accelerated
         ! Q = Ka*fQp(E). The normalization constant Ka is obtained by the condition:
         ! integral Q(E')*E'dE' = Linj. [Linj]=erg/s, [Ka]=erg/s, [Q]=1/erg/seg
         ! The integral is solved with the trapezoid method, up to 10*Emax.

         n_E = count(E_p < 10.d0*Emax_p0(l_min))
         if (n_E < 2) cycle
         Q_E(1:n_E) = fQp(E_p(1:n_E), E_break, Emax_p0(l_min), E_corr(l_min), alpha_sh(l_min), alpha_1)
         integral = sum( (E_p(1:n_E-1)*Q_E(1:n_E-1) + E_p(2:n_E)*Q_E(2:n_E))/2.d0 * dE(1:n_E-1) )
         Ka = L_inj_p(l_min)/integral ! Injected at the starting point of the line

         !------------------------------------------------------------------------

         do i = 1, n_E
            ! N(E) = Q(E)*dt, with dt = min(t_cell, t_cool) and diffusion treated
            ! as an effective energy loss (as in CD_diste); [N(E)] = 1/erg
            Np = Ka*Q_E(i) * min(t_cell(l_min), 1.d0/(1.d0/t_coolp(l_min,i) + 1.d0/t_diff(l_min,i)))
            N_ptot(l_min,i) = N_ptot(l_min,i) + Np

            ! Evolved version of the injected distribution at each position "l":
            ! N_l = N_l-1*(Edot_l-1/Edot_l)*(t_cell(l)/t_cell(l-1)), where the energy
            ! evolves as E0(l) given by t_cell(l) = int_E0^El dE/Edot

            index1 = i            ! For the cell "l"
            index2 = i            ! For the previous cell, "l-1"
            l = l_min + 1
            n_count = 0           ! Counter for when to change the index if losses are small
            cool = 0.d0
            new_count = .true.

            do while( l <= m_theta .and. index1 > 1) ! if index1 = 1, the p+ are completely cooled
               ! The cooling parameter "cool" determines whether particles cool efficiently or not
               if (l == l_min+1) then ! The total cooling parameter from l_min to m_theta
                  cool = sum( t_cell(l:m_theta)/t_coolp(l:m_theta,index1) )
               else
                  cool = cool - t_cell(l-1)/t_coolp(l-1,index1)
               end if
               dE_aprox = (1.d0-t_cell(l)/t_coolp(l,index1))

               if ( 1.d0-cool > 0.99d0/Eint ) then ! negligible cooling
                  continue

               elseif ( 1.d0-cool <= 0.99d0/Eint .and. dE_aprox > 1.d0/Eint ) then
                  if (new_count) then
                     suma = 1.d0
                     steps = 0 ! number of steps until it changes E(index1) -> E(index1-1)
                     l2 = l
                     do while (suma > 1.d0/Eint .and. l2 <= m_theta)
                        steps = steps + 1
                        suma = suma*(1.d0-t_cell(l2)/t_coolp(l2,index1))
                        l2 = l2 + 1
                        if (suma > 1.d0/Eint .and. l2 == m_theta) steps = m_theta ! NOT COOLING
                     end do
                     new_count = .false.
                  end if
                  n_count = n_count + 1

                  if (n_count >= steps) then
                     n_count = 0
                     index1 = index1 - 1
                     new_count = .true.
                  end if

               elseif ( dE_aprox <= 1.d0/Eint ) then
                  suma = 0.d0
                  do while (suma < t_cell(l) .and. index1 > 1)
                     index1 = index1 - 1
                     suma = suma + dE(index1)/E_dot(l,index1)
                  end do
               end if

               ! The following takes into account energy losses and variable cell size and velocity
               Np = Np * ( E_dot(l,index2)/E_dot(l,index1) ) * (t_cell(l)/t_cell(l-1))

               N_ptot(l,index1) = N_ptot(l,index1) + Np  ! The energies match by construction

               index2 = index1  ! index2 keeps the value of the previous step

               l = l + 1
            end do              ! End loop in shock position (l)
         end do                 ! End loop in proton starting energies (i)

      end do                    ! End loop in starting position (l_min = # linear-emitters)
      !$omp end parallel do

      !----------------------------------------------------------------------
      ! Save the proton distribution at each position
      open(newunit=iu, file=output_file('distp_'//system//'_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(3es25.16e3)') (theta(l), E_p(i)/eV, N_ptot(l,i), i=1,mE_p)
         write(iu,*)
      end do
      close(iu)

      ! Non-thermal energy in each BS "ring" (4 quadrants), its density and ratio to U_B
      E_NT = matmul(N_ptot, E_p*dE) * 4.d0*mphi
      U_NT = E_NT/dV_ring
      U_B  = B**2/(8.d0*pi)
      open(newunit=iu, file=output_file('ENTp_'//system//'_'//shock//'.dat'))
      write(iu,'(4es17.8e3)') (theta(l), E_NT(l), U_NT(l), U_NT(l)/U_B(l), l=1,m_theta)
      close(iu)

      !------------------------------------------------------------
      ! Group in bins according to position angle
      if (calc_N_bin) call Np_bin(E_p,N_ptot,theta,shock)

      !------------------------------------------------------------

      deallocate( theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, R_acel, dV_ring )
      deallocate( efi_ac_sh, alpha_sh, E_corr, t_cell, t_adv, t_adi, t_adi_pp )
      deallocate( log10Emax, Emax_p0, rg_p_max, E_NT, U_NT, U_B )
      deallocate( E_p, dE, N_ptot, E_dot, t_pp, t_coolp, t_diff, t_esc, t_loss, t_ac )

      !call CD_distp_approx_run(shock) ! JUST FOR A CHECK

   END subroutine CD_distp_run

   !============================================================================
   ! CD_distp_run (END)
   !============================================================================

   !============================================================================
   ! CD_distp_approx_run (START)
   !============================================================================
   !  Same as CD_distp_run, but NEGLECTING ENERGY LOSSES (JUST FOR A CHECK).
   !  Emax is limited by advection (t_adv = t_ac) and by the Hillas criterion.
   !  Without losses, N(l,E) = t_cell(l) * sum_{l_min <= l} Q_lmin(E).
   !  Outputs: distp_approx_<shock>.dat, ENTp_approx_<shock>.dat and
   !  distp_approx_<shock>_total.dat.
   !============================================================================

   subroutine CD_distp_approx_run(shock)
      !$ use omp_lib
      character(LEN=*), intent(in) :: shock
      integer  :: l, i, l_min, n_E, iu
      real(dp) :: Eint, integral
      real(dp), allocatable :: theta(:), rho(:), B(:), v_t(:), n_a(:)
      real(dp), allocatable :: L_inj_p(:), dR(:), H_theta(:), R_acel(:), dV_ring(:)
      real(dp), allocatable :: efi_ac_sh(:), alpha_sh(:), E_corr(:), t_adv(:), t_cell(:)
      real(dp), allocatable :: Emax_p0(:), E_NT(:), U_NT(:), U_B(:)
      real(dp), allocatable :: E_p(:), dE(:), Q_inj(:,:), N_ptot(:,:)

      allocate( theta(m_theta), rho(m_theta), B(m_theta), v_t(m_theta), n_a(m_theta) )
      allocate( L_inj_p(m_theta), dR(m_theta), H_theta(m_theta), R_acel(m_theta), dV_ring(m_theta) )
      allocate( efi_ac_sh(m_theta), alpha_sh(m_theta), E_corr(m_theta), t_adv(m_theta), t_cell(m_theta) )
      allocate( Emax_p0(m_theta), E_NT(m_theta), U_NT(m_theta), U_B(m_theta) )
      allocate( E_p(mE_p), dE(mE_p), Q_inj(m_theta,mE_p), N_ptot(m_theta,mE_p) )

      if (verbose) print *, 'Calculate the non-thermal proton distribution in the '//shock//' (no energy losses)'

      ! set omp threads number
      !$ call omp_set_num_threads(nth)

      call read_shock_p(shock, theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, R_acel, dV_ring, &
                        efi_ac_sh, alpha_sh, t_adv)

      call vector_log(Emin_p,Emax_p,Eint,E_p)
      dE = E_p*(Eint - 1.d0)

      t_cell = dR/v_t                                  ! Cell advection time

      ! Maximum energy: advection (t_adv = t_ac) or Hillas
      Emax_p0 = MIN( t_adv*(B*c*qe)/(lambda_diff*efi_ac_sh), E_Hillas(B, R_acel, efi_ac_sh) )
      if (any(Emax_p0 >= 0.5d0*Emax_p)) print *, 'Emax_p in global.f90 is too small! max(Emax_p0) [eV] =', maxval(Emax_p0)/eV

      ! Injection Q_inj(l_min,E) = Ka*fQp(E) [1/erg/s], normalized with integral Q(E')*E'dE' = Linj
      E_corr = E_break**(alpha_1-alpha_sh)
      Q_inj  = 0.d0

      !$omp parallel do default(shared) schedule(dynamic) private(l_min, n_E, integral)
      do l_min = 1, m_theta
         n_E = count(E_p < 10.d0*Emax_p0(l_min))
         Q_inj(l_min,1:n_E) = fQp(E_p(1:n_E), E_break, Emax_p0(l_min), E_corr(l_min), alpha_sh(l_min), alpha_1)
         integral = sum( (E_p(1:n_E-1)*Q_inj(l_min,1:n_E-1) + E_p(2:n_E)*Q_inj(l_min,2:n_E))/2.d0 * dE(1:n_E-1) )
         Q_inj(l_min,:) = L_inj_p(l_min)/integral * Q_inj(l_min,:)
      end do
      !$omp end parallel do

      ! Protons injected upstream (l_min <= l) are advected without losses
      do l = 2, m_theta
         Q_inj(l,:) = Q_inj(l,:) + Q_inj(l-1,:)
      end do
      N_ptot = spread(t_cell,2,mE_p)*Q_inj

      open(newunit=iu, file=output_file('distp_approx_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(3es25.16e3)') (theta(l), E_p(i)/eV, N_ptot(l,i), i=1,mE_p)
         write(iu,*)
      end do
      close(iu)

      E_NT = matmul(N_ptot, E_p*dE) * 4.d0*mphi      ! Total NT energy in each BS "ring"
      U_NT = E_NT/dV_ring
      U_B  = B**2/(8.d0*pi)
      open(newunit=iu, file=output_file('ENTp_approx_'//shock//'.dat'))
      write(iu,'(i6,3es17.8e3)') (l, E_NT(l), U_NT(l), U_NT(l)/U_B(l), l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('distp_approx_'//shock//'_total.dat'))
      write(iu,'(2es17.8e3)') (E_p(i)/eV, sum(N_ptot(:,i))*mphi*1.d0, i=1,mE_p) ! Total distribution for 1 of the 4 quadrants
      close(iu)

      deallocate( theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, R_acel, dV_ring )
      deallocate( efi_ac_sh, alpha_sh, E_corr, t_adv, t_cell, Emax_p0, E_NT, U_NT, U_B )
      deallocate( E_p, dE, Q_inj, N_ptot )

   END subroutine CD_distp_approx_run

   !============================================================================
   ! CD_distp_approx_run (END)
   !============================================================================


   !--------------------------------------------------------------------------!
   !                                                                          !
   !                            SUBROUTINES                                   !
   !                                                                          !
   !--------------------------------------------------------------------------!

   subroutine read_shock_p(shock, theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, R_acel, dV_ring, &
                           efi_ac_sh, alpha_sh, t_adv)
      ! Reads, at each cell: the cell size dR and injected luminosity (Linj), the thermodynamical
      ! quantities (termo), the shock width, R_acel and volume (Hnorm), the acceleration efficiency
      ! (efi_ac), the injection index (alpha) and the advection time (t_conv).
      ! One read per line, so that the number of columns of each file does not matter.
      character(LEN=*), intent(in) :: shock
      real(dp), dimension(m_theta), intent(out) :: theta, rho, B, v_t, n_a, L_inj_p, dR, H_theta, &
                                                    R_acel, dV_ring, efi_ac_sh, alpha_sh, t_adv
      integer  :: l, l2, iu_L, iu_t, iu_H, iu_e, iu_a, iu_c
      real(dp) :: dum, L_inj(m_theta)

      open(newunit=iu_L, file=output_file('Linj_'//shock//'.dat'))
      open(newunit=iu_t, file=output_file('termo_'//shock//'.dat'))
      open(newunit=iu_H, file=output_file('Hnorm_'//shock//'.dat'))
      open(newunit=iu_e, file=output_file('efi_ac_'//shock//'.dat'))
      open(newunit=iu_a, file=output_file('alpha_'//shock//'.dat'))
      open(newunit=iu_c, file=output_file('t_conv_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu_L,*) dum, dum, L_inj(l), dR(l)
         read(iu_t,*) l2, rho(l), dum, v_t(l), dum, B(l), dum, theta(l) ! all in CGS units
         read(iu_H,*) dum, H_theta(l), R_acel(l), dum, dV_ring(l)
         read(iu_e,*) dum, efi_ac_sh(l)
         read(iu_a,*) dum, alpha_sh(l)
         read(iu_c,*) dum, t_adv(l)   ! Characteristic advection time
      end do
      close(iu_L)
      close(iu_t)
      close(iu_H)
      close(iu_e)
      close(iu_a)
      close(iu_c)
      dR = R0*dR               ! The file contained dR/R0, not dR [cm]
      L_inj_p = inj_p*L_inj    ! The cell size is already taken into account in L_inj
      if (shock == 'RS') then
         n_a = rho/mp/mu_sw    ! Ambient number density
      else
         n_a = rho/mp/mu_SISM
      end if
   end subroutine read_shock_p

   !=========================================================================

   subroutine Np_bin(Ep,Np_tot,ang,shock)
      ! Total distribution and distributions binned in N_bins intervals of the angle
      character(LEN=*), intent(in) :: shock
      real(dp), intent(in) :: Ep(mE_p),Np_tot(m_theta,mE_p),ang(m_theta)
      integer, parameter :: N_bins = 5
      integer  :: i, j, iu
      real(dp) :: ang_step, ang_lo(N_bins), ang_hi(N_bins), N_bin_i(N_bins,mE_p)

      ang_step = theta_max/N_bins
      ang_lo = [ -huge(1.d0), (j*ang_step, j=1,N_bins-1) ]   ! the first bin has no lower limit
      ang_hi = [ (j*ang_step, j=1,N_bins-1), theta_max ]
      ! Sum over the positions l with ang_lo(j) <= ang(l) < ang_hi(j)
      ! (do not build this mask with spread(): gfortran 13 -O3 miscompiles it)
      do j = 1, N_bins
         N_bin_i(j,:) = matmul( merge(1.d0, 0.d0, ang >= ang_lo(j) .and. ang < ang_hi(j)), Np_tot )
      end do

      open(newunit=iu, file=output_file('distp_'//shock//'_total.dat'))
      write(iu,'(2es17.8e3)') (Ep(i)/eV, sum(Np_tot(:,i))*mphi*1.d0, i=1,mE_p) ! Total distribution for 1 of the 4 quadrants
      close(iu)

      open(newunit=iu, file=output_file('distp_'//shock//'_bin.dat'))
      do j = 1, N_bins
         write(iu,'(i3,2es17.8e3)') (j, Ep(i)/eV, N_bin_i(j,i)*mphi*1.d0, i=1,mE_p)
         write(iu,*)
      end do
      close(iu)

   end subroutine Np_bin

   !=========================================================================

   elemental function E_Hillas(B, R_acel, efi)
      ! Maximum energy from t_ac = t_diff (Bohm), with t_ac = lambda_diff*efi*E/(qe*B*c)
      ! and t_diff = R_acel**2/(2*D), D = lambda_diff*E*c/(3*qe*B)
      real(dp) :: E_Hillas
      real(dp), intent(in) :: B, R_acel, efi
      E_Hillas = qe*B*R_acel*sqrt(1.5d0/efi)/lambda_diff
   end function E_Hillas

   !=========================================================================

   subroutine Biseccion_p(Bl,x3,t_coolp,t_adv,xi,Rdiff)
      ! Bisection in log10(E) for the maximum energy: t_loss(E) = t_ac(E)
      real(dp), intent(in) :: Bl,t_coolp,t_adv,xi,Rdiff
      real(dp), intent(out) :: x3
      real(dp), parameter :: err = 1.d-3  ! error in log10(Emax) < err
      real(dp) :: a_,b_,c_,fa,fc

      a_ = log10(1.2d0*Emin_p)
      b_ = log10(Emax_p)
      fa = resta_p(a_,Bl,t_coolp,t_adv,xi,Rdiff)
      if (fa < 0.1d0*err) then   ! t_ac > t_loss already at the lowest energy: no acceleration
         x3 = a_ - 1.d0          ! log10Emax < log10(1.5*Emin_p) -> cut in CD_distp_run
         return
      end if
      do while ( abs(b_ - a_) > err )
         c_ = (a_+b_)/2.d0
         fc = resta_p(c_,Bl,t_coolp,t_adv,xi,Rdiff)
         if (fc < 0.1d0*err ) then
            b_ = c_
         elseif (fa*fc < 0.d0) then
            b_ = c_
         elseif (fa*fc > 0.d0) then
            a_ = c_
            fa = fc
         end if
      end do
      x3 = b_
   end subroutine Biseccion_p

   !=========================================================================

   function resta_p(log10E,Bl,t_coolp,t_adv,xi,Rdiff)
      ! t_loss - t_ac at E = 10**log10E
      real(dp), intent(in) :: log10E,Bl,t_coolp,t_adv,xi,Rdiff
      real(dp) :: resta_p
      real(dp) :: E_l,D_diff,t_diff,t_esc,t_loss,t_ac

      E_l = 10.d0**log10E

      D_diff = Diff(E_l,Bl,'B')
      t_diff = Rdiff**2/(2.d0*D_diff)          ! Diffusion losses
      t_esc = 1.d0/(1.d0/t_adv + 1.d0/t_diff)  ! escape

      t_loss = 1.d0/(1.d0/t_coolp + 1.d0/t_esc) ! general losses

      t_ac = lambda_diff*xi*E_l/(Bl*c*qe)
      resta_p = t_loss - t_ac
   end function resta_p

   !=========================================================================

   elemental function fQp(E, E_tr, Emax, E_c, alpha, alpha_0)
      ! Injection function: E**(-alpha_0)*E_c for E < E_tr, and
      ! E**(-alpha)*exp(-E/Emax) above E_tr (E_c makes it continuous at E_tr)
      real(dp) :: fQp
      real(dp), intent(in) :: E, E_tr, Emax, E_c, alpha, alpha_0
      if (E < E_tr) then
         fQp = E**(-alpha_0)*E_c
      else
         fQp = E**(-alpha)*exp(-E/Emax)
      endif
   end function fQp

   !=====================================================================

END module CD_distp