module CD_diste
   use global
   use controls, only: verbose, calc_N_bin
   use numerical_utils, only: vector_log
   use diffusion, only: Diff
   implicit none
   private

   public :: CD_diste_run, CD_diste_approx_run

   contains

   !============================================================================
   ! CD_diste_run (START)
   !============================================================================
   !  This program calculates the electron energy distribution at each position
   !  of the shock: diste_<system>_<shock>.dat, with Ne(l,E(i)).
   !============================================================================

   subroutine CD_diste_run(shock)
      character(LEN=*), intent(in) :: shock
      integer  :: l, i, l2, l_min, l_cut, m_theta_inj_max, n_E, iu
      integer  :: index1, index2, n_count, steps
      logical  :: new_count, fast
      real(dp) :: Eint, coef_IC, C_ic_st, kC_ic_st, kC_ic_IR, coef_sy, gamm, D_diff
      real(dp) :: log_coul, coef_coul, gdot_coul, Kappa_st
      real(dp) :: Ka, integral, Ne, cool, dE_approx, accum, flux, dt_bin
      real(dp) :: Q_E(mE_e)
      real(dp), allocatable :: theta(:), rho(:), T(:), B(:), v_t(:), n_a(:)
      real(dp), allocatable :: L_inj_e(:), dR(:), H_theta(:), R_acel(:), dV_ring(:)
      real(dp), allocatable :: eff_ac_sh(:), alpha_sh(:), E_corr(:)
      real(dp), allocatable :: T_IR(:), T_IR_ad(:), Kappa_IR(:), C_ic_IR(:)
      real(dp), allocatable :: t_cell(:), t_adv(:), t_adi(:), t_br(:)
      real(dp), allocatable :: log10Emax(:), Emax_e0(:), rg_e_max(:)
      real(dp), allocatable :: E_NT(:), U_NT(:), U_B(:)
      real(dp), allocatable :: E_e(:), dE(:), N_etot(:,:), E_dot(:,:)
      real(dp), allocatable :: t_coole(:,:), t_loss(:,:), t_sy(:,:), t_ic(:,:), t_ic_st(:,:), t_ic_IR(:,:)
      real(dp), allocatable :: t_diff(:,:), t_esc(:,:), t_ac(:,:), t_coul(:,:), t_ion(:,:)

      allocate( theta(m_theta), rho(m_theta), T(m_theta), B(m_theta), v_t(m_theta), n_a(m_theta) )
      allocate( L_inj_e(m_theta), dR(m_theta), H_theta(m_theta), R_acel(m_theta), dV_ring(m_theta) )
      allocate( eff_ac_sh(m_theta), alpha_sh(m_theta), E_corr(m_theta) )
      allocate( T_IR(m_theta), T_IR_ad(m_theta), Kappa_IR(m_theta), C_ic_IR(m_theta) )
      allocate( t_cell(m_theta), t_adv(m_theta), t_adi(m_theta), t_br(m_theta) )
      allocate( log10Emax(m_theta), Emax_e0(m_theta), rg_e_max(m_theta) )
      allocate( E_NT(m_theta), U_NT(m_theta), U_B(m_theta) )
      allocate( E_e(mE_e), dE(mE_e), N_etot(m_theta,mE_e), E_dot(m_theta,mE_e) )
      allocate( t_coole(m_theta,mE_e), t_loss(m_theta,mE_e), t_sy(m_theta,mE_e), t_ic(m_theta,mE_e) )
      allocate( t_ic_st(m_theta,mE_e), t_ic_IR(m_theta,mE_e), t_diff(m_theta,mE_e), t_esc(m_theta,mE_e) )
      allocate( t_ac(m_theta,mE_e), t_coul(m_theta,mE_e), t_ion(m_theta,mE_e) )

      call validate_shock(shock, 'CD_diste_run')
      if (verbose) print *, 'Calculate the non-thermal electron distribution in the '//shock


      !---------------------------------------------------------------------------------
      !  Read the cell size, injected luminosity, thermodynamical quantities, shock
      !  width, acceleration parameters and t_adv at each cell

      call read_shock(shock, theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, R_acel, dV_ring, &
                      eff_ac_sh, alpha_sh, t_adv)

      open(newunit=iu, file=output_file('L_inj_e_'//shock//'.dat'))
      write(iu,'(es17.8e3)') sum(L_inj_e)
      close(iu)

      ! Dust (IR) photon field; T_IR.dat is also read by CD_radic
      T_IR     = 27.d0*(0.2d0)**(-1.d0/6.d0)*(L_st/1.d38)**(1.d0/6.d0)*(R_theta/pc)**(-1.d0/3.d0) ! Eq.10 del Valle 2012
      !T_IR    = 100.d0                ! K2017
      T_IR_ad  = k*T_IR/mec2           ! Adimensional dust temperature
      Kappa_IR = L_IR/(16.d0*pi*sigma*T_IR**4*R_theta**2) ! assume the dust is in a shell outside the BS
      open(newunit=iu, file=output_file('T_IR.dat'))
      write(iu,'(4es25.16e3)') (theta(l), T_IR(l), T_IR_ad(l), Kappa_IR(l), l=1,m_theta)
      close(iu)

      !--------------------------------------------------------------------
      ! Build an energy vector that is common for all linear-emitters:

      call vector_log(Emin_e,Emax_e,Eint,E_e)
      dE = E_e*(Eint - 1.d0)

      !========================================================================
      ! Calculate the cooling time at all energies and at each position
      ! (and also the maximum energy)

      coef_IC  = pi*h_bar**3*c**2/(2.d0*re**2*mec2**3)
      C_ic_st = coef_IC/(T_st_ad**2)
      C_ic_IR = coef_IC/(T_IR_ad**2)

      t_cell = dR/v_t                ! Cell advection time
      t_br   = 1.d15/n_a/z_br        ! added z_br according to Padovani+2018

      ! Characteristic adiabatic losses time
      t_adi(1:m_theta-1) = 3.d0*dR(1:m_theta-1)/(v_t(1:m_theta-1)*log(rho(1:m_theta-1)/rho(2:m_theta)))
      t_adi(m_theta) = t_adi(m_theta-1)
      where (t_adi <= 0.d0 .or. t_adi > 1.d20) t_adi = 1.d20   ! To avoid numerical errors (also Inf)

      !$omp parallel do default(shared) schedule(dynamic) &
      !$omp private(l, i, Kappa_st, kC_ic_st, kC_ic_IR, coef_sy, log_coul, coef_coul, gdot_coul, gamm, D_diff)
      do l = 1, m_theta
         Kappa_st = ( R_st/(2.d0*R_theta(l)) )**2
         kC_ic_st = C_ic_st/Kappa_st
         kC_ic_IR = C_ic_IR(l)/Kappa_IR(l)
         coef_sy   = 6.d0*pi*mec2**2/(c*sigma_T*B(l)**2)   ! Valid for isotropic e- distribution

         ! Maximum energy (the cut in l is applied below, after the loop)
         call bisection(coef_sy,B(l),log10Emax(l),kC_ic_st,kC_ic_IR,T_IR_ad(l),t_br(l),t_adv(l),eff_ac_sh(l),R_acel(l))

         log_coul = 29.7d0 + log( (T(l)/1.d6)/sqrt(n_a(l)) )
         coef_coul = 3.d-14*n_a(l)*log_coul

         do i = 1, mE_e
            D_diff = Diff(E_e(i),B(l),'B')
            t_diff(l,i) = R_acel(l)**2/(2.d0*D_diff) ! Diffusion losses
            t_esc(l,i)  = 1.d0/( 1.d0/t_adv(l) + 1.d0/t_diff(l,i) ) ! escape

            t_sy(l,i)    = coef_sy/E_e(i)
            gamm         = E_e(i)/mec2
            t_ic_st(l,i) = kC_ic_st*gamm/F_iso(4.d0*gamm*T_st_ad)
            t_ic_IR(l,i) = kC_ic_IR*gamm/F_iso(4.d0*gamm*T_IR_ad(l))
            t_ic(l,i)    = 1.d0/(1.d0/t_ic_st(l,i) + 1.d0/t_ic_IR(l,i))

            gdot_coul   = coef_coul/sqrt(gamm**2-1.d0) ! Pittard2006
            t_coul(l,i) = gamm/gdot_coul

            t_ion(l,i) = (4.8d7/n_a(l)) * (E_e(i)/(1.d9*eV)) * yr ! Kornecki+2020
         end do
      end do
      !$omp end parallel do

      t_coole = 1.d0/(1.d0/t_sy + 1.d0/t_ic + 1.d0/spread(t_br,2,mE_e) + 1.d0/spread(t_adi,2,mE_e) &
                      + 1.d0/t_coul + 1.d0/t_ion)
      t_loss  = 1.d0/(1.d0/t_coole + 1.d0/t_esc)                          ! general losses
      t_ac    = lambda_diff * spread(eff_ac_sh,2,mE_e)*spread(E_e,1,m_theta)/(spread(B,2,mE_e)*c*qe)
      E_dot   = spread(E_e,1,m_theta)/t_coole

      ! bisection returns log10Emax < log10(Emin_e_inj) when Emax < Emin_e_inj: particle
      ! acceleration stops at the first such position, l_cut
      l_cut = findloc(log10Emax < log10(Emin_e_inj), .true., dim=1)
      m_theta_inj_max = m_theta              ! Default value is to inject up to m_theta
      if (l_cut > 0) then
         m_theta_inj_max = l_cut
         log10Emax(l_cut) = log10(1.001d0*Emin_e)
      end if
      Emax_e0 = MIN( 10.d0**log10Emax, E_Hillas(B, R_acel, eff_ac_sh) )
      if (l_cut > 0) Emax_e0(l_cut+1:) = 1.001d0*Emin_e   ! Dummy value, no particle acceleration
      if (any(Emax_e0 >= 0.5d0*Emax_e)) print *, 'Emax_e in global.f90 is too small! max(Emax_e0) [eV] =', maxval(Emax_e0)/eV

      rg_e_max = Emax_e0/qe/B     ! Max gyroradius

      open(newunit=iu, file=output_file('CD_tcoole_'//system//'_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(i6,16es17.8e3)') (l, E_e(i)/eV, t_ic(l,i), t_sy(l,i), t_br(l), t_diff(l,i), t_ac(l,i), &
              t_coole(l,i), t_adv(l), t_cell(l), t_ic_st(l,i), t_ic_IR(l,i), t_esc(l,i), t_loss(l,i), &
              t_adi(l), t_coul(l,i), t_ion(l,i), i=1,mE_e)
      end do
      close(iu)

      open(newunit=iu, file=output_file('Emaxe'//shock//'.dat'))
      write(iu,'(i6,es17.8e3)') (l, Emax_e0(l)/eV, l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('Gyroradius_e_'//shock//'.dat'))
      write(iu,'(i6,3es17.8e3)') (l, theta(l), rg_e_max(l), H_theta(l)/rg_e_max(l), l=1,m_theta)
      close(iu)

      !========================================================================
      ! Transport along the linear emitters
      !
      ! The electrons injected in the cell l_min with energy E_e(i) are followed
      ! cell by cell (l = l_min, l_min+1, ...) as a packet in one energy bin,
      ! index1. The flow is stationary: if F particles per unit time enter a cell,
      ! the cell contains F*(time spent there). N_etot(l,:) is the sum of the
      ! packets of all the starting cells l_min and energies E_e(i).
      !
      ! Logarithmic grid: E(i+1) = Eint*E(i), dE(i) = (Eint-1)*E(i). A particle
      ! crosses the bin i in dt_bin = dE(i)/E_dot(l,i). In each cell:
      !  - Negligible cooling: if the energy lost from l to the end of the BS is
      !    less than one bin (1 - cool > 0.99/Eint), the packet keeps its bin.
      !  - Slow cooling (t_cell < dt_bin): the packet stays one or more cells in
      !    the same bin and moves to the next lower one after "steps" cells, when
      !    prod(1 - t_cell/t_coole) < 1/Eint. N*dE is conserved when the bin changes
      !    and N scales with t_cell from cell to cell (F = N*dE/t_cell constant).
      !  - Fast cooling (t_cell > dt_bin): the particles cross several bins within
      !    the cell, so the packet is spread over them: each bin gets F*dt_bin/dE,
      !    and the last one F*(remaining time of the cell)/dE. If the lowest bin is
      !    reached before the end of the cell, the particles have cooled below
      !    Emin_e and are removed. The packet passed to the next cell is in the last
      !    bin, with N = F*t_cell/dE (same F). The position within that bin is not
      !    kept: in the next cell the cascade starts in the next lower bin (error of
      !    at most one dt_bin per cell).
      ! Diffusive escape (t_diff) is only taken into account in the injection cell.
      !========================================================================

      E_corr = E_break**(alpha_1-alpha_sh) ! continuity of Q at E_break

      N_etot = 0.d0

      !$omp parallel do default(shared) schedule(dynamic) reduction(+:N_etot) &
      !$omp private(l_min, n_E, Q_E, integral, Ka, i, Ne, index1, index2, l, l2, &
      !$omp         n_count, steps, new_count, cool, dE_approx, accum, flux, dt_bin, fast)
      Do l_min = 1, m_theta_inj_max ! After m_theta_inj_max, particles are no longer accelerated
         ! Q = Ka*fQ(E). The normalization constant Ka is obtained by the condition:
         ! integral Q(E')*E'dE' = Linj. [Linj]=erg/s, [Ka]=erg/s, [Q]=1/erg/seg
         ! The integral is solved with the trapezoid method, up to 10*Emax.

         n_E = count(E_e < 10.d0*Emax_e0(l_min))
         if (n_E < 2) cycle
         Q_E(1:n_E) = fQ(E_e(1:n_E), Emin_e_inj, E_break, Emax_e0(l_min), E_corr(l_min), alpha_sh(l_min), alpha_1)
         integral = sum( (E_e(1:n_E-1)*Q_E(1:n_E-1) + E_e(2:n_E)*Q_E(2:n_E))/2.d0 * dE(1:n_E-1) )
         Ka = L_inj_e(l_min)/integral ! Injected at the starting point of the line

         !------------------------------------------------------------------------

         do i = 1, n_E
            ! Injection cell: N(E) = Q(E)*dt [1/erg], with dt the time the particles
            ! injected at E spend in the bin (t_cell, or less if they cool or escape)
            index1 = i
            if (dE(i)/E_dot(l_min,i) >= t_cell(l_min)) then
               ! Slow cooling: the particles stay in the bin i during t_cell
               Ne = Ka*Q_E(i) * min(t_cell(l_min), 1.d0/(1.d0/t_coole(l_min,i) + 1.d0/t_diff(l_min,i)))
               N_etot(l_min,i) = N_etot(l_min,i) + Ne
               ! Packet advected to the next cell: the particles that do not escape by diffusion
               Ne = Ka*Q_E(i) * min(t_cell(l_min), t_diff(l_min,i))
            else
               ! Fast cooling: flux F of particles that do not escape by diffusion
               ! (a fraction t_coole/(t_diff + t_coole) escapes before cooling), which
               ! cross the bins i, i-1, ... staying dt_bin in each one
               flux = Ka*Q_E(i)*dE(i) * t_diff(l_min,i)/(t_diff(l_min,i) + t_coole(l_min,i))
               accum = 0.d0
               index1 = i + 1
               do while (accum < t_cell(l_min) .and. index1 > 1)
                  index1 = index1 - 1
                  dt_bin = dE(index1)/E_dot(l_min,index1)
                  N_etot(l_min,index1) = N_etot(l_min,index1) + flux*min(dt_bin, t_cell(l_min)-accum)/dE(index1)
                  accum = accum + dt_bin
               end do
               if (accum < t_cell(l_min)) cycle ! cooled below Emin_e within the cell
               Ne = flux*t_cell(l_min)/dE(index1) ! packet advected to the next cell
            end if

            ! Following cells: the packet has energy E_e(index1) <= E_e(i) in the cell l

            index2 = index1 ! bin of the packet in the previous cell, l-1
            l = l_min + 1
            n_count = 0 ! cells spent in the bin index1 (slow cooling)
            cool = 0.d0
            new_count = .true.

            do while( l <= m_theta .and. index1 > 1) ! if index1 = 1, the e- are completely cooled
               ! cool = sum of t_cell/t_coole from l to the end of the BS (relative energy
               ! lost until the end), at the energy of the packet
               if (l == l_min+1) then
                  cool = sum( t_cell(l:m_theta)/t_coole(l:m_theta,index1) )
               else
                  cool = cool - t_cell(l-1)/t_coole(l-1,index1)
               end if
               dE_approx = (1.d0-t_cell(l)/t_coole(l,index1)) ! E(after the cell)/E(before)
               fast = .false.

               if ( 1.d0-cool > 0.99d0/Eint ) then ! negligible cooling: less than one bin until the end
                  continue

               elseif ( 1.d0-cool <= 0.99d0/Eint .and. dE_approx > 1.d0/Eint ) then
                  ! Slow cooling: less than one bin per cell
                  if (new_count) then
                     accum = 1.d0
                     steps = 0 ! number of cells until E(index1) -> E(index1-1)
                     l2 = l
                     do while (accum > 1.d0/Eint .and. l2 <= m_theta)
                        steps = steps + 1
                        accum = accum*(1.d0-t_cell(l2)/t_coole(l2,index1))
                        l2 = l2 + 1
                        if (accum > 1.d0/Eint .and. l2 == m_theta) steps = m_theta ! NOT COOLING
                     end do
                     new_count = .false.
                  end if
                  n_count = n_count + 1

                  if (n_count >= steps) then
                     n_count = 0
                     index1 = index1 - 1
                     new_count = .true.
                  end if

               elseif ( dE_approx <= 1.d0/Eint ) then
                  ! Fast cooling: one bin or more per cell. The particles cross the bins
                  ! index2-1, index2-2, ... staying dt_bin in each one, until the end of
                  ! the cell (t_cell(l)) or until the lowest bin (removed: index1 = 1 ends
                  ! the loop in l)
                  fast = .true.
                  flux = Ne*dE(index2)/t_cell(l-1) ! F, particles per unit time
                  accum = 0.d0 ! time spent in the cell
                  do while (accum < t_cell(l) .and. index1 > 1)
                     index1 = index1 - 1
                     dt_bin = dE(index1)/E_dot(l,index1)
                     N_etot(l,index1) = N_etot(l,index1) + flux*min(dt_bin, t_cell(l)-accum)/dE(index1)
                     accum = accum + dt_bin
                  end do
                  Ne = flux*t_cell(l)/dE(index1) ! packet advected, in the last bin
               end if

               if (.not. fast) then
                  ! The whole packet is in the bin index1: N*dE conserved if the bin
                  ! changed, and N proportional to the time spent in the cell
                  Ne = Ne * ( dE(index2)/dE(index1) ) * (t_cell(l)/t_cell(l-1))
                  N_etot(l,index1) = N_etot(l,index1) + Ne
               end if

               index2 = index1  ! index2 keeps the value of the previous step

               l = l + 1
            end do              ! End loop in shock position
         end do                 ! End loop in electron starting energies

      end do ! End loop in starting position (l_min = # linear emitters)
      !$omp end parallel do

      !----------------------------------------------------------------------
      ! Save the electron distribution at each position
      open(newunit=iu, file=output_file('diste_'//system//'_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(3es25.16e3)') (theta(l), E_e(i)/eV, N_etot(l,i), i=1,mE_e)
         write(iu,*)
      end do
      close(iu)

      ! Non-thermal energy in each BS "ring" (4 quadrants), its density and ratio to U_B
      E_NT = matmul(N_etot, E_e*dE) * n_emit
      U_NT = E_NT/dV_ring
      U_B  = B**2/(8.d0*pi)
      open(newunit=iu, file=output_file('ENTe_'//system//'_'//shock//'.dat'))
      write(iu,'(5es17.8e3)') (theta(l), E_NT(l), U_NT(l), U_NT(l)/U_B(l),  t_cell(l), l=1,m_theta)
      close(iu)
      open(newunit=iu, file=output_file('PNTe_'//system//'_'//shock//'.dat'))
      write(iu,'(2es17.8e3)') (theta(l), (U_NT(l)/3.d0)/U_B(l), l=1,m_theta)   ! Ratio P_e,rel/P_mag
      close(iu)

      !------------------------------------------------------------
      ! Group in bins according to position angle
      if (calc_N_bin) call Ne_bin(E_e,N_etot,theta,shock)

      !------------------------------------------------------------

      deallocate( theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, R_acel, dV_ring )
      deallocate( eff_ac_sh, alpha_sh, E_corr, T_IR, T_IR_ad, Kappa_IR, C_ic_IR )
      deallocate( t_cell, t_adv, t_adi, t_br, log10Emax, Emax_e0, rg_e_max, E_NT, U_NT, U_B )
      deallocate( E_e, dE, N_etot, E_dot, t_coole, t_loss, t_sy, t_ic, t_ic_st, t_ic_IR )
      deallocate( t_diff, t_esc, t_ac, t_coul, t_ion )

      !call CD_diste_approx_run(shock) ! JUST FOR A CHECK

   END subroutine CD_diste_run

   !============================================================================
   ! CD_diste_run (END)
   !============================================================================

   !============================================================================
   ! CD_diste_approx_run (START)
   !============================================================================
   !  Same as CD_diste_run, but NEGLECTING ENERGY LOSSES (JUST FOR A CHECK).
   !  N(l,E) = t_cell(l) * sum_{l_min <= l} Q_lmin(E)
   !============================================================================

   subroutine CD_diste_approx_run(shock)
      character(LEN=*), intent(in) :: shock
      integer  :: l, i, l_min, n_E, iu
      real(dp) :: Eint, integral
      real(dp), allocatable :: theta(:), rho(:), T(:), B(:), v_t(:), n_a(:)
      real(dp), allocatable :: L_inj_e(:), dR(:), H_theta(:), R_acel(:), dV_ring(:)
      real(dp), allocatable :: eff_ac_sh(:), alpha_sh(:), E_corr(:), t_adv(:), t_cell(:)
      real(dp), allocatable :: coef_sy(:), Emax_e0(:), E_NT(:), U_NT(:), U_B(:)
      real(dp), allocatable :: E_e(:), dE(:), Q_inj(:,:), N_etot(:,:)

      allocate( theta(m_theta), rho(m_theta), T(m_theta), B(m_theta), v_t(m_theta), n_a(m_theta) )
      allocate( L_inj_e(m_theta), dR(m_theta), H_theta(m_theta), R_acel(m_theta), dV_ring(m_theta) )
      allocate( eff_ac_sh(m_theta), alpha_sh(m_theta), E_corr(m_theta), t_adv(m_theta), t_cell(m_theta) )
      allocate( coef_sy(m_theta), Emax_e0(m_theta), E_NT(m_theta), U_NT(m_theta), U_B(m_theta) )
      allocate( E_e(mE_e), dE(mE_e), Q_inj(m_theta,mE_e), N_etot(m_theta,mE_e) )

      if (verbose) print *, 'Calculate the non-thermal electron distribution in the '//shock//' (no energy losses)'


      call read_shock(shock, theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, R_acel, dV_ring, &
                      eff_ac_sh, alpha_sh, t_adv)

      call vector_log(Emin_e,Emax_e,Eint,E_e)
      dE = E_e*(Eint - 1.d0)

      t_cell = dR/v_t                                  ! Cell advection time

      ! Maximum energy: t_sy = coef_sy/E = t_ac = lambda_diff*eff_ac*E/(qe*B*c), and Hillas
      coef_sy  = 6.d0*pi*mec2**2/(c*sigma_T*B**2)        ! Valid for isotropic e- distribution
      Emax_e0 = MIN( sqrt(coef_sy*c*qe*B/(lambda_diff*eff_ac_sh)), E_Hillas(B, R_acel, eff_ac_sh) )
      if (any(Emax_e0 >= 0.5d0*Emax_e)) print *, 'Emax_e in global.f90 is too small! max(Emax_e0) [eV] =', maxval(Emax_e0)/eV

      ! Injection Q_inj(l_min,E) = Ka*fQ(E) [1/erg/s], normalized with integral Q(E')*E'dE' = Linj
      ! (trapezoid method, up to 10*Emax)
      E_corr = E_break**(alpha_1-alpha_sh)
      Q_inj  = 0.d0

      !$omp parallel do default(shared) schedule(dynamic) private(l_min, n_E, integral)
      do l_min = 1, m_theta
         n_E = count(E_e < 10.d0*Emax_e0(l_min))
         Q_inj(l_min,1:n_E) = fQ(E_e(1:n_E), Emin_e_inj, E_break, Emax_e0(l_min), E_corr(l_min), alpha_sh(l_min), alpha_1)
         integral = sum( (E_e(1:n_E-1)*Q_inj(l_min,1:n_E-1) + E_e(2:n_E)*Q_inj(l_min,2:n_E))/2.d0 * dE(1:n_E-1) )
         Q_inj(l_min,:) = L_inj_e(l_min)/integral * Q_inj(l_min,:)
      end do
      !$omp end parallel do

      ! Electrons injected upstream (l_min <= l) are advected without losses
      do l = 2, m_theta
         Q_inj(l,:) = Q_inj(l,:) + Q_inj(l-1,:)
      end do
      N_etot = spread(t_cell,2,mE_e)*Q_inj

      !----------------------------------------------------------------------
      ! Save the electron distribution at each position
      open(newunit=iu, file=output_file('diste_approx_'//shock//'.dat'))
      do l = 1, m_theta
         write(iu,'(3es25.16e3)') (theta(l), E_e(i)/eV, N_etot(l,i), i=1,mE_e)
         write(iu,*)
      end do
      close(iu)

      E_NT = matmul(N_etot, E_e*dE) * n_emit      ! Total NT energy in each BS "ring"
      U_NT = E_NT/dV_ring
      U_B  = B**2/(8.d0*pi)
      open(newunit=iu, file=output_file('ENTe_approx_'//shock//'.dat'))
      write(iu,'(i6,3es17.8e3)') (l, E_NT(l), U_NT(l), U_NT(l)/U_B(l), l=1,m_theta)
      close(iu)

      open(newunit=iu, file=output_file('diste_approx_'//shock//'_total.dat'))
      write(iu,'(2es17.8e3)') (E_e(i)/eV, sum(N_etot(:,i))*n_emit/4.d0, i=1,mE_e) ! Total distribution for 1 of the 4 quadrants
      close(iu)

      deallocate( theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, R_acel, dV_ring )
      deallocate( eff_ac_sh, alpha_sh, E_corr, t_adv, t_cell, coef_sy, Emax_e0, E_NT, U_NT, U_B )
      deallocate( E_e, dE, Q_inj, N_etot )

   END subroutine CD_diste_approx_run

   !============================================================================
   ! CD_diste_approx_run (END)
   !============================================================================


   !--------------------------------------------------------------------------!
   !                                                                          !
   !                            SUBROUTINES                                   !
   !                                                                          !
   !--------------------------------------------------------------------------!

   subroutine read_shock(shock, theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, R_acel, dV_ring, &
                         eff_ac_sh, alpha_sh, t_adv)
      ! Reads the shock quantities at each cell
      character(LEN=*), intent(in) :: shock
      real(dp), dimension(m_theta), intent(out) :: theta, rho, T, B, v_t, n_a, L_inj_e, dR, H_theta, &
                                                    R_acel, dV_ring, eff_ac_sh, alpha_sh, t_adv
      integer  :: l, l2, iu_L, iu_t, iu_H, iu_e, iu_a, iu_c
      real(dp) :: dum, L_inj(m_theta)

      open(newunit=iu_L, file=output_file('Linj_'//shock//'.dat'))
      open(newunit=iu_t, file=output_file('thermo_'//shock//'.dat'))
      open(newunit=iu_H, file=output_file('Hnorm_'//shock//'.dat'))
      open(newunit=iu_e, file=output_file('eff_ac_'//shock//'.dat'))
      open(newunit=iu_a, file=output_file('alpha_'//shock//'.dat'))
      open(newunit=iu_c, file=output_file('t_conv_'//shock//'.dat'))
      do l = 1, m_theta
         read(iu_L,*) dum, dum, L_inj(l), dR(l)
         read(iu_t,*) l2, rho(l), dum, v_t(l), T(l), B(l), dum, theta(l) ! all in CGS units
         read(iu_H,*) dum, H_theta(l), R_acel(l), dum, dV_ring(l)
         read(iu_e,*) dum, eff_ac_sh(l)
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
      L_inj_e = inj_e*L_inj    ! The cell size is already taken into account in L_inj
      if (shock == 'RS') then
         n_a = rho/mp/mu_sw    ! Ambient number density
      else
         n_a = rho/mp/mu_SISM
      end if
   end subroutine read_shock

   !=========================================================================

   subroutine Ne_bin(Ee,Ne_tot,ang,shock)
      ! Total distribution and distributions binned in N_bins intervals of the angle
      character(LEN=*), intent(in) :: shock
      real(dp), intent(in) :: Ee(mE_e),Ne_tot(m_theta,mE_e),ang(m_theta)
      integer, parameter :: N_bins = 5 ! For the plot diste_sw_bins
      integer  :: i, j, iu
      real(dp) :: ang_step, ang_lo(N_bins), ang_hi(N_bins), N_bin_i(N_bins,mE_e)

      ang_step = theta_max/N_bins
      ang_lo = [ -huge(1.d0), (j*ang_step, j=1,N_bins-1) ]   ! the first bin has no lower limit
      ang_hi = [ (j*ang_step, j=1,N_bins-1), theta_max ]
      ! Sum over the positions l with ang_lo(j) <= ang(l) < ang_hi(j)
      ! (not with spread(): gfortran -O3 bug)
      do j = 1, N_bins
         N_bin_i(j,:) = matmul( merge(1.d0, 0.d0, ang >= ang_lo(j) .and. ang < ang_hi(j)), Ne_tot )
      end do

      open(newunit=iu, file=output_file('diste_'//shock//'_'//system//'_total.dat'))
      write(iu,'(2es17.8e3)') (Ee(i)/eV, sum(Ne_tot(:,i))*n_emit/4.d0, i=1,mE_e) ! Total distribution for 1 of the 4 quadrants
      close(iu)

      open(newunit=iu, file=output_file('diste_'//shock//'_'//system//'_bin.dat'))
      do j = 1, N_bins
         write(iu,'(i3,2es17.8e3)') (j, Ee(i)/eV, N_bin_i(j,i)*n_emit/4.d0, i=1,mE_e)
         write(iu,*)
      end do
      close(iu)

   end subroutine Ne_bin

   !=========================================================================

   elemental function E_Hillas(B, R_acel, efi)
      ! Maximum energy from t_ac = t_diff (Bohm), with t_ac = lambda_diff*efi*E/(qe*B*c)
      ! and t_diff = R_acel**2/(2*D), D = lambda_diff*E*c/(3*qe*B)
      real(dp) :: E_Hillas
      real(dp), intent(in) :: B, R_acel, efi
      E_Hillas = qe*B*R_acel*sqrt(1.5d0/efi)/lambda_diff
   end function E_Hillas

   !=========================================================================

   subroutine bisection(csy,Bl,x3,c1,c2,T1,t_br,t_adv,xi,Rdiff)
      ! Bisection in log10(E) for the maximum energy: t_loss(E) = t_ac(E)
      real(dp), intent(in) :: csy,Bl,c1,c2,T1,t_br,t_adv,xi,Rdiff
      real(dp), intent(out) :: x3
      real(dp), parameter :: err = 1.d-3  ! error in log10(Emax) < err
      real(dp) :: a_,b_,c_,fa,fc

      ! Error in nth step < |a-b|/2**n
      a_ = log10(Emin_e_inj)
      b_ = log10(Emax_e)
      fa = t_balance(a_,csy,Bl,c1,c2,T1,t_br,t_adv,xi,Rdiff)
      if (fa < 0.1d0*err) then   ! t_ac > t_loss already at Emin_e_inj: no acceleration
         x3 = a_ - 1.d0          ! log10Emax < log10(Emin_e_inj) -> cut in CD_diste_run
         return
      end if
      do while ( abs(b_ - a_) > err )
         c_ = (a_+b_)/2.d0
         fc = t_balance(c_,csy,Bl,c1,c2,T1,t_br,t_adv,xi,Rdiff)
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
   end subroutine bisection

   !=========================================================================

   function t_balance(log10E,csy,Bl,c1,c2,T1,t_br,t_adv,xi,Rdiff)
      ! t_loss - t_ac at E = 10**log10E
      real(dp), intent(in) :: log10E,csy,Bl,c1,c2,T1,t_br,t_adv,xi,Rdiff
      real(dp) :: t_balance
      real(dp) :: E_l,gamm,t_sy,t_ic_st,t_ic_IR,t_ac,t_coole
      real(dp) :: D_diff,t_diff,t_loss,t_esc

      E_l = 10.d0**log10E

      D_diff = Diff(E_l,Bl,'B')
      t_diff = Rdiff**2/(2.d0*D_diff)     ! Diffusion losses
      t_esc = 1.d0/(1.d0/t_adv + 1.d0/t_diff) ! escape

      gamm = E_l/mec2
      t_sy = csy/E_l
      t_ic_st = c1*gamm/F_iso(4.d0*gamm*T_st_ad)
      t_ic_IR = c2*gamm/F_iso(4.d0*gamm*T1)
      t_coole = 1.d0/(1.d0/t_sy + 1.d0/t_ic_st + 1.d0/t_ic_IR + 1.d0/t_br) ! radiative cooling

      t_loss = 1.d0/(1.d0/t_coole + 1.d0/t_esc) ! general losses

      t_ac = lambda_diff*xi*E_l/(Bl*c*qe)  ! t_ac = lambda_diff*eff_ac*E/(qe*B*c)
      t_balance = t_loss - t_ac
   end function t_balance

   !=========================================================================

   elemental function fQ(E, Emin, E_tr, Emax, E_c, alpha, alpha_0)
      ! Injection function: E**(-alpha_0)*E_c for Emin <= E < E_tr, and
      ! E**(-alpha)*exp(-E/Emax) above E_tr (E_c makes it continuous at E_tr)
      real(dp) :: fQ
      real(dp), intent(in) :: E, Emin, E_tr, Emax, alpha, alpha_0, E_c
      if (E < Emin) then
         fQ = 0.d0 ! no injection
      else if (E < E_tr) then
         fQ = E**(-alpha_0)*E_c
      else
         fQ = E**(-alpha)*exp(-E/Emax)
      endif
   end function fQ

   !=====================================================================

END module CD_diste