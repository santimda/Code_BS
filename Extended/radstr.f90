module radstr
    use global
    use diffusion
    implicit none

    public :: radstr_run

    contains

    !========================================================================================!
    ! This progrma estimates the luminosity emitted via p-p and NT Br. inside the H II region!
    ! We concentrate on particles with E = E_break                                           !
    !========================================================================================!

    subroutine radstr_run(shock)
        character(LEN=*), intent(in) :: shock
        integer :: i, i_E
        real(dp) :: R_str
        real(dp) :: L_inj_p_BS, L_inj_p, Lpp_sh_tot
        real(dp) :: L_inj_e_BS, L_inj_e, Lbr_sh_tot, LIC_sh_tot, Lsy_sh_tot
        real(dp) :: L_pp_Etran, L_br_Etran, L_IC_Etran, L_sy_Etran
        real(dp) :: E_pp(mE_f), E_br(mE_f)
        real(dp) :: Lpp_sh(mE_f), Lbr_sh(mE_f)
        real(dp) :: t_diff,t_adv,t_esc,t_pp,t_br,t_IC_CMB,t_sy
        real(dp) :: n_p
        real(dp) :: D_diff
        real(dp), parameter :: B_str = 10.d-6
        character(len=*), parameter :: diff_regime = 'G'

        90    format(A8,X,Es10.2)
        91    format(A7,X,Es10.2)
        81    format(A13,X,Es10.2,X,Es10.2)
        82    format(A15,X,Es10.2)
        83    format(A13,X,Es10.2)
        84    format(A14,X,Es10.2)

        open(100,file=output_file('Stromgren.dat'))
            read(100,*) R_str, dummy
        close(100)

        t_diff = R_str**2/(6.d0*Diff(E_break,B_str,diff_regime))
        t_adv = R_str/V_st
        t_esc = 1./( 1./t_adv + 1.d0/t_diff )
        write(*,90) 'D_diff =', Diff(E_break,B_str,diff_regime)
        write(*,90) 't_diff =', t_diff
        write(*,91) 't_adv =', t_adv
        write(*,91) 't_esc =', t_esc

        n_p = (mu/mu_i)*n_ISM
        print*, 'n_p=', n_p

        !---------------------------------------------------
        !---------------------------------------------------

        
        open(30,file=output_file'L_inj_p_'//shock//'.dat'))
            read(30,*) L_inj_p_BS
        close(30)

        !---------------------------------------------------
        open(31,file='radpp_'//shock//'_tot.dat')
        do i_E = 1, mE_f
            read(31,*) E_pp(i_E), Lpp_sh(i_E)
        end do
        close(31)

        Lpp_sh_tot = 0.d0
        do i_E = 1, mE_f
            Lpp_sh_tot = Lpp_sh_tot + Lpp_sh(i_E)
        end do
        !---------------------------------------------------

        L_inj_p = 4.d0*mphi*L_inj_p_BS - Lpp_sh_tot

        write(*,81) 'L_inj_p_HII =', L_inj_p, 4.d0*mphi*L_inj_p_BS

        !---------------------------------------------------
        
        t_pp = 1.d15/n_p

        ! We estimate that approx 0.1*L_pp is radiated by E_break
        L_pp_Etran = 1.d-1*(t_esc/t_pp)*L_inj_p

        ! Protons with E = E_break radiate mostly at sim 0.1E_break
        write(*,82) 'L_pp(100 MeV) =', L_pp_Etran

        open(50,file='rad_pp_str.dat')
            write(50,*) 0.1d0*E_break/eV, L_pp_Etran
        close(50) 

        !---------------------------------------------------
        !---------------------------------------------------


        open(40,file=output_file('L_inj_e_'//shock//'.dat'))
            read(40,*) L_inj_e_BS
        close(40)

        !---------------------------------------------------
        open(41,file='radbr_'//shock//'_tot.dat')
        do i_E = 1, mE_f
            read(41,*) E_br(i_E), Lbr_sh(mE_f)
        end do
        close(41)

        Lbr_sh_tot = 0.d0
        do i = 1, mE_f
            Lbr_sh_tot = Lbr_sh_tot + Lbr_sh(mE_f)
        end do

        !---------------------------------------------------
        open(42,file='rad_IC_st_tot_'//shock//'.dat')
            read(42,*) LIC_sh_tot
        close(42)

        !---------------------------------------------------
        open(43,file=output_file('radsyn_tot_'//shock//'.dat'))
            read(43,*) Lsy_sh_tot
        close(43)

        L_inj_e = 4.d0*mphi*L_inj_e_BS - (Lbr_sh_tot + LIC_sh_tot + Lsy_sh_tot)

        write(*,81) 'L_inj_e_HII =', L_inj_e, 4.d0*mphi*L_inj_e_BS

        !---------------------------------------------------

        !t_sy = E_break*(4.d0/3.d0*c*sigma_T*(B_str**2/8.d0/pi)*(E_break/mec2)**2)**(-1.)

        ! We estimate that approx 0.1*L_sy is radiated by E_break
        !L_sy_Etran = 0.1d0*(t_esc/t_sy)*L_inj_e

        ! Electrons with E = E_break radiate mostly at sim E_break
        !write(*,84) 'L_sy(75 MHz) =', L_sy_Etran

        !---------------------------------------------------
        
        t_br = 1.d15/n_p
        write(*,91) 't_Br. =', t_br

        ! We estimate that approx 0.1*L_Br is radiated by E_break
        L_br_Etran = (t_esc/t_br)*L_inj_e

        ! Electrons with E = E_break radiate mostly at sim E_break
        write(*,83) 'L_br(1 GeV) =', L_br_Etran

        open(60,file='rad_Br_str.dat')
            write(60,*) E_break/eV, L_br_Etran
        close(60)

        !---------------------------------------------------
        !t_IC_CMB = E_break*(4.d0/3.d0*c*sigma_T*U_CMB*(E_break/mec2)**2)**(-1.)

        ! We estimate that approx 0.1*L_IC_CMB is radiated by E_break
        !L_IC_Etran = 0.1d0*(t_esc/t_IC_CMB)*L_inj_e

        ! Electrons with E = E_break radiate mostly at sim E_cmb*(E_break/mc2)**2
        !write(*,84) 'L_IC(40 keV) =', L_IC_Etran

    end subroutine radstr_run


END module radstr