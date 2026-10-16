MODULE initialize
    use global
    use numerical_utils, only: vector_log
    implicit none

    contains

    !============================================================================
    ! initialize_allocate: allocates the CD vectors and writes the frequency
    ! grid nu.dat, shared by CD_radsyn, opac_ff and CD_abssyn
    !============================================================================
    subroutine initialize_allocate()
        real(dp) :: nu_grid(m_nu), nu_int
        integer  :: iu

        allocate( theta(m_theta), r_norm(m_theta), R_theta(m_theta) )

        call vector_log(nu_min_syn, nu_max_syn, nu_int, nu_grid)
        open(newunit=iu, file=output_file('nu.dat'))
            write(iu,'(es25.16e3)') nu_grid
        close(iu)

    end subroutine initialize_allocate

    !============================================================================
    ! initialize_deallocate: deallocates the CD vectors
    !============================================================================
    subroutine initialize_deallocate()

        deallocate( theta, r_norm, R_theta )

    end subroutine initialize_deallocate

END MODULE initialize