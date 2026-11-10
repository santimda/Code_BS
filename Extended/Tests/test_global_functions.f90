program test_global_functions
   ! Tests of the physical functions of global.f90 (compiled with the BD+43 test case)
   use global
   use numerical_utils, only: vector_log
   implicit none

   real(dp) :: E_ph(4001), ratio, T_bb, L_bb, x0, T_grid(2001), r_max
   integer :: ii

   ! Black body: L = 4*pi*R**2 * int c*nBB*E dE = 4*pi*R**2*sigma*T**4
   T_bb = 4.07d4
   call vector_log(1.d-4*k*T_bb, 1.d2*k*T_bb, ratio, E_ph)
   L_bb = sum( (/ (c*nBB(E_ph(ii),T_bb)*E_ph(ii)**2, ii=1,4001) /) )*log(ratio)
   call assert_close('nBB normalization', L_bb/(sigma*T_bb**4), 1.d0, 2.d-3)

   ! Khangulyan et al. 2014: F_1..F_4 -> pi**2/6 for x0 -> 0, and decrease with x0
   x0 = 1.d-14
   call assert_close('F_1(0) = pi**2/6', F_1(x0), pi**2/6.d0, 1.d-6)
   call assert_close('F_2(0) = pi**2/6', F_2(x0), pi**2/6.d0, 1.d-6)
   call assert_close('F_3(0) = pi**2/6', F_3(x0), pi**2/6.d0, 1.d-6)
   call assert_close('F_4(0) = pi**2/6', F_4(x0), pi**2/6.d0, 1.d-5)
   call assert_true('F_1..F_4 decrease with x0', F_1(1.d0) < F_1(0.1d0) .and. F_2(1.d0) < F_2(0.1d0) &
                    .and. F_3(1.d0) < F_3(0.1d0) .and. F_4(1.d0) < F_4(0.1d0))

   ! p-p cross section: zero below the threshold, ~34.3 mb at 1 TeV (Kelner et al. 2006)
   call assert_close('sigma_pp below the threshold', sigma_pp(log(0.9d0*E_th/TeV), 0.9d0*E_th), 0.d0, 0.d0)
   call assert_close('sigma_pp(1 TeV) [mb]', sigma_pp(0.d0, TeV)*1.d27, 34.3d0, 1.d-3)

   ! Cooling function: continuous at 1e5 K and 4e7 K
   call assert_close('lambda_cool continuous at 1e5 K', lambda_cool(1.00001d5)/lambda_cool(1.d5), 1.d0, 1.d-3)
   call assert_close('lambda_cool continuous at 4e7 K', lambda_cool(4.00001d7)/lambda_cool(4.d7), 1.d0, 1.5d-2)

   ! Hot layer of a radiative shock (shock_thermal): the f-f power is
   ! (lambda_ff/lambda_cool)*(1-T_ion/T)*L_inj, which must not exceed L_inj
   T_grid = (/ (T_ion*10.d0**(6.d0*(ii-1)/2000.d0), ii=1,2001) /)
   r_max = maxval( lambda_ff(T_grid)/lambda_cool(T_grid)*(1.d0 - T_ion/T_grid) )
   call assert_true('lambda_ff/lambda_cool*(1-T_ion/T) <= 1 for T >= T_ion', r_max <= 1.d0)

   print *, 'Global function tests passed.'

contains

   subroutine assert_close(name, actual, expected, allowed_error)
      ! relative error if expected /= 0, absolute otherwise
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, allowed_error
      real(dp) :: error

      error = abs(actual-expected)
      if (expected /= 0.d0) error = error/abs(expected)
      if (error > allowed_error) then
         print *, 'FAILED:', trim(name)
         print *, 'actual, expected, error =', actual, expected, error
         stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(name, condition)
      character(len=*), intent(in) :: name
      logical, intent(in) :: condition
      if (.not. condition) then
         print *, 'FAILED:', trim(name)
         stop 1
      end if
   end subroutine assert_true

end program test_global_functions
