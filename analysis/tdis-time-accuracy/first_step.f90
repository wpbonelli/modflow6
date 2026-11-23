! First-step forms from tdis_set_delt, for comparison with first_step_cmp.py.
! Run with no argument for literal inputs, or 'vol' to read them from stdin.
program first_step
  use iso_c_binding, only: c_double
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  interface
    real(c_double) function expm1(x) bind(c, name='expm1')
      import c_double
      real(c_double), value :: x
    end function
    real(c_double) function log1p(x) bind(c, name='log1p')
      import c_double
      real(c_double), value :: x
    end function
  end interface
  real(DP) :: perlen(5), tsmult(5), a, b, c
  integer :: nstp(5), i
  character(len=8) :: arg
  perlen = [1.0_DP, 1.0_DP, 1.0_DP, 100.0_DP, 10.0_DP]
  nstp = [10, 10, 3, 50, 100]
  tsmult = [1.0001_DP, 1.000001_DP, 1.0_DP + 1.0e-9_DP, 0.999_DP, 1.005_DP]
  call get_command_argument(1, arg)
  if (arg == 'vol') then ! defeat constant folding
    do i = 1, 5
      read (*, *) perlen(i), nstp(i), tsmult(i)
    end do
  end if
  do i = 1, 5
    a = perlen(i) * (1.0_DP - tsmult(i)) / (1.0_DP - tsmult(i)**nstp(i))
    b = perlen(i) * (1.0_DP - tsmult(i)) / &
        (1.0_DP - tsmult(i)**real(nstp(i), DP))
    c = perlen(i) * (tsmult(i) - 1.0_DP) / &
        expm1(nstp(i) * log1p(tsmult(i) - 1.0_DP))
    write (*, '(3(es24.16e3,1x))') a, b, c
  end do
end program
