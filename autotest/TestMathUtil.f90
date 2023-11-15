module TestMathUtil
  use KindModule, only: I4B, DP
  use testdrive, only: check, error_type, new_unittest, test_failed, &
                       to_string, unittest_type
  use MathUtilModule, only: is_same, is_close
  use ConstantsModule, only: LINELENGTH
  implicit none
  private
  public :: collect_mathutil

contains

  subroutine collect_mathutil(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                ! old version (is_same)
                new_unittest("is_same", test_is_same), &
                new_unittest("is_same_both_near_0", &
                             test_is_same_both_near_0, &
                             should_fail=.true.), &
                new_unittest("is_not_same", test_is_not_same), &
                ! new version (is_close)
                new_unittest("is_close", test_is_close), &
                new_unittest("is_close_both_near_0", &
                             test_is_close_both_near_0) &
                ]
  end subroutine collect_mathutil

  ! old version (is_same)

  subroutine test_is_same(error)
    type(error_type), allocatable, intent(out) :: error

    ! exact
    call check(error, is_same(1.0_DP, 1.0_DP))
    if (allocated(error)) return

    ! inexact (within tolerance)
    call check(error, is_same(1.0000_DP, 1.0001_DP, eps=0.01_DP))
    if (allocated(error)) return
  end subroutine test_is_same

  subroutine test_is_same_both_near_0(error)
    type(error_type), allocatable, intent(out) :: error

    ! relative comparison mode fails when a and b are close to 0
    call check(error, is_same(0.0000_DP, 0.0001_DP, eps=0.01_DP))
    if (allocated(error)) return
  end subroutine test_is_same_both_near_0

  subroutine test_is_not_same(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, (.not. (is_same(1.0_DP, 1.0001_DP))))
    if (allocated(error)) return

    ! with tolerance
    call check(error, (.not. is_same(1.0_DP, 1.0001_DP, eps=0.00005_DP)))
    if (allocated(error)) return
  end subroutine test_is_not_same

  ! new version (is_close)

  subroutine test_is_close(error)
    type(error_type), allocatable, intent(out) :: error

    ! exact match
    call check(error, is_close(1.0_DP, 1.0_DP), &
               "exact match failed")
    if (allocated(error)) return

    ! inexact match (within default relative tolerance)
    call check(error, is_close(1.0000_DP, 1.000001_DP), &
               "inexact match failed")
    if (allocated(error)) return

    ! mismatch (beyond default relative tolerance)
    call check(error, (.not. is_close(1.0000_DP, 1.0001_DP)), &
               "default rtol mismatch failed")
    if (allocated(error)) return

    ! mismatch (beyond absolute tolerance)
    call check(error, (.not. is_close(1.0000_DP, 1.01_DP, atol=1d-3)), &
               "atol mismatch failed")
    if (allocated(error)) return
  end subroutine test_is_close

  subroutine test_is_close_both_near_0(error)
    type(error_type), allocatable, intent(out) :: error

    ! exact match
    call check(error, is_close(0.0_DP, 0.0_DP), &
               "exact match failed")
    if (allocated(error)) return

    ! inexact mismatch (inappropriate atol)
    call check(error, (.not. is_close(0.0_DP, 0.000001_DP)), &
               "inexact mismatch failed")
    if (allocated(error)) return

    ! inexact mismatch (inappropriate atol)
    call check(error, is_close(0.0_DP, 0.000001_DP, atol=1d-5), &
               "inexact match failed")
    if (allocated(error)) return
  end subroutine test_is_close_both_near_0

end module TestMathUtil
