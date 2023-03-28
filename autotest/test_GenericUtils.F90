module test_genericutils
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use KindModule, only: DP
    use GenericUtilitiesModule, only: is_same

    implicit none
    private
    public :: collect_genericutils

  contains

  subroutine collect_genericutils(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
        new_unittest("is_same", test_is_same), &
        new_unittest("is_not_same", test_is_not_same) &
    ]
  end subroutine collect_genericutils

  subroutine test_is_same(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    real(DP) :: r1
    real(DP) :: r2

    r1 = 1.0
    r2 = 1.0
    call check(error, is_same(r1, r2))
    if (allocated(error)) return
  end subroutine test_is_same

  subroutine test_is_not_same(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    real(DP) :: r1
    real(DP) :: r2

    r1 = 0.0
    r2 = 1.0
    call check(error, (.not.(is_same(r1, r2))))
    if (allocated(error)) return
  end subroutine test_is_not_same

end module test_genericutils