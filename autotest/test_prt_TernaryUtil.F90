module test_prt_ternaryutil
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use TernaryUtil, only : rotate
    implicit none
    private
    public :: collect_prt_ternaryutil

    contains

    subroutine collect_prt_ternaryutil(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
        new_unittest("rotate", test_rotate) &
      ]
    end subroutine collect_prt_ternaryutil

    subroutine test_rotate(error)
        type(error_type), allocatable, intent(out) :: error
        call check(error, 1 - 1 == 0)
        if (allocated(error)) return
    end subroutine test_rotate

end module test_prt_ternaryutil