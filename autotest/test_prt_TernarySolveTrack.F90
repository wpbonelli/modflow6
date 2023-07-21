module test_prt_ternarysolvetrack
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use TernarySolveTrack, only : zeroin
    implicit none
    private
    public :: collect_prt_ternarysolvetrack

    contains

    subroutine collect_prt_ternarysolvetrack(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
        new_unittest("zeroin", test_zeroin) &
      ]
    end subroutine collect_prt_ternarysolvetrack

    subroutine test_zeroin(error)
        type(error_type), allocatable, intent(out) :: error
        ! todo
        call check(error, 1 - 1 == 0)
        if (allocated(error)) return
    end subroutine test_zeroin

end module test_prt_ternarysolvetrack