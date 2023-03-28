module test_demo
    use testdrive, only : error_type, unittest_type, new_unittest, check
    implicit none
    private
    public :: collect_demo
  contains
  
  subroutine collect_demo(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [new_unittest("demo", test_something)]
  end subroutine collect_demo

  subroutine test_something(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, 1 + 1 == 2, "Pass!")
    if (allocated(error)) return
  end subroutine test_something
end module test_demo