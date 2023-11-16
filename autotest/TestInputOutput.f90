module TestInputOutput
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use ConstantsModule, only: LINELENGTH
  ! use InputOutputModule, only: ???
  implicit none
  private
  public :: collect_inputoutput

contains

  subroutine collect_inputoutput(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("something", &
                             test_something), &
                ]
  end subroutine collect_inputoutput

  subroutine test_something(error)
    type(error_type), allocatable, intent(out) :: error

  end subroutine test_something

end module TestInputOutput
