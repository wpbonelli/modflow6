module test_arrayhandlers
  use KindModule, only: I4B
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use UtilMiscModule, only: FirstNonBlank, GetLayerRowColumn
  use ArrayHandlersModule, only: ExpandArray2D
  use ConstantsModule, only: LINELENGTH
  implicit none
  private
  public :: collect_arrayhandlers

contains

  subroutine collect_arrayhandlers(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("ExpandArray2D_int", test_ExpandArray2D_int) &
                ]
  end subroutine collect_arrayhandlers

  subroutine test_ExpandArray2D_int(error)
    type(error_type), allocatable, intent(out) :: error
    integer, allocatable :: array(:, :)

    ! allocate array
    allocate (array(2, 2))

    ! check initial array size
    call check(error, size(array, 1) == 2)
    call check(error, size(array, 2) == 2)
    if (allocated(error)) return

    ! resize array
    call ExpandArray2D(array, 2, 2)

    ! check that arrays have been resized
    call check(error, size(array, 1) == 4)
    call check(error, size(array, 2) == 4)
    if (allocated(error)) return

  end subroutine test_ExpandArray2D_int
end module test_arrayhandlers
