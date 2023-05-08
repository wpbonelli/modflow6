module test_prt_utilmisc
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use UtilMiscModule, only : FirstNonBlank, GetLayerRowColumn
    use ConstantsModule, only : LINELENGTH
    implicit none
    private
    public :: collect_prt_utilmisc

    contains

    subroutine collect_prt_utilmisc(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
        new_unittest("FirstNonBlank", test_FirstNonBlank), &
        new_unittest("GetLayerRowColumn", test_GetLayerRowColumn) &
      ]
    end subroutine collect_prt_utilmisc

    subroutine test_FirstNonBlank(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=LINELENGTH) :: str

      str = "hello world"
      call check(error, FirstNonBlank(str) == 1)
      if (allocated(error)) return

      str = "    hello world"
      call check(error, FirstNonBlank(str) == 5)
      if (allocated(error)) return

      str = "    "
      call check(error, FirstNonBlank(str) == 0)
      if (allocated(error)) return

      str = ""
      call check(error, FirstNonBlank(str) == 0)
      if (allocated(error)) return
    end subroutine test_FirstNonBlank

    subroutine test_GetLayerRowColumn(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: cellNumber
      integer :: layer
      integer :: row
      integer :: column
      integer :: nlay
      integer :: nrow
      integer :: ncol

      cellNumber = 1
      nlay = 1
      nrow = 1
      ncol = 1

      call GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
      print *, layer, row, column
      call check(error, layer == 1)
      call check(error, row == 1)
      call check(error, column == 1)

      cellNumber = 3
      nlay = 3
      nrow = 1
      ncol = 1

      call GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
      print *, layer, row, column
      call check(error, layer == 3)
      call check(error, row == 1)
      call check(error, column == 1)

      cellNumber = 3
      nlay = 1
      nrow = 3
      ncol = 1

      call GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
      print *, layer, row, column
      call check(error, layer == 1)
      call check(error, row == 3)
      call check(error, column == 1)

      cellNumber = 3
      nlay = 1
      nrow = 1
      ncol = 3

      call GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
      print *, layer, row, column
      call check(error, layer == 1)
      call check(error, row == 1)
      call check(error, column == 3)

      cellNumber = 4
      nlay = 2
      nrow = 2
      ncol = 2

      call GetLayerRowColumn(cellNumber, layer, row, column, nlay, nrow, ncol)
      print *, layer, row, column
      call check(error, layer == 1)
      call check(error, row == 2)
      call check(error, column == 2)
      
    end subroutine test_GetLayerRowColumn

    ! subroutine test_transform_coords(error)
    !   type(error_type), allocatable, intent(out) :: error
    ! end subroutine test_transform_coords

end module test_prt_utilmisc