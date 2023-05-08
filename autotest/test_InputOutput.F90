module test_InputOutput
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use ConstantsModule, only : LINELENGTH
    use InputOutputModule, only : get_node, get_ijk, get_jk
    implicit none
    private
    public :: collect_InputOutput

    contains

    subroutine collect_InputOutput(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
        new_unittest("get_node", test_get_node), &
        new_unittest("get_ijk", test_get_ijk) &
      ]
    end subroutine collect_InputOutput

    subroutine test_get_node(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ilay
      integer :: irow
      integer :: icol
      integer :: nlay
      integer :: nrow
      integer :: ncol

      ilay = 1
      irow = 1
      icol = 1
      nlay = 1
      nrow = 1
      ncol = 1

      call check(error, get_node(ilay, irow, icol, nlay, nrow, ncol) == 1)
    end subroutine test_get_node

    subroutine test_get_ijk(error)
      type(error_type), allocatable, intent(out) :: error
      ! integer :: i
      ! integer :: j
      ! integer :: k

      ! todo
      call check(error, 1 + 1 == 2)
    end subroutine test_get_ijk

end module test_InputOutput