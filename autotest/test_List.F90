module test_List
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use ListModule, only: ListType, ListNodeType

    implicit none
    private
    public :: collect_list

    type :: IntNodeType
      integer :: value
    end type IntNodeType
  contains

  subroutine collect_list(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
        new_unittest("count", test_count), &
        new_unittest("add", test_add), &
        new_unittest("remove_node_by_index", test_remove_node_by_index) &
        ! new_unittest("remove_this_node", test_remove_this_node) &
    ]
  end subroutine collect_list

  subroutine test_count(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    type(ListType), pointer :: list => null()
    type(IntNodeType), pointer :: p1 => null()
    class(*), pointer :: p => null()

    ! allocate
    allocate(list)
    allocate(p1)

    ! empty case
    call check(error, list%Count() == 0)
    if (allocated(error)) return

    ! one node case
    p1%value = 1
    p => p1
    call list%Add(p)
    call check(error, list%Count() == 1)
    if (allocated(error)) return

    ! deallocate
    deallocate(list)
    deallocate(p1)
  end subroutine test_count

  subroutine test_add(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    type(ListType), pointer :: list => null()
    type(IntNodeType), pointer :: p1 => null()
    class(*), pointer :: p => null()

    ! allocate
    allocate(list)
    allocate(p1)

    ! add a node
    p1%value = 1
    p => p1
    call list%Add(p)
    call check(error, list%Count() == 1)
    if (allocated(error)) return

    ! deallocate
    deallocate(list)
  end subroutine test_add

  subroutine test_remove_node_by_index(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    type(ListType), pointer :: list => null()
    type(IntNodeType), pointer :: p1 => null()
    class(*), pointer :: p => null()

    ! allocate
    allocate(list)
    allocate(p1)

    ! add a node
    p1%value = 1
    p => p1
    call list%Add(p)
    call check(error, list%Count() == 1)
    if (allocated(error)) return

    ! remove a node
    call list%RemoveNode(1, destroyValue=.false.)
    call check(error, list%Count() == 0)

    ! deallocate
    deallocate(list)
    deallocate(p1)
  end subroutine test_remove_node_by_index

  ! subroutine test_remove_this_node(error)
  !   ! dummy
  !   type(error_type), allocatable, intent(out) :: error
  !   ! local
  !   type(ListType), pointer :: list => null()
  !   type(IntNodeType), pointer :: p1 => null()
  !   class(*), pointer :: p => null()

  !   ! allocate
  !   allocate(list)

  !   ! add a node
  !   p1%value = 1
  !   p => p1
  !   call list%Add(p)
  !   call check(error, list%Count() == 1)

  !   ! remove a node
  !   call list%RemoveNode(p, destroyValue=.false.)
  !   call check(error, list%Count() == 0)

  !   ! deallocate
  !   deallocate(list)
  !   
  ! end subroutine test_remove_this_node

end module test_List