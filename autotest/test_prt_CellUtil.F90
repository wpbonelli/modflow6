module test_prt_CellUtil
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use CellDefnModule, only: CellDefnType
    use CellRectModule, only: CellRectType, create_cellRect
    use CellUtilModule, only: MetricForPointInCellPolygon

    implicit none
    private
    public :: collect_prt_CellUtil

    contains

    subroutine collect_prt_CellUtil(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
          new_unittest("MetricForPointInCellPolygon", test_MetricForPointInCellPolygon) &
      ]
    end subroutine collect_prt_CellUtil

    subroutine test_MetricForPointInCellPolygon(error)
        ! dummy
        type(error_type), allocatable, intent(out) :: error
        ! local
        ! type(CellRectType), pointer :: cellRect
        ! double precision :: xpt, ypt

        call check(error, 1 + 1 == 2)

        ! todo debug
        ! call create_cellRect(cellRect)
        ! cellRect%xOrigin = 0
        ! cellRect%yOrigin = 0
        ! cellRect%yOrigin = 0
        ! cellRect%dx = 100
        ! cellRect%dy = 100
        ! cellRect%dz = 100

        ! xpt = 1.0
        ! ypt = 1.0
        ! call check(error, MetricForPointInCellPolygon(xpt, ypt, cellRect%cellDefn) .gt. 0)
        ! if (allocated(error)) return

        ! xpt = -1.0
        ! ypt = -1.0
        ! call check(error, MetricForPointInCellPolygon(xpt, ypt, cellRect%cellDefn) .lt. 0)
    
    end subroutine

end module test_prt_CellUtil