module test_prt_MethodSubcellPollock
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use MethodSubcellPollockModule, only: CalculateDT
    implicit none
    private
    public :: collect_prt_MethodSubcellPollock
  contains

  subroutine collect_prt_MethodSubcellPollock(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
        new_unittest("CalculateDT", test_CalculateDT) &
    ]
  end subroutine collect_prt_MethodSubcellPollock

  subroutine test_CalculateDT(error)
    ! dummy
    type(error_type), allocatable, intent(out) :: error
    ! local
    doubleprecision :: v1, v2, dx, xL, v, dvdx, dt
    ! integer :: status

    v1 = -0.029491455779586524
    v2 = -0
    dx = 117.74869374311351
    xL = 1.0000000000000004
    v = 0.0
    dvdx = 0.0
    dt = 0.0

    ! todo when vr very close to 0, expect status 2 and very large dt
    ! status = CalculateDT(v1, v2, dx, xL, v, dvdx, dt)
    ! call check(error, status == 2)
    ! call check(error, dt .gt. 1.0d+19)

  end subroutine test_CalculateDT

end module test_prt_MethodSubcellPollock