module TestPeriodBlockTiming
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use PeriodBlockTimingModule, only: PeriodBlockTimingType
  use ConstantsModule, only: LINELENGTH

  implicit none
  private
  public :: collect_period_block_timing

contains

  subroutine collect_period_block_timing(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("first", test_first), &
                new_unittest("last", test_last), &
                new_unittest("all", test_all), &
                new_unittest("freq", test_freq), &
                new_unittest("step", test_step) &
                ]
  end subroutine collect_period_block_timing

  subroutine test_first(error)
    type(error_type), allocatable, intent(out) :: error
    type(PeriodBlockTimingType) :: timing
    character(len=LINELENGTH) :: line

    line = "FIRST"

    call timing%init()
    call timing%read(line)

    call check(error, timing%is_active(0, .false.))
    if (allocated(error)) return

    call check(error, .not. timing%is_active(1, .false.))
    if (allocated(error)) return

  end subroutine test_first

  subroutine test_last(error)
    type(error_type), allocatable, intent(out) :: error
    type(PeriodBlockTimingType) :: timing
    character(len=LINELENGTH) :: line

    line = "LAST"

    call timing%init()
    call timing%read(line)

    call check(error, .not. timing%is_active(0, .false.))
    if (allocated(error)) return
    
    call check(error, timing%is_active(0, .true.))
    if (allocated(error)) return

  end subroutine test_last

  subroutine test_all(error)
    type(error_type), allocatable, intent(out) :: error
    type(PeriodBlockTimingType) :: timing
    character(len=LINELENGTH) :: line

    line = "ALL"

    call timing%init()
    call timing%read(line)

    call check(error, timing%is_active(0, .true.))
    if (allocated(error)) return

    call check(error, timing%is_active(0, .false.))
    if (allocated(error)) return

    call check(error, timing%is_active(1, .true.))
    if (allocated(error)) return

    call check(error, timing%is_active(1, .false.))
    if (allocated(error)) return

  end subroutine test_all

  subroutine test_freq(error)
    type(error_type), allocatable, intent(out) :: error
    type(PeriodBlockTimingType) :: timing
    character(len=LINELENGTH) :: line

    line = "FREQUENCY 2"

    call timing%init()
    call timing%read(line)

    call check(error, timing%is_active(0, .false.))
    if (allocated(error)) return

    call check(error, .not. timing%is_active(1, .false.))
    if (allocated(error)) return

    call check(error, timing%is_active(2, .false.))
    if (allocated(error)) return

    call check(error, .not. timing%is_active(3, .false.))
    if (allocated(error)) return

  end subroutine test_freq

  subroutine test_step(error)
    type(error_type), allocatable, intent(out) :: error
    type(PeriodBlockTimingType) :: timing
    character(len=LINELENGTH) :: line

    line = "STEP 1"

    call timing%init()
    call timing%read(line)

    call check(error, .not. timing%is_active(0, .false.))
    if (allocated(error)) return

    call check(error, timing%is_active(1, .false.))
    if (allocated(error)) return

    call check(error, .not. timing%is_active(2, .false.))
    if (allocated(error)) return

  end subroutine test_step

end module TestPeriodBlockTiming