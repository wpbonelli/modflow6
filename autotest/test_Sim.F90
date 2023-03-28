module test_Sim
    use testdrive, only : error_type, unittest_type, new_unittest, check
    use SimModule, only : dev_feature, store_error, store_warning, store_note, &
                          initial_message, count_errors, count_notes, count_warnings
    use ConstantsModule, only : LINELENGTH
    use VersionModule, only: IDEVELOPMODE

    implicit none
    private
    public :: collect_sim

    contains

    subroutine collect_sim(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
          new_unittest("dev_feature", test_dev_feature), &
          new_unittest("store_and_count", test_store_and_count) &
      ]
    end subroutine collect_sim

    subroutine test_dev_feature(error)
      ! dummy
      type(error_type), allocatable, intent(out) :: error
      ! local
      character(len=LINELENGTH) :: errmsg

      ! nothing to test if in release mode,
      ! calling dev_feature will terminate.
      if (IDEVELOPMODE == 0) then
        return
      end if

      ! if in develop mode check subroutine
      ! doesn't terminate (else test fails)
      call dev_feature(errmsg)

    end subroutine test_dev_feature

    subroutine test_store_and_count(error)
      ! dummy
      type(error_type), allocatable, intent(out) :: error
      ! local
      character(len=LINELENGTH) :: ntemsg
      character(len=LINELENGTH) :: wrnmsg
      character(len=LINELENGTH) :: errmsg

      ! define messages
      ntemsg = "NOTE"
      wrnmsg = "WARNING"
      errmsg = "ERROR"

      ! initialize message arrays
      call initial_message()

      ! check no messages stored
      call check(error, count_errors() == 0)
      call check(error, count_warnings() == 0)
      call check(error, count_notes() == 0)
      if (allocated(error)) return

      ! todo store a note and check that it's stored
      call store_note(ntemsg)
      call check(error, count_notes() == 1)
      if (allocated(error)) return

      ! todo store a warning and check that it's stored
      call store_warning(wrnmsg)
      call check(error, count_warnings() == 1)
      if (allocated(error)) return

      ! store an error and check that it's stored
      call store_error(errmsg, terminate=.false.)
      call check(error, count_errors() == 1)
      if (allocated(error)) return

    end subroutine test_store_and_count
end module test_Sim