program tester
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
    & select_suite, run_selected, get_argument
  use test_Demo, only: collect_Demo
  use test_List, only: collect_List
  use test_Sim, only: collect_Sim
  use test_InputOutput, only: collect_InputOutput
  use test_GenericUtils, only: collect_GenericUtils
  use test_ArrayHandlers, only: collect_ArrayHandlers
  use test_prt_CellUtil, only: collect_prt_CellUtil
  use test_prt_MethodSubcellPollock, only: collect_prt_MethodSubcellPollock
  use test_prt_particle, only: collect_prt_Particle
  use test_prt_TernaryUtil, only: collect_prt_TernaryUtil
  use test_prt_TernarySolveTrack, only: collect_prt_TernarySolveTrack
  use test_prt_UtilMisc, only: collect_prt_UtilMisc
  use test_prt_TrackData, only: collect_prt_TrackData
  implicit none
  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0
  testsuites = [ &
               new_testsuite("Demo", collect_Demo), &
               new_testsuite("List", collect_List), &
               new_testsuite("Sim", collect_Sim), &
               new_testsuite("GenericUtils", collect_GenericUtils), &
               new_testsuite("InputOutput", collect_InputOutput), &
               new_testsuite("ArrayHandlers", collect_ArrayHandlers), &
               new_testsuite("prt_CellUtil", collect_prt_CellUtil), &
               new_testsuite("prt_MethodSubcellPollock", &
                             collect_prt_MethodSubcellPollock), &
               new_testsuite("prt_Particle", collect_prt_Particle), &
               new_testsuite("prt_TernaryUtil", collect_prt_TernaryUtil), &
               new_testsuite("prt_TernarySolveTrack", &
                             collect_prt_TernarySolveTrack), &
               new_testsuite("prt_UtilMisc", collect_prt_UtilMisc), &
               new_testsuite("prt_TrackData", collect_prt_TrackData) &
               ]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        write (error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
        if (stat < 0) then
          error stop 1
        end if
      else
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end if
    else
      write (error_unit, fmt) "Available testsuites"
      do is = 1, size(testsuites)
        write (error_unit, fmt) "-", testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

end program tester
