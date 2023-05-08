tests = [
  'Demo',
  'List',
  'Sim',
  'InputOutput',
  'GenericUtils',
  'prt_UtilMisc',
  'prt_CellUtil',
  'prt_MethodSubcellPollock',
  'prt_TernaryUtil',
  'prt_TernarySolveTrack'
]

test_srcs = files(
  'testdrive.F90',
  'driver.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.F90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  link_with: mf6core,
)

test('Test source modules', tester)

foreach t : tests
  test(t, tester, args: t)
endforeach