module GlobalDataModule
  implicit none
  integer, save :: narealsp, issflg, nper
  integer, save :: mpbasUnit, disUnit, tdisUnit, gridMetaUnit, headUnit, & ! kluge note: cull unnecessary stuff
                   headuUnit, budgetUnit, traceModeUnit, binPathlineUnit
  integer, save :: inUnit, pathlineUnit, endpointUnit, timeseriesUnit, &
                   mplistUnit, mpsimUnit, traceUnit, budchkUnit, aobsUnit, logUnit
  integer, save :: particleGroupCount
  integer, save :: gridFileType
  integer, parameter :: niunit = 100
  character(len=200) :: mpnamFile, mpsimFile, mplistFile, &
                        mpbasFile, disFile, &
                        tdisFile, gridFile, headFile, &
                        budgetFile, traceFile, gridMetaFile
  integer, parameter :: levelMin = 0, levelMax = 4, &
                        levelRelMax = levelMax - levelMin
end module
