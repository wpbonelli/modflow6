@echo off
rem Builds MODFLOW 6 (serial) with Intel ifx and CMake/Ninja. Run from the repository root.
rem Needs: Intel oneAPI (ifx), CMake and Ninja (both come with Visual Studio 2022 "CMake tools").
setlocal
if exist "%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat" (
  call "%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat" intel64 vs2022
) else if exist "%ProgramFiles%\Intel\oneAPI\setvars.bat" (
  call "%ProgramFiles%\Intel\oneAPI\setvars.bat" intel64 vs2022
) else (
  echo Intel oneAPI setvars.bat not found. Install the Intel oneAPI HPC Toolkit.
  exit /b 1
)
cmake --preset ifx-release || exit /b 1
cmake --build --preset ifx-release || exit /b 1
echo.
echo Done: build\ifx-release\bin\mf6.exe
endlocal
