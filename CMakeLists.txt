# CMake project for MODFLOW 6 (serial build, with the FULL_3D DISU extension and VTU output).
#
# It mirrors the meson build: the source list is read from src/meson.build, so the two build systems stay in sync.
# Works with Visual Studio 2022 (open the folder, "CMake" project) using Intel ifx, and with gfortran on Linux/MSYS2.
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build --config Release
#
# Outputs: mf6 (executable) and libmf6 (shared library, BMI/XMI) in <build>/bin.
cmake_minimum_required(VERSION 3.20)
project(MODFLOW6 LANGUAGES Fortran)

option(MF6_BUILD_LIBRARY "Build the shared library libmf6 (BMI/XMI)" ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ---- source list taken from src/meson.build (block "modflow_sources = files( ... )") ----
file(READ ${CMAKE_SOURCE_DIR}/src/meson.build MESON_TEXT)
string(REGEX MATCH "modflow_sources = files\\(([^)]*)\\)" _m "${MESON_TEXT}")
set(MESON_BLOCK "${CMAKE_MATCH_1}")
string(REGEX MATCHALL "'[^']+'( */ *'[^']+')*" MESON_ENTRIES "${MESON_BLOCK}")
set(MF6_SOURCES "")
foreach(entry IN LISTS MESON_ENTRIES)
  string(REGEX REPLACE " */ *" "/" entry "${entry}")
  string(REPLACE "'" "" entry "${entry}")
  list(APPEND MF6_SOURCES ${CMAKE_SOURCE_DIR}/src/${entry})
endforeach()
list(LENGTH MF6_SOURCES MF6_NSRC)
if(MF6_NSRC LESS 400)
  message(FATAL_ERROR "Only ${MF6_NSRC} sources found in src/meson.build (parsing problem)")
endif()
message(STATUS "MODFLOW 6: ${MF6_NSRC} source files")

file(READ ${CMAKE_SOURCE_DIR}/src/meson.build _ext_text)
string(REGEX MATCH "external_libraries = files\\(([^)]*)\\)" _e "${_ext_text}")
string(REGEX MATCHALL "'[^']+'( */ *'[^']+')*" EXT_ENTRIES "${CMAKE_MATCH_1}")
set(MF6_EXTERNAL "")
foreach(entry IN LISTS EXT_ENTRIES)
  string(REGEX REPLACE " */ *" "/" entry "${entry}")
  string(REPLACE "'" "" entry "${entry}")
  list(APPEND MF6_EXTERNAL ${CMAKE_SOURCE_DIR}/src/${entry})
endforeach()

# version.f90 is generated from version.f90.in (meson does the same with vcs_tag)
set(MF6_VCS_TAG "" CACHE STRING "Suffix of the version string (for example the git commit)")
set(VCS_TAG "${MF6_VCS_TAG}")
configure_file(${CMAKE_SOURCE_DIR}/src/Utilities/version.f90.in ${CMAKE_BINARY_DIR}/generated/version.f90 @ONLY)
list(APPEND MF6_SOURCES ${CMAKE_BINARY_DIR}/generated/version.f90)

# ---- compiler options (same as meson.build) ----
set(MF6_OPTS "")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  list(APPEND MF6_OPTS -cpp -fall-intrinsics -pedantic -Wcharacter-truncation -Wno-unused-dummy-argument
       -Wno-intrinsic-shadow -Wno-maybe-uninitialized -Wno-uninitialized -std=f2018)
  if(WIN32)
    list(APPEND MF6_OPTS -D_WIN32)
  elseif(APPLE)
    list(APPEND MF6_OPTS -D__APPLE__)
  else()
    list(APPEND MF6_OPTS -D__linux__)
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if(WIN32)
    list(APPEND MF6_OPTS /fpe:0 /heap-arrays:0 /traceback /fpp /Qdiag-disable:7416 /Qdiag-disable:7025
         /Qdiag-disable:5268 /Qdiag-disable:10448)
    set(CMAKE_Fortran_FLAGS_RELEASE "/O2")
    # static runtime, as in the meson build (b_vscrt=static_from_buildtype)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  else()
    list(APPEND MF6_OPTS -fpe0 -no-heap-arrays -traceback -fpp -diag-disable:7416 -diag-disable:7025 -diag-disable:5268)
  endif()
endif()

add_library(mf6_external STATIC ${MF6_EXTERNAL})
target_compile_options(mf6_external PRIVATE ${MF6_OPTS})
set_target_properties(mf6_external PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(mf6core STATIC ${MF6_SOURCES})
target_compile_options(mf6core PRIVATE ${MF6_OPTS})
target_link_libraries(mf6core PUBLIC mf6_external)
set_target_properties(mf6core PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(mf6 src/mf6.f90)
target_compile_options(mf6 PRIVATE ${MF6_OPTS})
target_link_libraries(mf6 PRIVATE mf6core)

if(MF6_BUILD_LIBRARY)
  add_library(libmf6 SHARED
    srcbmi/bmi.f90 srcbmi/mf6bmi.f90 srcbmi/mf6bmiError.f90 srcbmi/mf6bmiGrid.f90 srcbmi/mf6bmiUtil.f90 srcbmi/mf6xmi.F90)
  target_compile_options(libmf6 PRIVATE ${MF6_OPTS})
  target_link_libraries(libmf6 PRIVATE mf6core)
  set_target_properties(libmf6 PROPERTIES OUTPUT_NAME mf6 PREFIX "lib")
endif()

install(TARGETS mf6 RUNTIME DESTINATION bin)
