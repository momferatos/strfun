# CMake build for strfun.
#
# Target:
#   strfun.exe   structure functions of turbulent velocity/magnetic fields
#                (reads an ALIAKMON-style HDF5 file; needs HDF5 Fortran)
#
# Quick start:
#   cmake -B build -DHDF5_ROOT=/path/to/hdf5
#   cmake --build build -j
#   ./build/strfun.exe n nord maxpoints nfile
#
# Common options:
#   -DCMAKE_Fortran_COMPILER=gfortran   (or nvfortran)
#   -DSTRFUN_OPENMP=ON                   enable OpenMP threading   [OFF]
#   -DCMAKE_BUILD_TYPE=Release|Debug|...                          [Release]
#
# Note: the HDF5 Fortran modules must have been built with the same Fortran
# compiler used here, otherwise hdf5.mod cannot be read. Point CMake at the
# right install with -DHDF5_ROOT=<prefix>.

cmake_minimum_required(VERSION 3.12)
project(strfun LANGUAGES Fortran)

# Default to an optimised build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
               STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(STRFUN_OPENMP "Enable OpenMP threading" OFF)

find_package(HDF5 REQUIRED COMPONENTS Fortran)

add_executable(strfun strfun.f90)
set_target_properties(strfun PROPERTIES OUTPUT_NAME strfun.exe)
target_include_directories(strfun PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(strfun PRIVATE ${HDF5_Fortran_LIBRARIES})

if(STRFUN_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
  target_link_libraries(strfun PRIVATE OpenMP::OpenMP_Fortran)
endif()

message(STATUS "strfun: OpenMP=${STRFUN_OPENMP} build=${CMAKE_BUILD_TYPE}")
