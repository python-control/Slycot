# CMake file for use in conjunction with scikit-build

cmake_minimum_required(VERSION 3.11.0)

if (CMAKE_VERSION VERSION_GREATER "3.11.99")
  cmake_policy(SET CMP0074 NEW)
endif()

project(slycot VERSION ${SLYCOT_VERSION} LANGUAGES NONE)

enable_language(C)
enable_language(Fortran)

find_package(PythonLibs REQUIRED)
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "NumPy included from: ${NumPy_INCLUDE_DIR}")
message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")
message(STATUS "LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Slycot version: ${SLYCOT_VERSION}")

add_subdirectory(slycot)
