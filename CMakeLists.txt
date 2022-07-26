# CMake file for use in conjunction with scikit-build

cmake_minimum_required(VERSION 3.14.0)

cmake_policy(SET CMP0074 NEW)

project(slycot LANGUAGES NONE)

enable_language(C)
enable_language(Fortran)

find_package(Python COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "Python headers included from: ${Python_INCLUDE_DIRS}")
message(STATUS "NumPy headers included from: ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "F2PY headers included from: ${F2PY_INCLUDE_DIRS}")
message(STATUS "LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "BLAS: ${BLAS_LIBRARIES}")

add_subdirectory(slycot)
