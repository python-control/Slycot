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

# https://github.com/python-control/Slycot/issues/193
if((EXISTS "${Python_INCLUDE_DIRS}/numpy")
    AND (NOT ("${Python_INCLUDE_DIRS}/numpy" EQUAL "${Python_NumPy_INCLUDE_DIRS}")))

  message(FATAL_ERROR
  "Python include directory has a numpy sub-directory,
      ${Python_INCLUDE_DIRS}/numpy,
  which is different from Numpy include directory
    ${Python_NumPy_INCLUDE_DIRS}.
  You're probably building in a virtual environment, in which case
  uninstall numpy from the base environment and try again.")

endif()

add_subdirectory(slycot)
