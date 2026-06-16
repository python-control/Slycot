cmake_minimum_required(VERSION 3.17.2...3.29)

if(SLYCOT_WINDOWS_CONDA_BUILD)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "")
  message(STATUS "Set CMAKE_MSVC_DEBUG_INFORMATION_FORMAT to empty string")
endif()

project(slycot LANGUAGES C Fortran)

set(SLYCOT_BUNDLE_OPENBLAS OFF CACHE BOOL "Bundle libraries from scipy-openblas32 in built wheel")

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

cmake_path(CONVERT ${F2PY_INCLUDE_DIR} TO_CMAKE_PATH_LIST F2PY_INCLUDE_DIR)

if(SLYCOT_BUNDLE_OPENBLAS)
  set(OPENBLAS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build-libs")

  find_library(
    OPENBLAS_LIBRARY_FILE
    NAMES libscipy_openblas.lib libscipy_openblas.dylib libscipy_openblas.so scipy_openblas.lib
    PATHS ${OPENBLAS_DIR}
    NO_DEFAULT_PATH
    REQUIRED)
else()
  find_package(LAPACK REQUIRED)
endif()

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory(slycot)
