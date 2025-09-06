set BLAS_ROOT=%PREFIX%
set LAPACK_ROOT=%PREFIX%

set "SKBUILD_CONFIGURE_OPTIONS=-DBLA_VENDOR=Generic"
set "CMAKE_GENERATOR=NMake Makefiles"
"%PYTHON%" -m pip install -v .

if errorlevel 1 exit 1
