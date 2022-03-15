set BLAS_ROOT=%PREFIX%
set LAPACK_ROOT=%PREFIX%

"%PYTHON%" setup.py install -G "NMake Makefiles" -DBLA_VENDOR=Generic

if errorlevel 1 exit 1
