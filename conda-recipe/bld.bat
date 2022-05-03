set BLAS_ROOT=%PREFIX%
set LAPACK_ROOT=%PREFIX%

# Keep deprecated setup.py install for now https://github.com/scikit-build/scikit-build/issues/705
"%PYTHON%" setup.py install -G "NMake Makefiles" -DBLA_VENDOR=Generic

if errorlevel 1 exit 1
