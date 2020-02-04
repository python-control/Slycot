:: Uncoment following two lines for local test build
cd %RECIPE_DIR%
cd ..

:: Clean old build attempts
RD /S /Q _skbuild

set FC=%BUILD_PREFIX%\Library\bin\flang.exe
set BLA_VENDOR=Intel10_64lp
:: Prefer f2py.exe, if it exists; this is provided by numpy 1.16 (and, we assume, later)
if EXIST "%PREFIX%\Scripts\f2py.exe" (
  set F2PY=%PREFIX%\Scripts\f2py.exe
) ELSE (
:: Otherwise use f2py.bat, which is provided by numpy 1.15 and earlier
  set F2PY=%PREFIX%\Scripts\f2py.bat
)

"%PYTHON%" -m pip install . --no-deps --ignore-installed -vv

if errorlevel 1 exit 1
