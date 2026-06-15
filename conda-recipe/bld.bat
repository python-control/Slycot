@echo on

set "SKBUILD_CMAKE_ARGS=-GNinja;-DVERBOSE=ON;-DSLYCOT_CI=ON"

echo SKBUILD_CMAKE_ARGS is - %SKBUILD_CMAKE_ARGS% -

%PYTHON% -m pip install --no-deps --no-build-isolation -vv .
if %ERRORLEVEL% neq 0 exit 1
