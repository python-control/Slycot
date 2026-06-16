@echo on

@rem -GNinja selects Ninja build generator
@rem SLYCOT_WINDOWS_CONDA_BUILD clears a debug build setting that
@rem otherwise prevents the build
set "SKBUILD_CMAKE_ARGS=-GNinja;-DSLYCOT_WINDOWS_CONDA_BUILD=ON"

%PYTHON% -m pip install --no-deps --no-build-isolation -vv .
if %ERRORLEVEL% neq 0 exit 1
