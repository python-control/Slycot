:: Uncoment following two lines for local test build
cd %RECIPE_DIR%
cd ..

:: indicating fortran compiler is essential
set FC=%BUILD_PREFIX%\Library\bin\flang.exe

:: The batch file created by conda-build sets a load of environment variables
:: Building worked fine without conda; apparently one or more of these 
:: variables produce test & link failures. Resetting most of these here
set ARCH=
set BUILD=
set BUILD_PREFIX=
set CMAKE_GENERATOR=
set CommandPromptType=
set CPU_COUNT=
set DISTUTILS_USE_SDK=
set folder=
set cpu_optimization_target=
set fortran_compiler=
set Framework40Version=
set FrameworkDir=
set FrameworkDIR64=
set FrameworkVersion=
set FrameworkVersion64=
set ignore_build_only_deps=
set CFLAGS=
set CXXFLAGS=
set cxx_compiler=
set c_compiler=
set INCLUDE=
set LDFLAGS_SHARED=
set LIBPATH=
set LIB=;%LIB%
set MSSdk=
set MSYS2_ARG_CONV_EXCL=
set MSYS2_ENV_CONV_EXCL=
set NETFSXDIR=
set PIP_IGNORE_INSTALLED=
set platform=
set WindowsLibPath=
set WindowsSdkDir=
set CYGWIN_PREFIX=
set SRC_DIR=
set STDLIB_DIR=
set SUBDIR=
set SYS_PREFIX=
set target_platform=
set UCRTVersion=
set UniversalCRTSdkDir=
set VCINSTALLDIR=
set vc=
set win=
set VisualStudioVersion=
set VSINSTALLDIR=
set VSREGKEY=
set VS_MAJOR=
set VS_VERSION=
set VS_YEAR=
set WindowsSDKLibVersion=
set WindowsSDKVersion=
set WindowsSDKExecutablePath_x64=
set WindowsSDKExecutablePath_x86=

:: information on remaining variables
set

set BLAS_ROOT=%CONDA_PREFIX%
set LAPACK_ROOT=%CONDA_PREFIX%

"%PYTHON%" setup.py install

if errorlevel 1 exit 1

:: Add more build steps here, if they are necessary.

:: See
:: https://conda.io/docs/user-guide/tasks/build-packages/environment-variables.html
:: for a list of environment variables that are set during the build process.
