:: Uncoment following two lines for local test build
cd %RECIPE_DIR%
cd ..

:: until scikit-build is updated to 0.7.1 or higher ...
"%PYTHON%" -m pip install \
	https://github.com/scikit-build/scikit-build/archive/0.7.1.zip

set F77=%BUILD_PREFIX%\Library\bin\flang.exe
set F90=%BUILD_PREFIX%\Library\bin\flang.exe
set LAPACKLIBS=lapack:blas

"%PYTHON%" setup.py build 
"%PYTHON%" setup.py install

:: remove scikit-build again, don't want to include that
"%PYTHON%" -m pip uninstall --yes scikit-build

if errorlevel 1 exit 1

:: Add more build steps here, if they are necessary.

:: See
:: http://docs.continuum.io/conda/build.html
:: for a list of environment variables that are set during the build process.
