:: Uncoment following two lines for local test build
cd %RECIPE_DIR%
cd ..

set F77=%BUILD_PREFIX%\Library\bin\flang.exe
set F90=%BUILD_PREFIX%\Library\bin\flang.exe

"%PYTHON%" setup.py build 
"%PYTHON%" setup.py install

if errorlevel 1 exit 1

:: Add more build steps here, if they are necessary.

:: See
:: http://docs.continuum.io/conda/build.html
:: for a list of environment variables that are set during the build process.
