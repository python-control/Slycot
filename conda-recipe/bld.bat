:: Uncoment following two lines for local test build
cd %RECIPE_DIR%
cd ..

"%PYTHON%" setup.py install --compiler=mingw32
if errorlevel 1 exit 1

:: Add more build steps here, if they are necessary.

:: See
:: http://docs.continuum.io/conda/build.html
:: for a list of environment variables that are set during the build process.
