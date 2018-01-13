cd $RECIPE_DIR/..
LDFLAGS="-shared" FFLAGS="-fPIC" $PYTHON setup.py install
