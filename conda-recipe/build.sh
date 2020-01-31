cd $RECIPE_DIR/..

# specify blas vendor should be MKL
export CMAKE_EXTRA_ARGS="-DBLA_VENDOR=MKL"

# ensure we are not building with old cmake files
rm -rf _skbuild

# do the build
$PYTHON setup.py build_ext install 
