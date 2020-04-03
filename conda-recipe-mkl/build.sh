cd $RECIPE_DIR/..

# specify blas vendor should be MKL
export BLA_VENDOR=Intel10_64lp

# ensure we are not building with old cmake files
rm -rf _skbuild

# do the build
$PYTHON -m pip install . --no-deps --ignore-installed -vv
