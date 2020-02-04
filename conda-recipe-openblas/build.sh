cd $RECIPE_DIR/..

# specify blas vendor should be OpenBLAS
export BLA_VENDOR=OpenBLAS

# ensure we are not building with old cmake files
rm -rf _skbuild

# do the build
$PYTHON -m pip install . --no-deps --ignore-installed -vv

