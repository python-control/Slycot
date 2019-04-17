cd $RECIPE_DIR/..

# specify where CMAKE will search for lapack and blas
# needs recent cmake (conda's 3.12) and policy CMP0074 NEW
# the ${PREFIX} points to conda-build's host environment
export BLAS_ROOT=${PREFIX}
export LAPACK_ROOT=${PREFIX}

# ensure we are not building with old cmake files
rm -rf _skbuild

env

# do the build
$PYTHON setup.py install
