
cd $RECIPE_DIR/..

# specify where CMAKE will search for lapack and blas
# needs recent cmake (conda's 3.12) and policy CMP0074 NEW
export BLAS_ROOT=${CONDA_PREFIX}
export LAPACK_ROOT=${CONDA_PREFIX}

$PYTHON setup.py install
