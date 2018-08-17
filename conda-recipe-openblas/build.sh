
cd $RECIPE_DIR/..
# HACK, until scikit-build is updated to 0.7.1 or higher ...
$PYTHON -m pip install \
	https://github.com/scikit-build/scikit-build/archive/0.7.1.zip

env

# specify where CMAKE will search for lapack and blas
# needs recent cmake (conda's 3.12) and policy CMP0074 NEW
export BLAS_ROOT=${CONDA_PREFIX}
export LAPACK_ROOT=${CONDA_PREFIX}

$PYTHON setup.py install

# same HACK, remove again
$PYTHON -m pip uninstall --yes scikit-build
$PYTHON -m pip uninstall --yes packaging
$PYTHON -m pip uninstall --yes pyparsing
$PYTHON -m pip uninstall --yes setuptools
$PYTHON -m pip uninstall --yes six
$PYTHON -m pip uninstall --yes wheel
