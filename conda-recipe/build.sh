cd $RECIPE_DIR/..
# HACK, until scikit-build is updated to 0.7.1 or higher ...
$PYTHON -m pip install \
	https://github.com/scikit-build/scikit-build/archive/0.7.1.zip
$PYTHON setup.py install
# same HACK
$PYTHON -m pip uninstall --yes scikit-build
