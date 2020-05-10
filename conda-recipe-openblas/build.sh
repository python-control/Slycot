cd $RECIPE_DIR/..

# ensure we are not building with old cmake files
rm -rf _skbuild
rm -rf _cmake_test_compile

# do the build
$PYTHON setup.py build_ext install -- \
	-DNumPy_INCLUDE_DIR=${SP_DIR}/numpy/core/include \
	-DBLA_VENDOR=OpenBLAS

