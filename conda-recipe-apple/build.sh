cd $RECIPE_DIR/..

# ensure we are not building with old cmake files
rm -rf _skbuild
rm -rf _cmake_test_compile

export LDFLAGS="$LDFLAGS -v"
if [[ "$target_platform" == osx-64 ]]; then
  export LDFLAGS="${LDFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
  export CFLAGS="${CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
fi

$PYTHON setup.py build_ext install -- \
	-DNumPy_INCLUDE_DIR=${SP_DIR}/numpy/core/include \
	-DCMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT} \
	-DBLA_VENDOR=Apple
