# ensure we are not building with old cmake files
rm -rf _skbuild
rm -rf _cmake_test_compile

export LDFLAGS="$LDFLAGS -v"
if [[ "$target_platform" == osx-64 ]]; then
  export LDFLAGS="${LDFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
  export CFLAGS="${CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
  export CMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT}
fi

# Always build against netlib implementation
# https://conda-forge.org/docs/maintainer/knowledge_base.html#blas
$PYTHON setup.py build_ext install -DBLA_VENDOR=Generic
