export LDFLAGS="$LDFLAGS -v"
if [[ "$target_platform" == osx-64 ]]; then
  export LDFLAGS="${LDFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
  export CFLAGS="${CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT}"
fi
$PYTHON setup.py build_ext install -- -DCMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT}
