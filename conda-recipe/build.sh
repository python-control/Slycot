# ensure we are not building with old cmake files
rm -rf _skbuild
rm -rf _cmake_test_compile

export SKBUILD_CONFIGURE_OPTIONS="-DBLA_VENDOR=Generic"
$PYTHON -m pip install -v .
