#!/bin/bash

# Demonstrate in-place, editable install
#
#  This re-uses the build directory, including CMake config, and is
#  fast when doing iterative development.
#
#  The downside is that the setup is fragile: although the build tools
#  (scikit-build-core, CMake, and the lower-level build tool like
#  ninja) are pretty good about detecting and handling changes to the
#  source, they're not perfect, and one can end up with apparently
#  impossible results that go away when doing a clean build.
#
#  For this script to work, the following must be installed:
#    - CMake,
#    - a CMake-compatible build tool, e.g., Ninja
#    - the following components that must all be discoverable by CMake:
#      - C compiler
#      - Fortran compiler
#      - BLAS and LAPACK libraries
#
#  Tested on Debian 13.

set -euo pipefail

if [ ! -f ./dev-tools/inplace-editable-build.bash ]; then
    echo "Run this from project root"
    exit 1;
fi

echo "--Build"
python -m venv venv-build
source venv-build/bin/activate

pip install --progress-bar off setuptools_scm scikit-build-core numpy pytest scipy

pip install --no-build-isolation --editable . --verbose --config-settings skbuild.build.verbose=true --config-settings skbuild.build-dir=build

pytest --pyargs slycot
