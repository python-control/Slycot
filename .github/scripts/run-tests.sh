#!/bin/bash

set -e

echo "::group::Slycot unit tests"
pytest -v --pyargs slycot \
       --cov=${slycot_libdir:=$(python -c "import slycot; print(slycot.__path__[0])")} \
       --cov-config=${slycot_srcdir:=$(realpath ./slycot-src)}/.coveragerc
mv .coverage ${slycot_srcdir}/.coverage.slycot
echo "::endgroup::"

echo "::group::python-control unit tests"
pushd ${python_control_srcdir:=./python-control}
# problems with the toolbar for MPL backends, not relevant to Slycot
donttest="test_root_locus_zoom or test_sisotool"
# don't care about deprecation warnings here
donttest="$donttest or test_default_deprecation"
pytest control/tests \
        -m slycot \
        --cov=$slycot_libdir \
        --cov-config=${slycot_srcdir}/.coveragerc \
        --ignore=control/tests/docstrings_test.py \
        -k "not ($donttest)"
mv .coverage ${slycot_srcdir}/.coverage.control
popd
echo "::endgroup::"

echo "::group::run slycot.test() inside interpreter"
echo 'import slycot; slycot.test()' > runtest.py
coverage run --source ${slycot_libdir} --rcfile ${slycot_srcdir}/.coveragerc runtest.py
mv .coverage ${slycot_srcdir}/.coverage.slycot-inline

echo "::group::Combine coverage"
# needs to run from within slycot source dir
cd ${slycot_srcdir}
echo "  ${slycot_libdir}" >> .coveragerc
coverage combine
coverage report
coverage xml
echo "::endgroup::"
