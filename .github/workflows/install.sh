. ./.github/workflows/set_os_env.sh

#
# Install miniconda to allow quicker installation of dependencies
# See https://conda.io/docs/user-guide/tasks/use-conda-with-travis-ci.html
#

if [[ $RUNNER_OS ==  "macOS" ]]; then
    . ./.github/workflows/wget_install_miniconda.sh
fi

hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q --all
if [[ $TEST_CONDA == 1 ]]; then
    conda install conda-build;
    conda install conda-verify;
fi
conda info -a

#
# Set up a test environment for testing everything out
conda create -q -n test-environment python="$SLYCOT_PYTHON_VERSION" pip coverage nose numpy openblas
source activate test-environment


# install scikit-build
if [[ $TEST_CONDA == 0 ]]; then
    conda config --append channels conda-forge;
    conda install -c conda-forge scikit-build >=0.8.0 ;
fi
#
# Install the slycot package (two ways, to improve robustness).  For the
# conda version, need to install lapack from conda-forge (no way to specify
# this in the recipe).
# add the conda-forge channel to the config, otherwise openblas or
# lapack cannot be found in the check
# with --override-channels to make sure the locally built slycot is installed
#
if [[ $TEST_CONDA == 1 ]]; then
    conda config --append channels conda-forge;
    conda build --python "$SLYCOT_PYTHON_VERSION" conda-recipe-openblas;
    conda install conda-forge::openblas>=0.3.0;
    conda install local::slycot;
else
    CMAKE_GENERATOR="Unix Makefiles" python setup.py install;
fi
