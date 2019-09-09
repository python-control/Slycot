# Local unit tests
# TODO: replace with nose?

echo "start xvfb ====================================="
Xvfb :99 &
export DISPLAY=:99

. ./.github/workflows/set_os_env.sh

echo "source activate test-environment ==============="
source activate test-environment
echo "change directory to .. ========================="
cd ..
echo "python Slycot/runtests.py ======================"
python Slycot/runtests.py --coverage --no-build

#
# As a deeper set of tests, get test against python-control as well
#
# Additional packages required for python-control

echo "conda install scipy matplotlib ================="
conda install scipy matplotlib
# Get python-control from source and install
echo "git clone python-control ======================="
git clone --depth 1 https://github.com/python-control/python-control.git control
cd control
echo "python python-control/setup.py test  ==========="
python setup.py test
