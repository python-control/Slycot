if [[ $TEST_CONDA == 0 && $RUNNER_OS != "Linux" ]]; then
    echo "Only Linux supported for non-Conda builds";
    exit 1;
fi
# from here on assume $TEST_CONDA == 0 implies $TRAVIS_OS_NAME == linux

if [[ $TEST_CONDA == 0 ]]; then
    sudo apt-get install liblapack-dev libblas-dev;
    sudo apt-get install gfortran;
    sudo apt-get install cmake;
fi
