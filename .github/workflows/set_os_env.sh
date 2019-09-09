echo "exporting env vars ============================="

if [[ $RUNNER_OS ==  "Linux" ]]; then
    export MINICONDA_PATH=$CONDA
    export CONDA_SCRIPT=bin         # l o
    export BASHRC=.bashrc           # l 
elif [[ $RUNNER_OS ==  "macOS" ]]; then
    export MINICONDA_PATH=$RUNNER_WORKSPACE/miniconda
    export CONDA_SCRIPT=bin         # l o
    export BASHRC=.bash_profile     # o w
elif [[ $RUNNER_OS ==  "Windows" ]]; then
    export MINICONDA_PATH=`cygpath --unix $CONDA`
    export CONDA_SCRIPT=Scripts     # w
    export BASHRC=.bash_profile     # o w
fi

export MINICONDA_SUB_PATH=$MINICONDA_PATH/$CONDA_SCRIPT
export MINICONDA_PYTEST=$MINICONDA_PATH/envs/test-environment/$CONDA_SCRIPT/py.test
#
# Make sure that fortran compiler can find conda libraries
#
export LIBRARY_PATH="$MINICONDA_PATH/envs/test-environment/lib";
export PATH="$MINICONDA_PATH:$MINICONDA_SUB_PATH:$PATH"
