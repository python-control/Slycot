#!/bin/bash
cd $RECIPE_DIR/..

# conda's activate-gcc_linux-64.sh sets LDFLAGS to have the directory
# entries for finding libpython, liblapack, etc., so we can't
# overwrite it
export LDFLAGS="-shared $LDFLAGS"

# conda's activate-gfortran_linux-64.sh sets F95, but
# numpy.distutils.fcompiler looks for F90
# this is hopefully a temporary hack
if [ ! -v F90 ]; then
    if [ ! -v F95 ]; then
        echo "build.sh:error: neither F90 nor F95 variable set"
        exit 1
    fi
    export F90="$F95"
fi

$PYTHON setup.py install
