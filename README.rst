Slycot
======

.. image:: https://img.shields.io/pypi/v/slycot.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://travis-ci.org/python-control/Slycot.svg?branch=master
   :target: https://travis-ci.org/python-control/Slycot

.. image:: https://coveralls.io/repos/python-control/slycot/badge.png
   :target: https://coveralls.io/r/python-control/slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov, and Sylvester equations.

Dependencies
------------

Supported Python versions are 2.7, and 3.5 and later.

Slycot depends on Numpy and, if you are installing a binary distribution,
Numpy should be the only prerequisite (though you may need LAPACK
libraries as well, depending on your particular system configuration).

If you are installing Slycot from source, you will need a FORTRAN
compiler, such as gfortran, and BLAS/LAPACK libraries. Openblas is
also supported. The build system uses skbuild (scikit-buildsystem >=
0.8.1) and cmake.

On Debian derivatives you should be able to install OpenBLAS using::

    sudo apt-get install libopenblas-dev

Additionally install cmake and install scikit-build with pip or conda.

On Mac, you will first need to install the `developer tools
<https://developer.apple.com/xcode/>`_.  You can then install gfortran using
`homebrew <http://brew.sh>`_ with::

    brew install gcc

On Windows, the BLAS and LAPACK libraries can be obtained from:

http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.4.1/Dynamic-MINGW/Win32/

Alternatively, use conda to install BLAS and LAPACK or OpenBLAS

Installing
-----------

Using pip
~~~~~~~~~

Slycot supports the pip packaging system. You must first have pip installed.

On Debian Linux based systems you can install pip with the command::

    sudo apt-get install pip

Pip can then be used to install Slycot with the command::

    pip install slycot

Note that installing with pip may or may not require having the build
dependencies installed.  There are some binary "wheels" available on PyPI,
so if those versions match with your system, you may be able to avoid
installing from source.

Using conda
~~~~~~~~~~~

Slycot can be installed via the conda package manager from the conda-forge
channel with the following command::

    conda install -c conda-forge slycot

From source
~~~~~~~~~~~

Unpack the course code to a directory of your choice,
e.g. ``/path/to/slycot_src/``, and execute::

    cd /path/to/slycot_src/
    python setup.py install

Where # is for commands that needs to be executed as root/administrator.

If you need to specify a specific compiler, set the environment
variable FC before running the install::

    # Linux/OSX:
    export FC=/path/to/my/fortran

    # Windows:
    set FC=D:\path\to\my\fortran.exe

You can also use conda to build and install slycot from source::

    conda build conda-recipe
    conda install --use-local slycot

If you prefer to use the OpenBLAS library, a conda recipe is available in
``conda-recipe-openblas``.

Additional tips for how to install slycot from source can be found in the
.travis.yml (commands used for Travis CI) and conda-recipe/ (conda
pre-requisities).

The hardest part about installing from source is getting
a working version of FORTRAN and LAPACK installed on your system and working
properly with Python. On Windows, the build system currently uses
flang, which can be installed from conda-forge. Note that flang is
incompatible with Python 2.7.

If you are using conda, you can also get working
(binary) copies of LAPACK from conda-forge using the command::

   conda install -c conda-forge lapack

Slycot will also work with the OpenBLAS libraries.

Note that in some cases you may need to set the LIBRARY_PATH environment
variable to pick up dependencies such as -lpythonN.m (where N.m is the
version of python you are using).
