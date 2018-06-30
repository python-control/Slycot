Slycot
======

.. image:: https://img.shields.io/pypi/v/slycot.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://travis-ci.org/python-control/slycot.svg?branch=master
   :target: https://travis-ci.org/python-control/slycot

.. image:: https://coveralls.io/repos/python-control/slycot/badge.png
   :target: https://coveralls.io/r/python-control/slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov, and Sylvester equations.


Dependencies
------------

Slycot depends on Numpy and, if you are installing a binary distribution,
Numpy should be the only prerequisite (though you may need the LAPACK
libraries as well, depending on your particular system configuration).

If you are installing Slycot from source, you will need a FORTRAN
compiler, such as gfortran, and BLAS/LAPACK libraries.

On Debian derivatives you should be able to install all the above with a
single command::

    sudo apt-get build-dep python-scipy

On Mac, you will first need to install the `developer tools
<https://developer.apple.com/xcode/>`_.  You can then install gfortran using
`homebrew <http://brew.sh>`_ with::

    brew install gcc

On Windows, the BLAS and LAPACK libraries can be obtained from: 

http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.4.1/Dynamic-MINGW/Win32/


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

If the build fails and you are on a 64bit OS you may want to try::

    python setup.py config_fc --arch="-march=x86-64" build
    python setup.py install

You can also use conda to build and install slycot from source::

    conda build conda-recipe
    conda install --use-local slycot

If you prefer to use the OpenBLAS library, a conda recipe is available in
``conda-recipe-openblas``.

Additional tips for how to install slycot from source can be found in the
.travis.yml (commands used for Travis CI) and conda-recipe/ (conda
pre-requisities).  The hardest part about installing from source is getting
a working version of FORTRAN and LAPACK installed on your system and working
properly with Python.  If you are using conda, you can also get working
(binary) copies of LAPACK from conda-forge using the command::

	conda install -c conda-forge lapack

Slycot will also work with the OpenBLAS libraries.

Note that in some cases you may need to set the LIBRARY_PATH environment
variable to pick up dependencies such as -lpythonN.m (where N.m is the
version of python you are using).
