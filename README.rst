Slycot
======

.. image:: https://img.shields.io/pypi/v/slycot.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://travis-ci.org/python-control/Slycot.svg?branch=master
   :target: https://travis-ci.org/python-control/Slycot

.. image:: https://coveralls.io/repos/github/python-control/Slycot/badge.svg?branch=master
   :target: https://coveralls.io/github/python-control/Slycot?branch=master

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov, and Sylvester equations.

Dependencies
------------

Slycot supports Python versions 2.7, and 3.5 or later.

To run the compiled Slycot package, the following must be installed as
dependencies:

- Python 2.7, 3.5+
- NumPy

If you are compiling and installing Slycot from source, you will need the
following dependencies:

- Python 2.7, 3.5+
- NumPy
- scikit-build >= 0.10.0
- CMake
- C compiler (e.g. gcc, MS Visual C++)
- FORTRAN compiler (e.g. gfortran, ifort, flang)
- BLAS/LAPACK (e.g. OpenBLAS, ATLAS, MKL)

To run the Slycot unit tests and examples, you'll also need scipy and
pytest.

There are a variety of ways to install these dependencies on different
operating systems. See the individual packages' documentation for options.

Installing
-----------

The easiest way to get started with Slycot is to install pre-compiled
binaries from conda-forge (see below); these are available for Linux,
OSX, and Windows.

Compiling the Slycot source is unfortunately a bit tricky, especially
on Windows, but we give some pointers further below for doing this.

Using conda and conda-forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First install Miniconda or Anaconda.  Slycot can then be installed
from the conda-forge channel with the following command::

    conda install -c conda-forge slycot

From source without conda (Linux, macOS, Windows)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unpack the source code to a directory of your choice,
e.g. ``/path/to/slycot_src/``

If you need to specify a specific compiler, set the environment variable FC
before running the install::

    # Linux/OSX:
    export FC=/path/to/my/fortran

    # Windows:
    set FC=D:\path\to\my\fortran.exe

To build and install, execute::

    cd /path/to/slycot_src/
    python setup.py install

From source using a conda recipe (Linux and macOS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also use conda to build and install Slycot from source, but
you'll have to choose the right recipe directory.

On Linux you can choose between ``conda-recipe-openblas`` and
``conda-recipe-mkl``

On macOS you should use ``conda-recipe-apple``. See the
`conda-build documentation`_ how to get the required macOS SDK.

.. _conda-build documentation: https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk

For example, to build with the OpenBLAS recipe::

    conda build -c conda-forge conda-recipe-openblas
    conda install -c conda-forge --use-local slycot

From source in a conda environment (Windows)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A similar method can be used for Linux and macOS, but is detailed here
for Windows.  This method uses conda and conda-forge to get most build
dependencies, *except* for the C compiler.

This procedure has been tested on Python 3.7 and 3.8.

First, install the `correct Visual Studio compiler for the Python
version`_ you wish to build for.

.. _correct Visual Studio compiler for the Python version: https://wiki.python.org/moin/WindowsCompilers

To build, you'll need a command shell setup for both conda and the
Visual Studio build tools.  See `conda activation`_ and `Microsoft
Visual Studio setup`_ for information on this.

.. _conda activation: https://docs.conda.io/projects/conda/en/latest/user-guide/troubleshooting.html#windows-environment-has-not-been-activated
.. _Microsoft Visual Studio setup: https://docs.microsoft.com/en-us/cpp/build/setting-the-path-and-environment-variables-for-command-line-builds?view=vs-2019

In such a command shell, run the following commands to build and
install Slycot (this example creates a Python 3.8 environment)::

    conda create --channel conda-forge --name build-slycot python=3.8 numpy scipy libblas=*=*netlib liblapack=*=*netlib scikit-build flang pytest
    conda activate build-slycot

    python setup.py install

Testing
~~~~~~~
To test if the installation was successful, you can run the slycot unit tests::

    pytest --pyargs slycot

Running ``pytest`` without ``--pyargs slycot`` from inside the source directory
will fail, unless either ``setup.cfg`` or the compiled wrapper library have
been installed into that directory.

General notes on compiling
~~~~~~~~~~~~~~~~~~~~~~~~~~

Additional tips for how to install Slycot from source can be found in the
``.travis.yml`` (commands used for Travis CI) and the ``conda-recipe-*/``
directories (conda pre-requisites) both which are included in the source
code repository.

The hardest part about installing from source is getting a working
version of FORTRAN and LAPACK (provided by OpenBLAS, MKL, etc.)
installed on your system, and working properly with Python.

Note that in some cases you may need to set the ``LIBRARY_PATH`` environment
variable to pick up dependencies such as ``-lpythonN.m`` (where N.m is the
version of python you are using).

Using pip
~~~~~~~~~

We publish Slycot to the Python package index, but only as a source
package, so to install using pip you'll first need to install the
build prerequisites (compilers, libraries, etc.)

If you have these build prerequisites, install in the standard way with:

    pip install slycot
