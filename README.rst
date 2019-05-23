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

Slycot supports Python versions 2.7 and >=3.5.

To run the compiled Slycot package, the following must be installed as
dependencies:

- Python 2.7, 3.5+
- NumPy

If you are compiling and installing Slycot from source, you will need the
following dependencies:

- Python 2.7, 3.5+
- NumPy
- scikit-build >=0.8.1
- cmake
- C compiler (e.g. gcc, MS Visual C++)
- FORTRAN compiler (e.g. gfortran, ifort, flang)
- BLAS/LAPACK (e.g. OpenBLAS, ATLAS, MKL)

There are a variety of ways to install these dependencies on different
operating systems. See the individual packages' documentation for options.

Installing
-----------

In general Slycot requires non-trivial compilation to install on a given
system. The easiest way to get started using Slycot is by installing
pre-compiled binaries. The Slycot team provides pre-compiled binaries via the
conda package manager and conda forge package hosting channel for Linux, OSX,
and Windows.

Using conda
~~~~~~~~~~~

Install Miniconda or Anaconda and then Slycot can be installed via the conda
package manager from the conda-forge channel with the following command::

    conda install -c conda-forge slycot

Using pip
~~~~~~~~~

Slycot can also be installed via the pip package manager. Install pip as per
recommendations in pip's documentation. At a minimum, Python and pip must be
installed. If a pre-complied binary (i.e. "wheel") is available it will be
installed with no need for compilation. If not, pip will attempt to compile the
package from source and thus the compilation dependencies will be required
(scikit-build, gfortran, BLAS, etc.).

Pip can then be used to install Slycot with the command::

    pip install slycot

From source
~~~~~~~~~~~

Unpack the course code to a directory of your choice,
e.g. ``/path/to/slycot_src/``

If you need to specify a specific compiler, set the environment variable FC
before running the install::

    # Linux/OSX:
    export FC=/path/to/my/fortran

    # Windows:
    set FC=D:\path\to\my\fortran.exe

To build and install execute::

    cd /path/to/slycot_src/
    python setup.py install

You can also use conda to build and install Slycot from source::

    conda build conda-recipe
    conda install --use-local slycot

If you prefer to use the OpenBLAS library, a conda recipe is available in
``conda-recipe-openblas``.

Additional tips for how to install Slycot from source can be found in the
``.travis.yml`` (commands used for Travis CI) and conda-recipe/ (conda
pre-requisites) both which are included in the source code repository.

The hardest part about installing from source is getting a working version of
FORTRAN and LAPACK installed on your system and working properly with Python.
On Windows, the build system currently uses flang, which can be installed from
conda-forge. Note that flang is incompatible with Python 2.7.

If you are using conda, you can also get working (binary) copies of LAPACK from
conda-forge using the command::

   conda install -c conda-forge lapack

Slycot will also work with the OpenBLAS libraries.

Note that in some cases you may need to set the ``LIBRARY_PATH`` environment
variable to pick up dependencies such as ``-lpythonN.m`` (where N.m is the
version of python you are using).
