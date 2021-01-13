Slycot
======

.. image:: https://img.shields.io/pypi/v/slycot.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://travis-ci.org/python-control/Slycot.svg?branch=master
   :target: https://travis-ci.org/python-control/Slycot

.. image:: https://github.com/python-control/Slycot/workflows/Build%20and%20Test%20Slycot/badge.svg
   :target: https://github.com/python-control/Slycot/actions

.. image:: https://coveralls.io/repos/github/python-control/Slycot/badge.svg?branch=master
   :target: https://coveralls.io/github/python-control/Slycot?branch=master

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov, and Sylvester equations.

Dependencies
------------

Slycot supports Python versions 3.6 or later.

To run the compiled Slycot package, the following must be installed as
dependencies:

- Python 3.6+
- NumPy

If you are compiling and installing Slycot from source, you will need the
following dependencies:

- 3.6+
- NumPy
- scikit-build >= 0.10.0
- CMake
- C compiler (e.g. gcc, MS Visual C++, clang)
- FORTRAN compiler (e.g. gfortran, ifort, flang)
- BLAS/LAPACK (e.g. OpenBLAS, ATLAS, MKL)

To run the Slycot unit tests and examples, you'll also need SciPy and
pytest.

There are a variety of ways to install these dependencies on different
operating systems. See the individual packages' documentation for options.

Installing
----------

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

Compiling from source
---------------------

The hardest part about installing from source is getting a working
version of FORTRAN and LAPACK (provided by OpenBLAS, MKL, etc.)
installed on your system. Depending on where you get your NumPy and SciPy
from, you will need to use a compatible LAPACK implementation. Make sure that
the correct header files are installed, and specify the environment variable
`BLA_VENDOR`_, if necessary.

.. _BLA_VENDOR: https://cmake.org/cmake/help/latest/module/FindBLAS.html#input-variables


From source without conda (Linux, macOS, Windows)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unpack the source code (or clone the git repository) to a directory of your choice,
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

From source using the conda recipe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use conda to compile and install Slycot from source. The recipe is
located in the folder ``conda-recipe`` and is intended to work for all
platforms.

The ``conda-forge`` channel provides almost all requirements to compile
Slycot with `conda-build`_, except:

- On macOS, you need the macOS SDK. See the
  `conda-build documentation for macOS`_ how to get it.
- On Windows, you need to install `Microsoft Visual C++ 14.x`_ provided e.g.
  by `Microsoft Visual Studio`_.  To build, you'll need a command shell setup
  for both conda and the Visual Studio build tools.  See `conda activation`_
  and `Microsoft Visual Studio setup`_ for information on this.

.. _conda-build: https://docs.conda.io/projects/conda-build/en/latest/resources/commands/conda-build.html
.. _conda-build documentation for macOS: https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk
.. _Microsoft Visual C++ 14.x: https://wiki.python.org/moin/WindowsCompilers
.. _Microsoft Visual Studio: https://visualstudio.microsoft.com/de/vs/
.. _conda activation: https://docs.conda.io/projects/conda/en/latest/user-guide/troubleshooting.html#windows-environment-has-not-been-activated
.. _Microsoft Visual Studio setup: https://docs.microsoft.com/en-us/cpp/build/setting-the-path-and-environment-variables-for-command-line-builds

To build and install::

    conda build -c conda-forge conda-recipe
    conda install -c conda-forge --use-local slycot

From source in a conda environment (Windows)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A similar method can be used for Linux and macOS, but is detailed here
for Windows.  This method uses conda and conda-forge to get most build
dependencies, *except* for the C compiler.

This procedure has been tested on Python 3.7 and 3.8.

1. Install `Microsoft Visual Studio`_.
2. Unpack the source code to a directory of your choice,
3. Create a command shell setup that can run the conda commands and the Visual
   Studio build tools (see above)
4. In such a command shell, within the Slycot source code directory, run the
   following commands to build and install Slycot (this example creates a
   Python 3.8 environment)::

        conda create --channel conda-forge --name build-slycot python=3.8 numpy scipy libblas=*=*netlib liblapack=*=*netlib scikit-build flang pytest
        conda activate build-slycot

        python setup.py install

Using pip
~~~~~~~~~

We publish Slycot to the Python package index, but only as a source
package, so to install using pip you'll first need to install the
build prerequisites (compilers, libraries, etc.)

If you have these build prerequisites, the command::

    pip install slycot

will download the latest release of the source code from `PyPI`_, compile, and
install Slycot into the currently configured location (virtual environment or
user site-packages).

.. _PyPI: https://pypi.org/project/slycot

Additional hints
~~~~~~~~~~~~~~~~

Additional hints for how to install Slycot from source can be found in the
`.github`_ directory , (commands used to build and test in the GitHub Actions
CI), the `logs from the GitHub Actions`_, and the ``conda-recipe`` directory
(conda pre-requisites, install and test commands) which is included
in the source code repository.

.. _.github: https://github.com/python-control/Slycot/tree/master/.github
.. _`logs from the GitHub Actions`: https://github.com/python-control/Slycot/actions

Testing
-------
To test if the installation was successful, you can run the slycot unit tests::

    pytest --pyargs slycot

You may also run the tests by calling ``slycot.test()`` from within the python
interpreter::

    import slycot
    slycot.test()

Importing ``slycot`` or running ``pytest`` without ``--pyargs slycot`` from
inside the source directory will fail, unless the compiled wrapper library has
been installed into that directory. Note that the ``[tool:pytest]`` section
in ``setup.cfg`` enforces the ``--pyargs slycot`` argument by default.
