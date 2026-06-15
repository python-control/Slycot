=================
 Building Slycot
=================

Slycot can be built on Linux, macOS, and Windows, with building on
Linux being easiest and Windows hardest.  It can *probably* be built
on other Unix-like systems that have Python, CMake, and a BLAS/LAPACK
library.

Building from the sdist
-----------------------

This section describes building Slycot using the latest `source
distribution (sdist)`_ on PyPi.

.. _`source distribution (sdist)`: https://packaging.python.org/en/latest/glossary/#term-Source-Distribution-or-sdist

Ensure you have the the following:

- Python 3.10 or later, including development files like Python.h
- CMake, and a CMake-compatible build tool like Ninja
- C compiler (e.g. gcc, MS Visual C++, clang)
- FORTRAN 77 compiler (e.g. gfortran, ifort, flang)
- BLAS/LAPACK (e.g. OpenBLAS, ATLAS, MKL)

Create and activate a virtual environment::

  python -m venv venv-slycot
  source ./venv-slycot/bin/activate # linux, macosx
  .\venv-slycot\Scripts\Activate.ps1 # Windows powershell

Install slycot with test dependencies::

  pip install --no-binary slycot[test]

Test it::

  pytest --pyargs slycot

Building from Github source
---------------------------

Instead of building from the sdist, you can get the source by cloning
the repository, with submodules::

   git clone --recurse-submodules https://github.com/python-control/Slycot.git
  
and build from the working tree::

  pip install .[test]

The test command is the same as building from the sdist::

  pytest --pyargs slycot

Non-isolated, editable build
----------------------------

When doing development once builds over and over; in that case it's faster to do a non-isolated editable build.  Run or adapt developer script `inplace-editable-build.bash`_ for that.

.. _`inplace-editable-build.bash`: ./dev-tools/inplace-editable-build.bash


Customizing the build
----------------------

To specify the C and Fortran compilers and compilation flags use the CMake and scikit-build-core customization mechanisms.  For example, on a Unix-like system one could set specific environment variables before running the `pip install` step::

  export CC=clang # C compiler
  export CFLAGS=-Os # C compilation flags
  export FC=ifort # Fortran compiler
  export FFLAGS=-O3 # Fortran compilation flags
  export BLA_VENDOR=Atlas # BLAS/LAPACK library to use
  export SKBUILD_BUILD_VERBOSE=true # verbose build

There are other ways to specify these; see CMake and scikit-build-core docs.

See `CMake`_ documentation for more in general, and `BLA_VENDOR`_ for
information about specifying the BLAS/LAPACK library.

.. _BLA_VENDOR: https://cmake.org/cmake/help/latest/module/FindBLAS.html#input-variables


Building wheels
---------------

Wheels are built with `cibuildwheel`_ via Github Actions.  A test build of wheels for the "manylinux" x86_64 target can be built locally using developer script `check-linux-cibw.bash`_ .  To use this script you'll need to have installed cibuildwheel and `Docker`_.

The wheels bundle the OpenBLAS libraries provided by `scipy-openblas-libs`_.

.. _`cibuildwheel`: https://cibuildwheel.pypa.io/
.. _`Docker`: https://www.docker.com
.. _`check-linux-cibw.bash`: ./tools/check-linux-cibw.bash
.. _`scipy-openblas-libs`: https://github.com/MacPython/openblas-libs


Building the conda recipe
-------------------------

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

Build design
------------

..
   Slycot is a multi-language package used on different operating systems
   and computer architectures.  Building Slycot for a *single* operating
   system and architecture combination is a challenge, and supporting
   mutiple platforms is even harder.

   To solve the multi-platform problem we use Python packaging tools.
   Slycot uses `scikit-build-core`_ for its `build backend`_.
   scikit-build-core, in turn, uses `CMake`_ to find necessary compilers,
   libraries, and other components, to generate build instructions, and
   to execute these build instructions.

   The sdist published to PyPI contains the Slycot and SLICOT source
   code, and the scikit-build-core and CMake build configuration needed
   to build the package.

The build uses `scikit-build-core`_ for the build backend;
scikit-build-core, in turn, uses CMake to configure and orchestrate
compilers, libraries, and other components needed for the build.

The build configuration must be:

- cross-platform, including at least Linux, macOS, and Windows 11
- able to build binary wheels for PyPI
- compatible with building conda-forge packages
- able to be built from source with user choice of compilers and BLAS/LAPACK vendor

Most scikit-build-core configuration is in ``pyproject.toml``; see section ``[tool.scikit-build]`` in that file.  Some directives are in the Github workflow files, especially those for building wheels.

The CMake configuration is in files ``CMakeLists.txt`` and
``slycot/CMakeLists.txt``.  The former is the top-level build file,
responsible for finding the necessary components for the build, and high-level configuration.

The only custom CMake build variable is ``SLYCOT_BUNDLE_OPENBLAS``,
which, if set, arranges for scipy-openblas-lib bundling (more on this
below).

``slycot/CMakeLists.txt`` is where the actual Slycot build is defined:
what source files are included, and what is linked in.

Wheel building
~~~~~~~~~~~~~~

SLICOT needs a BLAS/LAPACK library.  For wheels, this is provided by a
bundled copy of the OpenBLAS libraries.  This requires some fiddly
setup, which is handled by the ``[tool.cibuildwheel*]`` sections in
``pyproject.toml``, and also the Github Action workflow configured in
``.github/workflows/cibuildwheel.yml``.

.. _`build backend`: https://packaging.python.org/en/latest/glossary/#term-Build-Backend
.. _`scikit-build-core`: https://scikit-build-core.readthedocs.io/
.. _`CMake`: https://cmake.org/


Additional hints
----------------

Additional hints for how to install Slycot from source can be found in the
`.github`_ directory , (commands used to build and test in the GitHub Actions
CI), the `logs from the GitHub Actions`_, and the ``conda-recipe`` directory
(conda pre-requisites, install and test commands) which is included
in the source code repository.

.. _.github: https://github.com/python-control/Slycot/tree/master/.github
.. _`logs from the GitHub Actions`: https://github.com/python-control/Slycot/actions

SLICOT version
--------------

Slycot uses a patched version of SLICOT managed in this `python-control organisation repository`_.
This SLICOT source is included the PyPI sdist, and will be automatically checked out if you use the
``--recurse-submodules`` submodules command when cloning Slycot.

.. _`python-control organisation repository`: https://github.com/python-control/SLICOT-Reference
