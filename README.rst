Slycot
======

.. image:: https://img.shields.io/pypi/v/nine.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://travis-ci.org/python-control/Slycot.svg?branch=master
   :target: https://travis-ci.org/python-control/Slycot

.. image:: https://coveralls.io/repos/python-control/Slycot/badge.png
   :target: https://coveralls.io/r/python-control/Slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov and Sylvester equations.

Dependencies
------------

Slycot depends primarily on NumPy, and if you are installing a binary
distribution, NumPy is the only dependency.

If you are installing Slycot from source, you will need a Fortran compiler,
such as gfortran, and BLAS/LAPACK libraries.

On Debian Linux operating system derivates you can install all the necessary
build dependencies with a single command::

   sudo apt-get build-dep python-scipy

On Mac, you will first need to install the `developer tools
<https://developer.apple.com/xcode/>`_. You can then install gfortran using
`homebrew <http://brew.sh>`_ with::

   brew install gcc

On Windows, we suggest installing on top of the Python(x,y) distribution, and
grabbing BLAS and LAPACK libraries from:

http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.4.1/Dynamic-MINGW/Win32/

* install dll files in C:\Python27\DLLs
* install lib files in C:\Python27\libs

Installing
-----------

Using pip
~~~~~~~~~

Slycot supports the pip packaging system. You must first have pip installed.

On Debian Linux based systems you can install pip alongside your system Python
with the command::

        sudo apt-get install pip

Pip can then be used to install Slycot with the command::

        pip install slycot

Note that installing with pip may or may not require having the build
dependencies installed. There are some binary "wheels" available on PyPI, so if
those versions match with your system, you may be able to avoid installing from
source.

Using conda
~~~~~~~~~~~

Slycot can be installed for Linux or Mac via the conda package manager from the
Conda Forge channel with the following command::

  conda install -c conda-forge slycot

Note that there may be other versions for different operating systems available
on `other Anaconda channels <https://anaconda.org/search?q=slycot>`_, if the
Conda Forge version is not suitable for your needs.

From Source
~~~~~~~~~~~

Unpack the source code to a directory of your choice, e.g.
``/path/to/slycot_src/``, and execute::

   cd /path/to/slycot_src/
   python setup.py install

If the build fails and you are on a 64bit OS you may want to try::

   cd /path/to/slycot_src/
   python setup.py config_fc --arch="-march=x86-64" build
   python setup.py install

For Windows, and using Python(x,y), specify that you are using the
mingw compiler. Create a file::

   C:\\Python27\\Lib\\distutils\\distutils.cfg

with as contents::

        [build]
        compiler=mingw32
