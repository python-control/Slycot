Slycot
=============

.. image:: https://travis-ci.org/python-control/Slycot.svg?branch=master
        :target: https://travis-ci.org/python-control/Slycot
.. image:: https://coveralls.io/repos/python-control/Slycot/badge.png
        :target: https://coveralls.io/r/python-control/Slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov and Sylvester equations.


Prerequisite:
-------------

Slycot depends on Numpy, and if you are installing a binary distribution, Numpy
is the only prerequisite.

If you are installing Slycot from source, you will need a fortran
compiler such as gfortran, and BLAS/LAPACK libraries.

On Debian derivates you can install all the above with a single command::

        sudo apt-get build-dep python-scipy

On Mac, you will first need to install the `developer tools
<https://developer.apple.com/xcode/>`_.  You can then install gfortran using
`homebrew <http://brew.sh>`_ with::

        brew install gcc

On Windows, I suggest installing on top of the Python(x,y) distribution, and
grabbing BLAS and LAPACK libraries from: 

http://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.4.1/Dynamic-MINGW/Win32/

* install dll files in C:\Python27\DLLs
* install lib files in C:\Python27\libs


Installing
-----------

Using pip
~~~~~~~~~

Slycot supports the pip packaging system. You must first have
pip installed.

On debian linux based systems you can install pip with the command::

        sudo apt-get install pip

Pip can then be used to install Slycot wih the command::

        sudo pip install slycot

There are some binary "wheels" available on PyPI, so if those versions match
with your system, you may be able to avoid installing from source.

Using conda
~~~~~~~~~~~

If you use `Anaconda or conda <http://continuum.io/downloads>`_ on Linux or Mac,
it should be straighforward to install Slycot, without needing any compilers or
other prerequisites.  Slycot is not included in the standard conda package
repository, but there are packages available on http://binstar.org for Linux and
Mac.  You can install with the following command::

  conda install -c http://conda.binstar.org/cwrowley slycot


From Source
~~~~~~~~~~~

Unpack to a directory of your choice, say /path/to/slycot_src/, and execute::

        cd /path/to/slycot_src/
        # python setup.py install

Where # is for commands that needs to be executed as root/administrator. 

If the build fails and you are on a 64bit OS you may want to try::

        cd /path/to/slycot_src/
        python setup.py config_fc --arch="-march=x86-64" build
        # python setup.py install

For Windows, and using Python(x,y), specify that you are using the
mingw compiler. Create a file

C:\\Python27\\Lib\\distutils\\distutils.cfg

with as contents::

        [build]
        compiler=mingw32

To-Do
------
 
- write unit tests, already added test script, and simple test
- add examples in the doc-strings
