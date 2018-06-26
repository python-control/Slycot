Slycot
=============

.. image:: https://travis-ci.org/python-control/slycot.svg?branch=master
        :target: https://travis-ci.org/python-control/slycot
.. image:: https://coveralls.io/repos/python-control/slycot/badge.png
        :target: https://coveralls.io/r/python-control/slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov, and Sylvester equations.


Prerequisites:
--------------

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

On Windows, we suggest installing on top of the Python(x,y) distribution, and
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
repository, but there are packages available on conda-forge for Linux and
Mac.  You can install with the following command::

  conda install -c conda-forge slycot


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

Additional tips for how to install slycot from source can be found in the
.travis.yml (commands used for Travis CI) and conda-recipe/ (conda
pre-requisities).  The hardest part about installing from source is getting
a working version of FORTRAN and LAPACK installed on your system and working
properly with Python.  If you are using conda, you can also get working
(binary) copies of LAPACK from conda-forge using the command::

	conda install -c conda-forge lapack

Note that in some cases you may need to set the LIBRARY_PATH environment
variable to pick up dependencies such as -lpythonN.m (where N.m is the
version of python you are using).


To-Do
------
 
- write unit tests, already added test script, and simple test
- add examples in the doc-strings
