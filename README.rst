Slycot
=============

.. image:: https://travis-ci.org/jgoppert/Slycot.svg
        :target: https://travis-ci.org/jgoppert/Slycot
.. image:: https://coveralls.io/repos/jgoppert/Slycot/badge.png
        :target: https://coveralls.io/r/jgoppert/Slycot

Python wrapper for selected SLICOT routines, notably including solvers for
Riccati, Lyapunov and Sylvester equations.


Prerequisite:
-------------

You will need Numpy, a fortran compiler such as gfortran and BLAS/LAPACK 
libraries for building Slycot.

On Debian derivates you can install all the above with a single command::

        sudo apt-get build-dep python-scipy

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

C:\Python27\Lib\distutils\distutils.cfg

with as contents::

        [build]
        compiler=mingw32

To-Do
------
 
- write unit tests, already added test script, and simple test
- add examples in the doc-strings
