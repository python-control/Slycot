Fortran sources
---------------

This directory contains the f2py wrappers and some helper functions to work
with the SLICOT Library routines. SLICOT-Reference is a git submodule
forked from [SLICOT-Reference](https://github.com/SLICOT/SLICOT-Reference)
plus some backported improvements.

If your local copy of the SLICOT-Reference directory is empty, get the correct
version from python-control/SLICOT-Reference (see the Slycot toplevel directory
README for instructions).

The codes follow the Fortran 77 language conventions.  SLICOT routines make
calls to the state-of-the-art packages LAPACK (Linear Algebra Package) and BLAS
(Basic Linear Algebra Subprograms).
