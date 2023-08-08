************
Introduction
************

Welcome to the Slycot (`slycot`) User's Manual. 
This manual contains information on using the `slycot` package, 
including documentation for all functions in the package and examples illustrating their use.

Installation
============

The `slycot` package can be installed using conda or pip.  The
package requires `NumPy <http://www.numpy.org>`_.

For users with the Anaconda distribution of Python, the following
command can be used::

  conda install -c conda-forge slycot

This installs `slycot` from conda-forge, including the
`openblas` package.  NumPy will also be installed if
they are not already present.

.. note::
   Mixing packages from conda-forge and the default conda channel
   can sometimes cause problems with dependencies, so it is usually best to
   instally NumPy, SciPy, and Matplotlib from conda-forge as well.

To install using pip::

  pip install slycot

.. note::
   If you install Slycot using pip you'll need a development
   environment (e.g., Python development files, C and Fortran compilers).
   Pip installation can be particularly complicated for Windows.

Users can check to insure that slycot is installed 
correctly by running the command::

  python -c "import slycot"

and verifying that no error message appears. More information on the 
Slycot package can be obtained from the `Slycot project page
<https://github.com/python-control/Slycot>`_.

Alternatively, to install from source, first `download the source
<https://github.com/python-control/Slycot>`_ and unpack it.
To install in your home directory, use::
    
    pip install .

Getting started
===============

There are two different ways to use the package. For the default interface
described in :ref:`function-ref`, simply import the slycot package as follows::

    >>> import slycot

.. warning::
  In general, the second way is not recommended.

  To access the f2py interface described in :ref:`inner_function-ref`, simply import the slycot package as follows::

    >>> import slycot._wrapper