======
Slycot
======

.. image:: https://img.shields.io/pypi/v/slycot.svg
   :target: https://pypi.org/project/slycot/

.. image:: https://anaconda.org/conda-forge/slycot/badges/version.svg
   :target: https://anaconda.org/conda-forge/slycot

.. image:: https://github.com/python-control/Slycot/workflows/Build%20and%20Test%20Slycot/badge.svg
   :target: https://github.com/python-control/Slycot/actions

.. image:: https://coveralls.io/repos/github/python-control/Slycot/badge.svg?branch=master
   :target: https://coveralls.io/github/python-control/Slycot?branch=master

Slycot is a Python wrapper for selected routines from `SLICOT`_, the Subroutine Library in Systems and Control Theory, and is primarily a dependency for the `python-control`_ package.

Wrapped routines include those that:

 - solve for Ricatti, Lyapunov, and Sylvester equations
 - find multivariable zeros of systems
 - find reduced-order realizations of MIMO state-space systems
 - find H2 and L2 norms of MIMO systems
 - synthesize eigenvalue-placement controllers
 - synthesize H-infinity and H2 controllers
 - convert MIMO state-space representations to transfer functions, and vice versa

If you find a bug or have a question, please open an issue at `issues`_.

.. _`python-control`: https://github.com/python-control/python-control/
.. _SLICOT: https://www.slicot.org/
.. _issues: https://github.com/python-control/Slycot/issues

Installing
----------

We recommend using pre-compiled wheels from PyPI, or packages from
conda-forge.

If a binary package isn't available for your system, you can compile
from source, but you will need some software-building expertise.

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

If you are new to installing packages, read the `tutorial`_.

Use pip to install::

    pip install --only-binary :all: slycot

If this doesn't work there isn't a binary wheel for your platform.

To install with test dependencies, run::

    pip install --only-binary :all: slycot[test]

And you can then run the tests with ``pytest --pyargs slycot``

.. _`tutorial`: https://packaging.python.org/en/latest/tutorials/installing-packages/

Installing from conda-forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you haven't used conda before, read the `conda-forge getting started guide`_.

First install Miniforge_ or another conda package manager.  To install Slycot in a new environment ``project-1`` run::

    conda create -n project-1 slycot

.. _`conda-forge getting started guide`: https://conda-forge.org/docs/user/introduction/
.. _`Miniforge`: https://conda-forge.org/download/

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

See `Building Slycot`_.

.. _`Building Slycot`: https://github.com/python-control/Slycot/tree/master/BUILD.rst

Dependencies
------------

Slycot supports Python versions 3.10 or later, and Numpy 2.0 or later.

To run the Slycot unit tests and examples, you'll also need SciPy and
pytest.

License
-------
Up until version 0.4, Slycot used a version of SLICOT that was released under
the GPLv2 license. This requires Slycot to be released under the same license. In
December 2020, SLICOT 5.7 was released under BSD-3-Clause. However, as the
existing Slycot wrappers have been submitted by many contributors, we cannot
move away from GPLv2 unless we get the permission to do so by all authors.
Thus, Slycot remains licensed under GPLv2 until further notice.
