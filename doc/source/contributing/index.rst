.. this page is referenced from the front page but it's unnecessary as a navigation section for now.

:orphan:

Contributing to Slycot
======================

Development process and tools
-----------------------------

The development process is currently described on the `slycot github repo <https://github.com/python-control/Slycot>`_ and the `slycot github wiki <https://github.com/python-control/Slycot/wiki>`_.
You should be familiar with following topics:

- `git <https://git-scm.com/>`_
- `github <https://skills.github.com/>`_
- `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_
- `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
- `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`_
- `f2py <https://numpy.org/devdocs/f2py/index.html>`_

numpydoc
--------

Slycot uses numpydoc for the docstring style in order to provide support the Numpy docstring format in sphinx,
`see numpydoc example <https://numpydoc.readthedocs.io/en/latest/example.html>`_.

F2PY
----

Slycot heavily relias on `F2PY <https://numpy.org/devdocs/f2py/index.html>`_, which is currently a part of `NumPy <http://www.numpy.org>`_.