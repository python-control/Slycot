#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Slycot: a wrapper for the SLICOT control and systems library

Slycot wraps the SLICOT library which is used for control and systems analysis.

"""

import builtins
import os
import subprocess

try:
    from skbuild import setup
    from skbuild.command.sdist import sdist
except ImportError:
    raise ImportError('scikit-build must be installed before running setup.py')

try:
    from setuptools_scm import get_version
except ImportError:
    raise ImportError('setuptools_scm must be installed before running setup.py')

# This is a bit hackish: we are setting a global variable so that the main
# slycot __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__SLYCOT_SETUP__ = True


def check_submodules():
    """ verify that the submodules are checked out and clean
        use `git submodule update --init`; on failure
    """
    if not os.path.exists('.git'):
        return
    with open('.gitmodules') as f:
        for l in f:
            if 'path' in l:
                p = l.split('=')[-1].strip()
                if not os.path.exists(p):
                    raise ValueError('Submodule %s missing' % p)

    proc = subprocess.Popen(['git', 'submodule', 'status'],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    for line in status.splitlines():
        if line.startswith('-') or line.startswith('+'):
            raise ValueError('Submodule not clean: %s' % line)


class sdist_checked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        check_submodules()
        sdist.run(self)

# These need to stay in setup.py
# https://scikit-build.readthedocs.io/en/latest/usage.html#setuptools-options
setup(
    packages=['slycot', 'slycot.tests'],
    cmdclass={'sdist': sdist_checked},
    cmake_languages=('C', 'Fortran'),
    include_package_data = False,
)
