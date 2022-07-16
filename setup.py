#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Slycot: a wrapper for the SLICOT control and systems library

Slycot wraps the SLICOT library which is used for control and systems analysis.

"""

import builtins
import os
import sys
import subprocess
import re
import platform

try:
    from skbuild import setup
    from skbuild.command.sdist import sdist
except ImportError:
    raise ImportError('scikit-build must be installed before running setup.py')

try:
    from setuptools_scm import get_version
except ImportError:
    raise ImportError('setuptools_scm must be installed before running setup.py')

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
License :: OSI Approved :: GNU General Public License v2 (GPLv2)
Programming Language :: C
Programming Language :: Fortran
Programming Language :: Python
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

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

def setup_package():
    src_path = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, src_path)


    metadata = dict(
        name='slycot',
        packages=['slycot', 'slycot.tests'],
        cmake_languages=('C', 'Fortran'),
        use_scm_version=True,
        maintainer="Slycot developers",
        maintainer_email="python-control-discuss@lists.sourceforge.net",
        description=DOCLINES[0],
        long_description=open('README.rst').read(),
        url='https://github.com/python-control/Slycot',
        author='Enrico Avventi et al.',
        license='GPL-2.0',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        cmdclass={"sdist": sdist_checked},
        zip_safe=False,
        install_requires=["numpy"],
        python_requires=">=3.7"
    )

    try:
        setup(**metadata)
    finally:
        del sys.path[0]


if __name__ == '__main__':
    setup_package()
