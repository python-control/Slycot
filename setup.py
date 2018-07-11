#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Slycot: a wrapper for the SLICOT control and systems library

Slycot wraps the SLICOT library which is used for control and systems analysis.

"""
from skbuild import setup

DOCLINES = __doc__.split("\n")

import os
import sys
import subprocess


if sys.version_info[:2] < (2, 6) or (3, 0) <= sys.version_info[0:2] < (3, 2):
    raise RuntimeError("Python version 2.6, 2.7 or >= 3.2 required.")

if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins

# Fix a bug in python v3.4 installation
if (sys.version_info[0:2] == (3,4)):
    import importlib.machinery

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

MAJOR = 0
MINOR = 3
MICRO = 4
POST = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
if POST != 0:
    VERSION += '-post{:d}'.format(POST)

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# slycot __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__SLYCOT_SETUP__ = True


def get_version_info():
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of slycot.version messes up
    # the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('slycot/version.py'):
        # must be a source distribution, use existing version file
        try:
            from slycot.version import git_revision as GIT_REVISION
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing "
                              "slycot/version.py and the build directory "
                              "before building.")
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='slycot/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SLYCOT SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('slycot')
    config.get_version('slycot/version.py')  # sets config.version
    return config


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

from distutils.command.sdist import sdist


class sdist_checked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        # slycot had no submodules currently
        # check_submodules()
        sdist.run(self)


def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    write_version_py()

    metadata = dict(
        name='slycot',
	version=VERSION,
        maintainer="Slycot developers",
        maintainer_email="python-control-discuss@lists.sourceforge.net",
        description=DOCLINES[0],
        long_description="\n".join(DOCLINES[2:]),
        url='https://github.com/python-control/Slycot',
        author='Enrico Avventi et al.',
        license='GPLv2',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        cmdclass={"sdist": sdist_checked},
	cmake_args=[ '-DSLYCOT_VERSION=' + VERSION ],
	zip_safe=False,
    )

    # Windows builds use Flang.
    # Flang detection and configuration is not automatic yet; the CMAKE
    # settings below are to circumvent that; when scikit-build and cmake
    # tools have improved, most of this might be removed?
    import platform
    if platform.system() == 'Windows':
        pbase = r'/'.join(sys.executable.split(os.sep)[:-1])
        metadata['cmake_args'].extend([ 
	    '-DF2PY_EXECUTABLE=' + pbase + r'/Scripts/f2py.bat',
	    '-DCMAKE_Fortran_COMPILER=' + pbase + r'/Library/bin/flang.exe',
	    '-DCMAKE_Fortran_COMPILER_ID=Flang',
	    '-DCMAKE_C_COMPILER_ID=MSVC',
	    '-DCMAKE_C_COMPILER_VERSION=19.0.0', 
	    '-DNumPy_INCLUDE_DIR=' + pbase + r'/Include',
	    '-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON' ])
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()
