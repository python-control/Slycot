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
    import configparser
except ImportError:
    import ConfigParser as configparser

try:
    from skbuild import setup
    from skbuild.command.sdist import sdist
except ImportError:
    raise ImportError('scikit-build must be installed before running setup.py')

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

# defaults
ISRELEASED = True
# assume a version set by conda, next update with git,
# otherwise count on default
VERSION = 'Unknown'


class GitError(RuntimeError):
    """Exception for git errors occuring in in git_version"""
    pass


def git_version(srcdir=None):
    """Return the git version, revision and cycle

    Uses rev-parse to get the revision tag to get the version number from the
    latest tag and detects (approximate) revision cycles

    """
    def _minimal_ext_cmd(cmd, srcdir):
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
        proc = subprocess.Popen(
            cmd,
            cwd=srcdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env)
        out, err = proc.communicate()
        if proc.returncode:
            errmsg = err.decode('ascii', errors='ignore').strip()
            raise GitError("git err; return code %d, error message:\n  '%s'"
                           % (proc.returncode, errmsg))
        return out

    try:
        GIT_VERSION = VERSION
        GIT_REVISION = 'Unknown'
        GIT_CYCLE = 0
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'], srcdir)
        GIT_REVISION = out.strip().decode('ascii')
        out = _minimal_ext_cmd(['git', 'tag'], srcdir)
        GIT_VERSION = out.strip().decode('ascii').split('\n')[-1][1:]
        out = _minimal_ext_cmd(['git', 'describe', '--tags',
                                '--long', '--always'], srcdir)
        try:
            # don't get a good description with shallow clones
            GIT_CYCLE = out.strip().decode('ascii').split('-')[1]
        except IndexError:
            pass
    except OSError:
        pass

    return GIT_VERSION, GIT_REVISION, GIT_CYCLE

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# slycot __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__SLYCOT_SETUP__ = True


def rewrite_setup_cfg(version, gitrevision, release):
    toreplace = dict(locals())
    data = ''.join(open('setup.cfg.in', 'r').readlines()).split('@')
    for k, v in toreplace.items():
        idx = data.index(k)
        data[idx] = v
    cfg = open('setup.cfg', 'w')
    cfg.write(''.join(data))
    cfg.close()


def get_version_info(srcdir=None):
    global ISRELEASED
    GIT_CYCLE = 0

    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of slycot.version messes up
    # the build under Python 3.
    if os.environ.get('CONDA_BUILD', False):
        FULLVERSION = os.environ.get('PKG_VERSION', '???')
        GIT_REVISION = os.environ.get('GIT_DESCRIBE_HASH', '')
        ISRELEASED = True
        rewrite_setup_cfg(FULLVERSION, GIT_REVISION, 'yes')
    elif os.path.exists('.git'):
        FULLVERSION, GIT_REVISION, GIT_CYCLE = git_version(srcdir)
        ISRELEASED = (GIT_CYCLE == 0)
        rewrite_setup_cfg(FULLVERSION, GIT_REVISION,
                          (ISRELEASED and 'yes') or 'no')
    elif os.path.exists('setup.cfg'):
        # valid distribution
        setupcfg = configparser.ConfigParser(allow_no_value=True)
        setupcfg.read('setup.cfg')

        FULLVERSION = setupcfg.get(section='metadata', option='version')

        if FULLVERSION is None:
            FULLVERSION = "Unknown"

        GIT_REVISION = setupcfg.get(section='metadata', option='gitrevision')

        if GIT_REVISION is None:
            GIT_REVISION = ""

        return FULLVERSION, GIT_REVISION
    else:

        # try to find a version number from the dir name
        dname = os.getcwd().split(os.sep)[-1]

        m = re.search(r'[0-9.]+', dname)
        if m:
            FULLVERSION = m.group()
            GIT_REVISION = ''

        else:
            FULLVERSION = VERSION
            GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.' + str(GIT_CYCLE)

    return FULLVERSION, GIT_REVISION

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

    # Rewrite the version file everytime
    VERSION, gitrevision = get_version_info(src_path)

    metadata = dict(
        name='slycot',
        packages=['slycot', 'slycot.tests'],
        cmake_languages=('C', 'Fortran'),
        version=VERSION,
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
        cmake_args=['-DSLYCOT_VERSION:STRING=' + VERSION,
                    '-DGIT_REVISION:STRING=' + gitrevision,
                    '-DISRELEASE:STRING=' + str(ISRELEASED),
                    '-DFULL_VERSION=' + VERSION + '.git' + gitrevision[:7]],
        zip_safe=False,
        install_requires=['numpy'],
        python_requires=">=3.7"
    )

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
    return


if __name__ == '__main__':
    setup_package()
