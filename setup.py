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
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

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

# defaults
ISRELEASED = True
# assume a version set by conda, next update with git,
# otherwise count on default
VERSION = 'Unknown'

class GitError(RuntimeError):
    """Exception for git errors occuring in in git_version"""
    pass

# Return the git version, revision and cycle
#
# Uses rev-parse to get the revision
#      tag to get the version number from the latest tag
#      and detects (approximate) revision cycles
def git_version(srcdir=None):
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
            errmsg = err.decode('ascii',errors='ignore').strip()
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
        out = _minimal_ext_cmd(['git', 'describe', '--tags', '--long','--always'], srcdir)
        try:
            # don't get a good description with shallow clones, e.g., on Travis
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
        setupcfg = configparser.ConfigParser()
        setupcfg.read('setup.cfg')
        FULLVERSION = setupcfg['metadata'].get('version', 'Unknown')
        GIT_REVISION = setupcfg['metadata'].get('gitrevision', '')
        return FULLVERSION, GIT_REVISION
    else:

        # try to find a version number from the dir name
        dname = os.getcwd().split(os.sep)[-1]
        import re

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

from skbuild.command.sdist import sdist

class sdist_checked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        # slycot had no submodules currently
        # check_submodules()
        sdist.run(self)


def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    #os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    VERSION, gitrevision = get_version_info(src_path)
    
    metadata = dict(
        name='slycot',
        cmake_languages=('C', 'Fortran'),
        version=VERSION,
        maintainer="Slycot developers",
        maintainer_email="python-control-discuss@lists.sourceforge.net",
        description=DOCLINES[0],
        long_description="\n".join(DOCLINES[2:]),
        url='https://github.com/python-control/Slycot',
        author='Enrico Avventi et al.',
        license='GPL-2.0',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        cmdclass={"sdist": sdist_checked},
        cmake_args=[ '-DSLYCOT_VERSION:STRING=' + VERSION,
                     '-DGIT_REVISION:STRING=' + gitrevision,
                     '-DISRELEASE:STRING=' + str(ISRELEASED),
                     '-DFULL_VERSION=' + VERSION + '.git' + gitrevision[:7] ],
        #cmake_source_dir=src_path,
        zip_safe=False,
    )

    # Windows builds use Flang.
    # Flang detection and configuration is not automatic yet; the CMAKE
    # settings below are to circumvent that; when scikit-build and cmake
    # tools have improved, most of this might be removed?
    import platform
    if platform.system() == 'Windows':
        
        pbase = r'/'.join(sys.executable.split(os.sep)[:-1])
        env2cmakearg = {
            'FC': ('-DCMAKE_Fortran_COMPILER=',
                   pbase + r'/Library/bin/flang.exe'),
            'F2PY': ('-DF2PY_EXECUTABLE=',
                     pbase + r'/Scripts/f2py.exe'),
            'NUMPY_INCLUDE': ('-DNumPy_INCLUDE_DIR=',
                              pbase + r'/Include')
        }
            
        metadata['cmake_args'].extend([ 
	    '-GNMake Makefiles'])

        for k, v in env2cmakearg.items():
            print(k, v, os.environ.get(k, ''))
            envval = os.environ.get(k, None)
            if envval:
                # get from environment
                metadata['cmake_args'].append(
                    v[0] + envval.replace('\\', '/'))
            else:
                # default
                metadata['cmake_args'].append(v[0] + v[1])

        metadata['cmake_args'].extend([ 
            '-DCMAKE_Fortran_SIMULATE_VERSION=5.0.0',
	    '-DCMAKE_Fortran_COMPILER_ID=Flang',
	    '-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON' ])
        print(metadata['cmake_args'])
    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        #os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()
