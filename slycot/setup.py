#!/usr/bin/env python
from __future__ import division, print_function

import glob
import os
import sys


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('slycot', parent_package, top_path)

    # can disable building extension to test packaging
    build_fortran = False
    # to test c files from f2c
    build_c = True

    fortran_sources = []
    c_sources = []

    if build_fortran:
        fortran_sources = glob.glob(
            os.path.join('slycot', 'src', '*.f'))
    elif build_c:
        c_sources = glob.glob(
            os.path.join('slycot', 'src', '*.c'))
    else:
        print('WARNING FORTRAN BUILD DISABLED')
        fortran_sources = []

    f2py_sources = ['src/_wrapper.pyf']

    if sys.platform == 'win32':
        liblist = ['liblapack', 'libblas']
    else:
        liblist = ['lapack']

    if build_c:
        config.add_extension(
            name='_wrapper',
            libraries=liblist,
            sources=c_sources + f2py_sources)
    else:
        config.add_extension(
            name='_wrapper',
            libraries=liblist,
            sources=fortran_sources + f2py_sources)

    config.make_config_py()  # installs __config__.py

    config.add_subpackage('tests')

    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')
