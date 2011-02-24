#!/usr/bin/env python
#
#       setup.py
#       
#       Copyright 2010 Enrico Avventi <avventi@kth.se>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License version 2 as 
#       published by the Free Software Foundation.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


from glob import glob
import os

if __name__ == '__main__':
    from numpy.distutils.core import setup,Extension
    wrappers = Extension('_wrapper', 
              libraries=['lapack'],
              sources=['slycot/src/_wrapper.pyf'] + glob(os.path.join('slycot','src','*.f')))
    setup(name = 'Slycot',
    	  author = 'Enrico Avventi',
    	  license = 'GPLv2',
    	  version = '0.1.0',
    	  description = 'A wrapper for the Slicot library',
    	  packages = ['slycot'],
    	  py_modules  = ['slycot/examples'],
    	  ext_package = 'slycot',
    	  ext_modules = [wrappers])

