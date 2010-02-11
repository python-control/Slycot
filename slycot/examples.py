#!/usr/bin/env python
#
#       examples.py
#       
#       Copyright 2010 Enrico Avventi <avventi@kth.se>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
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
from numpy import array
import slycot

def sb02md_example():
	A = array([ [0, 1],
				[0, 0]])
	Q = array([ [1, 0],
				[0, 2]])
	G = array([ [0, 0],
				[0, 1]])
	out = slycot.sb02md('C',2,A,G,Q)
	print '--- Example for sb02md ---'		
	print 'The solution X is'
	print out[0]
	print 'rcond =', out[1]
	
def sb03md_example():
	from numpy import zeros
	A = array([ [3, 1, 1],
				[1, 3, 0],
				[0, 0, 3]])
	C = array([ [25, 24, 15],
				[24, 32,  8],
				[15,  8, 40]])
	U = zeros((3,3))
	out = slycot.sb03md('D',3,C,A,U)
	print '--- Example for sb03md ---'		
	print 'The solution X is'
	print out[0]
	print 'scaling factor:', out[1]
	
def ab08nd_example():
	from numpy import zeros, size
	from scipy.linalg import eigvals
	A = array([ [1, 0, 0, 0, 0, 0],
				[0, 1, 0, 0, 0, 0],
				[0, 0, 3, 0, 0, 0],
				[0, 0, 0,-4, 0, 0],
				[0, 0, 0, 0,-1, 0],
				[0, 0, 0, 0, 0, 3]])
	B = array([ [0,-1],
				[-1,0],
				[1,-1],
				[0, 0],
				[0, 1],
				[-1,-1]])
	C = array([ [1, 0, 0, 1, 0, 0],
				[0, 1, 0, 1, 0, 1],
				[0, 0, 1, 0, 0, 1]])
	D = zeros((3,2))
	out = slycot.ab08nd(6,2,3,A,B,C,D)
	nu = out[0]
	print '--- Example for ab08nd ---'
	print 'The finite invariant zeros are'
	print eigvals(out[8][0:nu,0:nu],out[9][0:nu,0:nu])
