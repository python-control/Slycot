#!/usr/bin/env python
#
#       glue.py
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

import slycot._raw_wrapper

def tb03ad(n,m,p,A,B,C,D,leri):
	""" Find a left/right polynomial matrix representation of a given state-space system."""
	if leri == 'L':
		out = slycot._raw_wrapper.tb03ad_l(n,m,p,A,B,C,D)
		return out
	if leri == 'R':
		out = slycot._raw_wrapper.tb03ad_r(n,m,p,A,B,C,D)
		return out
	raise ValueError('leri must be either L or R')

def tc04ad(m,p,index,pcoeff,qcoeff,leri):
	""" Find a state-space realization of a given left/right polynomial matrix representation."""
	n = sum(index)
	if leri == 'L':
		out = slycot._raw_wrapper.tc04ad_l(m,p,index,pcoeff,qcoeff,n)
		return out
	if leri == 'R':
		out = slycot._raw_wrapper.tc04ad_r(m,p,index,pcoeff,qcoeff,n)
		return out
	raise ValueError('leri must be either L or R')
	
def tc01od(m,p,indlin,pcoeff,qcoeff,leri):
	""" Find a left/right polynomial matrix representation for a given right/left polynomial matrix representation."""
	if leri == 'L':
		out = slycot._raw_wrapper.tc01od_l(m,p,indlin,pcoeff,qcoeff)
		return out
	if leri == 'R':
		out = slycot._raw_wrapper.tc01od_r(m,p,indlin,pcoeff,qcoeff)
		return out
	raise ValueError('leri must be either L or R')
