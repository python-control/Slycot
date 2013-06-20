#!/usr/bin/env python
#
#       examples.py
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
from numpy import array, ones
import slycot

def sb02md_example():
    A = array([ [0, 1],
                [0, 0]])
    Q = array([ [1, 0],
                [0, 2]])
    G = array([ [0, 0],
                [0, 1]])
    out = slycot.sb02md(2,A,G,Q,'C')
    print('--- Example for sb02md ---')      
    print('The solution X is')
    print(out[0])
    print('rcond =', out[1])
    
def sb03md_example():
    from numpy import zeros
    A = array([ [3, 1, 1],
                [1, 3, 0],
                [0, 0, 3]])
    C = array([ [25, 24, 15],
                [24, 32,  8],
                [15,  8, 40]])
    U = zeros((3,3))
    out = slycot.sb03md(3,C,A,U,'D')
    print('--- Example for sb03md ---')   
    print('The solution X is')
    print(out[0])
    print('scaling factor:', out[1])
    
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
    print('--- Example for ab08nd ---')
    print('The finite invariant zeros are')
    print(eigvals(out[8][0:nu,0:nu],out[9][0:nu,0:nu]))
    
def mc01td_example():
    p = array([2, 0, 1, -1, 1])
    out = slycot.mc01td('C',4,p)
    print('--- Example for mc01td ...')
    if out[1]:
        print('The polynomial is stable')
    else:
        print('The polynomial has', out[2], 'unstable zeros')
        
def sb02od_example():
    from numpy import zeros, shape, dot, ones
    A = array([ [0, 1],
                [0, 0]])
    B = array([ [0],
                [1]])
    C = array([ [1, 0],
                [0, 1],
                [0, 0]])
    Q = dot(C.T,C)
    R = ones((1,1))
    out = slycot.sb02od(2,1,A,B,Q,R,'C')
    print('--- Example for sb02od ...')
    print('The solution X is')
    print(out[0])
    print('rcond =', out[1])
    
def tb03ad_example():
    A = array([ [1, 2, 0],
                [4,-1, 0],
                [0, 0, 1]])
    B = array([ [1, 0],
                [0, 0],
                [1, 0]])
    C = array([ [0, 1,-1],
                [0, 0, 1]])
    D = array([ [0, 1],
                [0, 0]])
    n = 3
    m = 1
    p = 2
    out = slycot.tb03ad(n,m,p,A,B,C,D,'R')
    #out = slycot.tb03ad_l(n,m,p,A,B,C,D)
    print('--- Example for tb03ad ...')
    print('The right polynomial representation of' )
    print('    W(z) = C(zI-A)^-1B + D')
    print('is the following:' )
    print('index', out[4])
    k_max = max(out[4]) + 1
    for k in range(0,k_max):
        print('P_%d =' %(k))
        print(out[5][0:m,0:m,k])
    for k in range(0,k_max):
        print('Q_%d =' %(k))
        print(out[6][0:m,0:p,k])
        
def tc04ad_example():
    from numpy import shape,zeros
    A = array([ [1, 2, 0],
                [4,-1, 0],
                [0, 0, 1]])
    B = array([ [1, 0],
                [0, 0],
                [1, 0]])
    C = array([ [0, 1,-1],
                [0, 0, 1]])
    D = array([ [0, 1],
                [0, 0]])
    n = 3
    m = 1
    p = 2
    out = slycot.tb03ad(n,m,p,A,B,C,D,'R')
    qcoeff = zeros((max(m,p),max(m,p),shape(out[6])[2]))
    qcoeff[0:shape(out[6])[0],0:shape(out[6])[1],0:shape(out[6])[1]]
    out2 = slycot.tc04ad(m,p,out[4],out[5][0:m,0:m,:],qcoeff,'R')
    print('--- Example for tb04ad ...')
    print('The system has a state space realization (A,B,C,D) where')
    print('A =')
    print(out2[1])
    print('B =')
    print(out2[2])
    print('C =')
    print(out2[3])
    print('D =')
    print(out2[4])

def tb01pd_example():
    A = array([[-1, 0],[0,-1]])
    B = ones((2,1))
    C = array([[0,1]])
    out = slycot.tb01pd(2, 1, 1, A, B, C)
    print('--- Example for tb01pd ...')
    print('Minimal realization for A, B, C')
    print('reduced order', out[-2])
    print(out)

