#!/usr/bin/env python
#
#       analysis.py
#       
#       Copyright 2010 Enrico Avventi <avventi@Lonewolf>
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

from slycot import _wrapper

def ab01nd(n,m,A,B,jobz='N',tol=0,ldwork=None):
	""" Ac,Bc,ncont,indcon,nblk,Z,tau,info = ab01nd_i(n,m,A,B,[jobz,tol,ldwork]) 
	
    To find a controllable realization for the linear time-invariant
    multi-input system

          dX/dt = A * X + B * U,

    where A and B are N-by-N and N-by-M matrices, respectively,
    which are reduced by this routine to orthogonal canonical form
    using (and optionally accumulating) orthogonal similarity
    transformations.  Specifically, the pair (A, B) is reduced to
    the pair (Ac, Bc),  Ac = Z' * A * Z,  Bc = Z' * B,  given by

          [ Acont     *    ]         [ Bcont ]
     Ac = [                ],   Bc = [       ],
          [   0    Auncont ]         [   0   ]

     and

             [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
             [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
             [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
     Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],
             [  .   .     . .    .     .  ]         [ .  ]
             [  .   .       .    .     .  ]         [ .  ]
             [  0   0   . . .  Ap,p-1 App ]         [ 0  ]

    where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
    p is the controllability index of the pair.  The size of the
    block  Auncont is equal to the dimension of the uncontrollable
    subspace of the pair (A, B).
	
    Required arguments:
        n : input int
            The order of the original state-space representation, i.e. 
            the order of the matrix A.  N > 0.
        m : input int
            The number of system inputs, or of columns of B.  M > 0.
        A : rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the original 
            state dynamics matrix A.
        B : rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input 
            matrix B.
    Optional arguments:        
        jobz := 'N' input string(len=1)
            Indicates whether the user wishes to accumulate in a matrix Z 
            the orthogonal similarity transformations for reducing the system, 
            as follows:
            = 'N':  Do not form Z and do not store the orthogonal transformations;
            = 'F':  Do not form Z, but store the orthogonal transformations in 
                    the factored form;
            = 'I':  Z is initialized to the unit matrix and the orthogonal 
                    transformation matrix Z is returned.
        tol := 0 input float
            The tolerance to be used in rank determination when transforming 
            (A, B). If tol <= 0 a default value is used.
        ldwork := max(n,3*m) input int
            The length of the cache array. ldwork >= max(n, 3*m).
            For optimum performance it should be larger.
    Return objects:
        Ac : rank-2 array('d') with bounds (n,n)
            The leading ncont-by-ncont part contains the upper block 
            Hessenberg state dynamics matrix Acont in Ac, given by Z'*A*Z, 
            of a controllable realization for the original system. The 
            elements below the first block-subdiagonal are set to zero.
        Bc : rank-2 array('d') with bounds (n,m)
            The leading ncont-by-m part of this array contains the transformed 
            input matrix Bcont in Bc, given by Z'*B, with all elements but the 
            first block set to zero.
        ncont : int
            The order of the controllable state-space representation.
        indcon : int
            The controllability index of the controllable part of the system 
            representation.
        nblk : rank-1 array('i') with bounds (n)
            The leading indcon elements of this array contain the the orders of 
            the diagonal blocks of Acont.
        Z : rank-2 array('d') with bounds (n,n)
            If jobz = 'I', then the leading N-by-N part of this array contains 
            the matrix of accumulated orthogonal similarity transformations 
            which reduces the given system to orthogonal canonical form.
            If jobz = 'F', the elements below the diagonal, with the array tau, 
            represent the orthogonal transformation matrix as a product of 
            elementary reflectors. The transformation matrix can then be 
            obtained by calling the LAPACK Library routine DORGQR.
            If jobz = 'N', the array Z is None.
        tau : rank-1 array('d') with bounds (n)
            The elements of tau contain the scalar factors of the
            elementary reflectors used in the reduction of B and A.
        info : int
            = 0:  successful exit;
            < 0:  if info = -i, the i-th argument had an illegal
                value.
	"""
	if ldwork is None:
	    ldwork = max(n,3*m)
	if jobz == 'N':
		out = _wrapper.ab01nd_n(n,m,A,B,tol=tol,ldwork=ldwork)
		# sets Z to None
		out[5] = None
		return out
	if jobz == 'I':
		return _wrapper.ab01nd_i(n,m,A,B,tol=tol,ldwork=ldwork)
	if jobz == 'F':
		return _wrapper.ab01nd_f(n,m,A,B,tol=tol,ldwork=ldwork)
	raise ValueError('jobz must be either N, I or F')

def ab05md(n1,m1,p1,n2,p2,A1,B1,C1,D1,A2,B2,C2,D2,uplo='U'):
    """ n,a,b,c,d,info = ab05md(n1,m1,p1,n2,p2,a1,b1,c1,d1,a2,b2,c2,d2,[uplo])
    
    To obtain the state-space model (A,B,C,D) for the cascaded
    inter-connection of two systems, each given in state-space form.
    
    Required arguments:
        n1 : input int
            The number of state variables in the first system, i.e. the order 
            of the matrix A1.  n1 > 0.
        m1 : input int
            The number of input variables for the first system. m1 > 0.
        p1 : input int
            The number of output variables from the first system and the number 
            of input variables for the second system. p1 > 0.
        n2 : input int
            The number of state variables in the second system, i.e. the order 
            of the matrix A2.  n2 > 0.
        p2 : input int
            The number of output variables from the second system. p2 > 0.
        A1 : input rank-2 array('d') with bounds (n1,n1)
            The leading n1-by-n1 part of this array must contain the state 
            transition matrix A1 for the first system.
        B1 : input rank-2 array('d') with bounds (n1,m1)
            The leading n1-by-m1 part of this array must contain the input/state 
            matrix B1 for the first system.
        C1 : input rank-2 array('d') with bounds (p1,n1)
            The leading p1-by-n1 part of this array must contain the state/output 
            matrix C1 for the first system.
        D1 : input rank-2 array('d') with bounds (p1,m1)
            The leading p1-by-m1 part of this array must contain the input/output 
            matrix D1 for the first system.
        A2 : input rank-2 array('d') with bounds (n2,n2)
            The leading n2-by-n2 part of this array must contain the state 
            transition matrix A2 for the second system.
        B2 : input rank-2 array('d') with bounds (n2,p1)
            The leading n2-by-p1 part of this array must contain the input/state 
            matrix B2 for the second system.
        C2 : input rank-2 array('d') with bounds (p2,n2)
            The leading p2-by-n2 part of this array must contain the state/output 
            matrix C2 for the second system.
        D2 : input rank-2 array('d') with bounds (p2,p1)
            The leading p2-by-p1 part of this array must contain the input/output 
            matrix D2 for the second system.
    Optional arguments:
        uplo := 'U' input string(len=1)
            Indicates whether the user wishes to obtain the matrix A in 
            the upper or lower block diagonal form, as follows:
                = 'U':  Obtain A in the upper block diagonal form;
                = 'L':  Obtain A in the lower block diagonal form.
    Return objects:
        n : int
            The number of state variables (n1 + n2) in the resulting system, 
            i.e. the order of the matrix A, the number of rows of B and 
            the number of columns of C.
        A : rank-2 array('d') with bounds (n1+n2,n1+n2)
            The leading N-by-N part of this array contains the state transition 
            matrix A for the cascaded system.
        B : rank-2 array('d') with bounds (n1+n2,m1)
            The leading n-by-m1 part of this array contains the input/state 
            matrix B for the cascaded system.
        C : rank-2 array('d') with bounds (p2,n1+n2)
            The leading p2-by-n part of this array contains the state/output 
            matrix C for the cascaded system.
        D : rank-2 array('d') with bounds (p2,m1)
            The leading p2-by-m1 part of this array contains the input/output 
            matrix D for the cascaded system.
        info : int
            = 0:  successful exit;
            < 0:  if INFO = -i, the i-th argument had an illegal value.
    """
    return _wrapper.ab05md(n1,m1,p1,n2,p2,A1,B1,C1,D1,A2,B2,C2,D2,uplo=uplo)
	
def ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,alpha=1.0,ldwork=None):
    """  n,A,B,C,D,info = ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,[alpha,ldwork])
    
    To obtain the state-space model (A,B,C,D) for the feedback inter-connection 
    of two systems, each given in state-space form.

    Required arguments:
        n1 : input int
            The number of state variables in the first system, i.e. the order 
            of the matrix A1.  n1 > 0.
        m1 : input int
            The number of input variables for the first system and the number 
            of output variables from the second system. m1 > 0.
        p1 : input int
            The number of output variables from the first system and the number 
            of input variables for the second system. p1 > 0.
        n2 : input int
            The number of state variables in the second system, i.e. the order 
            of the matrix A2.  n2 > 0.
        A1 : input rank-2 array('d') with bounds (n1,n1)
            The leading n1-by-n1 part of this array must contain the state 
            transition matrix A1 for the first system.
        B1 : input rank-2 array('d') with bounds (n1,m1)
            The leading n1-by-m1 part of this array must contain the input/state 
            matrix B1 for the first system.
        C1 : input rank-2 array('d') with bounds (p1,n1)
            The leading p1-by-n1 part of this array must contain the state/output 
            matrix C1 for the first system.
        D1 : input rank-2 array('d') with bounds (p1,m1)
            The leading p1-by-m1 part of this array must contain the input/output 
            matrix D1 for the first system.
        A2 : input rank-2 array('d') with bounds (n2,n2)
            The leading n2-by-n2 part of this array must contain the state 
            transition matrix A2 for the second system.
        B2 : input rank-2 array('d') with bounds (n2,p1)
            The leading n2-by-p1 part of this array must contain the input/state 
            matrix B2 for the second system.
        C2 : input rank-2 array('d') with bounds (m1,n2)
            The leading m1-by-n2 part of this array must contain the state/output 
            matrix C2 for the second system.
        D2 : input rank-2 array('d') with bounds (m1,p1)
            The leading m1-by-p1 part of this array must contain the input/output 
            matrix D2 for the second system.
    Optional arguments:
        alpha := 1.0 input float
            A coefficient multiplying the transfer-function matrix (or the 
            output equation) of the second system. i.e alpha = +1 corresponds 
            to positive feedback, and alpha = -1 corresponds to negative 
            feedback.
        ldwork := max(p1*p1,m1*m1,n1*p1) input int
            The length of the cache array. ldwork >= max(p1*p1,m1*m1,n1*p1).
    Return objects:
        n : int
            The number of state variables (n1 + n2) in the connected system, i.e. 
            the order of the matrix A, the number of rows of B and the number of 
            columns of C.
        A : rank-2 array('d') with bounds (n1+n2,n1+n2)
            The leading n-by-n part of this array contains the state transition 
            matrix A for the connected system.
        B : rank-2 array('d') with bounds (n1+n2,m1)
            The leading n-by-m1 part of this array contains the input/state 
            matrix B for the connected system.
        C : rank-3 array('d') with bounds (p1,n1,n2)
            The leading p1-by-n part of this array contains the state/output 
            matrix C for the connected system.
        D : rank-2 array('d') with bounds (p1,m1)
            The leading p1-by-m1 part of this array contains the input/output 
            matrix D for the connected system.
        info : int
            = 0:  successful exit;
            < 0:  if info = -i, the i-th argument had an illegal value.
            > 0:  if info = i, 1 <= i <= p1, the system is not completely 
                controllable. That is, the matrix   (I + alpha*D1*D2) is 
                exactly singular (the element U(i,i) of the upper triangular 
                factor of LU factorization is exactly zero), possibly due to
                rounding errors.
    """
    if ldwork is None:
	    ldwork = max(p1*p1,m1*m1,n1*p1)
    return _wrapper.ab05nd(n1,m1,p1,n2,alpha,A1,B1,C1,D1,A2,B2,C2,D2,ldwork=ldwork)
	
# to be replaced by python wrappers
ab07nd = _wrapper.ab07nd
ab08nd = _wrapper.ab08nd
