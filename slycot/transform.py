#!/usr/bin/env python
#
#       transform.py
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

from slycot import _wrapper

def tb03ad(n,m,p,A,B,C,D,leri,equil='N',tol=0.0,ldwork=None):
	""" A_min,b_min,C_min,nr,index,pcoeff,qcoeff,vcoeff,info = tb03ad_l(n,m,p,A,B,C,D,leri,[equil,tol,ldwork])
	
    To find a relatively prime left or right polynomial matrix representation
    with the same transfer matrix as that of a given state-space representation, 
    i.e. if leri = 'L'

        inv(P(s))*Q(s) = C*inv(s*I-A)*B + D
        
    or, if leri = 'R'
    
        Q(s)*inv(P(s)) = C*inv(s*I-A)*B + D.
        
    Additionally a minimal realization (A_min,B_min,C_min) of the original 
    system (A,B,C) is returned.

    Required arguments:
        n : input int
            The order of the state-space representation, i.e. the order of 
            the original state dynamics matrix A.  n > 0.
        m : input int
            The number of system inputs.  m > 0.
        p : input int
            The number of system outputs.  p > 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the original 
            state dynamics matrix A.
        B : input rank-2 array('d') with bounds (n,max(m,p))
            The leading n-by-m part of this array must contain the original 
            input/state matrix B; the remainder of the leading n-by-max(m,p) 
            part is used as internal workspace.
        C : input rank-2 array('d') with bounds (max(m,p),n)
            The leading p-by-n part of this array must contain the original 
            state/output matrix C; the remainder of the leading max(m,p)-by-n 
            part is used as internal workspace.
        D : input rank-2 array('d') with bounds (max(m,p),max(m,p))
            The leading p-by-m part of this array must contain the original 
            direct transmission matrix D; the remainder of the leading 
            max(m,p)-by-max(m,p) part is used as internal workspace.
        leri : input string(len=1)
            Indicates whether the left polynomial matrix representation or 
            the right polynomial matrix representation is required.
    Optional arguments:
        equil := 'N' input string(len=1)
            Specifies whether the user wishes to balance the triplet (A,B,C), 
            before computing a minimal state-space representation, as follows:
            = 'S':  Perform balancing (scaling);
            = 'N':  Do not perform balancing.
        tol := 0.0 input float
            The tolerance to be used in rank determination when transforming 
            (A, B). If tol <= 0 a default value is used.
        ldwork := max(2*n+3*max(m,p), p*(p+2)) input int
            The length of the cache array.
            ldwork >= max( n + max(n, 3*m, 3*p), pm*(pm + 2))
            where pm = p, if leri = 'L';
                  pm = m, if leri = 'R'.
            For optimum performance it should be larger.
    Return objects:
        A_min : rank-2 array('d') with bounds (n,n)
            The leading nr-by-nr part of this array contains the upper block 
            Hessenberg state dynamics matrix A_min of a minimal realization for 
            the original system.
        B_min : rank-2 array('d') with bounds (n,max(m,p))
            The leading nr-by-m part of this array contains the transformed 
            input/state matrix B_min.
        C_min : rank-2 array('d') with bounds (max(m,p),n)
            The leading p-by-nr part of this array contains the transformed 
            state/output matrix C_min.
        nr : int
            The order of the minimal state-space representation 
            (A_min,B_min,C_min).
        index : rank-1 array('i') with bounds either (p) or (m)
            If leri = 'L', index(i), i = 1,2,...,p, contains the maximum degree 
            of the polynomials in the i-th row of the denominator matrix P(s) 
            of the left polynomial matrix representation. These elements are 
            ordered so that index(1) >= index(2) >= ... >= index(p).
            If leri = 'R', index(i), i = 1,2,...,m, contains the maximum degree 
            of the polynomials in the i-th column of the denominator matrix P(s) 
            of the right polynomial matrix representation. These elements are 
            ordered so that index(1) >= index(2) >= ... >= index(m). 
        pcoeff : rank-3 array('d') with bounds either (p,p,n+1) or (m,m,n+1)
            If leri = 'L' then porm = p, otherwise porm = m.
            The leading porm-by-porm-by-kpcoef part of this array contains 
            the coefficients of the denominator matrix P(s), where 
            kpcoef = max(index) + 1.
            pcoeff(i,j,k) is the coefficient in s**(index(iorj)-k+1) of 
            polynomial (i,j) of P(s), where k = 1,2,...,kpcoef; if leri = 'L' 
            then iorj = I, otherwise iorj = J. Thus for leri = 'L', 
            P(s) = diag(s**index)*(pcoeff(.,.,1)+pcoeff(.,.,2)/s+...).
        qcoeff : rank-3 array('d') with bounds (p,m,n + 1) or (max(m,p),max(m,p))
            If leri = 'L' then porp = m, otherwise porp = p.
            If leri = 'L', the leading porm-by-porp-by-kpcoef part of this array 
            contains the coefficients of the numerator matrix Q(s).
            If leri = 'R', the leading porp-by-porm-by-kpcoef part of this array 
            contains the coefficients of the numerator matrix Q(s).
            qcoeff(i,j,k) is defined as for pcoeff(i,j,k).
        vcoeff : rank-3 array('d') with bounds (p,n,n+1) or (m,n,n+1)
            The leading porm-by-nr-by-kpcoef part of this array contains 
            the coefficients of the intermediate matrix V(s). 
            vcoeff(i,j,k) is defined as for pcoeff(i,j,k).
        info : int
            = 0: successful exit;
            < 0: if info = -i, the i-th argument had an illegal value;
            = 1: if a singular matrix was encountered during the computation of V(s);
            = 2: if a singular matrix was encountered during the computation of P(s).
	"""
	if leri == 'L':
	    if ldwork is None:
		    out = _wrapper.tb03ad_l(n,m,p,A,B,C,D,equil=equil,tol=tol)
        else:
            out = _wrapper.tb03ad_l(n,m,p,A,B,C,D,equil=equil,tol=tol,ldwork=ldwork)
        return out
	if leri == 'R':
		if ldwork is None:
		    out = _wrapper.tb03ad_r(n,m,p,A,B,C,D,equil=equil,tol=tol)
        else:
            out = _wrapper.tb03ad_r(n,m,p,A,B,C,D,equil=equil,tol=tol,ldwork=ldwork)
        return out
	raise ValueError('leri must be either L or R')

def tc04ad(m,p,index,pcoeff,qcoeff,leri,ldwork=None):
    """ n,rcond,a,b,c,d,info = tc04ad_l(m,p,index,pcoeff,qcoeff,leri,[ldwork])

    To find a state-space representation (A,B,C,D) with the same
    transfer matrix as that of a given left or right polynomial
    matrix representation, i.e.

        C*inv(sI-A)*B + D = inv(P(s))*Q(s) 

    or

        C*inv(sI-A)*B + D = Q(s)*inv(P(s))

    respectively.
        
    Required arguments:
        m : input int
            The number of system inputs.  m > 0.
        p := len(index) input int
            The number of system outputs.  p > 0.
        index : input rank-1 array('i') with bounds (p) or (m)
            If leri = 'L', index(i), i = 1,2,...,p, must contain the maximum 
            degree of the polynomials in the I-th row of the denominator matrix 
            P(s) of the given left polynomial matrix representation.
            If leri = 'R', index(i), i = 1,2,...,m, must contain the maximum 
            degree of the polynomials in the I-th column of the denominator 
            matrix P(s) of the given right polynomial matrix representation.
        pcoeff : input rank-3 array('d') with bounds (p,p,*) or (m,m,*)
            If leri = 'L' then porm = p, otherwise porm = m. The leading 
            porm-by-porm-by-kpcoef part of this array must contain 
            the coefficients of the denominator matrix P(s). pcoeff(i,j,k) is 
            the coefficient in s**(index(iorj)-K+1) of polynomial (I,J) of P(s), 
            where k = 1,2,...,kpcoef and kpcoef = max(index) + 1; if leri = 'L' 
            then iorj = i, otherwise iorj = j. Thus for leri = 'L', 
            P(s) = diag(s**index)*(pcoeff(.,.,1)+pcoeff(.,.,2)/s+...).
            If leri = 'R', pcoeff is modified by the routine but restored on exit.
        qcoeff : input rank-3 array('d') with bounds (p,m,*) or (max(m,p),max(m,p),*)
            If leri = 'L' then porp = m, otherwise porp = p. The leading 
            porm-by-porp-by-kpcoef part of this array must contain 
            the coefficients of the numerator matrix Q(s). 
            qcoeff(i,j,k) is defined as for pcoeff(i,j,k).
            If leri = 'R', qcoeff is modified by the routine but restored on exit.
        leri : input string(len=1)
            Indicates whether a left polynomial matrix representation or a right 
            polynomial matrix representation is input as follows:
            = 'L':  A left matrix fraction is input;
            = 'R':  A right matrix fraction is input.
    Optional arguments:
        ldwork := max(m,p)*(max(m,p)+4) input int
            The length of the cache array. ldwork >= max(m,p)*(max(m,p)+4)
            For optimum performance it should be larger.
    Return objects:
        n : int
            The order of the resulting state-space representation.
            That is, n = sum(index).
        rcond : float
            The estimated reciprocal of the condition number of the leading row 
            (if leri = 'L') or the leading column (if leri = 'R') coefficient 
            matrix of P(s).
            If rcond is nearly zero, P(s) is nearly row or column non-proper.
        A : rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array contains the state dynamics matrix A.
        B : rank-2 array('d') with bounds (n,max(m,p))
            The leading n-by-n part of this array contains the input/state matrix B; 
            the remainder of the leading n-by-max(m,p) part is used as internal 
            workspace.
        C : rank-2 array('d') with bounds (max(m,p),n)
            The leading p-by-n part of this array contains the state/output matrix C; 
            the remainder of the leading max(m,p)-by-n part is used as internal 
            workspace.
        D : rank-2 array('d') with bounds (max(m,p),max(m,p))
            The leading p-by-m part of this array contains the direct transmission 
            matrix D; the remainder of the leading max(m,p)-by-max(m,p) part is 
            used as internal workspace.
        info : int
            = 0:  successful exit;
            < 0:  if info = -i, the i-th argument had an illegal value;
            = 1:  if P(s) is not row (if leri = 'L') or column  (if leri = 'R') 
                proper. Consequently, no state-space representation is calculated.
    """
    if ldwork is None:
        ldwork = max(m,p)*(max(m,p)+4)
    n = sum(index)
    if leri == 'L':
        out = _wrapper.tc04ad_l(m,p,index,pcoeff,qcoeff,n)
        return out
    if leri == 'R':
        out = _wrapper.tc04ad_r(m,p,index,pcoeff,qcoeff,n)
        return out
	raise ValueError('leri must be either L or R')
	
def tc01od(m,p,indlin,pcoeff,qcoeff,leri):
	""" pcoeff,qcoeff,info = tc01od_l(m,p,indlim,pcoeff,qcoeff,leri)
	
	To find the dual right (left) polynomial matrix representation of a given 
	left (right) polynomial matrix representation, where the right and left 
	polynomial matrix representations are of the form Q(s)*inv(P(s)) and 
	inv(P(s))*Q(s) respectively.
	
	Required arguments:
	    m : input int
	        The number of system inputs.  m > 0.
	    p : input int
	        The number of system outputs.  p > 0.
	    indlim : input int
	        The highest value of k for which pcoeff(.,.,k) and qcoeff(.,.,k) 
	        are to be transposed.
            k = kpcoef + 1, where kpcoef is the maximum degree of the polynomials 
            in P(s).  indlim > 0.
	    pcoeff : input rank-3 array('d') with bounds (p,p,indlim) or (m,m,indlim)
            If leri = 'L' then porm = p, otherwise porm = m. 
            On entry, the leading porm-by-porm-by-indlim part of this array 
            must contain the coefficients of the denominator matrix P(s).
            pcoeff(i,j,k) is the coefficient in s**(indlim-k) of polynomial 
            (i,j) of P(s), where k = 1,2,...,indlim.
	    qcoeff : input rank-3 array('d') with bounds (max(m,p),max(m,p),indlim)
	        On entry, the leading p-by-m-by-indlim part of this array must 
	        contain the coefficients of the numerator matrix Q(s).
            qcoeff(i,j,k) is the coefficient in s**(indlim-k) of polynomial 
            (i,j) of Q(s), where k = 1,2,...,indlim.
	    leri : input string(len=1)
    Return objects:
        pcoeff : rank-3 array('d') with bounds (p,p,indlim)
            On exit, the leading porm-by-porm-by-indlim part of this array 
            contains the coefficients of the denominator matrix P'(s) of 
            the dual system.
        qcoeff : rank-3 array('d') with bounds (max(m,p),max(m,p),indlim)
            On exit, the leading m-by-p-by-indlim part of the array contains 
            the coefficients of the numerator matrix Q'(s) of the dual system.
        info : int
            = 0:  successful exit;
            < 0:  if info = -i, the i-th argument had an illegal value.
    """
	if leri == 'L':
		out = _wrapper.tc01od_l(m,p,indlin,pcoeff,qcoeff)
		return out
	if leri == 'R':
		out = _wrapper.tc01od_r(m,p,indlin,pcoeff,qcoeff)
		return out
	raise ValueError('leri must be either L or R')
	
# to be replaced by python wrappers
tb01id = _wrapper.tb01id
