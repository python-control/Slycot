#!/usr/bin/env python
#
#       math.py
#
#       Copyright 2010 Enrico Avventi <avventi@Lonewolf>
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

from . import _wrapper

def mc01td(dico,dp,p):
    """ dp,stable,nz = mc01td(dico,dp,p)

    To determine whether or not a given polynomial P(x) with real
    coefficients is stable, either in the continuous-time or discrete-
    time case.

    A polynomial is said to be stable in the continuous-time case
    if all its zeros lie in the left half-plane, and stable in the
    discrete-time case if all its zeros lie inside the unit circle.


    Required arguments:
        dico : input string(len=1)
            Indicates whether the stability test to be applied to P(x) is in
            the continuous-time or discrete-time case as follows:
            = 'C':  continuous-time case;
            = 'D':  discrete-time case.
        dp : input int
            The degree of the polynomial P(x).  dp >= 0.
        p : input rank-1 array('d') with bounds (dp + 1)
            This array must contain the coefficients of P(x) in increasing
            powers of x.
    Return objects:
        dp : int
            If P(dp+1) = 0.0 on entry, then dp contains the index of the highest
            power of x for which P(dp+1) <> 0.0.
        stable : int
            Equal to 1 if P(x) if stable, 0 otherwise.
        nz : int
            The number of unstable zeros.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'dp', 'P', 'stable', 'nz', 'DWORK', 'IWARN'+hidden,
        'INFO'+hidden]
    out = _wrapper.mc01td(dico,dp,p)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        warings.warn('entry P(x) is the zero polynomial.')
    if out[-1] == 2:
        warings.warn('P(x) may have zeros very close to stability boundary.')
    if out[-2] > 0:
        warnings.warn('The degree of P(x) has been reduced to %i' %(dp-k))
    return out[:-2]


def mb05md(a, delta, balanc='N'):
    """Ar, Vr, Yr, VALRr, VALDr = mb05md(a, delta, balanc='N')

    Matrix exponential for a real non-defective matrix

    To compute exp(A*delta) where A is a real N-by-N non-defective
    matrix with real or complex eigenvalues and delta is a scalar
    value. The routine also returns the eigenvalues and eigenvectors
    of A as well as (if all eigenvalues are real) the matrix product
    exp(Lambda*delta) times the inverse of the eigenvector matrix of
    A, where Lambda is the diagonal matrix of eigenvalues.
    Optionally, the routine computes a balancing transformation to
    improve the conditioning of the eigenvalues and eigenvectors.

    Required arguments:
        A : input rank-2 array('d') with bounds (n,n)
            Square matrix
        delta : input 'd'
            The scalar value delta of the problem.
     
    Optional arguments:
        balanc : input char*1
            Indicates how the input matrix should be diagonally scaled
            to improve the conditioning of its eigenvalues as follows:
            = 'N':  Do not diagonally scale;
            = 'S':  Diagonally scale the matrix, i.e. replace A by
                    D*A*D**(-1), where D is a diagonal matrix chosen
                    to make the rows and columns of A more equal in
                    norm. Do not permute.

    Return objects:
        Ar : output rank-2 array('d') with bounds (n,n)
            Contains the solution matrix exp(A*delta)
        Vr : output rank-2 array('d') with bounds (n,n)
            Contains the eigenvector matrix for A.  If the k-th
            eigenvalue is real the k-th column of the eigenvector
            matrix holds the eigenvector corresponding to the k-th
            eigenvalue.  Otherwise, the k-th and (k+1)-th eigenvalues
            form a complex conjugate pair and the k-th and (k+1)-th
            columns of the eigenvector matrix hold the real and
            imaginary parts of the eigenvectors corresponding to these
            eigenvalues as follows.  If p and q denote the k-th and
            (k+1)-th columns of the eigenvector matrix, respectively,
            then the eigenvector corresponding to the complex
            eigenvalue with positive (negative) imaginary value is
            given by 
              p + q*j (p - q*j), where j^2  = -1.
        Yr : output rank-2 array('d') with bounds (n,n)
            contains an intermediate result for computing the matrix
            exponential.  Specifically, exp(A*delta) is obtained as the
            product V*Y, where V is the matrix stored in the leading
            N-by-N part of the array V. If all eigenvalues of A are
            real, then the leading N-by-N part of this array contains
            the matrix product exp(Lambda*delta) times the inverse of
            the (right) eigenvector matrix of A, where Lambda is the
            diagonal matrix of eigenvalues.

        VALr : output rank-1 array('c') with bounds (n)
            Contains the eigenvalues of the matrix A. The eigenvalues
            are unordered except that complex conjugate pairs of values
            appear consecutively with the eigenvalue having positive
            imaginary part first.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = [ 'balanc', 'n', 'delta', 'a', 'lda'+hidden, 'v', 'ldv'+hidden,
                 'y','ldy'+hidden,'valr','vali',
                 'iwork'+hidden,'dwork'+hidden,'ldwork'+hidden,'info'+hidden]
    out = _wrapper.mb05md(balanc=balanc,n=min(a.shape),delta=delta,a=a)
    if out[-1] == 0:
        return out[:-1]
    elif out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
    elif out[-1] > 0 and out[-1] <= n:
        error_text = "Incomplete eigenvalue calculation, missing %i eigenvalues" % out[-1]
    elif out[-1] == n+1:
        error_text = "Eigenvector matrix singular"
    elif out[-1] == n+2:
        error_text = "A matrix defective"
    e = ValueError(error_text)
    e.info = out[-1]
    raise e

"""
from slycot import mb05nd
import numpy as np
a = np.mat('[-2. 0; 0.1 -3.]')
mb05nd(a.shape[0], a, 0.1) 
"""

def mb05nd(a, delta, tol=1e-7):
    """F, H = mb05nd(n, a, delta, tol=1e-7)

    To compute

     (a)    F(delta) =  exp(A*delta) and

     (b)    H(delta) =  Int[F(s) ds] from s = 0 to s = delta,

    where A is a real N-by-N matrix and delta is a scalar value.

    Required arguments:
        A : input rank-2 array('d') with bounds (n,n)
            Square matrix
        delta : input 'd'
            The scalar value delta of the problem.
        tol : input 'd'
            Tolerance. A good value is sqrt(eps)

    Return objects:
        F : exp(A*delta)
        H : Int[F(s) ds] from s = 0 to s = delta,
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = [ 'n', 'delta', 'a', 'lda'+hidden, 'ex', 'ldex'+hidden,
                 'exint', 'ldexin'+hidden, 'tol', 'iwork'+hidden,
                 'dwork'+hidden, 'ldwork'+hidden]
    out = _wrapper.mb05nd(n=min(a.shape), delta=delta, a=a, tol=tol)
    if out[-1] == 0:
        return out[:-1]
    elif out[-1] < 0:
        error_text = "The following argument had an illegal value: " \
                     +arg_list[-out[-1]-1]
    elif out[-1] == n+1:
        error_text = "Delta too large"
    e = ValueError(error_text)
    e.info = out[-1]
    raise e
    


# to be replaced by python wrappers
