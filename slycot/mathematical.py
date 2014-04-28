#!/usr/bin/env python
#
#       mathematical.py
#
#       Copyright 2013 René van Paassen <rene.vanpaassen@gmail.com>
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

# current website with overview:
# http://www.icm.tu-bs.de/NICONET/

from . import _wrapper

def mb05md(n, A, delta, balanc='N')
    """Ar, Vr, Yr, VALRr, VALDr = mb05md(n, a, delta, balanc='N')

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
        n : input int
            The order of the matrix a.  n >= 0.
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
           The order of the matrix A.  N >= 0.

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
                                       2
              p + q*j (p - q*j), where j  = -1.
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
    out = _wrapper.mb05md(balanc=balanc,n=n,delta=delta,a=a,
    



def mb05nd(N, DELTA, A, LDA, EX, LDEX, EXINT, LDEXIN,
           TOL, IWORK, DWORK, LDWORK, INFO):
    """

    Matrix exponential and integral for a real matrix

    To compute

     (a)    F(delta) =  exp(A*delta) and

     (b)    H(delta) =  Int[F(s) ds] from s = 0 to s = delta,

     where A is a real N-by-N matrix and delta is a scalar value.

     ARGUMENTS

     Input/Output Parameters

     N       (input) INTEGER
             The order of the matrix A.  N >= 0.

     DELTA   (input) DOUBLE PRECISION
             The scalar value delta of the problem.

     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
             The leading N-by-N part of this array must contain the
             matrix A of the problem. (Array A need not be set if
             DELTA = 0.)

     LDA     INTEGER
             The leading dimension of array A.  LDA >= max(1,N).

     EX      (output) DOUBLE PRECISION array, dimension (LDEX,N)
             The leading N-by-N part of this array contains an
             approximation to F(delta).

     LDEX    INTEGER
             The leading dimension of array EX.  LDEX >= MAX(1,N).

     EXINT   (output) DOUBLE PRECISION array, dimension (LDEXIN,N)
             The leading N-by-N part of this array contains an
             approximation to H(delta).

     LDEXIN  INTEGER
             The leading dimension of array EXINT.  LDEXIN >= MAX(1,N).

     Tolerances

     TOL     DOUBLE PRECISION
             The tolerance to be used in determining the order of the
             Pade approximation to H(t), where t is a scale factor
             determined by the routine. A reasonable value for TOL may
             be SQRT(EPS), where EPS is the machine precision (see
             LAPACK Library routine DLAMCH).

     Workspace

     IWORK   INTEGER array, dimension (N)

     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
             On exit, if INFO = 0, DWORK(1) returns the optimal value
             of LDWORK.

     LDWORK  INTEGER
             The length of the array DWORK. LDWORK >= MAX(1,N*(N+1)).
             For optimum performance LDWORK should be larger (2*N*N).

     Error Indicator

     INFO    INTEGER
             = 0:  successful exit;
             < 0:  if INFO = -i, the i-th argument had an illegal
                   value;
             > 0:  if INFO = i, the (i,i) element of the denominator of
                   the Pade approximation is zero, so the denominator
                   is exactly singular;
             = N+1:  if DELTA = (delta * frobenius norm of matrix A) is
                   probably too large to permit meaningful computation.
                   That is, DELTA > SQRT(BIG), where BIG is a
                   representable number near the overflow threshold of
                   the machine (see LAPACK Library Routine DLAMCH)."""


