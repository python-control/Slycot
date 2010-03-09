#!/usr/bin/env python
#
#       synthesis.py
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

def sb02md(n,A,G,Q,dico,hinv='D',uplo='U',scal='N',sort='S',ldwork=None):
    """  X,rcond,wr,wi,S,U = sb02md(dico,n,A,G,Q,[hinv,uplo,scal,sort,ldwork])
    
    To solve for X either the continuous-time algebraic Riccati
    equation
                               -1
         Q + A'*X + X*A - X*B*R  B'*X = 0                            (1)

    or the discrete-time algebraic Riccati equation
                                         -1
         X = A'*X*A - A'*X*B*(R + B'*X*B)  B'*X*A + Q                (2)

    where A, B, Q and R are n-by-n, n-by-m, n-by-n and m-by-m matrices
    respectively, with Q symmetric and R symmetric nonsingular; X is
    an N-by-N symmetric matrix.
                      -1
    The matrix G = B*R  B' must be provided on input, instead of B and
    R, that is, for instance, the continuous-time equation

         Q + A'*X + X*A - X*G*X = 0                                  (3)

    is solved, where G is an N-by-N symmetric matrix. Slycot Library
    routine sb02mt should be used to compute G, given B and R. sb02mt
    also enables to solve Riccati equations corresponding to optimal
    problems with coupling terms.

    The routine also returns the computed values of the closed-loop
    spectrum of the optimal system, i.e., the stable eigenvalues
    lambda(1),...,lambda(n) of the corresponding Hamiltonian or
    symplectic matrix associated to the optimal problem.
    
    Required arguments:
        n : input int
            The order of the matrices A, Q, G and X.  n > 0.
        A : input rank-2 array('d') with bounds (n,n)
            On entry, the leading n-by-n part of this array must contain 
            the coefficient matrix A of the equation. On exit, if dico = 'D',
                                                                       -1
            the leading N-by-N part of this array contains the matrix A  .
            Otherwise, the array A is unchanged on exit.
        G : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n upper triangular part (if uplo = 'U')
            or lower triangular part (if uplo = 'L') of this array
            must contain the upper triangular part or lower triangular
            part, respectively, of the symmetric matrix G.
        Q : input rank-2 array('d') with bounds (n,n)
            On entry, the leading n-by-n upper triangular part (if
            uplo = 'U') or lower triangular part (if uplo = 'L') of
            this array must contain the upper triangular part or lower
            triangular part, respectively, of the symmetric matrix Q.
        dico : input string(len=1)
            Specifies the type of Riccati equation to be solved as follows:
            = 'C':  Equation (3), continuous-time case;
            = 'D':  Equation (2), discrete-time case.
    Optional arguments:
        hinv := 'D' input string(len=1)
            If dico = 'D', specifies which symplectic matrix is to be 
            constructed, as follows:
            = 'D':  The matrix H in (5) (see SLICOT reference) is constructed;
            = 'I':  The inverse of the matrix H in (5) is constructed.
            hinv is not used if DICO = 'C'.
        uplo := 'U' input string(len=1)
            Specifies which triangle of the matrices G and Q is stored, 
            as follows:
            = 'U':  Upper triangle is stored;
            = 'L':  Lower triangle is stored.
        scal := 'N' input string(len=1)
            Specifies whether or not a scaling strategy should be used, 
            as follows:
            = 'G':  General scaling should be used;
            = 'N':  No scaling should be used.
        sort := 'S' input string(len=1)
            Specifies which eigenvalues should be obtained in the top of 
            the Schur form, as follows:
            = 'S':  Stable   eigenvalues come first;
            = 'U':  Unstable eigenvalues come first.
        ldwork := max(3,6*n) input int
            The length of the cache array. ldwork >= max(3, 6*n).
            For optimum performance it should be larger.
    Return objects:
        X : rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array contains the solution matrix 
            of the problem.
        rcond : float
            An estimate of the reciprocal of the condition number (in
            the 1-norm) of the N-th order system of algebraic
            equations from which the solution matrix X is obtained.
        wr : rank-1 array('d') with bounds (2 * n)
        wi : rank-1 array('d') with bounds (2 * n)
            These arrays contain the real and imaginary parts, respectively, 
            of the eigenvalues of the 2n-by-2n matrix S, ordered as specified 
            by sort (except for the case hinv = 'D', when the order is opposite 
            to that specified by sort). The leading n elements of these arrays 
            contain the closed-loop spectrum of the system
                          -1
            matrix A - B*R  *B'*X, if dico = 'C', or of the matrix
                              -1
            A - B*(R + B'*X*B)  B'*X*A, if dico = 'D'. Specifically,
            lambda(k) = wr(k) + j*wi(k), for k = 1,2,...,n.
        S : rank-2 array('d') with bounds (2 * n,2 * n)
            The leading 2n-by-2n part of this array contains the ordered real 
            Schur form S of the Hamiltonian or symplectic matrix H. That is,

                    (S   S  )
                    ( 11  12)
                S = (       ),
                    (0   S  )
                    (     22)

            where S  , S   and S   are n-by-n matrices.
                   11   12      22

        U : rank-2 array('d') with bounds (2 * n,2 * n)
            The leading 2n-by-2n part of this array contains the transformation 
            matrix U which reduces the Hamiltonian or symplectic matrix H to 
            the ordered real Schur form S. That is,

                    (U   U  )
                    ( 11  12)
                U = (       ),
                    (U   U  )
                    ( 21  22)

            where U  , U  , U   and U   are n-by-n matrices.
                   11   12   21      22"""
    
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'hinv', 'uplo', 'scal', 'sort', 'n', 'A', 'LDA'+hidden, 
    'G', 'LDG'+hidden, 'Q', 'LDQ'+hidden, 'rcond', 'wr', 'wi', 'S', 
    'LDS'+hidden, 'U', 'LDU'+hidden, 'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 
    'BWORK'+hidden, 'INFO'+hidden]
    if ldwork is None:
	    ldwork = max(3,6*n)
    out = _wrapper.sb02md(dico,n,A,G,Q,hinv=hinv,uplo=uplo,scal=scal,sort=sort,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        raise ValueError(error_text)
    if out[-1] == 1:
        raise ArithmeticError('matrix A is (numerically) singular in discrete-time case')
    if out[-1] == 2:
        raise ArithmeticError('the Hamiltonian or symplectic matrix H cannot be reduced to real Schur form')
    if out[-1] == 3:
        raise ArithmeticError('the real Schur form of the Hamiltonian or symplectic matrix H cannot be appropriately ordered')
    if out[-1] == 4:
        raise ArithmeticError('the Hamiltonian or symplectic matrix H has less than n stable eigenvalues')
    if out[-1] == 5:
        raise ArithmeticError('if the N-th order system of linear algebraic equations is singular to working precision')
    return out[:-1]

# to be replaced by python wrappers
sb02od = _wrapper.sb02od
sb03md = _wrapper.sb03md
