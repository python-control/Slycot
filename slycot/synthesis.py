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
import numpy as _np

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
        ldwork := None input int
            The length of the cache array. The default value is max(3, 6*n),
            for optimum performance it should be larger.
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

def sb02od(n,m,A,B,Q,R,dico,p=None,L=None,fact='N',uplo='U',sort='S',tol=0.0,ldwork=None):
    """ rcond,x,alfar,alfai,beta,s,t = sb02od(n,m,A,B,Q,R,dico,[p,L,fact,uplo,sort,tol,ldwork])
    
    To solve for X either the continuous-time algebraic Riccati
    equation
                              -1
        Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1)

    or the discrete-time algebraic Riccati equation
                                     -1
        X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2)

    where A, B, Q, R, and L are n-by-n, n-by-m, n-by-n, m-by-m and
    N-by-M matrices, respectively, such that Q = C'C, R = D'D and
    L = C'D; X is an n-by-n symmetric matrix.
    The routine also returns the computed values of the closed-loop
    spectrum of the system, i.e., the stable eigenvalues lambda(1),
    ..., lambda(n) of the corresponding Hamiltonian or symplectic
    pencil, in the continuous-time case or discrete-time case,
    respectively.

    Optionally, Q and/or R may be given in a factored form, Q = C'C, 
    R = D'D, and L may be treated as a zero matrix.

    The routine uses the method of deflating subspaces, based on
    reordering the eigenvalues in a generalized Schur matrix pair.
    
    Required arguments:
        n : input int
            The actual state dimension, i.e. the order of the matrices A, Q, 
            and X, and the number of rows of the matrices B and L.  n > 0.
        m : input int
            The number of system inputs, the order of the matrix R, and the 
            number of columns of the matrix B.  m > 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state matrix 
            A of the system.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input matrix 
            B of the system.
        Q : input rank-2 array('d') with bounds (n,n) or (p,n)
            If fact = 'N' or 'D', the leading n-by-n upper/lower triangular part 
            (depending on uplo) of this array must contain the upper/lower 
            triangular part of the symmetric state weighting matrix Q. 
            If fact = 'C' or 'B', the leading p-by-n part of this array must 
            contain the output matrix C of the system.
        R : input rank-2 array('d') with bounds (m,m) or (p,m)
            If fact = 'N' or 'C', the leading m-by-m upper/lower triangular part 
            (depending on uplo) of this array must contain the upper/lower 
            triangular part of the symmetric input weighting matrix R.
            If FACT = 'D' or 'B', the leading P-by-M part of this array must 
            contain the direct transmission matrix D of the system.
        dico : input string(len=1)
            Specifies the type of Riccati equation to be solved as follows:
            = 'C':  Equation (1), continuous-time case;
            = 'D':  Equation (2), discrete-time case.
    Optional arguments:
        p : input int
            The number of system outputs. If fact = 'C' or 'D' or 'B',
            p is the number of rows of the matrices C and/or D. p > 0.
        L := None input rank-2 array('d') with bounds (n,m)
            If L is None it will considered as a zero matrix of the appropriate
            dimensions. Otherwise the leading n-by-m part of this array must 
            contain the cross weighting matrix L.
        fact := 'N' input string(len=1)
            Specifies whether or not the matrices Q and/or R are factored, 
            as follows:
            = 'N':  Not factored, Q and R are given;
            = 'C':  C is given, and Q = C'C;
            = 'D':  D is given, and R = D'D;
            = 'B':  Both factors C and D are given, Q = C'C, R = D'D.
        uplo := 'U' input string(len=1)
            If fact = 'N', specifies which triangle of Q and R is stored, 
            as follows:
            = 'U':  Upper triangle is stored;
            = 'L':  Lower triangle is stored.
        sort := 'S' input string(len=1)
            Specifies which eigenvalues should be obtained in the top of 
            the generalized Schur form, as follows:
            = 'S':  Stable   eigenvalues come first;
            = 'U':  Unstable eigenvalues come first.
        tol := 0 input float
            The tolerance to be used in rank determination of the original 
            matrix pencil, specifically of the triangular factor obtained during 
            the reduction process. If tol <= 0 a default value is used.
        ldwork := None input int
            The length of the cache array. The default value is 
            max(7*(2*n+1)+16,16*n,2*n+m,3*m), for optimum performance it should 
            be larger.
    Return objects:
        rcond : float
            An estimate of the reciprocal of the condition number (in 
            the 1-norm) of the N-th order system of algebraic equations 
            from which the solution matrix X is obtained.
        X : rank-2 array('d') with bounds (n,n)
            The leading N-by-N part of this array contains the solution matrix 
            of the problem.
        alfar : rank-1 array('d') with bounds (2 * n)
        alfai : rank-1 array('d') with bounds (2 * n)
        beta : rank-1 array('d') with bounds (2 * n)
            The generalized eigenvalues of the 2n-by-2n matrix pair, ordered as 
            specified by sort. For instance, if sort = 'S', the leading n 
            elements of these arrays contain the closed-loop spectrum of the 
            system matrix A - BF, where F is the optimal feedback matrix computed
            based on the solution matrix X. Specifically,
            lambda(k) = [alfar(k)+j*alfai(k)]/beta(k) for k = 1,2,...,n.
        S : rank-2 array('d') with bounds (2*n+m,2 * n)
            ???
        T : rank-2 array('d') with bounds (2*n+m+1,2 * n)
            The leading 2n-by-2n part of this array contains the ordered upper 
            triangular form T of the second matrix in the reduced matrix pencil 
            associated to the optimal problem. That is,

                    (T   T  )
                    ( 11  12)
                T = (       ),
                    (0   T  )
                    (     22)

            where T  , T   and T   are n-by-n matrices.
                   11   12      22
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'jobb'+hidden, 'fact', 'uplo', 'jobl', 'sort', 'n', 
        'm', 'p', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'Q', 'LDQ'+hidden, 
        'R', 'LDR'+hidden, 'L', 'LDL'+hidden, 'rcond', 'X', 'LDX'+hidden, 
        'alfar', 'alfai', 'beta', 'S', 'LDS'+hidden, 'T', 'LDT'+hidden, 
        'U', 'LDU'+hidden, 'tol', 'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 
        'BWORK'+hidden, 'INFO'+hidden]
    if ldwork is None:
        ldwork = max([7*(2*n+1)+16,16*n,2*n+m,3*m])
    jobl = 'N'
    if L is None:
            L = _np.zeros((n,m))
            jobl = 'Z'
    out = None
    if fact == 'N':
        out = _wrapper.sb02od_n(dico,n,m,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if fact == 'C':
        if p is None:
            p = shape(Q)[0]
        out = _wrapper.sb02od_n(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if fact == 'D':
        if p is None:
            p = shape(R)[0]
        out = _wrapper.sb02od_n(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if fact == 'B':
        if p is None:
            p = shape(Q)[0]
        out = _wrapper.sb02od_n(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        raise ValueError(error_text)
    if out[-1] == 1:
        raise ArithmeticError('the computed extended matrix pencil is singular, possibly due to rounding errors')
    if out[-1] == 2:
        raise ArithmeticError('the QZ (or QR) algorithm failed')  
    if out[-1] == 3:
        raise ArithmeticError('reordering of the (generalized) eigenvalues failed')
    if out[-1] == 4:
        raise ArithmeticError('stability condition failed due to roudoff errors')
    if out[-1] == 5:
        raise ArithmeticError('the computed dimension of the solution does not equal N') 
    if out[-1] == 6:
        raise ArithmeticError('a singular matrix was encountered during the computation')
    return out[:-1]

def sb03md(n,C,A,U,dico,job='X',fact='N',trana='N',ldwork=None):
    """  X,scale,sep,ferr,wr,wi = sb03md(dico,n,C,A,U,[job,fact,trana,ldwork])
    
    To solve for X either the real continuous-time Lyapunov equation

       op(A)'*X + X*op(A) = scale*C                             (1)

    or the real discrete-time Lyapunov equation

       op(A)'*X*op(A) - X = scale*C                             (2)

    and/or estimate an associated condition number, called separation,
    where op(A) = A or A' (A**T) and C is symmetric (C = C').
    (A' denotes the transpose of the matrix A.) A is n-by-n, the right
    hand side C and the solution X are n-by-n, and scale is an output
    scale factor, set less than or equal to 1 to avoid overflow in X.
    
    Required arguments:
        n : input int
            The order of the matrices A, X, and C.  n > 0.
        C : input rank-2 array('d') with bounds (n,n)
            If job = 'X' or 'B', the leading n-by-n part of this array must 
            contain the symmetric matrix C. If job = 'S', C is not referenced.
        A : input rank-2 array('d') with bounds (n,n)
            On entry, the leading n-by-n part of this array must contain the 
            matrix A. If fact = 'F', then A contains an upper quasi-triangular 
            matrix in Schur canonical form; the elements below the upper 
            Hessenberg part of the array A are not referenced.
            On exit, the leading n-by-n upper Hessenberg part of this array 
            contains the upper quasi-triangular matrix in Schur canonical form 
            from the Schur factorization of A. The contents of array A is not
            modified if fact = 'F'.
        U : input rank-2 array('d') with bounds (n,n)
            If fact = 'F', then U is an input argument and on entry the leading 
            n-by-n part of this array must contain the orthogonal matrix U of 
            the real Schur factorization of A.
            If fact = 'N', then U is an output argument and on exit, it contains 
            the orthogonal n-by-n matrix from the real Schur factorization of A.
        dico : input string(len=1)
            Specifies the equation from which X is to be determined as follows:
            = 'C':  Equation (1), continuous-time case;
            = 'D':  Equation (2), discrete-time case.
    Optional arguments:
        job := 'X' input string(len=1)
            Specifies the computation to be performed, as follows:
            = 'X':  Compute the solution only;
            = 'S':  Compute the separation only;
            = 'B':  Compute both the solution and the separation.
        fact := 'N' input string(len=1)
            Specifies whether or not the real Schur factorization of the matrix 
            A is supplied on entry, as follows:
            = 'F':  On entry, A and U contain the factors from the real Schur 
            factorization of the matrix A;
            = 'N':  The Schur factorization of A will be computed and the factors 
            will be stored in A and U.
        trana := 'N' input string(len=1)
            Specifies the form of op(A) to be used, as follows:
            = 'N':  op(A) = A    (No transpose);
            = 'T':  op(A) = A**T (Transpose);
            = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
        ldwork := None input int
            The length of the cache array. The default value is max(2*n*n,3*n),
            for optimum performance it should be larger.
    Return objects:
        X : rank-2 array('d') with bounds (n,n)
            If job = 'X' or 'B', the leading n-by-n part contains the symmetric 
            solution matrix.
        scale : float
            The scale factor, scale, set less than or equal to 1 to prevent 
            the solution from overflowing.
        sep : float
            If job = 'S' or 'B', sep contains the estimated separation of 
            the matrices op(A) and -op(A)', if dico = 'C' or of op(A) and op(A)', 
            if dico = 'D'.
        ferr : float
            If job = 'B', ferr contains an estimated forward error bound for 
            the solution X. If X_true is the true solution, ferr bounds the 
            relative error in the computed solution, measured in the Frobenius
            norm:  norm(X - X_true)/norm(X_true).
        wr : rank-1 array('d') with bounds (n)
        wi : rank-1 array('d') with bounds (n)
            If fact = 'N', wr and wi contain the real and imaginary parts, 
            respectively, of the eigenvalues of A.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'fact', 'trana', 'n', 'A', 'LDA'+hidden, 'U', 
        'LDU'+hidden, 'C', 'LDC'+hidden, 'scale', 'sep', 'ferr', 'wr', 'wi', 
        'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 'INFO'+hidden]
    if ldwork is None:
	    ldwork = max(2*n*n,3*n)
    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be either D or C')
    out = _wrapper.sb03md(dico,n,C,A,U,job=job,fact=fact,trana=trana,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        raise ValueError(error_text)
    if out[-1] == n+1:
        if dico == 'D':
            warnings.warn('The matrix A has eigenvalues that are almost reciprocal.')
        else:
            warnings.warn('The matrix A and -A have common or very close eigenvalues.')
    else:
        if out[-1] > 0:
            warn_text = """The QR algorithm failed to compute all the eigenvalues 
(see LAPACK Library routine DGEES); elements %i:%i of wr and wi 
contains eigenvalues which have converged, A contains the partially 
converged Shur form'""" %(out[-1],n) # not sure about the indenting here
            warnings.warn(warn_text)
    return out[:-1]

# to be replaced by python wrappers
