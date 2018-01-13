#!/usr/bin/env python
#
#       synthesis.py
#
#       Copyright 2010-2011 Enrico Avventi <avventi@Lonewolf>
#       Copyright 2011 Jerker Nordh <jerker.nordh@control.lth.se>
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
import numpy as _np
import warnings

def sb01bd(n,m,np,alpha,A,B,w,dico,tol=0.0,ldwork=None):
    """ A_z,w,nfp,nap,nup,F,Z = sb01bd(n,m,np,alpha,A,B,w,dico,[tol,ldwork])

    To determine the state feedback matrix F for a given system (A,B) such that
    the closed-loop state matrix A+B*F has specified eigenvalues.

    Required arguments
    ------------------

        n : int
            The dimension of the state vector, n >= 0.
        m : int
            The dimension of input vector, m >= 0.
        np : int
            The number of given eigenvalues. At most n eigenvalues can be
            assigned.  0 <= np <= n.
        alpha : float
            Specifies the maximum admissible value, either for real parts,
            if dico = 'C', or for moduli, if dico = 'D', of the eigenvalues of
            A which will not be modified by the eigenvalue assignment algorithm.
            alpha >= 0 if dico = 'D'.
        A : rank-2 array('d'), shape (n,n)
            State dynamics matrix.
        B : rank-2 array('d'), shape (n,m)
            Input/state matrix.
        w : rank-1 array('c'), shape (np,)
            Array of the desired eigenvalues of the closed-loop system state-matrix
            A+B*F. The eigenvalues can be unordered, except that complex conjugate
            pairs must appear consecutively.
        dico : {'C', 'D'}
            Specifies the type of the original system as follows:
            = 'C':  continuous-time system;
            = 'D':  discrete-time system.

    Optional arguments
    ------------------

        tol : float
            The absolute tolerance level below which the elements of A or B are
            considered zero (used for controllability tests).
            If tol <= 0 the default value is used.
        ldwork : int
            The length of the cache array. The default value is
            max(1,5*m,5*n,2*n+4*m), for optimum performance it should be larger.

    Returns
    -------

        A_z : rank-2 array('d'), shape (n,n)
            This array contains the matrix Z'*(A+B*F)*Z in a real Schur form.
            The diagonal block A[:nfp,:nfp] corresponds to the fixed (unmodified)
            eigenvalues having real parts less than alpha, if dico = 'C', or moduli
            less than alpha if dico = 'D'.
            The diagonal block A[n-nup:,n-nup:] corresponds to the uncontrollable
            eigenvalues detected by the eigenvalue assignment algorithm.
            The elements under the first subdiagonal are set to zero.
        w : rank-1 array('c'), shape (np,)
            The first part w[:nap] contain the assigned eigenvalues.
            The rest w[np-nap:] contain the unassigned eigenvalues.
        nfp : int
            The number of eigenvalues of A having real parts less than alpha,
            if dico = 'C', or moduli less than alpha, if dico = 'D'. These
            eigenvalues are not modified by the eigenvalue assignment algorithm.
        nap : int
            The number of assigned eigenvalues.
        nup : int
            The number of uncontrollable eigenvalues detected by the eigenvalue
            assignment algorithm.
        F : rank-2 array('d'), shape (m,n)
            The state feedback F, which assigns nap closed-loop eigenvalues and
            keeps unaltered n-nap open-loop eigenvalues.
        Z : rank-2 array('d'), shape (n,n)
            The orthogonal matrix Z which reduces the closed-loop system state
            matrix A + B*F to upper real Schur form.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             = 1:  the reduction of A to a real Schur form failed;
             = 2:  a failure was detected during the ordering of the
                real Schur form of A, or in the iterative process
                for reordering the eigenvalues of Z'*(A + B*F)*Z
                along the diagonal.
             = 3:  the number of eigenvalues to be assigned is less
                than the number of possibly assignable eigenvalues;
                nap eigenvalues have been properly assigned,
                but some assignable eigenvalues remain unmodified.
             = 4:  an attempt is made to place a complex conjugate
                pair on the location of a real eigenvalue. This
                situation can only appear when n-nfp is odd,
                np > n-nfp-nup is even, and for the last real
                eigenvalue to be modified there exists no available
                real eigenvalue to be assigned. However, nap
                eigenvalues have been already properly assigned.

    Example
    -------

    >>> import numpy as np
    >>> import slycot
    >>> A = np.array([[0, 1, 0],[0, 0, 1],[-2, 1, 3]])
    >>> B = np.array([[0], [0], [1]])
    >>> np.linalg.eig(A)[0] # open loop eigenvalues
    array([ 3.11490754,  0.74589831, -0.86080585])
    >>> w = np.array([0.5,0.4,0.2])
    >>> out = slycot.sb01bd(3,1,3,1,A,B,w,'D')
    >>> A_fb = A + np.dot(B,out[5])
    >>> np.linalg.eig(A_fb)[0] # closed loop eigenvalues
    array([ 0.2       ,  0.40000001,  0.5       ])
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'n', 'm', 'np', 'alpha', 'A', 'LDA'+hidden, 'B',
        'LDB'+hidden, 'wr'+hidden, 'wi'+hidden, 'nfp', 'nap', 'nup', 'F',
        'LDF'+hidden, 'Z', 'LDZ'+hidden, 'tol', 'DWORK'+hidden, 'ldwork',
        'IWARN'+hidden, 'INFO'+hidden]
    if ldwork is None:
        ldwork = max(1,5*m,5*n,2*n+4*m)
    A_z,wr,wi,nfp,nap,nup,F,Z,warn,info = _wrapper.sb01bd(dico,n,m,np,alpha,A,B,w.real,w.imag,tol=tol,ldwork=ldwork)
    if info < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-info-1]
        e = ValueError(error_text)
        e.info = info
        raise e
    if info == 1:
        e = ValueError('the reduction of A to a real Schur form failed')
        e.info = info
        raise e
    if info == 2:
        e = ValueError('a failure was detected during the ordering of eigenvalues')
        e.info = info
        raise e
    if info == 3:
        e = ValueError('the number of eigenvalues to be assigned is less than the number of possibly assignable eigenvalues')
        e.info = info
        raise e
    if info == 4:
        e = ValueError('an attempt was made to place a complex conjugate pair on the location of a real eigenvalue')
        e.info = info
        raise e
    if warn != 0:
        warnings.warn('%i violations of the numerical stability condition occured during the assignment of eigenvalues' % warn)
    # put togheter wr and wi into a complex array of eigenvalues
    w = _np.zeros(np,'complex64')
    w.real = wr[0:np]
    w.imag = wi[0:np]
    return A_z,w,nfp,nap,nup,F,Z

def sb02md(n,A,G,Q,dico,hinv='D',uplo='U',scal='N',sort='S',ldwork=None):
    """  X,rcond,w,S,U,A_inv = sb02md(dico,n,A,G,Q,[hinv,uplo,scal,sort,ldwork])

    To solve for X either the continuous-time algebraic Riccati
    equation
                               -1
         Q + A'*X + X*A - X*B*R  B'*X = 0                            (1)

    or the discrete-time algebraic Riccati equation
                                         -1
         X = A'*X*A - A'*X*B*(R + B'*X*B)  B'*X*A + Q                (2)

    where A, B, Q and R are n-by-n, n-by-m, n-by-n and m-by-m matrices
    respectively, with Q symmetric and R symmetric nonsingular; X is
    an n-by-n symmetric matrix.
                      -1
    The matrix G = B*R  B' must be provided on input, instead of B and
    R, that is, for instance, the continuous-time equation

         Q + A'*X + X*A - X*G*X = 0                                  (3)

    is solved, where G is an n-by-n symmetric matrix. Slycot Library
    routine sb02mt should be used to compute G, given B and R. sb02mt
    also enables to solve Riccati equations corresponding to optimal
    problems with coupling terms.

    The routine also returns the computed values of the closed-loop
    spectrum of the optimal system, i.e., the stable eigenvalues
    lambda(1),...,lambda(n) of the corresponding Hamiltonian or
    symplectic matrix associated to the optimal problem.

    Required arguments
    ------------------

        n : int
            The order of the matrices A, Q, G and X.  n > 0.
        A : rank-2 array('d'), shape (n,n)
        G : rank-2 array('d'), shape (n,n)
            The upper triangular part (if uplo = 'U') or lower triangular
            part (if uplo = 'L') of this array must contain the upper
            triangular part or lower triangular part, respectively, of the
            symmetric matrix G.
        Q : rank-2 array('d'), shape (n,n)
            The upper triangular part (if uplo = 'U') or lower triangular
            part (if uplo = 'L') of this array must contain the upper
            triangular part or lower triangular part, respectively,
            of the symmetric matrix Q.
        dico : {'C', 'D'}
            Specifies the type of Riccati equation to be solved as follows:
            = 'C':  Equation (3), continuous-time case;
            = 'D':  Equation (2), discrete-time case.

    Optional arguments
    ------------------

        hinv : {'D', 'I'}
            If dico = 'D', specifies which symplectic matrix is to be
            constructed, as follows:
            = 'D':  The Hamiltonian or sympletic matrix H is constructed;
            = 'I':  The inverse of the matrix H is constructed.
            The default value is 'D'. hinv is not used if DICO = 'C'.
        uplo : {'U', 'L'}
            Specifies which triangle of the matrices G and Q is stored,
            as follows:
            = 'U':  Upper triangle is stored;
            = 'L':  Lower triangle is stored.
            The default value is 'U'.
        scal : {'N', 'G'}
            Specifies whether or not a scaling strategy should be used,
            as follows:
            = 'G':  General scaling should be used;
            = 'N':  No scaling should be used.
            The default value is 'N'.
        sort : {'S', 'U'}
            Specifies which eigenvalues should be obtained in the top of
            the Schur form, as follows:
            = 'S':  Stable   eigenvalues come first;
            = 'U':  Unstable eigenvalues come first.
            The default value is 'S'.
        ldwork : int
            The length of the cache array. The default value is max(3, 6*n),
            for optimum performance it should be larger.

    Returns
    -------

        X : rank-2 array('d'), shape (n,n)
            The solution matrix of the problem.
        rcond : float
            An estimate of the reciprocal of the condition number (in
            the 1-norm) of the n-th order system of algebraic
            equations from which the solution matrix X is obtained.
        w : rank-1 array('c'), shape (2 * n)
            This array contain the eigenvalues of the 2n-by-2n matrix S, ordered
            as specified by sort (except for the case hinv = 'D', when the order
            is opposite to that specified by sort). The leading n elements of
            this array contain the closed-loop spectrum of the system
                          -1
            matrix A - B*R  *B'*X, if dico = 'C', or of the matrix
                              -1
            A - B*(R + B'*X*B)  B'*X*A, if dico = 'D'.
        S : rank-2 array('d'), shape (2 * n,2 * n)
            This array contains the ordered real Schur form S of the Hamiltonian
            or symplectic matrix H.
        U : rank-2 array('d'), shape (2 * n,2 * n)
            This array contains the transformation matrix U which reduces the
            Hamiltonian or symplectic matrix H to the ordered real Schur form S.
        A_inv : rank-2 array('d'), shape (n,n)
            The inverse of A if dico = 'D' or a copy of A itself otherwise.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             = 1:  if matrix A is (numerically) singular in discrete-
                time case;
             = 2:  if the Hamiltonian or symplectic matrix H cannot be
                reduced to real Schur form;
             = 3:  if the real Schur form of the Hamiltonian or
                symplectic matrix H cannot be appropriately ordered;
             = 4:  if the Hamiltonian or symplectic matrix H has less
                than n stable eigenvalues;
             = 5:  if the n-th order system of linear algebraic
                equations, from which the solution matrix X would
                be obtained, is singular to working precision.


        Example
        -------

            >>> import numpy as np
            >>> import slycot
            >>> A = np.array([[0, 1],[0, 0]])
            >>> Q = np.array([[1, 0],[0, 2]])
            >>> G = np.array([[0, 0],[0, 1]])
            >>> out = slycot.sb02md(2,A,G,Q,'C')
            >>> out[0] # X
            array([[ 2.,  1.],
                   [ 1.,  2.]])
            >>> out[1] # rcond
            0.3090169943749475
        """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'hinv', 'uplo', 'scal', 'sort', 'n', 'A', 'LDA'+hidden,
    'G', 'LDG'+hidden, 'Q', 'LDQ'+hidden, 'rcond', 'wr'+hidden, 'wi'+hidden, 'S',
    'LDS'+hidden, 'U', 'LDU'+hidden, 'IWORK'+hidden, 'DWORK'+hidden, 'ldwork',
    'BWORK'+hidden, 'INFO'+hidden]
    if ldwork is None:
        ldwork = max(3,6*n)
    A_inv,X,rcond,wr,wi,S,U,info = _wrapper.sb02md(dico,n,A,G,Q,hinv=hinv,uplo=uplo,scal=scal,sort=sort,ldwork=ldwork)
    if info < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-info-1]
        e = ValueError(error_text)
        e.info = info
        raise e
    if info == 1:
        e = ValueError('matrix A is (numerically) singular in discrete-time case')
        e.info = info
        raise e
    if info == 2:
        e = ValueError('the Hamiltonian or symplectic matrix H cannot be reduced to real Schur form')
        e.info = info
        raise e
    if info == 3:
        e = ValueError('the real Schur form of the Hamiltonian or symplectic matrix H cannot be appropriately ordered')
        e.info = info
        raise e
    if info == 4:
        e = ValueError('the Hamiltonian or symplectic matrix H has less than n stable eigenvalues')
        e.info = info
        raise e
    if info == 5:
        e = ValueError('if the N-th order system of linear algebraic equations is singular to working precision')
        e.info = info
        raise e
    w = _np.zeros(2*n,'complex64')
    w.real = wr[0:2*n]
    w.imag = wi[0:2*n]
    return X,rcond,w,S,U,A_inv

def sb02mt(n,m,B,R,A=None,Q=None,L=None,fact='N',jobl='Z',uplo='U',ldwork=None):
    """ A_b,B_b,Q_b,R_b,L_b,ipiv,oufact,G = sb02mt(n,m,B,R,[A,Q,L,fact,jobl,uplo,ldwork])

    To compute the following matrices

               -1
        G = B*R  *B',

                       -1
        A_b = A - B*R  *L',

                       -1
        Q_b = Q - L*R  *L',

    where A, B, Q, R, L, and G are n-by-n, n-by-m, n-by-n, m-by-m, n-by-m,
    and n-by-n matrices, respectively, with Q, R and G symmetric matrices.

    When R is well-conditioned with respect to inversion, standard algorithms
    for solving linear-quadratic optimization problems will then also solve
    optimization problems with coupling weighting matrix L.

    Required arguments
    ------------------

        n : int
            The order of the matrices A, Q, and G, and the number of rows of
            the matrices B and L.  n >= 0.
        m : int
            The order of the matrix R, and the number of columns of the matrices
            B and L.  m >= 0.
        B : rank-2 array('d'), shape (n,m)
        R : rank-2 array('d'), shape (m,m)
            If fact = 'N', the upper/lower triangular part of this array must
            contain the upper/lower triangular part, of the symmetric input
            weighting matrix R.
            If fact = 'C', the upper/lower triangular part of this array must
            contain the Cholesky factor of the positive definite input weighting
            matrix R.

    Optional arguments
    ------------------

        A : rank-2 array('d'), shape (n,n)
            If jobl = 'Z', this matrix is not needed.
        Q : rank-2 array('d'), shape (n,n)
            If jobl = 'Z', this matrix is not needed. Otherwise the upper/lower
            triangular part of this array (depending on uplo) must contain the
            corresponding part of matrix Q.
        L : rank-2 array('d'), shape (n,m)
            If jobl = 'Z', this matrix is not needed.
        fact : {'N', 'C'}
            Specifies how the matrix R is given (factored or not), as follows:
            = 'N':  Array R contains the matrix R,
            = 'C':  Array R contains the Cholesky factor of R.
            The default value is 'N'.
        jobl : {'Z', 'N'}
            When equal to 'Z', L is considered as a zero matrix and A,Q and L are
            not needed. A,Q and L are required otherwise.The default value is 'Z'.
        uplo : {'U', 'L'}
            Specifies which triangle of the matrices R and Q (if jobl = 'N')
            is stored, as follows:
            = 'U':  Upper triangle is stored;
            = 'L':  Lower triangle is stored.
            The default value is 'U'.
        ldwork : int
            The length of the cache array. Whenever fact = 'N' the default value
            is max(2,3*m,n*m), for optimum performance it should be larger.
            No cache is needed if fact = 'C', defaults at 1.

    Returns
    -------

        A_b : rank-2 array('d'), shape (n,n)
            If jobl = 'Z', this is None.
        B_b : rank-2 array('d'), shape (n,m)
                                                                  -1
            If oufact = 1 this array contains the matrix B*chol(R)  . It is a copy
            of input B if oufact = 2.
        Q_b : rank-2 array('d'), shape (n,n)
            If jobl = 'Z', this is None. Otherwise the upper/lower triangular part
            of this array contain the corresponding triangular part of matrix Q_b
            (depending on uplo).
        R_b : rank-2 array('d'), shape (m,m)
            If oufact = 1, the upper/lower triangular part of this array contains
            the Cholesky factor of the given input weighting matrix.
            If oufact = 2, the upper/lower triangular part of this array contains
            the factors of the UdU' or LdL' factorization, respectively, of the given
            input weighting matrix.
            If fact = 'C' it is a copy of input R.
        L_b : rank-2 array('d'), shape (n,m)
            If jobl = 'Z', this is None. If oufact = 1, this array contains the matrix
                     -1
            L*chol(R)  .If oufact = 2 this contains a copy of input L.
        ipiv : rank-1 array('i'), shape (m,)
            If oufact = 2, this array contains details of the interchanges
            performed and the block structure of the d factor in the UdU' or
            LdL' factorization of matrix R, as produced by LAPACK routine DSYTRF.
            Otherwise it is None.
        oufact : int
            Information about the factorization finally used.
            oufact = 1:  Cholesky factorization of R has been used;
            oufact = 2:  UdU' (if uplo = 'U') or LdL' (if uplo = 'L')
                         factorization of R has been used.
        G : rank-2 array('d'), shape (n,n)
            The upper/lower triangular part of this array contains the corresponding
            triangular part of the matrix G.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             = i:  if the i-th element (1 <= i <= m) of the d factor is
                exactly zero; the UdU' (or LdL') factorization has
                been completed, but the block diagonal matrix d is
                exactly singular;
             = m+1:  if the matrix R is numerically singular.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['JOBG'+hidden, 'jobl', 'fact', 'uplo', 'n', 'm', 'A',
    'LDA'+hidden, 'B', 'LDB'+hidden, 'Q', 'LDQ'+hidden, 'R', 'LDR'+hidden, 'L',
    'LDL'+hidden, 'ipiv', 'oufact', 'G', 'LDG'+hidden, 'IWORK'+hidden,
    'DWORK'+hidden, 'ldwork', 'INFO'+hidden]
    out = None
    if fact == 'N' and ldwork is None:
        ldwork = max(2,3*m,n*m)
    if jobl == 'Z':
        if fact == 'C':
            out = _wrapper.sb02mt_c(n,m,B,R,uplo=uplo)
        if fact == 'N':
            out = _wrapper.sb02mt_n(n,m,B,R,uplo=uplo,ldwork=ldwork)
        if out is None:
            raise ValueError('fact must be either C or N.')
    else:
        if A is None or Q is None or L is None:
            raise ValueError('matrices A,Q and L are required if jobl is not Z.')
        if fact == 'C':
            out = _wrapper.sb02mt_cl(n,m,A,B,Q,R,L,uplo=uplo)
        if fact == 'N':
            out = _wrapper.sb02mt_nl(n,m,A,B,Q,R,L,uplo=uplo,ldwork=ldwork)
        if out is None:
            raise ValueError('fact must be either C or N.')
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] > 0 and out[-1] <= m:
        e = ValueError('the %i-th elemend of d in the UdU (LdL) factorization is zero.'%out[-1])
        e.info = out[-1]
        raise e
    if out[-1] == m+1:
        e = ValueError('matrix R is numerically singular.')
        e.info = out[-1]
        raise e
    return out[:-1]

def sb02od(n,m,A,B,Q,R,dico,p=None,L=None,fact='N',uplo='U',sort='S',tol=0.0,ldwork=None):
    """ X,rcond,w,S,T = sb02od(n,m,A,B,Q,R,dico,[p,L,fact,uplo,sort,tol,ldwork])

    To solve for X either the continuous-time algebraic Riccati
    equation
                              -1
        Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1)

    or the discrete-time algebraic Riccati equation
                                     -1
        X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2)

    where A, B, Q, R, and L are n-by-n, n-by-m, n-by-n, m-by-m and
    n-by-m matrices, respectively, such that Q = C'C, R = D'D and
    L = C'D; X is an n-by-n symmetric matrix.
    The routine also returns the computed values of the closed-loop
    spectrum of the system, i.e., the stable eigenvalues w(1),
    ..., w(n) of the corresponding Hamiltonian or symplectic
    pencil, in the continuous-time case or discrete-time case,
    respectively.

    Optionally, Q and/or R may be given in a factored form, Q = C'C,
    R = D'D, and L may be treated as a zero matrix.

    The routine uses the method of deflating subspaces, based on
    reordering the eigenvalues in a generalized Schur matrix pair.

    Required arguments
    ------------------

        n : int
            The actual state dimension, i.e. the order of the matrices A, Q,
            and X, and the number of rows of the matrices B and L.  n > 0.
        m : int
            The number of system inputs, the order of the matrix R, and the
            number of columns of the matrix B.  m > 0.
        A : rank-2 array('d'), shape (n,n)
            The state matrix of the system.
        B : rank-2 array('d'), shape (n,m)
            The input matrix of the system.
        Q : rank-2 array('d'), shape (n,n) or (p,n)
            If fact = 'N' or 'D', the shape must be (n,n) and the upper/lower
            triangular part (depending on uplo) of this array must contain
            the corresponding triangular part of the symmetric state weighting
            matrix Q.
            If fact = 'C' or 'B', the shape must be (p,n) and of this array must
            contain the output matrix C of the system.
        R : rank-2 array('d'), shape (m,m) or (p,m)
            If fact = 'N' or 'C', the shape must be (m,m) and the upper/lower
            triangular part (depending on uplo) of this array must contain the
            corresponding triangular part of the symmetric input weighting matrix R.
            If FACT = 'D' or 'B', the shape must be (p,m) and this array must
            contain the direct transmission matrix D of the system.
        dico : {'C', 'D'}
            Specifies the type of Riccati equation to be solved as follows:
            = 'C':  Equation (1), continuous-time case;
            = 'D':  Equation (2), discrete-time case.

    Optional arguments
    ------------------

        p : int
            The number of system outputs. If fact = 'C' or 'D' or 'B',
            p is the number of rows of the matrices C and/or D. p > 0.
            If fact = 'N' it is not needed.
        L : rank-2 array('d'), shape (n,m)
            If L is not specified it will considered as a zero matrix of the
            appropriate dimensions.
        fact : {'N', 'C', 'D', 'B'}
            Specifies whether or not the matrices Q and/or R are factored,
            as follows:
            = 'N':  Not factored, Q and R are given;
            = 'C':  C is given, and Q = C'C;
            = 'D':  D is given, and R = D'D;
            = 'B':  Both factors C and D are given, Q = C'C, R = D'D.
            The default value is 'N'.
        uplo : {'U', 'L'}
            If fact = 'N', specifies which triangle of Q and R is stored,
            as follows:
            = 'U':  Upper triangle is stored;
            = 'L':  Lower triangle is stored.
            The default value is 'U'.
        sort : {'S', 'U'}
            Specifies which eigenvalues should be obtained in the top of
            the generalized Schur form, as follows:
            = 'S':  Stable   eigenvalues come first;
            = 'U':  Unstable eigenvalues come first.
            The default value is 'S'.
        tol : float
            The tolerance to be used in rank determination of the original
            matrix pencil, specifically of the triangular factor obtained during
            the reduction process. If tol <= 0 a default value is used.
        ldwork : int
            The length of the cache array. The default value is
            max(7*(2*n+1)+16,16*n,2*n+m,3*m), for optimum performance it should
            be larger.

    Returns
    -------

        X : rank-2 array('d'), shape (n,n)
            The solution matrix of the problem.
        rcond : float
            An estimate of the reciprocal of the condition number (in
            the 1-norm) of the n-th order system of algebraic equations
            from which the solution matrix X is obtained.
        w : rank-1 array('c'), shape (2 * n)
            The generalized eigenvalues of the 2n-by-2n matrix pair, ordered as
            specified by sort. For instance, if sort = 'S', the leading n
            elements of these arrays contain the closed-loop spectrum of the
            system matrix A - BF, where F is the optimal feedback matrix computed
            based on the solution matrix X.
        S : rank-2 array('d'), shape (2*n+m,2 * n)
            This array contains the ordered real Schur form S of the first matrix
            in the reduced matrix pencil associated to the optimal problem, or of
            the corresponding Hamiltonian matrix
        T : rank-2 array('d'), shape (2*n+m+1,2 * n)
            This array contains the ordered upper triangular form T of the second
            matrix in the reduced matrix pencil associated to the optimal problem.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             = 1:  if the computed extended matrix pencil is singular,
                possibly due to rounding errors;
             = 2:  if the QZ (or QR) algorithm failed;
             = 3:  if reordering of the (generalized) eigenvalues failed;
             = 4:  if after reordering, roundoff changed values of
                some complex eigenvalues so that leading eigenvalues
                in the (generalized) Schur form no longer satisfy
                the stability condition; this could also be caused
                due to scaling;
             = 5:  if the computed dimension of the solution does not
                equal n;
             = 6:  if a singular matrix was encountered during the
                computation of the solution matrix X.

    Example
    -------

    >>> import numpy as np
    >>> import slycot
    >>> A = np.array([[0, 1],[0, 0]])
    >>> B = np.array([[0],[1]])
    >>> C = np.array([[1, 0],[0, 1],[0, 0]])
    >>> Q = np.dot(C.T,C)
    >>> R = np.ones((1,1))
    >>> out = slycot.sb02od(2,1,A,B,Q,R,'C')
    >>> out[0] # X
    array([[ 1.73205081,  1.        ],
           [ 1.        ,  1.73205081]])
    >>> out[1] # rcond
    0.4713846564431681
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'jobb'+hidden, 'fact', 'uplo', 'jobl', 'sort', 'n',
        'm', 'p', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'Q', 'LDQ'+hidden,
        'R', 'LDR'+hidden, 'L', 'LDL'+hidden, 'rcond', 'X', 'LDX'+hidden,
        'alfar'+hidden, 'alfai'+hidden, 'beta'+hidden, 'S', 'LDS'+hidden, 'T',
        'LDT'+hidden, 'U', 'LDU'+hidden, 'tol', 'IWORK'+hidden, 'DWORK'+hidden,
        'ldwork', 'BWORK'+hidden, 'INFO'+hidden]
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
        out = _wrapper.sb02od_c(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if fact == 'D':
        if p is None:
            p = shape(R)[0]
        out = _wrapper.sb02od_d(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if fact == 'B':
        if p is None:
            p = shape(Q)[0]
        out = _wrapper.sb02od_b(dico,n,m,p,A,B,Q,R,L,uplo=uplo,jobl=jobl,sort=sort,tol=tol,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = info
        raise e
    if out[-1] == 1:
        e = ValueError('the computed extended matrix pencil is singular, possibly due to rounding errors')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ValueError('the QZ (or QR) algorithm failed')
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        e = ValueError('reordering of the (generalized) eigenvalues failed')
        e.info = out[-1]
        raise e
    if out[-1] == 4:
        e = ValueError('stability condition failed due to roudoff errors')
        e.info = out[-1]
        raise e
    if out[-1] == 5:
        e = ValueError('the computed dimension of the solution does not equal N')
        e.info = out[-1]
        raise e
    if out[-1] == 6:
        e = ValueError('a singular matrix was encountered during the computation')
        e.info = out[-1]
        raise e
    rcond,X,alphar,alphai,beta,S,T = out[:-1]
    w = _np.zeros(2*n,'complex64')
    w.real = alphar[0:2*n]/beta[0:2*n]
    w.imag = alphai[0:2*n]/beta[0:2*n]
    return X,rcond,w,S,T

def sb03md(n,C,A,U,dico,job='X',fact='N',trana='N',ldwork=None):
    """  X,scale,sep,ferr,w = sb03md(dico,n,C,A,U,[job,fact,trana,ldwork])

    To solve for X either the real continuous-time Lyapunov equation

       op(A)'*X + X*op(A) = scale*C                             (1)

    or the real discrete-time Lyapunov equation

       op(A)'*X*op(A) - X = scale*C                             (2)

    and/or estimate an associated condition number, called separation,
    where op(A) = A or A' (A**T) and C is symmetric (C = C').
    (A' denotes the transpose of the matrix A.) A is n-by-n, the right
    hand side C and the solution X are n-by-n, and scale is an output
    scale factor, set less than or equal to 1 to avoid overflow in X.

    Required arguments
    ------------------

        n : input int
            The order of the matrices A, X, and C.  n > 0.
        C : input rank-2 array('d'), shape  (n,n)
            If job = 'X' or 'B', the leading n-by-n part of this array must
            contain the symmetric matrix C. If job = 'S', C is not referenced.
        A : input rank-2 array('d'), shape  (n,n)
            On entry, the leading n-by-n part of this array must contain the
            matrix A. If fact = 'F', then A contains an upper quasi-triangular
            matrix in Schur canonical form; the elements below the upper
            Hessenberg part of the array A are not referenced.
            On exit, the leading n-by-n upper Hessenberg part of this array
            contains the upper quasi-triangular matrix in Schur canonical form
            from the Schur factorization of A. The contents of array A is not
            modified if fact = 'F'.
        U : input rank-2 array('d'), shape  (n,n)
            If fact = 'F', then U is an input argument and on entry the leading
            n-by-n part of this array must contain the orthogonal matrix U of
            the real Schur factorization of A.
            If fact = 'N', then U is an output argument and on exit, it contains
            the orthogonal n-by-n matrix from the real Schur factorization of A.
        dico : input string(len=1)
            Specifies the equation from which X is to be determined as follows:
            = 'C':  Equation (1), continuous-time case;
            = 'D':  Equation (2), discrete-time case.

    Optional arguments
    ------------------

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

    Return objects
    --------------

        X : rank-2 array('d'), shape  (n,n)
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
        w : rank-1 array('c'), shape  (n)
            If fact = 'N', this array contain the eigenvalues of A.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             > 0:  if info = i, the QR algorithm failed to compute all
                the eigenvalues (see LAPACK Library routine DGEES);
                elements i+1:n of w contain eigenvalues which have converged,
                and A contains the partially converged Schur form;
             = N+1:  if dico = 'C', and the matrices A and -A' have
                common or very close eigenvalues, or
                if dico = 'D', and matrix A has almost reciprocal
                eigenvalues (that is, lambda(i) = 1/lambda(j) for
                some i and j, where lambda(i) and lambda(j) are
                eigenvalues of A and i <> j); perturbed values were
                used to solve the equation (but the matrix A is
                unchanged).
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'fact', 'trana', 'n', 'A', 'LDA'+hidden, 'U',
        'LDU'+hidden, 'C', 'LDC'+hidden, 'scale', 'sep', 'ferr', 'wr'+hidden,
        'wi'+hidden, 'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 'INFO'+hidden]
    if ldwork is None:
        ldwork = max(2*n*n,3*n)
    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be either D or C')
    out = _wrapper.sb03md(dico,n,C,A,U,job=job,fact=fact,trana=trana,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == n+1:
        if dico == 'D':
            error_text = 'The matrix A has eigenvalues that are almost reciprocal.'
        else:
            error_text = 'The matrix A and -A have common or very close eigenvalues.'
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    else:
        if out[-1] > 0:
            error_text = """The QR algorithm failed to compute all the eigenvalues
(see LAPACK Library routine DGEES); elements %i:%i of w contains
eigenvalues which have converged, A contains the partially
converged Shur form'""" %(out[-1],n) # not sure about the indenting here
            e = ValueError(error_text)
            e.info = out[-1]
            raise e
    X,scale,sep,ferr,wr,wi = out[:-1]
    w = _np.zeros(n,'complex64')
    w.real = wr[0:n]
    w.imag = wi[0:n]
    return X,scale,sep,ferr,w

def sb03od(n,m,A,Q,B,dico,fact='N',trans='N',ldwork=None):
    """  U,scale,w = sb03od(dico,n,m,A,Q,B,[fact,trans,ldwork])

    To solve for X = op(U)'*op(U) either the stable non-negative
    definite continuous-time Lyapunov equation
                                 2
      op(A)'*X + X*op(A) = -scale *op(B)'*op(B)                   (1)

    or the convergent non-negative definite discrete-time Lyapunov
    equation
                                 2
      op(A)'*X*op(A) - X = -scale *op(B)'*op(B)                   (2)

    where op(K) = K or K' (i.e., the transpose of the matrix K), A is
    an N-by-N matrix, op(B) is an M-by-N matrix, U is an upper
    triangular matrix containing the Cholesky factor of the solution
    matrix X, X = op(U)'*op(U), and scale is an output scale factor,
    set less than or equal to 1 to avoid overflow in X. If matrix B
    has full rank then the solution matrix X will be positive-definite
    and hence the Cholesky factor U will be nonsingular, but if B is
    rank deficient then X may be only positive semi-definite and U
    will be singular.

    In the case of equation (1) the matrix A must be stable (that
    is, all the eigenvalues of A must have negative real parts),
    and for equation (2) the matrix A must be convergent (that is,
    all the eigenvalues of A must lie inside the unit circle).

    Required arguments
    ------------------

        n : input int
            The order of the matrix A and the number of columns in
            matrix op(B).  n >= 0.
        m : input int
            The number of rows in matrix op(B).  m >= 0.
        A : input rank-2 array('d'), shape  (n,n)
            On entry, the leading n-by-n part of this array must
            contain the matrix A. If fact = 'F', then A contains
            an upper quasi-triangular matrix S in Schur canonical
            form; the elements below the upper Hessenberg part of the
            array A are not referenced.
            On exit, the leading n-by-n upper Hessenberg part of this
            array contains the upper quasi-triangular matrix S in
            Schur canonical form from the Shur factorization of A.
            The contents of array A is not modified if fact = 'F'.
        Q : input rank-2 array('d'), shape  (n,n)
            On entry, if fact = 'F', then the leading n-by-n part of
            this array must contain the orthogonal matrix Q of the
            Schur factorization of A.
            Otherwise, Q need not be set on entry.
            On exit, the leading n-by-n part of this array contains
            the orthogonal matrix Q of the Schur factorization of A.
            The contents of array Q is not modified if fact = 'F'.
        B : input rank-2 array('d'), shape  (m,n)
            On entry, if trans = 'N', the leading m-by-n part of this
            array must contain the coefficient matrix B of the
            equation.
            On entry, if trans = 'T', the leading N-by-m part of this
            array must contain the coefficient matrix B of the
            equation.
            On exit, the leading n-by-n part of this array contains
            the upper triangular Cholesky factor U of the solution
            matrix X of the problem, X = op(U)'*op(U).
            If m = 0 and n > 0, then U is set to zero.
        dico : input string(len=1)
            Specifies the type of Lyapunov equation to be solved as
            follows:
            = 'C':  Equation (1), continuous-time case;
            = 'D':  Equation (2), discrete-time case.

    Optional arguments
    ------------------

        fact := 'N' input string(len=1)
            Specifies whether or not the real Schur factorization
            of the matrix A is supplied on entry, as follows:
            = 'F':  On entry, A and Q contain the factors from the
            real Schur factorization of the matrix A;
            = 'N':  The Schur factorization of A will be computed
            and the factors will be stored in A and Q.
        trans := 'N' input string(len=1)
            Specifies the form of op(K) to be used, as follows:
            = 'N':  op(K) = K    (No transpose);
            = 'T':  op(K) = K**T (Transpose).
        ldwork := None input int
            The length of the array DWORK.
            If m > 0, ldwork >= max(1,4*n + min(m,n));
            If m = 0, ldwork >= 1.
            For optimum performance ldwork should sometimes be larger.

    Return objects
    ______________

        U : rank-2 array('d'), shape  (n,n)
            The leading n-by-n part of this array contains
            the upper triangular Cholesky factor U of the solution
            matrix X of the problem, X = op(U)'*op(U).
        scale : float
            The scale factor, scale, set less than or equal to 1 to
            prevent the solution overflowing.
        w : rank-1 array('c'), shape (n)
            If fact = 'N', this array contains the eigenvalues of A.

    Raises
    ______

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             = 1:  if the Lyapunov equation is (nearly) singular
                (warning indicator);
                if DICO = 'C' this means that while the matrix A
                (or the factor S) has computed eigenvalues with
                negative real parts, it is only just stable in the
                sense that small perturbations in A can make one or
                more of the eigenvalues have a non-negative real
                part;
                if DICO = 'D' this means that while the matrix A
                (or the factor S) has computed eigenvalues inside
                the unit circle, it is nevertheless only just
                convergent, in the sense that small perturbations
                in A can make one or more of the eigenvalues lie
                outside the unit circle;
                perturbed values were used to solve the equation;
             = 2:  if FACT = 'N' and DICO = 'C', but the matrix A is
                not stable (that is, one or more of the eigenvalues
                of A has a non-negative real part), or DICO = 'D',
                but the matrix A is not convergent (that is, one or
                more of the eigenvalues of A lies outside the unit
                circle); however, A will still have been factored
                and the eigenvalues of A returned in WR and WI.
             = 3:  if FACT = 'F' and DICO = 'C', but the Schur factor S
                supplied in the array A is not stable (that is, one
                or more of the eigenvalues of S has a non-negative
                real part), or DICO = 'D', but the Schur factor S
                supplied in the array A is not convergent (that is,
                one or more of the eigenvalues of S lies outside the
                unit circle);
             = 4:  if FACT = 'F' and the Schur factor S supplied in
                the array A has two or more consecutive non-zero
                elements on the first sub-diagonal, so that there is
                a block larger than 2-by-2 on the diagonal;
             = 5:  if FACT = 'F' and the Schur factor S supplied in
                the array A has a 2-by-2 diagonal block with real
                eigenvalues instead of a complex conjugate pair;
             = 6:  if FACT = 'N' and the LAPACK Library routine DGEES
                has failed to converge. This failure is not likely
                to occur. The matrix B will be unaltered but A will
                be destroyed.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico','fact', 'trans', 'n', 'm', 'a', 'lda'+hidden, 'q',
        'ldq'+hidden, 'b', 'ldb'+hidden, 'scale', 'wr'+hidden,
        'wi'+hidden, 'dwork'+hidden, 'ldwork', 'info'+hidden]
    if ldwork is None:
        if m > 0:
            ldwork = max(1,4*n + min(m,n))
        elif m == 0:
            ldwork = 1
    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be either D or C')
    out = _wrapper.sb03od(dico,n,m,A,Q,B,fact=fact,trans=trans,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        if dico == 'D':
            error_text = """this means that while the matrix A
                (or the factor S) has computed eigenvalues inside
                the unit circle, it is nevertheless only just
                convergent, in the sense that small perturbations
                in A can make one or more of the eigenvalues lie
                outside the unit circle;
                perturbed values were used to solve the equation;"""
        else:
            error_text = """this means that while the matrix A
                (or the factor S) has computed eigenvalues with
                negative real parts, it is only just stable in the
                sense that small perturbations in A can make one or
                more of the eigenvalues have a non-negative real
                part;"""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        if dico == 'D':
            error_text = """the matrix A is not convergent (that is, one or
                more of the eigenvalues of A lies outside the unit
                circle); however, A will still have been factored
                and the eigenvalues of A returned in WR and WI."""
        else:
            error_text = """the matrix A is
                not stable (that is, one or more of the eigenvalues
                of A has a non-negative real part)."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        if dico == 'D':
            error_text = """the Schur factor S
                supplied in the array A is not convergent (that is,
                one or more of the eigenvalues of S lies outside the
                unit circle)."""
        else:
            error_text = """the Schur factor S
                supplied in the array A is not stable (that is, one
                or more of the eigenvalues of S has a non-negative
                real part)."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 4:
        if fact == 'F':
            error_text = """the Schur factor S supplied in
                the array A has two or more consecutive non-zero
                elements on the first sub-diagonal, so that there is
                a block larger than 2-by-2 on the diagonal."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 5:
        if fact == 'F':
            error_text = """the Schur factor S supplied in
                the array A has a 2-by-2 diagonal block with real
                eigenvalues instead of a complex conjugate pair."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 6:
        if fact == 'N':
            error_text = """the LAPACK Library routine DGEES
                has failed to converge. This failure is not likely
                to occur. The matrix B will be unaltered but A will
                be destroyed."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    U,scale,wr,wi = out[:-1]
    w = _np.zeros(n,'complex64')
    w.real = wr[0:n]
    w.imag = wi[0:n]
    return U,scale,w

def sb04md(n,m,A,B,C,ldwork=None):
    """X = sb04md(n,m,A,B,C[,ldwork])

    To solve for X the continuous-time Sylvester equation

     AX + XB = C

    where A, B, C and X are general n-by-n, m-by-m, n-by-m and
    n-by-m matrices respectively.

    Required arguments
    ------------------

        n : input int
        m : input int
        A : input rank-2 array('d'), shape  (n,n)
        B : input rank-2 array('d'), shape  (m,m)
        C : input rank-2 array('d'), shape  (n,m)

    Return objects
    --------------

        X : rank-2 array('d'), shape  (n,m)

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value;
             > 0:  if info = i, 1 <= i <= m, the QR algorithm failed to
                compute all the eigenvalues of B (see LAPACK Library
                routine DGEES)
             > m:  if a singular matrix was encountered whilst solving
                for the (info-m)-th column of matrix X.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'A', 'LDA'+hidden,  'B', 'LDB'+hidden, 'C',
        'LDC'+hidden,  'Z', 'LDZ'+hidden, 'IWORK'+hidden, 'DWORK'+hidden,
        'ldwork', 'INFO'+hidden]
    if ldwork is None:
        out = _wrapper.sb04md(n,m,A,B,C)
    else:
        out = _wrapper.sb04md(n,m,A,B,C,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value:"+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] > 0 and out[-1] <= m:
        error_text = """The QR algorithm failed to compute all the eigenvalues
(see LAPACK Library routine DGEES)"""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] > m:
        error_text = """a singular matrix was encountered whilst solving
for the %i-th column of matrix X.""" % (out[-1]-m)
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    return out[2]

def sb04qd(n,m,A,B,C,ldwork=None):
    """X = sb04qd(n,m,A,B,C[,ldwork])

    To solve for X the discrete-time Sylvester equation

        AXB + X + C = 0,

    where A, B, C and X are general n-by-n, m-by-m, n-by-m and
    n-by-m matrices respectively. A Hessenberg-Schur method, which
    reduces A to upper Hessenberg form, H = U'AU, and B' to real
    Schur form, S = Z'B'Z (with U, Z orthogonal matrices), is used.

    Required arguments
    ------------------

        n : input int
        m : input int
        A : input rank-2 array('d'), shape  (n,n)
        B : input rank-2 array('d'), shape  (m,m)
        C : input rank-2 array('d'), shape  (n,m)

    Return objects
    --------------

        X : rank-2 array('d'), shape  (n,m)

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal value
             > 0:  if info = i, 1 <= i <= m, the QR algorithm failed to
                compute all the eigenvalues of B (see LAPACK Library
                routine DGEES)
             > m:  if a singular matrix was encountered whilst solving
                for the (info-m)-th column of matrix X.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'A', 'LDA'+hidden,  'B', 'LDB'+hidden, 'C',
        'LDC'+hidden,  'Z', 'LDZ'+hidden, 'IWORK'+hidden, 'DWORK'+hidden,
        'ldwork', 'INFO'+hidden]
    if ldwork is None:
        out = _wrapper.sb04qd(n,m,A,B,C)
    else:
        out = _wrapper.sb04qd(n,m,A,B,C,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value:"+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] > 0 and out[-1] <= m:
        error_text = """The QR algorithm failed to compute all the eigenvalues
(see LAPACK Library routine DGEES)"""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] > m:
        error_text = """a singular matrix was encountered whilst solving
for the %i-th column of matrix X.""" % (out[-1]-m)
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    return out[2]

def sb10ad(n,m,np,ncon,nmeas,gamma,A,B,C,D,job=3,gtol=0.0,actol=0.0,liwork=None,ldwork=None):
    """ gamma_est, Ak, Bk, Ck, Dk, Ac, Bc, Cc, Dc, rcond = sb10ad(n,m,np,ncon,nmeas,gamma,A,B,C,D,[job,gtol,actol,liwork,ldwork])

    To compute the matrices of an H-infinity optimal n-state
    controller

           | Ak | Bk |
       K = |----|----|,
           | Ck | Dk |

    using modified Glover's and Doyle's 1988 formulas, for the system

           | A  | B1  B2  |   | A | B |
       P = |----|---------| = |---|---|
           | C1 | D11 D12 |   | C | D |
           | C2 | D21 D22 |

    and for the estimated minimal possible value of gamma with respect
    to gtol, where B2 has as column size the number of control inputs
    (ncon) and C2 has as row size the number of measurements (nmeas)
    being provided to the controller, and then to compute the matrices
    of the closed-loop system

           | AC | BC |
       G = |----|----|,
           | CC | DC |

    if the stabilizing controller exists.

    It is assumed that

    (A1) (A,B2) is stabilizable and (C2,A) is detectable,

    (A2) D12 is full column rank and D21 is full row rank,

    (A3) | A-j*omega*I  B2  | has full column rank for all omega,
         |    C1        D12 |

    (A4) | A-j*omega*I  B1  |  has full row rank for all omega.
         |    C2        D21 |


    Required arguments
    ------------------

        n : int
            The order of the system. (size of matrix A).
        m : int
            The column size of the matrix B.
        np : int
            The row size of the matrix C.
        ncon : int
            The number of control inputs.  m >= ncon >= 0, np-nmeas >= ncon.
        nmeas : int
            The number of measurements.  np >= nmeas >= 0, m-ncon >= nmeas.
        gamma : double
            The initial value of gamma on input.  It is assumed that gamma
            is sufficiently large so that the controller is admissible.  gamma >= 0.
        A : rank-2 array('d'), shape (n,n)
        B : rank-2 array('d'), shape (n,m)
        C : rank-2 array('d'), shape (np,n)
        D : rank-2 array('d'), shape (np,m)

    Optional arguments
    ------------------

        job := 3 int
            Specifies the computation to be performed, as follows:
            = 1: Use bisection method for decreasing gamma until the closed-loop
                system leaves stability
            = 2: Scan from gamma to 0 trying to find the minimal gamma for which
                the closed-loop system retains stability
            = 3: First bisection, then scanning.
            = 4: Find suboptimal controller only.
        gtol : double
            Tolerance used for controlling the accuracy of gamma
            and its distance to the estimated minimal possible
            value of gamma.
            If gtol <= 0, then a default value equal to sqrt(eps)
            is used, where eps is the relative machine precision.
        actol : double
            Upper bound for the poles of the closed-loop system used for determining
            if it is stable.  actol <= 0 for stable systems
        liwork : int
            The dimension of the integer cache array.
        ldwork : int
            The dimension of the double cache array.

    Return objects
    --------------

        gamma_est : double
            The minimal estimated gamma.
        Ak : rank-2 array('d'), shape (n,n)
            The controller state matrix Ak.
        Bk : rank-2 array('d') with bound s(n,nmeas)
            The controller input matrix Bk.
        Ck : rank-2 array('d'), shape  (ncon,n)
            The controller output matrix Ck.
        Dk : rank-2 array('d'), shape  (ncon,nmeas)
            The controller input/output matrix DK.
        Ac : rank-2 array('d'), shape  (2n,2n)
            The closed-loop system state matrix AC.
        Bc : rank-2 array('d'), shape  (2n,m-ncon)
            The closed-loop system input matrix BC.
        Cc : rank-2 array('d'), shape  (np-nmeas,2n)
            The closed-loop system output matrix CC.
        Dc : rank-2 array('d'), shape  (np-nmeas,m-ncon)
            The the closed-loop system input/output matrix DC.
        rcond : rank-1 array('d'), shape  (4)
           For the last successful step:
           rcond(1) contains the reciprocal condition number of the
                control transformation matrix;
           rcond(2) contains the reciprocal condition number of the
                measurement transformation matrix;
           rcond(3) contains an estimate of the reciprocal condition
                number of the X-Riccati equation;
           rcond(4) contains an estimate of the reciprocal condition
                number of the Y-Riccati equation.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal
                   value;
             = 1:  if the matrix | A-j*omega*I  B2  | had not full
                                 |    C1        D12 |
                   column rank in respect to the tolerance eps;
             = 2:  if the matrix | A-j*omega*I  B1  |  had not full row
                                 |    C2        D21 |
                   rank in respect to the tolerance eps;
             = 3:  if the matrix D12 had not full column rank in
                   respect to the tolerance SQRT(eps);
             = 4:  if the matrix D21 had not full row rank in respect
                   to the tolerance SQRT(eps);
             = 5:  if the singular value decomposition (SVD) algorithm
                   did not converge (when computing the SVD of one of
                   the matrices |A   B2 |, |A   B1 |, D12 or D21);
                                |C1  D12|  |C2  D21|
             = 6:  if the controller is not admissible (too small value
                   of gamma);
             = 7:  if the X-Riccati equation was not solved
                   successfully (the controller is not admissible or
                   there are numerical difficulties);
             = 8:  if the Y-Riccati equation was not solved
                   successfully (the controller is not admissible or
                   there are numerical difficulties);
             = 9:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is
                   zero [3];
             = 10: if there are numerical problems when estimating
                   singular values of D1111, D1112, D1111', D1121';
             = 11: if the matrices Inp2 - D22*DK or Im2 - DK*D22
                   are singular to working precision;
             = 12: if a stabilizing controller cannot be found.

    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['job', 'n', 'm', 'np', 'ncon', 'nmeas', 'gamma',
        'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'C', 'LDC'+hidden,
        'D', 'LDD'+hidden, 'AK', 'LDAK'+hidden, 'BK', 'LDBK'+hidden,
        'CK', 'LDCK'+hidden, 'DK', 'LDDK'+hidden, 'AC', 'LDAC'+hidden,
        'BC', 'LDBC'+hidden, 'CC', 'LDCC'+hidden, 'DC', 'LDDC'+hidden,
        'rcond', 'gtol', 'actol', 'IWORK'+hidden, 'liwork',
        'DWORK'+hidden, 'ldwork', 'BWORK'+hidden, 'LBWORK'+hidden, 'info']
    if liwork is None:
        liwork = max(2*max(n,m-ncon,np-nmeas,ncon,nmeas),n*n)
    if ldwork is None:
        m2 = ncon
        np2 = nmeas
        m1 = m - m2
        np1 = np - np2
        nd1 = np1 - m2
        nd2 = m1 - np2
        LW1 = n*m + np*n + np*m + m2*m2 + np2*np2
        LW2 = max((n + np1 + 1)*(n + m2) + max(3*(n + m2) + n + np1, 5*(n + m2)),\
            (n + np2)*(n + m1 + 1) + max(3*(n + np2) + n + m1, 5*(n + np2)),\
            m2 + np1*np1 + max(np1*max(n, m1), 3*m2 + np1, 5*m2),\
            np2 + m1*m1 + max(max(n, np1)*m1, 3*np2 + m1, 5*np2))
        LW3 = max(nd1*m1 + max(4*min(nd1, m1) + max(nd1, m1), 6*min(nd1, m1)),\
            np1*nd2 + max(4*min(np1, nd2) + max(np1, nd2), 6*min(np1, nd2)))
        LW4 = 2*m*m + np*np + 2*m*n + m*np + 2*n*np
        LW5 = 2*n*n + m*n + n*np
        LW6 = max(m*m + max(2*m1, 3*n*n + max(n*m, 10*n*n + 12*n + 5)),\
            np*np + max(2*np1, 3*n*n + max(n*np, 10*n*n + 12*n +5)))
        LW7 = m2*np2 + np2*np2 + m2*m2 + max(nd1*nd1 + max(2*nd1, (nd1 + nd2)*np2),\
            nd2*nd2 + max(2*nd2, nd2*m2), 3*n, n*(2*np2 + m2) +\
            max(2*n*m2, m2*np2 + max(m2*m2 + 3*m2, np2*(2*np2 + m2 + max(np2, n)))))
        ldwork = LW1 + max(1, LW2, LW3, LW4, LW5 + max(LW6,LW7))
    out = _wrapper.sb10ad(job,n,m,np,ncon,nmeas,gamma,A,B,C,D,gtol,actol,liwork,ldwork)

    if out[-1] != 0:
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "\
                +arg_list[-out[-1]-1]
        if out[-1] == 1:
            error_text = "The matrix [A-j*omega*I,  B2 ; C1,  D12] does not\
                have full column rank with respect to the tolerance eps."
        if out[-1] == 2:
            error_text = "The matrix [A-j*omega*I,  B1 ; C2,  D21] does not\
                have full row rank with respect to the tolerance eps."
        if out[-1] == 3:
            error_text = "The matrix D12 does not have full column rank with\
                respect to tolerance sqrt(eps)."
        if out[-1] == 4:
            error_text = "The matrix D21 does not have full column rank with\
                respect to tolerance sqrt(eps)."
        if out[-1] == 5:
            error_text = "The singular value decomposition (SVD) algorithm did\
                not converge (when computing the SVD of one of the matrices\
                [A,  B2; C1,  D12], [A,  B1; C2,  D21] , D12 or D21."
        if out[-1] == 6:
            error_text = "The controller is not admissible (too small value of\
                gamma)."
        if out[-1] == 7:
            error_text = "The X-Riccati equation was not solved successfully\
                (the controller is not admissible or there are numerical\
                difficulties)."
        if out[-1] == 8:
            error_text = "The Y-Riccati equation was not solved successfully\
                (the controller is not admissible or there are numerical\
                difficulties)."
        if out[-1] == 9:
            error_text = "The determinant of Im2 + Tu*D11Hat*Ty*D22 is zero,\
                see ref [3] in SLICOT doc."
        if out[-1] == 10:
            error_text = "There are numerical problems when estimating singular\
                values of D1111, D1112, D1111', D1121'."
        if out[-1] == 11:
            error_text = "The matrices Inp2 - D22*DK or Im2 - DK*D22 are singular\
                to working precision."
        if out[-1] == 12:
            error_text = "A stabilizing controller cannot be found."
        e = ValueError(error_text)
        e.info = out[-1]
        raise e

    return out[:-1]

def sb10dd(n,m,np,ncon,nmeas,gamma,A,B,C,D,tol=0.0,ldwork=None):
    """ gamma_est, Ak, Bk, Ck, Dk, rcond = sb10dd(n,m,np,ncon,nmeas,gamma,A,B,C,D,[tol,ldwork])

    To compute the matrices of an H-infinity (sub)optimal n-state
    controller

            | AK | BK |
        K = |----|----|,
            | CK | DK |

     for the discrete-time system

            | A  | B1  B2  |   | A | B |
        P = |----|---------| = |---|---|
            | C1 | D11 D12 |   | C | D |
            | C2 | D21 D22 |

     and for a given value of gamma, where B2 has as column size the
     number of control inputs (NCON) and C2 has as row size the number
     of measurements (NMEAS) being provided to the controller.

     It is assumed that

     (A1) (A,B2) is stabilizable and (C2,A) is detectable,

     (A2) D12 is full column rank and D21 is full row rank,

               j*Theta
     (A3) | A-e       *I  B2  | has full column rank for all
          |    C1         D12 |

          0 <= Theta < 2*Pi ,

               j*Theta
     (A4) | A-e       *I  B1  |  has full row rank for all
          |    C2         D21 |

          0 <= Theta < 2*Pi .

     Required arguments
     ------------------

        n : int
            The order of the system. (size of matrix A).
        m : int
            The column size of the matrix B.
        np : int
            The row size of the matrix C.
        ncon : int
            The number of control inputs.  m >= ncon >= 0, np-nmeas >= ncon.
        nmeas : int
            The number of measurements.  np >= nmeas >= 0, m-ncon >= nmeas.
        gamma : double
            The initial value of gamma on input.  It is assumed that gamma
            is sufficiently large so that the controller is admissible.  gamma >= 0.
        A : rank-2 array('d'), shape (n,n)
        B : rank-2 array('d'), shape (n,m)
        C : rank-2 array('d'), shape (np,n)
        D : rank-2 array('d'), shape (np,m)

    Optional arguments
    ------------------

        tol : double
          Tolerance used in neglecting the small singular values
          in rank determination. If tol <= 0, then a default value
          equal to 1000*eps is used, where eps is the relative
          machine precision.

        ldwork : int
          The dimension of the array dwork.

    Return objects
    --------------

        gamma_est : double
            The minimal estimated gamma.
        Ak : rank-2 array('d'), shape (n,n)
            The controller state matrix Ak.
        Bk : rank-2 array('d') with bound s(n,nmeas)
            The controller input matrix Bk.
        Ck : rank-2 array('d'), shape  (ncon,n)
            The controller output matrix Ck.
        Dk : rank-2 array('d'), shape  (ncon,nmeas)
            The controller input/output matrix DK.
        X  : rank-2 array('d'), shape (n,n)
            The matrix X, solution of the X-Riccati equation.
        Z  : rank-2 array('d'), shape (n,n)
            The matrix Z, solution of the Z-Riccati equation.
        rcond : rank-1 array('d'), shape  (8)
            rcond contains estimates of the reciprocal condition
            numbers of the matrices which are to be inverted and
            estimates of the reciprocal condition numbers of the
            Riccati equations which have to be solved during the
            computation of the controller. (See the description of
            the algorithm in [2].)
            rcond(1) contains the reciprocal condition number of the
                   matrix R3;
            rcond(2) contains the reciprocal condition number of the
                   matrix R1 - R2'*inv(R3)*R2;
            rcond(3) contains the reciprocal condition number of the
                   matrix V21;
            rcond(4) contains the reciprocal condition number of the
                   matrix St3;
            rcond(5) contains the reciprocal condition number of the
                   matrix V12;
            rcond(6) contains the reciprocal condition number of the
                   matrix Im2 + dkhat*D22
            rcond(7) contains the reciprocal condition number of the
                   X-Riccati equation;
            rcond(8) contains the reciprocal condition number of the
                   Z-Riccati equation.


    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if INFO = -i, the i-th argument had an illegal
                   value;
                                      j*Theta
             = 1:  if the matrix | A-e       *I  B2  | had not full
                                 |      C1       D12 |
                   column rank;
                                      j*Theta
             = 2:  if the matrix | A-e       *I  B1  | had not full
                                 |      C2       D21 |
                   row rank;
             = 3:  if the matrix D12 had not full column rank;
             = 4:  if the matrix D21 had not full row rank;
             = 5:  if the controller is not admissible (too small value
                   of gamma);
             = 6:  if the X-Riccati equation was not solved
                   successfully (the controller is not admissible or
                   there are numerical difficulties);
             = 7:  if the Z-Riccati equation was not solved
                   successfully (the controller is not admissible or
                   there are numerical difficulties);
             = 8:  if the matrix Im2 + DKHAT*D22 is singular.
             = 9:  if the singular value decomposition (SVD) algorithm
                   did not converge (when computing the SVD of one of
                   the matrices |A   B2 |, |A   B1 |, D12 or D21).
                                |C1  D12|  |C2  D21|

    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'np', 'ncon', 'nmeas', 'gamma',
        'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'C', 'LDC'+hidden,
        'D', 'LDD'+hidden, 'AK', 'LDAK'+hidden, 'BK', 'LDBK'+hidden,
        'CK', 'LDCK'+hidden, 'DK', 'LDDK'+hidden, 'X', 'LDX'+hidden, 'Z', 'LDZ'+hidden,
        'rcond', 'tol', 'IWORK'+hidden,
        'DWORK'+hidden, 'ldwork', 'BWORK'+hidden, 'info']
    if ldwork is None:
        m1 = m - ncon
        m2 = ncon
        np1 = np - nmeas
        np2 = nmeas
        LW1 = (n+np1+1)*(n+m2) + max(3*(n+m2)+n+np1,5*(n+m2))
        LW2 = (n+np2)*(n+m1+1) + max(3*(n+np2)+n+m1,5*(n+np2))
        LW3 = 13*n*n + 2*m*m + n*(8*m+np2) + m1*(m2+np2) + 6*n + max(14*n+23,16*n,2*n+m,3*m)
        LW4 = 13*n*n + m*m + (8*n+m+m2+2*np2)*(m2+np2) + 6*n + n*(m+np2) + max(14*n+23,16*n,2*n+m2+np2,3*(m2+np2))
        ldwork = max(LW1,LW2,LW3,LW4)
    out = _wrapper.sb10dd(n,m,np,ncon,nmeas,gamma,A,B,C,D,tol,ldwork)
    
    if out[-1] != 0:
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "\
                +arg_list[-out[-1]-1]
        if out[-1] == 1:
            error_text = "  j*Theta\
            The matrix | A-e       *I  B2  | had not full column rank.\
                       |      C1       D12 |"
        if out[-1] == 2:
            error_text = "  j*Theta\
            The matrix | A-e       *I  B2  | had not full row rank.\
                       |      C1       D12 |"
        if out[-1] == 3:
            error_text = "The matrix D12 had not full column rank."
        if out[-1] == 4:
            error_text = "The matrix D21 had not full row rank."
        if out[-1] == 5:
            error_text = "The controller is not admissible (too small value of gamma)"
        if out[-1] == 6:
            error_text = "The X-Riccati equation was not solved\
                successfully (the controller is not admissible or\
                there are numerical difficulties)."
        if out[-1] == 7:
            error_text = "The Z-Riccati equation was not solved\
                successfully (the controller is not admissible or\
                there are numerical difficulties)."
        if out[-1] == 8:
            error_text = "The matrix Im2 + DKHAT*D22 is singular."
        if out[-1] == 9:
            error_text = "the singular value decomposition (SVD) algorithm\
                did not converge (when computing the SVD of one of\
                the matrices |A   B2 |, |A   B1 |, D12 or D21).\
                             |C1  D12|  |C2  D21|"
        e = ValueError(error_text)
        e.info = out[-1]
        raise e

    return out[:-1]

def sb10hd(n,m,np,ncon,nmeas,A,B,C,D,tol=0.0,ldwork=None):
    """ Ak, Bk, Ck, Dk, rcond = sb10hd(n,m,np,ncon,nmeas,a,b,c,d,[tol,ldwork])

    To compute the matrices of the H2 optimal n-state controller

           | AK | BK |
       K = |----|----|
           | CK | DK |

    for the system

                | A  | B1  B2  |   | A | B |
            P = |----|---------| = |---|---| ,
                | C1 |  0  D12 |   | C | D |
                | C2 | D21 D22 |

    where B2 has as column size the number of control inputs (ncon)
    and C2 has as row size the number of measurements (nmeas) being
    provided to the controller.

    It is assumed that

    (A1) (A,B2) is stabilizable and (C2,A) is detectable,

    (A2) The block D11 of D is zero,

    (A3) D12 is full column rank and D21 is full row rank.

    Required arguments
    ------------------

        n : int
            The order of the system. (size of matrix A).
        m : int
            The column size of the matrix B
        np : int
            The row size of the matrix C
        ncon : int
            The number of control inputs.  m >= ncon >= 0, np-nmeas >= ncon.
        nmeas : int
            The number of measurements.  np >= nmeas >= 0, m-ncon >= nmeas.
        A : rank-2 array('d'), shape (n,n)
        B : rank-2 array('d'), shape (n,m)
        C : rank-2 array('d'), shape (np,n)
        D : rank-2 array('d'), shape (np,m)

    Optional arguments
    ------------------

        tol : double
            Tolerance used for controlling the accuracy of the applied
            transformations for computing the normalized form in
            SLICOT Library routine SB10UD. Transformation matrices
            whose reciprocal condition numbers are less than tol are
            not allowed. If tol <= 0, then a default value equal to
            sqrt(eps) is used, where eps is the relative machine
            precision.
        ldwork : int
            The dimension of the cache array.

    Return objects
    --------------

        Ak : rank-2 array('d'), shape (n,n)
            The controller state matrix Ak.
        Bk : rank-2 array('d'), shape (n,nmeas)
            The controller input matrix Bk.
        Ck : rank-2 array('d'), shape (ncon,n)
            The controller output matrix Ck.
        Dk : rank-2 array('d'), shape (ncon,nmeas)
            The controller input/output matrix Dk.
        rcond : rank-1 array('d'), shape (4)
           For the last successful step:
           rcond(1) contains the reciprocal condition number of the
                control transformation matrix;
           rcond(2) contains the reciprocal condition number of the
                measurement transformation matrix;
           rcond(3) contains an estimate of the reciprocal condition
                number of the X-Riccati equation;
           rcond(4) contains an estimate of the reciprocal condition
                number of the Y-Riccati equation.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal
                   value;
             = 1:  if the matrix D12 had not full column rank in
                   respect to the tolerance tol;
             = 2:  if the matrix D21 had not full row rank in respect
                   to the tolerance tol;
             = 3:  if the singular value decomposition (SVD) algorithm
                   did not converge (when computing the SVD of one of
                   the matrices D12 or D21).
             = 4:  if the X-Riccati equation was not solved
                   successfully;
             = 5:  if the Y-Riccati equation was not solved
                   successfully.
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'np', 'ncon', 'nmeas', 'A', 'LDA'+hidden,
        'B', 'LDB'+hidden, 'C', 'LDC'+hidden, 'D', 'LDD'+hidden, 'Ak',
        'LDAK'+hidden, 'Bk', 'LDBK'+hidden, 'Ck', 'LDCK'+hidden, 'Dk',
        'LDDK'+hidden, 'rcond', 'tol', 'IWORK'+hidden, 'DWORK'+hidden,
        'LDWORK', 'BWORK'+hidden, 'INFO']
    if ldwork is None:
        Q = max(m-ncon,ncon,np-nmeas,nmeas)
        ldwork = 2*Q*(3*Q+2*n)+max(1,Q*(Q+max(n,5)+1),n*(14*n+12+2*Q)+5)
    out = _wrapper.sb10hd(n,m,np,ncon,nmeas,A,B,C,D,tol,ldwork)

    if out[-1] != 0:
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "\
                +arg_list[-out[-1]-1]
        if out[-1] == 1:
            error_text = "The matrix D12 does not have full column rank with\
                respect to the tolerance tol."
        if out[-1] == 2:
            error_text = "The matrix D21 does not have full row rank with\
                respect to the tolerance tol."
        if out[-1] == 3:
            error_text = "The singular value decomposition (SVD) algorithm\
                did not converge (when computing the SVD of one of the matrices\
                D12 or D21.)"
        if out[-1] == 4:
            error_text = "The X-Riccati equation was not solved successfully\
                (the controller is not admissible or there are numerical difficulties)."
        if out[-1] == 5:
            error_text = "The Y-Riccati equation was not solved successfully\
                (the controller is not admissible or there are numerical\
                difficulties)."
        e = ValueError(error_text)
        e.info = out[-1]
        raise e

    return out[:-1]

def sg03ad(dico,job,fact,trans,uplo,N,A,E,Q,Z,X,ldwork=None):
    """ A,E,Q,Z,X,scale,sep,ferr,alphar,alphai,beta =
            sg03ad(dico,job,fact,trans,uplo,N,A,E,Q,Z,X,[ldwork])

    To solve for X either the generalized continuous-time Lyapunov
    equation

          T                T
     op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1)

    or the generalized discrete-time Lyapunov equation

          T                T
     op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2)

    where op(M) is either M or M**T for M = A, E and the right hand
    side Y is symmetric. A, E, Y, and the solution X are N-by-N
    matrices. SCALE is an output scale factor, set to avoid overflow
    in X.

    Estimates of the separation and the relative forward error norm
    are provided.

    Required arguments
    ------------------

        dico : input string(len=1)
            Specifies which type of the equation is considered:
              = 'C':  Continuous-time equation (1);
              = 'D':  Discrete-time equation (2).

        job : input string(len=1)
            Specifies if the solution is to be computed and if the
            separation is to be estimated:
              = 'X':  Compute the solution only;
              = 'S':  Estimate the separation only;
              = 'B':  Compute the solution and estimate the separation.

        fact : input string(len=1)
            Specifies whether the generalized real Schur
            factorization of the pencil A - lambda * E is supplied
            on entry or not:
              = 'N':  Factorization is not supplied;
              = 'F':  Factorization is supplied.

        trans : input string(len=1)
            Specifies whether the transposed equation is to be solved
            or not:
              = 'N':  op(A) = A,    op(E) = E;
              = 'T':  op(A) = A**T, op(E) = E**T.

        uplo : input string(len=1)
            Specifies whether the lower or the upper triangle of the
            array X is needed on input:
              = 'L':  Only the lower triangle is needed on input;
              = 'U':  Only the upper triangle is needed on input.

        A : input rank-2 array('d') with bounds (n,n)
            On entry, if FACT = 'F', then the leading N-by-N upper
            Hessenberg part of this array must contain the
            generalized Schur factor A_s of the matrix A (see
            definition (3) in section METHOD). A_s must be an upper
            quasitriangular matrix. The elements below the upper
            Hessenberg part of the array A are not referenced.
            If FACT = 'N', then the leading N-by-N part of this
            array must contain the matrix A.
            On exit, the leading N-by-N part of this array contains
            the generalized Schur factor A_s of the matrix A. (A_s is
            an upper quasitriangular matrix.)

        E : input rank-2 array('d') with bounds (n,n)
            On entry, if FACT = 'F', then the leading N-by-N upper
            triangular part of this array must contain the
            generalized Schur factor E_s of the matrix E (see
            definition (4) in section METHOD). The elements below the
            upper triangular part of the array E are not referenced.
            If FACT = 'N', then the leading N-by-N part of this
            array must contain the coefficient matrix E of the
            equation.
            On exit, the leading N-by-N part of this array contains
            the generalized Schur factor E_s of the matrix E. (E_s is
            an upper triangular matrix.)

        Q : input rank-2 array('d') with bounds (n,n)
            On entry, if FACT = 'F', then the leading N-by-N part of
            this array must contain the orthogonal matrix Q from
            the generalized Schur factorization (see definitions (3)
            and (4) in section METHOD).
            If FACT = 'N', Q need not be set on entry.
            On exit, the leading N-by-N part of this array contains
            the orthogonal matrix Q from the generalized Schur
            factorization.

        Z : input rank-2 array('d') with bounds (n,n)
            On entry, if FACT = 'F', then the leading N-by-N part of
            this array must contain the orthogonal matrix Z from
            the generalized Schur factorization (see definitions (3)
            and (4) in section METHOD).
            If FACT = 'N', Z need not be set on entry.
            On exit, the leading N-by-N part of this array contains
            the orthogonal matrix Z from the generalized Schur
            factorization.

        X : input rank-2 array('d') with bounds (n,n)
            On entry, if JOB = 'B' or 'X', then the leading N-by-N
            part of this array must contain the right hand side matrix
            Y of the equation. Either the lower or the upper
            triangular part of this array is needed (see mode
            parameter UPLO).
            If JOB = 'S', X is not referenced.
            On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then
            the leading N-by-N part of this array contains the
            solution matrix X of the equation.
            If JOB = 'S', X is not referenced.


    Optional arguments
    ------------------

        ldwork := max(1,max(2*n*n,4*n)) input int
            The length of the array DWORK. The following table
              contains the minimal work space requirements depending
              on the choice of JOB and FACT.

                     JOB        FACT    |  LDWORK
                     -------------------+-------------------
                     'X'        'F'     |  MAX(1,N)
                     'X'        'N'     |  MAX(1,4*N)
                     'B', 'S'   'F'     |  MAX(1,2*N**2)
                     'B', 'S'   'N'     |  MAX(1,2*N**2,4*N)

              For optimum performance, LDWORK should be larger.


    Return objects
    --------------

        A : rank-2 array('d') with bounds (n,n)
            On exit, the leading N-by-N part of this array contains
            the generalized Schur factor A_s of the matrix A. (A_s is
            an upper quasitriangular matrix.)

        E : rank-2 array('d') with bounds (n,n)
            On exit, the leading N-by-N part of this array contains
            the generalized Schur factor E_s of the matrix E. (E_s is
            an upper triangular matrix.)

        Q : rank-2 array('d') with bounds (n,n)
            On exit, the leading N-by-N part of this array contains
            the orthogonal matrix Q from the generalized Schur
            factorization.

        Z : rank-2 array('d') with bounds (n,n)
            On exit, the leading N-by-N part of this array contains
            the orthogonal matrix Z from the generalized Schur
            factorization.

        X : rank-2 array('d') with bounds (n,n)
            On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then
            the leading N-by-N part of this array contains the
            solution matrix X of the equation.
            If JOB = 'S', X is not referenced.

        scale : float
            The scale factor set to avoid overflow in X.
            (0 < SCALE <= 1)
        sep : float
            If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then
            SEP contains an estimate of the separation of the
            Lyapunov operator.
        ferr : float
            If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an
            estimated forward error bound for the solution X. If XTRUE
            is the true solution, FERR estimates the relative error
            in the computed solution, measured in the Frobenius norm:
            norm(X - XTRUE) / norm(XTRUE)

        alphar : rank-1 array('d') with bounds (n)
        alphai : rank-1 array('d') with bounds (n)
        beta : rank-1 array('d') with bounds (n)
            If FACT = 'N' and INFO = 0, 3, or 4, then
            (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the
            eigenvalues of the matrix pencil A - lambda * E.
            If FACT = 'F', ALPHAR, ALPHAI, and BETA are not
            referenced.

    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if info = -i, the i-th argument had an illegal
                   value;
             = 1:  FACT = 'F' and the matrix contained in the upper
                   Hessenberg part of the array A is not in upper
                   quasitriangular form;
             = 2:  FACT = 'N' and the pencil A - lambda * E cannot be
                   reduced to generalized Schur form: LAPACK routine
                   DGEGS has failed to converge;
             = 3:  DICO = 'D' and the pencil A - lambda * E has a
                   pair of reciprocal eigenvalues. That is, lambda_i =
                   1/lambda_j for some i and j, where lambda_i and
                   lambda_j are eigenvalues of A - lambda * E. Hence,
                   equation (2) is singular;  perturbed values were
                   used to solve the equation (but the matrices A and
                   E are unchanged);
             = 4:  DICO = 'C' and the pencil A - lambda * E has a
                   degenerate pair of eigenvalues. That is, lambda_i =
                   -lambda_j for some i and j, where lambda_i and
                   lambda_j are eigenvalues of A - lambda * E. Hence,
                   equation (1) is singular;  perturbed values were
                   used to solve the equation (but the matrices A and
                   E are unchanged).
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'fact', 'trans', 'uplo', 'N'+hidden, 'A', 'LDA'+hidden, 'E',
                'LDE'+hidden, 'Q', 'LDQ'+hidden, 'Z', 'LDZ'+hidden, 'X', 'LDX'+hidden,
                'scale', 'sep', 'ferr', 'alphar', 'alphai', 'beta', 'IWORK'+hidden,
                'DWORK'+hidden, 'ldwork', 'info' ]

    if ldwork is None:
        if job == 'X' and fact == 'F':
            ldwork = max(1,N)
        elif job == 'X' and fact == 'N':
            ldwork = max(1,8*N+16)
        elif (job == 'B' or job == 'S') and fact == 'F':
            ldwork = max(1,2*N**2)
        elif (job == 'B' or job == 'S') and fact == 'N':
            ldwork = max(1,2*N**2,8*N+16)

    out = _wrapper.sg03ad(dico,job,fact,trans,uplo,N,A,E,Q,Z,X,ldwork)

    if out[-1] != 0:
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "\
                +arg_list[-out[-1]-1]
        elif out[-1] == 1:
            error_text = "FACT = 'F' and the matrix contained in the upper \
                        Hessenberg part of the array A is not in upper \
                        quasitriangular form"
        elif out[-1] == 2:
            error_text = "FACT = 'N' and the pencil A - lambda * E cannot be \
                        reduced to generalized Schur form: LAPACK routine \
                        DGEGS has failed to converge"
        elif out[-1] == 3:
            error_text = "DICO = 'D' and the pencil A - lambda * E has a \
                        pair of reciprocal eigenvalues. That is, lambda_i = \
                        1/lambda_j for some i and j, where lambda_i and \
                        lambda_j are eigenvalues of A - lambda * E. Hence, \
                        equation (2) is singular;  perturbed values were \
                        used to solve the equation (but the matrices A and \
                        E are unchanged)"
        elif out[-1] == 4:
            error_text = "TDICO = 'C' and the pencil A - lambda * E has a \
                        degenerate pair of eigenvalues. That is, lambda_i = \
                        -lambda_j for some i and j, where lambda_i and \
                        lambda_j are eigenvalues of A - lambda * E. Hence, \
                        equation (1) is singular;  perturbed values were \
                        used to solve the equation (but the matrices A and \
                        E are unchanged)"

        e = ValueError(error_text)
        e.info = out[-1]
        raise e

    return out[:-1]

def sg02ad(dico,jobb,fact,uplo,jobl,scal,sort,acc,N,M,P,A,E,B,Q,R,L,ldwork=None,tol=-1):
    """rcondu,x,alfar,alfai,beta,s,t,u,iwarn,info =
    sg02ad(dico,jobb,fact,uplo,jobl,scal,sort,acc,N,M,P,A,E,B,Q,R,L,[ldwork,tol=-1])

    To solve for X either the continuous-time algebraic Riccati
    equation
                                -1
     Q + A'XE + E'XA - (L+E'XB)R  (L+E'XB)' = 0 ,              (1)

    or the discrete-time algebraic Riccati equation
                                     -1
     E'XE = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q ,         (2)

    where A, E, B, Q, R, and L are N-by-N, N-by-N, N-by-M, N-by-N,
    M-by-M and N-by-M matrices, respectively, such that Q = C'C,
    R = D'D and L = C'D; X is an N-by-N symmetric matrix.
    The routine also returns the computed values of the closed-loop
    spectrum of the system, i.e., the stable eigenvalues
    lambda(1),...,lambda(N) of the pencil (A - BF,E), where F is
    the optimal gain matrix,
          -1
     F = R  (L+E'XB)' ,        for (1),

    and
                 -1
     F = (R+B'XB)  (L+A'XB)' , for (2).
                           -1
    Optionally, matrix G = BR  B' may be given instead of B and R.
    Other options include the case with Q and/or R given in a
    factored form, Q = C'C, R = D'D, and with L a zero matrix.

    The routine uses the method of deflating subspaces, based on
    reordering the eigenvalues in a generalized Schur matrix pair.

    It is assumed that E is nonsingular, but this condition is not
    checked. Note that the definition (1) of the continuous-time
    algebraic Riccati equation, and the formula for the corresponding
    optimal gain matrix, require R to be nonsingular, but the
    associated linear quadratic optimal problem could have a unique
    solution even when matrix R is singular, under mild assumptions
    (see METHOD). The routine SG02AD works accordingly in this case.




    Required arguments
    ------------------

        dico : input string(len=1)
            Specifies the type of Riccati equation to be solved as
            follows:
              = 'C':  Equation (1), continuous-time case;
              = 'D':  Equation (2), discrete-time case.

        jobb : input string(len=1)
            Specifies whether or not the matrix G is given, instead
            of the matrices B and R, as follows:
              = 'B':  B and R are given;
              = 'G':  G is given.

        fact : input string(len=1)
            Specifies whether or not the matrices Q and/or R (if
            JOBB = 'B') are factored, as follows:
              = 'N':  Not factored, Q and R are given;
              = 'C':  C is given, and Q = C'C;
              = 'D':  D is given, and R = D'D;
              = 'B':  Both factors C and D are given, Q = C'C, R = D'D.

        uplo : input string(len=1)
            If JOBB = 'G', or FACT = 'N', specifies which triangle of
            the matrices G, or Q and R, is stored, as follows:
              = 'U':  Upper triangle is stored;
              = 'L':  Lower triangle is stored.

        jobl : input string(len=1)
            Specifies whether or not the matrix L is zero, as follows:
              = 'Z':  L is zero;
              = 'N':  L is nonzero.
              JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed.
              SLICOT Library routine SB02MT should be called just before
              SG02AD, for obtaining the results when JOBB = 'G' and
              JOBL = 'N'.

        scal : input string(len=1)
            If JOBB = 'B', specifies whether or not a scaling strategy
            should be used to scale Q, R, and L, as follows:
              = 'G':  General scaling should be used;
              = 'N':  No scaling should be used.
              SCAL is not used if JOBB = 'G'.

        sort : input string(len=1)
            Specifies which eigenvalues should be obtained in the top
            of the generalized Schur form, as follows:
              = 'S':  Stable   eigenvalues come first;
              = 'U':  Unstable eigenvalues come first.

        acc : input string(len=1)
            Specifies whether or not iterative refinement should be
            used to solve the system of algebraic equations giving
            the solution matrix X, as follows:
              = 'R':  Use iterative refinement;
              = 'N':  Do not use iterative refinement.

        N : input int
            The actual state dimension, i.e., the order of the
            matrices A, E, Q, and X, and the number of rows of the
            matrices B and L.  N > 0.
        M : input int
            The number of system inputs. If JOBB = 'B', M is the
            order of the matrix R, and the number of columns of the
            matrix B.  M >= 0.
            M is not used if JOBB = 'G'.

        P : input int            out = _wrapper.sg02ad_bn(dico,uplo,jobl,scal,sort,acc,N,M,A,E,B,Q,R,L,tol,ldwork)
            The number of system outputs. If FACT = 'C' or 'D' or 'B',
            P is the number of rows of the matrices C and/or D.
            P >= 0.
            Otherwise, P is not used.

        A : input rank-2 array('d') with bounds (max(1,N),N)
            The leading N-by-N part of this array must contain the
            state matrix A of the descriptor system.

        E : input rank-2 array('d') with bounds (max(1,N),N)
            The leading N-by-N part of this array must contain the
            matrix E of the descriptor system.

        B : input rank-2 array('d') with bounds (max(1,N),*)
            If JOBB = 'B', the leading N-by-M part of this array must
            contain the input matrix B of the system.
            If JOBB = 'G', the leading N-by-N upper triangular part
            (if UPLO = 'U') or lower triangular part (if UPLO = 'L')
            of this array must contain the upper triangular part or
            lower triangular part, respectively, of the matrix
                -1
            G = BR  B'. The stricly lower triangular part (if
            UPLO = 'U') or stricly upper triangular part (if
            UPLO = 'L') is not referenced.


        Q : input rank-2 array('d') with bounds (ldq,N)
            If FACT = 'N' or 'D', the leading N-by-N upper triangular
            part (if UPLO = 'U') or lower triangular part (if UPLO =
            'L') of this array must contain the upper triangular part
            or lower triangular part, respectively, of the symmetric
            state weighting matrix Q. The stricly lower triangular
            part (if UPLO = 'U') or stricly upper triangular part (if
            UPLO = 'L') is not referenced.
            If FACT = 'C' or 'B', the leading P-by-N part of this
            array must contain the output matrix C of the system.
            If JOBB = 'B' and SCAL = 'G', then Q is modified
            internally, but is restored on exit.

            The leading dimension of array Q.
            LDQ >= MAX(1,N) if FACT = 'N' or 'D';
            LDQ >= MAX(1,P) if FACT = 'C' or 'B'.

        R : input rank-2 array('d') with bounds (ldr,M)
            If FACT = 'N' or 'C', the leading M-by-M upper triangular
            part (if UPLO = 'U') or lower triangular part (if UPLO =
            'L') of this array must contain the upper triangular part
            or lower triangular part, respectively, of the symmetric
            input weighting matrix R. The stricly lower triangular
            part (if UPLO = 'U') or stricly upper triangular part (if
            UPLO = 'L') is not referenced.
            If FACT = 'D' or 'B', the leading P-by-M part of this
            array must contain the direct transmission matrix D of the
            system.
            If JOBB = 'B' and SCAL = 'G', then R is modified
            internally, but is restored on exit.
            If JOBB = 'G', this array is not referenced.

            The leading dimension of array R.
            LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C';
            LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B';
            LDR >= 1        if JOBB = 'G'.


        L : input rank-2 array('d') with bounds (n,M)
            If JOBL = 'N' and JOBB = 'B', the leading N-by-M part of
            this array must contain the cross weighting matrix L.
            If JOBB = 'B' and SCAL = 'G', then L is modified
            internally, but is restored on exit.
            If JOBL = 'Z' or JOBB = 'G', this array is not referenced.


    Optional arguments
    ------------------

        ldwork := max(7*(2*n+1)+16,16*n) input int
            The length of the array DWORK.
              LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G';
              LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'.
            For optimum performance LDWORK should be larger.

        tol := -1 input float
            The tolerance to be used to test for near singularity of
            the original matrix pencil, specifically of the triangular
            M-by-M factor obtained during the reduction process. If
            the user sets TOL > 0, then the given value of TOL is used
            as a lower bound for the reciprocal condition number of
            that matrix; a matrix whose estimated condition number is
            less than 1/TOL is considered to be nonsingular. If the
            user sets TOL <= 0, then a default tolerance, defined by
            TOLDEF = EPS, is used instead, where EPS is the machine
            precision (see LAPACK Library routine DLAMCH).
            This parameter is not referenced if JOBB = 'G'.




    Return objects
    --------------
        rcondu : float
            If N > 0 and INFO = 0 or INFO = 7, an estimate of the
            reciprocal of the condition number (in the 1-norm) of
            the N-th order system of algebraic equations from which
            the solution matrix X is obtained.

        X : rank-2 array('d') with bounds (n,n)
            If INFO = 0, the leading N-by-N part of this array
            contains the solution matrix X of the problem.

        alfar : rank-1 array('d') with bounds (2 * n)

        alfai : rank-1 array('d') with bounds (2 * n)

        beta : rank-1 array('d') with bounds (2 * n)
            The generalized eigenvalues of the 2N-by-2N matrix pair,
            ordered as specified by SORT (if INFO = 0, or INFO >= 5).
            For instance, if SORT = 'S', the leading N elements of
            these arrays contain the closed-loop spectrum of the
            system. Specifically,
             lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for
            k = 1,2,...,N.

        S : rank-2 array('d') with bounds (2 * n,2 * n)
            The leading 2N-by-2N part of this array contains the
            ordered real Schur form S of the first matrix in the
            reduced matrix pencil associated to the optimal problem,
            corresponding to the scaled Q, R, and L, if JOBB = 'B'
            and SCAL = 'G'. That is,

                 (S   S  )
                 ( 11  12)
             S = (       ),
                 (0   S  )
                 (     22)

            where S  , S   and S   are N-by-N matrices.
                 11   12      22
            Array S must have 2*N+M columns if JOBB = 'B', and 2*N
            columns, otherwise.


        T : rank-2 array('d') with bounds (2 * n,2 * n)
            The leading 2N-by-2N part of this array contains the
            ordered upper triangular form T of the second matrix in
            the reduced matrix pencil associated to the optimal
            problem, corresponding to the scaled Q, R, and L, if
            JOBB = 'B' and SCAL = 'G'. That is,

                 (T   T  )
                 ( 11  12)
             T = (       ),
                 (0   T  )
                 (     22)

            where T  , T   and T   are N-by-N matrices.
                 11   12      22

        U : rank-2 array('d') with bounds (2 * n,2 * n)
            The leading 2N-by-2N part of this array contains the right
            transformation matrix U which reduces the 2N-by-2N matrix
            pencil to the ordered generalized real Schur form (S,T).
            That is,

                 (U   U  )
                 ( 11  12)
             U = (       ),
                 (U   U  )
                 ( 21  22)

            where U  , U  , U   and U   are N-by-N matrices.
                 11   12   21      22
            If JOBB = 'B' and SCAL = 'G', then U corresponds to the
            scaled pencil. If a basis for the stable deflating
            subspace of the original problem is needed, then the
            submatrix U   must be multiplied by the scaling factor
                     21
            contained in DWORK(4).

        iwarn : int
            = 0:  no warning;
            = 1:  the computed solution may be inaccurate due to poor
                  scaling or eigenvalues too close to the boundary of
                  the stability domain (the imaginary axis, if
                  DICO = 'C', or the unit circle, if DICO = 'D').


    Raises
    ------

        ValueError : e
            e.info contains information about the exact type of exception
              < 0:  if INFO = -i, the i-th argument had an illegal
                    value;
              = 1:  if the computed extended matrix pencil is singular,
                    possibly due to rounding errors;
              = 2:  if the QZ algorithm failed;
              = 3:  if reordering of the generalized eigenvalues failed;
              = 4:  if after reordering, roundoff changed values of
                    some complex eigenvalues so that leading eigenvalues
                    in the generalized Schur form no longer satisfy the
                    stability condition; this could also be caused due
                    to scaling;
              = 5:  if the computed dimension of the solution does not
                    equal N;
              = 6:  if the spectrum is too close to the boundary of
                    the stability domain;
              = 7:  if a singular matrix was encountered during the
                    computation of the solution matrix X.
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'jobb', 'fact', 'uplo', 'jobl',
                'scal', 'sort', 'acc', 'N', 'M', 'P',
                'A', 'LDA'+hidden, 'E', 'LDE'+hidden, 'B', 'LDB'+hidden,
                'Q', 'LDQ'+hidden, 'R', 'LDR'+hidden, 'L', 'LDL'+hidden,
                'rcondu', 'X', 'LDX'+hidden, 'alfar', 'alfai',
                'beta', 'S', 'LDS'+hidden, 'T', 'LDT'+hidden, 'U',
                'LDU'+hidden, 'tol', 'IWORK'+hidden, 'DWORK'+hidden,
                'ldwork', 'BWORK'+hidden, 'iwarn', 'INFO'+hidden ]

    if ldwork is None:
        if (jobb == 'G'):
            ldwork = max(7*(2*N+1)+16,16*N)
        elif (jobb == 'B'):
            ldwork = max(7*(2*N+1)+16,16*N,2*N+M,3*M)

    if (jobb == 'G'):
        out = _wrapper.sg02ad_g(dico,uplo,sort,acc,N,A,E,B,Q,ldwork)
    elif (jobb == 'B'):
        if (fact == 'N'):
            out = _wrapper.sg02ad_bn(dico,uplo,jobl,scal,sort,acc,N,M,A,E,B,Q,R,L,tol,ldwork)
        elif (fact == 'C'):
            out = _wrapper.sg02ad_bc(dico,jobl,scal,sort,acc,N,M,P,A,E,B,Q,R,L,tol,ldwork)
        elif (fact == 'D'):
            out = _wrapper.sg02ad_bc(dico,jobl,scal,sort,acc,N,M,P,A,E,B,Q,R,L,tol,ldwork)
        elif (fact == 'B'):
            out = _wrapper.sg02ad_bb(dico,jobl,scal,sort,acc,N,M,P,A,E,B,Q,R,L,tol,ldwork)

    if out[-1] != 0:
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "\
                +arg_list[-out[-1]-1]
        elif out[-1] == 1:
            error_text = "The computed extended matrix pencil is singular,\
                possibly due to rounding errors"
        elif out[-1] == 2:
            error_text = "The QZ algorithm failed"
        elif out[-1] == 3:
            error_text = "Reordering of the generalized eigenvalues failed"
        elif out[-1] == 4:
            error_text = "After reordering, roundoff changed values of\
                some complex eigenvalues so that leading eigenvalues\
                in the generalized Schur form no longer satisfy the\
                stability condition; this could also be caused due\
                to scaling"
        elif out[-1] == 5:
            error_text = "The computed dimension of the solution does not\
                equal N"
        elif out[-1] == 6:
            error_text = "The spectrum is too close to the boundary of\
                the stability domain"
        elif out[-1] == 7:
            error_text = "A singular matrix was encountered during the\
                computation of the solution matrix X"
        e = ValueError(error_text)
        e.info = out[-1]
        raise e

    return out[:-1]

def sg03bd(n,m,A,E,Q,Z,B,dico,fact='N',trans='N',ldwork=None):
    """U,scale,alpha = sg03bd(dico,n,m,A,E,Q,Z,B,[fact,trans,ldwork])

     To compute the Cholesky factor U of the matrix X,

                 T
        X = op(U)  * op(U),

     which is the solution of either the generalized
     c-stable continuous-time Lyapunov equation

             T                    T
        op(A)  * X * op(E) + op(E)  * X * op(A)

                 2        T
        = - SCALE  * op(B)  * op(B),                                (1)

     or the generalized d-stable discrete-time Lyapunov equation

             T                    T
        op(A)  * X * op(A) - op(E)  * X * op(E)

                 2        T
        = - SCALE  * op(B)  * op(B),                                (2)

     without first finding X and without the need to form the matrix
     op(B)**T * op(B).

     op(K) is either K or K**T for K = A, B, E, U. A and E are N-by-N
     matrices, op(B) is an M-by-N matrix. The resulting matrix U is an
     N-by-N upper triangular matrix with non-negative entries on its
     main diagonal. SCALE is an output scale factor set to avoid
     overflow in U.

     In the continuous-time case (1) the pencil A - lambda * E must be
     c-stable (that is, all eigenvalues must have negative real parts).
     In the discrete-time case (2) the pencil A - lambda * E must be
     d-stable (that is, the moduli of all eigenvalues must be smaller
     than one).

     Required arguments
     __________________

         n : input int
             The order of the matrix A.  n >= 0.
         m : input int
             The number of rows in the matrix op(B).  m >= 0.
         A : input rank-2 array('d'), shape  (n,n)
             On entry, if fact = 'F', then the leading n-by-n upper
             Hessenberg part of this array must contain the
             generalized Schur factor A_s of the matrix A (see
             definition (3) in section METHOD). A_s must be an upper
             quasitriangular matrix. The elements below the upper
             Hessenberg part of the array A are not referenced.
             If fact = 'N', then the leading n-by-n part of this
             array must contain the matrix A.
             On exit, the leading n-by-n part of this array contains
             the generalized Schur factor A_s of the matrix A. (A_s is
             an upper quasitriangular matrix.)
         E : input rank-2 array('d'), shape  (n,n)
             On entry, if fact = 'F', then the leading n-by-n upper
             triangular part of this array must contain the
             generalized Schur factor E_s of the matrix E (see
             definition (4) in section METHOD). The elements below the
             upper triangular part of the array E are not referenced.
             If fact = 'N', then the leading n-by-n part of this
             array must contain the coefficient matrix E of the
             equation.
             On exit, the leading n-by-n part of this array contains
             the generalized Schur factor E_s of the matrix E. (E_s is
             an upper triangular matrix.)
         Q : input rank-2 array('d'), shape  (n,n)
             On entry, if fact = 'F', then the leading n-by-n part of
             this array must contain the orthogonal matrix Q from
             the generalized Schur factorization (see definitions (3)
             and (4) in section METHOD).
             If fact = 'N', Q need not be set on entry.
             On exit, the leading n-by-n part of this array contains
             the orthogonal matrix Q from the generalized Schur
             factorization.
         Z : input rank-2 array('d'), shape  (n,n)
             On entry, if fact = 'F', then the leading n-by-n part of
             this array must contain the orthogonal matrix Z from
             the generalized Schur factorization (see definitions (3)
             and (4) in section METHOD).
             If fact = 'N', Z need not be set on entry.
             On exit, the leading n-by-n part of this array contains
             the orthogonal matrix Z from the generalized Schur
             factorization.
         B : input rank-2 array('d'), shape  (n,n1)
             On entry, if trans = 'T', the leading n-by-m part of this
             array must contain the matrix B and n1 >= max(m,n).
             If trans = 'N', the leading m-by-n part of this array
             must contain the matrix B and n1 >= n.
             On exit, the leading n-by-n part of this array contains
             the Cholesky factor U of the solution matrix X of the
             problem, X = op(U)**T * op(U).
             If m = 0 and n > 0, then U is set to zero.
         dico : input string(len=1)
             Specifies which type of the equation is considered:
             = 'C':  Continuous-time equation (1);
             = 'D':  Discrete-time equation (2).

     Optional arguments
     __________________

         fact := 'N' input string(len=1)
             Specifies whether the generalized real Schur
             factorization of the pencil A - lambda * E is supplied
             on entry or not:
             = 'N':  Factorization is not supplied;
             = 'F':  Factorization is supplied.
         trans := 'N' input string(len=1)
             Specifies whether the transposed equation is to be solved
             or not:
             = 'N':  op(A) = A,    op(E) = E;
             = 'T':  op(A) = A**T, op(E) = E**T.
         ldwork := None input int
             The dimension of the array dwork.
             ldwork >= max(1,4*n,6*n-6),  if fact = 'N';
             ldwork >= max(1,2*n,6*n-6),  if fact = 'F'.
             For good performance, ldwork should be larger.

     Return objects
     ______________
    
         U : rank-2 array('d'), shape  (n,n)
             The leading n-by-b part of this array contains
             the Cholesky factor U of the solution matrix X of the
             problem, X = op(U)**T * op(U).
             If m = 0 and m > 0, then U is set to zero.
         scale : float
             The scale factor set to avoid overflow in U.
             0 < scale <= 1.
         alpha : rank-1 array('c'), shape (n)
             If INFO = 0, 3, 5, 6, or 7, then
             (alpha(j), j=1,...,n, are the
             eigenvalues of the matrix pencil A - lambda * E.

     Raises
     ______

         ValueError : e
             e.info contains information about the exact type of exception
              = 0:  successful exit;
              < 0:  if info = -i, the i-th argument had an illegal
                 value;
              = 1:  the pencil A - lambda * E is (nearly) singular;
                 perturbed values were used to solve the equation
                 (but the reduced (quasi)triangular matrices A and E
                 are unchanged);
              = 2:  fact = 'F' and the matrix contained in the upper
                 Hessenberg part of the array A is not in upper
                 quasitriangular form;
              = 3:  fact = 'F' and there is a 2-by-2 block on the main
                 diagonal of the pencil A_s - lambda * E_s whose
                 eigenvalues are not conjugate complex;
              = 4:  fact = 'N' and the pencil A - lambda * E cannot be
                 reduced to generalized Schur form: LAPACK routine
                 DGEGS (or DGGES) has failed to converge;
              = 5:  dico = 'C' and the pencil A - lambda * E is not
                 c-stable;
              = 6:  dico = 'D' and the pencil A - lambda * E is not
                 d-stable;
              = 7:  the LAPACK routine DSYEVX utilized to factorize M3
                 failed to converge in the discrete-time case (see
                 section METHOD for SLICOT Library routine SG03BU).
                 This error is unlikely to occur.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'fact', 'trans', 'n', 'm', 'A', 'LDA'+hidden, 'E',
        'LDE'+hidden, 'Q', 'LDQ'+hidden, 'Z', 'LDZ'+hidden, 'B', 'LDB'+hidden,
        'scale', 'alphar'+hidden, 'alphai'+hidden, 'IWORK'+hidden, 'DWORK'+hidden,
        'ldwork', 'INFO'+hidden]
    if ldwork is None:
        ldwork = max(1,4*n,6*n-6)
    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be either D or C')
    out = _wrapper.sg03bd(dico,n,m,A,E,Q,Z,B,fact=fact,trans=trans,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        error_text = """the pencil A - lambda * E is (nearly) singular;
                 perturbed values were used to solve the equation
                 (but the reduced (quasi)triangular matrices A and E
                 are unchanged)."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        error_text = """the matrix contained in the upper
                 Hessenberg part of the array A is not in upper
                 quasitriangular form."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        error_text = """there is a 2-by-2 block on the main
                 diagonal of the pencil A_s - lambda * E_s whose
                 eigenvalues are not conjugate complex."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 4:
        error_text = """the pencil A - lambda * E cannot be
                 reduced to generalized Schur form: LAPACK routine
                 DGEGS (or DGGES) has failed to converge."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 5:
        error_text = """the pencil A - lambda * E is not
                 c-stable."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 6:
        error_text = """the pencil A - lambda * E is not
                 d-stable."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 7:
        error_text = """the LAPACK routine DSYEVX utilized to factorize M3
                 failed to converge in the discrete-time case (see
                 section METHOD for SLICOT Library routine SG03BU).
                 This error is unlikely to occur."""
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    U,scale,alphar,alphai,beta = out[:-1]
    alpha = _np.zeros(n,'complex64')
    alpha.real = alphar[0:n]
    alpha.imag = alphai[0:n]
    return U,scale,alpha/beta
