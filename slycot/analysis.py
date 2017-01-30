#!/usr/bin/env python
#
#       analysis.py
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
import warnings

def ab01nd(n,m,A,B,jobz='N',tol=0,ldwork=None):
    """ Ac,Bc,ncont,indcon,nblk,Z,tau = ab01nd_i(n,m,A,B,[jobz,tol,ldwork])

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
            elementary reflectors used in the reduction of B and A."""

    hidden = ' (hidden by the wrapper)'
    arg_list = ['jobz', 'n', 'm', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden,
    'ncont', 'indcon', 'nblk', 'Z', 'LDZ'+hidden, 'tau', 'tol',
    'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 'info'+hidden]
    if ldwork is None:
        ldwork = max(n,3*m)
    if jobz == 'N':
        out = _wrapper.ab01nd_n(n,m,A,B,tol=tol,ldwork=ldwork)
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
            raise ValueError(error_text)
        # sets Z to None
        out[5] = None
        return out[:-1]
    if jobz == 'I':
        out = _wrapper.ab01nd_i(n,m,A,B,tol=tol,ldwork=ldwork)
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
            e = ValueError(error_text)
            e.info = out[-1]
            raise e
        return out[:-1]
    if jobz == 'F':
        out = _wrapper.ab01nd_f(n,m,A,B,tol=tol,ldwork=ldwork)
        if out[-1] < 0:
            error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
            e = ValueError(error_text)
            e.info = out[-1]
            raise e
        return out[:-1]
    raise ValueError('jobz must be either N, I or F')

def ab05md(n1,m1,p1,n2,p2,A1,B1,C1,D1,A2,B2,C2,D2,uplo='U'):
    """ n,a,b,c,d = ab05md(n1,m1,p1,n2,p2,a1,b1,c1,d1,a2,b2,c2,d2,[uplo])

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

    Notes:
        The implemented methods rely on accuracy enhancing square-root or
        balancing-free square-root techniques.
        The algorithms require less than 30N^3  floating point operations.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['uplo', 'OVER'+hidden, 'n1', 'm1', 'p1', 'n2', 'p2', 'A1',
        'LDA1'+hidden, 'B1', 'LDB1'+hidden, 'C1', 'LDC1'+hidden, 'D1',
        'LDD1'+hidden, 'A2', 'LDA2'+hidden, 'B2', 'LDB2'+hidden, 'C2',
        'LDC2'+hidden, 'D2', 'LDD2'+hidden, 'n', 'A', 'LDA'+hidden, 'B',
        'LDB'+hidden, 'C', 'LDC'+hidden, 'D', 'LDD'+hidden, 'DWORK'+hidden,
        'ldwork', 'info'+hidden ]
    out = _wrapper.ab05md(n1,m1,p1,n2,p2,A1,B1,C1,D1,A2,B2,C2,D2,uplo=uplo)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    return out[:-1]

def ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,alpha=1.0,ldwork=None):
    """  n,A,B,C,D = ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,[alpha,ldwork])

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
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['over'+hidden, 'n1', 'm1', 'p1', 'n2', 'alpha', 'A1', 'LDA1'+hidden,
        'B1', 'LDB1'+hidden, 'C1', 'LDC1'+hidden, 'D1', 'LDD1'+hidden, 'A2',
        'LDA2'+hidden, 'B2', 'LDB2'+hidden, 'C2', 'LDC2'+hidden, 'D2',
        'LDD2'+hidden, 'n', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'C',
        'LDC'+hidden, 'D', 'LDD'+hidden, 'IWORK'+hidden, 'DWORK'+hidden,
        'ldwork', 'info'+hidden]
    if ldwork is None:
        ldwork = max(p1*p1,m1*m1,n1*p1)
    out = _wrapper.ab05nd(n1,m1,p1,n2,alpha,A1,B1,C1,D1,A2,B2,C2,D2,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] > 0:
        e = ArithmeticError('The resulting system is not completely controllable.')
        e.info = out[-1]
        raise e
    return out[:-1]

def ab07nd(n,m,A,B,C,D,ldwork=None):
    """ A_i,B_i,C_i,D_i,rcond = ab07nd(n,m,A,B,C,D,[ldwork])

    To compute the inverse (A_i,B_i,C_i,D_i) of a given system (A,B,C,D).

    Required arguments:
        n : input int
            The order of the state matrix A.  n >= 0.
        m : input int
            The number of system inputs and outputs.  m >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state matrix
            A of the original system.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input matrix
            B of the original system.
        C : input rank-2 array('d') with bounds (m,n)
            The leading m-by-n part of this array must contain the output matrix
            C of the original system.
        D : input rank-2 array('d') with bounds (m,m)
            The leading m-by-m part of this array must contain the feedthrough
            matrix D of the original system.
    Optional arguments:
        ldwork := None input int
            The length of the cache array. The default value is max(1,4*m),
            for better performance should be larger.
    Return objects:
        A_i : rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array contains the state matrix A_i
            of the inverse system.
        B_i : rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array contains the input matrix B_i
            of the inverse system.
        C_i : rank-2 array('d') with bounds (m,n)
            The leading m-by-n part of this array contains the output matrix C_i
            of the inverse system.
        D_i : rank-2 array('d') with bounds (m,m)
            The leading m-by-m part of this array contains the feedthrough
            matrix D_i of the inverse system.
        rcond : float
            The estimated reciprocal condition number of the feedthrough matrix
            D of the original system.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'C',
    'LDC'+hidden, 'D', 'LDD'+hidden, 'rcond', 'IWORK'+hidden, 'DWORK'+hidden,
    'ldwork', 'INFO'+hidden]
    if ldwork is None:
        ldwork = max(1,4*m)
    out = _wrapper.ab07nd(n,m,A,B,C,D,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == m+1:
        e = ValueError('Entry matrix D is numerically singular.')
        e.info = out[-1]
        raise e
    if out[-1] > 0:
        e = ValueError('Entry matrix D is exactly singular, the (%i,%i) diagonal element is zero.' %(out[-1],out[-1]))
        e.info = out[-1]
        raise e
    return out[:-1]

def ab08nd(n,m,p,A,B,C,D,equil='N',tol=0,ldwork=None):
    """ nu,rank,dinfz,nkror,nkrol,infz,kronr,kronl,Af,Bf = ab08nd(n,m,p,A,B,C,D,[equil,tol,ldwork])

    To construct for a linear multivariable system described by a state-space
    model (A,B,C,D) a regular pencil (Af - lambda*Bf ) which has the invariant
    zeros of the system as generalized eigenvalues.
    The routine also computes the orders of the infinite zeros and the
    right and left Kronecker indices of the system (A,B,C,D).

    Required arguments:
        n : input int
            The number of state variables.  n >= 0.
        m : input int
            The number of system inputs.  m >= 0.
        p : input int
            The number of system outputs.  p >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state
            dynamics matrix A of the system.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input/state
            matrix B of the system.
        C : input rank-2 array('d') with bounds (p,n)
            The leading p-by-n part of this array must contain the state/output
            matrix C of the system.
        D : input rank-2 array('d') with bounds (p,m)
            The leading p-by-m part of this array must contain the direct
            transmission matrix D of the system.
    Optional arguments:
        equil := 'N' input string(len=1)
            Specifies whether the user wishes to balance the compound matrix
            as follows:
            = 'S':  Perform balancing (scaling);
            = 'N':  Do not perform balancing.
        tol := 0.0 input float
            A tolerance used in rank decisions to determine the effective rank,
            which is defined as the order of the largest leading (or trailing)
            triangular submatrix in the QR (or RQ) factorization with column
            (or row) pivoting whose estimated condition number is less than 1/tol.
        ldwork := None input int
            The length of the cache array. The default value is n + 3*max(m,p),
            for better performance should be larger.
    Return objects:
        nu : int
            The number of (finite) invariant zeros.
        rank : int
            The normal rank of the transfer function matrix.
        dinfz : int
            The maximum degree of infinite elementary divisors.
        nkror : int
            The number of right Kronecker indices.
        nkrol : int
            The number of left Kronecker indices.
        infz : rank-1 array('i') with bounds (n)
            The leading dinfz elements of infz contain information on the
            infinite elementary divisors as follows: the system has infz(i)
            infinite elementary divisors of degree i, where i = 1,2,...,dinfz.
        kronr : rank-1 array('i') with bounds (max(n,m)+1)
            the leading nkror elements of this array contain the right kronecker
            (column) indices.
        kronl : rank-1 array('i') with bounds (max(n,p)+1)
            the leading nkrol elements of this array contain the left kronecker
            (row) indices.
        Af : rank-2 array('d') with bounds (max(1,n+m),n+min(p,m))
            the leading nu-by-nu part of this array contains the coefficient
            matrix Af of the reduced pencil. the remainder of the leading
            (n+m)-by-(n+min(p,m)) part is used as internal workspace.
        Bf : rank-2 array('d') with bounds (max(1,n+p),n+m)
            The leading nu-by-nu part of this array contains the coefficient
            matrix Bf of the reduced pencil. the remainder of the leading
            (n+p)-by-(n+m) part is used as internal workspace.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['equil', 'n', 'm', 'p', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden,
        'C', 'LDC'+hidden, 'D', 'LDD'+hidden, 'nu', 'rank', 'dinfz', 'nkror',
        'nkrol', 'infz', 'kronr', 'kronl', 'Af', 'LDAF'+hidden, 'Bf',
        'LDBF'+hidden, 'tol', 'IWORK'+hidden, 'DWORK'+hidden, 'ldwork',
        'INFO'+hidden]
    if ldwork is None:
        ldwork = n+3*max(m,p) #only an upper bound
    out = _wrapper.ab08nd(n,m,p,A,B,C,D,equil=equil,tol=tol,ldwork=ldwork)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    return out[:-1]

def ab09ad(dico,job,equil,n,m,p,A,B,C,nr=None,tol=0,ldwork=None):
    """ nr,Ar,Br,Cr,hsv = ab09ad(dico,job,equil,n,m,p,A,B,C,[nr,tol,ldwork])

    Compute reduced order State-Space-Model (Ar, Br, Cr) for a stable system
    (A, B, C) by using either the square-root or the balancing-free square-
    root Balance & truncate (B & T) model reduction method.

    Required arguments:
        dico : {'D', 'C'} input string(len=1)
            Indicate whether the system is discrete `D` or continuous `C`
        job : {'B', 'N'} input string(len=1)
            Balance `B` or not `N`
        equil : {'S', 'N'} input string(len=1)
            Scale `S` or not `N`
        n : input int
            The number of state variables.  n >= 0.
        m : input int
            The number of system inputs.  m >= 0.
        p : input int
            The number of system outputs.  p >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state
            dynamics matrix A of the system.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input/state
            matrix B of the system.
        C : input rank-2 array('d') with bounds (p,n)
            The leading p-by-n part of this array must contain the
            state/output matrix C of the system.

    Optional arguments:
        nr := None input int
            `nr` is the desired order of the resulting reduced order
            system.  ``0 <= nr <= n``. Automatically determined by `tol` if
            ``nr is None`` and returned. See return object `nr`.
        tol := 0 input double precision
            If ``nr is None``, `tol`contains the tolerance for determining the
            order of the reduced system. For model reduction, th recommended
            value is ``tol = c * HNORM(A, B, C)``, where `c` is a constan in the
            interval ``[0.00001, 0.001]`` and ``HNORM(A, B, C)`` is the
            Hankel-Norm of the given sysstem (computed in ``HSV(1)``). For
            computing a minimal realization, the recommended value is
            ``tol = n * eps * HNORM(A, B, C)``, where `eps` is the machine
            precision (see LAPACK Library Routine `DLAMCH`). This value is
            used by default if ``tol <= 0`` on entry. If `nr` is specified,
            the value of `tol` is ignored.
        ldwork := None input int
            The length of the cache array. The default value is
            ``n*(2*n+max(n,m,p)+5) + n*(n+1)/2 ~= 3.5*n**2 + 5*n``,
            a larger value should lead to better performance.

    Return objects :
        nr : output int
            `nr` is the order of the resulting reduced order model.
            `nr` is set as follows:
            If on input ``nr is not None``, `nr` is equal to ``MIN(nr,NMIN)``,
            where `nr` is the desired order on entry and `NMIN` is the order
            of a minimal realization of the given system; `NMIN` is
            determined as the number of Hankel singular values greater
            than ``n*eps*HNORM(A,B,C)``, where `eps` is the machine
            precision (see LAPACK Library Routine DLAMCH) and
            ``HNORM(A,B,C)`` is the Hankel norm of the system (computed
            in ``HSV(1)``);
            If on input ``nr is None``, `nr` is equal to the number of Hankel
            singular values greater than ``MAX(tol,n*eps*HNORM(A,B,C))``.
        Ar : rank-2 array('d') with bounds ``(nr,nr)``
            This array contains the state dynamics matrix `Ar` of the reduced
            order system.
        Br : rank-2 array('d') with bounds ``(nr,m)``
            Tthis array contains the input/state matrix `Br` of the reduced
            order system.
        Cr : rank-2 array('d') with bounds ``(p,nr)``
            This array contains the state/output matrix `Cr` of the reduced
            order system.
        hsv : output double precision array, dimension ``(n)``
            If ``INFO = 0``, it contains the Hankel singular values of
            the original system ordered decreasingly. ``HSV(1)`` is the
            Hankel norm of the system.

    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'A',
        'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'hsv', 'tol',
        'iwork'+hidden, 'dwork'+hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1,n*(2*n+max(n,max(m,p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0 #order will be computed by the routine
    else:
        ordsel = 'F'
    if dico != 'C' and dico != 'D':
        raise ValueError('Parameter dico had an illegal value')
    if job != 'B' and job != 'N':
        raise ValueError('Parameter job had an illegal value')
    if equil != 'S' and equil != 'N':
        raise ValueError('Parameter equil had an illegal value')
    out = _wrapper.ab09ad(dico,job,equil,ordsel,n,m,p,nr,A,B,C,tol,ldwork)
    if out[-2] == 1:
        warnings.warn("The selected order nr is greater\
                than the order of a minimal realization of the\
                given system. It was set automatically to a value\
                corresponding to the order of a minimal realization\
                of the system")
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        e = ArithmeticError('The reduction of A to the real Schur form failed')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ArithmeticError('The state matrix A is not stable or not convergent')
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        e = ArithmeticError('The computation of Hankel singular values failed')
        e.info = out[-1]
        raise e
    Nr,A,B,C,hsv = out[:-2]
    return Nr, A[:Nr,:Nr], B[:Nr,:], C[:,:Nr], hsv

def ab09ax(dico,job,n,m,p,A,B,C,nr=None,tol=0.0,ldwork=None):
    """``nr,Ar,Br,Cr,hsv,T,Ti = ab09ad(dico,job,equil,n,m,p,nr,A,B,C,[nr,tol,ldwork])``

    To compute a reduced order model ``(Ar,Br,Cr)`` for a stable original
    state-space representation ``(A,B,C)`` by using either the square-root
    or the balancing-free square-root Balance & Truncate model
    reduction method. The state dynamics matrix `A` of the original
    system is an upper quasi-triangular matrix in *real Schur canonical
    form.* The matrices of the reduced order system are computed using
    the truncation formulas:

        ``Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T`` .

    Required arguments :
        dico : {'D', 'C'} input string(len=1)
            Indicate whether the system is discrete `D` or continuous `C`
        job : {'B', 'N'} input string(len=1)
            Balance `B` or not `N`
        n : input int
            The number of state variables.  n >= 0.
        m : input int
            The number of system inputs.  m >= 0.
        p : input int
            The number of system outputs.  p >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state
            dynamics matrix A of the system *in real Schur form.*
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input/state
            matrix B of the system.
        C : input rank-2 array('d') with bounds (p,n)
            The leading p-by-n part of this array must contain the
            state/output matrix C of the system.

    Optional arguments:
        nr := None input int
            `nr` is the desired order of the resulting reduced order
            system.  ``0 <= nr <= n``. Automatically determined by `tol` if
            ``nr is None`` and returned. See return object `nr`.
        tol := 0 input double precision
            If ``nr is None``, `tol`contains the tolerance for determining the
            order of the reduced system. For model reduction, the recommended
            value is ``tol = c * HNORM(A, B, C)``, where `c` is a constant in
            the interval ``[0.00001, 0.001]`` and ``HNORM(A, B, C)`` is
            the Hankel-Norm of the given sysstem (computed in ``HSV(1)``). For
            computing a minimal realization, the recommended value is
            ``tol = n * eps * HNORM(A, B, C)``, where `eps` is the machine
            precision (see LAPACK Library Routine `DLAMCH`). This value is
            used by default if ``tol <= 0`` on entry. If `nr` is specified,
            the value of `tol` is ignored.
        ldwork := None input int
            The length of the cache array. The default value is
            ``n*(2*n+max(n,m,p)+5) + n*(n+1)/2 ~= 3.5*n**2 + 5*n``,
            a larger value should lead to better performance.

    Return objects :
        nr : output int
            `nr` is the order of the resulting reduced order model.
            `nr` is set as follows:
            If on input ``nr is not None``, `nr` is equal to ``MIN(nr,NMIN)``,
            where `nr` is the desired order on entry and `NMIN` is the order
            of a minimal realization of the given system; `NMIN` is
            determined as the number of Hankel singular values greater
            than ``n*eps*HNORM(A,B,C)``, where `eps` is the machine
            precision (see LAPACK Library Routine DLAMCH) and
            ``HNORM(A,B,C)`` is the Hankel norm of the system (computed
            in ``HSV(1)``);
            If on input ``nr is None``, `nr` is equal to the number of Hankel
            singular values greater than ``MAX(tol,n*eps*HNORM(A,B,C))``.
        Ar : rank-2 array('d') with bounds ``(nr,nr)``
            This array contains the state dynamics matrix `Ar` of the reduced
            order system.
        Br : rank-2 array('d') with bounds ``(nr,m)``
            Tthis array contains the input/state matrix `Br` of the reduced
            order system.
        Cr : rank-2 array('d') with bounds ``(p,nr)``
            This array contains the state/output matrix `Cr` of the reduced
            order system.
        hsv : output double precision array, dimension ``(n)``
            If ``INFO = 0``, it contains the Hankel singular values of
            the original system ordered decreasingly. ``HSV(1)`` is the
            Hankel norm of the system.
        T : rank-2 array('d') with bounds ``(n,nr)``
            This array contains the right truncation matrix `T` of the reduced
            order system.
        Ti : rank-2 array('d') with bounds ``(nr,n)``
            This array contains the left truncation matrix `Ti` of the reduced
            order system.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'ordsel', 'n', 'm', 'p', 'nr', 'A',
        'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'hsv', 'T',
        'ldt'+hidden, 'Ti', 'ldti'+hidden, 'tol', 'iwork'+hidden,
        'dwork'+hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1,n*(2*n+max(n,max(m,p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0 #order will be computed by the routine
    else:
        ordsel = 'F'
    if dico != 'C' and dico != 'D':
        raise ValueError('Parameter dico had an illegal value')
    if job != 'B' and job != 'N':
        raise ValueError('Parameter job had an illegal value')
    out = _wrapper.ab09ax(dico,job,ordsel,n,m,p,nr,A,B,C,tol,ldwork)
    if out[-2] == 1:
        warnings.warn("The selected order nr is greater\
                than the order of a minimal realization of the\
                given system. It was set automatically to a value\
                corresponding to the order of a minimal realization\
                of the system")
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        e = ArithmeticError('The state matrix A is not stable or not convergent')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ArithmeticError('The computation of Hankel singular values failed')
        e.info = out[-1]
        raise e
    nr,A,B,C,hsv,T,Ti = out[:-2]
    return nr, A[:nr,:nr], B[:nr,:], C[:,:nr], hsv, T[:,:nr], Ti[:nr,:]

def ab09bd(dico,job,equil,n,m,p,A,B,C,D,nr=None,tol1=0,tol2=0,ldwork=None):
    """ nr,Ar,Br,Cr,Dr,hsv = ab09bd(dico,job,equil,n,m,p,A,B,C,D,[nr,tol1,tol2,ldwork])
    
    To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
    original state-space representation (A,B,C,D) by using either the
    square-root or the balancing-free square-root Singular
    Perturbation Approximation (SPA) model reduction method.
    Must supply either nr or tolerance values.

    Arguments
        Mode Parameters
            dico
                Specifies the type of the original system as follows:
                = 'C':  continuous-time system;
                = 'D':  discrete-time system.
            job
                Specifies the model reduction approach to be used
                as follows:
                = 'B':  use the square-root SPA method;
                = 'N':  use the balancing-free square-root SPA method.
            equil
                Specifies whether the user wishes to preliminarily
                equilibrate the triplet (A,B,C) as follows:
                = 'S':  perform equilibration (scaling);
                = 'N':  do not perform equilibration.

        Required arguments
            n : input int
                The order of the original state-space representation, i.e.
                the order of the matrix A.  n >= 0.
            m : input int
                The number of system inputs.  m >= 0.
            p : input int
                The number of system outputs.  p >= 0.
            A : input rank-2 array('d') with bounds (n,n)
                On entry, the leading n-by-n part of this array must
                contain the state dynamics matrix A.
            B : input rank-2 array('d') with bounds (n,m)
                On entry, the leading n-by-m part of this array must
                contain the original input/state matrix B.
            C : input rank-2 array('d') with bounds (p,n)
                On entry, the leading p-by-n part of this array must
                contain the original state/output matrix C.
            D : input rank-2 array('d') with bounds (p,m)
                On entry, the leading p-by-m part of this array must
                contain the original input/output matrix D.

        Optional arguments :
            nr :=None input int
                nr is the desired order of
                the resulting reduced order system.  0 <= nr <= n.
            tol1 :=0 input double precision
                If ordsel = 'A', tol1 contains the tolerance for
                determining the order of reduced system.
                For model reduction, the recommended value is
                tol1 = c*hnorm(A,B,C), where c is a constant in the
                interval [0.00001,0.001], and hnorm(A,B,C) is the
                Hankel-norm of the given system (computed in hsv(1)).
                For computing a minimal realization, the recommended
                value is tol1 = n*eps*hnorm(A,B,C), where eps is the
                machine precision (see LAPACK Library Routine DLAMCH).
                This value is used by default if tol1 <= 0 on entry.
                If ordsel = 'F', the value of tol1 is ignored.
            tol2 :=0 input double precision
                The tolerance for determining the order of a minimal
                realization of the given system. The recommended value is
                tol2 = n*eps*hnorm(A,B,C). This value is used by default
                if tol2 <= 0 on entry.
                If tol2 > 0, then tol2 <= tol1.
            ldwork := None input int
                The length of the cache array. The default value is n + 3*max(m,p),
                for better performance should be larger.

        Return objects
            nr : output int
                nr is the order of the resulting reduced order model.
                nr is set as follows:
                if ordsel = 'F', nr is equal to min(nr,nmin), where nr
                is the desired order on entry and nmin is the order of a
                minimal realization of the given system; nmin is
                determined as the number of Hankel singular values greater
                than n*eps*hnorm(A,B,C), where eps is the machine
                precision (see LAPACK Library Routine DLAMCH) and
                hnorm(A,B,C) is the Hankel norm of the system (computed
                in hsv(1));
                if ordsel = 'A', nr is equal to the number of Hankel
                singular values greater than max(tol1,n*eps*hnorm(A,B,C)).
            Ar : rank-2 array('d') with bounds (nr,nr)
                the leading nr-by-nr part of this array contains the
                state dynamics matrix Ar of the reduced order system.
            Br : rank-2 array('d') with bounds (nr,m)
                the leading nr-by-m part of this array contains the
                input/state matrix Br of the reduced order system.
            Cr : rank-2 array('d') with bounds (p,nr)
                the leading p-by-nr part of this array contains the
                state/output matrix Cr of the reduced order system.
            Dr : rank-2 array('d') with bounds (p,m)
                the leading p-by-m part of this array contains the
                input/output matrix Dr of the reduced order system.
            hsv : output double precision array, dimension (n)
                If info = 0, it contains the Hankel singular values of
                the original system ordered decreasingly. hsv(1) is the
                Hankel norm of the system.
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'A',
        'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'D', 'ldd'+hidden, 'hsv', 'tol1', 'tol2',
        'iwork'+hidden, 'dwork'+hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1,n*(2*n+max(n,max(m,p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0 #order will be computed by the routine
    else:
        ordsel = 'F'
    if dico != 'C' and dico != 'D':
        raise ValueError('Parameter dico had an illegal value')
    if job != 'B' and job != 'N':
        raise ValueError('Parameter job had an illegal value')
    if equil != 'S' and equil != 'N':
        raise ValueError('Parameter equil had an illegal value')
    out = _wrapper.ab09bd(dico,job,equil,ordsel,n,m,p,nr,A,B,C,D,tol1,tol2,ldwork)
    if out[-2] == 1:
        warnings.warn("The selected order nr is greater\
                than the order of a minimal realization of the\
                given system. It was set automatically to a value\
                corresponding to the order of a minimal realization\
                of the system")
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        e = ArithmeticError('The reduction of A to the real Schur form failed')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ArithmeticError('The state matrix A is not stable (if dico = C) or not convergent (if dico = D)')
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        e = ArithmeticError('The computation of Hankel singular values failed')
        e.info = out[-1]
        raise e
    Nr,A,B,C,D,hsv = out[:-2]
    return Nr, A[:Nr,:Nr], B[:Nr,:], C[:,:Nr],D[:,:], hsv

def ab09md(dico,job,equil,n,m,p,A,B,C,alpha=None,nr=None,tol=0,ldwork=None):
    """ nr,Ar,Br,Cr,ns,hsv = ab09md(dico,job,equil,n,m,p,A,B,C,[alpha,nr,tol,ldwork])

    To compute a reduced order model (Ar,Br,Cr) for an original
    state-space representation (A,B,C) by using either the square-root
    or the balancing-free square-root Balance & Truncate (B & T)
    model reduction method for the ALPHA-stable part of the system.

    Arguments
        Mode Parameters
            dico
                Specifies the type of the original system as follows:
                = 'C':  continuous-time system;
                = 'D':  discrete-time system.
            job
                Specifies the model reduction approach to be used
                as follows:
                = 'B':  use the square-root Balance & Truncate method;
                = 'N':  use the balancing-free square-root
                        Balance & Truncate method.
            equil
                Specifies whether the user wishes to preliminarily
                equilibrate the triplet (A,B,C) as follows:
                = 'S':  perform equilibration (scaling);
                = 'N':  do not perform equilibration.

        Required arguments
            n : input int
                The order of the original state-space representation, i.e.
                the order of the matrix A.  n >= 0.
            m : input int
                The number of system inputs.  m >= 0.
            p : input int
                The number of system outputs.  p >= 0.
            A : input rank-2 array('d'), dimension (n,n)
                On entry, the leading N-by-N part of this array must
                contain the state dynamics matrix A.                
            B : input rank-2 array('d'), dimension (n,m)
                On entry, the leading N-by-M part of this array must
                contain the original input/state matrix B.
            C : input rank-2 array('d'), dimension (p,n)
                On entry, the leading P-by-N part of this array must
                contain the original state/output matrix C.

         Optional arguments
            alpha :=None input double precision
                Specifies the alpha-stability boundary for the eigenvalues
                of the state dynamics matrix A. For a continuous-time
                system (dico = 'C'), alpha <= 0 is the boundary value for
                the real parts of eigenvalues, while for a discrete-time
                system (dico = 'D'), 0 <= alpha <= 1 represents the
                boundary value for the moduli of eigenvalues.
                The alpha-stability domain does not include the boundary.
            nr := None input int
                On entry with ordsel = 'F', nr is the desired order of the
                resulting reduced order system.  0 <= nr <= n.
            tol :=0 input double precision
                If ordsel = 'A', tol contains the tolerance for
                determining the order of reduced system.
                For model reduction, the recommended value is
                tol = c*hnorm(As,Bs,Cs), where c is a constant in the
                interval [0.00001,0.001], and hnorm(As,Bs,Cs) is the
                Hankel-norm of the alpha-stable part of the given system
                (computed in hsv(1)).
                If tol <= 0 on entry, the used default value is
                tol = ns*eps*hnorm(As,Bs,Cs), where ns is the number of
                alpha-stable eigenvalues of A and eps is the machine
                precision (see LAPACK Library Routine DLAMCH).
                This value is appropriate to compute a minimal realization
                of the alpha-stable part.
                If ordsel = 'F', the value of tol is ignored.
            ldwork :=None input int
                The length of the array dwork.
                ldwork >= max(1,n*(2*n+max(n,m,p)+5) + n*(n+1)/2).
                For optimum performance ldwork should be larger.

         Return objects
            nr : output int
                On exit, if info = 0, nr is the order of the resulting
                reduced order model. For a system with nu alpha-unstable
                eigenvalues and ns alpha-stable eigenvalues (nu+ns = n),
                nr is set as follows: if ordsel = 'F', nr is equal to
                nu+min(max(0,nr-nu),nmin), where nr is the desired order
                on entry, and nmin is the order of a minimal realization
                of the alpha-stable part of the given system; nmin is
                determined as the number of Hankel singular values greater
                than ns*eps*hnorm(As,Bs,Cs), where eps is the machine
                precision (see LAPACK Library Routine DLAMCH) and
                hnorm(As,Bs,Cs) is the Hankel norm of the alpha-stable
                part of the given system (computed in hsv(1));
                if ordsel = 'A', nr is the sum of nu and the number of
                Hankel singular values greater than
                max(tol,ns*eps*hnorm(As,Bs,Cs)).
            Ar : rank-2 array('d') with bounds (nr,nr)
                On exit, if info = 0, the leading nr-by-nr part of this
                array contains the state dynamics matrix Ar of the reduced
                order system.
                The resulting A has a block-diagonal form with two blocks.
                For a system with nu alpha-unstable eigenvalues and
                ns alpha-stable eigenvalues (nu+ns = n), the leading
                nu-by-nu block contains the unreduced part of A
                corresponding to alpha-unstable eigenvalues in an
                upper real Schur form.
                The trailing (nr+ns-n)-by-(nr+ns-n) block contains
                the reduced part of A corresponding to alpha-stable
                eigenvalues.
            Br : rank-2 array('d') with bounds (nr,m)
                On exit, if info = 0, the leading nr-by-m part of this
                array contains the input/state matrix Br of the reduced
                order system.
            Cr : rank-2 array('d') with bounds (p,nr)
                On exit, if info = 0, the leading p-by-nr part of this
                array contains the state/output matrix Cr of the reduced
                order system.
            ns : output int
                The dimension of the alpha-stable subsystem.
            hsv : output double precision array, dimension (n)
                If info = 0, the leading ns elements of hsv contain the
                Hankel singular values of the alpha-stable part of the
                original system ordered decreasingly.
                hsv(1) is the Hankel norm of the alpha-stable subsystem.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'alpha', 'A',
        'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'ns', 'hsv', 'tol1', 'tol2',
        'iwork'+hidden, 'dwork'+hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1,n*(2*n+max(n,max(m,p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0 #order will be computed by the routine
    else:
        ordsel = 'F'
    if dico != 'C' and dico != 'D':
        raise ValueError('Parameter dico had an illegal value')
    if alpha is None:
        if dico == 'C':
            alpha = 0.
        elif dico == 'D':
            alpha = 1.
    if job != 'B' and job != 'N':
        raise ValueError('Parameter job had an illegal value')
    if equil != 'S' and equil != 'N':
        raise ValueError('Parameter equil had an illegal value')
    out = _wrapper.ab09md(dico,job,equil,ordsel,n,m,p,nr,alpha,A,B,C,tol,ldwork)
    if out[-2] == 1:
        warnings.warn("with ordsel = 'F', the selected order nr is greater\
                than nsmin, the sum of the order of the\
                alpha-unstable part and the order of a minimal\
                realization of the alpha-stable part of the given\
                system. In this case, the resulting nr is set equal\
                to nsmin.")
    if out[-2] == 2:
        warnings.warn("with ordsel = 'F', the selected order nr is less\
                than the order of the alpha-unstable part of the\
                given system. In this case nr is set equal to the\
                order of the alpha-unstable part.")
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        e = ArithmeticError('The reduction of A to the real Schur form failed')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ArithmeticError('the separation of the alpha-stable/unstable diagonal\
                blocks failed because of very close eigenvalues;')
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        e = ArithmeticError('The computation of Hankel singular values failed')
        e.info = out[-1]
        raise e
    Nr,A,B,C,Ns,hsv = out[:-2]
    return Nr, A[:Nr,:Nr], B[:Nr,:], C[:,:Nr], Ns, hsv

def ab09nd(dico,job,equil,n,m,p,A,B,C,D,alpha=None,nr=None,tol1=0,tol2=0,ldwork=None):
    """ nr,Ar,Br,Cr,Dr,ns,hsv = ab09nd(dico,job,equil,n,m,p,A,B,C,D,[alpha,nr,tol1,tol2,ldwork])
    
    To compute a reduced order model (Ar,Br,Cr,Dr) for an original
    state-space representation (A,B,C,D) by using either the
    square-root or the balancing-free square-root Singular
    Perturbation Approximation (SPA) model reduction method for the
    alpha-stable part of the system.

    Arguments
        Mode Parameters
            dico
                Specifies the type of the original system as follows:
                = 'C':  continuous-time system;
                = 'D':  discrete-time system.
            job
                Specifies the model reduction approach to be used
                as follows:
                = 'B':  use the square-root SPA method;
                = 'N':  use the balancing-free square-root SPA method.
            equil
                Specifies whether the user wishes to preliminarily
                equilibrate the triplet (A,B,C) as follows:
                = 'S':  perform equilibration (scaling);
                = 'N':  do not perform equilibration.

        Required arguments
            n : input int
                The order of the original state-space representation, i.e.
                the order of the matrix A.  n >= 0.
            m : input int
                The number of system inputs.  m >= 0.
            p : input int
                The number of system outputs.  p >= 0.
            A : input rank-2 array('d') with bounds (n,n)
                On entry, the leading n-by-n part of this array must
                contain the state dynamics matrix A.
            B : input rank-2 array('d') with bounds (n,m)
                On entry, the leading n-by-m part of this array must
                contain the original input/state matrix B.
            C : input rank-2 array('d') with bounds (p,n)
                On entry, the leading p-by-n part of this array must
                contain the original state/output matrix C.
            D : input rank-2 array('d') with bounds (p,m)
                On entry, the leading p-by-m part of this array must
                contain the original input/output matrix D.

        Optional arguments
            alpha :=None input double precision
                Specifies the alpha-stability boundary for the eigenvalues
                of the state dynamics matrix A. For a continuous-time
                system (dico = 'C'), alpha <= 0 is the boundary value for
                the real parts of eigenvalues, while for a discrete-time
                system (dico = 'D'), 0 <= alpha <= 1 represents the
                boundary value for the moduli of eigenvalues.
                The alpha-stability domain does not include the boundary.
            nr :=None input int
                nr is the desired order of
                the resulting reduced order system.  0 <= nr <= n.
            tol1 :=0 input double precision
                If ordsel = 'A', tol1 contains the tolerance for
                determining the order of reduced system.
                For model reduction, the recommended value is
                tol1 = c*hnorm(As,Bs,Cs), where c is a constant in the
                interval [0.00001,0.001], and hnorm(As,Bs,Cs) is the
                Hankel-norm of the alpha-stable part of the given system
                (computed in hsv(1)).
                If tol1 <= 0 on entry, the used default value is
                tol1 = ns*eps*hnorm(As,Bs,Cs), where NS is the number of
                alpha-stable eigenvalues of A and eps is the machine
                precision (see LAPACK Library Routine DLAMCH).
                This value is appropriate to compute a minimal realization
                of the alpha-stable part.
                If ordsel = 'F', the value of tol1 is ignored.
            tol2 :=0 input double precision
                The tolerance for determining the order of a minimal
                realization of the alpha-stable part of the given system.
                The recommended value is tol2 = ns*eps*hnorm(As,Bs,Cs).
                This value is used by default if tol2 <= 0 on entry.
                If tol2 > 0, then tol2 <= tol1.
            ldwork := None input int
                The length of the array dwork.
                ldwork >= max(1,n*(2*n+max(n,m,p)+5) + n*(n+1)/2).
                For optimum performance ldwork should be larger.

        Return objects
            nr : output int
                nr is the order of the resulting reduced order model.
                nr is set as follows:
                if ordsel = 'F', nr is equal to min(nr,nmin), where nr
                is the desired order on entry and nmin is the order of a
                minimal realization of the given system; nmin is
                determined as the number of Hankel singular values greater
                than n*eps*hnorm(A,B,C), where eps is the machine
                precision (see LAPACK Library Routine DLAMCH) and
                hnorm(A,B,C) is the Hankel norm of the system (computed
                in hsv(1));
                if ordsel = 'A', nr is equal to the number of Hankel
                singular values greater than max(TOL1,n*eps*hnorm(A,B,C)).
            Ar : rank-2 array('d') with bounds (nr,nr)
                the leading nr-by-nr part of this array contains the
                state dynamics matrix Ar of the reduced order system.
            Br : rank-2 array('d') with bounds (nr,m)
                the leading nr-by-m part of this array contains the
                input/state matrix Br of the reduced order system.
            Cr : rank-2 array('d') with bounds (p,nr)
                the leading p-by-nr part of this array contains the
                state/output matrix Cr of the reduced order system.
            Dr : rank-2 array('d') with bounds (p,m)
                the leading p-by-m part of this array contains the
                input/output matrix Dr of the reduced order system.
            ns : output int
                The dimension of the alpha-stable subsystem.
            hsv : output double precision array, dimension (n)
                If info = 0, it contains the Hankel singular values of
                the original system ordered decreasingly. hsv(1) is the
                Hankel norm of the system.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'alpha', 'A',
        'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'D', 'ldd'+hidden, 'ns', 'hsv', 'tol1', 'tol2',
        'iwork'+hidden, 'dwork'+hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1,n*(2*n+max(n,max(m,p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0 #order will be computed by the routine
    else:
        ordsel = 'F'
    if dico != 'C' and dico != 'D':
        raise ValueError('Parameter dico had an illegal value')
    if alpha is None:
        if dico == 'C':
            alpha = 0.
        elif dico == 'D':
            alpha = 1.
    if job != 'B' and job != 'N':
        raise ValueError('Parameter job had an illegal value')
    if equil != 'S' and equil != 'N':
        raise ValueError('Parameter equil had an illegal value')
    out = _wrapper.ab09nd(dico,job,equil,ordsel,n,m,p,nr,alpha,A,B,C,D,tol1,tol2,ldwork)
    if out[-2] == 1:
        warnings.warn("with ordsel = 'F', the selected order nr is greater\
                than nsmin, the sum of the order of the\
                alpha-unstable part and the order of a minimal\
                realization of the alpha-stable part of the given\
                system. In this case, the resulting nr is set equal\
                to nsmin.")
    if out[-2] == 2:
        warnings.warn("with ordsel = 'F', the selected order nr is less\
                than the order of the alpha-unstable part of the\
                given system. In this case nr is set equal to the\
                order of the alpha-unstable part.")
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        e = ArithmeticError('The reduction of A to the real Schur form failed')
        e.info = out[-1]
        raise e
    if out[-1] == 2:
        e = ArithmeticError('the separation of the alpha-stable/unstable diagonal\
                 blocks failed because of very close eigenvalues;')
        e.info = out[-1]
        raise e
    if out[-1] == 3:
        e = ArithmeticError('The computation of Hankel singular values failed')
        e.info = out[-1]
        raise e
    Nr,A,B,C,D,Ns,hsv = out[:-2]
    return Nr, A[:Nr,:Nr], B[:Nr,:], C[:,:Nr], D, Ns, hsv

def ab13bd(dico, jobn, n, m, p, A, B, C, D, tol = 0.0):
    """norm = ab13bd(dico, jobn, n, m, p, A, B, C, D, [tol])

    To compute the H2 or L2 norm of the transfer-function matrix G
    of the system (A,B,C,D). G must not have poles on the imaginary
    axis, for a continuous-time system, or on the unit circle, for
    a discrete-time system. If the H2-norm is computed, the system
    must be stable.

    Required arguments:
        dico : {'D', 'C'} input string(len=1)
               Indicate whether the system is discrete 'D' or continuous 'C'.
        jobn : {'H', 'L'} input string(len=1)
               H2-norm 'H' or L2-norm 'L' to be computed.
        n : input int
            The number of state variables.  n >= 0.
        m : input int
            The number of system inputs.  m >= 0.
        p : input int
            The number of system outputs.  p >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state
            dynamics matrix A of the system.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input/state
            matrix B of the system.
        C : input rank-2 array('d') with bounds (p,n)
            The leading p-by-n part of this array must contain the state/output
            matrix C of the system.
        D : input rank-2 array('d') with bounds (p,m)
            The leading p-by-m part of this array must contain the direct
            transmission matrix D of the system.

    Optional arguments:
        tol : The absolute tolerance level below which the elements of
              B are considered zero (used for controllability tests).
              If the user sets tol <= 0, then an implicitly computed,
              default tolerance, defined by  toldef = n*eps*norm(B),
              is used instead, where eps is the machine precision
              (see LAPACK Library routine DLAMCH) and norm(B) denotes
              the 1-norm of B.
    """

    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be "C" or "D"')
    if jobn != 'H' and jobn != 'L':
        raise ValueError('jobn must be "H" or "L"')
    out = _wrapper.ab13bd(dico, jobn, n, m, p, A, B, C, D, tol)
    if out[-1] == 0:
        # success
        return out[0]
    elif out[-1] < 0:
        hidden = ' (hidden by the wrapper)'
        arg_list = ['dico', 'jobn', 'n', 'm', 'p',
                    'A', 'lda'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'D', 'ldd'+hidden,
                    'nq'+hidden,'tol', 'dwork'+hidden, 'ldwork'+hidden, 'iwarn'+hidden, 'info'+hidden]
        error_text = "The following argument had an illegal value: " + arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] == 1:
        e = ArithmeticError("the reduction of A to a real Schur form failed")
        e.info = out[-1]
        raise e
    elif out[-1] == 2:
        e = ArithmeticError("a failure was detected during the reordering of the real Schur form of A, "
                            "or in the iterative process for reordering the eigenvalues of "
                            "Z'*(A + B*F)*Z along the diagonal")
        e.info = out[-1]
        raise e
    elif out[-1] == 3:
        e = ArithmeticError("the matrix A has a controllable eigenvalue on the " +
                            ( "imaginary axis" if dico == 'C' else "unit circle" ))
        e.info = out[-1]
        raise e
    elif out[-1] == 4:
        e = ArithmeticError("the solution of Lyapunov equation failed because the equation is singular")
        e.info = out[-1]
        raise e
    elif out[-1] == 5:
        e = ArithmeticError("DICO = 'C' and D is a nonzero matrix")
        e.info = out[-1]
        raise e
    elif out[-1] == 6:
        e = ArithmeticError("JOBN = 'H' and the system is unstable")
        e.info = out[-1]
        raise e
    else:
        raise RuntimeError("unknown error code %r" % out[-1])

def ab13dd(dico, jobe, equil, jobd, n, m, p, A, E, B, C, D, tol = 1e-10):
    """gpeak, fpeak = ab13dd(dico, jobe, equil, jobd, n, m, p, A, E, B, C, D, [tol])

    To compute the L-infinity norm of a continuous-time or
    discrete-time system, either standard or in the descriptor form,

                                  -1
     G(lambda) = C*( lambda*E - A ) *B + D .

    The norm is finite if and only if the matrix pair (A,E) has no
    eigenvalue on the boundary of the stability domain, i.e., the
    imaginary axis, or the unit circle, respectively. It is assumed
    that the matrix E is nonsingular.

    Required arguments:
        dico : {'D', 'C'} input string(len=1)
               Indicate whether the system is discrete 'D' or continuous 'C'.
        jobe : {'G', 'I'} input string(len=1)
               Specifies whether E is a general square or an identity
               matrix, as follows:
               = 'G':  E is a general square matrix;
               = 'I':  E is the identity matrix.
        equil : {'S', 'N'} input string(len=1)
                Specifies whether the user wishes to preliminarily
                equilibrate the system (A,E,B,C) or (A,B,C), as follows:
                = 'S':  perform equilibration (scaling);
                = 'N':  do not perform equilibration.
        jobd : {'D', 'Z'} input string(len=1)
               Specifies whether or not a non-zero matrix D appears in
               the given state space model:
               = 'D':  D is present;
               = 'Z':  D is assumed a zero matrix.
        n : input int
            The number of state variables.  n >= 0.
        m : input int
            The number of system inputs.  m >= 0.
        p : input int
            The number of system outputs.  p >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the state
            dynamics matrix A of the system.
        E : input rank-2 array('d') with bounds (n,n)
            If jobe = 'G', the leading N-by-N part of this array must
            contain the descriptor matrix E of the system.
            If jobe = 'I', then E is assumed to be the identity
            matrix and is not referenced.
        B : input rank-2 array('d') with bounds (n,m)
            The leading n-by-m part of this array must contain the input/state
            matrix B of the system.
        C : input rank-2 array('d') with bounds (p,n)
            The leading p-by-n part of this array must contain the state/output
            matrix C of the system.
        D : input rank-2 array('d') with bounds (p,m)
            The leading p-by-m part of this array must contain the direct
            transmission matrix D of the system.

    Optional arguments:
        tol : Tolerance used to set the accuracy in determining the
              norm.  0 <= tol < 1.

    Return objects:
        gpeak : float
                The L-infinity norm of the system, i.e., the peak gain
                of the frequency response (as measured by the largest
                singular value in the MIMO case).
        fpeak : float
                The frequency where the gain of the frequency response
                achieves its peak value gpeak, i.e.,

                  || G ( j*fpeak ) || = gpeak ,  if dico = 'C', or

                          j*fpeak
                  || G ( e       ) || = gpeak ,  if dico = 'D'.
    """

    if dico != 'C' and dico != 'D':
        raise ValueError('dico must be "C" or "D"')
    if jobe != 'G' and jobe != 'I':
        raise ValueError('jobe must be "G" or "I"')
    if equil != 'S' and equil != 'N':
        raise ValueError('equil must be "S" or "N"')
    if jobd != 'D' and jobd != 'Z':
        raise ValueError('jobd must be "D" or "Z"')
    out = _wrapper.ab13dd(dico, jobe, equil, jobd, n, m, p, [0.0, 1.0], A, E, B, C, D, tol)
    if out[-1] == 0:
        # success
        fpeak = out[0][0] if out[0][1] > 0 else float('inf')
        gpeak = out[1][0] if out[1][1] > 0 else float('inf')
        return gpeak, fpeak
    elif out[-1] < 0:
        hidden = ' (hidden by the wrapper)'
        arg_list = ['dico', 'jobe', 'equil', 'jobd', 'n', 'm', 'p', 'fpeak'+hidden,
                    'A', 'lda'+hidden, 'E', 'lde'+hidden, 'B', 'ldb'+hidden, 'C', 'ldc'+hidden, 'D', 'ldd'+hidden,
                    'gpeak'+hidden, 'tol', 'iwork'+hidden, 'dwork'+hidden, 'ldwork'+hidden,
                    'cwork'+hidden, 'lcwork'+hidden, 'info'+hidden]
        error_text = "The following argument had an illegal value: " + arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] == 1:
        e = ArithmeticError("the matrix E is (numerically) singular")
        e.info = out[-1]
        raise e
    elif out[-1] == 2:
        e = ArithmeticError("the (periodic) QR (or QZ) algorithm for computing eigenvalues did not converge")
        e.info = out[-1]
        raise e
    elif out[-1] == 3:
        e = ArithmeticError("the SVD algorithm for computing singular values did not converge")
        e.info = out[-1]
        raise e
    elif out[-1] == 4:
        e = ArithmeticError("the tolerance is too small and the algorithm did not converge")
        e.info = out[-1]
        raise e
    else:
        raise RuntimeError("unknown error code %r" % out[-1])

def ab13ed(n, A, tol = 9.0):
    """low, high = ab13ed(n, A, [tol])

    To estimate beta(A), the 2-norm distance from a real matrix A to
    the nearest complex matrix with an eigenvalue on the imaginary
    axis. The estimate is given as

         low <= beta(A) <= high,

    where either

         (1 + tol) * low >= high,

    or

         low = 0   and   high = delta,

    and delta is a small number approximately equal to the square root
    of machine precision times the Frobenius norm (Euclidean norm)
    of A. If A is stable in the sense that all eigenvalues of A lie
    in the open left half complex plane, then beta(A) is the distance
    to the nearest unstable complex matrix, i.e., the complex
    stability radius.

    Required arguments:
        n : input int
            The order of the matrix A.  n >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the matrix A.

    Optional arguments:
        tol : Specifies the accuracy with which low and high approximate
              beta(A). If the user sets tol to be less than sqrt(eps),
              where eps is the machine precision (see LAPACK Library
              Routine DLAMCH), then the tolerance is taken to be
              sqrt(eps).
              The recommended value is tol = 9, which gives an estimate
              of beta(A) correct to within an order of magnitude.

    Return objects:
        low : float
              A lower bound for beta(A).
        high : float
               An upper bound for beta(A).
    """
    out = _wrapper.ab13ed(n, A, tol)
    if out[-1] == 0:
        # success
        return out[0], out[1]
    elif out[-1] < 0:
        hidden = ' (hidden by the wrapper)'
        arg_list = ['n', 'A', 'lda'+hidden, 'low'+hidden, 'high'+hidden, 'tol',
                    'dwork'+hidden, 'ldwork'+hidden, 'info'+hidden]
        error_text = "The following argument had an illegal value: " + arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] == 1:
        e = ArithmeticError("the QR algorithm fails to converge")
        e.info = out[-1]
        raise e
    else:
        raise RuntimeError("unknown error code %r" % out[-1])

def ab13fd(n, A, tol = 0.0):
    """beta, omega = ab13fd(n, A, [tol])

    To compute beta(A), the 2-norm distance from a real matrix A to
    the nearest complex matrix with an eigenvalue on the imaginary
    axis. If A is stable in the sense that all eigenvalues of A lie
    in the open left half complex plane, then beta(A) is the complex
    stability radius, i.e., the distance to the nearest unstable
    complex matrix. The value of beta(A) is the minimum of the
    smallest singular value of (A - jwI), taken over all real w.
    The value of w corresponding to the minimum is also computed.

    Required arguments:
        n : input int
            The order of the matrix A.  n >= 0.
        A : input rank-2 array('d') with bounds (n,n)
            The leading n-by-n part of this array must contain the matrix A.

    Optional arguments:
        tol : Specifies the accuracy with which beta(A) is to be
              calculated. (See the Numerical Aspects section below.)
              If the user sets tol to be less than eps, where eps is the
              machine precision (see LAPACK Library Routine DLAMCH),
              then the tolerance is taken to be eps.

    Return objects:
        beta : float
               The computed value of beta(A), which actually is an upper
               bound.
        omega : float
                The value of w such that the smallest singular value of
                (A - jwI) equals beta(A).

    Numerical Aspects:
        In the presence of rounding errors, the computed function value
        beta satisfies
              beta(A) <= beta + epsilon,
              beta/(1+tol) - delta <= max(beta(A), sqrt(2*n*eps)*norm(A)),
        where norm(A) is the Frobenius norm of A,
              epsilon = p(n) * eps * norm(A),
        and
              delta   = p(n) * sqrt(eps) * norm(A),
        and p(n) is a low degree polynomial. It is recommended to choose
        tol greater than sqrt(eps). Although rounding errors can cause
        AB13FD to fail for smaller values of tol, nevertheless, it usually
        succeeds. Regardless of success or failure, the first inequality
        holds.
    """
    out = _wrapper.ab13fd(n, A, tol)
    if out[-1] == 0:
        # success
        return out[0], out[1]
    elif out[-1] < 0:
        hidden = ' (hidden by the wrapper)'
        arg_list = ['n', 'A', 'lda'+hidden, 'beta'+hidden, 'omega'+hidden, 'tol',
                    'dwork'+hidden, 'ldwork'+hidden, 'cwork'+hidden, 'lcwork'+hidden, 'info'+hidden]
        error_text = "The following argument had an illegal value: " + arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] == 1:
        warnings.warn("the routine fails to compute beta(A) within the specified tolerance")
        return out[0], out[1]  # the returned value is an upper bound on beta(A)
    elif out[-1] == 2:
        e = ArithmeticError("either the QR or SVD algorithm fails to converge")
        e.info = out[-1]
        raise e
    else:
        raise RuntimeError("unknown error code %r" % out[-1])

# to be replaced by python wrappers
