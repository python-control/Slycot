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

import numpy as np

from . import _wrapper
from .exceptions import raise_if_slycot_error, SlycotParameterError


def ab01nd(n, m, A, B, jobz='N', tol=0, ldwork=None):
    """ Ac,Bc,ncont,indcon,nblk,Z,tau = ab01nd_i(n,m,A,B,[jobz,tol,ldwork])

    To find a controllable realization for the linear time-invariant
    multi-input system

    ::

          dX/dt = A * X + B * U,

    where A and B are N-by-N and N-by-M matrices, respectively,
    which are reduced by this routine to orthogonal canonical form
    using (and optionally accumulating) orthogonal similarity
    transformations.  Specifically, the pair ``(A, B)`` is reduced to
    the pair ``(Ac, Bc),  Ac = Z' * A * Z,  Bc = Z' * B``,  given by

    ::

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

    where the blocks ``B1, A21, ..., Ap,p-1`` have full row ranks and
    `p` is the controllability index of the pair.  The size of the
    block `Auncont` is equal to the dimension of the uncontrollable
    subspace of the pair ``(A, B)``.

    Parameters
    ----------
    n : int
        The order of the original state-space representation, i.e.
        the order of the matrix A.  ``n > 0``.
    m : int
        The number of system inputs, or of columns of B.  ``m > 0``.
    A : (n, n) array_like
        The original state dynamics matrix A.
    B : (n, m) array_like
        The input matrix B.
    jobz : {'N', 'F', 'I'}, optional
        Indicates whether the user wishes to accumulate in a matrix Z
        the orthogonal similarity transformations for reducing the system,
        as follows:
        := 'N':  Do not form Z and do not store the orthogonal transformations;
                 (default)
        := 'F':  Do not form Z, but store the orthogonal transformations in
                 the factored form;
        := 'I':  Z is initialized to the unit matrix and the orthogonal
                 transformation matrix Z is returned.
    tol : float, optional
        The tolerance to be used in rank determination when transforming
        ``(A, B)``. If ``tol <= 0`` a default value is used.
    ldwork : int, optional
        The length of the cache array. ``ldwork >= max(n, 3*m)``.
        For optimum performance it should be larger.
        default: ``ldwork = max(n, 3*m)``

    Returns
    -------
    Ac : (n, n) ndarray
        The leading ncont-by-ncont part contains the upper block
        Hessenberg state dynamics matrix Acont in Ac, given by Z'*A*Z,
        of a controllable realization for the original system. The
        elements below the first block-subdiagonal are set to zero.
    Bc : (n, m) ndarray
        The leading ncont-by-m part of this array contains the transformed
        input matrix Bcont in Bc, given by ``Z'*B``, with all elements but the
        first block set to zero.
    ncont : int
        The order of the controllable state-space representation.
    indcon : int
        The controllability index of the controllable part of the system
        representation.
    nblk : (n, ) int ndarray
        The leading indcon elements of this array contain the the orders of
        the diagonal blocks of Acont.
    Z : (n, n) ndarray
        - If jobz = 'I', then the leading N-by-N part of this array contains
          the matrix of accumulated orthogonal similarity transformations
          which reduces the given system to orthogonal canonical form.
        - If jobz = 'F', the elements below the diagonal, with the array tau,
          represent the orthogonal transformation matrix as a product of
          elementary reflectors. The transformation matrix can then be
          obtained by calling the LAPACK Library routine DORGQR.
        - If jobz = 'N', the array Z is `None`.
    tau : (n, ) ndarray
        The elements of tau contain the scalar factors of the
        elementary reflectors used in the reduction of B and A.
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['jobz', 'n', 'm', 'A', 'LDA'+hidden, 'B', 'LDB'+hidden,
                'ncont', 'indcon', 'nblk', 'Z', 'LDZ'+hidden, 'tau', 'tol',
                'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 'info'+hidden]

    if ldwork is None:
        ldwork = max(n, 3*m)

    Ac, Bc, ncont, indcon, nblk, Z, tau, info = _wrapper.ab01nd(
        jobz, n, m, A, B, tol=tol, ldwork=ldwork)
    raise_if_slycot_error(info, arg_list)

    if jobz == "N":
        Z = None
    return Ac, Bc, ncont, indcon, nblk, Z, tau

def ab04md(type_t, n, m, p, A, B, C, D, alpha=1.0, beta=1.0, ldwork=None):
    """ At,Bt,Ct,Dt = ab04md(type_t, n, m, p, A, B, C, D, [alpha, beta,ldwork])

    Parameters
    ----------
    type_t : {'D','C'}
            Indicates the type of the original system and the
            transformation to be performed as follows:
            = 'D':  discrete-time   -> continuous-time;
            = 'C':  continuous-time -> discrete-time.
    n : int
        The order of the matrix A, the number of rows of matrix B and
        the number of columns of matrix C. It represents the dimension of
        the state vector.  n > 0.        
    m : int
        The number of columns of matrix B. It represents the dimension of
        the input vector.  m > 0.
    p : int
        The number of rows of matrix C. It represents the dimension of
        the output vector.  p > 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the system state
        matrix A.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the system input
        matrix B.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the system output
        matrix C.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the system direct
        transmission matrix D.
    alpha : float, optional
        Parameter specifying the bilinear transformation.
        Recommended values for stable systems: alpha = 1, alpha != 0,
        Default is 1.0.
    beta : float, optional
        Parameter specifying the bilinear transformation.
        Recommended values for stable systems: beta = 1, beta != 0,
        Default is 1.0.
    ldwork : int, optional
        The length of the cache array.
        ldwork >= max(1, n), default is max(1, n)
    Returns
    -------
    At : (n, n) ndarray
        The state matrix At of the transformed system.
    Bt : (n, m) ndarray
        The input matrix Bt of the transformed system.
    Ct : (p, n) ndarray
        The output matrix Ct of the transformed system.
    Dt : (p, m) ndarray
        The transmission matrix Dt of the transformed system.
    Raises
    ------
    SlycotArithmeticError
        :info == 1: 
            If the matrix (ALPHA*I + A) is exactly singular
        :info == 2: 
            If the matrix  (BETA*I - A) is exactly singular.
    """

    hidden = ' (hidden by the wrapper)'  
    arg_list = ['type_t', 'n', 'm', 'p', 'alpha', 'beta', 
                'A', 'LDA'+hidden, 'B', 'LDB'+hidden, 'C', 'LDC'+hidden, 'D', 'LDD'+hidden,
                'IWORK'+hidden, 'DWORK'+hidden, 'ldwork', 'info'+hidden]

    if ldwork is None:
        ldwork = max(1, n)

    out = _wrapper.ab04md(type_t, n, m, p, alpha, beta, A, B, C, D, ldwork=ldwork)
    info=out[-1]
    raise_if_slycot_error(info, arg_list)

    return out[:-1]

def ab05md(n1,m1,p1,n2,p2,A1,B1,C1,D1,A2,B2,C2,D2,uplo='U'):
    """ n,a,b,c,d = ab05md(n1,m1,p1,n2,p2,a1,b1,c1,d1,a2,b2,c2,d2,[uplo])

    To obtain the state-space model (A,B,C,D) for the cascaded
    inter-connection of two systems, each given in state-space form.

    Parameters
    ----------
    n1 : int
        The number of state variables in the first system, i.e. the order
        of the matrix A1.  n1 > 0.
    m1 : int
        The number of input variables for the first system. m1 > 0.
    p1 : int
        The number of output variables from the first system and the number
        of input variables for the second system. p1 > 0.
    n2 : int
        The number of state variables in the second system, i.e. the order
        of the matrix A2.  n2 > 0.
    p2 : int
        The number of output variables from the second system. p2 > 0.
    A1 : (n1, n1) array_like
        The leading n1-by-n1 part of this array must contain the state
        transition matrix A1 for the first system.
    B1 : (n1, m1) array_like
        The leading n1-by-m1 part of this array must contain the input/state
        matrix B1 for the first system.
    C1 : (p1, n1) array_like
        The leading p1-by-n1 part of this array must contain the state/output
        matrix C1 for the first system.
    D1 : (p1, m1) array_like
        The leading p1-by-m1 part of this array must contain the input/output
        matrix D1 for the first system.
    A2 : (n2, n2) array_like
        The leading n2-by-n2 part of this array must contain the state
        transition matrix A2 for the second system.
    B2 : (n2, p1) array_like
        The leading n2-by-p1 part of this array must contain the input/state
        matrix B2 for the second system.
    C2 : (p2, n2) array_like
        The leading p2-by-n2 part of this array must contain the state/output
        matrix C2 for the second system.
    D2 : (p2, p1) array_like
        The leading p2-by-p1 part of this array must contain the input/output
        matrix D2 for the second system.
    uplo :  {'U', 'L'}, optional
        Indicates whether the user wishes to obtain the matrix A in
        the upper or lower block diagonal form, as follows:
            = 'U':  Obtain A in the upper block diagonal form;
            = 'L':  Obtain A in the lower block diagonal form.
        Default is `U`.

    Returns
    -------
    n : int
        The number of state variables (n1 + n2) in the resulting system,
        i.e. the order of the matrix A, the number of rows of B and
        the number of columns of C.
    A : (n1+n2, n1+n2) ndarray
        The leading N-by-N part of this array contains the state transition
        matrix A for the cascaded system.
    B : (n1+n2, m1) ndarray
        The leading n-by-m1 part of this array contains the input/state
        matrix B for the cascaded system.
    C : (p2, n1+n2) ndarray
        The leading p2-by-n part of this array contains the state/output
        matrix C for the cascaded system.
    D : (p2, m1) ndarray
        The leading p2-by-m1 part of this array contains the input/output
        matrix D for the cascaded system.

    Notes
    -----
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
    raise_if_slycot_error(out[-1], arg_list)
    return out[:-1]

def ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,alpha=1.0,ldwork=None):
    """  n,A,B,C,D = ab05nd(n1,m1,p1,n2,A1,B1,C1,D1,A2,B2,C2,D2,[alpha,ldwork])

    To obtain the state-space model (A,B,C,D) for the feedback inter-connection
    of two systems, each given in state-space form.

    Parameters
    ----------
    n1 : int
        The number of state variables in the first system, i.e. the order
        of the matrix A1.  n1 > 0.
    m1 : int
        The number of input variables for the first system and the number
        of output variables from the second system. m1 > 0.
    p1 : int
        The number of output variables from the first system and the number
        of input variables for the second system. p1 > 0.
    n2 : int
        The number of state variables in the second system, i.e. the order
        of the matrix A2.  n2 > 0.
    A1 : (n1, n1) array_like
        The leading n1-by-n1 part of this array must contain the state
        transition matrix A1 for the first system.
    B1 : (n1, m1) array_like
        The leading n1-by-m1 part of this array must contain the input/state
        matrix B1 for the first system.
    C1 : (p1, n1) array_like
        The leading p1-by-n1 part of this array must contain the state/output
        matrix C1 for the first system.
    D1 : (p1, m1) array_like
        The leading p1-by-m1 part of this array must contain the input/output
        matrix D1 for the first system.
    A2 : (n2, n2) array_like
        The leading n2-by-n2 part of this array must contain the state
        transition matrix A2 for the second system.
    B2 : (n2, p1) array_like
        The leading n2-by-p1 part of this array must contain the input/state
        matrix B2 for the second system.
    C2 : (m1, n2) array_like
        The leading m1-by-n2 part of this array must contain the state/output
        matrix C2 for the second system.
    D2 : (m1, p1) array_like
        The leading m1-by-p1 part of this array must contain the input/output
        matrix D2 for the second system.
    alpha : float, optional
        A coefficient multiplying the transfer-function matrix (or the
        output equation) of the second system. i.e alpha = +1 corresponds
        to positive feedback, and alpha = -1 corresponds to negative
        feedback.
        Default is `1.0`.
    ldwork : int, optional
        The length of the cache array. ldwork >= max(p1*p1,m1*m1,n1*p1).
        Default is max(p1*p1,m1*m1,n1*p1).

    Returns
    -------
    n : int
        The number of state variables (n1 + n2) in the connected system, i.e.
        the order of the matrix A, the number of rows of B and the number of
        columns of C.
    A : (n1+n2, n1+n2) ndarray
        The leading n-by-n part of this array contains the state transition
        matrix A for the connected system.
    B : (n1+n2, m1) ndarray
        The leading n-by-m1 part of this array contains the input/state
        matrix B for the connected system.
    C : (p1, n1, n2) ndarray
        The leading p1-by-n part of this array contains the state/output
        matrix C for the connected system.
    D : (p1, m1) ndarray
        The leading p1-by-m1 part of this array contains the input/output
        matrix D for the connected system.

    Raises
    ------
    SlycotArithmeticError
        :1 <= info <= p1:
            the system is not completely controllable. That is, the matrix
            ``(I + ALPHA*D1*D2)`` is exactly singular (the element
            ``U(i,i)```` of the upper triangular factor of ``LU```
            factorization is exactly zero), possibly due to
            rounding errors.
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
    raise_if_slycot_error(out[-1], arg_list, ab05nd, locals())
    return out[:-1]

def ab07nd(n,m,A,B,C,D,ldwork=None):
    """ Ai,Bi,Ci,Di,rcond = ab07nd(n,m,A,B,C,D,[ldwork])

    To compute the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).

    Parameters
    ----------
    n : int
        The order of the state matrix A. n >= 0.
    m : int
        The number of system inputs and outputs. m >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state matrix
        A of the original system.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input matrix
        B of the original system.
    C : (m, n) array_like
        The leading m-by-n part of this array must contain the output matrix
        C of the original system.
    D : (m, m) array_like
        The leading m-by-m part of this array must contain the feedthrough
        matrix D of the original system.
    ldwork : int, optional
        The length of the cache array. The default value is max(1,4*m),
        for better performance should be larger.

    Returns
    -------
    Ai : (n, n) ndarray
        The leading n-by-n part of this array contains the state matrix Ai
        of the inverse system.
    Bi : (n, m) ndarray
        The leading n-by-m part of this array contains the input matrix Bi
        of the inverse system.
    Ci : (m, n) ndarray
        The leading m-by-n part of this array contains the output matrix Ci
        of the inverse system.
    Di : (m, m) ndarray
        The leading m-by-m part of this array contains the feedthrough
        matrix Di of the inverse system.
    rcond : float
        The estimated reciprocal condition number of the feedthrough matrix
        D of the original system.

    Warns
    -----
    SlycotResultWarning
        :1 <= info <= m:
            the matrix `D` is exactly singular; the ({info},{info})
            diagonal element is zero, `RCOND` was set to zero;
        :info == m+1:
            the matrix `D` is numerically singular, i.e., `RCOND`
            is less than the relative machine precision, `EPS`
            (see LAPACK Library routine DLAMCH). The
            calculations have been completed, but the results
            could be very inaccurate.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'm', 'A', 'LDA' + hidden, 'B', 'LDB' + hidden,
                'C', 'LDC' + hidden, 'D', 'LDD' + hidden, 'rcond',
                'IWORK' + hidden, 'DWORK' + hidden, 'ldwork', 'INFO' + hidden]
    if ldwork is None:
        ldwork = max(1, 4*m)
    out = _wrapper.ab07nd(n, m, A, B, C, D, ldwork=ldwork)
    raise_if_slycot_error(out[-1], arg_list, ab07nd.__doc__, locals())
    return out[:-1]


def ab08nd(n,m,p,A,B,C,D,equil='N',tol=0,ldwork=None):
    """ nu,rank,dinfz,nkror,nkrol,infz,kronr,kronl,Af,Bf = ab08nd(n,m,p,A,B,C,D,[equil,tol,ldwork])

    To construct for a linear multivariable system described by a state-space
    model (A,B,C,D) a regular pencil (Af - lambda*Bf ) which has the invariant
    zeros of the system as generalized eigenvalues.
    The routine also computes the orders of the infinite zeros and the
    right and left Kronecker indices of the system (A,B,C,D).

    Parameters
    ----------
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the state/output
        matrix C of the system.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the direct
        transmission matrix D of the system.
    equil : {'S', 'N'}, optional
        Specifies whether the user wishes to balance the compound matrix
        as follows:
        = 'S':  Perform balancing (scaling);
        = 'N':  Do not perform balancing.
        Default is `N`.
    tol : float, optional
        A tolerance used in rank decisions to determine the effective rank,
        which is defined as the order of the largest leading (or trailing)
        triangular submatrix in the QR (or RQ) factorization with column
        (or row) pivoting whose estimated condition number is less than 1/tol.
        Default is `0.0`.
    ldwork : int, optional
        The length of the cache array. The default value is n + 3*max(m,p),
        for better performance should be larger.
        Default is None.

    Returns
    -------
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
    infz : (n, ) ndarray
        The leading dinfz elements of infz contain information on the
        infinite elementary divisors as follows: the system has infz(i)
        infinite elementary divisors of degree i, where i = 1,2,...,dinfz.
    kronr :(max(n,m)+1, ) ndarray
        the leading nkror elements of this array contain the right kronecker
        (column) indices.
    kronl : (max(n,p)+1, ) ndarray
        the leading nkrol elements of this array contain the left kronecker
        (row) indices.
    Af : (max(1,n+m), n+min(p,m)) ndarray
        the leading nu-by-nu part of this array contains the coefficient
        matrix Af of the reduced pencil. the remainder of the leading
        (n+m)-by-(n+min(p,m)) part is used as internal workspace.
    Bf : (max(1,n+p), n+m) ndarray
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
    raise_if_slycot_error(out[-1], arg_list)
    return out[:-1]

def ab08nz(n, m, p, A, B, C, D, equil='N', tol=0., lzwork=None):
    """ nu,rank,dinfz,nkror,nkrol,infz,kronr,kronl,Af,Bf = ab08nz(n,m,p,A,B,C,D,[equil,tol,lzwork])

    To construct for a linear multivariable system described by a state-space
    model (A,B,C,D) a regular pencil (Af - lambda*Bf ) which has the invariant
    zeros of the system as generalized eigenvalues.
    The routine also computes the orders of the infinite zeros and the
    right and left Kronecker indices of the system (A,B,C,D).

    Parameters
    ----------
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the state/output
        matrix C of the system.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the direct
        transmission matrix D of the system.
    equil : {'S', 'N'}, optional
        Specifies whether the user wishes to balance the compound matrix
        as follows:
        = 'S':  Perform balancing (scaling);
        = 'N':  Do not perform balancing.
        Default is `N`.
    tol : float, optional
        A tolerance used in rank decisions to determine the effective rank,
        which is defined as the order of the largest leading (or trailing)
        triangular submatrix in the QR (or RQ) factorization with column
        (or row) pivoting whose estimated condition number is less than 1/tol.
        If tol is set to less than SQRT((N+P)*(N+M))*EPS
        then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
        where EPS is the machine precision (see LAPACK Library
        Routine DLAMCH).
        Default is 0.0.
    lzwork : int, optional
        The length of the internal cache array ZWORK. The default value is
        calculated to
            MAX( 1,
                MIN(P,M) + MAX(3*M-1,N),
                MIN(P,N) + MAX(3*P-1,N+P,N+M),
                MIN(M,N) + MAX(3*M-1,N+M) )
        For optimum performance lzwork should be larger.
        If lzwork = -1, then a workspace query is assumed;
        the routine only calculates the optimal size of the
        ZWORK array, and returns this value in lzwork_opt
        Default is None.

    Returns
    -------
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
    infz : (n, ) ndarray
        The leading dinfz elements of infz contain information on the
        infinite elementary divisors as follows: the system has infz(i)
        infinite elementary divisors of degree i, where i = 1,2,...,dinfz.
    kronr : (max(n,m)+1, ) ndarray
        the leading nkror elements of this array contain the right kronecker
        (column) indices.
    kronl : (max(n,p)+1, ) ndarray
        the leading nkrol elements of this array contain the left kronecker
        (row) indices.
    Af : (max(1,n+m), n+min(p,m)) ndarray
        the leading nu-by-nu part of this array contains the coefficient
        matrix Af of the reduced pencil. the remainder of the leading
        (n+m)-by-(n+min(p,m)) part is used as internal workspace.
    Bf : (max(1,n+p), n+m) ndarray
        The leading nu-by-nu part of this array contains the coefficient
        matrix Bf of the reduced pencil. the remainder of the leading
        (n+p)-by-(n+m) part is used as internal workspace.
    lzwork_opt : int
        The optimal value of lzwork.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['equil', 'n', 'm', 'p',
                'a', 'lda' + hidden, 'b', 'ldb' + hidden,
                'c', 'ldc' + hidden, 'd', 'ldd' + hidden,
                'nu', 'rank', 'dinfz', 'nkror', 'nkrol', 'infz', 'kronr',
                'kronl', 'af', 'ldaf' + hidden, 'bf', 'ldbf' + hidden,
                'tol', 'iwork' + hidden, 'dwork' + hidden, 'zwork',
                'lzwork', 'info']
    if lzwork is None:
        lzwork = max(min(p, m) + max(3*m-1, n),
                     min(p, n) + max(3*p-1, n+p, n+m),
                     min(m, n) + max(3*m-1, n+m))

    nu, rank, dinfz, nkror, nkrol, infz, kronr, kronl, Af, Bf, zwork, info \
        = _wrapper.ab08nz(n, m, p, A, B, C, D,
                          equil=equil, tol=tol, lzwork=lzwork)

    raise_if_slycot_error(info, arg_list)
    return (nu, rank, dinfz, nkror, nkrol, infz, kronr, kronl, Af, Bf,
            int(zwork[0].real))


def ab09ad(dico,job,equil,n,m,p,A,B,C,nr=None,tol=0,ldwork=None):
    """ nr,Ar,Br,Cr,hsv = ab09ad(dico,job,equil,n,m,p,A,B,C,[nr,tol,ldwork])

    Compute reduced order State-Space-Model (Ar, Br, Cr) for a stable system
    (A, B, C) by using either the square-root or the balancing-free square-
    root Balance & truncate (B & T) model reduction method.

    Parameters
    ----------
    dico : {'D', 'C'}
        Indicate whether the system is discrete `D` or continuous `C`
    job : {'B', 'N'}
        Balance `B` or not `N`
    equil : {'S', 'N'}
        Scale `S` or not `N`
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the
        state/output matrix C of the system.
    nr : int, optional
        `nr` is the desired order of the resulting reduced order
        system.  ``0 <= nr <= n``. Automatically determined by `tol` if
        ``nr is None`` and returned. See return object `nr`.
        Default is None.
    tol : float, optional
        If ``nr is None``, `tol`contains the tolerance for determining the
        order of the reduced system. For model reduction, th recommended
        value is ``tol = c * HNORM(A, B, C)``, where `c` is a constan in the
        interval ``[0.00001, 0.001]`` and ``HNORM(A, B, C)`` is the
        Hankel-Norm of the given sysstem (computed in ``HSV(1)``). For
        computing a minimal realization, the recommended value is
        ``tol = n * eps * HNORM(A, B, C)``, where `eps` is the machine
        precision (see LAPACK Library Routine `DLAMCH`). This value is
        used by default if ``tol <= 0`` on entry. If `nr` is specified,
        the value of `tol` is ignored. Default is `0.0`.
    ldwork : int, optional
        The length of the cache array. The default value is
        ``n*(2*n+max(n,m,p)+5) + n*(n+1)/2 ~= 3.5*n**2 + 5*n``,
        a larger value should lead to better performance.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) ndarray
        This array contains the state dynamics matrix `Ar` of the reduced
        order system.
    Br : (nr, m) ndarray
        This array contains the input/state matrix `Br` of the reduced
        order system.
    Cr : (p, nr) ndarray
        This array contains the state/output matrix `Cr` of the reduced
        order system.
    hsv : (n, ) ndarray
        If ``INFO = 0``, it contains the Hankel singular values of
        the original system ordered decreasingly. ``HSV(1)`` is the
        Hankel norm of the system.

    Raises
    ------
    SlycotArithmeticError
        :info == 1:
            The reduction of A to the real Schur form failed
        :info == 2 and dico == 'C':
            The state matrix A is not stable
        :info == 2 and dico == 'D':
            The state matrix A is not convergent
        :info == 3:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
                The selected order {nr} is greater
                than the order of a minimal realization of the
                given system. `nr` was set automatically to {Nr}
                corresponding to the order of a minimal realization
                of the system
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'hsv', 'tol', 'iwork' + hidden, 'dwork ' + hidden, 'ldwork',
                'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n+max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    out = _wrapper.ab09ad(dico, job, equil, ordsel,
                          n, m, p, nr, A, B, C, tol, ldwork)
    Nr, A, B, C, hsv = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09ad.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], hsv


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

    Parameters
    ----------
    dico : {'D', 'C'}
        Indicate whether the system is discrete `D` or continuous `C`
    job : {'B', 'N'}
        Balance `B` or not `N`
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system *in real Schur form.*
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the
        state/output matrix C of the system.
    nr : int, optional
        `nr` is the desired order of the resulting reduced order
        system.  ``0 <= nr <= n``. Automatically determined by `tol` if
        ``nr is None`` and returned. See return object `nr`.
        Default is None.
    tol : float, optional
        If ``nr is None``, `tol`contains the tolerance for determining the
        order of the reduced system. For model reduction, the recommended
        value is ``tol = c * HNORM(A, B, C)``, where `c` is a constant in
        the interval ``[0.00001, 0.001]`` and ``HNORM(A, B, C)`` is
        the Hankel-Norm of the given sysstem (computed in ``HSV(1)``). For
        computing a minimal realization, the recommended value is
        ``tol = n * eps * HNORM(A, B, C)``, where `eps` is the machine
        precision (see LAPACK Library Routine `DLAMCH`). This value is
        used by default if ``tol <= 0`` on entry. If `nr` is specified,
        the value of `tol` is ignored. Default is `0.0`.
    ldwork : int, optional
        The length of the cache array. The default value is
        ``n*(2*n+max(n,m,p)+5) + n*(n+1)/2 ~= 3.5*n**2 + 5*n``,
        a larger value should lead to better performance.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) ndarray
        This array contains the state dynamics matrix `Ar` of the reduced
        order system.
    Br : (nr, m) ndarray
        Tthis array contains the input/state matrix `Br` of the reduced
        order system.
    Cr : (p, nr) ndarray
        This array contains the state/output matrix `Cr` of the reduced
        order system.
    hsv : (n, ) ndarray
        If ``INFO = 0``, it contains the Hankel singular values of
        the original system ordered decreasingly. ``HSV(1)`` is the
        Hankel norm of the system.
    T : (n, nr) ndarray
        This array contains the right truncation matrix `T` of the reduced
        order system.
    Ti : (nr, n) ndarray
        This array contains the left truncation matrix `Ti` of the reduced
        order system.

    Raises
    ------
    SlycotArithmeticError
        :info == 1 and dico == 'C':
            The state matrix A is not stable
        :info == 1 and dico == 'D':
            The state matrix A is not convergent
        :info == 2:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
                The selected order {nr} is greater
                than the order of a minimal realization of the
                given system. `nr` was set automatically to {Nr}
                corresponding to the order of a minimal realization
                of the system
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'ordsel', 'n', 'm', 'p', 'nr',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'hsv', 'T', 'ldt' + hidden, 'Ti', 'ldti' + hidden, 'tol',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n + max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    out = _wrapper.ab09ax(dico, job, ordsel, n, m, p, nr, A, B, C, tol, ldwork)
    Nr, A, B, C, hsv, T, Ti = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09ax.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], hsv, T[:, :Nr], Ti[:Nr, :]


def ab09bx(dico,job,n,m,p,A,B,C,D,nr=None,tol1=0.0,tol2=0.0,ldwork=None):
    """``nr,Ar,Br,Cr,Dr,hsv,T,Ti = ab09bx(dico,job,equil,n,m,p,nr,A,B,C,D,[nr,tol1,tol2,ldwork])``

    TODO update docs for ab09bx

    To compute a reduced order model ``(Ar,Br,Cr,Dr)`` for a stable original
    state-space representation ``(A,B,C,D)`` by using either the square-root
    or the balancing-free square-root Singular Perturbation model
    reduction method. The state dynamics matrix `A` of the original
    system is an upper quasi-triangular matrix in *real Schur canonical
    form.* The matrices of the reduced order system are computed using
    the truncation formulas:

        ``Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T`` .

    Parameters
    ----------
    dico : {'D', 'C'}
        Indicate whether the system is discrete `D` or continuous `C`
    job : {'B', 'N'}
        Balance `B` or not `N`
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system *in real Schur form.*
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the
        state/output matrix C of the system.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the
        state/output matrix D of the system.
    nr : int, optional
        `nr` is the desired order of the resulting reduced order
        system.  ``0 <= nr <= n``. Automatically determined by `tol` if
        ``nr is None`` and returned. See return object `nr`.
        Default is None.
    tol2 : float, optional
        If ``nr is None``, `tol`contains the tolerance for determining the
        order of the reduced system. For model reduction, the recommended
        value is ``tol = c * HNORM(A, B, C)``, where `c` is a constant in
        the interval ``[0.00001, 0.001]`` and ``HNORM(A, B, C)`` is
        the Hankel-Norm of the given sysstem (computed in ``HSV(1)``). For
        computing a minimal realization, the recommended value is
        ``tol = n * eps * HNORM(A, B, C)``, where `eps` is the machine
        precision (see LAPACK Library Routine `DLAMCH`). This value is
        used by default if ``tol <= 0`` on entry. If `nr` is specified,
        the value of `tol` is ignored. Default is `0.0`.
    tol2 : float, optional
        additional tolerance TODO update but slicot docs not illuminating
    ldwork : int, optional
        The length of the cache array. The default value is
        ``n*(2*n+max(n,m,p)+5) + n*(n+1)/2 ~= 3.5*n**2 + 5*n``,
        a larger value should lead to better performance.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) ndarray
        This array contains the state dynamics matrix `Ar` of the reduced
        order system.
    Br : (nr, m) ndarray
        Tthis array contains the input/state matrix `Br` of the reduced
        order system.
    Cr : (p, nr) ndarray
        This array contains the state/output matrix `Cr` of the reduced
        order system.
    Dr : (p, m) ndarray
        This array contains the state/output matrix `Dr` of the reduced
        order system.
    hsv : (n, ) ndarray
        If ``INFO = 0``, it contains the Hankel singular values of
        the original system ordered decreasingly. ``HSV(1)`` is the
        Hankel norm of the system.
    T : (n, nr) ndarray
        This array contains the right truncation matrix `T` of the reduced
        order system.
    Ti : (nr, n) ndarray
        This array contains the left truncation matrix `Ti` of the reduced
        order system.

    Raises
    ------
    SlycotArithmeticError
        :info == 1 and dico == 'C':
            The state matrix A is not stable
        :info == 1 and dico == 'D':
            The state matrix A is not convergent
        :info == 2:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
                The selected order {nr} is greater
                than the order of a minimal realization of the
                given system. `nr` was set automatically to {Nr}
                corresponding to the order of a minimal realization
                of the system
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'ordsel', 'n', 'm', 'p', 'nr',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'hsv', 'T', 'ldt' + hidden, 'Ti', 'ldti' + hidden, 'tol',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n + max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    out = _wrapper.ab09bx(dico, job, ordsel, n, m, p, nr, A, B, C, D, tol1, tol2, ldwork)
    Nr, A, B, C, D, hsv, T, Ti = out[:-2]  # 
    raise_if_slycot_error(out[-2:], arg_list, ab09bx.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], D, hsv, T[:, :Nr], Ti[:Nr, :]


# TODO wrap ab09hx, the interface is similar to ab09bx

def ab09bd(dico,job,equil,n,m,p,A,B,C,D,nr=None,tol1=0,tol2=0,ldwork=None):
    """ nr,Ar,Br,Cr,Dr,hsv = ab09bd(dico,job,equil,n,m,p,A,B,C,D,[nr,tol1,tol2,ldwork])

    To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
    original state-space representation (A,B,C,D) by using either the
    square-root or the balancing-free square-root Singular
    Perturbation Approximation (SPA) model reduction method.
    Must supply either nr or tolerance values.

    Parameters
    ----------
    dico : {'C', 'D'}
        Specifies the type of the original system as follows:
        = 'C':  continuous-time system;
        = 'D':  discrete-time system.
    job : {'B', 'N'}
        Specifies the model reduction approach to be used
        as follows:
        = 'B':  use the square-root SPA method;
        = 'N':  use the balancing-free square-root SPA method.
    equil : {'S', 'N'}
        Specifies whether the user wishes to preliminarily
        equilibrate the triplet (A,B,C) as follows:
        = 'S':  perform equilibration (scaling);
        = 'N':  do not perform equilibration.
    n : int
        The order of the original state-space representation, i.e.
        the order of the matrix A. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        On entry, the leading n-by-n part of this array must
        contain the state dynamics matrix A.
    B : (n, m) array_like
        On entry, the leading n-by-m part of this array must
        contain the original input/state matrix B.
    C : (p, n) array_like
        On entry, the leading p-by-n part of this array must
        contain the original state/output matrix C.
    D : (p, m) array_like
        On entry, the leading p-by-m part of this array must
        contain the original input/output matrix D.
    nr : int, optional
        nr is the desired order of
        the resulting reduced order system.  0 <= nr <= n.
        Default is None.
    tol1 : float, optional
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
        Default is `0.0`.
    tol2 : float, optional
        The tolerance for determining the order of a minimal
        realization of the given system. The recommended value is
        tol2 = n*eps*hnorm(A,B,C). This value is used by default
        if tol2 <= 0 on entry.
        If tol2 > 0, then tol2 <= tol1.
        Default is `0.0`.
    ldwork : int, optional
        The length of the cache array. The default value is n + 3*max(m,p),
        for better performance should be larger.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) ndarray
        the leading nr-by-nr part of this array contains the
        state dynamics matrix Ar of the reduced order system.
    Br : (nr, m) ndarray
        the leading nr-by-m part of this array contains the
        input/state matrix Br of the reduced order system.
    Cr : (p, nr) ndarray
        the leading p-by-nr part of this array contains the
        state/output matrix Cr of the reduced order system.
    Dr : (p, m) ndarray
        the leading p-by-m part of this array contains the
        input/output matrix Dr of the reduced order system.
    hsv : (n, ) ndarray
        If info = 0, it contains the Hankel singular values of
        the original system ordered decreasingly. hsv(1) is the
        Hankel norm of the system.

    Raises
    ------
    SlycotArithmeticError
        :info == 1:
            The reduction of A to the real Schur form failed
        :info == 2 and dico == 'C':
            The state matrix A is not stable
        :info == 2 and dico == 'D':
            The state matrix A is not convergent
        :info == 3:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
                The selected order {nr} is greater
                than the order of a minimal realization of the
                given system. `nr` was set automatically to {Nr}
                corresponding to the order of a minimal realization
                of the system
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'D', 'ldd' + hidden, 'hsv', 'tol1', 'tol2',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n+max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    out = _wrapper.ab09bd(dico, job, equil, ordsel,
                          n, m, p, nr, A, B, C, D, tol1, tol2, ldwork)
    Nr, A, B, C, D, hsv = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09bd.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[ :,:Nr], D[:, :], hsv

def ab09md(dico,job,equil,n,m,p,A,B,C,alpha=None,nr=None,tol=0,ldwork=None):
    """ nr,Ar,Br,Cr,ns,hsv = ab09md(dico,job,equil,n,m,p,A,B,C,[alpha,nr,tol,ldwork])

    To compute a reduced order model (Ar,Br,Cr) for an original
    state-space representation (A,B,C) by using either the square-root
    or the balancing-free square-root Balance & Truncate (B & T)
    model reduction method for the ALPHA-stable part of the system.

    Parameters
    ----------
    dico : {'C', 'D'}
        Specifies the type of the original system as follows:
        = 'C':  continuous-time system;
        = 'D':  discrete-time system.
    job : {'B', 'N'}
        Specifies the model reduction approach to be used
        as follows:
        = 'B':  use the square-root Balance & Truncate method;
        = 'N':  use the balancing-free square-root
                Balance & Truncate method.
    equil : {'S', 'N'}
        Specifies whether the user wishes to preliminarily
        equilibrate the triplet (A,B,C) as follows:
        = 'S':  perform equilibration (scaling);
        = 'N':  do not perform equilibration.
    n : int
        The order of the original state-space representation, i.e.
        the order of the matrix A. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        On entry, the leading N-by-N part of this array must
        contain the state dynamics matrix A.
    B : (n, m) array_like
        On entry, the leading N-by-M part of this array must
        contain the original input/state matrix B.
    C : (p, n) array_like
        On entry, the leading P-by-N part of this array must
        contain the original state/output matrix C.
    alpha : float, optional
        Specifies the alpha-stability boundary for the eigenvalues
        of the state dynamics matrix A. For a continuous-time
        system (dico = 'C'), alpha <= 0 is the boundary value for
        the real parts of eigenvalues, while for a discrete-time
        system (dico = 'D'), 0 <= alpha <= 1 represents the
        boundary value for the moduli of eigenvalues.
        The alpha-stability domain does not include the boundary.
        Default is None.
    nr : int, optional
        On entry with ordsel = 'F', nr is the desired order of the
        resulting reduced order system.  0 <= nr <= n.
        Default is None.
    tol : float, optional
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
        Default is `0.0`.
    ldwork : int, optional
        The length of the array dwork.
        ldwork >= max(1,n*(2*n+max(n,m,p)+5) + n*(n+1)/2).
        For optimum performance ldwork should be larger.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) array_like
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
    Br : (nr, m) array_like
        On exit, if info = 0, the leading nr-by-m part of this
        array contains the input/state matrix Br of the reduced
        order system.
    Cr : (p, nr) array_like
        On exit, if info = 0, the leading p-by-nr part of this
        array contains the state/output matrix Cr of the reduced
        order system.
    ns : int
        The dimension of the alpha-stable subsystem.
    hsv : (n, ) array_like
        If info = 0, the leading ns elements of hsv contain the
        Hankel singular values of the alpha-stable part of the
        original system ordered decreasingly.
        hsv(1) is the Hankel norm of the alpha-stable subsystem.

    Raises
    ------
    SlycotArithmeticError : e
        :info == 1:
            The computation of the ordered real Schur form of A failed
        :info == 2:
            The separation of the {alpha}-stable/unstable diagonal
            blocks failed because of very close eigenvalues
        :info == 3:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning : e
        :iwarn == 1:
            The selected order {nr} is greater
            than `nsmin`, the sum of the order of the
            {alpha}-unstable part and the order of a minimal
            realization of the {alpha}-stable part of the given
            system. The resulting `nr`  is set to `nsmin` = {Nr}
        :iwarn == 2:
            The selected order {nr} is less
            than the order of the {alpha}-unstable part of the
            given system. In this case `nr` is set equal to the
            order of the {alpha}-unstable part {Nr}.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'alpha',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'ns', 'hsv', 'tol',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n+max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    if alpha is None:
        alpha = {'C': 0, 'D': 1.}[dico]
    out = _wrapper.ab09md(dico, job, equil, ordsel,
                          n, m, p, nr, alpha, A, B, C, tol, ldwork)
    Nr, A, B, C, Ns, hsv = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09md.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], Ns, hsv

def ab09nd(dico,job,equil,n,m,p,A,B,C,D,alpha=None,nr=None,tol1=0,tol2=0,ldwork=None):
    """ nr,Ar,Br,Cr,Dr,ns,hsv = ab09nd(dico,job,equil,n,m,p,A,B,C,D,[alpha,nr,tol1,tol2,ldwork])

    To compute a reduced order model (Ar,Br,Cr,Dr) for an original
    state-space representation (A,B,C,D) by using either the
    square-root or the balancing-free square-root Singular
    Perturbation Approximation (SPA) model reduction method for the
    alpha-stable part of the system.

    Parameters
    ----------
    dico : {'C', 'D'}
        Specifies the type of the original system as follows:
        = 'C':  continuous-time system;
        = 'D':  discrete-time system.
    job : {'B', 'N'}
        Specifies the model reduction approach to be used
        as follows:
        = 'B':  use the square-root SPA method;
        = 'N':  use the balancing-free square-root SPA method.
    equil : {'S', 'N'}
        Specifies whether the user wishes to preliminarily
        equilibrate the triplet (A,B,C) as follows:
        = 'S':  perform equilibration (scaling);
        = 'N':  do not perform equilibration.
    n : int
        The order of the original state-space representation, i.e.
        the order of the matrix A. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        On entry, the leading n-by-n part of this array must
        contain the state dynamics matrix A.
    B : (n, m) array_like
        On entry, the leading n-by-m part of this array must
        contain the original input/state matrix B.
    C : (p, n) array_like
        On entry, the leading p-by-n part of this array must
        contain the original state/output matrix C.
    D : (p, m) array_like
        On entry, the leading p-by-m part of this array must
        contain the original input/output matrix D.
    alpha : float, optional
        Specifies the alpha-stability boundary for the eigenvalues
        of the state dynamics matrix A. For a continuous-time
        system (dico = 'C'), alpha <= 0 is the boundary value for
        the real parts of eigenvalues, while for a discrete-time
        system (dico = 'D'), 0 <= alpha <= 1 represents the
        boundary value for the moduli of eigenvalues.
        The alpha-stability domain does not include the boundary.
        Default is None.
    nr : int, optional
        nr is the desired order of
        the resulting reduced order system.  0 <= nr <= n.
        Default is None.
    tol1 : float, optional
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
        Default is `0.0`.
    tol2 : float, optional
        The tolerance for determining the order of a minimal
        realization of the alpha-stable part of the given system.
        The recommended value is tol2 = ns*eps*hnorm(As,Bs,Cs).
        This value is used by default if tol2 <= 0 on entry.
        If tol2 > 0, then tol2 <= tol1.
        Default is `0.0`.
    ldwork : int, optional
        The length of the array dwork.
        ldwork >= max(1,n*(2*n+max(n,m,p)+5) + n*(n+1)/2).
        For optimum performance ldwork should be larger.
        Default is None.

    Returns
    -------
    nr : int
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
    Ar : (nr, nr) ndarray
        the leading nr-by-nr part of this array contains the
        state dynamics matrix Ar of the reduced order system.
    Br : (nr, m) ndarray
        the leading nr-by-m part of this array contains the
        input/state matrix Br of the reduced order system.
    Cr : (p, nr) ndarray
        the leading p-by-nr part of this array contains the
        state/output matrix Cr of the reduced order system.
    Dr : (p, m) ndarray
        the leading p-by-m part of this array contains the
        input/output matrix Dr of the reduced order system.
    ns : int
        The dimension of the alpha-stable subsystem.
    hsv : (n, ) ndarray
        If info = 0, it contains the Hankel singular values of
        the original system ordered decreasingly. hsv(1) is the
        Hankel norm of the system.

    Raises
    ------
    SlycotArithmeticError
        :info == 1:
            The computation of the ordered real Schur form of A failed
        :info == 2:
            The separation of the {alpha}-stable/unstable diagonal
            blocks failed because of very close eigenvalues
        :info == 3:
            The computation of Hankel singular values failed

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
            The selected order {nr} is greater
            than `nsmin`, the sum of the order of the
            {alpha}-unstable part and the order of a minimal
            realization of the {alpha}-stable part of the given
            system. The resulting `nr`  is set to `nsmin` = {Nr}
        :iwarn == 2:
            The selected order {nr} is less
            than the order of the {alpha}-unstable part of the
            given system. In this case `nr` is set equal to the
            order of the {alpha}-unstable part {Nr}.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'alpha',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'D', 'ldc' + hidden, 'ns', 'hsv', 'tol1', 'tol2',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if ldwork is None:
        ldwork = max(1, n*(2*n+max(n, max(m, p))+5)+n*(n+1)/2)
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    if alpha is None:
        alpha = {'C': 0, 'D': 1.}[dico]
    out = _wrapper.ab09nd(dico, job, equil, ordsel,
                          n, m, p, nr, alpha, A, B, C, D, tol1, tol2, ldwork)
    #out = _wrapper.ab09nd_alt(dico, job, equil, ordsel, nr, alpha, A, B, C, D, tol1, tol2, ldwork)
    Nr, A, B, C, D, Ns, hsv = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09nd.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], D, Ns, hsv


def ab09hd(dico,job,equil,n,m,p,A,B,C,D,alpha=None,beta=0,nr=None,tol1=0,tol2=0,ldwork=None):
    """ nr,Ar,Br,Cr,Dr,ns,hsv = ab09hd(dico,job,equil,n,m,p,A,B,C,D,[alpha,nr,tol1,tol2,ldwork])

    Purpose

    To compute a reduced order model (Ar,Br,Cr,Dr) for an original
    state-space representation (A,B,C,D) by using the stochastic
    balancing approach in conjunction with the square-root or
    the balancing-free square-root Balance & Truncate (B&T)
    or Singular Perturbation Approximation (SPA) model reduction
    methods for the ALPHA-stable part of the system.

    Arguments

    Mode Parameters

    DICO    CHARACTER*1
            Specifies the type of the original system as follows:
            = 'C':  continuous-time system;
            = 'D':  discrete-time system.

    JOB     CHARACTER*1
            Specifies the model reduction approach to be used
            as follows:
            = 'B':  use the square-root Balance & Truncate method;
            = 'F':  use the balancing-free square-root
                    Balance & Truncate method;
            = 'S':  use the square-root Singular Perturbation
                    Approximation method;
            = 'P':  use the balancing-free square-root
                    Singular Perturbation Approximation method.

    EQUIL   CHARACTER*1
            Specifies whether the user wishes to preliminarily
            equilibrate the triplet (A,B,C) as follows:
            = 'S':  perform equilibration (scaling);
            = 'N':  do not perform equilibration.

    ORDSEL  CHARACTER*1
            Specifies the order selection method as follows:
            = 'F':  the resulting order NR is fixed;
            = 'A':  the resulting order NR is automatically determined
                    on basis of the given tolerance TOL1.

    Input/Output Parameters

    N       (input) INTEGER
            The order of the original state-space representation,
            i.e., the order of the matrix A.  N >= 0.

    M       (input) INTEGER
            The number of system inputs.  M >= 0.

    P       (input) INTEGER
            The number of system outputs.  P >= 0.
            P <= M if BETA = 0.

    NR      (input/output) INTEGER
            On entry with ORDSEL = 'F', NR is the desired order of the
            resulting reduced order system.  0 <= NR <= N.
            On exit, if INFO = 0, NR is the order of the resulting
            reduced order model. For a system with NU ALPHA-unstable
            eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N),
            NR is set as follows: if ORDSEL = 'F', NR is equal to
            NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order
            on entry, and NMIN is the order of a minimal realization
            of the ALPHA-stable part of the given system; NMIN is
            determined as the number of Hankel singular values greater
            than NS*EPS, where EPS is the machine precision
            (see LAPACK Library Routine DLAMCH);
            if ORDSEL = 'A', NR is the sum of NU and the number of
            Hankel singular values greater than MAX(TOL1,NS*EPS);
            NR can be further reduced to ensure that
            HSV(NR-NU) > HSV(NR+1-NU).

    ALPHA   (input) DOUBLE PRECISION
            Specifies the ALPHA-stability boundary for the eigenvalues
            of the state dynamics matrix A. For a continuous-time
            system (DICO = 'C'), ALPHA <= 0 is the boundary value for
            the real parts of eigenvalues, while for a discrete-time
            system (DICO = 'D'), 0 <= ALPHA <= 1 represents the
            boundary value for the moduli of eigenvalues.
            The ALPHA-stability domain does not include the boundary.

    BETA    (input) DOUBLE PRECISION
            BETA > 0 specifies the absolute/relative error weighting
            parameter. A large positive value of BETA favours the
            minimization of the absolute approximation error, while a
            small value of BETA is appropriate for the minimization
            of the relative error.
            BETA = 0 means a pure relative error method and can be
            used only if rank(D) = P.

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the leading N-by-N part of this array must
            contain the state dynamics matrix A.
            On exit, if INFO = 0, the leading NR-by-NR part of this
            array contains the state dynamics matrix Ar of the reduced
            order system.
            The resulting A has a block-diagonal form with two blocks.
            For a system with NU ALPHA-unstable eigenvalues and
            NS ALPHA-stable eigenvalues (NU+NS = N), the leading
            NU-by-NU block contains the unreduced part of A
            corresponding to ALPHA-unstable eigenvalues in an
            upper real Schur form.
            The trailing (NR+NS-N)-by-(NR+NS-N) block contains
            the reduced part of A corresponding to ALPHA-stable
            eigenvalues.

    LDA     INTEGER
            The leading dimension of array A.  LDA >= MAX(1,N).

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
            On entry, the leading N-by-M part of this array must
            contain the original input/state matrix B.
            On exit, if INFO = 0, the leading NR-by-M part of this
            array contains the input/state matrix Br of the reduced
            order system.

    LDB     INTEGER
            The leading dimension of array B.  LDB >= MAX(1,N).

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
            On entry, the leading P-by-N part of this array must
            contain the original state/output matrix C.
            On exit, if INFO = 0, the leading P-by-NR part of this
            array contains the state/output matrix Cr of the reduced
            order system.

    LDC     INTEGER
            The leading dimension of array C.  LDC >= MAX(1,P).

    D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
            On entry, the leading P-by-M part of this array must
            contain the original input/output matrix D.
            On exit, if INFO = 0, the leading P-by-M part of this
            array contains the input/output matrix Dr of the reduced
            order system.

    LDD     INTEGER
            The leading dimension of array D.  LDD >= MAX(1,P).

    NS      (output) INTEGER
            The dimension of the ALPHA-stable subsystem.

    HSV     (output) DOUBLE PRECISION array, dimension (N)
            If INFO = 0, the leading NS elements of HSV contain the
            Hankel singular values of the phase system corresponding
            to the ALPHA-stable part of the original system.
            The Hankel singular values are ordered decreasingly.

    Tolerances

    TOL1    DOUBLE PRECISION
            If ORDSEL = 'A', TOL1 contains the tolerance for
            determining the order of reduced system.
            For model reduction, the recommended value of TOL1 lies
            in the interval [0.00001,0.001].
            If TOL1 <= 0 on entry, the used default value is
            TOL1 = NS*EPS, where NS is the number of
            ALPHA-stable eigenvalues of A and EPS is the machine
            precision (see LAPACK Library Routine DLAMCH).
            If ORDSEL = 'F', the value of TOL1 is ignored.
            TOL1 < 1.

    TOL2    DOUBLE PRECISION
            The tolerance for determining the order of a minimal
            realization of the phase system (see METHOD) corresponding
            to the ALPHA-stable part of the given system.
            The recommended value is TOL2 = NS*EPS.
            This value is used by default if TOL2 <= 0 on entry.
            If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
            TOL2 < 1.

    Workspace

    IWORK   INTEGER array, dimension (MAX(1,2*N))
            On exit with INFO = 0, IWORK(1) contains the order of the
            minimal realization of the system.

    DWORK   DOUBLE PRECISION array, dimension (LDWORK)
            On exit, if INFO = 0, DWORK(1) returns the optimal value
            of LDWORK and DWORK(2) contains RCOND, the reciprocal
            condition number of the U11 matrix from the expression
            used to compute the solution X = U21*inv(U11) of the
            Riccati equation for spectral factorization.
            A small value RCOND indicates possible ill-conditioning
            of the respective Riccati equation.

    LDWORK  INTEGER
            The length of the array DWORK.
            LDWORK >= 2*N*N + MB*(N+P) + MAX( 2, N*(MAX(N,MB,P)+5),
                                    2*N*P+MAX(P*(MB+2),10*N*(N+1) ) ),
            where MB = M if BETA = 0 and MB = M+P if BETA > 0.
            For optimum performance LDWORK should be larger.

    BWORK   LOGICAL array, dimension 2*N

    Raises
    ------
    SlycotArithmeticError
        :info == 1:
            The computation of the ordered real Schur form of A failed
        :info == 2:
            The separation of the {alpha}-stable/unstable diagonal
            blocks failed because of very close eigenvalues
        :info == 3:
            The computation of Hankel singular values failed
        :info == 1:
            the computation of the ordered real Schur form of A
            failed;
        :info == 2:
            the reduction of the Hamiltonian matrix to real
            Schur form failed;
        :info == 3:
            the reordering of the real Schur form of the
            Hamiltonian matrix failed;
        :info == 4:
            the Hamiltonian matrix has less than N stable
            eigenvalues;
        :info == 5:
            the coefficient matrix U11 in the linear system
            X*U11 = U21 to determine X is singular to working
            precision;
        :info == 6:
            BETA = 0 and D has not a maximal row rank;
        :info == 7:
            the computation of Hankel singular values failed;
        :info == 8:
            the separation of the ALPHA-stable/unstable diagonal
            blocks failed because of very close eigenvalues;
        :info == 9:
            the resulting order of reduced stable part is less
            than the number of unstable zeros of the stable
            part.

    Warns
    -----
    SlycotResultWarning
        :iwarn == 1:
            with ORDSEL = 'F', the selected order NR is greater
            than NSMIN, the sum of the order of the
            ALPHA-unstable part and the order of a minimal
            realization of the ALPHA-stable part of the given
            system; in this case, the resulting NR is set equal
            to NSMIN;
        :iwarn == 2:
            with ORDSEL = 'F', the selected order NR corresponds
            to repeated singular values for the ALPHA-stable
            part, which are neither all included nor all
            excluded from the reduced model; in this case, the
            resulting NR is automatically decreased to exclude
            all repeated singular values;
        :iwarn == 3:
              with ORDSEL = 'F', the selected order NR is less
              than the order of the ALPHA-unstable part of the
              given system; in this case NR is set equal to the
              order of the ALPHA-unstable part.


    Method

    Let be the following linear system

        d[x(t)] = Ax(t) + Bu(t)
        y(t)    = Cx(t) + Du(t),                      (1)

    where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
    for a discrete-time system. The subroutine AB09HD determines for
    the given system (1), the matrices of a reduced order system

        d[z(t)] = Ar*z(t) + Br*u(t)
        yr(t)   = Cr*z(t) + Dr*u(t),                  (2)

    such that

        INFNORM[inv(conj(W))*(G-Gr)] <=
                        (1+HSV(NR+NS-N+1)) / (1-HSV(NR+NS-N+1)) + ...
                        + (1+HSV(NS)) / (1-HSV(NS)) - 1,

    where G and Gr are transfer-function matrices of the systems
    (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, W is the right, minimum
    phase spectral factor satisfying

        G1*conj(G1) = conj(W)* W,                      (3)

    G1 is the NS-order ALPHA-stable part of G, and INFNORM(G) is the
    infinity-norm of G. HSV(1), ... , HSV(NS) are the Hankel-singular
    values of the stable part of the phase system (Ap,Bp,Cp)
    with the transfer-function matrix

        P = inv(conj(W))*G1.

    If BETA > 0, then the model reduction is performed on [G BETA*I]
    instead of G. This is the recommended approach to be used when D
    has not a maximal row rank or when a certain balance between
    relative and absolute approximation errors is desired. For
    increasingly large values of BETA, the obtained reduced system
    assymptotically approaches that computed by using the
    Balance & Truncate or Singular Perturbation Approximation methods.

    Note: conj(G)  denotes either G'(-s) for a continuous-time system
            or G'(1/z) for a discrete-time system.
            inv(G) is the inverse of G.

    The following procedure is used to reduce a given G:

    1) Decompose additively G as

        G = G1 + G2,

        such that G1 = (As,Bs,Cs,D) has only ALPHA-stable poles and
        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles.

    2) Determine G1r, a reduced order approximation of the
        ALPHA-stable part G1 using the balancing stochastic method
        in conjunction with either the B&T [1,2] or SPA methods [3].

    3) Assemble the reduced model Gr as

            Gr = G1r + G2.

    Note: The employed stochastic truncation algorithm [2,3] has the
    property that right half plane zeros of G1 remain as right half
    plane zeros of G1r. Thus, the order can not be chosen smaller than
    the sum of the number of unstable poles of G and the number of
    unstable zeros of G1.

    The reduction of the ALPHA-stable part G1 is done as follows.

    If JOB = 'B', the square-root stochastic Balance & Truncate
    method of [1] is used.
    For an ALPHA-stable continuous-time system (DICO = 'C'),
    the resulting reduced model is stochastically balanced.

    If JOB = 'F', the balancing-free square-root version of the
    stochastic Balance & Truncate method [1] is used to reduce
    the ALPHA-stable part G1.

    If JOB = 'S', the stochastic balancing method is used to reduce
    the ALPHA-stable part G1, in conjunction with the square-root
    version of the Singular Perturbation Approximation method [3,4].

    If JOB = 'P', the stochastic balancing method is used to reduce
    the ALPHA-stable part G1, in conjunction with the balancing-free
    square-root version of the Singular Perturbation Approximation
    method [3,4].

    References

    [1] Varga A. and Fasol K.H.
        A new square-root balancing-free stochastic truncation model
        reduction algorithm.
        Proc. 12th IFAC World Congress, Sydney, 1993.

    [2] Safonov M. G. and Chiang R. Y.
        Model reduction for robust control: a Schur relative error
        method.
        Int. J. Adapt. Contr. Sign. Proc., vol. 2, pp. 259-272, 1988.

    [3] Green M. and Anderson B. D. O.
        Generalized balanced stochastic truncation.
        Proc. 29-th CDC, Honolulu, Hawaii, pp. 476-481, 1990.

    [4] Varga A.
        Balancing-free square-root algorithm for computing
        singular perturbation approximations.
        Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
        Vol. 2, pp. 1062-1065.

    Numerical Aspects

    The implemented methods rely on accuracy enhancing square-root or
    balancing-free square-root techniques. The effectiveness of the
    accuracy enhancing technique depends on the accuracy of the
    solution of a Riccati equation. An ill-conditioned Riccati
    solution typically results when [D BETA*I] is nearly
    rank deficient.
    The algorithm requires about 100N^3  floating point operations.

    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['dico', 'job', 'equil', 'ordsel', 'n', 'm', 'p', 'nr', 'alpha',
                'beta',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'D', 'ldc' + hidden, 'ns', 'hsv', 'tol1', 'tol2',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'iwarn', 'info']
    if nr is None:
        ordsel = 'A'
        nr = 0  # order will be computed by the routine
    else:
        ordsel = 'F'
    if alpha is None:
        alpha = {'C': 0, 'D': 1.}[dico]
    out = _wrapper.ab09hd(
        dico, job, equil, ordsel,
        n, m, p, nr, alpha, beta, A, B, C, D, tol1, tol2,
    )
    Nr, A, B, C, D, Ns, hsv = out[:-2]
    raise_if_slycot_error(out[-2:], arg_list, ab09hd.__doc__, locals())
    return Nr, A[:Nr, :Nr], B[:Nr, :], C[:, :Nr], D, Ns, hsv


def ab13bd(dico, jobn, n, m, p, A, B, C, D, tol = 0.0):
    """norm = ab13bd(dico, jobn, n, m, p, A, B, C, D, [tol])

    To compute the H2 or L2 norm of the transfer-function matrix G
    of the system (A,B,C,D). G must not have poles on the imaginary
    axis, for a continuous-time system, or on the unit circle, for
    a discrete-time system. If the H2-norm is computed, the system
    must be stable.

    Parameters
    ----------
    dico : {'D', 'C'}
        Indicate whether the system is discrete 'D' or continuous 'C'.
    jobn : {'H', 'L'}
        H2-norm 'H' or L2-norm 'L' to be computed.
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) ndarray
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system.
    B : (n, m) ndarray
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) ndarray
        The leading p-by-n part of this array must contain the state/output
        matrix C of the system.
    D : (p, m) ndarray
        The leading p-by-m part of this array must contain the direct
        transmission matrix D of the system.
    tol : float, optional
        The absolute tolerance level below which the elements of
        B are considered zero (used for controllability tests).
        If the user sets tol <= 0, then an implicitly computed,
        default tolerance, defined by  toldef = n*eps*norm(B),
        is used instead, where eps is the machine precision
        (see LAPACK Library routine DLAMCH) and norm(B) denotes
        the 1-norm of B.

    Returns
    -------
    norm:  H2 or L2 norm of the system (A,B,C,D)

    Raises
    ------
    SlycotArithmeticError
        :info == 1:
            The reduction of A to a real Schur form failed
        :info == 2:
            A failure was detected during the reordering of the
            real Schur form of A, or in the iterative process for
            reordering the eigenvalues of `` Z'*(A + B*F)*Z`` along the
            diagonal (see SLICOT routine SB08DD)
        :info == 3 and dico == 'C':
            The matrix A has a controllable eigenvalue on the imaginary axis
        :info == 3 and dico == 'D':
            The matrix A has a controllable eigenvalue on the unit circle
        :info == 4:
            The solution of Lyapunov equation failed because the
            equation is singular
        :info == 5:
            D is a nonzero matrix
        :info == 6:
            The system is unstable

    Warns
    -----
    SlycotResultWarning
        :iwarn > 0:
            {iwarn} violations of the numerical stability condition
            occured during the assignment of eigenvalues in
            computing the right coprime factorization with inner
            denominator of `G` (see the SLICOT subroutine SB08DD).
    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ('dico', 'jobn', 'n', 'm', 'p',
                'A', 'lda' + hidden, 'B', 'ldb' + hidden, 'C', 'ldc' + hidden,
                'D', 'ldd' + hidden, 'nq' + hidden,'tol', 'dwork' + hidden,
                'ldwork' + hidden, 'iwarn', 'info')
    
    out = _wrapper.ab13bd(dico, jobn, n, m, p, A, B, C, D, tol)

    raise_if_slycot_error(out[-2:], arg_list, ab13bd.__doc__, locals())
    return out[0]

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

    Parameters
    ----------
    dico : {'D', 'C'}
            Indicate whether the system is discrete 'D' or continuous 'C'.
    jobe : {'G', 'I'}
            Specifies whether E is a general square or an identity
            matrix, as follows:
            = 'G':  E is a general square matrix;
            = 'I':  E is the identity matrix.
    equil : {'S', 'N'}
            Specifies whether the user wishes to preliminarily
            equilibrate the system (A,E,B,C) or (A,B,C), as follows:
            = 'S':  perform equilibration (scaling);
            = 'N':  do not perform equilibration.
    jobd : {'D', 'Z'}
            Specifies whether or not a non-zero matrix D appears in
            the given state space model:
            = 'D':  D is present;
            = 'Z':  D is assumed a zero matrix.
    n : int
        The number of state variables. n >= 0.
    m : int
        The number of system inputs. m >= 0.
    p : int
        The number of system outputs. p >= 0.
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the state
        dynamics matrix A of the system.
    E : (n, n) array_like
        If jobe = 'G', the leading N-by-N part of this array must
        contain the descriptor matrix E of the system.
        If jobe = 'I', then E is assumed to be the identity
        matrix and is not referenced.
    B : (n, m) array_like
        The leading n-by-m part of this array must contain the input/state
        matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must contain the state/output
        matrix C of the system.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the direct
        transmission matrix D of the system.
    tol : float
        Tolerance used to set the accuracy in determining the norm.
        0 <= tol < 1. Default tol=1e-10.

    Returns
    -------
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

    Raises
    ------
    SlycotArithmeticError
        :info = 1:
            The matrix E is (numerically) singular
        :info = 2:
            The (periodic) QR (or QZ) algorithm for computing
            eigenvalues did not converge
        :info = 3:
            The SVD algorithm for computing singular values did
            not converge
        :info = 4:
            The tolerance is too small and the algorithm did not converge
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ('dico', 'jobe', 'equil', 'jobd', 'n', 'm', 'p',
                'fpeak' + hidden,
                'A', 'lda' + hidden, 'E', 'lde' + hidden, 'B', 'ldb' + hidden,
                'C', 'ldc' + hidden, 'D', 'ldd' + hidden,
                'gpeak' + hidden, 'tol', 'iwork' + hidden, 'dwork' + hidden,
                'ldwork' + hidden, 'cwork' + hidden, 'lcwork' + hidden,
                'info' + hidden)
    if dico != 'C' and dico != 'D':
        raise SlycotParameterError('dico must be "C" or "D"', -1)
    if jobe != 'G' and jobe != 'I':
        raise SlycotParameterError('jobe must be "G" or "I"', -2)
    if equil != 'S' and equil != 'N':
        raise SlycotParameterError('equil must be "S" or "N"', -3)
    if jobd != 'D' and jobd != 'Z':
        raise SlycotParameterError('jobd must be "D" or "Z"', -4)
    out = _wrapper.ab13dd(dico, jobe, equil, jobd,
                          n, m, p, [0.0, 1.0], A, E, B, C, D, tol)
    raise_if_slycot_error(out[-1], arg_list, ab13dd.__doc__)

    fpeak = out[0][0] if out[0][1] > 0 else float('inf')
    gpeak = out[1][0] if out[1][1] > 0 else float('inf')
    return gpeak, fpeak



def ab13ed(n, A, tol = 9.0):
    """low, high = ab13ed(n, A, [tol])

    To estimate beta(A), the 2-norm distance from a real matrix A to
    the nearest complex matrix with an eigenvalue on the imaginary
    axis. The estimate is given as

         ``low <= beta(A) <= high,``

    where either

         ``(1 + tol) * low >= high,``

    or

         ``low = 0   and   high = delta,``

    and delta is a small number approximately equal to the square root
    of machine precision times the Frobenius norm (Euclidean norm)
    of A. If A is stable in the sense that all eigenvalues of A lie
    in the open left half complex plane, then beta(A) is the distance
    to the nearest unstable complex matrix, i.e., the complex
    stability radius.

    Parameters
    ----------
    n : int
        The order of the matrix A.  ``n >= 0.``
    A : (n, n) array_like
        The leading n-by-n part of this array must contain the matrix A.
    tol : float, optional
        Specifies the accuracy with which low and high approximate
        beta(A). If the user sets tol to be less than sqrt(eps),
        where eps is the machine precision (see LAPACK Library
        Routine DLAMCH), then the tolerance is taken to be
        sqrt(eps).
        The recommended value is tol = 9, which gives an estimate
        of beta(A) correct to within an order of magnitude.

    Returns
    -------
    low : float
          A lower bound for beta(A).
    high : float
           An upper bound for beta(A).

    Raises
    ------
    SlycotArithmeticError
        :info = 1:
            The QR algorithm fails to converge
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'A', 'lda' + hidden, 'low' + hidden, 'high' + hidden, 'tol',
                'dwork' + hidden, 'ldwork' + hidden, 'info' + hidden]
    out = _wrapper.ab13ed(n, A, tol)
    raise_if_slycot_error(out[-1], arg_list, ab13ed.__doc__)
    return out[:-1]

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

    Parameters
    ----------
    n : int
        The order of the matrix A.  n >= 0.
    A : (n, n), array_like
        The leading n-by-n part of this array must contain the matrix A.
    tol : float, optional
        Specifies the accuracy with which beta(A) is to be
        calculated. (See the Numerical Aspects section below.)
        If the user sets tol to be less than eps, where eps is the
        machine precision (see LAPACK Library Routine DLAMCH),
        then the tolerance is taken to be eps.

    Returns
    -------
    beta : float
            The computed value of beta(A), which actually is an upper
            bound.
    omega : float
            The value of w such that the smallest singular value of
            (A - jwI) equals beta(A).

    Raises
    ------
    SlycotArithmeticError
        :info = 2:
            Either the QR or SVD algorithm fails to converge

    Warns
    -----
    SlycotResultWarning
        :info = 1:
            Failed to compute beta(A) within the specified tolerance.
            Nevertheless, the returned value is an upper bound on beta(A);

    Notes
    -----
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
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'A', 'lda' + hidden, 'beta' + hidden, 'omega' + hidden, 'tol',
                'dwork' + hidden, 'ldwork' + hidden, 'cwork' + hidden,
                'lcwork' + hidden, 'info' + hidden]
    out = _wrapper.ab13fd(n, A, tol)
    raise_if_slycot_error(out[-1], arg_list, ab13fd.__doc__)
    return out[0], out[1]


def ab13md(Z, nblock, itype, x=None):
    """mubound, d, g, xout = ab13md(Z, nblock, itype, [x])

    Find an upper bound for the structured singular value of complex
    matrix Z and given block diagonal structure.

    Parameters
    ----------
    Z : (n, n) array_like
        Matrix to find structured singular value upper bound of
    nblock : (m, ) array_like
        The size of the block diagonals of the uncertainty structure;
        i.e., nblock(i)=p means that the ith block is pxp.
    itype : (m, ) array_like
        The type of each block diagonal uncertainty defined in nblock.
        itype(i)==1 means that the ith block is real, while itype(i)==2
        means the the ith block is complex.  Real blocks must be 1x1,
        i.e., if itype(i)==1, nblock(i) must be 1.
    x : (q, ) array_like, optional
        If not None, must be the output of a previous call to ab13md.
        The previous call must have been with the same values of n,
        nblock, and itype; and the previous call's Z should be "close"
        to the current call's Z.
        q is determined by the block structure; see SLICOT AB13MD for
        details. Default is None.

    Returns
    -------
    mubound : non-negative real scalar
        Upper bound on structure singular value for given arguments
    d, g : (n, ) ndarray
        Real arrays such that if D=np.diag(g), G=np.diag(G), and ZH = Z.T.conj(), then
          ZH @ D**2 @ Z + 1j * (G@Z - ZH@G) - mu**2 * D**2
        will be negative semi-definite.
    xout : (q, ) ndarray
        For use as ``x`` argument in subsequent call to ``ab13md``.

    For scalar Z and real uncertainty (ntype=1, itype=1), returns 0
    instead of abs(Z).

    Raises
    ------
    SlycotArithmeticError
        :info = 1: Block sizes must be positive
        :info = 2: Block sizes must sum to n
        :info = 3: Real blocks must be of size 1
        :info = 4: Block types must be 1 or 2
        :info = 5: Error in linear equation solution
        :info = 6: Error in eigenvalue or singular value computation

    Notes
    -----
    This wraps SLICOT routine AB13MD, which implements the upper bound
    of [1].

    References
    ----------
    .. [1] Fan, M.K.H., Tits, A.L., and Doyle, J.C., "Robustness in
       the presence of mixed parametric uncertainty and unmodeled
       dynamics," IEEE Trans. Automatic Control, vol. AC-36, 1991,
       pp. 25-38.

    """
    hidden = ' (hidden by the wrapper)'

    arg_list = ['fact', 'n' + hidden, 'z', 'ldz' + hidden, 'm' + hidden,
                'nblock', 'itype', 'x', 'bound', 'd', 'g',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork' + hidden,
                'zwork' + hidden, 'lzwork' + hidden, 'info']

    # prepare the "x" input and output

    # x, in SLICOT, needs to be length m+mr-1. m is the length of
    # nblock (and itype), and mr is the number of real blocks.

    # In analysis.pyf x is specified as length 2*m-1, since I couldn't
    # figure out how to express the m+mr-1 constraint there.

    # The code here is to arrange for the user-visible part of x,
    # which is length m+mr-1, to be packed into the 2*m-1-length array
    # to pass the SLICOT routine.

    m = len(nblock)
    mr = np.count_nonzero(1==itype)

    if x is None:
        fact='N'
        x = np.empty(2*m-1)
    else:
        fact='F'
        if len(x) != m+mr-1:
            raise ValueError(f'Require len(x)==m+mr-1, but'
                             + f' len(x)={len(x)}, m={m}, mr={mr}')
        x = np.concatenate([x,np.zeros(2*m-1-len(x))])

    x, bound, d, g, info = _wrapper.ab13md(fact, Z, nblock, itype, x)

    raise_if_slycot_error(info, arg_list, ab13md.__doc__)

    return bound, d, g, x[:m+mr-1]


def ag08bd(l,n,m,p,A,E,B,C,D,equil='N',tol=0.0,ldwork=None):
    """ Af,Ef,nrank,niz,infz,kronr,infe,kronl = ag08bd(l,n,m,p,A,E,B,C,D,[equil,tol,ldwork])

    To extract from the system pencil

                        ( A-lambda*E B )
            S(lambda) = (              )
                        (      C     D )

    a regular pencil Af-lambda*Ef which has the finite Smith zeros of
    S(lambda) as generalized eigenvalues. The routine also computes
    the orders of the infinite Smith zeros and determines the singular
    and infinite Kronecker structure of system pencil, i.e., the right
    and left Kronecker indices, and the multiplicities of infinite
    eigenvalues.

    Parameters
    ----------
    l : int
        The number of rows of matrices A, B, and E.  l >= 0.
    n : int
        The number of columns of matrices A, E, and C.  n >= 0.
    m : int
        The number of columns of matrix B.  m >= 0.
    p : int
        The number of rows of matrix C.  p >= 0.
    A : (l, n) array_like
        The leading l-by-n part of this array must
        contain the state dynamics matrix A of the system.
    E : (l, n) array_like
        The leading l-by-n part of this array must
        contain the descriptor matrix E of the system.
    B : (l, m) array_like
        The leading l-by-m part of this array must
        contain the input/state matrix B of the system.
    C : (p, n) array_like
        The leading p-by-n part of this array must
        contain the state/output matrix C of the system.
    D : (p, m) array_like
        The leading p-by-m part of this array must contain the
        direct transmission matrix D of the system.
    equil : {'S', 'N'}, optional
        Specifies whether the user wishes to balance the system
        matrix as follows:
        = 'S':  Perform balancing (scaling);
        = 'N':  Do not perform balancing.
        Default is `N`
    tol : float, optional
        A tolerance used in rank decisions to determine the
        effective rank, which is defined as the order of the
        largest leading (or trailing) triangular submatrix in the
        QR (or RQ) factorization with column (or row) pivoting
        whose estimated condition number is less than 1/TOL.
        If the user sets TOL <= 0, then default tolerances are
        used instead, as follows: TOLDEF = L*N*EPS in TG01FD
        (to determine the rank of E) and TOLDEF = (L+P)*(N+M)*EPS
        in the rest, where EPS is the machine precision
        (see LAPACK Library routine DLAMCH).  TOL < 1.
        Default is 0.
    ldwork : int, optional
        The length of the cache array.
        ldwork >= max( 4*(l,n), ldw ), if equil = 'S',
        ldwork >= ldw,                 if equil = 'N', where
        ldw = max(l+p,m+n)*(m+n) + max(1,5*max(l+p,m+n)).
        For optimum performance ldwork should be larger.
        Default is None.

    Returns
    -------
    Af : (nfz, nfz) ndarray
        the leading NFZ-by-NFZ part of this array
        contains the matrix Af of the reduced pencil.
    Ef : (nfz, nfz) ndarray
        the leading NFZ-by-NFZ part of this array
        contains the matrix Ef of the reduced pencil.
    nrank : int
        The normal rank of the system pencil.
    niz : int
        The number of infinite zeros.
    infz : (n+1, ) ndarray
        Contains information on the infinite elementary divisors.
    kronr : (n+m+1, ) ndarray
        Contains the right Kronecker (column) indices.
    infe : (1+min(l+p,n+m), ) ndarray
        Contains the multiplicities of infinite eigenvalues.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['equil', 'l', 'n', 'm', 'p',
                'A', 'lda' + hidden, 'E', 'lde' + hidden, 'B', 'ldb' + hidden,
                'C', 'ldc' + hidden, 'D', 'ldd' + hidden,
                'nfz', 'nrank', 'niz', 'dinfz', 'nkror', 'ninfe', 'nkrol',
                'infz', 'kronr', 'infe', 'kronl', 'tol',
                'iwork' + hidden, 'dwork' + hidden, 'ldwork', 'info']

    if ldwork is None:
        ldw = max(l+p, m+n)*(m+n) + max(1, 5*max(l+p, m+n))
        if equil == 'S':
            ldwork = max(4*(l+n), ldw)
        else:  # equil == 'N'
            ldwork = ldw

    out = _wrapper.ag08bd(equil, l, n, m, p, A, E, B, C, D, tol, ldwork)
    [Af, Ef, nfz, nrank, niz,
     dinfz, nkror, ninfe, nkrol, infz, kronr, infe, kronl, info] = out
    raise_if_slycot_error(info, arg_list)

    return (Af[:nfz, :nfz], Ef[:nfz, :nfz], nrank, niz,
            infz[:dinfz], kronr[:nkror], infe[:ninfe], kronl[:nkrol])
