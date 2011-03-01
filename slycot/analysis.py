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

from slycot import _wrapper

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
    
def ab09ad(dico,job,equil,n,m,p,a,b,c,nr=None,tol=0,ldwork=None):
    """ nr,A_r,B_r,C_r,hsv = ab09ad(dico,job,equil,n,m,p,nr,A,B,C,[nr,tol,ldwork,overwrite_a,overwrite_b,overwrite_c]) """
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
    out = _wrapper.ab09ad(dico,job,equil,ordsel,n,m,p,nr,a,b,c,tol,ldwork)
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
    
# to be replaced by python wrappers
