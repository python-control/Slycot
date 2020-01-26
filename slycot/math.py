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
import warnings

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
        warnings.warn('entry P(x) is the zero polynomial.')
    if out[-1] == 2:
        warnings.warn('P(x) may have zeros very close to stability boundary.')
    if out[-2] > 0:
        warnings.warn('The degree of P(x) has been reduced to %i' %(dp-out[-2]))
    return out[:-2]


def mb03vd(n, ilo, ihi, A):
    """ HQ, Tau = mb03vd(n, ilo, ihi, A)

    To reduce a product of p real general matrices A = A_1*A_2*...*A_p
    to upper Hessenberg form, H = H_1*H_2*...*H_p, where H_1 is
    upper Hessenberg, and H_2, ..., H_p are upper triangular, by using
    orthogonal similarity transformations on A,

            Q_1' * A_1 * Q_2 = H_1,
            Q_2' * A_2 * Q_3 = H_2,
                   ...
            Q_p' * A_p * Q_1 = H_p.

    Parameters
    ----------

    n : int
            The order of the square matrices A_1, A_2, ..., A_p.
            n >= 0.

    ilo, ihi : int
            It is assumed that all matrices A_j, j = 2, ..., p, are
            already upper triangular in rows and columns [:ilo] and
            [ihi:n], and A_1 is upper Hessenberg in rows and columns
            [:ilo] and [ihi:n], with A_1[ilo-1,ilo] = 0 (unless
            ilo = 1), and A_1[ihi,ihi-1] = 0 (unless ihi = n).
            If this is not the case, ilo and ihi should be set to 1
            and n, respectively.
            1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.

    A : ndarray
            A[:n,:n,:p] must contain the matrices of factors to be reduced;
            specifically, A[:,:,j-1] must contain A_j, j = 1, ..., p.


    Returns
    -------

    HQ : ndarray
            The upper triangle and the first
            subdiagonal of HQ[:n,:n,0] contain the upper Hessenberg
            matrix H_1, and the elements below the first subdiagonal,
            with the first column of the array Tau represent the
            orthogonal matrix Q_1 as a product of elementary
            reflectors. See FURTHER COMMENTS.
            For j > 1, the upper triangle of HQ[:n,_n,j-1]
            contains the upper triangular matrix H_j, and the elements
            below the diagonal, with the j-th column of the array TAU
            represent the orthogonal matrix Q_j as a product of
            elementary reflectors. See FURTHER COMMENTS.

    Tau : ndarray
            The leading n-1 elements in the j-th column contain the
            scalar factors of the elementary reflectors used to form
            the matrix Q_j, j = 1, ..., p. See FURTHER COMMENTS.

    Raises
    ------

    ValueError : e
            e.info contains information about the exact type of exception

    Further Comments
    ----------------

    Each matrix Q_j is represented as a product of (ihi-ilo)
    elementary reflectors,

       Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1).

    Each H_j(i), i = ilo, ..., ihi-1, has the form

       H_j(i) = I - tau_j * v_j * v_j',

    where tau_j is a real scalar, and v_j is a real vector with
    v_j(1:i) = 0, v_j(i+1) = 1 and v_j(ihi+1:n) = 0; v_j(i+2:ihi)
    is stored on exit in A_j(i+2:ihi,i), and tau_j in TAU(i,j).

    The contents of A_1 are illustrated by the following example
    for n = 7, ilo = 2, and ihi = 6:

    on entry                         on exit

    ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a )
    ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   a )
    ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h )
    ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a )

    where a denotes an element of the original matrix A_1, h denotes
    a modified element of the upper Hessenberg matrix H_1, and vi
    denotes an element of the vector defining H_1(i).

    The contents of A_j, j > 1, are illustrated by the following
    example for n = 7, ilo = 2, and ihi = 6:

    on entry                         on exit

    ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a )
    ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h )
    ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  v5  h   h )
    ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a )

    where a denotes an element of the original matrix A_j, h denotes
    a modified element of the upper triangular matrix H_j, and vi
    denotes an element of the vector defining H_j(i). (The element
    (1,2) in A_p is also unchanged for this example.)

    Note that for P = 1, the LAPACK Library routine DGEHRD could be
    more efficient on some computer architectures than this routine
    (a BLAS 2 version).
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'p' + hidden,
                'ilo', 'ihi', 'a',
                'lda1' + hidden, 'lda2' + hidden, 'tau',
                'ldtau' + hidden, 'dwork' + hidden, 'info']

    HQ, Tau, info = _wrapper.mb03vd(n, ilo, ihi, A)

    if info != 0:
        e = ValueError(
                "Argument '{}' had an illegal value".format(arg_list[-info-1]))
        e.info = info
        raise e
    return (HQ, Tau)


def mb03vy(n, ilo, ihi, A, Tau, ldwork=None):
    """ Q = mb03vy(n, ilo, ihi, A, Tau, [ldwork])

    To generate the real orthogonal matrices Q_1, Q_2, ..., Q_p,
    which are defined as the product of ihi-ilo elementary reflectors
    of order n, as returned by SLICOT Library routine MB03VD:

        Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1).

    Parameters
    ----------

    n : int
            The order of the matrices Q_1, Q_2, ..., Q_p.  N >= 0.

    ilo, ihi : int
            The values of the indices ilo and ihi, respectively, used
            in the previous call of the SLICOT Library routine MB03VD.
            1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.

    A : ndarray
            A[:n,:n,j-1] must contain the vectors which define the
            elementary reflectors used for reducing A_j, as returned
            by SLICOT Library routine MB03VD, j = 1, ..., p.

    Tau : ndarray
            The leading N-1 elements in the j-th column must contain
            the scalar factors of the elementary reflectors used to
            form the matrix Q_j, as returned by SLICOT Library routine
            MB03VD.

    ldwork : int, optional
            The length of the array DWORK.  LDWORK >= MAX(1,N).
            For optimum performance LDWORK should be larger.


    Returns
    -------

    Q : ndarray
           Q[:n,:n,j-1] contains the
             N-by-N orthogonal matrix Q_j, j = 1, ..., p.

    Raises
    ------

    ValueError :
            e.info contains the number of the argument that was invalid

    """

    hidden = ' (hidden by the wrapper)'
    arg_list = ['n', 'p' + hidden,
                'ilo', 'ihi', 'a',
                'lda1' + hidden, 'lda2' + hidden, 'tau',
                'ldtau' + hidden, 'dwork' + hidden, 'info']

    if not ldwork:
        ldwork = max(1, 2 * n)

    Q, info = _wrapper.mb03vy(n, ilo, ihi, A, Tau, ldwork)

    if info != 0:
        e = ValueError(
                "Argument '{}' had an illegal value".format(arg_list[-info-1]))
        e.info = info
        raise e

    return Q


def mb03wd(job, compz, n, ilo, ihi, iloz, ihiz, H, Q, ldwork=None):
    """ T, Z, Wr = mb03wd(job, compz, n, ilo, ihi, iloz, ihiz, H, Q, [ldwork])

    To compute the Schur decomposition and the eigenvalues of a
    product of matrices, H = H_1*H_2*...*H_p, with H_1 an upper
    Hessenberg matrix and H_2, ..., H_p upper triangular matrices,
    without evaluating the product. Specifically, the matrices Z_i
    are computed, such that

            Z_1' * H_1 * Z_2 = T_1,
            Z_2' * H_2 * Z_3 = T_2,
                   ...
            Z_p' * H_p * Z_1 = T_p,

    where T_1 is in real Schur form, and T_2, ..., T_p are upper
    triangular.

    The routine works primarily with the Hessenberg and triangular
    submatrices in rows and columns ilo to ihi, but optionally applies
    the transformations to all the rows and columns of the matrices
    H_i, i = 1,...,p. The transformations can be optionally
    accumulated.

    Parameters
    ----------

    job : {'E', 'S'}
            Indicates whether the user wishes to compute the full
            Schur form or the eigenvalues only, as follows:
            = 'E':  Compute the eigenvalues only;
            = 'S':  Compute the factors T_1, ..., T_p of the full
                    Schur form, T = T_1*T_2*...*T_p.

   compz : {'N', 'I', 'V'}
            Indicates whether or not the user wishes to accumulate
            the matrices Z_1, ..., Z_p, as follows:
            = 'N':  The matrices Z_1, ..., Z_p are not required;
            = 'I':  Z_i is initialized to the unit matrix and the
                    orthogonal transformation matrix Z_i is returned,
                    i = 1, ..., p;
            = 'V':  Z_i must contain an orthogonal matrix Q_i on
                    entry, and the product Q_i*Z_i is returned,
                    i = 1, ..., p.

    n : int
            The order of the matrix H.  n >= 0

    ilo, ihi : int
            It is assumed that all matrices H_j, j = 2, ..., p, are
            already upper triangular in rows and columns [:ilo] and
            [ihi+1:n], and H_1 is upper quasi-triangular in rows and
            columns [:ilo] and [ihi+1:n], with H_1[ilo-1,ilo] = 0
            (unless ilo = 1), and H_1[ihi,ihi-1] = 0 (unless ihi = n).
            The routine works primarily with the Hessenberg submatrix
            in rows and columns ilo to ihi, but applies the
            transformations to all the rows and columns of the
            matrices H_i, i = 1,...,p, if JOB = 'S'.
            1 <= ilo <= max(1,n); min(ilo,n) <= ihi <= n.

    iloz, ihiz : int
            Specify the rows of Z to which the transformations must be
            applied if compz = 'I' or compz = 'V'.
            1 <= iloz <= ilo; ihi <= ihiz <= n.

    H : ndarray
            H[:n,:n,0] must contain the upper Hessenberg matrix H_1 and
            H[:n,:n,j-1] for j > 1 must contain the upper triangular matrix
            H_j, j = 2, ..., p.

    Q : ndarray
            If compz = 'V', Q[:n,:n,:p] must contain the current matrix Q of
            transformations accumulated by SLICOT Library routine
            MB03VY.
            If compz = 'I', Q is ignored

    ldwork : int, optinal
            The length of the cache array. The default value is
            ihi-ilo+p-1



    Returns
    -------

    T : ndarray
            If JOB = 'S', T[:n,:n,0] s upper quasi-triangular in rows
            and columns [ilo-1:ihi], with any 2-by-2 diagonal blocks
            corresponding to a pair of complex conjugated eigenvalues, and
            T[:n,:n,j-1] for j > 1 contains the resulting upper
            triangular matrix T_j.
            If job = 'E', T is  None

    Z : ndarray
            If compz = 'V', or compz = 'I', the leading
            N-by-N-by-P part of this array contains the transformation
            matrices which produced the Schur form; the
            transformations are applied only to the submatrices
            Z[iloz-1:ihiz,ilo-1:ihi,j-1], j = 1, ..., P.
            If compz = 'N', Z is None


    W : ndarray (dtype=complex)
            The computed eigenvalues ilo to ihi. If two eigenvalues
            are computed as a complex conjugate pair, they are stored
            in consecutive elements of Wr say the i-th and
            (i+1)th, with imag(W][i]) > 0 and imag(W[i+1]) < 0.
            If JOB = 'S', the eigenvalues are stored in the same order
            as on the diagonal of the Schur form returned in H.

    Raises
    ------

    ValueError : e
            e.info contains information about the exact type of exception

    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['job', 'compz', 'n', 'p' + hidden,
                'ilo', 'ihi', 'iloz', 'ihiz',
                'h', 'ldh1' + hidden, 'ldh2' + hidden,
                'z', 'ldz1' + hidden, 'ldz2' + hidden,
                'wr', 'wi',
                'dwork' + hidden, 'ldwork', 'info' + hidden]

    if not ldwork:
        p = H.shape[2]
        ldwork = max(1, ihi - ilo + p - 1)

    T, Z, Wr, Wi, info = _wrapper.mb03wd(
            job, compz, n, ilo, ihi, iloz, ihiz, H, Q, ldwork)

    if info < 0:
        e = ValueError(
                "Argument '{}' had an illegal value".format(arg_list[-info-1]))
        e.info = info
        raise e
    elif info > 0:
        warnings.warn(("failed to compute all the eigenvalues {ilo} to {ihi} "
                       "in a total of 30*({ihi}-{ilo}+1) iterations "
                       "the elements {i}:{ihi} of Wr contain those "
                       "eigenvalues which have been successfully computed."
                       ).format(i=info, ilo=ilo, ihi=ihi))
    if job == 'E':
        T = None
    if compz == 'N':
        Z = None

    W = Wr + Wi*1J
    return (T, Z, W)


def mb05md(a, delta, balanc='N'):
    """Ar, Vr, Yr, VAL = mb05md(a, delta, balanc='N')

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

        VAL : output rank-1 array('c') with bounds (n)
            Contains the eigenvalues of the matrix A. The eigenvalues
            are unordered except that complex conjugate pairs of values
            appear consecutively with the eigenvalue having positive
            imaginary part first.
    """
    hidden = ' (hidden by the wrapper)'
    arg_list = ['balanc', 'n', 'delta', 'a', 'lda'+hidden, 'v', 'ldv'+hidden,
                'y', 'ldy'+hidden, 'valr', 'vali',
                'iwork'+hidden, 'dwork'+hidden, 'ldwork'+hidden,
                'info'+hidden]
    n = min(a.shape)
    (Ar, Vr, Yr, VALr, VALi, INFO) = _wrapper.mb05md(balanc=balanc,
                                                     n=n,
                                                     delta=delta,
                                                     a=a)
    if INFO == 0:
        if not all(VALi == 0):
            VAL = VALr + 1J*VALi
        else:
            VAL = VALr
        return (Ar, Vr, Yr, VAL)
    elif INFO < 0:
        error_text = "The following argument had an illegal value: " \
                     + arg_list[-INFO-1]
    elif INFO > 0 and INFO <= n:
        error_text = "Incomplete eigenvalue calculation, missing %i eigenvalues" % INFO
    elif INFO == n+1:
        error_text = "Eigenvector matrix singular"
    elif INFO == n+2:
        error_text = "A matrix defective"
    e = ValueError(error_text)
    e.info = INFO
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
    arg_list = ['n', 'delta', 'a', 'lda'+hidden, 'ex', 'ldex'+hidden,
                'exint', 'ldexin'+hidden, 'tol', 'iwork'+hidden,
                'dwork'+hidden, 'ldwork'+hidden]
    n = min(a.shape)
    out = _wrapper.mb05nd(n=n, delta=delta, a=a, tol=tol)
    if out[-1] == 0:
        return out[:-1]
    elif out[-1] < 0:
        error_text = "The following argument had an illegal value: " \
                     + arg_list[-out[-1]-1]
    elif out[-1] == n+1:
        error_text = "Delta too large"
    e = ValueError(error_text)
    e.info = out[-1]
    raise e


# to be replaced by python wrappers
