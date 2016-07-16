#!/usr/bin/env python
#
#       identify.py
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


from slycot import _wrapper
import numpy as _np
import warnings


def ib01ad(meth, alg, jobd, conct, ctrl, nobr, data, rcond=-1, tol=-1):
    """To preprocess the input-output data for estimating the matrices of a
    linear time-invariant dynamical system and to find an
    estimate of the system order. The input-output data can,
    optionally, be processed sequentially.

    Required arguments
    ------------------

        meth : {'M', 'N'}
             Specifies the subspace identification method to be used,
             as follows:
             = 'M':  MOESP  algorithm with past inputs and outputs;
             = 'N':  N4SID  algorithm.
        alg  : {'C', 'F', 'Q'}
             Specifies the algorithm for computing the triangular
             factor R, as follows:
             = 'C':  Cholesky algorithm applied to the correlation
                     matrix of the input-output data;
             = 'F':  Fast QR algorithm;
             = 'Q':  QR algorithm applied to the concatenated block
                     Hankel matrices.
        jobd : {'M', 'N'}
             Specifies whether or not the matrices B and D should later
             be computed using the MOESP approach, as follows:
             = 'M':  the matrices B and D should later be computed
                     using the MOESP approach;
             = 'N':  the matrices B and D should not be computed using
                     the MOESP approach.
             This parameter is not relevant for METH = 'N'.
        conct : {'C', 'N'}
             Specifies whether or not the successive data blocks in
             sequential data processing belong to a single experiment,
             as follows:
             = 'C':  the current data block is a continuation of the
                     previous data block and/or it will be continued
                     by the next data block;
             = 'N':  there is no connection between the current data
                     block and the previous and/or the next ones.
             This parameter is not used if BATCH = 'O'.

        ctrl : {'C', 'N'}
             Specifies whether or not the user's confirmation of the
             system order estimate is desired, as follows:
             = 'C':  user's confirmation;
             = 'N':  no confirmation.
             If  CTRL = 'C',  a reverse communication routine,  IB01OY,
             is indirectly called (by SLICOT Library routine IB01OD),
             and, after inspecting the singular values and system order
             estimate,  n,  the user may accept  n  or set a new value.
             IB01OY  is not called if CTRL = 'N'.
        nobr : int
              The number of block rows,  s,  in the input and output
              block Hankel matrices to be processed.  NOBR > 0.
              (In the MOESP theory,  NOBR  should be larger than  n,
              the estimated dimension of state vector.)
        data : [(u1, y1), (u2, y2)...]. A list of tuples containing experiment
               data of inputs and outputs. The number of inputs and number of
               outputs should be the same in every tuple. u1 has dimension
               (n1 x nu). y1 has dimension (n1 x ny). u2 has dimension (n2 x nu)

    Optional:

      rcond   DOUBLE
              The tolerance to be used for estimating the rank of
              matrices. If the user sets  RCOND > 0,  the given value
              of  RCOND  is used as a lower bound for the reciprocal
              condition number;  an m-by-n matrix whose estimated
              condition number is less than  1/RCOND  is considered to
              be of full rank.  If the user sets  RCOND <= 0,  then an
              implicitly computed, default tolerance, defined by
              RCONDEF = m*n*EPS,  is used instead, where  EPS  is the
              relative machine precision (see LAPACK Library routine
              DLAMCH).
              This parameter is not used for  METH = 'M'.

      tol     DOUBLE
              Absolute tolerance used for determining an estimate of
              the system order. If  TOL >= 0,  the estimate is
              indicated by the index of the last singular value greater
              than or equal to  TOL.  (Singular values less than  TOL
              are considered as zero.) When  TOL = 0,  an internally
              computed default value,  TOL = NOBR*EPS*SV(1),  is used,
              where  SV(1)  is the maximal singular value, and  EPS  is
              the relative machine precision (see LAPACK Library routine
              DLAMCH). When  TOL < 0,  the estimate is indicated by the
              index of the singular value that has the largest
              logarithmic gap to its successor.
    Outputs:

        R       (output) rank-2 array('d'), dimension
                ( LDR,2*(M+L)*NOBR )
                On exit, if ALG = 'C' and BATCH = 'F' or 'I', the leading
                2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this
                array contains the current upper triangular part of the
                correlation matrix in sequential data processing.
                If ALG = 'F' and BATCH = 'F' or 'I', the array R is not
                referenced.
                On exit, if INFO = 0, ALG = 'Q', and BATCH = 'F' or 'I',
                the leading 2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular
                part of this array contains the current upper triangular
                factor R from the QR factorization of the concatenated
                block Hankel matrices. Denote  R_ij, i,j = 1:4,  the
                ij submatrix of  R,  partitioned by M*NOBR,  M*NOBR,
                L*NOBR,  and  L*NOBR  rows and columns.
                On exit, if INFO = 0 and BATCH = 'L' or 'O', the leading
                2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of
                this array contains the matrix S, the processed upper
                triangular factor R from the QR factorization of the
                concatenated block Hankel matrices, as required by other
                subroutines. Specifically, let  S_ij, i,j = 1:4,  be the
                ij submatrix of  S,  partitioned by M*NOBR,  L*NOBR,
                M*NOBR,  and  L*NOBR  rows and columns. The submatrix
                S_22  contains the matrix of left singular vectors needed
                subsequently. Useful information is stored in  S_11  and
                in the block-column  S_14 : S_44.  For METH = 'M' and
                JOBD = 'M', the upper triangular part of  S_31  contains
                the upper triangular factor in the QR factorization of the
                matrix  R_1c = [ R_12'  R_22'  R_11' ]',  and  S_12
                contains the corresponding leading part of the transformed
                matrix  R_2c = [ R_13'  R_23'  R_14' ]'.  For  METH = 'N',
                the subarray  S_41 : S_43  contains the transpose of the
                matrix contained in  S_14 : S_34.
                The details of the contents of R need not be known if this
                routine is followed by SLICOT Library routine IB01BD.
                On entry, if ALG = 'C', or ALG = 'Q', and BATCH = 'I' or
                'L', the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  upper
                triangular part of this array must contain the upper
                triangular matrix R computed at the previous call of this
                routine in sequential data processing. The array R need
                not be set on entry if ALG = 'F' or if BATCH = 'F' or 'O'.
         N       (output) int
                 The estimated order of the system.
                 If  CTRL = 'C',  the estimated order has been reset to a
                 value specified by the user.
       sv      (output) rank-1 rray('d'), dimension ( L*NOBR )
               The singular values used to estimate the system order.
       iwarn :  int
               = 0:  no warning;
               = 1:  the number of 100 cycles in sequential data
                     processing has been exhausted without signaling
                     that the last block of data was get; the cycle
                     counter was reinitialized;
               = 2:  a fast algorithm was requested (ALG = 'C' or 'F'),
                     but it failed, and the QR algorithm was then used
                     (non-sequential data processing);
               = 3:  all singular values were exactly zero, hence  N = 0
                     (both input and output were identically zero);
               = 4:  the least squares problems with coefficient matrix
                     U_f,  used for computing the weighted oblique
                     projection (for METH = 'N'), have a rank-deficient
                     coefficient matrix;
               = 5:  the least squares problem with coefficient matrix
                     r_1  [6], used for computing the weighted oblique
                     projection (for METH = 'N'), has a rank-deficient
                     coefficient matrix.
               NOTE: the values 4 and 5 of IWARN have no significance
                     for the identification problem.

    Raises
    ------
        ValueError : e
            e.info contains information about the exact type of exception
             < 0:  if INFO = -i, the i-th argument had an illegal
                   value;
             = 1:  a fast algorithm was requested (ALG = 'C', or 'F')
                   in sequential data processing, but it failed; the
                   routine can be repeatedly called again using the
                   standard QR algorithm;
             = 2:  the singular value decomposition (SVD) algorithm did
                   not converge.
 """

    n_exp = len(data)     # number of experiments
    m = data[0][0].shape[1]   # m: number of inputs
    l = data[0][1].shape[1]   # l: number of outputs
    nsmpl = 0                 # total number of samples

    if meth == 'M' and jobd == 'M':
        ldr = max(2*(m+l)*nobr, 3*m*nobr)
    elif meth == 'N' or (meth == 'M' and jobd == 'N'):
        ldr = 2*(m+l)*nobr
    else:
        raise ValueError("ib01ad: could not handle 'ldr' case")

    r = _np.zeros((ldr, 2*(m+l)*nobr), order='F')
    for i in range(0, n_exp):

        if n_exp == 1:
            batch = 'O'        # one block only
        elif i == 0:
            batch = 'F'        # first block
        elif i == n_exp-1:
            batch = 'L'        # last block
        else:
            batch = 'I'       # intermediate block

        y = data[i][1]
        u = data[i][0]

        # y.rows == u.rows  is checked by iddata class
        # octave_idx_type m = u.columns ();   // m: number of inputs
        # octave_idx_type l = y.columns ();   // l: number of outputs
        nsmp = y.shape[0]   # nsmp: number of samples in the current experiment
        nsmpl += nsmp      # nsmpl: total number of samples of all experiments

        # minimal nsmp size checked by __slicot_identification__.m
        if batch == 'O':
            if nsmp < 2*(m+l+1)*nobr - 1:
                raise ValueError("identify: require NSMP >= 2*(M+L+1)*NOBR - 1")
        else:
            if nsmp < 2*nobr:
                raise ValueError("identify: require NSMP >= 2*NOBR")
        if batch == 'L':
            if nsmpl < 2*(m+l+1)*nobr - 1:
                raise ValueError('''identify: total number of samples of all
                                 experiments should be >= 2*(M+L+1)*NOBR -1''')

        # if m == 0:
        #    ldu = 1
        # else:                    # m > 0
        #    ldu = nsmp

        # ldy = nsmp

        # workspace

        if meth == 'N':            # if METH = 'N'
            liwork = (m+l)*nobr
        elif alg == 'F':       # if METH = 'M' and ALG = 'F'
            liwork = m+l
        else:                        # if METH = 'M' and ALG = 'C' or 'Q'
            liwork = 0

        # TODO: Handle 'k' for DWORK


#             The length of the array DWORK.
        ns = nsmp - 2*nobr + 1

        if (alg == 'C'):
            if (batch == 'F' or batch == 'I'):
                if (conct == 'C'):
                    ldwork = (4*nobr-2)*(m+l)
                else:    # (conct == 'N')
                    ldwork = 1
            elif (meth == 'M'):   # and (batch == 'L' or batch == 'O')
                if(conct == 'C' and batch == 'L'):
                    ldwork = max((4*nobr-2)*(m+l), 5*l*nobr)
                elif(jobd == 'M'):
                    ldwork = max((2*m-1)*nobr, (m+l)*nobr, 5*l*nobr)
                else:    # (jobd == 'N')
                    ldwork = 5*l*nobr
            else:    # meth_b == 'N' and (batch == 'L' or batch == 'O')
                ldwork = 5*(m+l)*nobr + 1
        elif alg == 'F':
            if batch != 'O' and conct == 'C':
                ldwork = (m+l)*2*nobr*(m+l+3)
            elif batch == 'F' or batch == 'I':  # and conct == 'N'
                ldwork = (m+l)*2*nobr*(m+l+1)
            else:    # (batch == 'L' or '0' and conct == 'N')
                ldwork = (m+l)*4*nobr*(m+l+1)+(m+l)*2*nobr
        else:    # (alg == 'Q')
            # octave_idx_type ns = nsmp - 2*nobr + 1

            if ldr >= ns and batch == 'F':
                ldwork = 4*(m+l)*nobr
            elif ldr >= ns and batch == 'O':
                if meth == 'M':
                    ldwork = max(4*(m+l)*nobr, 5*l*nobr)
                else:    # (meth == 'N')
                    ldwork = 5*(m+l)*nobr + 1
            elif(conct == 'C' and (batch == 'I' or batch == 'L')):
                ldwork = 4*(nobr+1)*(m+l)*nobr
            else:
                # if ALG = 'Q', (BATCH = 'F' or 'O', and LDR < NS),
                # or (BATCH = 'I' or 'L' and CONCT = 'N')
                ldwork = 6*(m+l)*nobr

        '''
        IB01AD.f Lines 438-445
        C     FURTHER COMMENTS
        C
        C     For ALG = 'Q', BATCH = 'O' and LDR < NS, or BATCH <> 'O', the
        C     calculations could be rather inefficient if only minimal workspace
        C     (see argument LDWORK) is provided. It is advisable to provide as
        C     much workspace as possible. Almost optimal efficiency can be
        C     obtained for  LDWORK = (NS+2)*(2*(M+L)*NOBR),  assuming that the
        C     cache size is large enough to accommodate R, U, Y, and DWORK.
        */'''

        ldwork = max(ldwork, (ns+2)*(2*(m+l)*nobr))

        '''
        IB01AD.f Lines 291-195:
        c             the workspace used for alg = 'q' is
        c                       ldrwrk*2*(m+l)*nobr + 4*(m+l)*nobr,
        c             where ldrwrk = ldwork/(2*(m+l)*nobr) - 2; recommended
        c             value ldrwrk = ns, assuming a large enough cache size.
        c             for good performance,  ldwork  should be larger.

        somehow ldrwrk and ldwork must have been mixed up here

            '''

        iwork = _np.zeros((liwork,), dtype=_np.int)
        dwork = _np.zeros((ldwork,))

        args_list = ['meth', 'alg', 'jobd', 'batch', 'conct', 'ctrl', 'nobr',
                     'm', 'l', 'nsmp', 'u', 'ldu', 'y', 'ldy', 'n', 'r', 'ldr',
                     'sv', 'rcond', 'tol', 'iwork', 'dwork', 'ldwork', 'iwarn',
                     'info']

        out = _wrapper.ib01ad(meth, alg, jobd, batch, conct, ctrl, nobr,
                              nsmp, u, y, r, rcond, tol, iwork, dwork)
        out = (r,) + out

        if out[-1] < 0:
            error_text = '''The following argument had an illegal
            value: ''' + args_list[-out[-1]-1]
            e = ValueError(error_text)
            e.info = out[-1]
            raise e
        if out[-1] == 1:
            e = ValueError('''a fast algorithm was requested (ALG = 'C', or 'F')
                   in sequential data processing, but it failed; the
                   routine can be repeatedly called again using the
                   standard QR algorithm''')
            e.info = out[-1]
            raise e
        elif out[-1] == 2:
            e = ValueError('''the singular value decomposition (SVD) algorithm
                           did not converge.''')
            e.info = out[-1]
            raise e

    return out[:-1]


def ib01bd(meth, job, jobck, nobr, n, m, l, nsmpl, r, A, C, tol=0):
    '''
     PURPOSE

     To estimate the system matrices A, C, B, and D, the noise
     covariance matrices Q, Ry, and S, and the Kalman gain matrix K
     of a linear time-invariant state space model, using the
     processed triangular factor R of the concatenated block Hankel
     matrices, provided by SLICOT Library routine IB01AD.

     METHOD

     In the MOESP approach, the matrices  A  and  C  are first
     computed from an estimated extended observability matrix [1],
     and then, the matrices  B  and  D  are obtained by solving an
     extended linear system in a least squares sense.
     In the N4SID approach, besides the estimated extended
     observability matrix, the solutions of two least squares problems
     are used to build another least squares problem, whose solution
     is needed to compute the system matrices  A,  C,  B,  and  D.  The
     solutions of the two least squares problems are also optionally
     used by both approaches to find the covariance matrices.
     The Kalman gain matrix is obtained by solving a discrete-time
     algebraic Riccati equation.

     INPUTS
     -------

     Mode Parameters

     meth    CHARACTER*1
             Specifies the subspace identification method to be used,
             as follows:
             = 'M':  MOESP  algorithm with past inputs and outputs;
             = 'N':  N4SID  algorithm;
             = 'C':  combined method:  MOESP  algorithm for finding the
                     matrices A and C, and  N4SID  algorithm for
                     finding the matrices B and D.

     job     CHARACTER*1
             Specifies which matrices should be computed, as follows:
             = 'A':  compute all system matrices, A, B, C, and D;
             = 'C':  compute the matrices A and C only;
             = 'B':  compute the matrix B only;
             = 'D':  compute the matrices B and D only.

     jobck   CHARACTER*1
             Specifies whether or not the covariance matrices and the
             Kalman gain matrix are to be computed, as follows:
             = 'C':  the covariance matrices only should be computed;
             = 'K':  the covariance matrices and the Kalman gain
                     matrix should be computed;
             = 'N':  the covariance matrices and the Kalman gain matrix
                     should not be computed.

     nobr    INTEGER
             The number of block rows,  s,  in the input and output
             Hankel matrices processed by other routines.  NOBR > 1.

     n       INTEGER
             The order of the system.  NOBR > N > 0.

     m       INTEGER
             The number of system inputs.  M >= 0.

     l       INTEGER
             The number of system outputs.  L > 0.

     nsmpl   INTEGER
             If  JOBCK = 'C' or 'K',  the total number of samples used
             for calculating the covariance matrices.
             NSMPL >= 2*(M+L)*NOBR.
             This parameter is not meaningful if  JOBCK = 'N'.

     r       rank-2 array type-d, dimension
             ( LDR,2*(M+L)*NOBR )
             On entry, the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  part
             of this array must contain the relevant data for the MOESP
             or N4SID algorithms, as constructed by SLICOT Library
             routine IB01AD. Let  R_ij,  i,j = 1:4,  be the
             ij submatrix of  R  (denoted  S  in IB01AD),  partitioned
             by  M*NOBR,  L*NOBR,  M*NOBR,  and  L*NOBR  rows and
             columns. The submatrix  R_22  contains the matrix of left
             singular vectors used. Also needed, for  METH = 'N'  or
             JOBCK <> 'N',  are the submatrices  R_11,  R_14 : R_44,
             and, for  METH = 'M' or 'C'  and  JOB <> 'C', the
             submatrices  R_31  and  R_12,  containing the processed
             matrices  R_1c  and  R_2c,  respectively, as returned by
             SLICOT Library routine IB01AD.
             Moreover, if  METH = 'N'  and  JOB = 'A' or 'C',  the
             block-row  R_41 : R_43  must contain the transpose of the
             block-column  R_14 : R_34  as returned by SLICOT Library
             routine IB01AD.
             The remaining part of  R  is used as workspace.
             On exit, part of this array is overwritten. Specifically,
             if  METH = 'M',  R_22  and  R_31  are overwritten if
                 JOB = 'B' or 'D',  and  R_12,  R_22,  R_14 : R_34,
                 and possibly  R_11  are overwritten if  JOBCK <> 'N';
             if  METH = 'N',  all needed submatrices are overwritten.
             The details of the contents of  R  need not be known if
             this routine is called once just after calling the SLICOT
             Library routine IB01AD.

     A       (input or output) DOUBLE PRECISION array, dimension
             (LDA,N)
             On entry, if  METH = 'N' or 'C'  and  JOB = 'B' or 'D',
             the leading N-by-N part of this array must contain the
             system state matrix.
             If  METH = 'M'  or  (METH = 'N' or 'C'  and JOB = 'A'
             or 'C'),  this array need not be set on input.
             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  the
             leading N-by-N part of this array contains the system
             state matrix.

     C       (input or output) DOUBLE PRECISION array, dimension
             (LDC,N)
             On entry, if  METH = 'N' or 'C'  and  JOB = 'B' or 'D',
             the leading L-by-N part of this array must contain the
             system output matrix.
             If  METH = 'M'  or  (METH = 'N' or 'C'  and JOB = 'A'
             or 'C'),  this array need not be set on input.
             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  or
             INFO = 3  (or  INFO >= 0,  for  METH = 'M'),  the leading
             L-by-N part of this array contains the system output
             matrix.

     Tolerances

     TOL     DOUBLE
             The tolerance to be used for estimating the rank of
             matrices. If the user sets  TOL > 0,  then the given value
             of  TOL  is used as a lower bound for the reciprocal
             condition number;  an m-by-n matrix whose estimated
             condition number is less than  1/TOL  is considered to
             be of full rank.  If the user sets  TOL <= 0,  then an
             implicitly computed, default tolerance, defined by
             TOLDEF = m*n*EPS,  is used instead, where  EPS  is the
             relative machine precision (see LAPACK Library routine
             DLAMCH).

     OUTPUTS
     --------
     B       (output) rank-2 type('d') array, dimension (LDB,M)
             If  M > 0,  JOB = 'A', 'B', or 'D'  and  INFO = 0,  the
             leading N-by-M part of this array contains the system
             input matrix. If  M = 0  or  JOB = 'C',  this array is
             not referenced.

     D       (output) rank-2 type('d') array, dimension (LDD,M)
             If  M > 0,  JOB = 'A' or 'D'  and  INFO = 0,  the leading
             L-by-M part of this array contains the system input-output
             matrix. If  M = 0  or  JOB = 'C' or 'B',  this array is
             not referenced.

     Q       (output) rank-2 type('d') array, dimension (LDQ,N)
             If  JOBCK = 'C' or 'K',  the leading N-by-N part of this
             array contains the positive semidefinite state covariance
             matrix. If  JOBCK = 'K',  this matrix has been used as
             state weighting matrix for computing the Kalman gain.
             This parameter is not referenced if JOBCK = 'N'.

     RY      (output) rank-2 type('d') array, dimension (LDRY,L)
             If  JOBCK = 'C' or 'K',  the leading L-by-L part of this
             array contains the positive (semi)definite output
             covariance matrix. If  JOBCK = 'K',  this matrix has been
             used as output weighting matrix for computing the Kalman
             gain.
             This parameter is not referenced if JOBCK = 'N'.

     S       (output) rank-2 type('d') array, dimension (LDS,L)
             If  JOBCK = 'C' or 'K',  the leading N-by-L part of this
             array contains the state-output cross-covariance matrix.
             If  JOBCK = 'K',  this matrix has been used as state-
             output weighting matrix for computing the Kalman gain.
             This parameter is not referenced if JOBCK = 'N'.

     K       (output) rank-2 type('d') array, dimension ( LDK,L )
             If  JOBCK = 'K',  the leading  N-by-L  part of this array
             contains the estimated Kalman gain matrix.
             If  JOBCK = 'C' or 'N',  this array is not referenced.

     Warning Indicator

     iwarn   INTEGER
             = 0:  no warning;
             = 4:  a least squares problem to be solved has a
                   rank-deficient coefficient matrix;
             = 5:  the computed covariance matrices are too small.
                   The problem seems to be a deterministic one; the
                   gain matrix is set to zero.

     Error Indicator

     info    INTEGER
             = 0:  successful exit;
             < 0:  if INFO = -i, the i-th argument had an illegal
                   value;
             = 2:  the singular value decomposition (SVD) algorithm did
                   not converge;
             = 3:  a singular upper triangular matrix was found;
             = 3+i:  if  JOBCK = 'K'  and the associated Riccati
                   equation could not be solved, where i = 1,...,6;
                   (see the description of the parameter INFO for the
                   SLICOT Library routine SB02RD for the meaning of
                   the i values);
             = 10: the QR algorithm did not converge.

    '''
    if (nsmpl < 2*(m+l)*nobr):
        raise ValueError("identify : nsmpl (%d) < 2*(m+l)*nobr (%d)",
                         nsmpl, nobr)

    # arguments out

    lda = max(1, n)
    ldc = max(1, l)
    ldb = max(1, n)
    ldd = max(1, l)
    if jobck == 'C' or jobck == 'K':
        ldq = n            # if JOBCK = 'C' or 'K'
        ldry = l           # if JOBCK = 'C' or 'K'
        lds = n            # if JOBCK = 'C' or 'K'
    else:
        ldq = 1
        ldry = 1
        lds = 1
    if jobck == 'K':
        ldk = n            # if JOBCK = 'K'
    else:
        ldk = 1

    if (meth == 'N' or meth == 'C') and (job == 'B' or job == 'D'):
        if A.shape[0] < n or A.shape[1] < n:
            raise ValueError('''identify_ib01bd : Matrix A input shape not
                             correct''')
        if C.shape[0] < l or C.shape[1] < n:
            raise ValueError('''identify_ib01bd : Matrix C input shape not
                             correct''')
    else:
        A = _np.zeros((lda, n))
        C = _np.zeros((ldc, n))

    # workspace

    liw1 = max(n, m*nobr+n, l*nobr, m*(n+l))
    liw2 = n*n
    liwork = max(liw1, liw2)

    if (meth == 'M'):

        ldw1a = max(2*(l*nobr-l)*n+2*n, (l*nobr-l)*n+n*n+7*n)
        ldw1b = max(2*(l*nobr-l)*n+n*n+7*n, (l*nobr-l)*n+n+6*m*nobr,
                    (l*nobr-l)*n+n+max(l+m*nobr, l*nobr + max(3*l*nobr+1, m)))
        ldw1 = max(ldw1a, ldw1b)
        #THE FOLLOWING LINE IS NOT THERE IN THE DOC OF IB01BD BUT CODE DOES
        #NOT WORK WITOUT THIS
        ldw1 = ldw1 + l*nobr*n
        if (m == 0 or job == 'C'):
            aw = n + n*n
        else:
            aw = 0

        ldw2 = l*nobr*n + max((l*nobr-l)*n+aw+2*n+max(5*n, (2*m+l)*nobr+l),
                              4*(m*nobr+n)+1, m*nobr+2*n+l)

    elif(meth == 'N'):
        ldw1 = l*nobr*n + max((l*nobr-l)*n+2*n+(2*m+l)*nobr+l,
                              2*(l*nobr-l)*n+n*n+8*n, n+4*(m*nobr+n)+1,
                              m*nobr+3*n+l)

        if (m == 0 or job == 'C'):
            ldw2 = 0
        else:
            ldw2 = l*nobr*n+m*nobr*(n+l)*(m*(n+l)+1) + max((n+l)*(n+l),
                                                           4*m*(n+l)+1)

    else:    # (meth == 'C')
        ldw1a = max(2*(l*nobr-l)*n+2*n, (l*nobr-l)*n+n*n+7*n)
        ldw1b = l*nobr*n + max((l*nobr-l)*n+2*n+(2*m+l)*nobr+l,
                               2*(l*nobr-l)*n+n*n+8*n,
                               n+4*(m*nobr+n)+1,
                               m*nobr+3*n+l)

        ldw1 = max(ldw1a, ldw1b)

        ldw2 = l*nobr*n+m*nobr*(n+l)*(m*(n+l)+1) + max((n+l)*(n+l), 4*m*(n+l)+1)

    ldw3 = max(4*n*n + 2*n*l + l*l + max(3*l, n*l), 14*n*n + 12*n + 5)
    ldwork = max(ldw1, ldw2, ldw3)

    iwork = _np.zeros((liwork, ), dtype=_np.int)
    dwork = _np.zeros((ldwork, ))
    bwork = _np.zeros((2*n, ), dtype=_np.int32)

    # error indicators

    # SLICOT routine IB01BD
    out = _wrapper.ib01bd(meth, job, jobck, nobr, m, l, nsmpl, r, A, lda, C,
                          ldc, ldb, ldd, ldq, ldry, lds, ldk, iwork, dwork,
                          bwork)

    err_msg = (
        "0: OK",
        "1: error message not specified",
        "2: the singular value decomposition (SVD) algorithm did not converge",
        "3: a singular upper triangular matrix was found",
        "4: matrix A is (numerically) singular in discrete-time case",
        "5: the Hamiltonian or symplectic matrix H cannot be reduced to real"
        "Schur form",
        "6: the real Schur form of the Hamiltonian or "
        "symplectic matrix H cannot be appropriately ordered",
        "7: the Hamiltonian or symplectic matrix H has less "
        "than N stable eigenvalues",
        "8: the N-th order system of linear algebraic "
        "equations, from which the solution matrix X would "
        "be obtained, is singular to working precision",
        "9: the QR algorithm failed to complete the reduction "
        "of the matrix Ac to Schur canonical form, T",
        "10: the QR algorithm did not converge")

    '''warn_msg = (
        "0: OK",
        "1: warning message not specified",
        "2: warning message not specified",
        "3: warning message not specified",
        "4: a least squares problem to be solved has a "
            "rank-deficient coefficient matrix",
        "5: the computed covariance matrices are too small. "
            "The problem seems to be a deterministic one; the "
            "gain matrix is set to zero")
    '''
    args_list = ['meth', 'job', 'jobck', 'nobr', 'n', 'm', 'l', 'nsmpl', 'r',
                 'ldr', 'a', 'lda', 'c', 'ldc', 'b', 'ldb', 'd', 'ldd', 'q',
                 'ldq', 'ry', 'ldry', 's', 'lds', 'k', 'ldk', 'tol', 'iwork',
                 'dwork', 'ldwork', 'bwork', 'iwarn', 'info']

    if out[-1] < 0:
        error_text = '''The following argument had an illegal
        value: ''' + args_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    elif out[-1] > 0:
        e = ValueError(err_msg[out[-1]])
        raise e

    # error_msg ("ident: IB01BD", info_b, 10, err_msg_b);
    # warning_msg ("ident: IB01BD", iwarn_b, 5, warn_msg_b);

    return out
