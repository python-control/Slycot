      SUBROUTINE MB04XD( JOBU, JOBV, M, N, RANK, THETA, A, LDA, U, LDU,
     $                   V, LDV, Q, INUL, TOL, RELTOL, DWORK, LDWORK,
     $                   IWARN, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To compute a basis for the left and/or right singular subspace of
C     an M-by-N matrix A corresponding to its smallest singular values.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBU    CHARACTER*1
C             Specifies whether to compute the left singular subspace
C             as follows:
C             = 'N':  Do not compute the left singular subspace;
C             = 'A':  Return the (M - RANK) base vectors of the desired
C                     left singular subspace in U;
C             = 'S':  Return the first (min(M,N) - RANK) base vectors
C                     of the desired left singular subspace in U.
C
C     JOBV    CHARACTER*1
C             Specifies whether to compute the right singular subspace
C             as follows:
C             = 'N':  Do not compute the right singular subspace;
C             = 'A':  Return the (N - RANK) base vectors of the desired
C                     right singular subspace in V;
C             = 'S':  Return the first (min(M,N) - RANK) base vectors
C                     of the desired right singular subspace in V.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns in matrix A.  N >= 0.
C
C     RANK    (input/output) INTEGER
C             On entry, if RANK < 0, then the rank of matrix A is
C             computed by the routine as the number of singular values
C             greater than THETA.
C             Otherwise, RANK must specify the rank of matrix A.
C             RANK <= min(M,N).
C             On exit, if RANK < 0 on entry, then RANK contains the
C             computed rank of matrix A. That is, the number of singular
C             values of A greater than THETA.
C             Otherwise, the user-supplied value of RANK may be changed
C             by the routine on exit if the RANK-th and the (RANK+1)-th
C             singular values of A are considered to be equal.
C             See also the description of parameter TOL below.
C
C     THETA   (input/output) DOUBLE PRECISION
C             On entry, if RANK < 0, then THETA must specify an upper
C             bound on the smallest singular values of A corresponding
C             to the singular subspace to be computed.  THETA >= 0.0.
C             Otherwise, THETA must specify an initial estimate (t say)
C             for computing an upper bound on the (min(M,N) - RANK)
C             smallest singular values of A. If THETA < 0.0, then t is
C             computed by the routine.
C             On exit, if RANK >= 0 on entry, then THETA contains the
C             computed upper bound such that precisely RANK singular
C             values of A are greater than THETA + TOL.
C             Otherwise, THETA is unchanged.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading M-by-N part of this array must contain the
C             matrix A from which the basis of a desired singular
C             subspace is to be computed.
C             NOTE that this array is destroyed.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,M).
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,*)
C             If JOBU = 'A', then the leading M-by-M part of this array
C             contains the (M - RANK) M-dimensional base vectors of the
C             desired left singular subspace of A corresponding to its
C             singular values less than or equal to THETA. These vectors
C             are stored in the i-th column(s) of U for which
C             INUL(i) = .TRUE., where i = 1,2,...,M.
C
C             If JOBU = 'S', then the leading M-by-min(M,N) part of this
C             array contains the first (min(M,N) - RANK) M-dimensional
C             base vectors of the desired left singular subspace of A
C             corresponding to its singular values less than or equal to
C             THETA. These vectors are stored in the i-th column(s) of U
C             for which INUL(i) = .TRUE., where i = 1,2,..., min(M,N).
C
C             Otherwise, U is not referenced (since JOBU = 'N') and can
C             be supplied as a dummy array (i.e. set parameter LDU = 1
C             and declare this array to be U(1,1) in the calling
C             program).
C
C     LDU     INTEGER
C             The leading dimension of array U.
C             LDU >= max(1,M) if JOBU = 'A' or JOBU = 'S',
C             LDU >= 1        if JOBU = 'N'.
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,*)
C             If JOBV = 'A', then the leading N-by-N part of this array
C             contains the (N - RANK) N-dimensional base vectors of the
C             desired right singular subspace of A corresponding to its
C             singular values less than or equal to THETA. These vectors
C             are stored in the i-th column(s) of V for which
C             INUL(i) = .TRUE., where i = 1,2,...,N.
C
C             If JOBV = 'S', then the leading N-by-min(M,N) part of this
C             array contains the first (min(M,N) - RANK) N-dimensional
C             base vectors of the desired right singular subspace of A
C             corresponding to its singular values less than or equal to
C             THETA. These vectors are stored in the i-th column(s) of V
C             for which INUL(i) = .TRUE., where i = 1,2,...,MIN( M,N).
C
C             Otherwise, V is not referenced (since JOBV = 'N') and can
C             be supplied as a dummy array (i.e. set parameter LDV = 1
C             and declare this array to be V(1,1) in the calling
C             program).
C
C     LDV     INTEGER
C             The leading dimension of array V.
C             LDV >= max(1,N) if JOBV = 'A' or JOBV = 'S',
C             LDV >= 1        if JOBV = 'N'.
C
C     Q       (output) DOUBLE PRECISION array, dimension (2*min(M,N)-1)
C             This array contains the partially diagonalized bidiagonal
C             matrix J computed from A, at the moment that the desired
C             singular subspace has been found. Specifically, the
C             leading p = min(M,N) entries of Q contain the diagonal
C             elements q(1),q(2),...,q(p) and the entries Q(p+1),
C             Q(p+2),...,Q(2*p-1) contain the superdiagonal elements
C             e(1),e(2),...,e(p-1) of J.
C
C     INUL    (output) LOGICAL array, dimension (max(M,N))
C             If JOBU <> 'N' or JOBV <> 'N', then the indices of the
C             elements of this array with value .TRUE. indicate the
C             columns in U and/or V containing the base vectors of the
C             desired left and/or right singular subspace of A. They
C             also equal the indices of the diagonal elements of the
C             bidiagonal submatrices in the array Q, which correspond
C             to the computed singular subspaces.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             This parameter defines the multiplicity of singular values
C             by considering all singular values within an interval of
C             length TOL as coinciding. TOL is used in checking how many
C             singular values are less than or equal to THETA. Also in
C             computing an appropriate upper bound THETA by a bisection
C             method, TOL is used as a stopping criterion defining the
C             minimum (absolute) subinterval width. TOL is also taken
C             as an absolute tolerance for negligible elements in the
C             QR/QL iterations. If the user sets TOL to be less than or
C             equal to 0, then the tolerance is taken as specified in
C             SLICOT Library routine MB04YD document.
C
C     RELTOL  DOUBLE PRECISION
C             This parameter specifies the minimum relative width of an
C             interval. When an interval is narrower than TOL, or than
C             RELTOL times the larger (in magnitude) endpoint, then it
C             is considered to be sufficiently small and bisection has
C             converged. If the user sets RELTOL to be less than
C             BASE * EPS, where BASE is machine radix and EPS is machine
C             precision (see LAPACK Library routine DLAMCH), then the
C             tolerance is taken as BASE * EPS.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK = max(1, LDW + max(2*P + max(M,N), LDY)), where
C                  P = min(M,N);
C                LDW = max(2*N, N*(N+1)/2), if JOBU <> 'N' and M large
C                                                        enough than N;
C                LDW = 0,                   otherwise;
C                LDY = 8*P - 5, if JOBU <> 'N' or  JOBV <> 'N';
C                LDY = 6*P - 3, if JOBU =  'N' and JOBV =  'N'.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  if the rank of matrix A (as specified by the user)
C                   has been lowered because a singular value of
C                   multiplicity greater than 1 was found.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the maximum number of QR/QL iteration steps
C                   (30*MIN(M,N)) has been exceeded.
C
C     METHOD
C
C     The method used is the Partial Singular Value Decomposition (PSVD)
C     approach proposed by Van Huffel, Vandewalle and Haegemans, which
C     is an efficient technique (see [1]) for computing the singular
C     subspace of a matrix corresponding to its smallest singular
C     values. It differs from the classical SVD algorithm [3] at three
C     points, which results in high efficiency. Firstly, the Householder
C     transformations of the bidiagonalization need only to be applied
C     on the base vectors of the desired singular subspaces; secondly,
C     the bidiagonal matrix need only be partially diagonalized; and
C     thirdly, the convergence rate of the iterative diagonalization can
C     be improved by an appropriate choice between QL and QR iterations.
C     (Note, however, that LAPACK Library routine DGESVD, for computing
C     SVD, also uses either QL and QR iterations.) Depending on the gap,
C     the desired numerical accuracy and the dimension of the desired
C     singular subspace, the PSVD can be up to three times faster than
C     the classical SVD algorithm.
C
C     The PSVD algorithm [1-2] for an M-by-N matrix A proceeds as
C     follows:
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C      (a) If M is large enough than N, transform A into upper
C          triangular form R.
C
C      (b) Transform A (or R) into bidiagonal form:
C
C                |q(1) e(1)  0   ...  0   |
C           (0)  | 0   q(2) e(2)      .   |
C          J   = | .                  .   |
C                | .                e(N-1)|
C                | 0            ...  q(N) |
C
C     if M >= N, or
C
C                |q(1)  0    0   ...  0     0   |
C           (0)  |e(1) q(2)  0        .     .   |
C          J   = | .                  .     .   |
C                | .                 q(M-1) .   |
C                | 0             ... e(M-1) q(M)|
C
C     if M < N, using Householder transformations.
C     In the second case, transform the matrix to the upper bidiagonal
C     form by applying Givens rotations.
C
C      (c) If U is requested, initialize U with the identity matrix.
C          If V is requested, initialize V with the identity matrix.
C
C     Step 2: Partial diagonalization phase
C             -----------------------------
C     If the upper bound THETA is not given, then compute THETA such
C     that precisely (min(M,N) - RANK) singular values of the bidiagonal
C     matrix are less than or equal to THETA, using a bisection method
C     [4]. Diagonalize the given bidiagonal matrix J partially, using
C     either QR iterations (if the upper left diagonal element of the
C     considered bidiagonal submatrix is larger than the lower right
C     diagonal element) or QL iterations, such that J is split into
C     unreduced bidiagonal submatrices whose singular values are either
C     all larger than THETA or all less than or equal to THETA.
C     Accumulate the Givens rotations in U and/or V (if desired).
C
C     Step 3: Back transformation phase
C             -------------------------
C      (a) Apply the Householder transformations of Step 1(b) onto the
C          columns of U and/or V associated with the bidiagonal
C          submatrices with all singular values less than or equal to
C          THETA (if U and/or V is desired).
C
C      (b) If M is large enough than N, and U is desired, then apply the
C          Householder transformations of Step 1(a) onto each computed
C          column of U in Step 3(a).
C
C     REFERENCES
C
C     [1] Van Huffel, S., Vandewalle, J. and Haegemans, A.
C         An efficient and reliable algorithm for computing the singular
C         subspace of a matrix associated with its smallest singular
C         values.
C         J. Comput. and Appl. Math., 19, pp. 313-330, 1987.
C
C     [2] Van Huffel, S.
C         Analysis of the total least squares problem and its use in
C         parameter estimation.
C         Doctoral dissertation, Dept. of Electr. Eng., Katholieke
C         Universiteit Leuven, Belgium, June 1987.
C
C     [3] Chan, T.F.
C         An improved algorithm for computing the singular value
C         decomposition.
C         ACM TOMS, 8, pp. 72-83, 1982.
C
C     [4] Van Huffel, S. and Vandewalle, J.
C         The partial total least squares algorithm.
C         J. Comput. and Appl. Math., 21, pp. 333-341, 1988.
C
C     NUMERICAL ASPECTS
C
C     Using the PSVD a large reduction in computation time can be
C     gained in total least squares applications (cf [2 - 4]), in the
C     computation of the null space of a matrix and in solving
C     (non)homogeneous linear equations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB04PD by S. Van Huffel, Katholieke
C     University Leuven, Belgium.
C
C     REVISIONS
C
C     July 10, 1997.
C
C     KEYWORDS
C
C     Bidiagonalization, singular subspace, singular value
C     decomposition, singular values.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBU, JOBV
      INTEGER           INFO, IWARN, LDA, LDU, LDV, LDWORK, M, N, RANK
      DOUBLE PRECISION  RELTOL, THETA, TOL
C     .. Array Arguments ..
      LOGICAL           INUL(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), Q(*), U(LDU,*), V(LDV,*)
C     .. Local Scalars ..
      CHARACTER*1       JOBUY, JOBVY
      LOGICAL           ALL, LJOBUA, LJOBUS, LJOBVA, LJOBVS, QR, WANTU,
     $                  WANTV
      INTEGER           I, IHOUSH, IJ, ITAU, ITAUP, ITAUQ, J, JU, JV,
     $                  JWORK, K, LDW, LDY, MA, P, PP1, WRKOPT
      DOUBLE PRECISION  CS, SN, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBRD, DGEQRF, DLARTG, DLASET, DLASR,
     $                  MB04XY, MB04YD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      IWARN = 0
      INFO = 0
      P = MIN( M, N )
      K = MAX( M, N )
C
C     Determine whether U and/or V are/is to be computed.
C
      LJOBUA = LSAME( JOBU, 'A' )
      LJOBUS = LSAME( JOBU, 'S' )
      LJOBVA = LSAME( JOBV, 'A' )
      LJOBVS = LSAME( JOBV, 'S' )
      WANTU  = LJOBUA.OR.LJOBUS
      WANTV  = LJOBVA.OR.LJOBVS
      ALL = ( LJOBUA .AND. M.GT.N ) .OR. ( LJOBVA .AND. M.LT.N )
      QR  = M.GE.ILAENV( 6, 'DGESVD', 'N' // 'N', M, N, 0, 0 )
      IF ( QR.AND.WANTU ) THEN
         LDW = MAX( 2*N, N*( N + 1 )/2 )
      ELSE
         LDW = 0
      END IF
      IF ( WANTU.OR.WANTV ) THEN
         LDY = 8*P - 5
      ELSE
         LDY = 6*P - 3
      END IF
C
C     Test the input scalar arguments.
C
      IF( .NOT.WANTU .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTV .AND. .NOT.LSAME( JOBV, 'N' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( RANK.GT.P ) THEN
         INFO = -5
      ELSE IF( RANK.LT.0 .AND. THETA.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( ( .NOT.WANTU .AND. LDU.LT.1 )             .OR.
     $         (      WANTU .AND. LDU.LT.MAX( 1, M ) ) ) THEN
         INFO = -10
      ELSE IF( ( .NOT.WANTV .AND. LDV.LT.1 )             .OR.
     $         (      WANTV .AND. LDV.LT.MAX( 1, N ) ) ) THEN
         INFO = -12
      ELSE IF( LDWORK.LT.MAX( 1, LDW + MAX( 2*P + K, LDY ) ) ) THEN
         INFO = -18
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB04XD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( P.EQ.0 ) THEN
         IF ( RANK.GE.0 )
     $      THETA = ZERO
         RANK = 0
         RETURN
      END IF
C
C     Initializations.
C
      PP1 = P + 1
C
      IF ( ALL .AND. ( .NOT.QR ) ) THEN
C
         DO 20 I = 1, P
            INUL(I) = .FALSE.
   20    CONTINUE
C
         DO 40 I = PP1, K
            INUL(I) = .TRUE.
   40    CONTINUE
C
      ELSE
C
         DO 60 I = 1, K
            INUL(I) = .FALSE.
   60    CONTINUE
C
      END IF
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( QR ) THEN
C
C        1.a.: M is large enough than N; transform A into upper
C              triangular form R by Householder transformations.
C
C        Workspace: need 2*N;  prefer N + N*NB.
C
         ITAU = 1
         JWORK = ITAU + N
         CALL DGEQRF( M, N, A, LDA, DWORK(ITAU), DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
         WRKOPT = INT( DWORK(JWORK) )+JWORK-1
C
C        If (WANTU), store information on the Householder
C        transformations performed on the columns of A in N*(N+1)/2
C        extra storage locations DWORK(K), for K = 1,2,...,N*(N+1)/2.
C        (The first N locations store the scalar factors of Householder
C        transformations.)
C
C        Workspace: LDW = max(2*N, N*(N+1)/2).
C
         IF ( WANTU ) THEN
            IHOUSH = JWORK
            K = IHOUSH
            I = N
         ELSE
            K = 1
         END IF
C
         DO 100 J = 1, N - 1
            IF ( WANTU ) THEN
               I = I - 1
               CALL DCOPY( I, A(J+1,J), 1, DWORK(K), 1 )
               K = K + I
            END IF
C
            DO 80 IJ = J + 1, N
               A(IJ,J) = ZERO
   80       CONTINUE
C
  100    CONTINUE
C
         MA = N
         WRKOPT = MAX( WRKOPT, K )
      ELSE
C
C        Workspace: LDW = 0.
C
         K  = 1
         MA = M
         WRKOPT = 1
      END IF
C
C     1.b.: Transform A (or R) into bidiagonal form Q using Householder
C           transformations.
C
C     Workspace: need   LDW + 2*min(M,N) + max(M,N);
C                prefer LDW + 2*min(M,N) + (M+N)*NB.
C
      ITAUQ = K
      ITAUP = ITAUQ + P
      JWORK = ITAUP + P
      CALL DGEBRD( MA, N, A, LDA, Q, Q(PP1), DWORK(ITAUQ),
     $             DWORK(ITAUP), DWORK(JWORK), LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     1.c.: Initialize U (if WANTU) and V (if WANTV) with the identity
C           matrix.
C
      IF ( WANTU ) THEN
         IF ( ALL ) THEN
            JU = M
         ELSE
            JU = P
         END IF
         CALL DLASET( 'Full', M, JU, ZERO, ONE, U, LDU )
         JOBUY = 'U'
      ELSE
         JOBUY = 'N'
      END IF
      IF ( WANTV ) THEN
         IF ( ALL ) THEN
            JV = N
         ELSE
            JV = P
         END IF
         CALL DLASET( 'Full', N, JV, ZERO, ONE, V, LDV )
         JOBVY = 'U'
      ELSE
         JOBVY = 'N'
      END IF
C
C     If the matrix is lower bidiagonal, rotate to be upper bidiagonal
C     by applying Givens rotations on the left.
C
      IF ( M.LT.N ) THEN
C
         DO 120 I = 1, P - 1
            CALL DLARTG( Q(I), Q(P+I), CS, SN, TEMP )
            Q(I)   = TEMP
            Q(P+I) = SN*Q(I+1)
            Q(I+1) = CS*Q(I+1)
            IF ( WANTU ) THEN
C
C              Workspace: LDW + 4*min(M,N) - 2.
C
               DWORK(JWORK+I-1) = CS
               DWORK(JWORK+P+I-2) = SN
            END IF
  120    CONTINUE
C
C        Update left singular vectors if desired.
C
         IF( WANTU )
     $      CALL DLASR( 'Right', 'Variable pivot', 'Forward', M, JU,
     $                  DWORK(JWORK), DWORK(JWORK+P-1), U, LDU )
C
      END IF
C
C     Step 2: Partial diagonalization phase.
C             -----------------------------
C             Diagonalize the bidiagonal Q partially until convergence
C             to  the desired left and/or right singular subspace.
C
C              Workspace: LDW + 8*min(M,N) - 5, if WANTU or WANTV;
C              Workspace: LDW + 6*min(M,N) - 3, if JOBU = JOBV = 'N'.
C
      CALL MB04YD( JOBUY, JOBVY, M, N, RANK, THETA, Q, Q(PP1), U, LDU,
     $             V, LDV, INUL, TOL, RELTOL, DWORK(JWORK),
     $             LDWORK-JWORK+1, IWARN, INFO )
      IF ( WANTU.OR.WANTV ) THEN
         WRKOPT = MAX( WRKOPT, JWORK - 6 + 8*P )
      ELSE
         WRKOPT = MAX( WRKOPT, JWORK - 4 + 6*P )
      END IF
      IF ( INFO.GT.0 )
     $   RETURN
C
C     Step 3: Back transformation phase.
C             -------------------------
C     3.a.: Apply the Householder transformations of the bidiagonaliza-
C           tion onto the base vectors associated with the desired
C           bidiagonal submatrices.
C
C           Workspace: LDW + 2*min(M,N).
C
      CALL MB04XY( JOBU, JOBV, MA, N, A, LDA, DWORK(ITAUQ),
     $             DWORK(ITAUP), U, LDU, V, LDV, INUL, INFO )
C
C     3.b.: If A was reduced to upper triangular form R and JOBU = 'A'
C           or JOBU = 'S' apply the Householder transformations of the
C           triangularization of A onto the desired base vectors.
C
      IF ( QR.AND.WANTU ) THEN
         IF ( ALL ) THEN
C
            DO 140 I = PP1, M
               INUL(I) = .TRUE.
  140       CONTINUE
C
         END IF
         K = IHOUSH
         I = N
C
         DO 160 J = 1, N - 1
            I = I - 1
            CALL DCOPY( I, DWORK(K), 1, A(J+1,J), 1 )
            K = K + I
  160    CONTINUE
C
C        Workspace: MIN(M,N) + 1.
C
         JWORK = PP1
         CALL MB04XY( JOBU, 'No V', M, N, A, LDA, DWORK(ITAU),
     $                DWORK(ITAU), U, LDU, DWORK(JWORK), 1, INUL, INFO )
         WRKOPT = MAX( WRKOPT, PP1 )
      END IF
C
C     Set the optimal workspace.
C
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of MB04XD ***
      END
