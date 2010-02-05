      SUBROUTINE MB02ND( M, N, L, RANK, THETA, C, LDC, X, LDX, Q, INUL,
     $                   TOL, RELTOL, IWORK, DWORK, LDWORK, BWORK,
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
C     To solve the Total Least Squares (TLS) problem using a Partial
C     Singular Value Decomposition (PSVD) approach.
C     The TLS problem assumes an overdetermined set of linear equations
C     AX = B, where both the data matrix A as well as the observation
C     matrix B are inaccurate. The routine also solves determined and
C     underdetermined sets of equations by computing the minimum norm
C     solution.
C     It is assumed that all preprocessing measures (scaling, coordinate
C     transformations, whitening, ... ) of the data have been performed
C     in advance.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in the data matrix A and the
C             observation matrix B.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns in the data matrix A.  N >= 0.
C
C     L       (input) INTEGER
C             The number of columns in the observation matrix B.
C             L >= 0.
C
C     RANK    (input/output) INTEGER
C             On entry, if RANK < 0, then the rank of the TLS
C             approximation [A+DA|B+DB] (r say) is computed by the
C             routine.
C             Otherwise, RANK must specify the value of r.
C             RANK <= min(M,N).
C             On exit, if RANK < 0 on entry and INFO = 0, then RANK
C             contains the computed rank of the TLS approximation
C             [A+DA|B+DB].
C             Otherwise, the user-supplied value of RANK may be
C             changed by the routine on exit if the RANK-th and the
C             (RANK+1)-th singular values of C = [A|B] are considered
C             to be equal, or if the upper triangular matrix F (as
C             defined in METHOD) is (numerically) singular.
C
C     THETA   (input/output) DOUBLE PRECISION
C             On entry, if RANK < 0, then the rank of the TLS
C             approximation [A+DA|B+DB] is computed using THETA as
C             (min(M,N+L) - d), where d is the number of singular
C             values of [A|B] <= THETA. THETA >= 0.0.
C             Otherwise, THETA is an initial estimate (t say) for
C             computing a lower bound on the RANK largest singular
C             values of [A|B]. If THETA < 0.0 on entry however, then
C             t is computed by the routine.
C             On exit, if RANK >= 0 on entry, then THETA contains the
C             computed bound such that precisely RANK singular values
C             of C = [A|B] are greater than THETA + TOL.
C             Otherwise, THETA is unchanged.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N+L)
C             On entry, the leading M-by-(N+L) part of this array must
C             contain the matrices A and B. Specifically, the first N
C             columns must contain the data matrix A and the last L
C             columns the observation matrix B (right-hand sides).
C             On exit, if INFO = 0, the first N+L components of the
C             columns of this array whose index i corresponds with
C             INUL(i) = .TRUE., are the possibly transformed (N+L-RANK)
C             base vectors of the right singular subspace corresponding
C             to the singular values of C = [A|B] which are less than or
C             equal to THETA. Specifically, if L = 0, or if RANK = 0 and
C             IWARN <> 2, these vectors are indeed the base vectors
C             above. Otherwise, these vectors form the matrix V2,
C             transformed as described in Step 4 of the PTLS algorithm
C             (see METHOD). The TLS solution is computed from these
C             vectors. The other columns of array C contain no useful
C             information.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= max(1,M,N+L).
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,L)
C             If INFO = 0, the leading N-by-L part of this array
C             contains the solution X to the TLS problem specified by
C             A and B.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= max(1,N).
C
C     Q       (output) DOUBLE PRECISION array, dimension
C             (max(1,2*min(M,N+L)-1))
C             This array contains the partially diagonalized bidiagonal
C             matrix J computed from C, at the moment that the desired
C             singular subspace has been found. Specifically, the
C             leading p = min(M,N+L) entries of Q contain the diagonal
C             elements q(1),q(2),...,q(p) and the entries Q(p+1),Q(p+2),
C             ...,Q(2*p-1) contain the superdiagonal elements e(1),e(2),
C             ...,e(p-1) of J.
C
C     INUL    (output) LOGICAL array, dimension (N+L)
C             The indices of the elements of this array with value
C             .TRUE. indicate the columns in C containing the base
C             vectors of the right singular subspace of C from which
C             the TLS solution has been computed.
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
C     IWORK   INTEGER array, dimension (N+2*L)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, and DWORK(2) returns the reciprocal of the
C             condition number of the matrix F.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK = max(2, max(M,N+L) + 2*min(M,N+L),
C                          min(M,N+L) + LW + max(6*(N+L)-5,
C                                                L*L+max(N+L,3*L)),
C             where
C             LW = (N+L)*(N+L-1)/2,  if M >= N+L,
C             LW = M*(N+L-(M-1)/2),  if M <  N+L.
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension (N+L)
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warnings;
C             = 1:  if the rank of matrix C has been lowered because a
C                   singular value of multiplicity greater than 1 was
C                   found;
C             = 2:  if the rank of matrix C has been lowered because the
C                   upper triangular matrix F is (numerically) singular.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the maximum number of QR/QL iteration steps
C                   (30*MIN(M,N)) has been exceeded;
C             = 2:  if the computed rank of the TLS approximation
C                   [A+DA|B+DB] exceeds MIN(M,N). Try increasing the
C                   value of THETA or set the value of RANK to min(M,N).
C
C     METHOD
C
C     The method used is the Partial Total Least Squares (PTLS) approach
C     proposed by Van Huffel and Vandewalle [5].
C
C     Let C = [A|B] denote the matrix formed by adjoining the columns of
C     B to the columns of A on the right.
C
C     Total Least Squares (TLS) definition:
C     -------------------------------------
C
C       Given matrices A and B, find a matrix X satisfying
C
C            (A + DA) X = B + DB,
C
C       where A and DA are M-by-N matrices, B and DB are M-by-L matrices
C       and X is an N-by-L matrix.
C       The solution X must be such that the Frobenius norm of [DA|DB]
C       is a minimum and each column of B + DB is in the range of
C       A + DA. Whenever the solution is not unique, the routine singles
C       out the minimum norm solution X.
C
C     Let V denote the right singular subspace of C. Since the TLS
C     solution can be computed from any orthogonal basis of the subspace
C     of V corresponding to the smallest singular values of C, the
C     Partial Singular Value Decomposition (PSVD) can be used instead of
C     the classical SVD. The dimension of this subspace of V may be
C     determined by the rank of C or by an upper bound for those
C     smallest singular values.
C
C     The PTLS algorithm proceeds as follows (see [2 - 5]):
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C      (a) If M is large enough than N + L, transform C into upper
C          triangular form R by Householder transformations.
C      (b) Transform C (or R) into upper bidiagonal form
C          (p = min(M,N+L)):
C
C                     |q(1) e(1)  0   ...  0   |
C                (0)  | 0   q(2) e(2)      .   |
C               J   = | .                  .   |
C                     | .                e(p-1)|
C                     | 0             ... q(p) |
C
C          if M >= N + L, or lower bidiagonal form:
C
C                     |q(1)  0    0   ...  0     0   |
C                (0)  |e(1) q(2)  0        .     .   |
C               J   = | .                  .     .   |
C                     | .                 q(p)   .   |
C                     | 0             ... e(p-1) q(p)|
C
C          if M < N + L, using Householder transformations.
C          In the second case, transform the matrix to the upper
C          bidiagonal form by applying Givens rotations.
C      (c) Initialize the right singular base matrix with the identity
C          matrix.
C
C     Step 2: Partial diagonalization phase
C             -----------------------------
C     If the upper bound THETA is not given, then compute THETA such
C     that precisely p - RANK singular values (p=min(M,N+L)) of the
C     bidiagonal matrix are less than or equal to THETA, using a
C     bisection method [5]. Diagonalize the given bidiagonal matrix J
C     partially, using either QL iterations (if the upper left diagonal
C     element of the considered bidiagonal submatrix is smaller than the
C     lower right diagonal element) or QR iterations, such that J is
C     split into unreduced bidiagonal submatrices whose singular values
C     are either all larger than THETA or are all less than or equal
C     to THETA. Accumulate the Givens rotations in V.
C
C     Step 3: Back transformation phase
C             -------------------------
C     Apply the Householder transformations of Step 1(b) onto the base
C     vectors of V associated with the bidiagonal submatrices with all
C     singular values less than or equal to THETA.
C
C     Step 4: Computation of F and Y
C             ----------------------
C     Let V2 be the matrix of the columns of V corresponding to the
C     (N + L - RANK) smallest singular values of C.
C     Compute with Householder transformations the matrices F and Y
C     such that:
C
C                       |VH   Y|
C              V2 x Q = |      |
C                       |0    F|
C
C     where Q is an orthogonal matrix, VH is an N-by-(N-RANK) matrix,
C     Y is an N-by-L matrix and F is an L-by-L upper triangular matrix.
C     If F is singular, then reduce the value of RANK by one and repeat
C     Steps 2, 3 and 4.
C
C     Step 5: Computation of the TLS solution
C             -------------------------------
C     If F is non-singular then the solution X is obtained by solving
C     the following equations by forward elimination:
C
C              X F = -Y.
C
C     Notes:
C     If RANK is lowered in Step 4, some additional base vectors must
C     be computed in Step 2. The additional computations are kept to
C     a minimum.
C     If RANK is lowered in Step 4 but the multiplicity of the RANK-th
C     singular value is larger than 1, then the value of RANK is further
C     lowered with its multiplicity defined by the parameter TOL. This
C     is done at the beginning of Step 2 by calling SLICOT Library
C     routine MB03MD (from MB04YD), which estimates THETA using a
C     bisection method. If F in Step 4 is singular, then the computed
C     solution is infinite and hence does not satisfy the second TLS
C     criterion (see TLS definition). For these cases, Golub and
C     Van Loan [1] claim that the TLS problem has no solution. The
C     properties of these so-called nongeneric problems are described
C     in [6] and the TLS computations are generalized in order to solve
C     them. As proven in [6], the proposed generalization satisfies the
C     TLS criteria for any number L of observation vectors in B provided
C     that, in addition, the solution | X| is constrained to be
C                                     |-I|
C     orthogonal to all vectors of the form |w| which belong to the
C                                           |0|
C     space generated by the columns of the submatrix |Y|.
C                                                     |F|
C
C     REFERENCES
C
C     [1] Golub, G.H. and Van Loan, C.F.
C         An Analysis of the Total Least-Squares Problem.
C         SIAM J. Numer. Anal., 17, pp. 883-893, 1980.
C
C     [2] Van Huffel, S., Vandewalle, J. and Haegemans, A.
C         An Efficient and Reliable Algorithm for Computing the
C         Singular Subspace of a Matrix Associated with its Smallest
C         Singular Values.
C         J. Comput. and Appl. Math., 19, pp. 313-330, 1987.
C
C     [3] Van Huffel, S.
C         Analysis of the Total Least Squares Problem and its Use in
C         Parameter Estimation.
C         Doctoral dissertation, Dept. of Electr. Eng., Katholieke
C         Universiteit Leuven, Belgium, June 1987.
C
C     [4] Chan, T.F.
C         An Improved Algorithm for Computing the Singular Value
C         Decomposition.
C         ACM TOMS, 8, pp. 72-83, 1982.
C
C     [5] Van Huffel, S. and Vandewalle, J.
C         The Partial Total Least Squares Algorithm.
C         J. Comput. Appl. Math., 21, pp. 333-341, 1988.
C
C     [6] Van Huffel, S. and Vandewalle, J.
C         Analysis and Solution of the Nongeneric Total Least Squares
C         Problem.
C         SIAM J. Matr. Anal. and Appl., 9, pp. 360-372, 1988.
C
C     NUMERICAL ASPECTS
C
C     The computational efficiency of the PTLS algorithm compared with
C     the classical TLS algorithm (see [2 - 5]) is obtained by making
C     use of PSVD (see [1]) instead of performing the entire SVD.
C     Depending on the gap between the RANK-th and the (RANK+1)-th
C     singular values of C, the number (N + L - RANK) of base vectors to
C     be computed with respect to the column dimension (N + L) of C and
C     the desired accuracy RELTOL, the algorithm used by this routine is
C     approximately twice as fast as the classical TLS algorithm at the
C     expense of extra storage requirements, namely:
C       (N + L) x (N + L - 1)/2  if M >= N + L or
C       M x (N + L - (M - 1)/2)  if M <  N + L.
C     This is because the Householder transformations performed on the
C     rows of C in the bidiagonalization phase (see Step 1) must be kept
C     until the end (Step 5).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB02BD by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     June 30, 1997, Oct. 19, 2003, Feb. 15, 2004.
C
C     KEYWORDS
C
C     Least-squares approximation, singular subspace, singular value
C     decomposition, singular values, total least-squares.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, IWARN, L, LDC, LDWORK, LDX, M, N, RANK
      DOUBLE PRECISION  RELTOL, THETA, TOL
C     .. Array Arguments ..
      LOGICAL           BWORK(*), INUL(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  C(LDC,*), DWORK(*), Q(*), X(LDX,*)
C     .. Local Scalars ..
      LOGICAL           LFIRST, SUFWRK
      INTEGER           I, I1, IFAIL, IHOUSH, IJ, IOFF, ITAUP, ITAUQ,
     $                  IWARM, J, J1, JF, JV, JWORK, K, KF, KJ, LDF, LW,
     $                  MC, MJ, MNL, N1, NJ, NL, P, WRKOPT
      DOUBLE PRECISION  CS, EPS, FIRST, FNORM, HH, INPROD, RCOND, SN,
     $                  TEMP
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANTR
      EXTERNAL          DLAMCH, DLANGE, DLANTR, ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBRD, DGEQRF, DGERQF, DLARF, DLARFG,
     $                  DLARTG, DLASET, DORMBR, DORMRQ, DTRCON, DTRSM,
     $                  MB04YD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      IWARN = 0
      INFO  = 0
      NL = N + L
      K  = MAX( M, NL )
      P  = MIN( M, NL )
      IF ( M.GE.NL ) THEN
         LW = ( NL*( NL - 1 ) )/2
      ELSE
         LW = M*NL - ( M*( M - 1 ) )/2
      END IF
      JV = P + LW + MAX( 6*NL - 5, L*L + MAX( NL, 3*L ) )
C
C     Test the input scalar arguments.
C
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 ) THEN
         INFO = -3
      ELSE IF( RANK.GT.MIN( M, N ) ) THEN
         INFO = -4
      ELSE IF( ( RANK.LT.0 ) .AND. ( THETA.LT.ZERO ) ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDWORK.LT.MAX( 2, K + 2*P, JV ) ) THEN
         INFO = -16
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB02ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( M, NL ).EQ.0 ) THEN
         IF ( M.EQ.0 ) THEN
            CALL DLASET( 'Full', NL, NL, ZERO, ONE, C, LDC )
            CALL DLASET( 'Full', N, L, ZERO, ZERO, X, LDX )
C
            DO 10 I = 1, NL
               INUL(I) = .TRUE.
   10       CONTINUE
C
         END IF
         IF ( RANK.GE.0 )
     $      THETA = ZERO
         RANK = 0
         DWORK(1) = TWO
         DWORK(2) = ONE
         RETURN
      END IF
C
      WRKOPT = 2
      N1 = N + 1
C
      EPS = DLAMCH( 'Precision' )
      LFIRST = .TRUE.
C
C     Initializations.
C
      DO 20 I = 1, P
         INUL(I)  = .FALSE.
         BWORK(I) = .FALSE.
   20 CONTINUE
C
      DO 40 I = P + 1, NL
         INUL(I)  = .TRUE.
         BWORK(I) = .FALSE.
   40 CONTINUE
C
C     Subroutine MB02ND solves a set of linear equations by a Total
C     Least Squares Approximation, based on the Partial SVD.
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C     1.a): If M is large enough than N+L, transform C into upper
C           triangular form R by Householder transformations.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( M.GE.MAX( NL,
     $               ILAENV( 6, 'DGESVD', 'N' // 'N', M, NL, 0, 0 ) ) )
     $      THEN
C
C        Workspace: need   2*(N+L),
C                   prefer N+L + (N+L)*NB.
C
         ITAUQ = 1
         JWORK = ITAUQ + NL
         CALL DGEQRF( M, NL, C, LDC, DWORK(ITAUQ), DWORK(JWORK),
     $                LDWORK-JWORK+1, IFAIL )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
         IF ( NL.GT.1 )
     $      CALL DLASET( 'Lower', NL-1, NL-1, ZERO, ZERO, C(2,1), LDC )
         MNL = NL
      ELSE
         MNL = M
      END IF
C
C     1.b): Transform C (or R) into bidiagonal form Q using Householder
C           transformations.
C     Workspace: need   2*min(M,N+L) + max(M,N+L),
C                prefer 2*min(M,N+L) + (M+N+L)*NB.
C
      ITAUP = 1
      ITAUQ = ITAUP + P
      JWORK = ITAUQ + P
      CALL DGEBRD( MNL, NL, C, LDC, Q, Q(P+1), DWORK(ITAUQ),
     $             DWORK(ITAUP), DWORK(JWORK), LDWORK-JWORK+1, IFAIL )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     If the matrix is lower bidiagonal, rotate to be upper bidiagonal
C     by applying Givens rotations on the left.
C
      IF ( M.LT.NL ) THEN
         IOFF = 0
C
         DO 60 I = 1, P - 1
            CALL DLARTG( Q(I), Q(P+I), CS, SN, TEMP )
            Q(I)   = TEMP
            Q(P+I) = SN*Q(I+1)
            Q(I+1) = CS*Q(I+1)
   60    CONTINUE
C
      ELSE
         IOFF = 1
      END IF
C
C     Store the Householder transformations performed onto the rows of C
C     in the extra storage locations DWORK(IHOUSH).
C     Workspace: need   LDW = min(M,N+L) + (N+L)*(N+L-1)/2, if M >= N+L,
C                       LDW = min(M,N+L) + M*(N+L-(M-1)/2), if M <  N+L;
C                prefer LDW = min(M,N+L) + (N+L)**2,        if M >= N+L,
C                       LDW = min(M,N+L) + M*(N+L),         if M <  N+L.
C
      IHOUSH = ITAUQ
      MC = NL - IOFF
      KF = IHOUSH + P*NL
      SUFWRK = LDWORK.GE.( KF + MAX( 6*(N+L)-5,
     $                               NL**2 + MAX( NL, 3*L ) - 1 ) )
      IF ( SUFWRK ) THEN
C
C        Enough workspace for a fast algorithm.
C
         CALL DLACPY( 'Upper', P, NL, C, LDC, DWORK(IHOUSH), P )
         KJ = KF
         WRKOPT = MAX( WRKOPT, KF - 1 )
      ELSE
C
C        Not enough workspace for a fast algorithm.
C
         KJ = IHOUSH
C
         DO 80 NJ = 1, MIN( P, MC )
            J = MC - NJ + 1
            CALL DCOPY( J, C(NJ,NJ+IOFF), LDC, DWORK(KJ), 1 )
            KJ = KJ + J
   80    CONTINUE
C
      END IF
C
C     1.c): Initialize the right singular base matrix V with the
C           identity matrix (V overwrites C).
C
      CALL DLASET( 'Full', NL, NL, ZERO, ONE, C, LDC )
      JV = KJ
      IWARM = 0
C
C     REPEAT
C
C     Compute the Householder matrix Q and matrices F and Y such that
C     F is nonsingular.
C
C     Step 2: Partial diagonalization phase.
C             -----------------------------
C     Diagonalize the bidiagonal Q partially until convergence to
C     the desired right singular subspace.
C     Workspace: LDW + 6*(N+L)-5.
C
  100 CONTINUE
      JWORK = JV
      CALL MB04YD( 'No U', 'Update V', P, NL, RANK, THETA, Q, Q(P+1),
     $             DUMMY, 1, C, LDC, INUL, TOL, RELTOL, DWORK(JWORK),
     $             LDWORK-JWORK+1, IWARN, INFO )
      WRKOPT = MAX( WRKOPT, JWORK + 6*NL - 6 )
C
      IWARN = MAX( IWARN, IWARM )
      IF ( INFO.GT.0 )
     $   RETURN
C
C     Set pointers to the selected base vectors in the right singular
C     matrix of C.
C
      K = 0
C
      DO 120 I = 1, NL
         IF ( INUL(I) ) THEN
            K = K + 1
            IWORK(K) = I
         END IF
  120 CONTINUE
C
      IF ( K.LT.L ) THEN
C
C        Rank of the TLS approximation is larger than min(M,N).
C
         INFO = 2
         RETURN
      END IF
C
C     Step 3: Back transformation phase.
C             -------------------------
C     Apply in backward order the Householder transformations (stored
C     in DWORK(IHOUSH)) performed onto the rows of C during the
C     bidiagonalization phase, to the selected base vectors (specified
C     by INUL(I) = .TRUE.). Already transformed vectors are those for
C     which BWORK(I) = .TRUE..
C
      KF = K
      IF ( SUFWRK.AND.LFIRST ) THEN
C
C        Enough workspace for a fast algorithm and first pass.
C
         IJ = JV
C
         DO 140 J = 1, K
            CALL DCOPY (NL, C(1,IWORK(J)), 1, DWORK(IJ), 1 )
            IJ = IJ + NL
  140    CONTINUE
C
C        Workspace: need   LDW + (N+L)*K + K,
C                   prefer LDW + (N+L)*K + K*NB.
C
         IJ = JV
         JWORK = IJ + NL*K
         CALL DORMBR( 'P vectors', 'Left', 'No transpose', NL, K,
     $                MNL, DWORK(IHOUSH), P, DWORK(ITAUP), DWORK(IJ),
     $                NL, DWORK(JWORK), LDWORK-JWORK+1, IFAIL )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
         DO 160 I = 1, NL
            IF ( INUL(I) .AND. ( .NOT. BWORK(I) ) )
     $         BWORK(I) = .TRUE.
  160    CONTINUE
C
      ELSE
C
C        Not enough workspace for a fast algorithm or subsequent passes.
C
         DO 180 I = 1, NL
            IF ( INUL(I) .AND. ( .NOT. BWORK(I) ) ) THEN
               KJ = JV
C
               DO 170 NJ = MIN( P, MC ), 1, -1
                  J = MC - NJ + 1
                  KJ = KJ - J
                  FIRST = DWORK(KJ)
                  DWORK(KJ) = ONE
                  CALL DLARF( 'Left', J, 1, DWORK(KJ), 1,
     $                        DWORK(ITAUP+NJ-1), C(NJ+IOFF,I), LDC,
     $                        DWORK(JWORK) )
                  DWORK(KJ) = FIRST
  170          CONTINUE
C
               BWORK(I) = .TRUE.
            END IF
  180    CONTINUE
      END IF
C
      IF ( RANK.LE.0 )
     $   RANK = 0
      IF ( MIN( RANK, L ).EQ.0 ) THEN
         IF ( SUFWRK.AND.LFIRST )
     $      CALL DLACPY( 'Full', NL, K, DWORK(JV), NL, C, LDC )
         DWORK(1) = WRKOPT
         DWORK(2) = ONE
         RETURN
      END IF
C
C     Step 4: Compute matrices F and Y
C             ------------------------
C             using Householder transformation Q.
C
C     Compute the orthogonal matrix Q (in factorized form) and the
C     matrices F and Y using RQ factorization. It is assumed that,
C     generically, the last L rows of V2 matrix have full rank.
C     The code could not be the most efficient when RANK has been
C     lowered, because the already created zero pattern of the last
C     L rows of V2 matrix is not exploited.
C
      IF ( SUFWRK.AND.LFIRST ) THEN
C
C        Enough workspace for a fast algorithm and first pass.
C        Workspace: need   LDW1 + 2*L,
C                   prefer LDW1 + L + L*NB, where
C                          LDW1 = LDW + (N+L)*K;
C
         ITAUQ = JWORK
         JWORK = ITAUQ + L
         CALL DGERQF( L, K, DWORK(JV+N), NL, DWORK(ITAUQ), DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Workspace: need   LDW1 + N+L,
C                   prefer LDW1 + L + N*NB.
C
         CALL DORMRQ( 'Right', 'Transpose', N, K, L, DWORK(JV+N), NL,
     $                DWORK(ITAUQ), DWORK(JV), NL, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
         JF = JV + NL*(K-L) + N
         LDF = NL
         JWORK = JF + LDF*L - N
         CALL DLASET( 'Full', L, K-L, ZERO, ZERO, DWORK(JV+N), LDF )
         IF ( L.GT.1 )
     $      CALL DLASET( 'Lower', L-1, L-1, ZERO, ZERO, DWORK(JF+1),
     $                   LDF )
         IJ = JV
C
         DO 200 J = 1, K
            CALL DCOPY( NL, DWORK(IJ), 1, C(1,IWORK(J)), 1 )
            IJ = IJ + NL
  200    CONTINUE
C
      ELSE
C
C        Not enough workspace for a fast algorithm or subsequent passes.
C        Workspace: LDW2 + N+L, where LDW2 = LDW + L*L.
C
         I  = NL
         JF = JV
         LDF = L
         JWORK  = JF + LDF*L
         WRKOPT = MAX( WRKOPT, JWORK+NL-1 )
C
C        WHILE ( ( K >= 1 ) .AND. ( I > N ) ) DO
  220    CONTINUE
         IF ( ( K.GE.1 ) .AND. ( I.GT.N ) ) THEN
C
            DO 240 J = 1, K
               DWORK(JWORK+J-1) = C(I,IWORK(J))
  240       CONTINUE
C
C           Compute Householder transformation.
C
            CALL DLARFG( K, DWORK(JWORK+K-1), DWORK(JWORK), 1, TEMP )
            C(I,IWORK(K)) = DWORK(JWORK+K-1)
            IF ( TEMP.NE.ZERO ) THEN
C
C              Apply Householder transformation onto the selected base
C              vectors.
C
               DO 300 I1 = 1, I - 1
                  INPROD = C(I1,IWORK(K))
C
                  DO 260 J = 1, K - 1
                     INPROD = INPROD + DWORK(JWORK+J-1)*C(I1,IWORK(J))
  260             CONTINUE
C
                  HH = INPROD*TEMP
                  C(I1,IWORK(K)) = C(I1,IWORK(K)) - HH
C
                  DO 280 J = 1, K - 1
                     J1 = IWORK(J)
                     C(I1,J1) = C(I1,J1) - DWORK(JWORK+J-1)*HH
                     C(I,J1)  = ZERO
  280             CONTINUE
C
  300          CONTINUE
C
            END IF
            CALL DCOPY( I-N, C(N1,IWORK(K)), 1, DWORK(JF+(I-N-1)*L), 1 )
            K = K - 1
            I = I - 1
            GO TO 220
         END IF
C        END WHILE 220
      END IF
C
C     Estimate the reciprocal condition number of the matrix F.
C     If F singular, lower the rank of the TLS approximation.
C     Workspace: LDW1 + 3*L or
C                LDW2 + 3*L.
C
      CALL DTRCON( '1-norm', 'Upper', 'Non-unit', L, DWORK(JF), LDF,
     $             RCOND, DWORK(JWORK), IWORK(KF+1), INFO )
      WRKOPT = MAX( WRKOPT, JWORK + 3*L - 1 )
C
      DO 320 J = 1, L
         CALL DCOPY( N, C(1,IWORK(KF-L+J)), 1, X(1,J), 1 )
  320 CONTINUE
C
      FNORM = DLANTR( '1-norm', 'Upper', 'Non-unit', L, L, DWORK(JF),
     $                LDF, DWORK(JWORK) )
      IF ( RCOND.LE.EPS*FNORM ) THEN
         RANK = RANK - 1
         GO TO 340
      END IF
      IF ( FNORM.LE.EPS*DLANGE( '1-norm', N, L, X, LDX,
     $                             DWORK(JWORK) ) ) THEN
         RANK = RANK - L
         GO TO 340
      ELSE
         GO TO 400
      END IF
C
  340    CONTINUE
         IWARM = 2
         THETA = -ONE
         IF ( SUFWRK.AND.LFIRST ) THEN
C
C           Rearrange the stored Householder transformations for
C           subsequent passes, taking care to avoid overwriting.
C
            IF ( P.LT.NL ) THEN
               KJ = IHOUSH + NL*(NL - 1)
               MJ = IHOUSH +  P*(NL - 1)
C
               DO 360 NJ = 1, NL
                  DO 350 J = P - 1, 0, -1
                     DWORK(KJ+J) = DWORK(MJ+J)
  350             CONTINUE
                  KJ = KJ - NL
                  MJ = MJ - P
  360          CONTINUE
C
            END IF
            KJ = IHOUSH
            MJ = IHOUSH + NL*IOFF
C
            DO 380 NJ = 1, MIN( P, MC )
               DO 370 J = 0, MC - NJ
                  DWORK(KJ) = DWORK(MJ+J*P)
                  KJ = KJ + 1
  370          CONTINUE
               MJ = MJ + NL + 1
  380       CONTINUE
C
            JV = KJ
            LFIRST = .FALSE.
         END IF
         GO TO 100
C     UNTIL ( F nonsingular, i.e., RCOND.GT.EPS*FNORM or
C                                  FNORM.GT.EPS*norm(Y) )
  400 CONTINUE
C
C     Step 5: Compute TLS solution.
C             --------------------
C     Solve X F = -Y  by forward elimination  (F is upper triangular).
C
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', N, L,
     $            -ONE, DWORK(JF), LDF, X, LDX )
C
C     Set the optimal workspace and reciprocal condition number of F.
C
      DWORK(1) = WRKOPT
      DWORK(2) = RCOND
C
      RETURN
C *** Last line of MB02ND ***
      END
