      SUBROUTINE FB01RD( JOBK, MULTBQ, N, M, P, S, LDS, A, LDA, B,
     $                   LDB, Q, LDQ, C, LDC, R, LDR, K, LDK, TOL,
     $                   IWORK, DWORK, LDWORK, INFO )
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
C     To calculate a combined measurement and time update of one
C     iteration of the time-invariant Kalman filter. This update is
C     given for the square root covariance filter, using the condensed
C     observer Hessenberg form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBK    CHARACTER*1
C             Indicates whether the user wishes to compute the Kalman
C             filter gain matrix K  as follows:
C                                 i
C             = 'K':  K  is computed and stored in array K;
C                      i
C             = 'N':  K  is not required.
C                      i
C
C     MULTBQ  CHARACTER*1                    1/2
C             Indicates how matrices B  and Q    are to be passed to
C                                     i      i
C             the routine as follows:
C             = 'P':  Array Q is not used and the array B must contain
C                                    1/2
C                     the product B Q   ;
C                                  i i
C             = 'N':  Arrays B and Q must contain the matrices as
C                     described below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices S    and A.  N >= 0.
C                       i-1
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C              1/2
C             Q   .  M >= 0.
C              i
C
C     P       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C              1/2
C             R   .  P >= 0.
C              i
C
C     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N)
C             On entry, the leading N-by-N lower triangular part of this
C             array must contain S   , the square root (left Cholesky
C                                 i-1
C             factor) of the state covariance matrix at instant (i-1).
C             On exit, the leading N-by-N lower triangular part of this
C             array contains S , the square root (left Cholesky factor)
C                             i
C             of the state covariance matrix at instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,N).
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain A,
C             the state transition matrix of the discrete system in
C             lower observer Hessenberg form (e.g., as produced by
C             SLICOT Library Routine TB01ND).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain B ,
C                                                        1/2      i
C             the input weight matrix (or the product B Q    if
C                                                      i i
C             MULTBQ = 'P') of the discrete system at instant i.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,*)
C             If MULTBQ = 'N', then the leading M-by-M lower triangular
C                                              1/2
C             part of this array must contain Q   , the square root
C                                              i
C             (left Cholesky factor) of the input (process) noise
C             covariance matrix at instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C             Otherwise, Q is not referenced and can be supplied as a
C             dummy array (i.e., set parameter LDQ = 1 and declare this
C             array to be Q(1,1) in the calling program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= MAX(1,M) if MULTBQ = 'N';
C             LDQ >= 1        if MULTBQ = 'P'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain C,
C             the output weight matrix of the discrete system in lower
C             observer Hessenberg form (e.g., as produced by SLICOT
C             Library routine TB01ND).
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,P)
C             On entry, the leading P-by-P lower triangular part of this
C                                 1/2
C             array must contain R   , the square root (left Cholesky
C                                 i
C             factor) of the output (measurement) noise covariance
C             matrix at instant i.
C             On exit, the leading P-by-P lower triangular part of this
C                                    1/2
C             array contains (RINOV )   , the square root (left Cholesky
C                                  i
C             factor) of the covariance matrix of the innovations at
C             instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,P).
C
C     K       (output) DOUBLE PRECISION array, dimension (LDK,P)
C             If JOBK = 'K', and INFO = 0, then the leading N-by-P part
C             of this array contains K , the Kalman filter gain matrix
C                                     i
C             at instant i.
C             If JOBK = 'N', or JOBK = 'K' and INFO = 1, then the
C             leading N-by-P part of this array contains AK ,  a matrix
C                                                          i
C             related to the Kalman filter gain matrix at instant i (see
C                                                            -1/2
C             METHOD). Specifically, AK  = A P     C'(RINOV')    .
C                                      i      i|i-1        i
C
C     LDK     INTEGER
C             The leading dimension of array K.  LDK >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If JOBK = 'K', then TOL is used to test for near
C                                               1/2
C             singularity of the matrix (RINOV )   . If the user sets
C                                             i
C             TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then an implicitly computed, default
C             tolerance, defined by TOLDEF = P*P*EPS, is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C             Otherwise, TOL is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             where LIWORK = P if JOBK = 'K',
C             and   LIWORK = 1 otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.  If INFO = 0 and JOBK = 'K', DWORK(2) returns
C             an estimate of the reciprocal of the condition number
C                                        1/2
C             (in the 1-norm) of (RINOV )   .
C                                      i
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,N*(P+N+1),N*(P+N)+2*P,N*(N+M+2)),
C                           if JOBK = 'N';
C             LDWORK >= MAX(2,N*(P+N+1),N*(P+N)+2*P,N*(N+M+2),3*P),
C                           if JOBK = 'K'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C                                                        1/2
C             = 1:  if JOBK = 'K' and the matrix (RINOV )   is singular,
C                                                      i           1/2
C                   i.e., the condition number estimate of (RINOV )
C                                                                i
C                   (in the 1-norm) exceeds 1/TOL.  The matrices S, AK ,
C                               1/2                                   i
C                   and (RINOV )    have been computed.
C                             i
C
C     METHOD
C
C     The routine performs one recursion of the square root covariance
C     filter algorithm, summarized as follows:
C
C      |  1/2                     |     |         1/2          |
C      | R      0        C x S    |     | (RINOV )     0     0 |
C      |  i                   i-1 |     |       i              |
C      |             1/2          | T = |                      |
C      | 0      B x Q    A x S    |     |     AK       S     0 |
C      |         i   i        i-1 |     |       i       i      |
C
C           (Pre-array)                      (Post-array)
C
C     where T is unitary and (A,C) is in lower observer Hessenberg form.
C
C     An example of the pre-array is given below (where N = 6, P = 2
C     and M = 3):
C
C          |x   |      | x          |
C          |x x |      | x x        |
C          |____|______|____________|
C          |    | x x x| x x x      |
C          |    | x x x| x x x x    |
C          |    | x x x| x x x x x  |
C          |    | x x x| x x x x x x|
C          |    | x x x| x x x x x x|
C          |    | x x x| x x x x x x|
C
C     The corresponding state covariance matrix P      is then
C                                                i|i-1
C     factorized as
C
C         P     = S  S'
C          i|i-1   i  i
C
C     and one combined time and measurement update for the state X
C                                                                 i|i-1
C     is given by
C
C         X     = A X      + K (Y - C X     )
C          i+1|i     i|i-1    i  i     i|i-1
C
C                          -1/2
C     where K = AK (RINOV )     is the Kalman filter gain matrix and Y
C            i    i      i                                            i
C     is the observed output of the system.
C
C     The triangularization is done entirely via Householder
C     transformations exploiting the zero pattern of the pre-array.
C
C     REFERENCES
C
C     [1] Anderson, B.D.O. and Moore, J.B.
C         Optimal Filtering.
C         Prentice Hall, Englewood Cliffs, New Jersey, 1979.
C
C     [2] Van Dooren, P. and Verhaegen, M.H.G.
C         Condensed Forms for Efficient Time-Invariant Kalman Filtering.
C         SIAM J. Sci. Stat. Comp., 9. pp. 516-530, 1988.
C
C     [3] Verhaegen, M.H.G. and Van Dooren, P.
C         Numerical Aspects of Different Kalman Filter Implementations.
C         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986.
C
C     [4] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G.
C         Algorithm 675: FORTRAN Subroutines for Computing the Square
C         Root Covariance Filter and Square Root Information Filter in
C         Dense or Hessenberg Forms.
C         ACM Trans. Math. Software, 15, pp. 243-256, 1989.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires
C
C            3    2                           2         3
C     1/6 x N  + N x (3/2 x P + M) + 2 x N x P + 2/3 x P
C
C     operations and is backward stable (see [3]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01FD by M. Vanbegin,
C     P. Van Dooren, and M.H.G. Verhaegen.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, February 14, 2004.
C
C     KEYWORDS
C
C     Kalman filtering, observer Hessenberg form, optimal filtering,
C     orthogonal transformation, recursive estimation, square-root
C     covariance filtering, square-root filtering.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBK, MULTBQ
      INTEGER           INFO, LDA, LDB, LDC, LDK, LDQ, LDR, LDS, LDWORK,
     $                  M, N, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*),
     $                  K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
C     .. Local Scalars ..
      LOGICAL           LJOBK, LMULTB
      INTEGER           I, II, ITAU, JWORK, N1, PL, PN, WRKOPT
      DOUBLE PRECISION  RCOND
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLACPY, DTRMM, DTRMV, MB02OD, MB04JD,
     $                  MB04LD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      PN = P + N
      N1 = MAX( 1, N )
      INFO = 0
      LJOBK  = LSAME( JOBK, 'K' )
      LMULTB = LSAME( MULTBQ, 'P' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBK .AND. .NOT.LSAME( JOBK, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LMULTB .AND. .NOT.LSAME( MULTBQ, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDS.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -11
      ELSE IF( LDQ.LT.1 .OR. ( .NOT.LMULTB .AND. LDQ.LT.M ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDR.LT.MAX( 1, P ) ) THEN
         INFO = -17
      ELSE IF( LDK.LT.N1 ) THEN
         INFO = -19
      ELSE IF( ( LJOBK .AND. LDWORK.LT.MAX( 2, PN*N + N, PN*N + 2*P,
     $                                      N*(N + M + 2), 3*P ) ) .OR.
     $    ( .NOT.LJOBK .AND. LDWORK.LT.MAX( 1, PN*N + N, PN*N + 2*P,
     $                                      N*(N + M + 2) ) ) ) THEN
         INFO = -23
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'FB01RD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         IF ( LJOBK ) THEN
            DWORK(1) = TWO
            DWORK(2) = ONE
         ELSE
            DWORK(1) = ONE
         END IF
         RETURN
      END IF
C
C     Construction of the needed part of the pre-array in DWORK.
C     To save workspace, only the blocks (1,3), (2,2), and (2,3) will be
C     constructed as shown below.
C
C     Storing C x S and A x S in the (1,1) and (2,1) blocks of DWORK,
C     respectively. The lower trapezoidal structure of [ C' A' ]' is
C     fully exploited. Specifically, if P <= N, the following partition
C     is used:
C
C        [ C1  0  ] [ S1  0  ]
C        [ A1  A3 ] [ S2  S3 ],
C        [ A2  A4 ]
C
C     where C1, S1, and A2 are P-by-P matrices, A1 and S2 are
C     (N-P)-by-P, A3 and S3 are (N-P)-by-(N-P), A4 is P-by-(N-P), and
C     C1, S1, A3, and S3 are lower triangular. The left hand side
C     matrix above is stored in the workspace. If P > N, the partition
C     is:
C
C        [ C1 ]
C        [ C2 ] [ S ],
C        [ A  ]
C
C     where C1 and C2 are N-by-N and (P-N)-by-N matrices, respectively,
C     and C1 and S are lower triangular.
C
C     Workspace: need (P+N)*N.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      CALL DLACPY( 'Lower', P, MIN( N, P ), C, LDC, DWORK, PN )
      CALL DLACPY( 'Full',  N, MIN( N, P ), A, LDA, DWORK(P+1), PN )
      IF ( N.GT.P )
     $   CALL DLACPY( 'Lower', N, N-P, A(1,P+1), LDA, DWORK(P*PN+P+1),
     $                PN )
C
C             [ C1  0  ]
C     Compute [        ] x S or C1 x S as a product of lower triangular
C             [ A1  A3 ]
C     matrices.
C     Workspace: need (P+N+1)*N.
C
      II = 1
      PL = N*PN + 1
      WRKOPT = PL + N - 1
C
      DO 10 I = 1, N
         CALL DCOPY( N-I+1, S(I,I), 1, DWORK(PL), 1 )
         CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', N-I+1,
     $               DWORK(II), PN, DWORK(PL), 1 )
         CALL DCOPY( N-I+1, DWORK(PL), 1, DWORK(II), 1 )
         II = II + PN + 1
   10 CONTINUE
C
C     Compute [ A2  A4 ] x S.
C
      CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Non-unit', P, N,
     $            ONE, S, LDS, DWORK(N+1), PN )
C
C     Triangularization (2 steps).
C
C     Step 1: annihilate the matrix C x S (hence C1 x S1, if P <= N).
C     Workspace: need (N+P)*N + 2*P.
C
      ITAU  = PL
      JWORK = ITAU + P
C
      CALL MB04LD( 'Lower', P, N, N, R, LDR, DWORK, PN, DWORK(P+1), PN,
     $             K, LDK, DWORK(ITAU), DWORK(JWORK) )
      WRKOPT = MAX( WRKOPT, PN*N + 2*P )
C
C     Now, the workspace for C x S is no longer needed.
C     Adjust the leading dimension of DWORK, to save space for the
C     following computations, and make room for B x Q.
C
      CALL DLACPY( 'Full', N, N, DWORK(P+1), PN, DWORK, N )
C
      DO 20 I = N*( N - 1 ) + 1, 1, -N
         CALL DCOPY( N, DWORK(I), 1, DWORK(I+N*M), 1 )
   20 CONTINUE
C
C     Storing B x Q in the (1,1) block of DWORK.
C     Workspace: need N*(M+N).
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
      IF ( .NOT.LMULTB )
     $   CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Non-unit', N, M,
     $               ONE, Q, LDQ, DWORK, N )
C
C     Step 2: LQ triangularization of the matrix [ B x Q  A x S ], where
C     A x S was modified at Step 1.
C     Workspace: need   N*(N+M+2);
C                prefer N*(N+M+1)+(P+1)*NB, where NB is the optimal
C                       block size for DGELQF (called in MB04JD).
C
      ITAU  = N*( M + N ) + 1
      JWORK = ITAU + N
C
      CALL MB04JD( N, M+N, MAX( N-P-1, 0 ), 0, DWORK, N, DWORK, N,
     $             DWORK(ITAU), DWORK(JWORK), LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     Output S and K (if needed) and set the optimal workspace
C     dimension (and the reciprocal of the condition number estimate).
C
      CALL DLACPY( 'Lower', N, N, DWORK, N, S, LDS )
C
      IF ( LJOBK ) THEN
C
C        Compute K.
C        Workspace: need 3*P.
C
         CALL MB02OD( 'Right', 'Lower', 'No transpose', 'Non-unit',
     $                '1-norm', N, P, ONE, R, LDR, K, LDK, RCOND, TOL,
     $                IWORK, DWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            WRKOPT = MAX( WRKOPT, 3*P )
            DWORK(2) = RCOND
         END IF
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of FB01RD ***
      END
