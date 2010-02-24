      SUBROUTINE FB01TD( JOBX, MULTRC, N, M, P, SINV, LDSINV, AINV,
     $                   LDAINV, AINVB, LDAINB, RINV, LDRINV, C, LDC,
     $                   QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,
     $                   DWORK, LDWORK, INFO )
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
C     given for the square root information filter, using the condensed
C     controller Hessenberg form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBX    CHARACTER*1
C             Indicates whether X    is to be computed as follows:
C                                i+1
C             = 'X':  X    is computed and stored in array X;
C                      i+1
C             = 'N':  X    is not required.
C                      i+1
C
C     MULTRC  CHARACTER*1             -1/2
C             Indicates how matrices R     and C    are to be passed to
C                                     i+1       i+1
C             the routine as follows:
C             = 'P':  Array RINV is not used and the array C must
C                                          -1/2
C                     contain the product R    C   ;
C                                          i+1  i+1
C             = 'N':  Arrays RINV and C must contain the matrices
C                     as described below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C                       -1      -1
C             matrices S   and A  .  N >= 0.
C                       i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C              -1/2
C             Q    .  M >= 0.
C              i
C
C     P       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C              -1/2
C             R    .  P >= 0.
C              i+1
C
C     SINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDSINV,N)
C             On entry, the leading N-by-N upper triangular part of this
C                                 -1
C             array must contain S  , the inverse of the square root
C                                 i
C             (right Cholesky factor) of the state covariance matrix
C             P    (hence the information square root) at instant i.
C              i|i
C             On exit, the leading N-by-N upper triangular part of this
C                             -1
C             array contains S   , the inverse of the square root (right
C                             i+1
C             Cholesky factor) of the state covariance matrix P
C                                                              i+1|i+1
C             (hence the information square root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDSINV  INTEGER
C             The leading dimension of array SINV.  LDSINV >= MAX(1,N).
C
C     AINV    (input) DOUBLE PRECISION array, dimension (LDAINV,N)
C                                                                 -1
C             The leading N-by-N part of this array must contain A  ,
C             the inverse of the state transition matrix of the discrete
C             system in controller Hessenberg form (e.g., as produced by
C             SLICOT Library Routine TB01MD).
C
C     LDAINV  INTEGER
C             The leading dimension of array AINV.  LDAINV >= MAX(1,N).
C
C     AINVB   (input) DOUBLE PRECISION array, dimension (LDAINB,M)
C                                                                  -1
C             The leading N-by-M part of this array must contain  A  B,
C                             -1
C             the product of A   and the input weight matrix B of the
C             discrete system, in upper controller Hessenberg form
C             (e.g., as produced by SLICOT Library Routine TB01MD).
C
C     LDAINB  INTEGER
C             The leading dimension of array AINVB.  LDAINB >= MAX(1,N).
C
C     RINV    (input) DOUBLE PRECISION array, dimension (LDRINV,*)
C             If MULTRC = 'N', then the leading P-by-P upper triangular
C                                              -1/2
C             part of this array must contain R    , the inverse of the
C                                              i+1
C             covariance square root (right Cholesky factor) of the
C             output (measurement) noise (hence the information square
C             root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C             Otherwise, RINV is not referenced and can be supplied as a
C             dummy array (i.e., set parameter LDRINV = 1 and declare
C             this array to be RINV(1,1) in the calling program).
C
C     LDRINV  INTEGER
C             The leading dimension of array RINV.
C             LDRINV >= MAX(1,P) if MULTRC = 'N';
C             LDRINV >= 1        if MULTRC = 'P'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain C   ,
C                                                       -1/2      i+1
C             the output weight matrix (or the product R    C    if
C                                                       i+1  i+1
C             MULTRC = 'P') of the discrete system at instant i+1.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     QINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDQINV,M)
C             On entry, the leading M-by-M upper triangular part of this
C                                 -1/2
C             array must contain Q    , the inverse of the covariance
C                                 i
C             square root (right Cholesky factor) of the input (process)
C             noise (hence the information square root) at instant i.
C             On exit, the leading M-by-M upper triangular part of this
C                                    -1/2
C             array contains (QINOV )    , the inverse of the covariance
C                                  i
C             square root (right Cholesky factor) of the process noise
C             innovation (hence the information square root) at
C             instant i.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDQINV  INTEGER
C             The leading dimension of array QINV.  LDQINV >= MAX(1,M).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain X , the estimated
C                                                i
C             filtered state at instant i.
C             On exit, if JOBX = 'X', and INFO = 0, then this array
C             contains X   , the estimated filtered state at
C                       i+1
C             instant i+1.
C             On exit, if JOBX = 'N', or JOBX = 'X' and INFO = 1, then
C                                  -1
C             this array contains S   X   .
C                                  i+1 i+1
C
C     RINVY   (input) DOUBLE PRECISION array, dimension (P)
C                                      -1/2
C             This array must contain R    Y   , the product of the
C                                      i+1  i+1
C                                      -1/2
C             upper triangular matrix R     and the measured output
C                                      i+1
C             vector Y    at instant i+1.
C                     i+1
C
C     Z       (input) DOUBLE PRECISION array, dimension (M)
C             This array must contain Z , the mean value of the state
C                                      i
C             process noise at instant i.
C
C     E       (output) DOUBLE PRECISION array, dimension (P)
C             This array contains E   , the estimated error at instant
C                                  i+1
C             i+1.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If JOBX = 'X', then TOL is used to test for near
C                                        -1
C             singularity of the matrix S   . If the user sets
C                                        i+1
C             TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then an implicitly computed, default
C             tolerance, defined by TOLDEF = N*N*EPS, is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C             Otherwise, TOL is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             where LIWORK = N if JOBX = 'X',
C             and   LIWORK = 1 otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.  If INFO = 0 and JOBX = 'X', DWORK(2) returns
C             an estimate of the reciprocal of the condition number
C                                 -1
C             (in the 1-norm) of S   .
C                                 i+1
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1)),
C                                 if JOBX = 'N';
C             LDWORK >= MAX(2,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1),
C                           3*N), if JOBX = 'X'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;                        -1
C             = 1:  if JOBX = 'X' and the matrix S    is singular,
C                                                 i+1       -1
C                   i.e., the condition number estimate of S    (in the
C                                                           i+1
C                                                         -1    -1/2
C                   1-norm) exceeds 1/TOL.  The matrices S   , Q
C                                                         i+1   i
C                   and E have been computed.
C
C     METHOD
C
C     The routine performs one recursion of the square root information
C     filter algorithm, summarized as follows:
C
C       |    -1/2             -1/2    |     |         -1/2             |
C       |   Q         0      Q    Z   |     | (QINOV )     *     *     |
C       |    i                i    i  |     |       i                  |
C       |                             |     |                          |
C       |           -1/2      -1/2    |     |             -1    -1     |
C     T |    0     R    C    R    Y   |  =  |    0       S     S   X   |
C       |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1|
C       |                             |     |                          |
C       |  -1 -1     -1 -1    -1      |     |                          |
C       | S  A  B   S  A     S  X     |     |    0         0     E     |
C       |  i         i        i  i    |     |                     i+1  |
C
C                   (Pre-array)                      (Post-array)
C
C     where T is an orthogonal transformation triangularizing the
C                        -1/2
C     pre-array, (QINOV )     is the inverse of the covariance square
C                      i
C     root (right Cholesky factor) of the process noise innovation
C                                                            -1  -1
C     (hence the information square root) at instant i and (A  ,A  B) is
C     in upper controller Hessenberg form.
C
C     An example of the pre-array is given below (where N = 6, M = 2,
C     and P = 3):
C
C         |x x |             | x|
C         |  x |             | x|
C         _______________________
C         |    | x x x x x x | x|
C         |    | x x x x x x | x|
C         |    | x x x x x x | x|
C         _______________________
C         |x x | x x x x x x | x|
C         |  x | x x x x x x | x|
C         |    | x x x x x x | x|
C         |    |   x x x x x | x|
C         |    |     x x x x | x|
C         |    |       x x x | x|
C
C     The inverse of the corresponding state covariance matrix P
C                                                               i+1|i+1
C     (hence the information matrix I) is then factorized as
C
C                    -1         -1     -1
C         I       = P       = (S   )' S
C          i+1|i+1   i+1|i+1    i+1    i+1
C
C     and one combined time and measurement update for the state is
C     given by X   .
C               i+1
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
C     The algorithm requires approximately
C
C           3    2                           2          3
C     (1/6)N  + N x (3/2 x M + P) + 2 x N x M  + 2/3 x M
C
C     operations and is backward stable (see [3]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01HD by M. Vanbegin,
C     P. Van Dooren, and M.H.G. Verhaegen.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, February 14, 2004.
C
C     KEYWORDS
C
C     Controller Hessenberg form, Kalman filtering, optimal filtering,
C     orthogonal transformation, recursive estimation, square-root
C     filtering, square-root information filtering.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBX, MULTRC
      INTEGER           INFO, LDAINB, LDAINV, LDC, LDQINV, LDRINV,
     $                  LDSINV, LDWORK, M, N, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  AINV(LDAINV,*), AINVB(LDAINB,*), C(LDC,*),
     $                  DWORK(*), E(*), QINV(LDQINV,*), RINV(LDRINV,*),
     $                  RINVY(*), SINV(LDSINV,*), X(*), Z(*)
C     .. Local Scalars ..
      LOGICAL           LJOBX, LMULTR
      INTEGER           I, I12, I13, I23, I32, I33, II, IJ, ITAU, JWORK,
     $                  LDW, M1, MP1, N1, NM, NP, WRKOPT
      DOUBLE PRECISION  RCOND
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DLACPY, DTRMM, DTRMV, MB02OD,
     $                  MB04ID, MB04KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      NP = N + P
      NM = N + M
      N1 = MAX( 1, N )
      M1 = MAX( 1, M )
      MP1 = M + 1
      INFO = 0
      LJOBX  = LSAME( JOBX, 'X' )
      LMULTR = LSAME( MULTRC, 'P' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBX .AND. .NOT.LSAME( JOBX, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LMULTR .AND. .NOT.LSAME( MULTRC, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDSINV.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDAINV.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDAINB.LT.N1 ) THEN
         INFO = -11
      ELSE IF( LDRINV.LT.1 .OR. ( .NOT.LMULTR .AND. LDRINV.LT.P ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDQINV.LT.M1 ) THEN
         INFO = -17
      ELSE IF( ( LJOBX .AND. LDWORK.LT.MAX( 2, N*(NM + M) + 3*M,
     $                                      NP*(N + 1) + N +
     $                                      MAX( N - 1, MP1 ), 3*N ) )
     $                                                             .OR.
     $    ( .NOT.LJOBX .AND. LDWORK.LT.MAX( 1, N*(NM + M) + 3*M,
     $                                      NP*(N + 1) + N +
     $                                      MAX( N - 1, MP1 ) ) ) ) THEN
         INFO = -25
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'FB01TD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, P ).EQ.0 ) THEN
         IF ( LJOBX ) THEN
            DWORK(1) = TWO
            DWORK(2) = ONE
         ELSE
            DWORK(1) = ONE
         END IF
         RETURN
      END IF
C
C     Construction of the needed part of the pre-array in DWORK.
C     To save workspace, only the blocks (1,3), (3,1)-(3,3), (2,2), and
C     (2,3) will be constructed when needed as shown below.
C
C     Storing SINV x AINVB and SINV x AINV in the (1,1) and (1,2)
C     blocks of DWORK, respectively. The upper trapezoidal structure of
C     [ AINVB AINV ] is fully exploited. Specifically, if M <= N, the
C     following partition is used:
C
C       [ S1  S2 ] [ B1  A1 A3 ]
C       [ 0   S3 ] [ 0   A2 A4 ],
C
C     where B1, A3, and S1 are M-by-M matrices, A1 and S2 are
C     M-by-(N-M), A2 and S3 are (N-M)-by-(N-M), A4 is (N-M)-by-M, and
C     B1, S1, A2, and S3 are upper triangular. The right hand side
C     matrix above is stored in the workspace. If M > N, the partition
C     is [ SINV ] [ B1 B2  A ], where B1 is N-by-N, B2 is N-by-(M-N),
C     and B1 and SINV are upper triangular.
C     The variables called Ixy define the starting positions where the
C     (x,y) blocks of the pre-array are initially stored in DWORK.
C     Workspace: need N*(M+N).
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      LDW = N1
      I32 = N*M + 1
C
      CALL DLACPY( 'Upper', N, M, AINVB, LDAINB, DWORK, LDW )
      CALL DLACPY( 'Full',  MIN( M, N ), N, AINV, LDAINV, DWORK(I32),
     $             LDW )
      IF ( N.GT.M )
     $   CALL DLACPY( 'Upper', N-M, N, AINV(MP1,1), LDAINV,
     $                DWORK(I32+M), LDW )
C
C                    [ B1  A1 ]
C     Compute SINV x [ 0   A2 ] or SINV x B1 as a product of upper
C     triangular matrices.
C     Workspace: need N*(M+N+1).
C
      II  = 1
      I13 = N*NM + 1
      WRKOPT = MAX( 1, N*NM + N )
C
      DO 10 I = 1, N
         CALL DCOPY( I, DWORK(II), 1, DWORK(I13), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I, SINV,
     $               LDSINV, DWORK(I13), 1 )
         CALL DCOPY( I, DWORK(I13), 1, DWORK(II), 1 )
         II = II + N
   10 CONTINUE
C
C                    [ A3 ]
C     Compute SINV x [ A4 ] or SINV x [ B2 A ].
C
      CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, M,
     $            ONE, SINV, LDSINV, DWORK(II), LDW )
C
C     Storing the process noise mean value in (1,3) block of DWORK.
C     Workspace: need N*(M+N) + M.
C
      CALL DCOPY( M, Z, 1, DWORK(I13), 1 )
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', M, QINV, LDQINV,
     $            DWORK(I13), 1 )
C
C     Computing SINV x X in X.
C
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', N, SINV, LDSINV,
     $            X, 1 )
C
C     Triangularization (2 steps).
C
C     Step 1: annihilate the matrix SINV x AINVB.
C     Workspace: need N*(N+2*M) + 3*M.
C
      I12   = I13  + M
      ITAU  = I12  + M*N
      JWORK = ITAU + M
C
      CALL MB04KD( 'Upper', M, N, N, QINV, LDQINV, DWORK, LDW,
     $             DWORK(I32), LDW, DWORK(I12), M1, DWORK(ITAU),
     $             DWORK(JWORK) )
      WRKOPT = MAX( WRKOPT, N*( NM + M ) + 3*M )
C
      IF ( N.EQ.0 ) THEN
         CALL DCOPY( P, RINVY, 1, E, 1 )
         IF ( LJOBX )
     $      DWORK(2) = ONE
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Apply the transformations to the last column of the pre-array.
C     (Only the updated (3,3) block is now needed.)
C
      IJ = 1
C
      DO 20 I = 1, M
         CALL DAXPY( MIN( I, N ), -DWORK(ITAU+I-1)*( DWORK(I13+I-1) +
     $               DDOT( MIN( I, N ), DWORK(IJ), 1, X, 1 ) ),
     $                                  DWORK(IJ), 1, X, 1 )
         IJ = IJ + N
   20 CONTINUE
C
C     Now, the workspace for SINV x AINVB, as well as for the updated
C     (1,2) block of the pre-array, are no longer needed.
C     Move the computed (3,2) and (3,3) blocks of the pre-array in the
C     (1,1) and (1,2) block positions of DWORK, to save space for the
C     following computations.
C     Then, adjust the implicitly defined leading dimension of DWORK,
C     to make space for storing the (2,2) and (2,3) blocks of the
C     pre-array.
C     Workspace: need (P+N)*(N+1).
C
      CALL DLACPY( 'Full', MIN( M, N ), N, DWORK(I32), LDW, DWORK, LDW )
      IF ( N.GT.M )
     $   CALL DLACPY( 'Upper', N-M, N, DWORK(I32+M), LDW, DWORK(MP1),
     $                LDW )
      LDW = MAX( 1, NP )
C
      DO 40 I = N, 1, -1
         DO 30 IJ = MIN( N, I+M ), 1, -1
            DWORK(NP*(I-1)+P+IJ) = DWORK(N*(I-1)+IJ)
   30    CONTINUE
   40 CONTINUE
C
C     Copy of RINV x C in the (1,1) block of DWORK.
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK, LDW )
      IF ( .NOT.LMULTR )
     $   CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', P, N,
     $               ONE, RINV, LDRINV, DWORK, LDW )
C
C     Copy the inclusion measurement in the (1,2) block and the updated
C     X in the (2,2) block of DWORK.
C
      I23 = NP*N + 1
      I33 = I23  + P
      CALL DCOPY( P, RINVY, 1, DWORK(I23), 1 )
      CALL DCOPY( N, X, 1, DWORK(I33), 1 )
      WRKOPT = MAX( WRKOPT, NP*( N + 1 ) )
C
C     Step 2: QR factorization of the first block column of the matrix
C
C        [ RINV x C     RINV x Y ],
C        [ SINV x AINV  SINV x X ]
C
C     where the second block row was modified at Step 1.
C     Workspace: need   (P+N)*(N+1) + N + MAX(N-1,M+1);
C                prefer (P+N)*(N+1) + N + (M+1)*NB, where NB is the
C                       optimal block size for DGEQRF called in MB04ID.
C
      ITAU  = I23  + NP
      JWORK = ITAU + N
C
      CALL MB04ID( NP, N, MAX( N-MP1, 0 ), 1, DWORK, LDW, DWORK(I23),
     $             LDW, DWORK(ITAU), DWORK(JWORK), LDWORK-JWORK+1,
     $             INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     Output SINV, X, and E and set the optimal workspace dimension
C     (and the reciprocal of the condition number estimate).
C
      CALL DLACPY( 'Upper', N, N, DWORK, LDW, SINV, LDSINV )
      CALL DCOPY( N, DWORK(I23), 1, X, 1 )
      IF( P.GT.0 )
     $   CALL DCOPY( P, DWORK(I23+N), 1, E, 1 )
C
      IF ( LJOBX ) THEN
C
C        Compute X.
C        Workspace: need 3*N.
C
         CALL MB02OD( 'Left', 'Upper', 'No transpose', 'Non-unit',
     $                '1-norm', N, 1, ONE, SINV, LDSINV, X, N, RCOND,
     $                TOL, IWORK, DWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            WRKOPT = MAX( WRKOPT, 3*N )
            DWORK(2) = RCOND
         END IF
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of FB01TD***
      END
