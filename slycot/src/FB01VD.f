      SUBROUTINE FB01VD( N, M, L, P, LDP, A, LDA, B, LDB, C, LDC, Q,
     $                   LDQ, R, LDR, K, LDK, TOL, IWORK, DWORK, LDWORK,
     $                   INFO )
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
C     To compute one recursion of the conventional Kalman filter
C     equations. This is one update of the Riccati difference equation
C     and the Kalman filter gain.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices P      and A .  N >= 0.
C                       i|i-1      i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C             Q .  M >= 0.
C              i
C
C     L       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C             R .  L >= 0.
C              i
C
C     P       (input/output) DOUBLE PRECISION array, dimension (LDP,N)
C             On entry, the leading N-by-N part of this array must
C             contain P     , the state covariance matrix at instant
C                      i|i-1
C             (i-1). The upper triangular part only is needed.
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains P     , the state covariance matrix at
C                             i+1|i
C             instant i. The strictly lower triangular part is not set.
C             Otherwise, the leading N-by-N part of this array contains
C             P     , its input value.
C              i|i-1
C
C     LDP     INTEGER
C             The leading dimension of array P.  LDP >= MAX(1,N).
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain A ,
C                                                                 i
C             the state transition matrix of the discrete system at
C             instant i.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain B ,
C                                                                 i
C             the input weight matrix of the discrete system at
C             instant i.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading L-by-N part of this array must contain C ,
C                                                                 i
C             the output weight matrix of the discrete system at
C             instant i.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,L).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,M)
C             The leading M-by-M part of this array must contain Q ,
C                                                                 i
C             the input (process) noise covariance matrix at instant i.
C             The diagonal elements of this array are modified by the
C             routine, but are restored on exit.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,M).
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,L)
C             On entry, the leading L-by-L part of this array must
C             contain R , the output (measurement) noise covariance
C                      i
C             matrix at instant i.
C             On exit, if INFO = 0, or INFO = L+1, the leading L-by-L
C                                                                  1/2
C             upper triangular part of this array contains (RINOV )   ,
C                                                                i
C             the square root (left Cholesky factor) of the covariance
C             matrix of the innovations at instant i.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,L).
C
C     K       (output) DOUBLE PRECISION array, dimension (LDK,L)
C             If INFO = 0, the leading N-by-L part of this array
C             contains K , the Kalman filter gain matrix at instant i.
C                       i
C             If INFO > 0, the leading N-by-L part of this array
C             contains the matrix product P     C'.
C                                          i|i-1 i
C
C     LDK     INTEGER
C             The leading dimension of array K.  LDK >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the matrix RINOV . If the user sets TOL > 0, then the
C                             i
C             given value of TOL is used as a lower bound for the
C             reciprocal condition number of that matrix; a matrix whose
C             estimated condition number is less than 1/TOL is
C             considered to be nonsingular. If the user sets TOL <= 0,
C             then an implicitly computed, default tolerance, defined by
C             TOLDEF = L*L*EPS, is used instead, where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (L)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, or INFO = L+1, DWORK(1) returns an
C             estimate of the reciprocal of the condition number (in the
C             1-norm) of the matrix RINOV .
C                                        i
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,L*N+3*L,N*N,N*M).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -k, the k-th argument had an illegal
C                   value;
C             = k:  if INFO = k, 1 <= k <= L, the leading minor of order
C                   k of the matrix RINOV  is not positive-definite, and
C                                        i
C                   its Cholesky factorization could not be completed;
C             = L+1: the matrix RINOV  is singular, i.e., the condition
C                                    i
C                   number estimate of RINOV  (in the 1-norm) exceeds
C                                           i
C                   1/TOL.
C
C     METHOD
C
C     The conventional Kalman filter gain used at the i-th recursion
C     step is of the form
C
C                            -1
C        K  = P     C'  RINOV  ,
C         i    i|i-1 i       i
C
C     where RINOV  = C P     C' + R , and the state covariance matrix
C                i    i i|i-1 i    i
C
C     P      is updated by the discrete-time difference Riccati equation
C      i|i-1
C
C        P      = A  (P      - K C P     ) A'  + B Q B'.
C         i+1|i    i   i|i-1    i i i|i-1   i     i i i
C
C     Using these two updates, the combined time and measurement update
C     of the state X      is given by
C                   i|i-1
C
C        X      = A X      + A K (Y  - C X     ),
C         i+1|i    i i|i-1    i i  i    i i|i-1
C
C     where Y  is the new observation at step i.
C            i
C
C     REFERENCES
C
C     [1] Anderson, B.D.O. and Moore, J.B.
C         Optimal Filtering,
C         Prentice Hall, Englewood Cliffs, New Jersey, 1979.
C
C     [2] Verhaegen, M.H.G. and Van Dooren, P.
C         Numerical Aspects of Different Kalman Filter Implementations.
C         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, 1986.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C
C             3   2
C      3/2 x N + N  x (3 x L + M/2)
C
C     operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01JD by M.H.G. Verhaegen,
C     M. Vanbegin, and P. Van Dooren.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, April 20, 2004.
C
C     KEYWORDS
C
C     Kalman filtering, optimal filtering, recursive estimation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDC, LDK, LDP, LDQ, LDR,
     $                  LDWORK, M, N
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*),
     $                  K(LDK,*), P(LDP,*), Q(LDQ,*), R(LDR,*)
C     .. Local Scalars ..
      INTEGER           J, JWORK, LDW, N1
      DOUBLE PRECISION  RCOND, RNORM, TOLDEF
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANSY
      EXTERNAL          DLAMCH, DLANSY
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DLACPY, DLASET, DPOCON,
     $                  DPOTRF, DSCAL, DTRMM, DTRSM, MB01RD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      N1 = MAX( 1, N )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDP.LT.N1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, L ) ) THEN
         INFO = -11
      ELSE IF( LDQ.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDR.LT.MAX( 1, L ) ) THEN
         INFO = -15
      ELSE IF( LDK.LT.N1 ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.MAX( 1, L*N + 3*L, N*N, N*M ) ) THEN
         INFO = -21
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'FB01VD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, L ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Efficiently compute RINOV = CPC' + R in R and put CP in DWORK and
C     PC' in K. (The content of DWORK on exit from MB01RD is used.)
C     Workspace: need L*N.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code.)
C
      CALL MB01RD( 'Upper', 'No transpose', L, N, ONE, ONE, R, LDR, C,
     $             LDC, P, LDP, DWORK, LDWORK, INFO )
      LDW = MAX( 1, L )
C
      DO 10 J = 1, L
         CALL DCOPY( N, DWORK(J), LDW, K(1,J), 1 )
   10 CONTINUE
C
      CALL DLACPY( 'Full', L, N, C, LDC, DWORK, LDW )
      CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', L, N, ONE,
     $            P, LDP, DWORK, LDW )
      CALL DSCAL( N, TWO, P, LDP+1 )
C
      DO 20 J = 1, L
         CALL DAXPY( N, ONE, K(1,J), 1, DWORK(J), LDW )
         CALL DCOPY( N, DWORK(J), LDW, K(1,J), 1 )
   20 CONTINUE
C
C     Calculate the Cholesky decomposition U'U of the innovation
C     covariance matrix RINOV, and its reciprocal condition number.
C     Workspace: need L*N + 3*L.
C
      JWORK = L*N + 1
      RNORM = DLANSY( '1-norm', 'Upper', L, R, LDR, DWORK(JWORK) )
C
      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO )
     $   TOLDEF = DBLE( L*L )*DLAMCH( 'Epsilon' )
      CALL DPOTRF( 'Upper', L, R, LDR, INFO )
      IF ( INFO.NE.0 )
     $   RETURN
C
      CALL DPOCON( 'Upper', L, R, LDR, RNORM, RCOND, DWORK(JWORK),
     $             IWORK, INFO )
C
      IF ( RCOND.LT.TOLDEF ) THEN
C
C        Error return: RINOV is numerically singular.
C
         INFO = L+1
         DWORK(1) = RCOND
         RETURN
      END IF
C
      IF ( L.GT.1 )
     $   CALL DLASET( 'Lower', L-1, L-1, ZERO, ZERO, R(2,1),LDR )
C                                                          -1
C     Calculate the Kalman filter gain matrix  K = PC'RINOV .
C     Workspace: need L*N.
C
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', N, L,
     $            ONE, R, LDR, K, LDK )
      CALL DTRSM( 'Right', 'Upper', 'Transpose', 'Non-unit', N, L,
     $            ONE, R, LDR, K, LDK )
C
C     First part of the Riccati equation update: compute A(P-KCP)A'.
C     The upper triangular part of the symmetric matrix P-KCP is formed.
C     Workspace: need max(L*N,N*N).
C
      JWORK = 1
C
      DO 30 J = 1, N
         CALL DGEMV( 'No transpose', J, L, -ONE, K, LDK, DWORK(JWORK),
     $               1, ONE, P(1,J), 1 )
         JWORK = JWORK + L
   30 CONTINUE
C
      CALL MB01RD( 'Upper', 'No transpose', N, N, ZERO, ONE, P, LDP, A,
     $             LDA, P, LDP, DWORK, LDWORK, INFO )
C
C     Second part of the Riccati equation update: add BQB'.
C     Workspace: need N*M.
C
      CALL MB01RD( 'Upper', 'No transpose', N, M, ONE, ONE, P, LDP, B,
     $             LDB, Q, LDQ, DWORK, LDWORK, INFO )
      CALL DSCAL( M, TWO, Q, LDQ+1 )
C
C     Set the reciprocal of the condition number estimate.
C
      DWORK(1) = RCOND
C
      RETURN
C *** Last line of FB01VD ***
      END
