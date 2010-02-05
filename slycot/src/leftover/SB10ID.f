      SUBROUTINE SB10ID( N, M, NP, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   FACTOR, NK, AK, LDAK, BK, LDBK, CK, LDCK,
     $                   DK, LDDK, RCOND, IWORK, DWORK, LDWORK, BWORK,
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
C     To compute the matrices of the positive feedback controller
C
C              | Ak | Bk |
C          K = |----|----|
C              | Ck | Dk |
C
C     for the shaped plant
C
C              | A | B |
C          G = |---|---|
C              | C | D |
C
C     in the McFarlane/Glover Loop Shaping Design Procedure.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the plant.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     NP      (input) INTEGER
C             The row size of the matrix C.  NP >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A of the shaped plant.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B of the shaped plant.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading NP-by-N part of this array must contain the
C             system output matrix C of the shaped plant.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system matrix D of the shaped plant.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     FACTOR  (input) DOUBLE PRECISION
C             = 1 implies that an optimal controller is required;
C             > 1 implies that a suboptimal controller is required,
C                 achieving a performance FACTOR less than optimal.
C             FACTOR >= 1.
C
C     NK      (output) INTEGER
C             The order of the positive feedback controller.  NK <= N.
C
C     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N)
C             The leading NK-by-NK part of this array contains the
C             controller state matrix Ak.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NP)
C             The leading NK-by-NP part of this array contains the
C             controller input matrix Bk.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N)
C             The leading M-by-NK part of this array contains the
C             controller output matrix Ck.
C
C     LDCK    INTEGER
C             The leading dimension of the array CK.  LDCK >= max(1,M).
C
C     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NP)
C             The leading M-by-NP part of this array contains the
C             controller matrix Dk.
C
C     LDDK    INTEGER
C             The leading dimension of the array DK.  LDDK >= max(1,M).
C
C     RCOND   (output) DOUBLE PRECISION array, dimension (2)
C             RCOND(1) contains an estimate of the reciprocal condition
C                      number of the X-Riccati equation;
C             RCOND(2) contains an estimate of the reciprocal condition
C                      number of the Z-Riccati equation.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(2*N,N*N,M,NP)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 4*N*N + M*M + NP*NP + 2*M*N + N*NP + 4*N +
C                       max( 6*N*N + 5 + max(1,4*N*N+8*N), N*NP + 2*N ).
C             For good performance, LDWORK must generally be larger.
C             An upper bound of LDWORK in the above formula is
C             LDWORK >= 10*N*N + M*M + NP*NP + 2*M*N + 2*N*NP + 4*N +
C                       5 + max(1,4*N*N+8*N).
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the X-Riccati equation is not solved successfully;
C             = 2:  the Z-Riccati equation is not solved successfully;
C             = 3:  the iteration to compute eigenvalues or singular
C                   values failed to converge;
C             = 4:  the matrix Ip - D*Dk is singular;
C             = 5:  the matrix Im - Dk*D is singular;
C             = 6:  the closed-loop system is unstable.
C
C     METHOD
C
C     The routine implements the formulas given in [1].
C
C     REFERENCES
C
C     [1] McFarlane, D. and Glover, K.
C         A loop shaping design procedure using H_infinity synthesis.
C         IEEE Trans. Automat. Control, vol. AC-37, no. 6, pp. 759-769,
C         1992.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the results depends on the conditioning of the
C     two Riccati equations solved in the controller design (see the
C     output parameter RCOND).
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000,
C     Feb. 2001.
C
C     KEYWORDS
C
C     H_infinity control, Loop-shaping design, Robust control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                   LDDK, LDWORK, M, N, NK, NP
      DOUBLE PRECISION   FACTOR
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), DWORK( * ),
     $                   RCOND( 2 )
C     ..
C     .. Local Scalars ..
      CHARACTER*1        HINV
      INTEGER            I, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     $                   I11, I12, I13, INFO2, IWRK, J, LWA, LWAMAX,
     $                   MINWRK, N2, NS, SDIM
      DOUBLE PRECISION   SEP, FERR, GAMMA
C     ..
C     .. External Functions ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEES, DGEMM, DLACPY, DLASET, DPOTRF, DPOTRS,
     $                   DSYRK, DTRSM, MB02VD, SB02RD, SB10JD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( NP.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, NP ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -11
      ELSE IF( FACTOR.LT.ONE ) THEN
         INFO = -12
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDCK.LT.MAX( 1, M ) ) THEN
         INFO = -19
      ELSE IF( LDDK.LT.MAX( 1, M ) ) THEN
         INFO = -21
      END IF
C
C     Compute workspace.
C
      MINWRK = 4*N*N + M*M + NP*NP + 2*M*N + N*NP + 4*N +
     $         MAX( 6*N*N + 5 + MAX( 1, 4*N*N + 8*N ), N*NP + 2*N )
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -25
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10ID', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 ) THEN
         RCOND( 1 ) = ONE
         RCOND( 2 ) = ONE
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Workspace usage.
C
      I1  = N*N
      I2  = I1  + N*N
      I3  = I2  + M*N
      I4  = I3  + M*N
      I5  = I4  + M*M
      I6  = I5  + NP*NP
      I7  = I6  + NP*N
      I8  = I7  + N*N
      I9  = I8  + N*N
      I10 = I9  + N*N
      I11 = I10 + N*N
      I12 = I11 + 2*N
      I13 = I12 + 2*N
C
      IWRK = I13 + 4*N*N
C
C     Compute D'*C .
C
      CALL DGEMM( 'T', 'N', M, N, NP, ONE, D, LDD, C, LDC, ZERO,
     $            DWORK( I2+1 ), M )
C
C     Compute S = Im + D'*D .
C
      CALL DLASET( 'U', M, M, ZERO, ONE, DWORK( I4+1 ), M )
      CALL DSYRK( 'U', 'T', M, NP, ONE, D, LDD, ONE, DWORK( I4+1 ), M )
C
C     Factorize S, S = T'*T, with T upper triangular.
C
      CALL DPOTRF( 'U', M, DWORK( I4+1 ), M, INFO2 )
C
C              -1
C     Compute S  D'*C .
C
      CALL DPOTRS( 'U', M, N, DWORK( I4+1 ), M, DWORK( I2+1 ), M,
     $             INFO2 )
C
C                -1
C     Compute B*T  .
C
      CALL DLACPY( 'F', N, M, B, LDB, DWORK( I3+1 ), N )
      CALL DTRSM(  'R', 'U', 'N', 'N', N, M, ONE, DWORK( I4+1 ), M,
     $             DWORK( I3+1 ), N )
C
C     Compute R = Ip + D*D' .
C
      CALL DLASET( 'U', NP, NP, ZERO, ONE, DWORK( I5+1 ), NP )
      CALL DSYRK( 'U', 'N', NP, M, ONE, D, LDD, ONE, DWORK( I5+1 ), NP )
C
C     Factorize R, R = U'*U, with U upper triangular.
C
      CALL DPOTRF( 'U', NP, DWORK( I5+1 ), NP, INFO2 )
C
C              -T
C     Compute U  C .
C
      CALL DLACPY( 'F', NP, N, C, LDC, DWORK( I6+1 ), NP )
      CALL DTRSM(  'L', 'U', 'T', 'N', NP, N, ONE, DWORK( I5+1 ), NP,
     $             DWORK( I6+1 ), NP )
C
C                         -1
C     Compute Ar = A - B*S  D'*C .
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I7+1 ), N )
      CALL DGEMM( 'N', 'N', N, N, M, -ONE, B, LDB, DWORK( I2+1 ), M,
     $            ONE, DWORK( I7+1 ), N )
C
C                                            -1
C     Compute the upper triangle of Cr = C'*R  *C .
C
      CALL DSYRK( 'U', 'T', N, NP, ONE, DWORK( I6+1 ), NP, ZERO,
     $            DWORK( I8+1 ), N )
C
C                                           -1
C     Compute the upper triangle of Dr = B*S  B' .
C
      CALL DSYRK( 'U', 'N', N, M, ONE, DWORK( I3+1 ), N, ZERO,
     $            DWORK( I9+1 ), N )
C
C     Solution of the Riccati equation Ar'*X + X*Ar + Cr - X*Dr*X = 0 .
C     Workspace:    need   10*N*N + M*M + NP*NP + 2*M*N + N*NP + 4*N +
C                                   5 + max(1,4*N*N+8*N).
C                   prefer larger.
C                   AK is used as workspace.
C
      N2 = 2*N
      CALL SB02RD( 'A', 'C', HINV, 'N', 'U', 'G', 'S', 'N', 'O', N,
     $             DWORK( I7+1 ), N, DWORK( I10+1 ), N, AK, LDAK,
     $             DWORK( I9+1 ), N, DWORK( I8+1 ), N, DWORK, N, SEP,
     $             RCOND( 1 ), FERR, DWORK( I11+1 ), DWORK( I12+1 ),
     $             DWORK( I13+1 ), N2, IWORK, DWORK( IWRK+1 ),
     $             LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( MINWRK, LWA )
C
C     Solution of the Riccati equation Ar*Z + Z*Ar' + Dr - Z*Cr*Z = 0 .
C
      CALL SB02RD( 'A', 'C', HINV, 'T', 'U', 'G', 'S', 'N', 'O', N,
     $             DWORK( I7+1 ), N, DWORK( I10+1 ), N, AK, LDAK,
     $             DWORK( I8+1 ), N, DWORK( I9+1 ), N, DWORK( I1+1 ),
     $             N, SEP, RCOND( 2 ), FERR, DWORK( I11+1 ),
     $             DWORK( I12+1 ), DWORK( I13+1 ), N2, IWORK,
     $             DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 2
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C                      -1        -1
C     Compute F1 = -( S  D'*C + S  B'*X ) .
C
      CALL DTRSM(  'R', 'U', 'T', 'N', N, M, ONE, DWORK( I4+1 ), M,
     $             DWORK( I3+1 ), N )
      CALL DGEMM( 'T', 'N', M, N, N, -ONE, DWORK( I3+1 ), N, DWORK, N,
     $            -ONE, DWORK( I2+1 ), M )
C
C     Compute gamma .
C
      CALL DGEMM( 'N', 'N', N, N, N, ONE, DWORK, N, DWORK( I1+1 ), N,
     $            ZERO, DWORK( I7+1 ), N )
      CALL DGEES( 'N', 'N', SELECT, N, DWORK( I7+1 ), N, SDIM,
     $            DWORK( I11+1 ), DWORK( I12+1 ), DWORK( IWRK+1 ), N,
     $            DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
      GAMMA = ZERO
      DO 10 I = 1, N
         GAMMA = MAX( GAMMA, DWORK( I11+I ) )
   10 CONTINUE
      GAMMA = FACTOR*SQRT( ONE + GAMMA )
C
C     Workspace usage.
C     Workspace:    need   4*N*N + M*N + N*NP.
C
      I4 = I3 + N*N
      I5 = I4 + N*N
C
C     Compute Ac = A + B*F1 .
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I4+1 ), N )
      CALL DGEMM( 'N', 'N', N, N, M, ONE, B, LDB, DWORK( I2+1 ), M,
     $            ONE, DWORK( I4+1 ), N )
C
C     Compute W1' = (1-gamma^2)*In + Z*X .
C
      CALL DLASET( 'F', N, N, ZERO, ONE-GAMMA*GAMMA, DWORK( I3+1 ), N )
      CALL DGEMM( 'N', 'N', N, N, N, ONE, DWORK( I1+1 ), N, DWORK, N,
     $            ONE, DWORK( I3+1 ), N )
C
C     Compute Bcp = gamma^2*Z*C' .
C
      CALL DGEMM( 'N', 'T', N, NP, N, GAMMA*GAMMA, DWORK( I1+1 ), N, C,
     $            LDC, ZERO, BK, LDBK )
C
C     Compute C + D*F1 .
C
      CALL DLACPY( 'F', NP, N, C, LDC, DWORK( I5+1 ), NP )
      CALL DGEMM( 'N', 'N', NP, N, M, ONE, D, LDD, DWORK( I2+1 ), M,
     $            ONE, DWORK( I5+1 ), NP )
C
C     Compute Acp = W1'*Ac + gamma^2*Z*C'*(C+D*F1) .
C
      CALL DGEMM( 'N', 'N', N, N, N, ONE, DWORK( I3+1 ), N,
     $            DWORK( I4+1 ), N, ZERO, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, NP, ONE, BK, LDBK,
     $            DWORK( I5+1 ), NP, ONE, AK, LDAK )
C
C     Compute Ccp = B'*X .
C
      CALL DGEMM( 'T', 'N', M, N, N, ONE, B, LDB, DWORK, N, ZERO,
     $             CK, LDCK )
C
C     Set Dcp = -D' .
C
      DO 30 I = 1, M
         DO 20 J = 1, NP
            DK( I, J ) = -D( J, I )
   20    CONTINUE
   30 CONTINUE
C
      IWRK = I4
C
C     Reduce the generalized state-space description to a regular one.
C     Workspace:             need   3*N*N + M*N.
C     Additional workspace:  need   2*N*N + 2*N + N*MAX(5,N+M+NP).
C                            prefer larger.
C
      CALL SB10JD( N, NP, M, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $             DWORK( I3+1 ), N, NK, DWORK( IWRK+1 ), LDWORK-IWRK,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C     Workspace usage.
C     Workspace:    need   4*N*N + M*M + NP*NP + 2*M*N + 2*N*NP.
C                          (NK <= N.)
C
      I2 = NP*NP
      I3 = I2 + NK*NP
      I4 = I3 + M*M
      I5 = I4 + N*M
      I6 = I5 + NP*NK
      I7 = I6 + M*N
C
      IWRK = I7 + ( N + NK )*( N + NK )
C
C     Compute Ip - D*Dk .
C
      CALL DLASET( 'Full', NP, NP, ZERO, ONE, DWORK, NP )
      CALL DGEMM( 'N', 'N', NP, NP, M, -ONE, D, LDD, DK, LDDK, ONE,
     $             DWORK, NP )
C
C                         -1
C     Compute Bk*(Ip-D*Dk)  .
C
      CALL DLACPY( 'F', NK, NP, BK, LDBK, DWORK( I2+1 ), NK )
      CALL MB02VD( 'N', NK, NP, DWORK, NP, IWORK, DWORK( I2+1 ), NK,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 4
         RETURN
      END IF
C
C     Compute Im - Dk*D .
C
      CALL DLASET( 'Full', M, M, ZERO, ONE, DWORK( I3+1 ), M )
      CALL DGEMM( 'N', 'N', M, M, NP, -ONE, DK, LDDK, D, LDD, ONE,
     $             DWORK( I3+1 ), M )
C
C                        -1
C     Compute B*(Im-Dk*D)  .
C
      CALL DLACPY( 'F', N, M, B, LDB, DWORK( I4+1 ), N )
      CALL MB02VD( 'N', N, M, DWORK( I3+1 ), M, IWORK, DWORK( I4+1 ), N,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 5
         RETURN
      END IF
C
C     Compute D*Ck .
C
      CALL DGEMM( 'N', 'N', NP, NK, M, ONE, D, LDD, CK, LDCK, ZERO,
     $             DWORK( I5+1 ), NP )
C
C     Compute Dk*C .
C
      CALL DGEMM( 'N', 'N', M, N, NP, ONE, DK, LDDK, C, LDC, ZERO,
     $            DWORK( I6+1 ), M )
C
C     Compute the closed-loop state matrix.
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I7+1 ), N+NK )
      CALL DGEMM( 'N', 'N', N, N, M, ONE, DWORK( I4+1 ), N,
     $            DWORK( I6+1 ), M, ONE, DWORK( I7+1 ), N+NK )
      CALL DGEMM( 'N', 'N', NK, N, NP, ONE, DWORK( I2+1 ), NK, C, LDC,
     $            ZERO, DWORK( I7+N+1 ), N+NK )
      CALL DGEMM( 'N', 'N', N, NK, M, ONE, DWORK( I4+1 ), N, CK, LDCK,
     $            ZERO, DWORK( I7+(N+NK)*N+1 ), N+NK )
      CALL DLACPY( 'F', NK, NK, AK, LDAK, DWORK( I7+(N+NK)*N+N+1 ),
     $             N+NK )
      CALL DGEMM( 'N', 'N', NK, NK, NP, ONE, DWORK( I2+1 ), NK,
     $            DWORK( I5+1 ), NP, ONE, DWORK( I7+(N+NK)*N+N+1 ),
     $            N+NK )
C
C     Compute the closed-loop poles.
C     Additional workspace:  need 3*(N+NK);  prefer larger.
C     The fact that M > 0, NP > 0, and NK <= N is used here.
C
      CALL DGEES( 'N', 'N', SELECT, N+NK, DWORK( I7+1 ), N+NK, SDIM,
     $            DWORK, DWORK( N+NK+1 ), DWORK( IWRK+1 ), N,
     $            DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C     Check the stability of the closed-loop system.
C
      NS = 0
      DO 40 I = 1, N+NK
         IF( DWORK( I ).GE.ZERO ) NS = NS + 1
   40 CONTINUE
      IF( NS.GT.0 ) THEN
         INFO = 6
         RETURN
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10ID ***
      END
