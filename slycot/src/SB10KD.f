      SUBROUTINE SB10KD( N, M, NP, A, LDA, B, LDB, C, LDC, FACTOR,
     $                   AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, RCOND,
     $                   IWORK, DWORK, LDWORK, BWORK, INFO )
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
C              | C | 0 |
C
C     in the Discrete-Time Loop Shaping Design Procedure.
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
C     FACTOR  (input) DOUBLE PRECISION
C             = 1  implies that an optimal controller is required;
C             > 1  implies that a suboptimal controller is required
C                  achieving a performance FACTOR less than optimal.
C             FACTOR >= 1.
C
C     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N)
C             The leading N-by-N part of this array contains the
C             controller state matrix Ak.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NP)
C             The leading N-by-NP part of this array contains the
C             controller input matrix Bk.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N)
C             The leading M-by-N part of this array contains the
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
C     RCOND   (output) DOUBLE PRECISION array, dimension (4)
C             RCOND(1) contains an estimate of the reciprocal condition
C                      number of the linear system of equations from
C                      which the solution of the P-Riccati equation is
C                      obtained;
C             RCOND(2) contains an estimate of the reciprocal condition
C                      number of the linear system of equations from
C                      which the solution of the Q-Riccati equation is
C                      obtained;
C             RCOND(3) contains an estimate of the reciprocal condition
C                      number of the linear system of equations from
C                      which the solution of the X-Riccati equation is
C                      obtained;
C             RCOND(4) contains an estimate of the reciprocal condition
C                      number of the matrix Rx + Bx'*X*Bx (see the
C                      comments in the code).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension 2*max(N,NP+M)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 15*N*N + 6*N +
C                       max( 14*N+23, 16*N, 2*N+NP+M, 3*(NP+M) ) +
C                       max( N*N, 11*N*NP + 2*M*M + 8*NP*NP + 8*M*N +
C                                 4*M*NP + NP ).
C             For good performance, LDWORK must generally be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the P-Riccati equation is not solved successfully;
C             = 2:  the Q-Riccati equation is not solved successfully;
C             = 3:  the X-Riccati equation is not solved successfully;
C             = 4:  the iteration to compute eigenvalues failed to
C                   converge;
C             = 5:  the matrix Rx + Bx'*X*Bx is singular;
C             = 6:  the closed-loop system is unstable.
C
C     METHOD
C
C     The routine implements the method presented in [1].
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
C     two Riccati equations solved in the controller design. For
C     better conditioning it is advised to take FACTOR > 1.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 2000.
C
C     REVISIONS
C
C     V. Sima, Katholieke University Leuven, January 2001,
C     February 2001.
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
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDDK,
     $                   LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   DK( LDDK, * ), DWORK( * ), RCOND( 4 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     $                   I11, I12, I13, I14, I15, I16, I17, I18, I19,
     $                   I20, I21, I22, I23, I24, I25, I26, INFO2,
     $                   IWRK, J, LWA, LWAMAX, MINWRK, N2, NS, SDIM
      DOUBLE PRECISION   GAMMA, RNORM
C     ..
C     .. External Functions ..
      LOGICAL            SELECT
      DOUBLE PRECISION   DLANSY, DLAPY2
      EXTERNAL           DLANSY, DLAPY2, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEES, DLACPY, DLASET, DPOTRF, DPOTRS,
     $                   DSYCON, DSYEV, DSYRK, DSYTRF, DSYTRS, SB02OD,
     $                   XERBLA
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
      ELSE IF( FACTOR.LT.ONE ) THEN
         INFO = -10
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDCK.LT.MAX( 1, M ) ) THEN
         INFO = -16
      ELSE IF( LDDK.LT.MAX( 1, M ) ) THEN
         INFO = -18
      END IF
C
C     Compute workspace.
C
      MINWRK = 15*N*N + 6*N + MAX( 14*N+23, 16*N, 2*N+NP+M, 3*(NP+M) ) +
     $         MAX( N*N, 11*N*NP + 2*M*M + 8*NP*NP + 8*M*N +
     $                   4*M*NP + NP )
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10KD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 ) THEN
         RCOND( 1 ) = ONE
         RCOND( 2 ) = ONE
         RCOND( 3 ) = ONE
         RCOND( 4 ) = ONE
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Workspace usage.
C
      N2 = 2*N
      I1 = N*N
      I2 = I1 + N*N
      I3 = I2 + N*N
      I4 = I3 + N*N
      I5 = I4 + N2
      I6 = I5 + N2
      I7 = I6 + N2
      I8 = I7 + N2*N2
      I9 = I8 + N2*N2
C
      IWRK = I9 + N2*N2
      LWAMAX = 0
C
C     Compute Cr = C'*C .
C
      CALL DSYRK( 'U', 'T', N, NP, ONE, C, LDC, ZERO, DWORK( I2+1 ), N )
C
C     Compute Dr = B*B' .
C
      CALL DSYRK( 'U', 'N', N, M, ONE, B, LDB, ZERO, DWORK( I3+1 ), N )
C                                                     -1
C     Solution of the Riccati equation A'*P*(In + Dr*P) *A - P + Cr = 0.
C
      CALL SB02OD( 'D', 'G', 'N', 'U', 'Z', 'S', N, M, NP, A, LDA,
     $             DWORK( I3+1 ), N, DWORK( I2+1 ), N, DWORK, M, DWORK,
     $             N, RCOND( 1 ), DWORK, N, DWORK( I4+1 ),
     $             DWORK( I5+1 ), DWORK( I6+1 ), DWORK( I7+1 ), N2,
     $             DWORK( I8+1 ), N2, DWORK( I9+1 ), N2, -ONE, IWORK,
     $             DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C     Transpose A in AK (used as workspace).
C
      DO 40 J = 1, N
         DO 30 I = 1, N
            AK( I,J ) = A( J,I )
   30    CONTINUE
   40 CONTINUE
C                                                    -1
C     Solution of the Riccati equation A*Q*(In + Cr*Q) *A' - Q + Dr = 0.
C
      CALL SB02OD( 'D', 'G', 'N', 'U', 'Z', 'S', N, M, NP, AK, LDAK,
     $             DWORK( I2+1 ), N, DWORK( I3+1 ), N, DWORK, M, DWORK,
     $             N, RCOND( 2 ), DWORK( I1+1 ), N, DWORK( I4+1 ),
     $             DWORK( I5+1 ), DWORK( I6+1 ), DWORK( I7+1 ), N2,
     $             DWORK( I8+1 ), N2, DWORK( I9+1 ), N2, -ONE, IWORK,
     $             DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 2
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C     Compute gamma.
C
      CALL DGEMM( 'N', 'N', N, N, N, ONE, DWORK( I1+1 ), N, DWORK, N,
     $            ZERO, AK, LDAK )
      CALL DGEES( 'N', 'N', SELECT, N, AK, LDAK, SDIM, DWORK( I6+1 ),
     $            DWORK( I7+1 ), DWORK( IWRK+1 ), N, DWORK( IWRK+1 ),
     $            LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 4
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
      GAMMA = ZERO
      DO 50 I = 1, N
         GAMMA = MAX( GAMMA, DWORK( I6+I ) )
   50 CONTINUE
      GAMMA = FACTOR*SQRT( ONE + GAMMA )
C
C     Workspace usage.
C
      I3  = I2  + N*NP
      I4  = I3  + NP*NP
      I5  = I4  + NP*NP
      I6  = I5  + NP*NP
      I7  = I6  + NP
      I8  = I7  + NP*NP
      I9  = I8  + NP*NP
      I10 = I9  + NP*NP
      I11 = I10 + N*NP
      I12 = I11 + N*NP
      I13 = I12 + ( NP+M )*( NP+M )
      I14 = I13 + N*( NP+M )
      I15 = I14 + N*( NP+M )
      I16 = I15 + N*N
      I17 = I16 + N2
      I18 = I17 + N2
      I19 = I18 + N2
      I20 = I19 + ( N2+NP+M )*( N2+NP+M )
      I21 = I20 + ( N2+NP+M )*N2
C
      IWRK = I21 + N2*N2
C
C     Compute Q*C' .
C
      CALL DGEMM( 'N', 'T', N, NP, N, ONE, DWORK( I1+1 ), N, C, LDC,
     $            ZERO, DWORK( I2+1 ), N )
C
C     Compute Ip + C*Q*C' .
C
      CALL DLASET( 'Full', NP, NP, ZERO, ONE, DWORK( I3+1 ), NP )
      CALL DGEMM( 'N', 'N', NP, NP, N, ONE, C, LDC, DWORK( I2+1 ), N,
     $            ONE, DWORK( I3+1 ), NP )
C
C     Compute the eigenvalues and eigenvectors of Ip + C'*Q*C
C
      CALL DLACPY( 'U', NP, NP, DWORK( I3+1 ), NP, DWORK( I5+1 ), NP )
      CALL DSYEV( 'V', 'U', NP, DWORK( I5+1 ), NP, DWORK( I6+1 ),
     $            DWORK( IWRK+1 ), LDWORK-IWRK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 4
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C                            -1
C     Compute ( Ip + C'*Q*C )  .
C
      DO 70 J = 1, NP
         DO 60 I = 1, NP
            DWORK( I9+I+(J-1)*NP ) = DWORK( I5+J+(I-1)*NP ) /
     $                               DWORK( I6+I )
   60    CONTINUE
   70 CONTINUE
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I5+1 ), NP,
     $            DWORK( I9+1 ), NP, ZERO, DWORK( I4+1 ), NP )
C
C     Compute Z2 .
C
      DO 90 J = 1, NP
         DO 80 I = 1, NP
            DWORK( I9+I+(J-1)*NP ) = DWORK( I5+J+(I-1)*NP ) /
     $                               SQRT( DWORK( I6+I ) )
   80    CONTINUE
   90 CONTINUE
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I5+1 ), NP,
     $            DWORK( I9+1 ), NP, ZERO, DWORK( I7+1 ), NP )
C               -1
C     Compute Z2  .
C
      DO 110 J = 1, NP
         DO 100 I = 1, NP
            DWORK( I9+I+(J-1)*NP ) = DWORK( I5+J+(I-1)*NP )*
     $                               SQRT( DWORK( I6+I ) )
  100    CONTINUE
  110 CONTINUE
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I5+1 ), NP,
     $            DWORK( I9+1 ), NP, ZERO, DWORK( I8+1 ), NP )
C
C     Compute A*Q*C' .
C
      CALL DGEMM( 'N', 'N', N, NP, N, ONE, A, LDA, DWORK( I2+1 ), N,
     $            ZERO, DWORK( I10+1 ), N )
C                                        -1
C     Compute H = -A*Q*C'*( Ip + C*Q*C' )  .
C
      CALL DGEMM( 'N', 'N', N, NP, NP, -ONE, DWORK( I10+1 ), N,
     $            DWORK( I4+1 ), NP, ZERO, DWORK( I11+1 ), N )
C
C     Compute Rx .
C
      CALL DLASET( 'F', NP+M, NP+M, ZERO, ONE, DWORK( I12+1 ), NP+M )
      DO 130 J = 1, NP
         DO 120 I = 1, NP
            DWORK( I12+I+(J-1)*(NP+M) ) = DWORK( I3+I+(J-1)*NP )
  120    CONTINUE
         DWORK( I12+J+(J-1)*(NP+M) ) = DWORK( I3+J+(J-1)*NP ) -
     $                                 GAMMA*GAMMA
  130 CONTINUE
C
C     Compute Bx .
C
      CALL DGEMM( 'N', 'N', N, NP, NP, -ONE, DWORK( I11+1 ), N,
     $            DWORK( I8+1 ), NP, ZERO, DWORK( I13+1 ), N )
      DO 150 J = 1, M
         DO 140 I = 1, N
            DWORK( I13+N*NP+I+(J-1)*N ) = B( I, J )
  140    CONTINUE
  150 CONTINUE
C
C     Compute Sx .
C
      CALL DGEMM( 'T', 'N', N, NP, NP, ONE, C, LDC, DWORK( I8+1 ), NP,
     $            ZERO, DWORK( I14+1 ), N )
      CALL DLASET( 'F', N, M, ZERO, ZERO, DWORK( I14+N*NP+1 ), N )
C
C     Solve the Riccati equation
C                                                      -1
C       X = A'*X*A + Cx - (Sx + A'*X*Bx)*(Rx + Bx'*X*B ) *(Sx'+Bx'*X*A).
C
      CALL SB02OD( 'D', 'B', 'C', 'U', 'N', 'S', N, NP+M, NP, A, LDA,
     $             DWORK( I13+1 ), N, C, LDC, DWORK( I12+1 ), NP+M,
     $             DWORK( I14+1 ), N, RCOND( 3 ), DWORK( I15+1 ), N,
     $             DWORK( I16+1 ), DWORK( I17+1 ), DWORK( I18+1 ),
     $             DWORK( I19+1 ), N2+NP+M, DWORK( I20+1 ), N2+NP+M,
     $             DWORK( I21+1 ), N2, -ONE, IWORK, DWORK( IWRK+1 ),
     $             LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
      I22 = I16
      I23 = I22 + ( NP+M )*N
      I24 = I23 + ( NP+M )*( NP+M )
      I25 = I24 + ( NP+M )*N
      I26 = I25 + M*N
C
      IWRK = I25
C
C     Compute Bx'*X .
C
      CALL DGEMM( 'T', 'N', NP+M, N, N, ONE, DWORK( I13+1 ), N,
     $            DWORK( I15+1 ), N, ZERO, DWORK( I22+1 ), NP+M )
C
C     Compute Rx + Bx'*X*Bx .
C
      CALL DLACPY( 'F', NP+M, NP+M, DWORK( I12+1 ), NP+M,
     $             DWORK( I23+1 ), NP+M )
      CALL DGEMM( 'N', 'N', NP+M, NP+M, N, ONE, DWORK( I22+1 ), NP+M,
     $            DWORK( I13+1 ), N, ONE, DWORK( I23+1 ), NP+M )
C
C     Compute -( Sx' + Bx'*X*A ) .
C
      DO 170 J = 1, N
         DO 160 I = 1, NP+M
            DWORK( I24+I+(J-1)*(NP+M) ) = DWORK( I14+J+(I-1)*N )
 160     CONTINUE
 170  CONTINUE
      CALL DGEMM( 'N', 'N', NP+M, N, N, -ONE, DWORK( I22+1 ), NP+M,
     $            A, LDA, -ONE, DWORK( I24+1 ), NP+M )
C
C     Factorize Rx + Bx'*X*Bx .
C
      RNORM = DLANSY( '1', 'U', NP+M, DWORK( I23+1 ), NP+M,
     $                DWORK( IWRK+1 ) )
      CALL DSYTRF( 'U', NP+M, DWORK( I23+1 ), NP+M, IWORK,
     $             DWORK( IWRK+1 ), LDWORK-IWRK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 5
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
      CALL DSYCON( 'U', NP+M, DWORK( I23+1 ), NP+M, IWORK, RNORM,
     $             RCOND( 4 ), DWORK( IWRK+1 ), IWORK( NP+M+1), INFO2 )
C                                   -1
C     Compute F = -( Rx + Bx'*X*Bx )  ( Sx' + Bx'*X*A ) .
C
      CALL DSYTRS( 'U', NP+M, N, DWORK( I23+1 ), NP+M, IWORK,
     $             DWORK( I24+1 ), NP+M, INFO2 )
C
C     Compute B'*X .
C
      CALL DGEMM( 'T', 'N', M, N, N, ONE, B, LDB, DWORK( I15+1 ), N,
     $            ZERO, DWORK( I25+1 ), M )
C
C     Compute Im + B'*X*B .
C
      CALL DLASET( 'F', M, M, ZERO, ONE, DWORK( I23+1 ), M )
      CALL DGEMM( 'N', 'N', M, M, N, ONE, DWORK( I25+1 ), M, B, LDB,
     $            ONE, DWORK( I23+1 ), M )
C
C     Factorize Im + B'*X*B .
C
      CALL DPOTRF( 'U', M, DWORK( I23+1 ), M, INFO2 )
C                            -1
C     Compute ( Im + B'*X*B )  B'*X .
C
      CALL DPOTRS( 'U', M, N, DWORK( I23+1 ), M, DWORK( I25+1 ), M,
     $             INFO2 )
C                                 -1
C     Compute Dk = ( Im + B'*X*B )  B'*X*H .
C
      CALL DGEMM( 'N', 'N', M, NP, N, ONE, DWORK( I25+1 ), M,
     $            DWORK( I11+1 ), N, ZERO, DK, LDDK )
C
C     Compute Bk = -H + B*Dk .
C
      CALL DLACPY( 'F', N, NP, DWORK( I11+1 ), N, BK, LDBK )
      CALL DGEMM( 'N', 'N', N, NP, M, ONE, B, LDB, DK, LDDK, -ONE,
     $            BK, LDBK )
C                  -1
C     Compute Dk*Z2  .
C
      CALL DGEMM( 'N', 'N', M, NP, NP, ONE, DK, LDDK, DWORK( I8+1 ),
     $            NP, ZERO, DWORK( I26+1 ), M )
C
C     Compute F1 + Z2*C .
C
      CALL DLACPY( 'F', NP, N, DWORK( I24+1 ), NP+M, DWORK( I12+1 ),
     $             NP )
      CALL DGEMM( 'N', 'N', NP, N, NP, ONE, DWORK( I7+1 ), NP, C, LDC,
     $            ONE, DWORK( I12+1 ), NP )
C                            -1
C     Compute Ck = F2 - Dk*Z2  *( F1 + Z2*C ) .
C
      CALL DLACPY( 'F', M, N, DWORK( I24+NP+1 ), NP+M, CK, LDCK )
      CALL DGEMM( 'N', 'N', M, N, NP, -ONE, DWORK( I26+1 ), M,
     $            DWORK( I12+1 ), NP, ONE, CK, LDCK )
C
C     Compute Ak = A + H*C + B*Ck .
C
      CALL DLACPY( 'F', N, N, A, LDA, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, NP, ONE, DWORK( I11+1 ), N, C, LDC,
     $            ONE, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, M, ONE, B, LDB, CK, LDCK, ONE, AK,
     $            LDAK )
C
C     Workspace usage.
C
      I1 = M*N
      I2 = I1 + N2*N2
      I3 = I2 + N2
C
      IWRK = I3 + N2
C
C     Compute Dk*C .
C
      CALL DGEMM( 'N', 'N', M, N, NP, ONE, DK, LDDK, C, LDC, ZERO,
     $            DWORK, M )
C
C     Compute the closed-loop state matrix.
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I1+1 ), N2 )
      CALL DGEMM( 'N', 'N', N, N, M, -ONE, B, LDB, DWORK, M, ONE,
     $            DWORK( I1+1 ), N2 )
      CALL DGEMM( 'N', 'N', N, N, NP, -ONE, BK, LDBK, C, LDC, ZERO,
     $            DWORK( I1+N+1 ), N2 )
      CALL DGEMM( 'N', 'N', N, N, M, ONE, B, LDB, CK, LDCK, ZERO,
     $            DWORK( I1+N2*N+1 ), N2 )
      CALL DLACPY( 'F', N, N, AK, LDAK, DWORK( I1+N2*N+N+1 ), N2 )
C
C     Compute the closed-loop poles.
C
      CALL DGEES( 'N', 'N', SELECT, N2, DWORK( I1+1 ), N2, SDIM,
     $            DWORK( I2+1 ), DWORK( I3+1 ), DWORK( IWRK+1 ), N,
     $            DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 4
         RETURN
      END IF
      LWA = INT( DWORK( IWRK+1 ) ) + IWRK
      LWAMAX = MAX( LWA, LWAMAX )
C
C     Check the stability of the closed-loop system.
C
      NS = 0
      DO 180 I = 1, N2
        IF( DLAPY2( DWORK( I2+I ), DWORK( I3+I ) ).GT.ONE ) NS = NS + 1
  180 CONTINUE
      IF( NS.GT.0 ) THEN
         INFO = 6
         RETURN
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10KD ***
      END
