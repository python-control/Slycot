      SUBROUTINE SB10ZD( N, M, NP, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   FACTOR, AK, LDAK, BK, LDBK, CK, LDCK, DK,
     $                   LDDK, RCOND, TOL, IWORK, DWORK, LDWORK, BWORK,
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
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D of the shaped plant.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     FACTOR  (input) DOUBLE PRECISION
C             = 1  implies that an optimal controller is required
C                  (not recommended);
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
C     RCOND   (output) DOUBLE PRECISION array, dimension (6)
C             RCOND(1) contains an estimate of the reciprocal condition
C                      number of the linear system of equations from
C                      which the solution of the P-Riccati equation is
C                      obtained;
C             RCOND(2) contains an estimate of the reciprocal condition
C                      number of the linear system of equations from
C                      which the solution of the Q-Riccati equation is
C                      obtained;
C             RCOND(3) contains an estimate of the reciprocal condition
C                      number of the matrix (gamma^2-1)*In - P*Q;
C             RCOND(4) contains an estimate of the reciprocal condition
C                      number of the matrix Rx + Bx'*X*Bx;
C             RCOND(5) contains an estimate of the reciprocal condition
C                                                  ^
C                      number of the matrix Ip + D*Dk;
C             RCOND(6) contains an estimate of the reciprocal condition
C                                                ^
C                      number of the matrix Im + Dk*D.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used for checking the nonsingularity of the
C             matrices to be inverted. If TOL <= 0, then a default value
C             equal to sqrt(EPS) is used, where EPS is the relative
C             machine precision.  TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension 2*max(N,M+NP)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 16*N*N + 5*M*M + 7*NP*NP + 6*M*N + 7*M*NP +
C                        7*N*NP + 6*N + 2*(M + NP) +
C                        max(14*N+23,16*N,2*M-1,2*NP-1).
C             For good performance, LDWORK must generally be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    (output) INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  the P-Riccati equation is not solved successfully;
C             =  2:  the Q-Riccati equation is not solved successfully;
C             =  3:  the iteration to compute eigenvalues or singular
C                    values failed to converge;
C             =  4:  the matrix (gamma^2-1)*In - P*Q is singular;
C             =  5:  the matrix Rx + Bx'*X*Bx is singular;
C                                      ^
C             =  6:  the matrix Ip + D*Dk is singular;
C                                    ^
C             =  7:  the matrix Im + Dk*D is singular;
C             =  8:  the matrix Ip - D*Dk is singular;
C             =  9:  the matrix Im - Dk*D is singular;
C             = 10:  the closed-loop system is unstable.
C
C     METHOD
C
C     The routine implements the formulas given in [1].
C
C     REFERENCES
C
C     [1] Gu, D.-W., Petkov, P.H., and Konstantinov, M.M.
C         On discrete H-infinity loop shaping design procedure routines.
C         Technical Report 00-6, Dept. of Engineering, Univ. of
C         Leicester, UK, 2000.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the results depends on the conditioning of the
C     two Riccati equations solved in the controller design. For
C     better conditioning it is advised to take FACTOR > 1.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2001.
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
     $                   LDDK, LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR, TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A ( LDA,  * ), AK( LDAK, * ), B ( LDB,  * ),
     $                   BK( LDBK, * ), C ( LDC,  * ), CK( LDCK, * ),
     $                   D ( LDD,  * ), DK( LDDK, * ), DWORK( * ),
     $                   RCOND( 6 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     $                   I11, I12, I13, I14, I15, I16, I17, I18, I19,
     $                   I20, I21, I22, I23, I24, I25, I26, INFO2, IWRK,
     $                   J, LWAMAX, MINWRK, N2, NS, SDIM
      DOUBLE PRECISION   ANORM, GAMMA, TOLL
C     ..
C     .. External Functions ..
      LOGICAL            SELECT
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY, DLAPY2
      EXTERNAL           DLAMCH, DLANGE, DLANSY, DLAPY2, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGECON, DGEES, DGEMM, DGETRF, DGETRS,
     $                   DLACPY, DLASCL, DLASET, DPOTRF, DPOTRS, DSWAP,
     $                   DSYCON, DSYEV, DSYRK, DSYTRF, DSYTRS, DTRSM,
     $                   DTRTRS, MA02AD, MB01RX, MB02VD, SB02OD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, SQRT
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
         INFO = -14
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF( LDCK.LT.MAX( 1, M ) ) THEN
         INFO = -18
      ELSE IF( LDDK.LT.MAX( 1, M ) ) THEN
         INFO = -20
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -22
      END IF
C
C     Compute workspace.
C
      MINWRK = 16*N*N + 5*M*M + 7*NP*NP + 6*M*N + 7*M*NP + 7*N*NP +
     $         6*N + 2*(M + NP) + MAX( 14*N+23, 16*N, 2*M-1, 2*NP-1 )
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -25
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10ZD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C     Note that some computation could be made if one or two of the
C     dimension parameters N, M, and P are zero, but the results are
C     not so meaningful.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 ) THEN
         RCOND( 1 ) = ONE
         RCOND( 2 ) = ONE
         RCOND( 3 ) = ONE
         RCOND( 4 ) = ONE
         RCOND( 5 ) = ONE
         RCOND( 6 ) = ONE
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Set the default tolerance, if needed.
C
      IF( TOL.LE.ZERO ) THEN
         TOLL = SQRT( DLAMCH( 'Epsilon' ) )
      ELSE
         TOLL = TOL
      END IF
C
C     Workspace usage.
C
      N2  = 2*N
      I1  = 1   + N*N
      I2  = I1  + N*N
      I3  = I2  + NP*NP
      I4  = I3  + M*M
      I5  = I4  + NP*NP
      I6  = I5  + M*M
      I7  = I6  + M*N
      I8  = I7  + M*N
      I9  = I8  + N*N
      I10 = I9  + N*N
      I11 = I10 + N2
      I12 = I11 + N2
      I13 = I12 + N2
      I14 = I13 + N2*N2
      I15 = I14 + N2*N2
C
      IWRK   = I15 + N2*N2
      LWAMAX = 0
C
C     Compute R1 = Ip + D*D' .
C
      CALL DLASET( 'U', NP, NP, ZERO, ONE, DWORK( I2 ), NP )
      CALL DSYRK(  'U', 'N', NP, M, ONE, D, LDD, ONE, DWORK( I2 ), NP )
      CALL DLACPY( 'U', NP, NP, DWORK( I2 ), NP, DWORK( I4 ), NP )
C
C     Factorize R1 = R'*R .
C
      CALL DPOTRF( 'U', NP, DWORK( I4 ), NP, INFO2 )
C                 -1
C     Compute C'*R   in BK .
C
      CALL MA02AD( 'F', NP, N, C, LDC, BK, LDBK )
      CALL DTRSM(  'R', 'U', 'N', 'N', N, NP, ONE, DWORK( I4 ), NP, BK,
     $             LDBK )
C
C     Compute R2 = Im + D'*D .
C
      CALL DLASET( 'U', M, M, ZERO, ONE, DWORK( I3 ), M )
      CALL DSYRK(  'U', 'T', M, NP, ONE, D, LDD, ONE, DWORK( I3 ), M )
      CALL DLACPY( 'U', M, M, DWORK( I3 ), M, DWORK( I5 ), M )
C
C     Factorize R2 = U'*U .
C
      CALL DPOTRF( 'U', M, DWORK( I5 ), M, INFO2 )
C               -1
C     Compute (U  )'*B' .
C
      CALL MA02AD( 'F', N, M, B, LDB, DWORK( I6 ), M )
      CALL DTRTRS( 'U', 'T', 'N', M, N, DWORK( I5 ), M, DWORK( I6 ), M,
     $             INFO2 )
C
C     Compute D'*C .
C
      CALL DGEMM( 'T', 'N', M, N, NP, ONE, D, LDD, C, LDC, ZERO,
     $            DWORK( I7 ), M )
C               -1
C     Compute (U  )'*D'*C .
C
      CALL DTRTRS( 'U', 'T', 'N', M, N, DWORK( I5 ), M, DWORK( I7 ), M,
     $             INFO2 )
C                          -1
C     Compute Ar = A - B*R2  D'*C .
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I8 ), N )
      CALL DGEMM(  'T', 'N', N, N, M, -ONE, DWORK( I6 ), M, DWORK( I7 ),
     $             M, ONE, DWORK( I8 ), N )
C                       -1
C     Compute Cr = C'*R1  *C .
C
      CALL DSYRK( 'U', 'N', N, NP, ONE, BK, LDBK, ZERO, DWORK( I9 ), N )
C                      -1
C     Compute Dr = B*R2  B' in AK .
C
      CALL DSYRK( 'U', 'T', N, M, ONE, DWORK( I6 ), M, ZERO, AK, LDAK )
C                                                       -1
C     Solution of the Riccati equation Ar'*P*(In + Dr*P)  Ar - P +
C                                              Cr = 0 .
      CALL SB02OD( 'D', 'G', 'N', 'U', 'Z', 'S', N, M, NP, DWORK( I8 ),
     $             N, AK, LDAK, DWORK( I9 ), N, DWORK, M, DWORK, N,
     $             RCOND( 1 ), DWORK, N, DWORK( I10 ), DWORK( I11 ),
     $             DWORK( I12 ), DWORK( I13 ), N2, DWORK( I14 ), N2,
     $             DWORK( I15 ), N2, -ONE, IWORK, DWORK( IWRK ),
     $             LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Transpose Ar .
C
      DO 10 J = 1, N - 1
         CALL DSWAP( J, DWORK( I8+J ), N, DWORK( I8+J*N ), 1 )
   10 CONTINUE
C                                                      -1
C     Solution of the Riccati equation Ar*Q*(In + Cr*Q)  *Ar' - Q +
C                                             Dr = 0 .
      CALL SB02OD( 'D', 'G', 'N', 'U', 'Z', 'S', N, M, NP, DWORK( I8 ),
     $             N, DWORK( I9 ), N, AK, LDAK, DWORK, M, DWORK, N,
     $             RCOND( 2 ), DWORK( I1 ), N, DWORK( I10 ),
     $             DWORK( I11 ), DWORK( I12 ), DWORK( I13 ), N2,
     $             DWORK( I14 ), N2, DWORK( I15 ), N2, -ONE, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 2
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Compute gamma.
C
      CALL DGEMM( 'N', 'N', N, N, N, ONE, DWORK( I1 ), N, DWORK, N,
     $            ZERO, DWORK( I8 ), N )
      CALL DGEES( 'N', 'N', SELECT, N, DWORK( I8 ), N, SDIM,
     $            DWORK( I10 ), DWORK( I11 ), DWORK( IWRK ), N,
     $            DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
      GAMMA  = ZERO
C
      DO 20 I = 0, N - 1
         GAMMA = MAX( GAMMA, DWORK( I10+I ) )
   20 CONTINUE
C
      GAMMA = FACTOR*SQRT( ONE + GAMMA )
C
C     Workspace usage.
C
      I5  = I4  + NP*NP
      I6  = I5  + M*M
      I7  = I6  + NP*NP
      I8  = I7  + NP*NP
      I9  = I8  + NP*NP
      I10 = I9  + NP
      I11 = I10 + NP*NP
      I12 = I11 + M*M
      I13 = I12 + M
C
      IWRK = I13 + M*M
C
C     Compute the eigenvalues and eigenvectors of R1 .
C
      CALL DLACPY( 'U', NP, NP, DWORK( I2 ), NP, DWORK( I8 ), NP )
      CALL DSYEV(  'V', 'U', NP, DWORK( I8 ), NP, DWORK( I9 ),
     $             DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C               -1/2
C     Compute R1     .
C
      DO 40 J = 1, NP
         DO 30 I = 1, NP
            DWORK( I10-1+I+(J-1)*NP ) = DWORK( I8-1+J+(I-1)*NP ) /
     $                                  SQRT( DWORK( I9+I-1 ) )
   30    CONTINUE
   40 CONTINUE
C
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I8 ), NP,
     $            DWORK( I10 ), NP, ZERO, DWORK( I4 ), NP )
C
C     Compute the eigenvalues and eigenvectors of R2 .
C
      CALL DLACPY( 'U', M, M, DWORK( I3 ), M, DWORK( I11 ), M )
      CALL DSYEV(  'V', 'U', M, DWORK( I11 ), M, DWORK( I12 ),
     $             DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C               -1/2
C     Compute R2     .
C
      DO 60 J = 1, M
         DO 50 I = 1, M
            DWORK( I13-1+I+(J-1)*M ) = DWORK( I11-1+J+(I-1)*M ) /
     $                                 SQRT( DWORK( I12+I-1 ) )
   50    CONTINUE
   60 CONTINUE
C
      CALL DGEMM( 'N', 'N', M, M, M, ONE, DWORK( I11 ), M, DWORK( I13 ),
     $            M, ZERO, DWORK( I5 ), M )
C
C     Compute R1 + C*Q*C' .
C
      CALL DGEMM(  'N', 'T', N, NP, N, ONE, DWORK( I1 ), N, C, LDC,
     $             ZERO, BK, LDBK )
      CALL MB01RX( 'L', 'U', 'N', NP, N, ONE, ONE, DWORK( I2 ), NP,
     $             C, LDC, BK, LDBK, INFO2 )
      CALL DLACPY( 'U', NP, NP, DWORK( I2 ), NP, DWORK( I8 ), NP )
C
C     Compute the eigenvalues and eigenvectors of R1 + C*Q*C' .
C
      CALL DSYEV( 'V', 'U', NP, DWORK( I8 ), NP, DWORK( I9 ),
     $            DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C                            -1
C     Compute ( R1 + C*Q*C' )   .
C
      DO 80 J = 1, NP
         DO 70 I = 1, NP
            DWORK( I10-1+I+(J-1)*NP ) = DWORK( I8-1+J+(I-1)*NP ) /
     $                                  DWORK( I9+I-1 )
   70    CONTINUE
   80 CONTINUE
C
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I8 ), NP,
     $            DWORK( I10 ), NP, ZERO, DWORK( I6 ), NP )
C               -1
C     Compute Z2  .
C
      DO 100 J = 1, NP
         DO 90 I = 1, NP
            DWORK( I10-1+I+(J-1)*NP ) = DWORK( I8-1+J+(I-1)*NP )*
     $                                  SQRT( DWORK( I9+I-1 ) )
   90    CONTINUE
  100 CONTINUE
C
      CALL DGEMM( 'N', 'N', NP, NP, NP, ONE, DWORK( I8 ), NP,
     $            DWORK( I10 ), NP, ZERO, DWORK( I7 ), NP )
C
C     Workspace usage.
C
      I9  = I8  + N*NP
      I10 = I9  + N*NP
      I11 = I10 + NP*M
      I12 = I11 + ( NP + M )*( NP + M )
      I13 = I12 + N*( NP + M )
      I14 = I13 + N*( NP + M )
      I15 = I14 + N*N
      I16 = I15 + N*N
      I17 = I16 + ( NP + M )*N
      I18 = I17 + ( NP + M )*( NP + M )
      I19 = I18 + ( NP + M )*N
      I20 = I19 + M*N
      I21 = I20 + M*NP
      I22 = I21 + NP*N
      I23 = I22 + N*N
      I24 = I23 + N*NP
      I25 = I24 + NP*NP
      I26 = I25 + M*M
C
      IWRK = I26 + N*M
C
C     Compute A*Q*C' + B*D' .
C
      CALL DGEMM( 'N', 'T', N, NP, M, ONE, B, LDB, D, LDD, ZERO,
     $            DWORK( I8 ), N )
      CALL DGEMM( 'N', 'N', N, NP, N, ONE, A, LDA, BK, LDBK,
     $            ONE, DWORK( I8 ), N )
C                                                 -1
C     Compute H = -( A*Q*C'+B*D' )*( R1 + C*Q*C' )   .
C
      CALL DGEMM( 'N', 'N', N, NP, NP, -ONE, DWORK( I8 ), N,
     $            DWORK( I6 ), NP, ZERO, DWORK( I9 ), N )
C               -1/2
C     Compute R1    D .
C
      CALL DGEMM( 'N', 'N', NP, M, NP, ONE, DWORK( I4 ), NP, D, LDD,
     $            ZERO, DWORK( I10 ), NP )
C
C     Compute Rx .
C
      DO 110 J = 1, NP
         CALL DCOPY( J, DWORK( I2+(J-1)*NP ), 1,
     $                  DWORK( I11+(J-1)*(NP+M) ), 1 )
         DWORK( I11-1+J+(J-1)*(NP+M) ) = DWORK( I2-1+J+(J-1)*NP ) -
     $                                   GAMMA*GAMMA
  110 CONTINUE
C
      CALL DGEMM(  'N', 'N', NP, M, NP, ONE, DWORK( I7 ), NP,
     $             DWORK( I10 ), NP, ZERO, DWORK( I11+(NP+M)*NP ),
     $             NP+M )
      CALL DLASET( 'U', M, M, ZERO, ONE, DWORK( I11+(NP+M)*NP+NP ),
     $             NP+M )
C
C     Compute Bx .
C
      CALL DGEMM( 'N', 'N', N, NP, NP, -ONE, DWORK( I9 ), N,
     $            DWORK( I7 ), NP, ZERO, DWORK( I12 ), N )
      CALL DGEMM( 'N', 'N', N, M, M, ONE, B, LDB, DWORK( I5 ), M,
     $            ZERO, DWORK( I12+N*NP ), N )
C
C     Compute Sx .
C
      CALL DGEMM( 'T', 'N', N, NP, NP, ONE, C, LDC, DWORK( I7 ), NP,
     $            ZERO, DWORK( I13 ), N )
      CALL DGEMM( 'T', 'N', N, M, NP, ONE, C, LDC, DWORK( I10 ), NP,
     $            ZERO, DWORK( I13+N*NP ), N )
C
C     Compute  (gamma^2 - 1)*In - P*Q .
C
      CALL DLASET( 'F', N, N, ZERO, GAMMA*GAMMA-ONE, DWORK( I14 ), N )
      CALL DGEMM(  'N', 'N', N, N, N, -ONE, DWORK, N, DWORK( I1 ), N,
     $             ONE, DWORK( I14 ), N )
C                                          -1
C     Compute X =  ((gamma^2 - 1)*In - P*Q)  *gamma^2*P .
C
      CALL DLACPY( 'F', N, N, DWORK, N, DWORK( I15 ), N )
      CALL DLASCL( 'G', 0, 0, ONE, GAMMA*GAMMA, N, N, DWORK( I15 ), N,
     $             INFO )
      ANORM = DLANGE( '1', N, N, DWORK( I14 ), N, DWORK( IWRK ) )
      CALL DGETRF( N, N, DWORK( I14 ), N, IWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 4
         RETURN
      END IF
      CALL DGECON( '1', N, DWORK( I14 ), N, ANORM, RCOND( 3 ),
     $             DWORK( IWRK ), IWORK( N+1 ), INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 3 ).LT.TOLL ) THEN
         INFO = 4
         RETURN
      END IF
      CALL DGETRS( 'N', N, N, DWORK( I14 ), N, IWORK, DWORK( I15 ),
     $             N, INFO2 )
C
C     Compute Bx'*X .
C
      CALL DGEMM( 'T', 'N', NP+M, N, N, ONE, DWORK( I12 ), N,
     $            DWORK( I15 ), N, ZERO, DWORK( I16 ), NP+M )
C
C     Compute Rx + Bx'*X*Bx .
C
      CALL DLACPY( 'U', NP+M, NP+M, DWORK( I11 ), NP+M, DWORK( I17 ),
     $             NP+M )
      CALL MB01RX( 'L', 'U', 'N', NP+M, N, ONE, ONE, DWORK( I17 ), NP+M,
     $             DWORK( I16 ), NP+M, DWORK( I12 ), N, INFO2 )
C
C     Compute  -( Sx' + Bx'*X*A ) .
C
      CALL MA02AD( 'F', N, NP+M, DWORK( I13 ), N, DWORK( I18 ), NP+M )
      CALL DGEMM(  'N', 'N', NP+M, N, N, -ONE, DWORK( I16 ), NP+M,
     $             A, LDA, -ONE, DWORK( I18 ), NP+M )
C
C     Factorize Rx + Bx'*X*Bx .
C
      ANORM = DLANSY( '1', 'U', NP+M, DWORK( I17 ), NP+M,
     $                DWORK( IWRK ) )
      CALL DSYTRF( 'U', NP+M, DWORK( I17 ), NP+M, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 5
         RETURN
      END IF
      CALL DSYCON( 'U', NP+M, DWORK( I17 ), NP+M, IWORK, ANORM,
     $             RCOND( 4 ), DWORK( IWRK ), IWORK( NP+M+1), INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 4 ).LT.TOLL ) THEN
         INFO = 5
         RETURN
      END IF
C                                   -1
C     Compute F = -( Rx + Bx'*X*Bx )  ( Sx' + Bx'*X*A ) .
C
      CALL DSYTRS( 'U', NP+M, N, DWORK( I17 ), NP+M, IWORK,
     $             DWORK( I18 ), NP+M, INFO2 )
C
C     Compute B'*X .
C
      CALL DGEMM( 'T', 'N', M, N, N, ONE, B, LDB, DWORK( I15 ), N,
     $            ZERO, DWORK( I19 ), M )
C
C     Compute  -( D' - B'*X*H ) .
C
      DO 130 J = 1, NP
         DO 120 I = 1, M
            DWORK( I20-1+I+(J-1)*M ) = -D( J, I )
  120    CONTINUE
  130 CONTINUE
C
      CALL DGEMM( 'N', 'N', M, NP, N, ONE, DWORK( I19 ), M,
     $            DWORK( I9 ), N, ONE, DWORK( I20 ), M )
C                   -1
C     Compute C + Z2  *F1 .
C
      CALL DLACPY( 'F', NP, N, C, LDC, DWORK( I21 ), NP )
      CALL DGEMM(  'N', 'N', NP, N, NP, ONE, DWORK( I7 ), NP,
     $             DWORK( I18 ), NP+M, ONE, DWORK( I21 ), NP )
C
C     Compute R2 + B'*X*B .
C
      CALL MB01RX( 'L', 'U', 'N', M, N, ONE, ONE, DWORK( I3 ), M,
     $             DWORK( I19 ), M, B, LDB, INFO2 )
C
C     Factorize R2 + B'*X*B .
C
      CALL DPOTRF( 'U', M, DWORK( I3 ), M, INFO2 )
C             ^                    -1
C     Compute Dk = -( R2 + B'*X*B )  (D' - B'*X*H) .
C
      CALL DLACPY( 'F', M, NP, DWORK( I20 ), M, DK, LDDK )
      CALL DPOTRS( 'U', M, NP, DWORK( I3 ), M, DK, LDDK, INFO2 )
C             ^           ^
C     Compute Bk = -H + B*Dk .
C
      CALL DLACPY( 'F', N, NP, DWORK( I9 ), N, DWORK( I23 ), N )
      CALL DGEMM(  'N', 'N', N, NP, M, ONE, B, LDB, DK, LDDK,
     $             -ONE, DWORK( I23 ), N )
C               -1/2
C     Compute R2    *F2  .
C
      CALL DGEMM( 'N', 'N', M, N, M, ONE, DWORK( I5 ), M,
     $            DWORK( I18+NP ), NP+M, ZERO, CK, LDCK )
C             ^      -1/2      ^          -1
C     Compute Ck = R2    *F2 - Dk*( C + Z2  *F1 ) .
C
      CALL DGEMM( 'N', 'N', M, N, NP, -ONE, DK, LDDK,
     $            DWORK( I21 ), NP, ONE, CK, LDCK )
C             ^                ^
C     Compute Ak = A + H*C + B*Ck .
C
      CALL DLACPY( 'F', N, N, A, LDA, AK, LDAK )
      CALL DGEMM(  'N', 'N', N, N, NP, ONE, DWORK( I9 ), N, C, LDC,
     $             ONE, AK, LDAK )
      CALL DGEMM(  'N', 'N', N, N, M, ONE, B, LDB, CK, LDCK,
     $             ONE, AK, LDAK )
C                    ^
C     Compute Ip + D*Dk .
C
      CALL DLASET( 'Full', NP, NP, ZERO, ONE, DWORK( I24 ), NP )
      CALL DGEMM(  'N', 'N', NP, NP, M, ONE, D, LDD, DK, LDDK,
     $             ONE, DWORK( I24 ), NP )
C                   ^
C     Compute  Im + Dk*D .
C
      CALL DLASET( 'Full', M, M, ZERO, ONE, DWORK( I25 ), M )
      CALL DGEMM(  'N', 'N', M, M, NP, ONE, DK, LDDK, D, LDD,
     $             ONE, DWORK( I25 ), M )
C                  ^ ^    ^         ^    -1
C     Compute Ck = M*Ck,  M = (Im + Dk*D)   .
C
      ANORM = DLANGE( '1', M, M, DWORK( I25 ), M, DWORK( IWRK ) )
      CALL DGETRF( M, M, DWORK( I25 ), M, IWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 7
         RETURN
      END IF
      CALL DGECON( '1', M, DWORK( I25 ), M, ANORM, RCOND( 6 ),
     $             DWORK( IWRK ), IWORK( M+1 ), INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 6 ).LT.TOLL ) THEN
         INFO = 7
         RETURN
      END IF
      CALL DGETRS( 'N', M, N, DWORK( I25 ), M, IWORK, CK, LDCK, INFO2 )
C                  ^ ^
C     Compute Dk = M*Dk .
C
      CALL DGETRS( 'N', M, NP, DWORK( I25 ), M, IWORK, DK, LDDK, INFO2 )
C             ^
C     Compute Bk*D .
C
      CALL DGEMM( 'N', 'N', N, M, NP, ONE, DWORK( I23 ), N, D, LDD,
     $            ZERO, DWORK( I26 ), N )
C                  ^    ^
C     Compute Ak = Ak - Bk*D*Ck.
C
      CALL DGEMM( 'N', 'N', N, N, M, -ONE, DWORK( I26 ), N, CK, LDCK,
     $            ONE, AK, LDAK )
C                  ^          ^  -1
C     Compute Bk = Bk*(Ip + D*Dk)   .
C
      ANORM = DLANGE( '1', NP, NP, DWORK( I24 ), NP, DWORK( IWRK ) )
      CALL DLACPY( 'Full', N, NP, DWORK( I23 ), N, BK, LDBK )
      CALL MB02VD( 'N', N, NP, DWORK( I24 ), NP, IWORK, BK, LDBK,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 6
         RETURN
      END IF
      CALL DGECON( '1', NP, DWORK( I24 ), NP, ANORM, RCOND( 5 ),
     $             DWORK( IWRK ), IWORK( NP+1 ), INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 5 ).LT.TOLL ) THEN
         INFO = 6
         RETURN
      END IF
C
C     Workspace usage.
C
      I2 = 1  + NP*NP
      I3 = I2 + N*NP
      I4 = I3 + M*M
      I5 = I4 + N*M
      I6 = I5 + NP*N
      I7 = I6 + M*N
      I8 = I7 + N2*N2
      I9 = I8 + N2
C
      IWRK = I9 + N2
C
C     Compute Ip - D*Dk .
C
      CALL DLASET( 'Full', NP, NP, ZERO, ONE, DWORK, NP )
      CALL DGEMM(  'N', 'N', NP, NP, M, -ONE, D, LDD, DK, LDDK, ONE,
     $              DWORK, NP )
C                         -1
C     Compute Bk*(Ip-D*Dk)   .
C
      CALL DLACPY( 'Full', N, NP, BK, LDBK, DWORK( I2 ), N )
      CALL MB02VD( 'N', N, NP, DWORK, NP, IWORK, DWORK( I2 ), N, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 8
         RETURN
      END IF
C
C     Compute Im - Dk*D .
C
      CALL DLASET( 'Full', M, M, ZERO, ONE, DWORK( I3 ), M )
      CALL DGEMM(  'N', 'N', M, M, NP, -ONE, DK, LDDK, D, LDD, ONE,
     $              DWORK( I3 ), M )
C                        -1
C     Compute B*(Im-Dk*D)    .
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK( I4 ), N )
      CALL MB02VD( 'N', N, M, DWORK( I3 ), M, IWORK, DWORK( I4 ), N,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 9
         RETURN
      END IF
C
C     Compute D*Ck .
C
      CALL DGEMM( 'N', 'N', NP, N, M, ONE, D, LDD, CK, LDCK, ZERO,
     $            DWORK( I5 ), NP )
C
C     Compute Dk*C .
C
      CALL DGEMM( 'N', 'N', M, N, NP, ONE, DK, LDDK, C, LDC, ZERO,
     $            DWORK( I6 ), M )
C
C     Compute the closed-loop state matrix.
C
      CALL DLACPY( 'F', N, N, A, LDA, DWORK( I7 ), N2 )
      CALL DGEMM(  'N', 'N', N, N, M, ONE, DWORK( I4 ), N,
     $             DWORK( I6 ), M, ONE, DWORK( I7 ), N2 )
      CALL DGEMM(  'N', 'N', N, N, M, ONE, DWORK( I4 ), N, CK, LDCK,
     $             ZERO, DWORK( I7+N2*N ), N2 )
      CALL DGEMM(  'N', 'N', N, N, NP, ONE, DWORK( I2 ), N, C, LDC,
     $             ZERO, DWORK( I7+N ), N2 )
      CALL DLACPY( 'F', N, N, AK, LDAK, DWORK( I7+N2*N+N ), N2 )
      CALL DGEMM(  'N', 'N', N, N, NP, ONE, DWORK( I2 ), N,
     $             DWORK( I5 ), NP, ONE, DWORK( I7+N2*N+N ), N2 )
C
C     Compute the closed-loop poles.
C
      CALL DGEES( 'N', 'N', SELECT, N2, DWORK( I7 ), N2, SDIM,
     $            DWORK( I8 ), DWORK( I9 ), DWORK( IWRK ), N,
     $            DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
C     Check the stability of the closed-loop system.
C
      NS = 0
C
      DO 140 I = 0, N2 - 1
         IF( DLAPY2( DWORK( I8+I ), DWORK( I9+I ) ).GT.ONE )
     $      NS = NS + 1
  140 CONTINUE
C
      IF( NS.GT.0 ) THEN
         INFO = 10
         RETURN
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10ZD ***
      END
