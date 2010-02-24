      SUBROUTINE SB10SD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $                   X, LDX, Y, LDY, RCOND, TOL, IWORK, DWORK,
     $                   LDWORK, BWORK, INFO )
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
C     To compute the matrices of the H2 optimal controller
C
C              | AK | BK |
C          K = |----|----|,
C              | CK | DK |
C
C     for the normalized discrete-time system
C
C                   | A  | B1  B2  |   | A | B |
C               P = |----|---------| = |---|---|
C                   | C1 | D11 D12 |   | C | D |
C                   | C2 | D21  0  |
C
C     where B2 has as column size the number of control inputs (NCON)
C     and C2 has as row size the number of measurements (NMEAS) being
C     provided to the controller.
C
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) D12 is full column rank with D12 = | 0 | and D21 is
C                                             | I |
C          full row rank with D21 = | 0 I | as obtained by the
C          SLICOT Library routine SB10PD,
C
C               j*Theta
C     (A3) | A-e       *I  B2  | has full column rank for all
C          |    C1         D12 |
C
C          0 <= Theta < 2*Pi ,
C
C
C               j*Theta
C     (A4) | A-e       *I  B1  | has full row rank for all
C          |    C2         D21 |
C
C          0 <= Theta < 2*Pi .
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     NP      (input) INTEGER
C             The row size of the matrix C.  NP >= 0.
C
C     NCON    (input) INTEGER
C             The number of control inputs (M2).  M >= NCON >= 0,
C             NP-NMEAS >= NCON.
C
C     NMEAS   (input) INTEGER
C             The number of measurements (NP2).  NP >= NMEAS >= 0,
C             M-NCON >= NMEAS.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading NP-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D. Only the leading
C             (NP-NP2)-by-(M-M2) submatrix D11 is used.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N)
C             The leading N-by-N part of this array contains the
C             controller state matrix AK.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS)
C             The leading N-by-NMEAS part of this array contains the
C             controller input matrix BK.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N)
C             The leading NCON-by-N part of this array contains the
C             controller output matrix CK.
C
C     LDCK    INTEGER
C             The leading dimension of the array CK.
C             LDCK >= max(1,NCON).
C
C     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS)
C             The leading NCON-by-NMEAS part of this array contains the
C             controller input/output matrix DK.
C
C     LDDK    INTEGER
C             The leading dimension of the array DK.
C             LDDK >= max(1,NCON).
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N)
C             The leading N-by-N part of this array contains the matrix
C             X, solution of the X-Riccati equation.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= max(1,N).
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY,N)
C             The leading N-by-N part of this array contains the matrix
C             Y, solution of the Y-Riccati equation.
C
C     LDY     INTEGER
C             The leading dimension of the array Y.  LDY >= max(1,N).
C
C     RCOND   (output) DOUBLE PRECISION array, dimension (4)
C             RCOND contains estimates of the reciprocal condition
C             numbers of the matrices which are to be inverted and the
C             reciprocal condition numbers of the Riccati equations
C             which have to be solved during the computation of the
C             controller. (See the description of the algorithm in [2].)
C             RCOND(1) contains the reciprocal condition number of the
C                      matrix Im2 + B2'*X2*B2;
C             RCOND(2) contains the reciprocal condition number of the
C                      matrix Ip2 + C2*Y2*C2';
C             RCOND(3) contains the reciprocal condition number of the
C                      X-Riccati equation;
C             RCOND(4) contains the reciprocal condition number of the
C                      Y-Riccati equation.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used in determining the nonsingularity of the
C             matrices which must be inverted. If TOL <= 0, then a
C             default value equal to sqrt(EPS) is used, where EPS is the
C             relative machine precision.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(M2,2*N,N*N,NP2)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal
C             LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= max(1, 14*N*N+6*N+max(14*N+23,16*N),
C                              M2*(N+M2+max(3,M1)), NP2*(N+NP2+3)),
C             where M1 = M - M2.
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
C             = 1:  if the X-Riccati equation was not solved
C                   successfully;
C             = 2:  if the matrix Im2 + B2'*X2*B2 is not positive
C                   definite, or it is numerically singular (with
C                   respect to the tolerance TOL);
C             = 3:  if the Y-Riccati equation was not solved
C                   successfully;
C             = 4:  if the matrix Ip2 + C2*Y2*C2' is not positive
C                   definite, or it is numerically singular (with
C                   respect to the tolerance TOL).
C
C     METHOD
C
C     The routine implements the formulas given in [1]. The X- and
C     Y-Riccati equations are solved with condition estimates.
C
C     REFERENCES
C
C     [1] Zhou, K., Doyle, J.C., and Glover, K.
C         Robust and Optimal Control.
C         Prentice-Hall, Upper Saddle River, NJ, 1996.
C
C     [2] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M.
C         Fortran 77 routines for Hinf and H2 design of linear
C         discrete-time control systems.
C         Report 99-8, Department of Engineering, Leicester University,
C         April 1999.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the result depends on the condition numbers of the
C     matrices which are to be inverted and on the condition numbers of
C     the matrix Riccati equations which are to be solved in the
C     computation of the controller. (The corresponding reciprocal
C     condition numbers are given in the output array RCOND.)
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999,
C     January 2003.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, H2 optimal control, LQG, LQR, optimal
C     regulator, robust control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                   LDDK, LDWORK, LDX, LDY, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), DWORK( * ),
     $                   RCOND( * ), X( LDX, * ), Y( LDY, * )
      LOGICAL            BWORK( * )
C     ..
C     .. Local Scalars ..
      INTEGER            INFO2, IW2, IWB, IWC, IWG, IWI, IWQ, IWR, IWRK,
     $                   IWS, IWT, IWU, IWV, J, LWAMAX, M1, M2, MINWRK,
     $                   ND1, ND2, NP1, NP2
      DOUBLE PRECISION   ANORM, FERR, RCOND2, SEPD, TOLL
C     ..
C     .. External functions ..
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           DLAMCH, DLANSY
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, DPOCON, DPOTRF, DPOTRS,
     $                   DSWAP, DSYRK, DTRSM, MB01RX, SB02OD, SB02SD,
     $                   XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      M1  = M - NCON
      M2  = NCON
      NP1 = NP - NMEAS
      NP2 = NMEAS
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( NP.LT.0 ) THEN
         INFO = -3
      ELSE IF( NCON.LT.0 .OR. M1.LT.0 .OR. M2.GT.NP1 ) THEN
         INFO = -4
      ELSE IF( NMEAS.LT.0 .OR. NP1.LT.0 .OR. NP2.GT.M1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, NP ) ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -13
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDCK.LT.MAX( 1, M2 ) ) THEN
         INFO = -19
      ELSE IF( LDDK.LT.MAX( 1, M2 ) ) THEN
         INFO = -21
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -23
      ELSE IF( LDY.LT.MAX( 1, N ) ) THEN
         INFO = -25
      ELSE
C
C        Compute workspace.
C
         MINWRK = MAX( 1, 14*N*N + 6*N + MAX( 14*N + 23, 16*N ),
     $            M2*( N + M2 + MAX( 3, M1 ) ), NP2*( N + NP2 + 3 ) )
         IF( LDWORK.LT.MINWRK )
     $      INFO = -30
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10SD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 .OR. M1.EQ.0 .OR. M2.EQ.0
     $    .OR. NP1.EQ.0 .OR. NP2.EQ.0 ) THEN
         RCOND( 1 ) = ONE
         RCOND( 2 ) = ONE
         RCOND( 3 ) = ONE
         RCOND( 4 ) = ONE
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
      ND1  = NP1 - M2
      ND2  = M1 - NP2
      TOLL = TOL
      IF( TOLL.LE.ZERO ) THEN
C
C        Set the default value of the tolerance for nonsingularity test.
C
         TOLL = SQRT( DLAMCH( 'Epsilon' )  )
      END IF
C
C     Workspace usage.
C
      IWQ  = 1
      IWG  = IWQ + N*N
      IWR  = IWG + N*N
      IWI  = IWR + 2*N
      IWB  = IWI + 2*N
      IWS  = IWB + 2*N
      IWT  = IWS + 4*N*N
      IWU  = IWT + 4*N*N
      IWRK = IWU + 4*N*N
      IWC  = IWR
      IWV  = IWC + N*N
C
C     Compute Ax = A - B2*D12'*C1 in AK .
C
      CALL DLACPY( 'Full', N, N, A, LDA, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, M2, -ONE, B( 1, M1+1 ), LDB,
     $            C( ND1+1, 1), LDC, ONE, AK, LDAK )
C
C     Compute Cx = C1'*C1 - C1'*D12*D12'*C1 .
C
      IF( ND1.GT.0 ) THEN
         CALL DSYRK( 'L', 'T', N, ND1, ONE, C, LDC, ZERO, DWORK( IWQ ),
     $               N )
      ELSE
         CALL DLASET( 'L', N, N, ZERO, ZERO, DWORK( IWQ ), N )
      END IF
C
C     Compute Dx = B2*B2' .
C
      CALL DSYRK( 'L', 'N', N, M2, ONE, B( 1, M1+1 ), LDB, ZERO,
     $            DWORK( IWG ), N )
C
C     Solution of the discrete-time Riccati equation
C        Ax'*inv(In + X2*Dx)*X2*Ax - X2 + Cx  = 0 .
C     Workspace:  need   14*N*N + 6*N + max(14*N+23,16*N);
C                 prefer larger.
C
      CALL SB02OD( 'D', 'G', 'N', 'L', 'Z', 'S', N, M2, NP1, AK, LDAK,
     $             DWORK( IWG ), N, DWORK( IWQ ), N, DWORK( IWRK ), M,
     $             DWORK( IWRK ), N, RCOND2, X, LDX, DWORK( IWR ),
     $             DWORK( IWI ), DWORK( IWB ), DWORK( IWS ), 2*N,
     $             DWORK( IWT ), 2*N, DWORK( IWU ), 2*N, TOLL, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
C     Condition estimation.
C     Workspace:  need   4*N*N + max(N*N+5*N,max(3,2*N*N)+N*N);
C                 prefer larger.
C
      IWRK = IWV + N*N
      CALL SB02SD( 'C', 'N', 'N', 'L', 'O', N, AK, LDAK, DWORK( IWC ),
     $             N, DWORK( IWV ), N, DWORK( IWG ), N, DWORK( IWQ ), N,
     $             X, LDX, SEPD, RCOND( 3 ), FERR, IWORK, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
      IF( INFO2.GT.0 ) RCOND( 3 ) = ZERO
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Workspace usage.
C
      IW2  = M2*N + 1
      IWRK = IW2 + M2*M2
C
C     Compute B2'*X2 .
C
      CALL DGEMM( 'T', 'N', M2, N, N, ONE, B( 1, M1+1 ), LDB, X, LDX,
     $            ZERO, DWORK, M2 )
C
C     Compute Im2 + B2'*X2*B2 .
C
      CALL DLASET( 'L', M2, M2, ZERO, ONE, DWORK( IW2 ), M2 )
      CALL MB01RX( 'Left', 'Lower', 'N', M2, N, ONE, ONE, DWORK( IW2 ),
     $            M2, DWORK, M2, B( 1, M1+1 ), LDB, INFO2 )
C
C     Compute the Cholesky factorization of Im2 + B2'*X2*B2 .
C     Workspace:  need   M2*N + M2*M2 + max(3*M2,M2*M1);
C                 prefer larger.
C
      ANORM = DLANSY( 'I', 'L', M2, DWORK( IW2 ), M2, DWORK( IWRK ) )
      CALL DPOTRF( 'L', M2, DWORK( IW2 ), M2, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
      CALL DPOCON( 'L', M2, DWORK( IW2 ), M2, ANORM, RCOND( 1 ),
     $             DWORK( IWRK ), IWORK, INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 1 ).LT.TOLL ) THEN
         INFO = 2
         RETURN
      END IF
C
C     Compute -( B2'*X2*A + D12'*C1 ) in CK .
C
      CALL DLACPY( 'Full', M2, N, C( ND1+1, 1 ), LDC, CK, LDCK )
      CALL DGEMM( 'N', 'N', M2, N, N, -ONE, DWORK, M2, A, LDA, -ONE, CK,
     $            LDCK )
C
C     Compute F2 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*A + D12'*C1 ) .
C
      CALL DPOTRS( 'L', M2, N, DWORK( IW2 ), M2, CK, LDCK, INFO2 )
C
C     Compute -( B2'*X2*B1 + D12'*D11 ) .
C
      CALL DLACPY( 'Full', M2, M1, D( ND1+1, 1 ), LDD, DWORK( IWRK ),
     $              M2 )
      CALL DGEMM( 'N', 'N', M2, M1, N, -ONE, DWORK, M2, B, LDB, -ONE,
     $            DWORK( IWRK ), M2 )
C
C     Compute F0 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*B1 + D12'*D11 ) .
C
      CALL DPOTRS( 'L', M2, M1, DWORK( IW2 ), M2, DWORK( IWRK ), M2,
     $             INFO2 )
C
C     Save F0*D21' in DK .
C
      CALL DLACPY( 'Full', M2, NP2, DWORK( IWRK+ND2*M2 ), M2, DK,
     $             LDDK )
C
C     Workspace usage.
C
      IWRK = IWU + 4*N*N
C
C     Compute Ay = A - B1*D21'*C2 in AK .
C
      CALL DLACPY( 'Full', N, N, A, LDA, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, NP2, -ONE, B( 1, ND2+1 ), LDB,
     $            C( NP1+1, 1 ), LDC, ONE, AK, LDAK )
C
C     Transpose Ay in-situ.
C
      DO 20 J = 1, N - 1
         CALL DSWAP( J, AK( J+1, 1 ), LDAK, AK( 1, J+1 ), 1 )
   20 CONTINUE
C
C     Compute Cy = B1*B1' - B1*D21'*D21*B1' .
C
      IF( ND2.GT.0 ) THEN
         CALL DSYRK( 'U', 'N', N, ND2, ONE, B, LDB, ZERO, DWORK( IWQ ),
     $               N )
      ELSE
         CALL DLASET( 'U', N, N, ZERO, ZERO, DWORK( IWQ ), N )
      END IF
C
C     Compute Dy = C2'*C2 .
C
      CALL DSYRK( 'U', 'T', N, NP2, ONE, C( NP1+1, 1 ), LDC, ZERO,
     $            DWORK( IWG ), N )
C
C     Solution of the discrete-time Riccati equation
C        Ay*inv( In + Y2*Dy )*Y2*Ay' - Y2 + Cy = 0 .
C
      CALL SB02OD( 'D', 'G', 'N', 'U', 'Z', 'S', N, NP2, M1, AK, LDAK,
     $             DWORK( IWG ), N, DWORK( IWQ ), N, DWORK( IWRK ), M,
     $             DWORK( IWRK ), N, RCOND2, Y, LDY, DWORK( IWR ),
     $             DWORK( IWI ), DWORK( IWB ), DWORK( IWS ), 2*N,
     $             DWORK( IWT ), 2*N, DWORK( IWU ), 2*N, TOLL, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Condition estimation.
C
      IWRK = IWV + N*N
      CALL SB02SD( 'C', 'N', 'N', 'U', 'O', N, AK, LDAK, DWORK( IWC ),
     $             N, DWORK( IWV ), N, DWORK( IWG ), N, DWORK( IWQ ), N,
     $             Y, LDY, SEPD, RCOND( 4 ), FERR, IWORK, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
      IF( INFO2.GT.0 ) RCOND( 4 ) = ZERO
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Workspace usage.
C
      IW2  = N*NP2 + 1
      IWRK = IW2 + NP2*NP2
C
C     Compute Y2*C2' .
C
      CALL DGEMM( 'N', 'T', N, NP2, N, ONE, Y, LDY, C( NP1+1, 1 ), LDC,
     $            ZERO, DWORK, N )
C
C     Compute Ip2 + C2*Y2*C2' .
C
      CALL DLASET( 'U', NP2, NP2, ZERO, ONE, DWORK( IW2 ), NP2 )
      CALL MB01RX( 'Left', 'Upper', 'N', NP2, N, ONE, ONE, DWORK( IW2 ),
     $            NP2, C( NP1+1, 1 ), LDC, DWORK, N, INFO2 )
C
C     Compute the Cholesky factorization of Ip2 + C2*Y2*C2' .
C
      ANORM = DLANSY( 'I', 'U', NP2, DWORK( IW2 ), NP2, DWORK( IWRK ) )
      CALL DPOTRF( 'U', NP2, DWORK( IW2 ), NP2, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 4
         RETURN
      END IF
      CALL DPOCON( 'U', NP2, DWORK( IW2 ), NP2, ANORM, RCOND( 2 ),
     $             DWORK( IWRK ), IWORK, INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND( 2 ).LT.TOLL ) THEN
         INFO = 4
         RETURN
      END IF
C
C     Compute A*Y2*C2' + B1*D21' in BK .
C
      CALL DLACPY ( 'Full', N, NP2, B( 1, ND2+1 ), LDB, BK, LDBK )
      CALL DGEMM( 'N', 'N', N, NP2, N, ONE, A, LDA, DWORK, N, ONE,
     $            BK, LDBK )
C
C     Compute L2 = -( A*Y2*C2' + B1*D21' )*inv( Ip2 + C2*Y2*C2' ) .
C
      CALL DTRSM( 'R', 'U', 'N', 'N', N, NP2, -ONE, DWORK( IW2 ), NP2,
     $            BK, LDBK )
      CALL DTRSM( 'R', 'U', 'T', 'N', N, NP2, ONE, DWORK( IW2 ), NP2,
     $            BK, LDBK )
C
C     Compute F2*Y2*C2' + F0*D21' .
C
      CALL DGEMM( 'N', 'N', M2, NP2, N, ONE, CK, LDCK, DWORK, N, ONE,
     $            DK, LDDK )
C
C     Compute DK = L0 = ( F2*Y2*C2' + F0*D21' )*inv( Ip2 + C2*Y2*C2' ) .
C
      CALL DTRSM( 'R', 'U', 'N', 'N', M2, NP2, ONE, DWORK( IW2 ), NP2,
     $            DK, LDDK )
      CALL DTRSM( 'R', 'U', 'T', 'N', M2, NP2, ONE, DWORK( IW2 ), NP2,
     $            DK, LDDK )
C
C     Compute CK = F2 - L0*C2 .
C
      CALL DGEMM( 'N', 'N', M2, N, NP2, -ONE, DK, LDDK, C( NP1+1, 1),
     $            LDC, ONE, CK, LDCK )
C
C     Find AK = A + B2*( F2 - L0*C2 ) + L2*C2 .
C
      CALL DLACPY( 'Full', N, N, A, LDA, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, M2, ONE, B(1, M1+1 ), LDB, CK, LDCK,
     $            ONE, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, NP2, ONE, BK, LDBK, C( NP1+1, 1),
     $            LDC, ONE, AK, LDAK )
C
C     Find BK = -L2 + B2*L0 .
C
      CALL DGEMM( 'N', 'N', N, NP2, M2, ONE, B( 1, M1+1 ), LDB, DK,
     $            LDDK, -ONE, BK, LDBK )
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10SD ***
      END
