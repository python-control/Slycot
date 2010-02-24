      SUBROUTINE SB10VD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   F, LDF, H, LDH, X, LDX, Y, LDY, XYCOND, IWORK,
     $                   DWORK, LDWORK, BWORK, INFO )
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
C     To compute the state feedback and the output injection
C     matrices for an H2 optimal n-state controller for the system
C
C                   | A  | B1  B2  |   | A | B |
C               P = |----|---------| = |---|---|
C                   | C1 |  0  D12 |   | C | D |
C                   | C2 | D21 D22 |
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
C          SLICOT Library routine SB10UD. Matrix D is not used
C          explicitly.
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
C     F       (output) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading NCON-by-N part of this array contains the
C             state feedback matrix F.
C
C     LDF     INTEGER
C             The leading dimension of the array F.  LDF >= max(1,NCON).
C
C     H       (output) DOUBLE PRECISION array, dimension (LDH,NMEAS)
C             The leading N-by-NMEAS part of this array contains the
C             output injection matrix H.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= max(1,N).
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
C     XYCOND  (output) DOUBLE PRECISION array, dimension (2)
C             XYCOND(1) contains an estimate of the reciprocal condition
C                       number of the X-Riccati equation;
C             XYCOND(2) contains an estimate of the reciprocal condition
C                       number of the Y-Riccati equation.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(2*N,N*N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal
C             LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 13*N*N + 12*N + 5.
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
C             = 2:  if the Y-Riccati equation was not solved
C                   successfully.
C
C     METHOD
C
C     The routine implements the formulas given in [1], [2]. The X-
C     and Y-Riccati equations are solved with condition and accuracy
C     estimates [3].
C
C     REFERENCES
C
C     [1] Zhou, K., Doyle, J.C., and Glover, K.
C         Robust and Optimal Control.
C         Prentice-Hall, Upper Saddle River, NJ, 1996.
C
C     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and
C         Smith, R.
C         mu-Analysis and Synthesis Toolbox.
C         The MathWorks Inc., Natick, Mass., 1995.
C
C     [3] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V.
C         DGRSVX and DMSRIC: Fortan 77 subroutines for solving
C         continuous-time matrix algebraic Riccati equations with
C         condition and accuracy estimates.
C         Preprint SFB393/98-16, Fak. f. Mathematik, Tech. Univ.
C         Chemnitz, May 1998.
C
C     NUMERICAL ASPECTS
C
C     The precision of the solution of the matrix Riccati equations
C     can be controlled by the values of the condition numbers
C     XYCOND(1) and XYCOND(2) of these equations.
C
C     FURTHER COMMENTS
C
C     The Riccati equations are solved by the Schur approach
C     implementing condition and accuracy estimates.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999.
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
      INTEGER            INFO, LDA, LDB, LDC, LDF, LDH, LDWORK, LDX,
     $                   LDY, M, N, NCON, NMEAS, NP
C     ..
C     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK( * ),  F( LDF, * ), H( LDH, * ),
     $                   X( LDX, * ), XYCOND( 2 ), Y( LDY, * )
C     ..
C     .. Local Scalars ..
      INTEGER            INFO2, IWG, IWI, IWQ, IWR, IWRK, IWS, IWT, IWV,
     $                   LWAMAX, M1, M2, MINWRK, N2, ND1, ND2, NP1, NP2
      DOUBLE PRECISION   FERR, SEP
C     ..
C     .. External Functions ..
C
      DOUBLE PRECISION   DLANSY
      EXTERNAL           DLANSY
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, DSYRK, SB02RD, XERBLA
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
      ELSE IF( LDF.LT.MAX( 1, NCON ) ) THEN
         INFO = -13
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDY.LT.MAX( 1, N ) ) THEN
         INFO = -19
      ELSE
C
C        Compute workspace.
C
         MINWRK = 13*N*N + 12*N + 5
         IF( LDWORK.LT.MINWRK )
     $      INFO = -23
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10VD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 .OR. M1.EQ.0 .OR. M2.EQ.0
     $    .OR. NP1.EQ.0 .OR. NP2.EQ.0 ) THEN
         DWORK( 1 )  = ONE
         XYCOND( 1 ) = ONE
         XYCOND( 2 ) = ONE
         RETURN
      END IF
C
      ND1 = NP1 - M2
      ND2 = M1 - NP2
      N2  = 2*N
C
C     Workspace usage.
C
      IWQ  = N*N + 1
      IWG  = IWQ + N*N
      IWT  = IWG + N*N
      IWV  = IWT + N*N
      IWR  = IWV + N*N
      IWI  = IWR + N2
      IWS  = IWI + N2
      IWRK = IWS + 4*N*N
C
C     Compute Ax = A - B2*D12'*C1 .
C
      CALL DLACPY ('Full', N, N, A, LDA, DWORK, N )
      CALL DGEMM( 'N', 'N', N, N, M2, -ONE, B( 1, M1+1 ), LDB,
     $            C( ND1+1, 1), LDC, ONE, DWORK, N )
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
C     Solution of the Riccati equation Ax'*X + X*Ax + Cx - X*Dx*X = 0 .
C     Workspace:  need   13*N*N + 12*N + 5;
C                 prefer larger.
C
      CALL SB02RD( 'All', 'Continuous', 'NotUsed', 'NoTranspose',
     $             'Lower', 'GeneralScaling', 'Stable', 'NotFactored',
     $             'Original', N, DWORK, N, DWORK( IWT ), N,
     $             DWORK( IWV ), N, DWORK( IWG ), N, DWORK( IWQ ), N,
     $             X, LDX, SEP, XYCOND( 1 ), FERR, DWORK( IWR ),
     $             DWORK( IWI ), DWORK( IWS ), N2, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
C     Compute F = -D12'*C1 - B2'*X .
C
      CALL DLACPY( 'Full', M2, N, C( ND1+1, 1 ), LDC, F, LDF )
      CALL DGEMM( 'T', 'N', M2, N, N, -ONE, B( 1, M1+1 ), LDB, X, LDX,
     $            -ONE, F, LDF )
C
C     Compute Ay = A - B1*D21'*C2 .
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK, N )
      CALL DGEMM( 'N', 'N', N, N, NP2, -ONE, B( 1, ND2+1 ), LDB,
     $            C( NP1+1, 1 ), LDC, ONE, DWORK, N )
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
C     Solution of the Riccati equation Ay*Y + Y*Ay' + Cy - Y*Dy*Y = 0 .
C     Workspace:  need   13*N*N + 12*N + 5;
C                 prefer larger.
C
      CALL SB02RD( 'All', 'Continuous', 'NotUsed', 'Transpose',
     $             'Upper', 'GeneralScaling', 'Stable', 'NotFactored',
     $             'Original', N, DWORK, N, DWORK( IWT ), N,
     $             DWORK( IWV ), N, DWORK( IWG ), N, DWORK( IWQ ), N,
     $             Y, LDY, SEP, XYCOND( 2 ), FERR, DWORK( IWR ),
     $             DWORK( IWI ), DWORK( IWS ), N2, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
C
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Compute H = -B1*D21' - Y*C2' .
C
      CALL DLACPY( 'Full', N, NP2, B( 1, ND2+1 ), LDB, H, LDH )
      CALL DGEMM( 'N', 'T', N, NP2, N, -ONE, Y, LDY, C( NP1+1, 1 ), LDC,
     $            -ONE, H, LDH )
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10VD ***
      END
