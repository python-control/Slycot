      SUBROUTINE SB10HD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $                   RCOND, TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
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
C     To compute the matrices of the H2 optimal n-state controller
C
C              | AK | BK |
C          K = |----|----|
C              | CK | DK |
C
C     for the system
C
C                   | A  | B1  B2  |   | A | B |
C               P = |----|---------| = |---|---| ,
C                   | C1 |  0  D12 |   | C | D |
C                   | C2 | D21 D22 |
C
C     where B2 has as column size the number of control inputs (NCON)
C     and C2 has as row size the number of measurements (NMEAS) being
C     provided to the controller.
c
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) The block D11 of D is zero,
C
C     (A3) D12 is full column rank and D21 is full row rank.
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
C             system input/output matrix D.
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
C     RCOND   (output) DOUBLE PRECISION array, dimension (4)
C             RCOND(1) contains the reciprocal condition number of the
C                      control transformation matrix;
C             RCOND(2) contains the reciprocal condition number of the
C                      measurement transformation matrix;
C             RCOND(3) contains an estimate of the reciprocal condition
C                      number of the X-Riccati equation;
C             RCOND(4) contains an estimate of the reciprocal condition
C                      number of the Y-Riccati equation.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used for controlling the accuracy of the applied
C             transformations for computing the normalized form in
C             SLICOT Library routine SB10UD. Transformation matrices
C             whose reciprocal condition numbers are less than TOL are
C             not allowed. If TOL <= 0, then a default value equal to
C             sqrt(EPS) is used, where EPS is the relative machine
C             precision.
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
C             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 +
C                       max(max(M2 + NP1*NP1 +
C                               max(NP1*N,3*M2+NP1,5*M2),
C                               NP2 + M1*M1 +
C                               max(M1*N,3*NP2+M1,5*NP2),
C                               N*M2,NP2*N,NP2*M2,1),
C                               N*(14*N+12+M2+NP2)+5),
C             where M1 = M - M2 and NP1 = NP - NP2.
C             For good performance, LDWORK must generally be larger.
C             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is
C             2*Q*(3*Q+2*N)+max(1,Q*(Q+max(N,5)+1),N*(14*N+12+2*Q)+5).
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix D12 had not full column rank in
C                   respect to the tolerance TOL;
C             = 2:  if the matrix D21 had not full row rank in respect
C                   to the tolerance TOL;
C             = 3:  if the singular value decomposition (SVD) algorithm
C                   did not converge (when computing the SVD of one of
C                   the matrices D12 or D21).
C             = 4:  if the X-Riccati equation was not solved
C                   successfully;
C             = 5:  if the Y-Riccati equation was not solved
C                   successfully.
C
C     METHOD
C
C     The routine implements the formulas given in [1], [2].
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
C     NUMERICAL ASPECTS
C
C     The accuracy of the result depends on the condition numbers of the
C     input and output transformations and on the condition numbers of
C     the two Riccati equations, as given by the values of RCOND(1),
C     RCOND(2), RCOND(3) and RCOND(4), respectively.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, Oct. 1998.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999,
C     Sept. 1999, Jan. 2000, Feb. 2000.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, H2 optimal control, optimal regulator,
C     robust control.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                   LDDK, LDWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), DWORK( * ),
     $                   RCOND( 4 )
C     ..
C     .. Local Scalars ..
      INTEGER            INFO2, IWC, IWD, IWF, IWH, IWRK, IWTU, IWTY,
     $                   IWY, LWAMAX, M1, M2, MINWRK, NP1, NP2
      DOUBLE PRECISION   TOLL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLACPY, SB10UD, SB10VD, SB10WD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
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
      ELSE
C
C        Compute workspace.
C
         MINWRK = N*M + NP*(N+M) + M2*M2 + NP2*NP2 +
     $            MAX( MAX( M2 + NP1*NP1 +
     $                      MAX( NP1*N, 3*M2 + NP1, 5*M2 ),
     $                      NP2 + M1*M1 +
     $                      MAX( M1*N, 3*NP2 + M1, 5*NP2 ),
     $                      N*M2, NP2*N, NP2*M2, 1 ),
     $                 N*( 14*N + 12 + M2 + NP2 ) + 5 )
         IF( LDWORK.LT.MINWRK )
     $      INFO = -26
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10HD', -INFO )
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
      TOLL = TOL
      IF( TOLL.LE.ZERO ) THEN
C
C        Set the default value of the tolerance for rank tests.
C
         TOLL = SQRT( DLAMCH( 'Epsilon' ) )
      END IF
C
C     Workspace usage.
C
      IWC  = N*M + 1
      IWD  = IWC + NP*N
      IWTU = IWD + NP*M
      IWTY = IWTU + M2*M2
      IWRK = IWTY + NP2*NP2
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
      CALL DLACPY( 'Full', NP, N, C, LDC, DWORK( IWC ), NP )
      CALL DLACPY( 'Full', NP, M, D, LDD, DWORK( IWD ), NP )
C
C     Transform the system so that D12 and D21 satisfy the formulas
C     in the computation of the H2 optimal controller.
C
      CALL SB10UD( N, M, NP, NCON, NMEAS, DWORK, N, DWORK( IWC ), NP,
     $             DWORK( IWD ), NP, DWORK( IWTU ), M2, DWORK( IWTY ),
     $             NP2, RCOND, TOLL, DWORK( IWRK ), LDWORK-IWRK+1,
     $             INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = INFO2
         RETURN
      END IF
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
      IWY  = IWRK
      IWF  = IWY + N*N
      IWH  = IWF + M2*N
      IWRK = IWH + N*NP2
C
C     Compute the optimal state feedback and output injection matrices.
C     AK is used to store X.
C
      CALL SB10VD( N, M, NP, NCON, NMEAS, A, LDA, DWORK, N,
     $             DWORK( IWC ), NP, DWORK( IWF ), M2, DWORK( IWH ), N,
     $             AK, LDAK, DWORK( IWY ), N, RCOND( 3 ), IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = INFO2 + 3
         RETURN
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Compute the H2 optimal controller.
C
      CALL SB10WD( N, M, NP, NCON, NMEAS, A, LDA, DWORK, N,
     $             DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWF ), M2,
     $             DWORK( IWH ), N, DWORK( IWTU ), M2, DWORK( IWTY ),
     $             NP2, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, INFO2 )
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10HD ***
      END
