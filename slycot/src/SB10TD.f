      SUBROUTINE SB10TD( N, M, NP, NCON, NMEAS, D, LDD, TU, LDTU, TY,
     $                   LDTY, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $                   RCOND, TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To compute the matrices of the H2 optimal discrete-time controller
C
C              | AK | BK |
C          K = |----|----|,
C              | CK | DK |
C
C     from the matrices of the controller for the normalized system,
C     as determined by the SLICOT Library routine SB10SD.
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
C             The number of control inputs (M2).  M >= NCON >= 0.
C             NP-NMEAS >= NCON.
C
C     NMEAS   (input) INTEGER
C             The number of measurements (NP2).  NP >= NMEAS >= 0.
C             M-NCON >= NMEAS.
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D. Only the trailing
C             NMEAS-by-NCON submatrix D22 is used.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     TU      (input) DOUBLE PRECISION array, dimension (LDTU,M2)
C             The leading M2-by-M2 part of this array must contain the
C             control transformation matrix TU, as obtained by the
C             SLICOT Library routine SB10PD.
C
C     LDTU    INTEGER
C             The leading dimension of the array TU.  LDTU >= max(1,M2).
C
C     TY      (input) DOUBLE PRECISION array, dimension (LDTY,NP2)
C             The leading NP2-by-NP2 part of this array must contain the
C             measurement transformation matrix TY, as obtained by the
C             SLICOT Library routine SB10PD.
C
C     LDTY    INTEGER
C             The leading dimension of the array TY.
C             LDTY >= max(1,NP2).
C
C     AK      (input/output) DOUBLE PRECISION array, dimension (LDAK,N)
C             On entry, the leading N-by-N part of this array must
C             contain controller state matrix for the normalized system
C             as obtained by the SLICOT Library routine SB10SD.
C             On exit, the leading N-by-N part of this array contains
C             controller state matrix AK.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (input/output) DOUBLE PRECISION array, dimension
C             (LDBK,NMEAS)
C             On entry, the leading N-by-NMEAS part of this array must
C             contain controller input matrix for the normalized system
C             as obtained by the SLICOT Library routine SB10SD.
C             On exit, the leading N-by-NMEAS part of this array
C             contains controller input matrix BK.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (input/output) DOUBLE PRECISION array, dimension (LDCK,N)
C             On entry, the leading NCON-by-N part of this array must
C             contain controller output matrix for the normalized
C             system as obtained by the SLICOT Library routine SB10SD.
C             On exit, the leading NCON-by-N part of this array contains
C             controller output matrix CK.
C
C     LDCK    INTEGER
C             The leading dimension of the array CK.
C             LDCK >= max(1,NCON).
C
C     DK      (input/output) DOUBLE PRECISION array, dimension
C             (LDDK,NMEAS)
C             On entry, the leading NCON-by-NMEAS part of this array
C             must contain controller matrix DK for the normalized
C             system as obtained by the SLICOT Library routine SB10SD.
C             On exit, the leading NCON-by-NMEAS part of this array
C             contains controller input/output matrix DK.
C
C     LDDK    INTEGER
C             The leading dimension of the array DK.
C             LDDK >= max(1,NCON).
C
C     RCOND   (output) DOUBLE PRECISION
C             RCOND contains an estimate of the reciprocal condition
C             number of the matrix Im2 + DKHAT*D22 which must be
C             inverted in the computation of the controller.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used in determining the nonsingularity of the
C             matrix which must be inverted. If TOL <= 0, then a default
C             value equal to sqrt(EPS) is used, where EPS is the
C             relative machine precision.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*M2)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= max(N*M2,N*NP2,M2*NP2,M2*M2+4*M2).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix Im2 + DKHAT*D22 is singular, or the
C                   estimated condition number is larger than or equal
C                   to 1/TOL.
C
C     METHOD
C
C     The routine implements the formulas given in [1].
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
C     input and output transformations and of the matrix Im2 +
C     DKHAT*D22.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999,
C     Jan. 2000.
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
      INTEGER            INFO, LDAK, LDBK, LDCK, LDD, LDDK, LDTU, LDTY,
     $                   LDWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   RCOND, TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AK( LDAK, * ), BK( LDBK, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), DWORK( * ),
     $                   TU( LDTU, * ), TY( LDTY, * )
C     ..
C     .. Local Scalars ..
      INTEGER            INFO2, IWRK, M1, M2, MINWRK, NP1, NP2
      DOUBLE PRECISION   ANORM, TOLL
C     ..
C     .. External Functions
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGECON, DGEMM, DGETRF, DGETRS, DLACPY, DLASET,
     $                   XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
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
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -7
      ELSE IF( LDTU.LT.MAX( 1, M2 ) ) THEN
         INFO = -9
      ELSE IF( LDTY.LT.MAX( 1, NP2 ) ) THEN
         INFO = -11
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDCK.LT.MAX( 1, M2 ) ) THEN
         INFO = -17
      ELSE IF( LDDK.LT.MAX( 1, M2 ) ) THEN
         INFO = -19
      ELSE
C
C        Compute workspace.
C
         MINWRK = MAX ( N*M2, N*NP2, M2*NP2, M2*( M2 + 4 ) )
         IF( LDWORK.LT.MINWRK )
     $      INFO = -24
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10TD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 .OR. M1.EQ.0 .OR. M2.EQ.0
     $    .OR. NP1.EQ.0 .OR. NP2.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
C
      TOLL = TOL
      IF( TOLL.LE.ZERO ) THEN
C
C        Set the default value of the tolerance for nonsingularity test.
C
         TOLL = SQRT( DLAMCH( 'Epsilon' )  )
      END IF
C
C     Find BKHAT .
C
      CALL DGEMM( 'N', 'N', N, NP2, NP2, ONE, BK, LDBK, TY, LDTY, ZERO,
     $            DWORK, N )
      CALL DLACPY ('Full', N, NP2, DWORK, N, BK, LDBK )
C
C     Find CKHAT .
C
      CALL DGEMM( 'N', 'N', M2, N, M2, ONE, TU, LDTU, CK, LDCK, ZERO,
     $            DWORK, M2 )
      CALL DLACPY ('Full', M2, N, DWORK, M2, CK, LDCK )
C
C     Compute DKHAT .
C
      CALL DGEMM( 'N', 'N', M2, NP2, M2, ONE, TU, LDTU, DK, LDDK, ZERO,
     $            DWORK, M2 )
      CALL DGEMM( 'N', 'N', M2, NP2, NP2, ONE, DWORK, M2, TY, LDTY,
     $            ZERO, DK, LDDK )
C
C     Compute Im2 + DKHAT*D22 .
C
      IWRK = M2*M2 + 1
      CALL DLASET( 'Full', M2, M2, ZERO, ONE, DWORK, M2 )
      CALL DGEMM( 'N', 'N', M2, M2, NP2, ONE, DK, LDDK,
     $            D( NP1+1, M1+1 ), LDD, ONE, DWORK, M2 )
      ANORM = DLANGE( '1', M2, M2, DWORK, M2, DWORK( IWRK ) )
      CALL DGETRF( M2, M2, DWORK, M2, IWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      CALL DGECON( '1', M2, DWORK, M2, ANORM, RCOND, DWORK( IWRK ),
     $             IWORK( M2+1 ), INFO2 )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND.LT.TOLL ) THEN
         INFO = 1
         RETURN
      END IF
C
C     Compute CK .
C
      CALL DGETRS( 'N', M2, N, DWORK, M2, IWORK, CK, LDCK, INFO2 )
C
C     Compute DK .
C
      CALL DGETRS( 'N', M2, NP2, DWORK, M2, IWORK, DK, LDDK, INFO2 )
C
C     Compute AK .
C
      CALL DGEMM( 'N', 'N', N, M2, NP2, ONE, BK, LDBK, D( NP1+1, M1+1 ),
     $            LDD, ZERO, DWORK, N )
      CALL DGEMM( 'N', 'N', N, N, M2, -ONE, DWORK, N, CK, LDCK, ONE, AK,
     $            LDAK )
C
C     Compute BK .
C
      CALL DGEMM( 'N', 'N', N, NP2, M2, -ONE, DWORK, N, DK, LDDK,
     $            ONE, BK, LDBK )
      RETURN
C *** Last line of SB10TD ***
      END
