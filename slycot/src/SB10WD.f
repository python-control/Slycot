      SUBROUTINE SB10WD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, F, LDF, H, LDH, TU, LDTU, TY, LDTY,
     $                   AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, INFO )
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
C     from the state feedback matrix F and output injection matrix H as
C     determined by the SLICOT Library routine SB10VD.
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
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B. Only the submatrix
C             B2 = B(:,M-M2+1:M) is used.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading NP-by-N part of this array must contain the
C             system output matrix C. Only the submatrix
C             C2 = C(NP-NP2+1:NP,:) is used.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D. Only the submatrix
C             D22 = D(NP-NP2+1:NP,M-M2+1:M) is used.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     F       (input) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading NCON-by-N part of this array must contain the
C             state feedback matrix F.
C
C     LDF     INTEGER
C             The leading dimension of the array F.  LDF >= max(1,NCON).
C
C     H       (input) DOUBLE PRECISION array, dimension (LDH,NMEAS)
C             The leading N-by-NMEAS part of this array must contain the
C             output injection matrix H.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= max(1,N).
C
C     TU      (input) DOUBLE PRECISION array, dimension (LDTU,M2)
C             The leading M2-by-M2 part of this array must contain the
C             control transformation matrix TU, as obtained by the
C             SLICOT Library routine SB10UD.
C
C     LDTU    INTEGER
C             The leading dimension of the array TU.  LDTU >= max(1,M2).
C
C     TY      (input) DOUBLE PRECISION array, dimension (LDTY,NP2)
C             The leading NP2-by-NP2 part of this array must contain the
C             measurement transformation matrix TY, as obtained by the
C             SLICOT Library routine SB10UD.
C
C     LDTY    INTEGER
C             The leading dimension of the array TY.
C             LDTY >= max(1,NP2).
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
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
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
C     input and output transformations.
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
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                   LDDK, LDF, LDH, LDTU, LDTY, M, N, NCON, NMEAS,
     $                   NP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), F( LDF, * ),
     $                   H( LDH, * ), TU( LDTU, * ), TY( LDTY, * )
C     ..
C     .. Local Scalars ..
      INTEGER            M1, M2, NP1, NP2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLASET, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
      ELSE IF( LDF.LT.MAX( 1, M2 ) ) THEN
         INFO = -15
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDTU.LT.MAX( 1, M2 ) ) THEN
         INFO = -19
      ELSE IF( LDTY.LT.MAX( 1, NP2 ) ) THEN
         INFO = -21
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -23
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -25
      ELSE IF( LDCK.LT.MAX( 1, M2 ) ) THEN
         INFO = -27
      ELSE IF( LDDK.LT.MAX( 1, M2 ) ) THEN
         INFO = -29
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10WD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 .OR. M1.EQ.0 .OR. M2.EQ.0
     $    .OR. NP1.EQ.0 .OR. NP2.EQ.0 ) RETURN
C
C     Compute the transpose of D22*F . BK is used as workspace.
C
      CALL DGEMM( 'T', 'T', N, NP2, M2, ONE, F, LDF, D( NP1+1, M1+1 ),
     $            LDD, ZERO, BK, LDBK )
C
C     Find AK = A + H*C2 + B2*F + H*D22*F .
C
      CALL DLACPY( 'Full', N, N, A, LDA, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, NP2, ONE, H, LDH, C( NP1+1, 1 ), LDC,
     $            ONE, AK, LDAK )
      CALL DGEMM( 'N', 'N', N, N, M2, ONE, B( 1, M1+1 ), LDB,
     $            F, LDF, ONE, AK, LDAK )
      CALL DGEMM( 'N', 'T', N, N, NP2, ONE, H, LDH, BK, LDBK, ONE, AK,
     $            LDAK )
C
C     Find BK = -H*Ty .
C
      CALL DGEMM( 'N', 'N', N, NP2, NP2, -ONE, H, LDH, TY, LDTY, ZERO,
     $            BK, LDBK )
C
C     Find CK = Tu*F .
C
      CALL DGEMM( 'N', 'N', M2, N, M2, ONE, TU, LDTU, F, LDF, ZERO, CK,
     $            LDCK )
C
C     Find DK .
C
      CALL DLASET( 'Full', M2, NP2, ZERO, ZERO, DK, LDDK )
C
      RETURN
C *** Last line of SB10WD ***
      END
