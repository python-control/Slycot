      SUBROUTINE SB10LD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $                   AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC, IWORK,
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
C     To compute the matrices of the closed-loop system
C
C              | AC | BC |
C          G = |----|----|,
C              | CC | DC |
C
C     from the matrices of the open-loop system
C
C               | A | B |
C           P = |---|---|
C               | C | D |
C
C     and the matrices of the controller
C
C              | AK | BK |
C          K = |----|----|.
C              | CK | DK |
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
C     AK      (input) DOUBLE PRECISION array, dimension (LDAK,N)
C             The leading N-by-N part of this array must contain the
C             controller state matrix AK.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (input) DOUBLE PRECISION array, dimension (LDBK,NMEAS)
C             The leading N-by-NMEAS part of this array must contain the
C             controller input matrix BK.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (input) DOUBLE PRECISION array, dimension (LDCK,N)
C             The leading NCON-by-N part of this array must contain the
C             controller output matrix CK.
C
C     LDCK    INTEGER
C             The leading dimension of the array CK.
C             LDCK >= max(1,NCON).
C
C     DK      (input) DOUBLE PRECISION array, dimension (LDDK,NMEAS)
C             The leading NCON-by-NMEAS part of this array must contain
C             the controller input/output matrix DK.
C
C     LDDK    INTEGER
C             The leading dimension of the array DK.
C             LDDK >= max(1,NCON).
C
C     AC      (output) DOUBLE PRECISION array, dimension (LDAC,2*N)
C             The leading 2*N-by-2*N part of this array contains the
C             closed-loop system state matrix AC.
C
C     LDAC    INTEGER
C             The leading dimension of the array AC.
C             LDAC >= max(1,2*N).
C
C     BC      (output) DOUBLE PRECISION array, dimension (LDBC,M-NCON)
C             The leading 2*N-by-(M-NCON) part of this array contains
C             the closed-loop system input matrix BC.
C
C     LDBC    INTEGER
C             The leading dimension of the array BC.
C             LDBC >= max(1,2*N).
C
C     CC      (output) DOUBLE PRECISION array, dimension (LDCC,2*N)
C             The leading (NP-NMEAS)-by-2*N part of this array contains
C             the closed-loop system output matrix CC.
C
C     LDCC    INTEGER
C             The leading dimension of the array CC.
C             LDCC >= max(1,NP-NMEAS).
C
C     DC      (output) DOUBLE PRECISION array, dimension (LDDC,M-NCON)
C             The leading (NP-NMEAS)-by-(M-NCON) part of this array
C             contains the closed-loop system input/output matrix DC.
C
C     LDDC    INTEGER
C             The leading dimension of the array DC.
C             LDDC >= max(1,NP-NMEAS).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension 2*max(NCON,NMEAS)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal
C             LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 2*M*M+NP*NP+2*M*N+M*NP+2*N*NP.
C             For good performance, LDWORK must generally be larger.
C
C     Error Indicactor
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix Inp2 - D22*DK is singular to working
C                   precision;
C             = 2:  if the matrix Im2 - DK*D22 is singular to working
C                   precision.
C
C     METHOD
C
C     The routine implements the formulas given in [1].
C
C     REFERENCES
C
C     [1] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and
C         Smith, R.
C         mu-Analysis and Synthesis Toolbox.
C         The MathWorks Inc., Natick, Mass., 1995.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the result depends on the condition numbers of the
C     matrices  Inp2 - D22*DK  and  Im2 - DK*D22.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999.
C     A. Markovski, Technical University, Sofia, April, 2003.
C
C     KEYWORDS
C
C     Closed loop systems, feedback control, robust control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAC, LDAK, LDB, LDBC, LDBK, LDC,
     $                   LDCC, LDCK, LDD, LDDC, LDDK, LDWORK, M, N,
     $                   NCON, NMEAS, NP
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ),
     $                   B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ),
     $                   C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ),
     $                   DWORK( * )
C     ..
C     .. Local Scalars ..
      INTEGER            INFO2, IW2, IW3, IW4, IW5, IW6, IW7, IW8, IWRK,
     $                   LWAMAX, M1, M2, MINWRK, N2, NP1, NP2
      DOUBLE PRECISION   ANORM, EPS, RCOND
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGECON, DGEMM, DGETRF, DGETRI, DLACPY, DLASET,
     $                   XERBLA
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      N2  = 2*N
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
      ELSE IF( LDAC.LT.MAX( 1, N2 ) ) THEN
         INFO = -23
      ELSE IF( LDBC.LT.MAX( 1, N2 ) ) THEN
         INFO = -25
      ELSE IF( LDCC.LT.MAX( 1, NP1 ) ) THEN
         INFO = -27
      ELSE IF( LDDC.LT.MAX( 1, NP1 ) ) THEN
         INFO = -29
      ELSE
C
C        Compute workspace.
C
         MINWRK = 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP
         IF( LDWORK.LT.MINWRK )
     $      INFO = -32
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10LD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 .OR. NP.EQ.0 .OR. M1.EQ.0 .OR. M2.EQ.0
     $    .OR. NP1.EQ.0 .OR. NP2.EQ.0 ) THEN
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Get the machine precision.
C
      EPS = DLAMCH( 'Epsilon' )
C
C     Workspace usage.
C
      IW2 = NP2*NP2 + 1
      IW3 = IW2 + M2*M2
      IW4 = IW3 + NP2*N
      IW5 = IW4 + M2*N
      IW6 = IW5 + NP2*M1
      IW7 = IW6 + M2*M1
      IW8 = IW7 + M2*N
      IWRK = IW8 + NP2*N
C
C     Compute inv(Inp2 - D22*DK) .
C
      CALL DLASET( 'Full', NP2, NP2, ZERO, ONE, DWORK, NP2 )
      CALL DGEMM( 'N', 'N', NP2, NP2, M2, -ONE, D( NP1+1, M1+1 ),
     $            LDD, DK, LDDK, ONE, DWORK, NP2 )
      ANORM = DLANGE( '1', NP2, NP2, DWORK, NP2, DWORK( IWRK ) )
      CALL DGETRF( NP2, NP2, DWORK, NP2, IWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      CALL DGECON( '1', NP2, DWORK, NP2, ANORM, RCOND, DWORK( IWRK ),
     $             IWORK( NP2+1 ), INFO )
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND.LT.EPS ) THEN
         INFO = 1
         RETURN
      END IF
      CALL DGETRI( NP2, DWORK, NP2, IWORK, DWORK( IWRK ), LDWORK-IWRK+1,
     $             INFO2 )
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Compute inv(Im2 - DK*D22) .
C
      CALL DLASET( 'Full', M2, M2, ZERO, ONE, DWORK( IW2 ), M2 )
      CALL DGEMM( 'N', 'N', M2, M2, NP2, -ONE, DK, LDDK,
     $            D( NP1+1, M1+1 ), LDD, ONE, DWORK( IW2 ), M2 )
      ANORM = DLANGE( '1', M2, M2, DWORK( IW2 ), M2, DWORK( IWRK ) )
      CALL DGETRF( M2, M2, DWORK( IW2 ), M2, IWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
      CALL DGECON( '1', M2, DWORK( IW2 ), M2, ANORM, RCOND,
     $             DWORK( IWRK ), IWORK( M2+1 ), INFO )
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND.LT.EPS ) THEN
         INFO = 2
         RETURN
      END IF
      CALL DGETRI( M2, DWORK( IW2 ), M2, IWORK, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
C     Compute inv(Inp2 - D22*DK)*C2 .
C
      CALL DGEMM( 'N', 'N', NP2, N, NP2, ONE, DWORK, NP2, C( NP1+1, 1 ),
     $            LDC, ZERO, DWORK( IW3 ), NP2 )
C
C     Compute DK*inv(Inp2 - D22*DK)*C2 .
C
      CALL DGEMM( 'N', 'N', M2, N, NP2, ONE, DK, LDDK,  DWORK( IW3 ),
     $            NP2, ZERO, DWORK( IW4 ), M2 )
C
C     Compute inv(Inp2 - D22*DK)*D21 .
C
      CALL DGEMM( 'N', 'N', NP2, M1, NP2, ONE, DWORK, NP2,
     $            D( NP1+1, 1 ), LDD, ZERO, DWORK( IW5 ), NP2 )
C
C     Compute DK*inv(Inp2 - D22*DK)*D21 .
C
      CALL DGEMM( 'N', 'N', M2, M1, NP2, ONE, DK, LDDK, DWORK( IW5 ),
     $            NP2, ZERO, DWORK( IW6 ), M2 )
C
C     Compute inv(Im2 - DK*D22)*CK .
C
      CALL DGEMM( 'N', 'N', M2, N, M2, ONE, DWORK( IW2 ), M2, CK, LDCK,
     $            ZERO, DWORK( IW7 ), M2 )
C
C     Compute D22*inv(Im2 - DK*D22)*CK .
C
      CALL DGEMM( 'N', 'N', NP2, N, M2, ONE, D( NP1+1, M1+1 ), LDD,
     $            DWORK( IW7 ), M2, ZERO, DWORK( IW8 ), NP2 )
C
C     Compute AC .
C
      CALL DLACPY( 'Full', N, N, A, LDA, AC, LDAC )
      CALL DGEMM( 'N', 'N', N, N, M2, ONE, B( 1, M1+1 ), LDB,
     $            DWORK( IW4 ), M2, ONE, AC, LDAC )
      CALL DGEMM( 'N', 'N', N, N, M2, ONE, B( 1, M1+1 ), LDB,
     $            DWORK( IW7 ), M2, ZERO, AC( 1, N+1 ), LDAC )
      CALL DGEMM( 'N', 'N', N, N, NP2, ONE, BK, LDBK, DWORK( IW3 ), NP2,
     $            ZERO, AC( N+1, 1 ), LDAC )
      CALL DLACPY( 'Full', N, N, AK, LDAK, AC( N+1, N+1 ), LDAC )
      CALL DGEMM( 'N', 'N', N, N, NP2, ONE, BK, LDBK, DWORK( IW8 ), NP2,
     $            ONE, AC( N+1, N+1 ), LDAC )
C
C     Compute BC .
C
      CALL DLACPY( 'Full', N, M1, B, LDB, BC, LDBC )
      CALL DGEMM( 'N', 'N', N, M1, M2, ONE, B( 1, M1+1 ), LDB,
     $            DWORK( IW6 ), M2, ONE, BC, LDBC )
      CALL DGEMM( 'N', 'N', N, M1, NP2, ONE, BK, LDBK, DWORK( IW5 ),
     $            NP2, ZERO, BC( N+1, 1 ), LDBC )
C
C     Compute CC .
C
      CALL DLACPY( 'Full', NP1, N, C, LDC, CC, LDCC )
      CALL DGEMM( 'N', 'N', NP1, N, M2, ONE, D( 1, M1+1 ), LDD,
     $            DWORK( IW4 ), M2, ONE, CC, LDCC )
      CALL DGEMM( 'N', 'N', NP1, N, M2, ONE, D( 1, M1+1 ), LDD,
     $            DWORK( IW7 ), M2, ZERO, CC( 1, N+1 ), LDCC )
C
C     Compute DC .
C
      CALL DLACPY( 'Full', NP1, M1, D, LDD, DC, LDDC )
      CALL DGEMM( 'N', 'N', NP1, M1, M2, ONE, D( 1, M1+1 ), LDD,
     $            DWORK( IW6 ), M2, ONE, DC, LDDC )
C
      RETURN
C *** Last line of SB10LD ***
      END
