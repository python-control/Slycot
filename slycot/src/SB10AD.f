      SUBROUTINE SB10AD( JOB, N, M, NP, NCON, NMEAS, GAMMA, A, LDA,
     $                   B, LDB, C, LDC, D, LDD, AK, LDAK, BK, LDBK, CK,
     $                   LDCK, DK, LDDK, AC, LDAC, BC, LDBC, CC, LDCC,
     $                   DC, LDDC, RCOND, GTOL, ACTOL, IWORK, LIWORK,
     $                   DWORK, LDWORK, BWORK, LBWORK, INFO )
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
C     To compute the matrices of an H-infinity optimal n-state
C     controller
C
C              | AK | BK |
C          K = |----|----|,
C              | CK | DK |
C
C     using modified Glover's and Doyle's 1988 formulas, for the system
C
C              | A  | B1  B2  |   | A | B |
C          P = |----|---------| = |---|---|
C              | C1 | D11 D12 |   | C | D |
C              | C2 | D21 D22 |
C
C     and for the estimated minimal possible value of gamma with respect
C     to GTOL, where B2 has as column size the number of control inputs
C     (NCON) and C2 has as row size the number of measurements (NMEAS)
C     being provided to the controller, and then to compute the matrices
C     of the closed-loop system
C
C              | AC | BC |
C          G = |----|----|,
C              | CC | DC |
C
C     if the stabilizing controller exists.
C
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) D12 is full column rank and D21 is full row rank,
C
C     (A3) | A-j*omega*I  B2  | has full column rank for all omega,
C          |    C1        D12 |
C
C     (A4) | A-j*omega*I  B1  |  has full row rank for all omega.
C          |    C2        D21 |
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     JOB     (input) INTEGER
C             Indicates the strategy for reducing the GAMMA value, as
C             follows:
C             = 1: Use bisection method for decreasing GAMMA from GAMMA
C                  to GAMMAMIN until the closed-loop system leaves
C                  stability.
C             = 2: Scan from GAMMA to 0 trying to find the minimal GAMMA
C                  for which the closed-loop system retains stability.
C             = 3: First bisection, then scanning.
C             = 4: Find suboptimal controller only.
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
C     GAMMA   (input/output) DOUBLE PRECISION
C             The initial value of gamma on input. It is assumed that
C             gamma is sufficiently large so that the controller is
C             admissible. GAMMA >= 0.
C             On output it contains the minimal estimated gamma.
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
C     RCOND   (output) DOUBLE PRECISION array, dimension (4)
C                      For the last successful step:
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
C     GTOL    DOUBLE PRECISION
C             Tolerance used for controlling the accuracy of GAMMA
C             and its distance to the estimated minimal possible
C             value of GAMMA.
C             If GTOL <= 0, then a default value equal to sqrt(EPS)
C             is used, where EPS is the relative machine precision.
C
C     ACTOL   DOUBLE PRECISION
C             Upper bound for the poles of the closed-loop system
C             used for determining if it is stable.
C             ACTOL <= 0 for stable systems.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C
C     LIWORK  INTEGER
C             The dimension of the array IWORK.
C             LIWORK >= max(2*max(N,M-NCON,NP-NMEAS,NCON,NMEAS),N*N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= LW1 + max(1,LW2,LW3,LW4,LW5 + MAX(LW6,LW7)),
C             where
C             LW1 = N*M + NP*N + NP*M + M2*M2 + NP2*NP2;
C             LW2 = max( ( N + NP1 + 1 )*( N + M2 ) +
C                          max( 3*( N + M2 ) + N + NP1, 5*( N + M2 ) ),
C                        ( N + NP2 )*( N + M1 + 1 ) +
C                          max( 3*( N + NP2 ) + N + M1, 5*( N + NP2 ) ),
C                        M2 + NP1*NP1 + max( NP1*max( N, M1 ),
C                                            3*M2 + NP1, 5*M2 ),
C                        NP2 + M1*M1 +  max( max( N, NP1 )*M1,
C                                            3*NP2 + M1, 5*NP2 ) );
C             LW3 = max( ND1*M1 + max( 4*min( ND1, M1 ) + max( ND1,M1 ),
C                                      6*min( ND1, M1 ) ),
C                        NP1*ND2 + max( 4*min( NP1, ND2 ) +
C                                                        max( NP1,ND2 ),
C                                       6*min( NP1, ND2 ) ) );
C             LW4 = 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP;
C             LW5 = 2*N*N + M*N + N*NP;
C             LW6 = max( M*M   + max( 2*M1, 3*N*N +
C                                     max( N*M, 10*N*N + 12*N + 5 ) ),
C                        NP*NP + max( 2*NP1, 3*N*N +
C                                     max( N*NP, 10*N*N + 12*N + 5 ) ));
C             LW7 = M2*NP2 + NP2*NP2 + M2*M2 +
C                   max( ND1*ND1 + max( 2*ND1, ( ND1 + ND2 )*NP2 ),
C                        ND2*ND2 + max( 2*ND2, ND2*M2 ), 3*N,
C                        N*( 2*NP2 + M2 ) +
C                        max( 2*N*M2, M2*NP2 +
C                                     max( M2*M2 + 3*M2, NP2*( 2*NP2 +
C                                          M2 + max( NP2, N ) ) ) ) );
C             M1  = M   - M2, NP1 = NP - NP2,
C             ND1 = NP1 - M2, ND2 = M1 - NP2.
C             For good performance, LDWORK must generally be larger.
C
C     BWORK   LOGICAL array, dimension (LBWORK)
C
C     LBWORK  INTEGER
C             The dimension of the array BWORK.  LBWORK >= 2*N.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix | A-j*omega*I  B2  | had not full
C                                 |    C1        D12 |
C                   column rank in respect to the tolerance EPS;
C             = 2:  if the matrix | A-j*omega*I  B1  |  had not full row
C                                 |    C2        D21 |
C                   rank in respect to the tolerance EPS;
C             = 3:  if the matrix D12 had not full column rank in
C                   respect to the tolerance SQRT(EPS);
C             = 4:  if the matrix D21 had not full row rank in respect
C                   to the tolerance SQRT(EPS);
C             = 5:  if the singular value decomposition (SVD) algorithm
C                   did not converge (when computing the SVD of one of
C                   the matrices |A   B2 |, |A   B1 |, D12 or D21);
C                                |C1  D12|  |C2  D21|
C             = 6:  if the controller is not admissible (too small value
C                   of gamma);
C             = 7:  if the X-Riccati equation was not solved
C                   successfully (the controller is not admissible or
C                   there are numerical difficulties);
C             = 8:  if the Y-Riccati equation was not solved
C                   successfully (the controller is not admissible or
C                   there are numerical difficulties);
C             = 9:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is
C                   zero [3];
C             = 10: if there are numerical problems when estimating
C                   singular values of D1111, D1112, D1111', D1121';
C             = 11: if the matrices Inp2 - D22*DK or Im2 - DK*D22
C                   are singular to working precision;
C             = 12: if a stabilizing controller cannot be found.
C
C     METHOD
C
C     The routine implements the Glover's and Doyle's 1988 formulas [1],
C     [2], modified to improve the efficiency as described in [3].
C
C     JOB = 1: It tries with a decreasing value of GAMMA, starting with
C     the given, and with the newly obtained controller estimates of the
C     closed-loop system. If it is stable, (i.e., max(eig(AC)) < ACTOL)
C     the iterations can be continued until the given tolerance between
C     GAMMA and the estimated GAMMAMIN is reached. Otherwise, in the
C     next step GAMMA is increased. The step in the all next iterations
C     is step = step/2. The closed-loop system is obtained by the
C     formulas given in [2].
C
C     JOB = 2: The same as for JOB = 1, but with non-varying step till
C     GAMMA = 0, step = max(0.1, GTOL).
C
C     JOB = 3: Combines the JOB = 1 and JOB = 2 cases for a quicker
C     procedure.
C
C     JOB = 4: Suboptimal controller for current GAMMA only.
C
C     REFERENCES
C
C     [1] Glover, K. and Doyle, J.C.
C         State-space formulae for all stabilizing controllers that
C         satisfy an Hinf norm bound and relations to risk sensitivity.
C         Systems and Control Letters, vol. 11, pp. 167-172, 1988.
C
C     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and
C         Smith, R.
C         mu-Analysis and Synthesis Toolbox.
C         The MathWorks Inc., Natick, MA, 1995.
C
C     [3] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M.
C         Fortran 77 routines for Hinf and H2 design of continuous-time
C         linear control systems.
C         Rep. 98-14, Department of Engineering, Leicester University,
C         Leicester, U.K., 1998.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the result depends on the condition numbers of the
C     input and output transformations and on the condition numbers of
C     the two Riccati equations, as given by the values of RCOND(1),
C     RCOND(2), RCOND(3) and RCOND(4), respectively.
C     This approach by estimating the closed-loop system and checking
C     its poles seems to be reliable.
C
C     CONTRIBUTORS
C
C     A. Markovski, P.Hr. Petkov, D.W. Gu and M.M. Konstantinov,
C     July 2003.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, H-infinity optimal control, robust
C     control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, P1, THOUS
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     P1 = 0.1D+0, THOUS = 1.0D+3 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, JOB, LBWORK, LDA, LDAC, LDAK, LDB, LDBC,
     $                   LDBK, LDC, LDCC, LDCK, LDD, LDDC, LDDK, LDWORK,
     $                   LIWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   ACTOL, GAMMA, GTOL
C     ..
C     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ),
     $                   B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ),
     $                   C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ),
     $                   DWORK( * ), RCOND( 4 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, INF, INFO2, INFO3, IWAC, IWC, IWD, IWD1,
     $                   IWF, IWH, IWRE, IWRK, IWS1, IWS2, IWTU, IWTY,
     $                   IWWI, IWWR, IWX, IWY, LW1, LW2, LW3, LW4, LW5,
     $                   LW6, LW7, LWAMAX, M1, M11, M2, MINWRK, MODE,
     $                   NP1, NP11, NP2
      DOUBLE PRECISION   GAMABS, GAMAMN, GAMAMX, GTOLL, MINEAC, STEPG,
     $                   TOL2
C     ..
C     .. External Functions ..
      LOGICAL            SELECT
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEES, DGESVD, DLACPY, SB10LD, SB10PD, SB10QD,
     $                   SB10RD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode and test input parameters.
C
      M1   = M - NCON
      M2   = NCON
      NP1  = NP - NMEAS
      NP2  = NMEAS
      NP11 = NP1 - M2
      M11  = M1 - NP2
C
      INFO = 0
      IF ( JOB.LT.1 .OR. JOB.GT.4 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( NP.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCON.LT.0 .OR. M1.LT.0 .OR. M2.GT.NP1 ) THEN
         INFO = -5
      ELSE IF( NMEAS.LT.0 .OR. NP1.LT.0 .OR. NP2.GT.M1 ) THEN
         INFO = -6
      ELSE IF( GAMMA.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDC.LT.MAX( 1, NP ) ) THEN
         INFO = -13
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -15
      ELSE IF( LDAK.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDBK.LT.MAX( 1, N ) ) THEN
         INFO = -19
      ELSE IF( LDCK.LT.MAX( 1, M2 ) ) THEN
         INFO = -21
      ELSE IF( LDDK.LT.MAX( 1, M2 ) ) THEN
         INFO = -23
      ELSE IF( LDAC.LT.MAX( 1, 2*N ) ) THEN
         INFO = -25
      ELSE IF( LDBC.LT.MAX( 1, 2*N ) ) THEN
         INFO = -27
      ELSE IF( LDCC.LT.MAX( 1, NP1 ) ) THEN
         INFO = -29
      ELSE IF( LDDC.LT.MAX( 1, NP1 ) ) THEN
         INFO = -31
      ELSE
C
C        Compute workspace.
C
         LW1 = N*M + NP*N + NP*M + M2*M2 + NP2*NP2
         LW2 = MAX( ( N + NP1 + 1 )*( N + M2 ) +
     $                MAX( 3*( N + M2 ) + N + NP1, 5*( N + M2 ) ),
     $              ( N + NP2 )*( N + M1 + 1 ) +
     $                MAX( 3*( N + NP2 ) + N + M1, 5*( N + NP2 ) ),
     $              M2 + NP1*NP1 + MAX( NP1*MAX( N, M1 ), 3*M2 + NP1,
     $                                  5*M2 ),
     $              NP2 + M1*M1 +  MAX( MAX( N, NP1 )*M1, 3*NP2 + M1,
     $                                  5*NP2 ) )
         LW3 = MAX( NP11*M1 + MAX( 4*MIN( NP11, M1 ) + MAX( NP11, M1 ),
     $                             6*MIN( NP11, M1 ) ),
     $              NP1*M11 + MAX( 4*MIN( NP1, M11 ) + MAX( NP1, M11 ),
     $                             6*MIN( NP1, M11 ) ) )
         LW4 = 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP
         LW5 = 2*N*N + M*N + N*NP
         LW6 = MAX( M*M   + MAX( 2*M1, 3*N*N +
     $                           MAX( N*M, 10*N*N + 12*N + 5 ) ),
     $              NP*NP + MAX( 2*NP1, 3*N*N +
     $                           MAX( N*NP, 10*N*N + 12*N + 5 ) ) )
         LW7 = M2*NP2 + NP2*NP2 + M2*M2 +
     $         MAX( NP11*NP11 + MAX( 2*NP11, ( NP11 + M11 )*NP2 ),
     $              M11*M11 + MAX( 2*M11, M11*M2 ), 3*N,
     $              N*( 2*NP2 + M2 ) +
     $              MAX( 2*N*M2, M2*NP2 +
     $                           MAX( M2*M2 + 3*M2, NP2*( 2*NP2 +
     $                                M2 + MAX( NP2, N ) ) ) ) )
         MINWRK = LW1 + MAX( 1, LW2, LW3, LW4, LW5 + MAX( LW6, LW7 ) )
         IF( LDWORK.LT.MINWRK ) THEN
            INFO = -38
         ELSE IF( LIWORK.LT.MAX( 2*MAX( N, M1, NP1, M2, NP2 ),
     $                           N*N ) ) THEN
            INFO = -36
         ELSE IF( LBWORK.LT.2*N ) THEN
            INFO = -40
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB10AD', -INFO )
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
      MODE = JOB
      IF ( MODE.GT.2 )
     $   MODE = 1
      GTOLL = GTOL
      IF( GTOLL.LE.ZERO ) THEN
C
C        Set the default value of the tolerance for GAMMA.
C
         GTOLL = SQRT( DLAMCH( 'Epsilon' ) )
      END IF
C
C     Workspace usage 1.
C
      IWC  = 1 + N*M
      IWD  = IWC + NP*N
      IWTU = IWD + NP*M
      IWTY = IWTU + M2*M2
      IWRK = IWTY + NP2*NP2
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
C
      CALL DLACPY( 'Full', NP, N, C, LDC, DWORK( IWC ), NP )
C
      CALL DLACPY( 'Full', NP, M, D, LDD, DWORK( IWD ), NP )
C
C     Transform the system so that D12 and D21 satisfy the formulas
C     in the computation of the Hinf optimal controller.
C     Workspace:  need   LW1 + MAX(1,LWP1,LWP2,LWP3,LWP4),
C                 prefer larger,
C                 where
C             LW1  = N*M + NP*N + NP*M + M2*M2 + NP2*NP2
C             LWP1 = (N+NP1+1)*(N+M2) + MAX(3*(N+M2)+N+NP1,5*(N+M2)),
C             LWP2 = (N+NP2)*(N+M1+1) + MAX(3*(N+NP2)+N+M1,5*(N+NP2)),
C             LWP3 = M2 + NP1*NP1 + MAX(NP1*MAX(N,M1),3*M2+NP1,5*M2),
C             LWP4 = NP2 + M1*M1 + MAX(MAX(N,NP1)*M1,3*NP2+M1,5*NP2),
C             with M1 = M - M2 and NP1 = NP - NP2.
C             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is
C             LW1 + MAX(1,(N+Q)*(N+Q+6),Q*(Q+MAX(N,Q,5)+1).
C
      TOL2 = -ONE
C
      CALL SB10PD( N, M, NP, NCON, NMEAS, A, LDA, DWORK, N,
     $             DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWTU ),
     $             M2, DWORK( IWTY ), NP2, RCOND, TOL2, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
C
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
      IF ( INFO2.NE.0 ) THEN
         INFO = INFO2
         RETURN
      END IF
C
C     Workspace usage 2.
C
      IWD1 = IWRK
      IWS1 = IWD1 + NP11*M1
C
C     Check if GAMMA < max(sigma[D1111,D1112],sigma[D1111',D1121']).
C     Workspace:  need   LW1 + MAX(1, LWS1, LWS2),
C                 prefer larger,
C                 where
C     LWS1 = NP11*M1 + MAX(4*MIN(NP11,M1)+MAX(NP11,M1),6*MIN(NP11,M1))
C     LWS2 = NP1*M11 + MAX(4*MIN(NP1,M11)+MAX(NP1,M11),6*MIN(NP1,M11))
C
      INFO2 = 0
      INFO3 = 0
C
      IF ( NP11.NE.0 .AND. M1.NE.0 ) THEN
         IWRK = IWS1 + MIN( NP11, M1 )
         CALL DLACPY( 'Full', NP11, M1, DWORK(IWD), LDD, DWORK(IWD1),
     $                NP11 )
         CALL DGESVD( 'N', 'N', NP11, M1, DWORK(IWD1), NP11,
     $                DWORK(IWS1), DWORK(IWS1), 1, DWORK(IWS1), 1,
     $                DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
         LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
      ELSE
         DWORK(IWS1) = ZERO
      END IF
C
      IWS2 = IWD1 + NP1*M11
      IF ( NP1.NE.0 .AND. M11.NE.0 ) THEN
         IWRK = IWS2 + MIN( NP1, M11 )
         CALL DLACPY( 'Full', NP1, M11, DWORK(IWD), LDD, DWORK(IWD1),
     $                NP1 )
         CALL DGESVD( 'N', 'N', NP1, M11, DWORK(IWD1), NP1, DWORK(IWS2),
     $                DWORK(IWS2), 1, DWORK(IWS2), 1, DWORK( IWRK ),
     $                LDWORK-IWRK+1, INFO3 )
         LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
      ELSE
         DWORK(IWS2) = ZERO
      END IF
C
      GAMAMN = MAX( DWORK(IWS1), DWORK(IWS2) )
C
      IF ( INFO2.GT.0 .OR. INFO3.GT.0 ) THEN
         INFO = 10
         RETURN
      ELSE IF ( GAMMA.LE.GAMAMN ) THEN
         INFO = 6
         RETURN
      END IF
C
C     Workspace usage 3.
C
      IWX  = IWD1
      IWY  = IWX + N*N
      IWF  = IWY + N*N
      IWH  = IWF + M*N
      IWRK = IWH + N*NP
      IWAC = IWD1
      IWWR = IWAC + 4*N*N
      IWWI = IWWR + 2*N
      IWRE = IWWI + 2*N
C
C     Prepare some auxiliary variables for the gamma iteration.
C
      STEPG  = GAMMA - GAMAMN
      GAMABS = GAMMA
      GAMAMX = GAMMA
      INF = 0
C
C     ###############################################################
C
C     Begin the gamma iteration.
C
   10 CONTINUE
         STEPG = STEPG/TWO
C
C        Try to compute the state feedback and output injection
C        matrices for the current GAMMA.
C
         CALL SB10QD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, DWORK, N,
     $                DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWF ),
     $                M, DWORK( IWH ), N, DWORK( IWX ), N, DWORK( IWY ),
     $                N, RCOND(3), IWORK, DWORK( IWRK ), LDWORK-IWRK+1,
     $                BWORK, INFO2 )
C
         IF ( INFO2.NE.0 ) GOTO 30
C
C        Try to compute the Hinf suboptimal (yet) controller.
C
         CALL SB10RD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, DWORK, N,
     $                DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWF ),
     $                M, DWORK( IWH ), N, DWORK( IWTU ), M2,
     $                DWORK( IWTY ), NP2, DWORK( IWX ), N, DWORK( IWY ),
     $                N, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, IWORK,
     $                DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
C
         IF ( INFO2.NE.0 ) GOTO 30
C
C        Compute the closed-loop system.
C        Workspace: need   LW1 + 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP;
C                   prefer larger.
C
         CALL SB10LD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, D,
     $                LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, AC,
     $                LDAC, BC, LDBC, CC, LDCC, DC, LDDC, IWORK,
     $                DWORK( IWD1 ), LDWORK-IWD1+1, INFO2 )
C
         IF ( INFO2.NE.0 ) GOTO 30
C
         LWAMAX = MAX( LWAMAX, INT( DWORK( IWD1 ) ) + IWD1 - 1 )
C
C        Compute the poles of the closed-loop system.
C        Workspace:  need   LW1 + 4*N*N + 4*N + max(1,6*N);
C                    prefer larger.
C
         CALL DLACPY( 'Full', 2*N, 2*N, AC, LDAC, DWORK(IWAC), 2*N )
C
         CALL DGEES( 'N', 'N', SELECT, 2*N, DWORK(IWAC), 2*N, IWORK,
     $               DWORK(IWWR), DWORK(IWWI), DWORK(IWRE), 1,
     $               DWORK(IWRE), LDWORK-IWRE+1, BWORK, INFO2 )
C
         LWAMAX = MAX( LWAMAX, INT( DWORK( IWRE ) ) + IWRE - 1 )
C
C        Now DWORK(IWWR+I)=Re(Lambda), DWORK(IWWI+I)=Im(Lambda),
C        for I=0,2*N-1.
C
         MINEAC = -THOUS
C
         DO 20 I = 0, 2*N - 1
            MINEAC = MAX( MINEAC, DWORK(IWWR+I) )
   20    CONTINUE
C
C        Check if the closed-loop system is stable.
C
   30    IF ( MODE.EQ.1 ) THEN
            IF ( INFO2.EQ.0 .AND. MINEAC.LT.ACTOL ) THEN
               GAMABS = GAMMA
               GAMMA  = GAMMA - STEPG
               INF = 1
            ELSE
               GAMMA = MIN( GAMMA + STEPG, GAMAMX )
            END IF
         ELSE IF ( MODE.EQ.2 ) THEN
            IF ( INFO2.EQ.0 .AND. MINEAC.LT.ACTOL ) THEN
               GAMABS = GAMMA
               INF = 1
            END IF
            GAMMA = GAMMA - MAX( P1, GTOLL )
         END IF
C
C        More iterations?
C
         IF ( MODE.EQ.1 .AND. JOB.EQ.3 .AND. TWO*STEPG.LT.GTOLL ) THEN
            MODE  = 2
            GAMMA = GAMABS
         END IF
C
         IF ( JOB.NE.4 .AND.
     $        ( MODE.EQ.1 .AND. TWO*STEPG.GE.GTOLL .OR.
     $          MODE.EQ.2 .AND. GAMMA.GT.ZERO ) ) THEN
            GOTO 10
         END IF
C
C     ###############################################################
C
C     End of the gamma iteration - Return if no stabilizing controller
C     was found.
C
      IF ( INF.EQ.0 ) THEN
         INFO = 12
         RETURN
      END IF
C
C     Now compute the state feedback and output injection matrices
C     using GAMABS.
C
      GAMMA = GAMABS
C
C     Integer workspace:  need   max(2*max(N,M-NCON,NP-NMEAS),N*N).
C     Workspace: need   LW1P +
C                       max(1,M*M + max(2*M1,3*N*N +
C                                       max(N*M,10*N*N+12*N+5)),
C                           NP*NP + max(2*NP1,3*N*N +
C                                       max(N*NP,10*N*N+12*N+5)));
C                prefer larger,
C             where LW1P = LW1 + 2*N*N + M*N + N*NP.
C             An upper bound of the second term after LW1P is
C             max(1,4*Q*Q+max(2*Q,3*N*N + max(2*N*Q,10*N*N+12*N+5))).
C
      CALL SB10QD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, DWORK, N,
     $             DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWF ),
     $             M, DWORK( IWH ), N, DWORK( IWX ), N, DWORK( IWY ),
     $             N, RCOND(3), IWORK, DWORK( IWRK ), LDWORK-IWRK+1,
     $             BWORK, INFO2 )
C
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
      IF ( INFO2.GT.0 ) THEN
         INFO = INFO2 + 5
         RETURN
      END IF
C
C     Compute the Hinf optimal controller.
C     Integer workspace:  need   max(2*(max(NP,M)-M2-NP2,M2,N),NP2).
C     Workspace: need   LW1P +
C                       max(1, M2*NP2 + NP2*NP2 + M2*M2 +
C                           max(D1*D1 + max(2*D1, (D1+D2)*NP2),
C                               D2*D2 + max(2*D2, D2*M2), 3*N,
C                               N*(2*NP2 + M2) +
C                               max(2*N*M2, M2*NP2 +
C                                           max(M2*M2+3*M2, NP2*(2*NP2+
C                                                  M2+max(NP2,N))))))
C                       where D1 = NP1 - M2 = NP11, D2 = M1 - NP2 = M11;
C                prefer larger.
C             An upper bound of the second term after LW1P is
C             max( 1, Q*(3*Q + 3*N + max(2*N, 4*Q + max(Q, N)))).
C
      CALL SB10RD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, DWORK, N,
     $             DWORK( IWC ), NP, DWORK( IWD ), NP, DWORK( IWF ),
     $             M, DWORK( IWH ), N, DWORK( IWTU ), M2, DWORK( IWTY ),
     $             NP2, DWORK( IWX ), N, DWORK( IWY ), N, AK, LDAK, BK,
     $             LDBK, CK, LDCK, DK, LDDK, IWORK, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
C
      LWAMAX = MAX( LWAMAX, INT( DWORK( IWRK ) ) + IWRK - 1 )
C
      IF( INFO2.EQ.1 ) THEN
         INFO = 6
         RETURN
      ELSE IF( INFO2.EQ.2 ) THEN
         INFO = 9
         RETURN
      END IF
C
C     Integer workspace:  need   2*max(NCON,NMEAS).
C     Workspace: need   2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP;
C                prefer larger.
C
      CALL SB10LD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, D,
     $             LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, AC,
     $             LDAC, BC, LDBC, CC, LDCC, DC, LDDC, IWORK, DWORK,
     $             LDWORK, INFO2 )
C
      IF( INFO2.GT.0 ) THEN
         INFO = 11
         RETURN
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10AD ***
      END
