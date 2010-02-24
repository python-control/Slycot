      SUBROUTINE AB09HY( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   SCALEC, SCALEO, S, LDS, R, LDR, IWORK,
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
C     To compute the Cholesky factors Su and Ru of the controllability
C     Grammian P = Su*Su' and observability Grammian Q = Ru'*Ru,
C     respectively, satisfying
C
C            A*P  + P*A' +  scalec^2*B*B'   = 0,       (1)
C
C            A'*Q + Q*A  +  scaleo^2*Cw'*Cw = 0,       (2)
C
C     where
C            Cw = Hw - Bw'*X,
C            Hw = inv(Dw)*C,
C            Bw = (B*D' + P*C')*inv(Dw'),
C            D*D' = Dw*Dw' (Dw upper triangular),
C
C     and, with Aw = A - Bw*Hw, X is the stabilizing solution of the
C     Riccati equation
C
C            Aw'*X + X*Aw + Hw'*Hw + X*Bw*Bw'*X = 0.   (3)
C
C     The P-by-M matrix D must have full row rank. Matrix A must be
C     stable and in a real Schur form.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of state-space representation, i.e.,
C             the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  M >= P >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             stable state dynamics matrix A in a real Schur canonical
C             form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B, corresponding to the Schur matrix A.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C, corresponding to the Schur
C             matrix A.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must
C             contain the full row rank input/output matrix D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     SCALEC  (output) DOUBLE PRECISION
C             Scaling factor for the controllability Grammian in (1).
C
C     SCALEO  (output) DOUBLE PRECISION
C             Scaling factor for the observability Grammian in (2).
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor Su of the cotrollability
C             Grammian P = Su*Su' satisfying (1).
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,N).
C
C     R       (output) DOUBLE PRECISION array, dimension (LDR,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor Ru of the observability
C             Grammian Q = Ru'*Ru satisfying (2).
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension 2*N
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) contains RCOND, the reciprocal
C             condition number of the U11 matrix from the expression
C             used to compute X = U21*inv(U11). A small value RCOND
C             indicates possible ill-conditioning of the Riccati
C             equation (3).
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 2, N*(MAX(N,M,P)+5),
C                            2*N*P+MAX(P*(M+2),10*N*(N+1) ) ).
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension 2*N
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the state matrix A is not stable or is not in a
C                   real Schur form;
C             = 2:  the reduction of Hamiltonian matrix to real Schur
C                   form failed;
C             = 3:  the reordering of the real Schur form of the
C                   Hamiltonian matrix failed;
C             = 4:  the Hamiltonian matrix has less than N stable
C                   eigenvalues;
C             = 5:  the coefficient matrix U11 in the linear system
C                   X*U11 = U21, used to determine X, is singular to
C                   working precision;
C             = 6:  the feedthrough matrix D has not a full row rank P.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000.
C     D. Sima, University of Bucharest, May 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2000.
C     Based on the RASP routines SRGRO and SRGRO1, by A. Varga, 1992.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001.
C
C     KEYWORDS
C
C     Minimal realization, model reduction, multivariable system,
C     state-space model, state-space representation,
C     stochastic balancing.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER          INFO, LDA, LDB, LDC, LDD, LDR, LDS, LDWORK, M, N,
     $                 P
      DOUBLE PRECISION SCALEC, SCALEO
C     .. Array Arguments ..
      INTEGER          IWORK(*)
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                 DWORK(*), R(LDR,*), S(LDS,*)
      LOGICAL          BWORK(*)
C     .. Local Scalars ..
      INTEGER          I, IERR, KBW, KCW, KD, KDW, KG, KQ, KS, KTAU, KU,
     $                 KW, KWI, KWR, LW, N2, WRKOPT
      DOUBLE PRECISION RCOND, RTOL
C     .. External Functions ..
      DOUBLE PRECISION DLANGE, DLAMCH
      EXTERNAL         DLANGE, DLAMCH
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DGERQF, DLACPY, DORGRQ, DSYRK, DTRMM,
     $                 DTRSM, SB02MD, SB03OU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        ABS, DBLE, INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LW   = MAX( 2, N*( MAX( N, M, P ) + 5 ),
     $            2*N*P + MAX( P*(M + 2), 10*N*(N + 1) ) )
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDS.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -20
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09HY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      SCALEC = ONE
      SCALEO = ONE
      IF( MIN( N, M, P ).EQ.0 ) THEN
         DWORK(1) = TWO
         DWORK(2) = ONE
         RETURN
      END IF
C
C     Solve for Su the Lyapunov equation
C                                      2
C     A*(Su*Su') + (Su*Su')*A' + scalec *B*B' = 0 .
C
C     Workspace:  need   N*(MAX(N,M) + 5);
C                 prefer larger.
C
      KU   = 1
      KTAU = KU + N*MAX( N, M )
      KW   = KTAU + N
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), N )
      CALL SB03OU( .FALSE., .TRUE., N, M, A, LDA, DWORK(KU), N,
     $             DWORK(KTAU), S, LDS, SCALEC, DWORK(KW),
     $             LDWORK - KW + 1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 1
         RETURN
      ENDIF
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Allocate workspace for Bw' (P*N), Cw (P*N), Q2 (P*M),
C     where Q2 = inv(Dw)*D.
C     Workspace:  need   2*N*P + P*M.
C
      KBW  = 1
      KCW  = KBW  + P*N
      KD   = KCW  + P*N
      KDW  = KD   + P*(M - P)
      KTAU = KD   + P*M
      KW   = KTAU + P
C
C     Compute an upper-triangular Dw such that D*D' = Dw*Dw', using
C     the RQ-decomposition of D: D = [0 Dw]*( Q1 ).
C                                           ( Q2 )
C     Additional workspace:  need 2*P; prefer P + P*NB.
C
      CALL DLACPY( 'F', P, M, D, LDD, DWORK(KD), P )
      CALL DGERQF( P, M, DWORK(KD), P, DWORK(KTAU), DWORK(KW),
     $             LDWORK-KW+1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Check the full row rank of D.
C
      RTOL = DBLE( M ) * DLAMCH( 'E' ) *
     $       DLANGE( '1', P, M, D, LDD, DWORK )
      DO 10 I = KDW, KDW+P*P-1, P+1
         IF( ABS( DWORK(I) ).LE.RTOL ) THEN
            INFO = 6
            RETURN
         END IF
   10 CONTINUE
C                    -1
C     Compute Hw = Dw  *C.
C
      CALL DLACPY( 'F', P, N, C, LDC, DWORK(KCW), P )
      CALL DTRSM( 'Left', 'Upper', 'No-transpose', 'Non-unit', P, N,
     $            ONE, DWORK(KDW), P, DWORK(KCW), P )
C
C     Compute Bw' = inv(Dw)*(D*B' + C*Su*Su').
C
C     Compute first Hw*Su*Su' in Bw'.
C
      CALL DLACPY( 'F', P, N, DWORK(KCW), P, DWORK(KBW), P )
      CALL DTRMM( 'Right', 'Upper', 'No-transpose', 'Non-unit', P, N,
     $            ONE, S, LDS, DWORK(KBW), P )
      CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', P, N,
     $            ONE, S, LDS, DWORK(KBW), P )
C
C     Compute Q2 = inv(Dw)*D, as the last P lines of the orthogonal
C     matrix ( Q1 ) from the RQ decomposition of D.
C            ( Q2 )
C     Additional workspace:  need P; prefer P*NB.
C
      CALL DORGRQ( P, M, P, DWORK(KD), P, DWORK(KTAU), DWORK(KW),
     $             LDWORK-KW+1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Compute Bw' <- Bw' + Q2*B'.
C
      CALL DGEMM( 'No-transpose', 'Transpose', P, N, M, ONE,
     $            DWORK(KD), P, B, LDB, ONE, DWORK(KBW), P )
C
C     Compute Aw = A - Bw*Hw in R.
C
      CALL DLACPY( 'F', N, N, A, LDA, R, LDR )
      CALL DGEMM( 'Transpose', 'No-transpose', N, N, P, -ONE,
     $            DWORK(KBW), P, DWORK(KCW), P, ONE, R, LDR )
C
C     Allocate storage to solve the Riccati equation (3) for
C     G(N*N), Q(N*N), WR(2N), WI(2N), S(2N*2N), U(2N*2N).
C
      N2  = N + N
      KG  = KD
      KQ  = KG  + N*N
      KWR = KQ  + N*N
      KWI = KWR + N2
      KS  = KWI + N2
      KU  = KS  + N2*N2
      KW  = KU  + N2*N2
C
C     Compute G = -Bw*Bw'.
C
      CALL DSYRK( 'Upper', 'Transpose', N, P, -ONE, DWORK(KBW), P, ZERO,
     $            DWORK(KG), N )
C
C     Compute Q = Hw'*Hw.
C
      CALL DSYRK( 'Upper', 'Transpose', N, P, ONE, DWORK(KCW), P, ZERO,
     $            DWORK(KQ), N )
C
C     Solve
C
C        Aw'*X + X*Aw + Q - X*G*X = 0,
C
C     with Q =  Hw'*Hw  and  G = -Bw*Bw'.
C     Additional workspace: need   6*N;
C                           prefer larger.
C
      CALL SB02MD( 'Continuous', 'None', 'Upper', 'General', 'Stable',
     $             N, R, LDR, DWORK(KG), N, DWORK(KQ), N, RCOND,
     $             DWORK(KWR), DWORK(KWI), DWORK(KS), N2,
     $             DWORK(KU), N2, IWORK, DWORK(KW), LDWORK-KW+1,
     $             BWORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Compute Cw = Hw - Bw'*X.
C
      CALL DGEMM ( 'No-transpose', 'No-transpose', P, N, N, -ONE,
     $              DWORK(KBW), P, DWORK(KQ), N, ONE, DWORK(KCW), P )
C
C     Solve for Ru the Lyapunov equation
C                                      2
C     A'*(Ru'*Ru) + (Ru'*Ru)*A + scaleo  * Cw'*Cw = 0 .
C
C     Workspace:  need   N*(MAX(N,P) + 5);
C                 prefer larger.
C
      KTAU = KCW  + N*MAX( N, P )
      KW   = KTAU + N
C
      CALL SB03OU( .FALSE., .FALSE., N, P, A, LDA, DWORK(KCW), P,
     $             DWORK(KTAU), R, LDR, SCALEO, DWORK(KW),
     $             LDWORK - KW + 1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Save optimal workspace and RCOND.
C
      DWORK(1) = WRKOPT
      DWORK(2) = RCOND
C
      RETURN
C *** Last line of AB09HY ***
      END
