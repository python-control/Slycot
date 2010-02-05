      SUBROUTINE SB16CY( DICO, JOBCF, N, M, P, A, LDA, B, LDB, C, LDC,
     $                   F, LDF, G, LDG, SCALEC, SCALEO, S, LDS, R, LDR,
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
C     To compute, for a given open-loop model (A,B,C,0), and for
C     given state feedback gain F and full observer gain G,
C     such that A+B*F and A+G*C are stable, the Cholesky factors
C     Su and Ru of a controllability Grammian P = Su*Su' and of
C     an observability Grammian Q = Ru'*Ru corresponding to a
C     frequency-weighted model reduction of the left or right coprime
C     factors of the state-feedback controller.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the open-loop system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOBCF   CHARACTER*1
C             Specifies whether a left or right coprime factorization
C             of the state-feedback controller is to be used as follows:
C             = 'L':  use a left coprime factorization;
C             = 'R':  use a right coprime factorization.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the open-loop state-space representation,
C             i.e., the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the open-loop system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B of the open-loop system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C of the open-loop system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     F       (input) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading M-by-N part of this array must contain a
C             stabilizing state feedback matrix.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     G       (input) DOUBLE PRECISION array, dimension (LDG,P)
C             The leading N-by-P part of this array must contain a
C             stabilizing observer gain matrix.
C
C     LDG     INTEGER
C             The leading dimension of array G.  LDG >= MAX(1,N).
C
C     SCALEC  (output) DOUBLE PRECISION
C             Scaling factor for the controllability Grammian.
C             See METHOD.
C
C     SCALEO  (output) DOUBLE PRECISION
C             Scaling factor for the observability Grammian.
C             See METHOD.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor Su of frequency-weighted
C             cotrollability Grammian P = Su*Su'. See METHOD.
C
C     LDS     INTEGER
C             The leading dimension of the array S.  LDS >= MAX(1,N).
C
C     R       (output) DOUBLE PRECISION array, dimension (LDR,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor Ru of the frequency-weighted
C             observability Grammian Q = Ru'*Ru. See METHOD.
C
C     LDR     INTEGER
C             The leading dimension of the array R.  LDR >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N*(N + MAX(N,M) + MIN(N,M) + 6)),
C                                                       if JOBCF = 'L';
C             LDWORK >= MAX(1, N*(N + MAX(N,P) + MIN(N,P) + 6)),
C                                                       if JOBCF = 'R'.
C             For optimum performance LDWORK should be larger.
C             An upper bound for both cases is
C             LDWORK >= MAX(1, N*(N + MAX(N,M,P) + 7)).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  eigenvalue computation failure;
C             = 2:  the matrix A+G*C is not stable;
C             = 3:  the matrix A+B*F is not stable;
C             = 4:  the Lyapunov equation for computing the
C                   observability Grammian is (nearly) singular;
C             = 5:  the Lyapunov equation for computing the
C                   controllability Grammian is (nearly) singular.
C
C     METHOD
C
C     In accordance with the type of the coprime factorization
C     of the controller (left or right), the Cholesky factors Su and Ru
C     of the frequency-weighted controllability Grammian P = Su*Su' and
C     of the frequency-weighted observability Grammian Q = Ru'*Ru are
C     computed by solving appropriate Lyapunov or Stein equations [1].
C
C     If JOBCF = 'L' and DICO = 'C', P and Q are computed as the
C     solutions of the following Lyapunov equations:
C
C            (A+B*F)*P + P*(A+B*F)' +  scalec^2*B*B' = 0,  (1)
C
C            (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*F'*F = 0.  (2)
C
C     If JOBCF = 'L' and DICO = 'D', P and Q are computed as the
C     solutions of the following Stein equations:
C
C            (A+B*F)*P*(A+B*F)' - P +  scalec^2*B*B' = 0,  (3)
C
C            (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*F'*F = 0.  (4)
C
C     If JOBCF = 'R' and DICO = 'C', P and Q are computed as the
C     solutions of the following Lyapunov equations:
C
C            (A+B*F)*P + P*(A+B*F)' +  scalec^2*G*G' = 0,  (5)
C
C            (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*C'*C = 0.  (6)
C
C     If JOBCF = 'R' and DICO = 'D', P and Q are computed as the
C     solutions of the following Stein equations:
C
C            (A+B*F)*P*(A+B*F)' - P +  scalec^2*G*G' = 0,  (7)
C
C            (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*C'*C = 0.  (8)
C
C     REFERENCES
C
C     [1] Liu, Y., Anderson, B.D.O. and Ly, O.L.
C         Coprime factorization controller reduction with Bezout
C         identity induced frequency weighting.
C         Automatica, vol. 26, pp. 233-249, 1990.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, October 2000.
C     D. Sima, University of Bucharest, October 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C
C     KEYWORDS
C
C     Controller reduction, frequency weighting, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBCF
      INTEGER          INFO, LDA, LDB, LDC, LDF, LDG, LDR, LDS, LDWORK,
     $                 M, N, P
      DOUBLE PRECISION SCALEC, SCALEO
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*),
     $                 F(LDF,*), G(LDG,*), R(LDR,*), S(LDS,*)
C     .. Local Scalars ..
      LOGICAL          DISCR, LEFTW
      INTEGER          IERR, KAW, KU, KW, KWI, KWR, LDU, LW, ME, MP,
     $                 WRKOPT
C     .. External Functions ..
      LOGICAL          LSAME
      EXTERNAL         LSAME
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DLACPY, SB03OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        INT, MAX, MIN
C     .. Executable Statements ..
C
      DISCR = LSAME( DICO,  'D' )
      LEFTW = LSAME( JOBCF, 'L' )
C
      INFO = 0
      IF( LEFTW ) THEN
         MP = M
      ELSE
         MP = P
      END IF
      LW = N*( N + MAX( N, MP ) + MIN( N, MP ) + 6 )
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LEFTW .OR. LSAME( JOBCF, 'R' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDS.LT.MAX( 1, N ) ) THEN
         INFO = -19
      ELSE IF( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -21
      ELSE IF( LDWORK.LT.MAX( 1, LW ) ) THEN
         INFO = -23
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB16CY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         SCALEC   = ONE
         SCALEO   = ONE
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Allocate storage for work arrays.
C
      KAW = 1
      KU  = KAW + N*N
      KWR = KU  + N*MAX( N, MP )
      KWI = KWR + N
      KW  = KWI + N
C
C     Form A+G*C.
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK(KAW), N )
      CALL DGEMM( 'No-transpose', 'No-transpose', N, N, P, ONE,
     $            G, LDG, C, LDC, ONE, DWORK(KAW), N )
C
C     Form the factor H of the free term.
C
      IF( LEFTW ) THEN
C
C        H = F.
C
         LDU = MAX( N, M )
         ME  = M
         CALL DLACPY( 'Full', M, N, F, LDF, DWORK(KU), LDU )
      ELSE
C
C        H = C.
C
         LDU = MAX( N, P )
         ME  = P
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), LDU )
      END IF
C
C     Solve for the Cholesky factor Ru of Q, Q = Ru'*Ru,
C     the continuous-time Lyapunov equation (if DICO = 'C')
C
C        (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*H'*H = 0,
C
C     or the discrete-time Lyapunov equation (if DICO = 'D')
C
C        (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*H'*H = 0.
C
C     Workspace:  need   N*(N + MAX(N,M) + MIN(N,M) + 6) if JOBCF = 'L';
C                        N*(N + MAX(N,P) + MIN(N,P) + 6) if JOBCF = 'R'.
C                 prefer larger.
C
      CALL SB03OD( DICO, 'NoFact', 'NoTransp', N, ME, DWORK(KAW), N,
     $             R, LDR, DWORK(KU), LDU, SCALEO, DWORK(KWR),
     $             DWORK(KWI), DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.EQ.2 ) THEN
            INFO = 2
         ELSE IF( IERR.EQ.1 ) THEN
            INFO = 4
         ELSE IF( IERR.EQ.6 ) THEN
            INFO = 1
         END IF
         RETURN
      END IF
C
      WRKOPT = INT( DWORK(KW) ) + KW - 1
      CALL DLACPY( 'Upper', N, N, DWORK(KU), LDU, R, LDR )
C
C     Form A+B*F.
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK(KAW), N )
      CALL DGEMM( 'No-transpose', 'No-transpose', N, N, M, ONE,
     $            B, LDB, F, LDF, ONE, DWORK(KAW), N )
C
C     Form the factor K of the free term.
C
      LDU = N
      IF( LEFTW ) THEN
C
C        K = B.
C
         ME = M
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), LDU )
      ELSE
C
C        K = G.
C
         ME = P
         CALL DLACPY( 'Full', N, P, G, LDG, DWORK(KU), LDU )
      END IF
C
C     Solve for the Cholesky factor Su of P, P = Su*Su',
C     the continuous-time Lyapunov equation (if DICO = 'C')
C
C         (A+B*F)*P + P*(A+B*F)' +  scalec^2*K*K' = 0,
C
C     or the discrete-time Lyapunov equation (if DICO = 'D')
C
C         (A+B*F)*P*(A+B*F)' - P +  scalec^2*K*K' = 0.
C
C     Workspace:  need   N*(N + MAX(N,M) + MIN(N,M) + 6) if JOBCF = 'L';
C                        N*(N + MAX(N,P) + MIN(N,P) + 6) if JOBCF = 'R'.
C                        prefer larger.
C
      CALL SB03OD( DICO, 'NoFact', 'Transp', N, ME, DWORK(KAW), N,
     $             S, LDS, DWORK(KU), LDU, SCALEC, DWORK(KWR),
     $             DWORK(KWI), DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.EQ.2 ) THEN
            INFO = 3
         ELSE IF( IERR.EQ.1 ) THEN
            INFO = 5
         ELSE IF( IERR.EQ.6 ) THEN
            INFO = 1
         END IF
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      CALL DLACPY( 'Upper', N, N, DWORK(KU), LDU, S, LDS )
C
C     Save the optimal workspace.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB16CY ***
      END
