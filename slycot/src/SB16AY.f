      SUBROUTINE SB16AY( DICO, JOBC, JOBO, WEIGHT, N, M, P, NC, NCS,
     $                   A, LDA, B, LDB, C, LDC, D, LDD,
     $                   AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC,
     $                   SCALEC, SCALEO, S, LDS, R, LDR,
     $                   IWORK, DWORK, LDWORK, INFO )
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
C     To compute for given state-space representations (A,B,C,D) and
C     (Ac,Bc,Cc,Dc) of the transfer-function matrices of the
C     open-loop system G and feedback controller K, respectively,
C     the Cholesky factors of the frequency-weighted
C     controllability and observability Grammians corresponding
C     to a frequency-weighted model reduction problem.
C     The controller must stabilize the closed-loop system.
C     The state matrix Ac must be in a block-diagonal real Schur form
C     Ac = diag(Ac1,Ac2), where Ac1 contains the unstable eigenvalues
C     of Ac and Ac2 contains the stable eigenvalues of Ac.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the systems as follows:
C             = 'C':  G and K are continuous-time systems;
C             = 'D':  G and K are discrete-time systems.
C
C     JOBC    CHARACTER*1
C             Specifies the choice of frequency-weighted controllability
C             Grammian as follows:
C             = 'S': choice corresponding to standard Enns' method [1];
C             = 'E': choice corresponding to the stability enhanced
C                    modified Enns' method of [2].
C
C     JOBO    CHARACTER*1
C             Specifies the choice of frequency-weighted observability
C             Grammian as follows:
C             = 'S': choice corresponding to standard Enns' method [1];
C             = 'E': choice corresponding to the stability enhanced
C                    modified combination method of [2].
C
C     WEIGHT  CHARACTER*1
C             Specifies the type of frequency-weighting, as follows:
C             = 'N':  no weightings are used (V = I, W = I);
C             = 'O':  stability enforcing left (output) weighting
C                               -1
C                     V = (I-G*K) *G is used (W = I);
C             = 'I':  stability enforcing right (input) weighting
C                               -1
C                     W = (I-G*K) *G is used (V = I);
C             = 'P':  stability and performance enforcing weightings
C                               -1                -1
C                     V = (I-G*K) *G ,  W = (I-G*K)  are used.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the open-loop system state-space
C             representation, i.e., the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NC      (input) INTEGER
C             The order of the controller state-space representation,
C             i.e., the order of the matrix AC.  NC >= 0.
C
C     NCS     (input) INTEGER
C             The dimension of the stable part of the controller, i.e.,
C             the order of matrix Ac2.  NC >= NCS >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the system with the transfer-function
C             matrix G.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             input/output matrix D of the open-loop system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     AC      (input) DOUBLE PRECISION array, dimension (LDAC,NC)
C             The leading NC-by-NC part of this array must contain
C             the state dynamics matrix Ac of the controller in a
C             block diagonal real Schur form Ac = diag(Ac1,Ac2), where
C             Ac1 is (NC-NCS)-by-(NC-NCS) and contains the unstable
C             eigenvalues of Ac, and Ac2 is NCS-by-NCS and contains
C             the stable eigenvalues of Ac.
C
C     LDAC    INTEGER
C             The leading dimension of array AC.  LDAC >= MAX(1,NC).
C
C     BC      (input) DOUBLE PRECISION array, dimension (LDBC,P)
C             The leading NC-by-P part of this array must contain
C             the input/state matrix Bc of the controller.
C
C     LDBC    INTEGER
C             The leading dimension of array BC.  LDBC >= MAX(1,NC).
C
C     CC      (input) DOUBLE PRECISION array, dimension (LDCC,NC)
C             The leading M-by-NC part of this array must contain
C             the state/output matrix Cc of the controller.
C
C     LDCC    INTEGER
C             The leading dimension of array CC.  LDCC >= MAX(1,M).
C
C     DC      (input) DOUBLE PRECISION array, dimension (LDDC,P)
C             The leading M-by-P part of this array must contain
C             the input/output matrix Dc of the controller.
C
C     LDDC    INTEGER
C             The leading dimension of array DC.  LDDC >= MAX(1,M).
C
C     SCALEC  (output) DOUBLE PRECISION
C             Scaling factor for the controllability Grammian.
C             See METHOD.
C
C     SCALEO  (output) DOUBLE PRECISION
C             Scaling factor for the observability Grammian. See METHOD.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,NCS)
C             The leading NCS-by-NCS upper triangular part of this array
C             contains the Cholesky factor S of the frequency-weighted
C             controllability Grammian P = S*S'. See METHOD.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,NCS).
C
C     R       (output) DOUBLE PRECISION array, dimension (LDR,NCS)
C             The leading NCS-by-NCS upper triangular part of this array
C             contains the Cholesky factor R of the frequency-weighted
C             observability Grammian Q = R'*R. See METHOD.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,NCS).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(LIWRK)
C             LIWRK = 0,       if WEIGHT = 'N';
C             LIWRK = 2(M+P),  if WEIGHT = 'O', 'I', or 'P'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, LFREQ ),
C             where
C             LFREQ = (N+NC)*(N+NC+2*M+2*P)+
C                     MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7), (M+P)*(M+P+4))
C                                      if WEIGHT = 'I' or 'O' or 'P';
C             LFREQ  = NCS*(MAX(M,P)+5) if WEIGHT = 'N'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the closed-loop system is not well-posed;
C                   its feedthrough matrix is (numerically) singular;
C             = 2:  the computation of the real Schur form of the
C                   closed-loop state matrix failed;
C             = 3:  the closed-loop state matrix is not stable;
C             = 4:  the solution of a symmetric eigenproblem failed;
C             = 5:  the NCS-by-NCS trailing part Ac2 of the state
C                   matrix Ac is not stable or not in a real Schur form.
C
C     METHOD
C
C     If JOBC = 'S', the controllability Grammian P is determined as
C     follows:
C
C     - if WEIGHT = 'O' or 'N', P satisfies for a continuous-time
C       controller the Lyapunov equation
C
C            Ac2*P + P*Ac2' +  scalec^2*Bc*Bc' = 0
C
C       and for a discrete-time controller
C
C            Ac2*P*Ac2' - P +  scalec^2*Bc*Bc' = 0;
C
C     - if WEIGHT = 'I' or 'P', let Pi be the solution of the
C       continuous-time Lyapunov equation
C
C            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0
C
C       or of the discrete-time Lyapunov equation
C
C            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0,
C
C       where Ai and Bi are the state and input matrices of a special
C       state-space realization of the input frequency weight (see [2]);
C       P results as the trailing NCS-by-NCS part of Pi partitioned as
C
C           Pi = ( *  * ).
C                ( *  P )
C
C     If JOBC = 'E', a modified controllability Grammian P1 >= P is
C     determined to guarantee stability for a modified Enns' method [2].
C
C     If JOBO = 'S', the observability Grammian Q is determined as
C     follows:
C
C     - if WEIGHT = 'I' or 'N', Q satisfies for a continuous-time
C       controller the Lyapunov equation
C
C            Ac2'*Q + Q*Ac2 +  scaleo^2*Cc'*Cc = 0
C
C       and for a discrete-time controller
C
C            Ac2'*Q*Ac2 - Q +  scaleo^2*Cc'*Cc = 0;
C
C     - if WEIGHT = 'O' or 'P', let Qo be the solution of the
C       continuous-time Lyapunov equation
C
C            Ao'*Qo + Qo*Ao +  scaleo^2*Co'*Co = 0
C
C       or of the discrete-time Lyapunov equation
C
C            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0,
C
C       where Ao and Co are the state and output matrices of a
C       special state-space realization of the output frequency weight
C       (see [2]); if WEIGHT = 'O', Q results as the leading NCS-by-NCS
C       part of Qo partitioned as
C
C           Qo = ( Q  * )
C                ( *  * )
C
C       while if WEIGHT = 'P', Q results as the trailing NCS-by-NCS
C       part of Qo partitioned as
C
C           Qo = ( *  * ).
C                ( *  Q )
C
C     If JOBO = 'E', a modified observability Grammian Q1 >= Q is
C     determined to guarantee stability for a modified Enns' method [2].
C
C     The routine computes directly the Cholesky factors S and R
C     such that P = S*S' and Q = R'*R according to formulas
C     developed in [2].
C
C     REFERENCES
C
C     [1] Enns, D.
C         Model reduction with balanced realizations: An error bound
C         and a frequency weighted generalization.
C         Proc. CDC, Las Vegas, pp. 127-132, 1984.
C
C     [2] Varga, A. and Anderson, B.D.O.
C         Frequency-weighted balancing related controller reduction.
C         Proceedings of the 15th IFAC World Congress, July 21-26, 2002,
C         Barcelona, Spain, Vol.15, Part 1, 2002-07-21.
C
C     CONTRIBUTORS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C     May 2009.
C     A. Varga, DLR Oberpfafenhofen, June 2001.
C
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
      CHARACTER        DICO, JOBC, JOBO, WEIGHT
      INTEGER          INFO, LDA, LDAC, LDB, LDBC, LDC, LDCC, LDD, LDDC,
     $                 LDR, LDS, LDWORK, M, N, NC, NCS, P
      DOUBLE PRECISION SCALEC, SCALEO
C     .. Array Arguments ..
      INTEGER          IWORK(*)
      DOUBLE PRECISION A(LDA,*), AC(LDAC,*), B(LDB,*), BC(LDBC,*),
     $                 C(LDC,*), CC(LDCC,*), D(LDD,*), DC(LDDC,*),
     $                 DWORK(*), R(LDR,*),   S(LDS,*)
C     .. Local Scalars ..
      CHARACTER        JOBFAC
      LOGICAL          DISCR, FRWGHT, LEFTW, PERF, RIGHTW
      INTEGER          I, IERR, J, JJ, KI, KL, KQ, KR, KTAU, KU, KW,
     $                 KWA, KWB, KWC, KWD, LDU, LW, MBBAR, ME, MP,
     $                 NCU, NCU1, NE, NNC, NNCU, PCBAR, PE, WRKOPT
      DOUBLE PRECISION RCOND, T, TOL
C     .. Local Arrays ..
      DOUBLE PRECISION DUM(1)
C     .. External Functions ..
      LOGICAL          LSAME
      DOUBLE PRECISION DLAMCH
      EXTERNAL         DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL         AB05PD, AB05QD, AB07ND, DCOPY, DLACPY, DLASET,
     $                 DSCAL, DSYEV, MB01WD, MB04OD, SB03OD, SB03OU,
     $                 XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        ABS, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      DISCR  = LSAME( DICO,   'D' )
      LEFTW  = LSAME( WEIGHT, 'O' )
      RIGHTW = LSAME( WEIGHT, 'I' )
      PERF   = LSAME( WEIGHT, 'P' )
      FRWGHT = LEFTW .OR. RIGHTW .OR. PERF
C
      INFO = 0
      NNC  = N + NC
      MP   = M + P
      IF( FRWGHT ) THEN
         LW = NNC*( NNC + 2*MP ) +
     $        MAX( NNC*( NNC + MAX( NNC, M, P ) + 7 ), MP*( MP + 4 ) )
      ELSE
         LW = NCS*( MAX( M, P ) + 5 )
      END IF
      LW = MAX( 1, LW )
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( JOBC, 'S' ) .OR. LSAME( JOBC, 'E' ) ) )
     $     THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( JOBO, 'S' ) .OR. LSAME( JOBO, 'E' ) ) )
     $     THEN
         INFO = -3
      ELSE IF( .NOT.( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( NC.LT.0 ) THEN
         INFO = -8
      ELSE IF( NCS.LT.0 .OR. NCS.GT.NC ) THEN
         INFO = -9
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -17
      ELSE IF( LDAC.LT.MAX( 1, NC ) ) THEN
         INFO = -19
      ELSE IF( LDBC.LT.MAX( 1, NC ) ) THEN
         INFO = -21
      ELSE IF( LDCC.LT.MAX( 1, M ) ) THEN
         INFO = -23
      ELSE IF( LDDC.LT.MAX( 1, M ) ) THEN
         INFO = -25
      ELSE IF( LDS.LT.MAX( 1, NCS ) ) THEN
         INFO = -29
      ELSE IF( LDR.LT.MAX( 1, NCS ) ) THEN
         INFO = -31
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -34
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB16AY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      SCALEC = ONE
      SCALEO = ONE
      IF( MIN( NCS, M, P ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      WRKOPT = 1
      NCU  = NC - NCS
      NCU1 = NCU + 1
C
      IF( .NOT.PERF ) THEN
C
C        Compute the Grammians in the case of no weighting or
C        one-sided weighting.
C
         IF( LEFTW .OR. LSAME( WEIGHT, 'N' ) ) THEN
C
C           Compute the standard controllability Grammian.
C
C           Solve for the Cholesky factor S of P, P = S*S',
C           the continuous-time Lyapunov equation (if DICO = 'C')
C
C               Ac2*P + P*Ac2' +  scalec^2*Bc2*Bc2' = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C
C               Ac2*P*Ac2' - P +  scalec^2*Bc2*Bc2' = 0,
C
C           where Bc2 is the matrix formed from the last NCS rows of Bc.
C
C           Workspace:  need   NCS*(P+5);
C                              prefer larger.
            KU   = 1
            KTAU = KU + NCS*P
            KW   = KTAU + NCS
C
            CALL DLACPY( 'Full', NCS, P, BC(NCU1,1), LDBC,
     $                   DWORK(KU), NCS )
            CALL SB03OU( DISCR, .TRUE., NCS, P, AC(NCU1,NCU1), LDAC,
     $                   DWORK(KU), NCS, DWORK(KTAU), S, LDS, SCALEC,
     $                   DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 5
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         END IF
C
         IF( RIGHTW .OR. LSAME( WEIGHT, 'N' ) ) THEN
C
C           Compute the standard observability Grammian.
C
C           Solve for the Cholesky factor R of Q, Q = R'*R,
C           the continuous-time Lyapunov equation (if DICO = 'C')
C
C               Ac2'*Q + Q*Ac2  +  scaleo^2*Cc2'*Cc2 = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C
C               Ac2'*Q*Ac2 - Q +  scaleo^2*Cc2'*Cc2 = 0,
C
C           where Cc2 is the matrix formed from the last NCS columns
C           of Cc.
C
C           Workspace:  need   NCS*(M + 5);
C                              prefer larger.
            KU   = 1
            KTAU = KU + M*NCS
            KW   = KTAU + NCS
C
            CALL DLACPY( 'Full', M, NCS, CC(1,NCU1), LDCC,
     $                   DWORK(KU), M )
            CALL SB03OU( DISCR, .FALSE., NCS, M, AC(NCU1,NCU1), LDAC,
     $                   DWORK(KU), M, DWORK(KTAU), R, LDR, SCALEO,
     $                   DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 5
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         END IF
C
C        Finish if there are no weights.
C
         IF( LSAME( WEIGHT, 'N' ) ) THEN
            DWORK(1) = WRKOPT
            RETURN
         END IF
      END IF
C
      IF( FRWGHT ) THEN
C
C        Allocate working storage for computing the weights.
C
C        Real workspace:    need MAX(1,NNC*NNC+2*NNC*MP+MP*(MP+4));
C        Integer workspace: need 2*MP.
C
         KWA = 1
         KWB = KWA + NNC*NNC
         KWC = KWB + NNC*MP
         KWD = KWC + NNC*MP
         KW  = KWD + MP*MP
         KL  = KWD
C
         IF( LEFTW ) THEN
C
C           Build the extended matrices
C
C           Ao = ( Ac+Bc*inv(R)*D*Cc   Bc*inv(R)*C   ),
C                (     B*inv(Rt)*Cc  A+B*Dc*inv(R)*C )
C
C           Co = ( -inv(R)*D*Cc  -inv(R)*C ) ,
C
C           where  R = I-D*Dc and Rt = I-Dc*D.
C                             -1
C           Method: Compute Ge  = ( Ge11 Ge12 ), where Ge = ( K   -Im ).
C                                 ( Ge21 Ge22 )             ( -Ip  G  )
C
C                               -1
C           Then  Ge11 = -(I-G*K) *G .
C
C           Construct first Ge = (  K  -Im ) such that the stable part
C                                ( -Ip  G  )
C           of K is in the leading position (to avoid updating of
C           QR factorization).
C
            CALL DLASET( 'Full', M, P, ZERO, ZERO, DWORK(KWD), MP )
            CALL AB05PD( 'N', NCS, P, M, NCU, ONE,
     $                   AC(NCU1,NCU1), LDAC, BC(NCU1,1), LDBC,
     $                   CC(1,NCU1), LDCC, DWORK(KWD), MP,
     $                   AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC,
     $                   NE, DWORK(KWA), NNC, DWORK(KWB), NNC,
     $                   DWORK(KWC), MP, DWORK(KWD), MP, IERR )
            CALL AB05QD( 'Over', NC, P, M, N, M, P, DWORK(KWA), NNC,
     $                   DWORK(KWB), NNC, DWORK(KWC), MP, DWORK(KWD),
     $                   MP, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   NE, ME, PE, DWORK(KWA), NNC, DWORK(KWB), NNC,
     $                   DWORK(KWC), MP, DWORK(KWD), MP, IERR )
            CALL DLASET( 'Full', M, M, ZERO, -ONE, DWORK(KWD+MP*P), MP )
            CALL DLASET( 'Full', P, P, ZERO, -ONE, DWORK(KWD+M), MP )
C
         ELSE
C
C           Build the extended matrices
C
C           Ai = ( A+B*Dc*inv(R)*C   B*inv(Rt)*Cc   ) ,
C                (   Bc*inv(R)*C  Ac+Bc*inv(R)*D*Cc )
C
C           Bi = ( B*Dc*inv(R)    B*inv(Rt)  ) ,
C                ( Bc*inv(R)    Bc*D*inv(Rt) )
C
C           Ci = (  -inv(R)*C   -inv(R)*D*Cc ) , where
C
C           R = I-D*Dc and Rt = I-Dc*D.
C
C                             -1
C           Method: Compute Ge  = ( Ge11 Ge12 ), where Ge = ( G   -Ip ).
C                                 ( Ge21 Ge22 )             ( -Im  K  )
C
C                              -1                     -1
C           Then Ge22 = -(I-G*K) *G and Ge21 = -(I-G*K) .
C
C           Construct first Ge = (  G  -Ip ).
C                                ( -Im  K  )
C
            CALL AB05QD( 'N', N, M, P, NC, P, M, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC,
     $                   NE, ME, PE, DWORK(KWA), NNC, DWORK(KWB), NNC,
     $                   DWORK(KWC), MP, DWORK(KWD), MP, IERR )
            CALL DLASET( 'Full', P, P, ZERO, -ONE, DWORK(KWD+MP*M), MP )
            CALL DLASET( 'Full', M, M, ZERO, -ONE, DWORK(KWD+P), MP )
         END IF
C                  -1
C        Compute Ge   = ( Ge11 Ge12 ).
C                       ( Ge21 Ge22 )
C
C        Additional real workspace: need 4*MP;
C        Integer workspace:         need 2*MP.
C
         CALL AB07ND( NNC, MP, DWORK(KWA), NNC, DWORK(KWB), NNC,
     $                DWORK(KWC), MP, DWORK(KWD), MP, RCOND,
     $                IWORK, DWORK(KW), LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C                     -1   ( A1 | B1  B2  )
C        Partition  Ge   = (--------------) and select appropriate
C                          ( C1 | D11 D12 )
C                          ( C2 | D21 D22 )
C
C        pointers to matrices and column dimensions to define weights.
C
         IF( RIGHTW ) THEN
C
C           Define B2 for Ge22.
C
            ME  = M
            KWB = KWB + NNC*P
         ELSE IF( PERF ) THEN
C
C           Define B1 and C2 for Ge21.
C
            ME  = P
            KWC = KWC + M
         END IF
      END IF
C
      IF( LEFTW .OR. PERF ) THEN
C
C        Compute the frequency-weighted observability Grammian.
C
C        Solve for the Cholesky factor Ro of Qo, Qo = Ro'*Ro,
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            Ao'*Qo + Qo*Ao  +  scaleo^2*Co'*Co = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0.
C
C        Additional workspace:  need   NNC*(NNC+MAX(NNC,P)+7);
C                               prefer larger.
C
         LDU = MAX( NNC, P )
         KU  = KL
         KQ  = KU + NNC*LDU
         KR  = KQ + NNC*NNC
         KI  = KR + NNC
         KW  = KI + NNC
C
         JOBFAC = 'N'
         CALL DLACPY( 'Full', P, NNC, DWORK(KWC), MP, DWORK(KU), LDU )
         CALL SB03OD( DICO, JOBFAC, 'No-transpose', NNC, P,
     $                DWORK(KWA), NNC, DWORK(KQ), NNC, DWORK(KU), LDU,
     $                SCALEO, DWORK(KR), DWORK(KI), DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            IF( IERR.EQ.6 ) THEN
               INFO = 2
            ELSE
               INFO = 3
            END IF
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C        Partition Ro as Ro = ( R11 R12 ).
C                             (  0  R22 )
C
         IF( LEFTW ) THEN
C
C           R = R11 (NCS-by-NCS).
C
            CALL DLACPY( 'Upper', NCS, NCS, DWORK(KU), LDU, R, LDR )
         ELSE
C
C           Compute R such that R'*R = R22'*R22 + R12'*R12, where
C           R22 is NCS-by-NCS and R12 is (N+NCU)-by-NCS.
C           R22 corresponds to the stable part of the controller.
C
            NNCU = N + NCU
            CALL DLACPY( 'Upper', NCS, NCS, DWORK(KU+(LDU+1)*NNCU), LDU,
     $                   R, LDR )
            KTAU = KU
            CALL MB04OD( 'Full', NCS, 0, NNCU, R, LDR,
     $                   DWORK(KU+LDU*NNCU), LDU, DUM, 1, DUM, 1,
     $                   DWORK(KTAU), DWORK(KW) )
C
            DO 10 J = 1, NCS
               IF( R(J,J).LT.ZERO )
     $            CALL DSCAL( NCS-J+1, -ONE, R(J,J), LDR )
   10       CONTINUE
         END IF
      END IF
C
      IF( RIGHTW .OR. PERF ) THEN
C
C        Compute the frequency-weighted controllability Grammian.
C
C        Solve for the Cholesky factor Si of Pi, Pi = Si*Si',
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0.
C
C        Additional workspace:  need   NNC*(NNC+MAX(NNC,P,M)+7);
C                               prefer larger.
C
         KU = KL
         KQ = KU + NNC*MAX( NNC, ME )
         KR = KQ + NNC*NNC
         KI = KR + NNC
         KW = KI + NNC
C
         CALL DLACPY( 'Full', NNC, ME, DWORK(KWB), NNC, DWORK(KU), NNC )
         JOBFAC = 'F'
         IF( RIGHTW ) JOBFAC = 'N'
         CALL SB03OD( DICO, JOBFAC, 'Transpose', NNC, ME,
     $                DWORK(KWA), NNC, DWORK(KQ), NNC, DWORK(KU), NNC,
     $                SCALEC, DWORK(KR), DWORK(KI), DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            IF( IERR.EQ.6 ) THEN
               INFO = 2
            ELSE
               INFO = 3
            END IF
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C        Partition Si as Si = ( S11 S12 ) with S22 NCS-by-NCS and
C                             (  0  S22 )
C        set S = S22.
C
         NNCU = N + NCU
         CALL DLACPY( 'Upper', NCS, NCS, DWORK(KU+(NNC+1)*NNCU), NNC,
     $                S, LDS )
      END IF
C
      KU = 1
      IF( LEFTW .OR. PERF ) THEN
         IF( LSAME( JOBO, 'E' ) ) THEN
C
C           Form Y = -Ac2'*(R'*R)-(R'*R)*Ac2 if DICO = 'C', or
C                Y = -Ac2'*(R'*R)*Ac2+(R'*R) if DICO = 'D'.
C
C           Workspace:  need   2*NCS*NCS.
C
            CALL DLACPY( 'Upper', NCS, NCS, R, LDR, DWORK(KU), NCS )
            CALL DLACPY( 'Full', NCS, NCS, AC(NCU1,NCU1), LDAC,
     $                   DWORK(KU+NCS*NCS), NCS )
            CALL MB01WD( DICO, 'Upper', 'No-transpose', 'Hessenberg',
     $                   NCS, -ONE, ZERO, R, LDR, DWORK(KU+NCS*NCS),
     $                   NCS, DWORK(KU), NCS, IERR )
C
C           Compute the eigendecomposition of Y as Y = Z*Sigma*Z'.
C
            KW = KU + NCS
            CALL DSYEV( 'Vectors', 'Upper', NCS, R, LDR, DWORK(KU),
     $                  DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 4
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C           Partition Sigma = (Sigma1,Sigma2), such that
C           Sigma1 <= 0, Sigma2 > 0.
C           Partition correspondingly Z = [Z1 Z2].
C
            TOL = MAX( ABS( DWORK(KU) ), ABS( DWORK(KU+NCS-1) ) )
     $            * DLAMCH( 'Epsilon')
C                _
C           Form Cc = [ sqrt(Sigma2)*Z2' ]
C
            PCBAR = 0
            JJ = KU
            DO 20 J = 1, NCS
               IF( DWORK(JJ).GT.TOL ) THEN
                  CALL DSCAL( NCS, SQRT( DWORK(JJ) ), R(1,J), 1 )
                  CALL DCOPY( NCS, R(1,J), 1, DWORK(KW+PCBAR), NCS )
                  PCBAR = PCBAR + 1
               END IF
               JJ = JJ + 1
   20       CONTINUE
C
C           Solve for the Cholesky factor R of Q, Q = R'*R,
C           the continuous-time Lyapunov equation (if DICO = 'C')
C                                               _   _
C                   Ac2'*Q + Q*Ac2  +  scaleo^2*Cc'*Cc = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C                                              _   _
C                   Ac2'*Q*Ac2 - Q +  scaleo^2*Cc'*Cc = 0.
C
C           Workspace:  need   NCS*(NCS + 6);
C                              prefer larger.
C
            KU   = KW
            KTAU = KU + NCS*NCS
            KW   = KTAU + NCS
C
            CALL SB03OU( DISCR, .FALSE., NCS, PCBAR, AC(NCU1,NCU1),
     $                   LDAC, DWORK(KU), NCS, DWORK(KTAU), R, LDR, T,
     $                   DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 5
               RETURN
            END IF
            SCALEO = SCALEO*T
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         END IF
C
      END IF
C
      IF( RIGHTW .OR. PERF ) THEN
         IF( LSAME( JOBC, 'E' ) ) THEN
C
C           Form X = -A2c*(S*S')-(S*S')*Ac2' if DICO = 'C', or
C                X = -Ac2*(S*S')*Ac2'+(S*S') if DICO = 'D'.
C
C           Workspace:  need   2*NCS*NCS.
C
            CALL DLACPY( 'Upper', NCS, NCS, S, LDS, DWORK(KU), NCS )
            CALL DLACPY( 'Full', NCS, NCS, AC(NCU1,NCU1), LDAC,
     $                   DWORK(KU+NCS*NCS), NCS )
            CALL MB01WD( DICO, 'Upper', 'Transpose', 'Hessenberg', NCS,
     $                   -ONE, ZERO, S, LDS, DWORK(KU+NCS*NCS), NCS,
     $                   DWORK(KU), NCS, IERR )
C
C           Compute the eigendecomposition of X as X = Z*Sigma*Z'.
C
            KW = KU + NCS
            CALL DSYEV( 'Vectors', 'Upper', NCS, S, LDS, DWORK(KU),
     $                  DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 4
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C           Partition Sigma = (Sigma1,Sigma2), such that
C           Sigma1 =< 0, Sigma2 > 0.
C           Partition correspondingly Z = [Z1 Z2].
C
            TOL = MAX( ABS( DWORK(KU) ), ABS( DWORK(KU+NCS-1) ) )
     $            * DLAMCH( 'Epsilon')
C                _
C           Form Bc = [ Z2*sqrt(Sigma2) ]
C
            MBBAR = 0
            I  = KW
            JJ = KU
            DO 30 J = 1, NCS
               IF( DWORK(JJ).GT.TOL ) THEN
                  MBBAR = MBBAR + 1
                  CALL DSCAL( NCS, SQRT( DWORK(JJ) ), S(1,J), 1 )
                  CALL DCOPY( NCS, S(1,J), 1, DWORK(I), 1 )
                  I = I + NCS
               END IF
               JJ = JJ + 1
   30       CONTINUE
C
C           Solve for the Cholesky factor S of P, P = S*S',
C           the continuous-time Lyapunov equation (if DICO = 'C')
C                                               _  _
C                   Ac2*P + P*Ac2'  +  scalec^2*Bc*Bc' = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C                                              _  _
C                   Ac2*P*Ac2' - P +  scalec^2*Bc*Bc' = 0.
C
C           Workspace:  need   maximum NCS*(NCS + 6);
C                       prefer larger.
C
            KU   = KW
            KTAU = KU + MBBAR*NCS
            KW   = KTAU + NCS
C
            CALL SB03OU( DISCR, .TRUE., NCS, MBBAR, AC(NCU1,NCU1), LDAC,
     $                   DWORK(KU), NCS, DWORK(KTAU), S, LDS, T,
     $                   DWORK(KW), LDWORK-KW+1, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 5
               RETURN
            END IF
            SCALEC = SCALEC*T
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         END IF
C
      END IF
C
C     Save optimal workspace.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB16AY ***
      END
