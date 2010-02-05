      SUBROUTINE SB16BD( DICO, JOBD, JOBMR, JOBCF, EQUIL, ORDSEL,
     $                   N, M, P, NCR, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   F, LDF, G, LDG, DC, LDDC, HSV, TOL1, TOL2,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     To compute, for a given open-loop model (A,B,C,D), and for
C     given state feedback gain F and full observer gain G,
C     such that A+B*F and A+G*C are stable, a reduced order
C     controller model (Ac,Bc,Cc,Dc) using a coprime factorization
C     based controller reduction approach. For reduction,
C     either the square-root or the balancing-free square-root
C     versions of the Balance & Truncate (B&T) or Singular Perturbation
C     Approximation (SPA) model reduction methods are used in
C     conjunction with stable coprime factorization techniques.
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
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears
C             in the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     JOBMR   CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root B&T method;
C             = 'F':  use the balancing-free square-root B&T method;
C             = 'S':  use the square-root SPA method;
C             = 'P':  use the balancing-free square-root SPA method.
C
C     JOBCF   CHARACTER*1
C             Specifies whether left or right coprime factorization is
C             to be used as follows:
C             = 'L':  use left coprime factorization;
C             = 'R':  use right coprime factorization.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to perform a
C             preliminary equilibration before performing
C             order reduction as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting controller order NCR is fixed;
C             = 'A':  the resulting controller order NCR is
C                     automatically determined on basis of the given
C                     tolerance TOL1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the open-loop state-space representation,
C             i.e., the order of the matrix A.  N >= 0.
C             N also represents the order of the original state-feedback
C             controller.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NCR     (input/output) INTEGER
C             On entry with ORDSEL = 'F', NCR is the desired order of
C             the resulting reduced order controller.  0 <= NCR <= N.
C             On exit, if INFO = 0, NCR is the order of the resulting
C             reduced order controller. NCR is set as follows:
C             if ORDSEL = 'F', NCR is equal to MIN(NCR,NMIN), where NCR
C             is the desired order on entry, and NMIN is the order of a
C             minimal realization of an extended system Ge (see METHOD);
C             NMIN is determined as the number of
C             Hankel singular values greater than N*EPS*HNORM(Ge),
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH) and HNORM(Ge) is the Hankel norm of the
C             extended system (computed in HSV(1));
C             if ORDSEL = 'A', NCR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS*HNORM(Ge)).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, if INFO = 0, the leading NCR-by-NCR part of this
C             array contains the state dynamics matrix Ac of the reduced
C             controller.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must
C             contain the original input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must
C             contain the original state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             If JOBD = 'D', the leading P-by-M part of this
C             array must contain the system direct input/output
C             transmission matrix D.
C             The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P), if JOBD = 'D';
C             LDD >= 1,        if JOBD = 'Z'.
C
C     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N)
C             On entry, the leading M-by-N part of this array must
C             contain a stabilizing state feedback matrix.
C             On exit, if INFO = 0, the leading M-by-NCR part of this
C             array contains the state/output matrix Cc of the reduced
C             controller.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,P)
C             On entry, the leading N-by-P part of this array must
C             contain a stabilizing observer gain matrix.
C             On exit, if INFO = 0, the leading NCR-by-P part of this
C             array contains the input/state matrix Bc of the reduced
C             controller.
C
C     LDG     INTEGER
C             The leading dimension of array G.  LDG >= MAX(1,N).
C
C     DC      (output) DOUBLE PRECISION array, dimension (LDDC,P)
C             If INFO = 0, the leading M-by-P part of this array
C             contains the input/output matrix Dc of the reduced
C             controller.
C
C     LDDC    INTEGER
C             The leading dimension of array DC.  LDDC >= MAX(1,M).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the N Hankel singular values
C             of the extended system ordered decreasingly (see METHOD).
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of the reduced extended system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(Ge), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(Ge) is the
C             Hankel norm of the extended system (computed in HSV(1)).
C             The value TOL1 = N*EPS*HNORM(Ge) is used by default if
C             TOL1 <= 0 on entry, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the coprime factorization controller
C             (see METHOD). The recommended value is
C             TOL2 = N*EPS*HNORM(Ge) (see METHOD).
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = 0,         if ORDSEL = 'F' and NCR = N.
C                                                 Otherwise,
C             LIWORK = MAX(PM,M), if JOBCF = 'L',
C             LIWORK = MAX(PM,P), if JOBCF = 'R', where
C             PM = 0,             if JOBMR = 'B',
C             PM = N,             if JOBMR = 'F',
C             PM = MAX(1,2*N),    if JOBMR = 'S' or 'P'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= P*N, if ORDSEL = 'F' and NCR = N. Otherwise,
C             LDWORK >= (N+M)*(M+P) + MAX(LWR,4*M), if JOBCF = 'L',
C             LDWORK >= (N+P)*(M+P) + MAX(LWR,4*P), if JOBCF = 'R',
C             where LWR = MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NCR is
C                   greater than the order of a minimal
C                   realization of the controller.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A+G*C to a real Schur form
C                   failed;
C             = 2:  the matrix A+G*C is not stable (if DICO = 'C'),
C                   or not convergent (if DICO = 'D');
C             = 3:  the computation of Hankel singular values failed;
C             = 4:  the reduction of A+B*F to a real Schur form
C                   failed;
C             = 5:  the matrix A+B*F is not stable (if DICO = 'C'),
C                   or not convergent (if DICO = 'D').
C
C     METHOD
C
C     Let be the linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t),                             (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system, and let Go(d) be the open-loop
C     transfer-function matrix
C                           -1
C          Go(d) = C*(d*I-A) *B + D .
C
C     Let F and G be the state feedback and observer gain matrices,
C     respectively, chosen so that A+B*F and A+G*C are stable matrices.
C     The controller has a transfer-function matrix K(d) given by
C                                        -1
C          K(d) = F*(d*I-A-B*F-G*C-G*D*F) *G .
C
C     The closed-loop transfer-function matrix is given by
C                                     -1
C          Gcl(d) = Go(d)(I+K(d)Go(d)) .
C
C     K(d) can be expressed as a left coprime factorization (LCF),
C                          -1
C          K(d) = M_left(d) *N_left(d) ,
C
C     or as a right coprime factorization (RCF),
C                                      -1
C          K(d) = N_right(d)*M_right(d) ,
C
C     where M_left(d), N_left(d), N_right(d), and M_right(d) are
C     stable transfer-function matrices.
C
C     The subroutine SB16BD determines the matrices of a reduced
C     controller
C
C          d[z(t)] = Ac*z(t) + Bc*y(t)
C          u(t)    = Cc*z(t) + Dc*y(t),                           (2)
C
C     with the transfer-function matrix Kr as follows:
C
C     (1) If JOBCF = 'L', the extended system
C         Ge(d)  = [ N_left(d) M_left(d) ] is reduced to
C         Ger(d) = [ N_leftr(d) M_leftr(d) ] by using either the
C         B&T or SPA methods. The reduced order controller Kr(d)
C         is computed as
C                           -1
C         Kr(d) = M_leftr(d) *N_leftr(d) ;
C
C     (2) If JOBCF = 'R', the extended system
C         Ge(d) = [ N_right(d) ] is reduced to
C                 [ M_right(d) ]
C         Ger(d) = [ N_rightr(d) ] by using either the
C                  [ M_rightr(d) ]
C         B&T or SPA methods. The reduced order controller Kr(d)
C         is computed as
C                                         -1
C         Kr(d) = N_rightr(d)* M_rightr(d) .
C
C     If ORDSEL = 'A', the order of the controller is determined by
C     computing the number of Hankel singular values greater than
C     the given tolerance TOL1. The Hankel singular values are
C     the square roots of the eigenvalues of the product of
C     the controllability and observability Grammians of the
C     extended system Ge.
C
C     If JOBMR = 'B', the square-root B&T method of [1] is used.
C
C     If JOBMR = 'F', the balancing-free square-root version of the
C     B&T method [1] is used.
C
C     If JOBMR = 'S', the square-root version of the SPA method [2,3]
C     is used.
C
C     If JOBMR = 'P', the balancing-free square-root version of the
C     SPA method [2,3] is used.
C
C     REFERENCES
C
C     [1] Tombs, M.S. and Postlethwaite, I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [2] Varga, A.
C         Efficient minimal realization procedure based on balancing.
C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), Vol. 2,
C         pp. 42-46, 1991.
C
C     [3] Varga, A.
C         Coprime factors model reduction method based on square-root
C         balancing-free techniques.
C         System Analysis, Modelling and Simulation, Vol. 11,
C         pp. 303-311, 1993.
C
C     [4] Liu, Y., Anderson, B.D.O. and Ly, O.L.
C         Coprime factorization controller reduction with Bezout
C         identity induced frequency weighting.
C         Automatica, vol. 26, pp. 233-249, 1990.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root or
C     balancing-free square-root techniques.
C                                         3
C     The algorithms require less than 30N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000.
C     D. Sima, University of Bucharest, August 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C              Aug. 2001.
C
C     KEYWORDS
C
C     Balancing, controller reduction, coprime factorization,
C     minimal realization, multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOBCF, JOBD, JOBMR, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDDC,
     $                  LDF, LDG, LDWORK, M, N, NCR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DC(LDDC,*), DWORK(*), F(LDF,*), G(LDG,*), HSV(*)
C     .. Local Scalars ..
      CHARACTER         JOB
      LOGICAL           BAL, BTA, DISCR, FIXORD, LEFT, LEQUIL, SPA,
     $                  WITHD
      INTEGER           KBE, KCE, KDE, KW, LDBE, LDCE, LDDE, LW1, LW2,
     $                  LWR, MAXMP, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09AD, AB09BD, DGEMM, DLACPY, DLASET, SB08GD,
     $                  SB08HD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      WITHD  = LSAME( JOBD,   'D' )
      BTA    = LSAME( JOBMR,  'B' ) .OR. LSAME( JOBMR, 'F' )
      SPA    = LSAME( JOBMR,  'S' ) .OR. LSAME( JOBMR, 'P' )
      BAL    = LSAME( JOBMR,  'B' ) .OR. LSAME( JOBMR, 'S' )
      LEFT   = LSAME( JOBCF,  'L' )
      LEQUIL = LSAME( EQUIL,  'S' )
      FIXORD = LSAME( ORDSEL, 'F' )
      MAXMP  = MAX( M, P )
C
      LWR = MAX( 1, N*( 2*N + MAX( N, M+P ) + 5 ) + ( N*(N+1) )/2 )
      LW1 = (N+M)*(M+P) + MAX( LWR, 4*M )
      LW2 = (N+P)*(M+P) + MAX( LWR, 4*P )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( WITHD .OR. LSAME( JOBD, 'Z' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LEFT .OR. LSAME( JOBCF, 'R' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( LEQUIL .OR. LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( M.LT.0 ) THEN
         INFO = -8
      ELSE IF( P.LT.0 ) THEN
         INFO = -9
      ELSE IF( FIXORD .AND. ( NCR.LT.0 .OR. NCR.GT.N ) ) THEN
         INFO = -10
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -20
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -22
      ELSE IF( LDDC.LT.MAX( 1, M ) ) THEN
         INFO = -24
      ELSE IF( .NOT.FIXORD .AND. TOL2.GT.ZERO .AND. TOL2.GT.TOL1 ) THEN
         INFO = -27
      ELSE IF( ( ( .NOT.FIXORD .OR. NCR.LT.N ) .AND.
     $         ( ( LEFT .AND. LDWORK.LT.LW1 ) ) .OR.
     $      ( .NOT.LEFT .AND. LDWORK.LT.LW2 ) ) .OR.
     $      ( FIXORD .AND. NCR.EQ.N .AND. LDWORK.LT.P*N ) ) THEN
         INFO = -30
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB16BD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 .OR.
     $    ( FIXORD .AND. BTA .AND. NCR.EQ.0 ) ) THEN
         NCR = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( NCR.EQ.N ) THEN
C
C        Form the controller state matrix,
C        Ac = A + B*F + G*C + G*D*F = A + B*F + G*(C+D*F) .
C        Real workspace:    need  P*N.
C        Integer workspace: need  0.
C
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK, P )
         IF( WITHD ) CALL DGEMM( 'NoTranspose', 'NoTranspose', P, N, M,
     $                          ONE, D, LDD, F, LDF, ONE,
     $                          DWORK, P )
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, P, ONE, G,
     $               LDG, DWORK, P, ONE, A, LDA )
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, M, ONE, B,
     $               LDB, F, LDF, ONE, A, LDA )
C
         DWORK(1) = P*N
         RETURN
      END IF
C
      IF( BAL ) THEN
         JOB = 'B'
      ELSE
         JOB = 'N'
      END IF
C
C     Reduce the coprime factors.
C
      IF( LEFT ) THEN
C
C        Form Ge(d) = [ N_left(d) M_left(d) ] as
C
C             ( A+G*C |  G  B+GD )
C             (------------------)
C             (   F   |  0   I   )
C
C        Real workspace:    need  (N+M)*(M+P).
C        Integer workspace: need  0.
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, P, ONE, G,
     $               LDG, C, LDC, ONE, A, LDA )
         KBE  = 1
         KDE  = KBE + N*(P+M)
         LDBE = MAX( 1, N )
         LDDE = M
         CALL DLACPY( 'Full', N, P, G, LDG, DWORK(KBE), LDBE )
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KBE+N*P), LDBE )
         IF( WITHD ) CALL DGEMM( 'NoTranspose', 'NoTranspose', N, M, P,
     $                          ONE, G, LDG, D, LDD, ONE,
     $                          DWORK(KBE+N*P), LDBE )
         CALL DLASET( 'Full', M, P, ZERO, ZERO, DWORK(KDE), LDDE )
         CALL DLASET( 'Full', M, M, ZERO, ONE, DWORK(KDE+M*P), LDDE )
C
C        Compute the reduced coprime factors,
C             Ger(d) = [ N_leftr(d) M_leftr(d) ] ,
C        by using either the B&T or SPA methods.
C
C        Real workspace:    need  (N+M)*(M+P) +
C                                 MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2).
C        Integer workspace: need  0,         if JOBMR = 'B',
C                                 N,         if JOBMR = 'F', and
C                                 MAX(1,2*N) if JOBMR = 'S' or 'P'.
C
         KW = KDE + M*(P+M)
         IF( BTA ) THEN
            CALL AB09AD( DICO, JOB, EQUIL, ORDSEL, N, M+P, M, NCR, A,
     $                   LDA, DWORK(KBE), LDBE, F, LDF, HSV, TOL1,
     $                   IWORK, DWORK(KW), LDWORK-KW+1, IWARN, INFO )
         ELSE
            CALL AB09BD( DICO, JOB, EQUIL, ORDSEL, N, M+P, M, NCR, A,
     $                   LDA, DWORK(KBE), LDBE, F, LDF, DWORK(KDE),
     $                   LDDE, HSV, TOL1, TOL2, IWORK, DWORK(KW),
     $                   LDWORK-KW+1, IWARN, INFO )
         END IF
         IF( INFO.NE.0 )
     $      RETURN
C
         WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C        Compute the reduced order controller,
C                             -1
C           Kr(d) = M_leftr(d)  *N_leftr(d).
C
C        Real workspace:    need  (N+M)*(M+P) + MAX(1,4*M).
C        Integer workspace: need  M.
C
         CALL SB08GD( NCR, P, M, A, LDA, DWORK(KBE), LDBE, F, LDF,
     $                DWORK(KDE), LDDE, DWORK(KBE+N*P), LDBE,
     $                DWORK(KDE+M*P), LDDE, IWORK, DWORK(KW), INFO )
C
C        Copy the reduced system matrices Bc and Dc.
C
         CALL DLACPY( 'Full', NCR, P, DWORK(KBE), LDBE, G, LDG )
         CALL DLACPY( 'Full', M, P, DWORK(KDE), LDDE, DC, LDDC )
C
      ELSE
C
C        Form Ge(d) = [ N_right(d) ]
C                     [ M_right(d) ] as
C
C             ( A+B*F | G )
C             (-----------)
C             (   F   | 0 )
C             ( C+D*F | I )
C
C        Real workspace:    need  (N+P)*(M+P).
C        Integer workspace: need  0.
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, M, ONE, B,
     $               LDB, F, LDF, ONE, A, LDA )
         KCE  = 1
         KDE  = KCE + N*(P+M)
         LDCE = M+P
         LDDE = LDCE
         CALL DLACPY( 'Full', M, N, F, LDF, DWORK(KCE), LDCE )
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KCE+M), LDCE )
         IF( WITHD ) CALL DGEMM( 'NoTranspose', 'NoTranspose', P, N, M,
     $                          ONE, D, LDD, F, LDF, ONE,
     $                          DWORK(KCE+M), LDCE )
         CALL DLASET( 'Full', M, P, ZERO, ZERO, DWORK(KDE), LDDE )
         CALL DLASET( 'Full', P, P, ZERO, ONE, DWORK(KDE+M), LDDE )
C
C        Compute the reduced coprime factors,
C             Ger(d) = [ N_rightr(d) ]
C                      [ M_rightr(d) ],
C        by using either the B&T or SPA methods.
C
C        Real workspace:    need  (N+P)*(M+P) +
C                                 MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2).
C        Integer workspace: need  0,         if JOBMR = 'B',
C                                 N,         if JOBMR = 'F', and
C                                 MAX(1,2*N) if JOBMR = 'S' or 'P'.
C
         KW = KDE + P*(P+M)
         IF( BTA ) THEN
            CALL AB09AD( DICO, JOB, EQUIL, ORDSEL, N, P, M+P, NCR, A,
     $                   LDA, G, LDG, DWORK(KCE), LDCE, HSV, TOL1,
     $                   IWORK, DWORK(KW), LDWORK-KW+1, IWARN, INFO )
         ELSE
            CALL AB09BD( DICO, JOB, EQUIL, ORDSEL, N, P, M+P, NCR, A,
     $                   LDA, G, LDG, DWORK(KCE), LDCE, DWORK(KDE),
     $                   LDDE, HSV, TOL1, TOL2, IWORK, DWORK(KW),
     $                   LDWORK-KW+1, IWARN, INFO )
         END IF
         IF( INFO.NE.0 ) THEN
            IF( INFO.NE.3 ) INFO = INFO + 3
            RETURN
         END IF
C
         WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C        Compute the reduced order controller,
C                                        -1
C           Kr(d) = N_rightr(d)*M_rightr(d) .
C
C        Real workspace:    need  (N+P)*(M+P) + MAX(1,4*P).
C        Integer workspace: need  P.
C
         CALL SB08HD( NCR, P, M, A, LDA, G, LDG, DWORK(KCE), LDCE,
     $                DWORK(KDE), LDDE, DWORK(KCE+M), LDCE,
     $                DWORK(KDE+M), LDDE, IWORK, DWORK(KW), INFO )
C
C        Copy the reduced system matrices Cc and Dc.
C
         CALL DLACPY( 'Full', M, NCR, DWORK(KCE), LDCE, F, LDF )
         CALL DLACPY( 'Full', M, P, DWORK(KDE), LDDE, DC, LDDC )
C
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB16BD ***
      END
