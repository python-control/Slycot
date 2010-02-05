      SUBROUTINE SB16CD( DICO, JOBD, JOBMR, JOBCF, ORDSEL, N, M, P, NCR,
     $                   A, LDA, B, LDB, C, LDC, D, LDD, F, LDF, G, LDG,
     $                   HSV, TOL, IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     controller model (Ac,Bc,Cc) using a coprime factorization
C     based controller reduction approach. For reduction of
C     coprime factors, a stability enforcing frequency-weighted
C     model reduction is performed using either the square-root or
C     the balancing-free square-root versions of the Balance & Truncate
C     (B&T) model reduction method.
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
C             in the given state space model, as follows:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     JOBMR   CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root B&T method;
C             = 'F':  use the balancing-free square-root B&T method.
C
C     JOBCF   CHARACTER*1
C             Specifies whether left or right coprime factorization
C             of the controller is to be used as follows:
C             = 'L':  use left coprime factorization;
C             = 'R':  use right coprime factorization.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting controller order NCR is fixed;
C             = 'A':  the resulting controller order NCR is
C                     automatically determined on basis of the given
C                     tolerance TOL.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation, i.e.
C             the order of the matrix A.  N >= 0.
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
C             if ORDSEL = 'F', NCR is equal to MIN(NCR,NCRMIN), where
C             NCR is the desired order on entry, and NCRMIN is the
C             number of Hankel-singular values greater than N*EPS*S1,
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH) and S1 is the largest Hankel singular
C             value (computed in HSV(1)); NCR can be further reduced
C             to ensure HSV(NCR) > HSV(NCR+1);
C             if ORDSEL = 'A', NCR is equal to the number of Hankel
C             singular values greater than MAX(TOL,N*EPS*S1).
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
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the open-loop system input/state matrix B.
C             On exit, this array is overwritten with a NCR-by-M
C             B&T approximation of the matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the open-loop system state/output matrix C.
C             On exit, this array is overwritten with a P-by-NCR
C             B&T approximation of the matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, if JOBD = 'D', the leading P-by-M part of this
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
C             array contains the output/state matrix Cc of the reduced
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
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, HSV contains the N frequency-weighted
C             Hankel singular values ordered decreasingly (see METHOD).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If ORDSEL = 'A', TOL contains the tolerance for
C             determining the order of reduced controller.
C             The recommended value is TOL = c*S1, where c is a constant
C             in the interval [0.00001,0.001], and S1 is the largest
C             Hankel singular value (computed in HSV(1)).
C             The value TOL = N*EPS*S1 is used by default if
C             TOL <= 0 on entry, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL is ignored.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension LIWORK, where
C             LIWORK = 0,   if JOBMR = 'B';
C             LIWORK = N,   if JOBMR = 'F'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 2*N*N + MAX( 1, 2*N*N + 5*N, N*MAX(M,P),
C                                    N*(N + MAX(N,MP) + MIN(N,MP) + 6)),
C             where     MP = M, if JOBCF = 'L';
C                       MP = P, if JOBCF = 'R'.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NCR is
C                   greater than the order of a minimal realization
C                   of the controller;
C             = 2:  with ORDSEL = 'F', the selected order NCR
C                   corresponds to repeated singular values, which are
C                   neither all included nor all excluded from the
C                   reduced controller. In this case, the resulting NCR
C                   is set automatically to the largest value such that
C                   HSV(NCR) > HSV(NCR+1).
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
C                   controllability Grammian is (nearly) singular;
C             = 6:  the computation of Hankel singular values failed.
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
C                          -1
C          Go(d) = C*(d*I-A) *B + D .
C
C     Let F and G be the state feedback and observer gain matrices,
C     respectively, chosen such that A+BF and A+GC are stable matrices.
C     The controller has a transfer-function matrix K(d) given by
C                                       -1
C          K(d) = F*(d*I-A-B*F-G*C-G*D*F) *G .
C
C     The closed-loop transfer function matrix is given by
C                                    -1
C          Gcl(d) = Go(d)(I+K(d)Go(d)) .
C
C     K(d) can be expressed as a left coprime factorization (LCF)
C                         -1
C          K(d) = M_left(d) *N_left(d),
C
C     or as a right coprime factorization (RCF)
C                                     -1
C          K(d) = N_right(d)*M_right(d) ,
C
C     where M_left(d), N_left(d), N_right(d), and M_right(d) are
C     stable transfer-function matrices.
C
C     The subroutine SB16CD determines the matrices of a reduced
C     controller
C
C          d[z(t)] = Ac*z(t) + Bc*y(t)
C          u(t)    = Cc*z(t),                                   (2)
C
C     with the transfer-function matrix Kr, using the following
C     stability enforcing approach proposed in [1]:
C
C     (1) If JOBCF = 'L', the frequency-weighted approximation problem
C         is solved
C
C         min||[M_left(d)-M_leftr(d)  N_left(d)-N_leftr(d)][-Y(d)]|| ,
C                                                          [ X(d)]
C         where
C                              -1
C               G(d) = Y(d)*X(d)
C
C         is a RCF of the open-loop system transfer-function matrix.
C         The B&T model reduction technique is used in conjunction
C         with the method proposed in [1].
C
C     (2) If JOBCF = 'R', the frequency-weighted approximation problem
C         is solved
C
C         min || [ -U(d) V(d) ] [ N_right(d)-N_rightr(d) ] || ,
C                               [ M_right(d)-M_rightr(d) ]
C         where
C                         -1
C               G(d) = V(d) *U(d)
C
C         is a LCF of the open-loop system transfer-function matrix.
C         The B&T model reduction technique is used in conjunction
C         with the method proposed in [1].
C
C     If ORDSEL = 'A', the order of the controller is determined by
C     computing the number of Hankel singular values greater than
C     the given tolerance TOL. The Hankel singular values are
C     the square roots of the eigenvalues of the product of
C     two frequency-weighted Grammians P and Q, defined as follows.
C
C     If JOBCF = 'L', then P is the controllability Grammian of a system
C     of the form (A+BF,B,*,*), and Q is the observability Grammian of a
C     system of the form (A+GC,*,F,*). This choice corresponds to an
C     input frequency-weighted order reduction of left coprime
C     factors [1].
C
C     If JOBCF = 'R', then P is the controllability Grammian of a system
C     of the form (A+BF,G,*,*), and Q is the observability Grammian of a
C     system of the form (A+GC,*,C,*). This choice corresponds to an
C     output frequency-weighted order reduction of right coprime
C     factors [1].
C
C     For the computation of truncation matrices, the B&T approach
C     is used in conjunction with accuracy enhancing techniques.
C     If JOBMR = 'B', the square-root B&T method of [2,4] is used.
C     If JOBMR = 'F', the balancing-free square-root version of the
C     B&T method [3,4] is used.
C
C     REFERENCES
C
C     [1] Liu, Y., Anderson, B.D.O. and Ly, O.L.
C         Coprime factorization controller reduction with Bezout
C         identity induced frequency weighting.
C         Automatica, vol. 26, pp. 233-249, 1990.
C
C     [2] Tombs, M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [3] Varga, A.
C         Efficient minimal realization procedure based on balancing.
C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), Vol. 2,
C         pp. 42-46, 1991.
C
C     [4] Varga, A.
C         Coprime factors model reduction method based on square-root
C         balancing-free techniques.
C         System Analysis, Modelling and Simulation, Vol. 11,
C         pp. 303-311, 1993.
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
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, October 2000.
C     D. Sima, University of Bucharest, October 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2001.
C
C     KEYWORDS
C
C     Controller reduction, coprime factorization, frequency weighting,
C     multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOBCF, JOBD, JOBMR, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD,
     $                  LDF, LDG, LDWORK, M, N, NCR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), F(LDF,*), G(LDG,*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           BAL, DISCR, FIXORD, LEFT, WITHD
      INTEGER           IERR, KT, KTI, KW, LW, MP, NMR, WRKOPT
      DOUBLE PRECISION  SCALEC, SCALEO
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09IX, DGEMM, DLACPY, SB16CY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      WITHD  = LSAME( JOBD,   'D' )
      BAL    = LSAME( JOBMR,  'B' )
      LEFT   = LSAME( JOBCF,  'L' )
      FIXORD = LSAME( ORDSEL, 'F' )
      IF( LEFT ) THEN
         MP = M
      ELSE
         MP = P
      END IF
      LW = 2*N*N + MAX( 1, 2*N*N + 5*N, N*MAX( M, P ),
     $                  N*( N + MAX( N, MP ) + MIN( N, MP ) + 6 ) )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( WITHD  .OR. LSAME( JOBD,   'Z' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( BAL    .OR. LSAME( JOBMR,  'F' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LEFT   .OR. LSAME( JOBCF,  'R' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( P.LT.0 ) THEN
         INFO = -8
      ELSE IF( FIXORD .AND. ( NCR.LT.0 .OR. NCR.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -17
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -19
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -21
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -26
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB16CD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 .OR.
     $    ( FIXORD .AND. NCR.EQ.0 ) ) THEN
         NCR = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Allocate working storage.
C
      KT  = 1
      KTI = KT  + N*N
      KW  = KTI + N*N
C
C     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors Su and Ru
C     of the frequency-weighted controllability and observability
C     Grammians, respectively.
C
C     Workspace:   need 2*N*N + MAX(1, N*(N + MAX(N,M) + MIN(N,M) + 6)),
C                                                        if JOBCF = 'L';
C                       2*N*N + MAX(1, N*(N + MAX(N,P) + MIN(N,P) + 6)),
C                                                        if JOBCF = 'R'.
C                  prefer larger.
C
      CALL SB16CY( DICO, JOBCF, N, M, P, A, LDA, B, LDB, C, LDC,
     $             F, LDF, G, LDG, SCALEC, SCALEO, DWORK(KTI), N,
     $             DWORK(KT), N, DWORK(KW), LDWORK-KW+1, INFO )
C
      IF( INFO.NE.0 )
     $   RETURN
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Compute a B&T approximation (Ar,Br,Cr) of (A,B,C) and
C     the corresponding truncation matrices TI and T.
C
C     Real workspace:  need   2*N*N + MAX( 1, 2*N*N+5*N, N*MAX(M,P) );
C                      prefer larger.
C     Integer workspace:  0,  if JOBMR = 'B';
C                         N,  if JOBMR = 'F'.
C
      CALL AB09IX( DICO, JOBMR, 'NotSchur', ORDSEL, N, M, P, NCR,
     $             SCALEC, SCALEO, A, LDA, B, LDB, C, LDC, D, LDD,
     $             DWORK(KTI), N, DWORK(KT), N, NMR, HSV, TOL, TOL,
     $             IWORK, DWORK(KW), LDWORK-KW+1, IWARN, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 6
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Compute reduced gains Bc = Gr = TI*G and Cc = Fr = F*T.
C     Workspace:  need   N*(2*N+MAX(M,P)).
C
      CALL DLACPY( 'Full', N, P, G, LDG, DWORK(KW), N )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NCR, P, N, ONE,
     $            DWORK(KTI), N, DWORK(KW), N, ZERO, G, LDG )
C
      CALL DLACPY( 'Full', M, N, F, LDF, DWORK(KW), M )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', M, NCR, N, ONE,
     $            DWORK(KW), M, DWORK(KT), N, ZERO, F, LDF )
C
C     Form the reduced controller state matrix,
C     Ac = Ar + Br*Fr + Gr*Cr + Gr*D*Fr = Ar + Br*Fr + Gr*(Cr+D*Fr) .
C
C     Workspace:    need  P*N.
C
      CALL DLACPY( 'Full', P, NCR, C, LDC, DWORK, P )
      IF( WITHD) CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NCR, M,
     $                       ONE, D, LDD, F, LDF, ONE, DWORK, P )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NCR, NCR, P, ONE, G,
     $            LDG, DWORK, P, ONE, A, LDA )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NCR, NCR, M, ONE, B,
     $            LDB, F, LDF, ONE, A, LDA )
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB16CD ***
      END
