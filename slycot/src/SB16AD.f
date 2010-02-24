      SUBROUTINE SB16AD( DICO, JOBC, JOBO, JOBMR, WEIGHT, EQUIL, ORDSEL,
     $                   N, M, P, NC, NCR, ALPHA, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, AC, LDAC, BC, LDBC, CC, LDCC,
     $                   DC, LDDC, NCS, HSVC, TOL1, TOL2, IWORK, DWORK,
     $                   LDWORK, IWARN, INFO )
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
C     To compute a reduced order controller (Acr,Bcr,Ccr,Dcr) for an
C     original state-space controller representation (Ac,Bc,Cc,Dc) by
C     using the frequency-weighted square-root or balancing-free
C     square-root Balance & Truncate (B&T) or Singular Perturbation
C     Approximation (SPA) model reduction methods. The algorithm tries
C     to minimize the norm of the frequency-weighted error
C
C           ||V*(K-Kr)*W||
C
C     where K and Kr are the transfer-function matrices of the original
C     and reduced order controllers, respectively. V and W are special
C     frequency-weighting transfer-function matrices constructed
C     to enforce closed-loop stability and/or closed-loop performance.
C     If G is the transfer-function matrix of the open-loop system, then
C     the following weightings V and W can be used:
C                      -1
C      (a)   V = (I-G*K) *G, W = I - to enforce closed-loop stability;
C                              -1
C      (b)   V = I,  W = (I-G*K) *G - to enforce closed-loop stability;
C                      -1              -1
C      (c)   V = (I-G*K) *G, W = (I-G*K)  - to enforce closed-loop
C            stability and performance.
C
C     G has the state space representation (A,B,C,D).
C     If K is unstable, only the ALPHA-stable part of K is reduced.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original controller as follows:
C             = 'C':  continuous-time controller;
C             = 'D':  discrete-time controller.
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
C     JOBMR   CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root B&T method;
C             = 'F':  use the balancing-free square-root B&T method;
C             = 'S':  use the square-root SPA method;
C             = 'P':  use the balancing-free square-root SPA method.
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
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplets (A,B,C) and (Ac,Bc,Cc) as
C             follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NCR is fixed;
C             = 'A':  the resulting order NCR is automatically
C                     determined on basis of the given tolerance TOL1.
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
C     NCR     (input/output) INTEGER
C             On entry with ORDSEL = 'F', NCR is the desired order of
C             the resulting reduced order controller.  0 <= NCR <= NC.
C             On exit, if INFO = 0, NCR is the order of the resulting
C             reduced order controller. For a controller with NCU
C             ALPHA-unstable eigenvalues and NCS ALPHA-stable
C             eigenvalues (NCU+NCS = NC), NCR is set as follows:
C             if ORDSEL = 'F', NCR is equal to
C             NCU+MIN(MAX(0,NCR-NCU),NCMIN), where NCR is the desired
C             order on entry, NCMIN is the number of frequency-weighted
C             Hankel singular values greater than NCS*EPS*S1, EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH) and
C             S1 is the largest Hankel singular value (computed in
C             HSVC(1)); NCR can be further reduced to ensure
C             HSVC(NCR-NCU) > HSVC(NCR+1-NCU);
C             if ORDSEL = 'A', NCR is the sum of NCU and the number of
C             Hankel singular values greater than MAX(TOL1,NCS*EPS*S1).
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the ALPHA-stability boundary for the eigenvalues
C             of the state dynamics matrix AC. For a continuous-time
C             controller (DICO = 'C'), ALPHA <= 0 is the boundary value
C             for the real parts of eigenvalues; for a discrete-time
C             controller (DICO = 'D'), 0 <= ALPHA <= 1 represents the
C             boundary value for the moduli of eigenvalues.
C             The ALPHA-stability domain does not include the boundary.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A of the open-loop
C             system.
C             On exit, if INFO = 0 and EQUIL = 'S', the leading N-by-N
C             part of this array contains the scaled state dynamics
C             matrix of the open-loop system.
C             If EQUIL = 'N', this array is unchanged on exit.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix B of the open-loop system.
C             On exit, if INFO = 0 and EQUIL = 'S', the leading N-by-M
C             part of this array contains the scaled input/state matrix
C             of the open-loop system.
C             If EQUIL = 'N', this array is unchanged on exit.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C of the open-loop system.
C             On exit, if INFO = 0 and EQUIL = 'S', the leading P-by-N
C             part of this array contains the scaled state/output matrix
C             of the open-loop system.
C             If EQUIL = 'N', this array is unchanged on exit.
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
C     AC      (input/output) DOUBLE PRECISION array, dimension (LDAC,NC)
C             On entry, the leading NC-by-NC part of this array must
C             contain the state dynamics matrix Ac of the original
C             controller.
C             On exit, if INFO = 0, the leading NCR-by-NCR part of this
C             array contains the state dynamics matrix Acr of the
C             reduced controller. The resulting Ac has a
C             block-diagonal form with two blocks.
C             For a system with NCU ALPHA-unstable eigenvalues and
C             NCS ALPHA-stable eigenvalues (NCU+NCS = NC), the leading
C             NCU-by-NCU block contains the unreduced part of Ac
C             corresponding to the ALPHA-unstable eigenvalues.
C             The trailing (NCR+NCS-NC)-by-(NCR+NCS-NC) block contains
C             the reduced part of Ac corresponding to ALPHA-stable
C             eigenvalues.
C
C     LDAC    INTEGER
C             The leading dimension of array AC.  LDAC >= MAX(1,NC).
C
C     BC      (input/output) DOUBLE PRECISION array, dimension (LDBC,P)
C             On entry, the leading NC-by-P part of this array must
C             contain the input/state matrix Bc of the original
C             controller.
C             On exit, if INFO = 0, the leading NCR-by-P part of this
C             array contains the input/state matrix Bcr of the reduced
C             controller.
C
C     LDBC    INTEGER
C             The leading dimension of array BC.  LDBC >= MAX(1,NC).
C
C     CC      (input/output) DOUBLE PRECISION array, dimension (LDCC,NC)
C             On entry, the leading M-by-NC part of this array must
C             contain the state/output matrix Cc of the original
C             controller.
C             On exit, if INFO = 0, the leading M-by-NCR part of this
C             array contains the state/output matrix Ccr of the reduced
C             controller.
C
C     LDCC    INTEGER
C             The leading dimension of array CC.  LDCC >= MAX(1,M).
C
C     DC      (input/output) DOUBLE PRECISION array, dimension (LDDC,P)
C             On entry, the leading M-by-P part of this array must
C             contain the input/output matrix Dc of the original
C             controller.
C             On exit, if INFO = 0, the leading M-by-P part of this
C             array contains the input/output matrix Dcr of the reduced
C             controller.
C
C     LDDC    INTEGER
C             The leading dimension of array DC.  LDDC >= MAX(1,M).
C
C     NCS     (output) INTEGER
C             The dimension of the ALPHA-stable part of the controller.
C
C     HSVC    (output) DOUBLE PRECISION array, dimension (NC)
C             If INFO = 0, the leading NCS elements of this array
C             contain the frequency-weighted Hankel singular values,
C             ordered decreasingly, of the ALPHA-stable part of the
C             controller.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of the reduced controller.
C             For model reduction, the recommended value is
C             TOL1 = c*S1, where c is a constant in the
C             interval [0.00001,0.001], and S1 is the largest
C             frequency-weighted Hankel singular value of the
C             ALPHA-stable part of the original controller
C             (computed in HSVC(1)).
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = NCS*EPS*S1, where NCS is the number of
C             ALPHA-stable eigenvalues of Ac and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the ALPHA-stable part of the given
C             controller. The recommended value is TOL2 = NCS*EPS*S1.
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(1,LIWRK1,LIWRK2)
C             LIWRK1 = 0,       if JOBMR  = 'B';
C             LIWRK1 = NC,      if JOBMR  = 'F';
C             LIWRK1 = 2*NC,    if JOBMR  = 'S' or 'P';
C             LIWRK2 = 0,       if WEIGHT = 'N';
C             LIWRK2 = 2*(M+P), if WEIGHT = 'O', 'I', or 'P'.
C             On exit, if INFO = 0, IWORK(1) contains NCMIN, the order
C             of the computed minimal realization of the stable part of
C             the controller.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 2*NC*NC + MAX( 1, LFREQ, LSQRED ),
C             where
C             LFREQ = (N+NC)*(N+NC+2*M+2*P)+
C                     MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7), (M+P)*(M+P+4))
C                                      if WEIGHT = 'I' or 'O' or 'P';
C             LFREQ  = NC*(MAX(M,P)+5) if WEIGHT = 'N' and EQUIL = 'N';
C             LFREQ  = MAX(N,NC*(MAX(M,P)+5)) if WEIGHT = 'N' and
C                                                EQUIL  = 'S';
C             LSQRED = MAX( 1, 2*NC*NC+5*NC );
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NCR is greater
C                   than NSMIN, the sum of the order of the
C                   ALPHA-unstable part and the order of a minimal
C                   realization of the ALPHA-stable part of the given
C                   controller; in this case, the resulting NCR is set
C                   equal to NSMIN;
C             = 2:  with ORDSEL = 'F', the selected order NCR
C                   corresponds to repeated singular values for the
C                   ALPHA-stable part of the controller, which are
C                   neither all included nor all excluded from the
C                   reduced model; in this case, the resulting NCR is
C                   automatically decreased to exclude all repeated
C                   singular values;
C             = 3:  with ORDSEL = 'F', the selected order NCR is less
C                   than the order of the ALPHA-unstable part of the
C                   given controller. In this case NCR is set equal to
C                   the order of the ALPHA-unstable part.
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
C             = 5:  the computation of the ordered real Schur form of Ac
C                   failed;
C             = 6:  the separation of the ALPHA-stable/unstable
C                   diagonal blocks failed because of very close
C                   eigenvalues;
C             = 7:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let K be the transfer-function matrix of the original linear
C     controller
C
C          d[xc(t)] = Ac*xc(t) + Bc*y(t)
C          u(t)     = Cc*xc(t) + Dc*y(t),                      (1)
C
C     where d[xc(t)] is dxc(t)/dt for a continuous-time system and
C     xc(t+1) for a discrete-time system. The subroutine SB16AD
C     determines the matrices of a reduced order controller
C
C          d[z(t)] = Acr*z(t) + Bcr*y(t)
C          u(t)    = Ccr*z(t) + Dcr*y(t),                      (2)
C
C     such that the corresponding transfer-function matrix Kr minimizes
C     the norm of the frequency-weighted error
C
C             V*(K-Kr)*W,                                      (3)
C
C     where V and W are special stable transfer-function matrices
C     chosen to enforce stability and/or performance of the closed-loop
C     system [3] (see description of the parameter WEIGHT).
C
C     The following procedure is used to reduce K in conjunction
C     with the frequency-weighted balancing approach of [2]
C     (see also [3]):
C
C     1) Decompose additively K, of order NC, as
C
C          K = K1 + K2,
C
C        such that K1 has only ALPHA-stable poles and K2, of order NCU,
C        has only ALPHA-unstable poles.
C
C     2) Compute for K1 a B&T or SPA frequency-weighted approximation
C        K1r of order NCR-NCU using the frequency-weighted balancing
C        approach of [1] in conjunction with accuracy enhancing
C        techniques specified by the parameter JOBMR.
C
C     3) Assemble the reduced model Kr as
C
C           Kr = K1r + K2.
C
C     For the reduction of the ALPHA-stable part, several accuracy
C     enhancing techniques can be employed (see [2] for details).
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
C     For each of these methods, two left and right truncation matrices
C     are determined using the Cholesky factors of an input
C     frequency-weighted controllability Grammian P and an output
C     frequency-weighted observability Grammian Q.
C     P and Q are determined as the leading NC-by-NC diagonal blocks
C     of the controllability Grammian of K*W and of the
C     observability Grammian of V*K. Special techniques developed in [2]
C     are used to compute the Cholesky factors of P and Q directly
C     (see also SLICOT Library routine SB16AY).
C     The frequency-weighted Hankel singular values HSVC(1), ....,
C     HSVC(NC) are computed as the square roots of the eigenvalues
C     of the product P*Q.
C
C     REFERENCES
C
C     [1] Enns, D.
C         Model reduction with balanced realizations: An error bound
C         and a frequency weighted generalization.
C         Proc. 23-th CDC, Las Vegas, pp. 127-132, 1984.
C
C     [2] Varga, A. and Anderson, B.D.O.
C         Square-root balancing-free methods for frequency-weighted
C         balancing related model reduction.
C         (report in preparation)
C
C     [3] Anderson, B.D.O and Liu, Y.
C         Controller reduction: concepts and approaches.
C         IEEE Trans. Autom. Control, Vol. 34, pp. 802-812, 1989.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root
C     techniques.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, Sept. 2000.
C     D. Sima, University of Bucharest, Sept. 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Sept.2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C              Sep. 2001.
C
C     KEYWORDS
C
C     Controller reduction, frequency weighting, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  C100, ONE, ZERO
      PARAMETER         ( C100 = 100.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOBC, JOBO, JOBMR, ORDSEL, WEIGHT
      INTEGER           INFO, IWARN, LDA, LDAC, LDB, LDBC, LDC, LDCC,
     $                  LDD, LDDC, LDWORK, M, N, NC, NCR, NCS, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AC(LDAC,*), B(LDB,*), BC(LDBC,*),
     $                  C(LDC,*), CC(LDCC,*), D(LDD,*), DC(LDDC,*),
     $                  DWORK(*), HSVC(*)
C     .. Local Scalars ..
      LOGICAL           BAL, BTA, DISCR, FIXORD, FRWGHT, ISTAB, LEFTW,
     $                  OSTAB, PERF, RIGHTW, SPA
      INTEGER           IERR, IWARNL, KI, KR, KT, KTI, KU, KW, LW, MP,
     $                  NCU, NCU1, NMR, NNC, NRA, WRKOPT
      DOUBLE PRECISION  ALPWRK, MAXRED, SCALEC, SCALEO
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09IX, SB16AY, TB01ID, TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      BTA    = LSAME( JOBMR,  'B' ) .OR. LSAME( JOBMR, 'F' )
      SPA    = LSAME( JOBMR,  'S' ) .OR. LSAME( JOBMR, 'P' )
      BAL    = LSAME( JOBMR,  'B' ) .OR. LSAME( JOBMR, 'S' )
      FIXORD = LSAME( ORDSEL, 'F' )
      ISTAB  = LSAME( WEIGHT, 'I' )
      OSTAB  = LSAME( WEIGHT, 'O' )
      PERF   = LSAME( WEIGHT, 'P' )
      LEFTW  = OSTAB .OR. PERF
      RIGHTW = ISTAB .OR. PERF
      FRWGHT = LEFTW .OR. RIGHTW
C
      LW  = 1
      NNC = N + NC
      MP  = M + P
      IF( FRWGHT ) THEN
         LW = NNC*( NNC + 2*MP ) +
     $        MAX( NNC*( NNC + MAX( NNC, M, P ) + 7 ), MP*( MP + 4 ) )
      ELSE
         LW = NC*( MAX( M, P ) + 5 )
         IF ( LSAME( EQUIL, 'S' ) )
     $      LW = MAX( N, LW )
      END IF
      LW = 2*NC*NC + MAX( 1, LW, NC*( 2*NC + 5 ) )
C
C     Check the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( JOBC, 'S' ) .OR. LSAME( JOBC, 'E' ) ) )
     $     THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( JOBO, 'S' ) .OR. LSAME( JOBO, 'E' ) ) )
     $     THEN
         INFO = -3
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -4
      ELSE IF( .NOT.( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -6
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -7
      ELSE IF( N.LT.0 ) THEN
         INFO = -8
      ELSE IF( M.LT.0 ) THEN
         INFO = -9
      ELSE IF( P.LT.0 ) THEN
         INFO = -10
      ELSE IF( NC.LT.0 ) THEN
         INFO = -11
      ELSE IF( FIXORD .AND. ( NCR.LT.0 .OR. NCR.GT.NC ) ) THEN
         INFO = -12
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -13
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -19
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -21
      ELSE IF( LDAC.LT.MAX( 1, NC ) ) THEN
         INFO = -23
      ELSE IF( LDBC.LT.MAX( 1, NC ) ) THEN
         INFO = -25
      ELSE IF( LDCC.LT.MAX( 1, M  ) ) THEN
         INFO = -27
      ELSE IF( LDDC.LT.MAX( 1, M  ) ) THEN
         INFO = -29
      ELSE IF( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 ) THEN
         INFO = -33
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -36
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB16AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( NC, M, P ).EQ.0 ) THEN
         NCR = 0
         NCS = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C and AC, BC and CC;
C        A <- inv(T1)*A*T1, B <- inv(T1)*B and C <- C*T1, where T1 is a
C        diagonal matrix;
C        AC <- inv(T2)*AC*T2, BC <- inv(T2)*BC and CC <- CC*T2, where T2
C        is a diagonal matrix.
C
C        Real workspace: need MAX(N,NC).
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
         MAXRED = C100
         CALL TB01ID( 'All', NC, P, M, MAXRED, AC, LDAC, BC, LDBC,
     $                CC, LDCC, DWORK, INFO )
      END IF
C
C     Correct the value of ALPHA to ensure stability.
C
      ALPWRK = ALPHA
      IF( DISCR ) THEN
         IF( ALPHA.EQ.ONE ) ALPWRK = ONE - SQRT( DLAMCH( 'E' ) )
      ELSE
         IF( ALPHA.EQ.ZERO ) ALPWRK = -SQRT( DLAMCH( 'E' ) )
      END IF
C
C     Reduce Ac to a block-diagonal real Schur form, with the
C     ALPHA-unstable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation, AC <- inv(T)*AC*T, and
C     apply the transformation to BC and CC:
C     BC <- inv(T)*BC and CC <- CC*T.
C
C     Workspace:  need   NC*(NC+5);
C                 prefer larger.
C
      WRKOPT = 1
      KU = 1
      KR = KU + NC*NC
      KI = KR + NC
      KW = KI + NC
C
      CALL TB01KD( DICO, 'Unstable', 'General', NC, P, M, ALPWRK,
     $             AC, LDAC, BC, LDBC, CC, LDCC, NCU, DWORK(KU), NC,
     $             DWORK(KR), DWORK(KI), DWORK(KW), LDWORK-KW+1, IERR )
C
      IF( IERR.NE.0 ) THEN
         IF( IERR.NE.3 ) THEN
            INFO = 5
         ELSE
            INFO = 6
         END IF
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
      IWARNL = 0
      NCS = NC - NCU
      IF( FIXORD ) THEN
         NRA = MAX( 0, NCR-NCU )
         IF( NCR.LT.NCU )
     $      IWARNL = 3
      ELSE
         NRA = 0
      END IF
C
C     Finish if only unstable part is present.
C
      IF( NCS.EQ.0 ) THEN
         NCR = NCU
         IWORK(1) = 0
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Allocate working storage.
C
      KT  = 1
      KTI = KT  + NC*NC
      KW  = KTI + NC*NC
C
C     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R
C     of the frequency-weighted controllability and observability
C     Grammians, respectively.
C
C     Real workspace:  need  2*NC*NC + MAX( 1, LFREQ ),
C                      where
C                      LFREQ = (N+NC)*(N+NC+2*M+2*P)+
C                              MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7),
C                                  (M+P)*(M+P+4))
C                                         if WEIGHT = 'I' or 'O' or 'P';
C                      LFREQ = NCS*(MAX(M,P)+5) if WEIGHT = 'N';
C                      prefer larger.
C     Integer workspace:      2*(M+P) if WEIGHT = 'I' or 'O' or 'P';
C                             0,      if WEIGHT = 'N'.
C
      CALL SB16AY( DICO, JOBC, JOBO, WEIGHT, N, M, P, NC, NCS,
     $             A, LDA, B, LDB, C, LDC, D, LDD,
     $             AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC,
     $             SCALEC, SCALEO, DWORK(KTI), NC, DWORK(KT), NC,
     $             IWORK, DWORK(KW), LDWORK-KW+1, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Compute a BTA or SPA of the stable part.
C     Real workspace:  need   2*NC*NC + MAX( 1, 2*NC*NC+5*NC,
C                                               NC*MAX(M,P) );
C                      prefer larger.
C     Integer workspace:      0,     if JOBMR = 'B';
C                             NC,    if JOBMR = 'F';
C                             2*NC,  if JOBMR = 'S' or 'P'.
C
      NCU1 = NCU + 1
      CALL AB09IX( DICO, JOBMR, 'Schur', ORDSEL, NCS, P, M, NRA, SCALEC,
     $             SCALEO, AC(NCU1,NCU1), LDAC, BC(NCU1,1), LDBC,
     $             CC(1,NCU1), LDCC, DC, LDDC, DWORK(KTI), NC,
     $             DWORK(KT), NC, NMR, HSVC, TOL1, TOL2, IWORK,
     $             DWORK(KW), LDWORK-KW+1, IWARN, IERR )
      IWARN = MAX( IWARN, IWARNL )
      IF( IERR.NE.0 ) THEN
         INFO = 7
         RETURN
      END IF
      NCR = NRA + NCU
      IWORK(1) = NMR
C
      DWORK(1) = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
      RETURN
C *** Last line of SB16AD ***
      END
