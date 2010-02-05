      SUBROUTINE AB09BD( DICO, JOB, EQUIL, ORDSEL, N, M, P, NR, A, LDA,
     $                   B, LDB, C, LDC, D, LDD, HSV, TOL1, TOL2, IWORK,
     $                   DWORK, LDWORK, IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
C     original state-space representation (A,B,C,D) by using either the
C     square-root or the balancing-free square-root Singular
C     Perturbation Approximation (SPA) model reduction method.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOB     CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root SPA method;
C             = 'N':  use the balancing-free square-root SPA method.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation, i.e.
C             the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of
C             the resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. NR is set as follows:
C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
C             is the desired order on entry and NMIN is the order of a
C             minimal realization of the given system; NMIN is
C             determined as the number of Hankel singular values greater
C             than N*EPS*HNORM(A,B,C), where EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH) and
C             HNORM(A,B,C) is the Hankel norm of the system (computed
C             in HSV(1));
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the
C             reduced order system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, if INFO = 0, the leading NR-by-M part of this
C             array contains the input/state matrix Br of the reduced
C             order system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, if INFO = 0, the leading P-by-NR part of this
C             array contains the state/output matrix Cr of the reduced
C             order system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the original input/output matrix D.
C             On exit, if INFO = 0, the leading P-by-M part of this
C             array contains the input/output matrix Dr of the reduced
C             order system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the Hankel singular values of
C             the original system ordered decreasingly. HSV(1) is the
C             Hankel norm of the system.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(A,B,C), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(A,B,C) is the
C             Hankel-norm of the given system (computed in HSV(1)).
C             For computing a minimal realization, the recommended
C             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH).
C             This value is used by default if TOL1 <= 0 on entry.
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the given system. The recommended value is
C             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default
C             if TOL2 <= 0 on entry.
C             If TOL2 > 0, then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(1,2*N)
C             On exit with INFO = 0, IWORK(1) contains the order of the
C             minimal realization of the system.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than the order of a minimal realization of the
C                   given system. In this case, the resulting NR is
C                   set automatically to a value corresponding to the
C                   order of a minimal realization of the system.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A to the real Schur form failed;
C             = 2:  the state matrix A is not stable (if DICO = 'C')
C                   or not convergent (if DICO = 'D');
C             = 3:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t)                           (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09BD determines for
C     the given system (1), the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t)                       (2)
C
C     such that
C
C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     If JOB = 'B', the balancing-based square-root SPA method of [1]
C     is used and the resulting model is balanced.
C
C     If JOB = 'N', the balancing-free square-root SPA method of [2]
C     is used.
C     By setting TOL1 = TOL2, the routine can be used to compute
C     Balance & Truncate approximations.
C
C     REFERENCES
C
C     [1] Liu Y. and Anderson B.D.O.
C         Singular Perturbation Approximation of Balanced Systems,
C         Int. J. Control, Vol. 50, pp. 1379-1405, 1989.
C
C     [2] Varga A.
C         Balancing-free square-root algorithm for computing singular
C         perturbation approximations.
C         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
C         Vol. 2, pp. 1062-1065.
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
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine SRBFSP.
C
C     REVISIONS
C
C     May 2, 1998.
C     November 11, 1998, V. Sima, Research Institute for Informatics,
C     Bucharest.
C
C     KEYWORDS
C
C     Balancing, minimal state-space representation, model reduction,
C     multivariable system, singular perturbation approximation,
C     state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, C100
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, C100 = 100.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK,
     $                  M, N, NR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           FIXORD
      INTEGER           IERR, KI, KR, KT, KTI, KW, NN
      DOUBLE PRECISION  MAXRED, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09BX, TB01ID, TB01WD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      FIXORD = LSAME( ORDSEL, 'F' )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. LSAME( DICO, 'D' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( JOB, 'B' ) .OR. LSAME( JOB, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( TOL2.GT.ZERO .AND. TOL2.GT.TOL1 ) THEN
         INFO = -19
      ELSE IF( LDWORK.LT.MAX( 1, N*( 2*N + MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -22
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09BD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Allocate working storage.
C
      NN = N*N
      KT = 1
      KR = KT + NN
      KI = KR + N
      KW = KI + N
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
      END IF
C
C     Reduce A to the real Schur form using an orthogonal similarity
C     transformation A <- T'*A*T and apply the transformation to
C     B and C: B <- T'*B and C <- C*T.
C
      CALL TB01WD( N, M, P, A, LDA, B, LDB, C, LDC, DWORK(KT), N,
     $             DWORK(KR), DWORK(KI), DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
      KTI = KT  + NN
      KW  = KTI + NN
      CALL AB09BX( DICO, JOB, ORDSEL, N, M, P, NR, A, LDA, B, LDB,
     $             C, LDC, D, LDD, HSV, DWORK(KT), N, DWORK(KTI), N,
     $             TOL1, TOL2, IWORK, DWORK(KW), LDWORK-KW+1, IWARN,
     $             IERR )
C
      IF( IERR.NE.0 ) THEN
         INFO = IERR + 1
         RETURN
      END IF
C
      DWORK(1) = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
      RETURN
C *** Last line of AB09BD ***
      END
