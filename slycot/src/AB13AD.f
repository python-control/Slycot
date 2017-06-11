      DOUBLE PRECISION FUNCTION AB13AD( DICO, EQUIL, N, M, P, ALPHA, A,
     $                                  LDA, B, LDB, C, LDC, NS, HSV,
     $                                  DWORK, LDWORK, INFO )
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
C     To compute the Hankel-norm of the ALPHA-stable projection of the
C     transfer-function matrix G of the state-space system (A,B,C).
C
C     FUNCTION VALUE
C
C     AB13AD  DOUBLE PRECISION
C             The Hankel-norm of the ALPHA-stable projection of G
C             (if INFO = 0).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state-space representation, i.e. the
C             order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the ALPHA-stability boundary for the eigenvalues
C             of the state dynamics matrix A. For a continuous-time
C             system (DICO = 'C'), ALPHA <= 0 is the boundary value for
C             the real parts of eigenvalues, while for a discrete-time
C             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the
C             boundary value for the moduli of eigenvalues.
C             The ALPHA-stability domain does not include the boundary
C             (see the Note below).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains the state dynamics matrix A in a block
C             diagonal real Schur form with its eigenvalues reordered
C             and separated. The resulting A has two diagonal blocks.
C             The leading NS-by-NS part of A has eigenvalues in the
C             ALPHA-stability domain and the trailing (N-NS) x (N-NS)
C             part has eigenvalues outside the ALPHA-stability domain.
C             Note: The ALPHA-stability domain is defined either
C                   as the open half complex plane left to ALPHA,
C                   for a continous-time system (DICO = 'C'), or the
C                   interior of the ALPHA-radius circle centered in the
C                   origin, for a discrete-time system (DICO = 'D').
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, if INFO = 0, the leading N-by-M part of this
C             array contains the input/state matrix B of the transformed
C             system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, if INFO = 0, the leading P-by-N part of this
C             array contains the state/output matrix C of the
C             transformed system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of HSV contain the
C             Hankel singular values of the ALPHA-stable part of the
C             original system ordered decreasingly.
C             HSV(1) is the Hankel norm of the ALPHA-stable subsystem.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*(MAX(N,M,P)+5)+N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the computation of the ordered real Schur form of A
C                   failed;
C             = 2:  the separation of the ALPHA-stable/unstable diagonal
C                   blocks failed because of very close eigenvalues;
C             = 3:  the computed ALPHA-stable part is just stable,
C                   having stable eigenvalues very near to the imaginary
C                   axis (if DICO = 'C') or to the unit circle
C                   (if DICO = 'D');
C             = 4:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the following linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t)                               (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system, and let G be the corresponding
C     transfer-function matrix. The following procedure is used to
C     compute the Hankel-norm of the ALPHA-stable projection of G:
C
C     1) Decompose additively G as
C
C          G = G1 + G2
C
C        such that G1 = (As,Bs,Cs) has only ALPHA-stable poles and
C        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles.
C        For the computation of the additive decomposition, the
C        algorithm presented in [1] is used.
C
C     2) Compute the Hankel-norm of ALPHA-stable projection G1 as the
C        the maximum Hankel singular value of the system (As,Bs,Cs).
C        The computation of the Hankel singular values is performed
C        by using the square-root method of [2].
C
C     REFERENCES
C
C     [1] Safonov, M.G., Jonckheere, E.A., Verma, M. and Limebeer, D.J.
C         Synthesis of positive real multivariable feedback systems,
C         Int. J. Control, Vol. 45, pp. 817-842, 1987.
C
C     [2] Tombs, M.S. and Postlethwaite, I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     NUMERICAL ASPECTS
C
C     The implemented method relies on a square-root technique.
C                                     3
C     The algorithms require about 17N  floating point operations.
C
C     CONTRIBUTOR
C
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine SHANRM.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest.
C     Jun. 2017, RvP, made 1st error return value zero
C
C     KEYWORDS
C
C     Additive spectral decomposition, model reduction,
C     multivariable system, state-space model, system norms.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  C100, ONE, ZERO
      PARAMETER         ( C100 = 100.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL
      INTEGER           INFO, LDA, LDB, LDC, LDWORK, M, N, NS, P
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           DISCR
      INTEGER           IERR, KT, KW, KW1, KW2
      DOUBLE PRECISION  ALPWRK, MAXRED, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  AB13AX, DLAMCH
      EXTERNAL          AB13AX, DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          TB01ID, TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO  = 0
      DISCR = LSAME( DICO, 'D' )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDWORK.LT.MAX( 1, N*( MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -16
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C     Error return.
C
C
         AB13AD = ZERO
         CALL XERBLA( 'AB13AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NS = 0
         AB13AD = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C        Workspace: N.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
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
C     Allocate working storage.
C
      KT  = 1
      KW1 = N*N + 1
      KW2 = KW1 + N
      KW  = KW2 + N
C
C     Reduce A to a block diagonal real Schur form, with the
C     ALPHA-stable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation A <- inv(T)*A*T and
C     apply the transformation to B and C: B <- inv(T)*B and C <- C*T.
C
C     Workspace needed:      N*(N+2);
C     Additional workspace:  need   3*N;
C                            prefer larger.
C
      CALL TB01KD( DICO, 'Stable', 'General', N, M, P, ALPWRK, A, LDA,
     $             B, LDB, C, LDC, NS, DWORK(KT), N, DWORK(KW1),
     $             DWORK(KW2), DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.NE.3 ) THEN
            INFO = 1
         ELSE
            INFO = 2
         END IF
         RETURN
      END IF
C
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
      IF( NS.EQ.0 ) THEN
         AB13AD = ZERO
      ELSE
C
C        Workspace:  need   N*(MAX(N,M,P)+5)+N*(N+1)/2;
C                    prefer larger.
C
         AB13AD = AB13AX( DICO, NS, M, P, A, LDA, B, LDB, C, LDC, HSV,
     $                    DWORK, LDWORK, IERR )
C
         IF( IERR.NE.0 ) THEN
            INFO = IERR + 2
            RETURN
         END IF
C
         DWORK(1) = MAX( WRKOPT, DWORK(1) )
      END IF
C
      RETURN
C *** Last line of AB13AD ***
      END
