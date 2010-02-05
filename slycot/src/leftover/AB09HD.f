      SUBROUTINE AB09HD( DICO, JOB, EQUIL, ORDSEL, N, M, P, NR, ALPHA,
     $                   BETA, A, LDA, B, LDB, C, LDC, D, LDD, NS, HSV,
     $                   TOL1, TOL2, IWORK, DWORK, LDWORK, BWORK, IWARN,
     $                   INFO )
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for an original
C     state-space representation (A,B,C,D) by using the stochastic
C     balancing approach in conjunction with the square-root or
C     the balancing-free square-root Balance & Truncate (B&T)
C     or Singular Perturbation Approximation (SPA) model reduction
C     methods for the ALPHA-stable part of the system.
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
C             = 'B':  use the square-root Balance & Truncate method;
C             = 'F':  use the balancing-free square-root
C                     Balance & Truncate method;
C             = 'S':  use the square-root Singular Perturbation
C                     Approximation method;
C             = 'P':  use the balancing-free square-root
C                     Singular Perturbation Approximation method.
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
C             The order of the original state-space representation,
C             i.e., the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C             P <= M if BETA = 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of the
C             resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. For a system with NU ALPHA-unstable
C             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N),
C             NR is set as follows: if ORDSEL = 'F', NR is equal to
C             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order
C             on entry, and NMIN is the order of a minimal realization
C             of the ALPHA-stable part of the given system; NMIN is
C             determined as the number of Hankel singular values greater
C             than NS*EPS, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH);
C             if ORDSEL = 'A', NR is the sum of NU and the number of
C             Hankel singular values greater than MAX(TOL1,NS*EPS);
C             NR can be further reduced to ensure that
C             HSV(NR-NU) > HSV(NR+1-NU).
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the ALPHA-stability boundary for the eigenvalues
C             of the state dynamics matrix A. For a continuous-time
C             system (DICO = 'C'), ALPHA <= 0 is the boundary value for
C             the real parts of eigenvalues, while for a discrete-time
C             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the
C             boundary value for the moduli of eigenvalues.
C             The ALPHA-stability domain does not include the boundary.
C
C     BETA    (input) DOUBLE PRECISION
C             BETA > 0 specifies the absolute/relative error weighting
C             parameter. A large positive value of BETA favours the
C             minimization of the absolute approximation error, while a
C             small value of BETA is appropriate for the minimization
C             of the relative error.
C             BETA = 0 means a pure relative error method and can be
C             used only if rank(D) = P.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the reduced
C             order system.
C             The resulting A has a block-diagonal form with two blocks.
C             For a system with NU ALPHA-unstable eigenvalues and
C             NS ALPHA-stable eigenvalues (NU+NS = N), the leading
C             NU-by-NU block contains the unreduced part of A
C             corresponding to ALPHA-unstable eigenvalues in an
C             upper real Schur form.
C             The trailing (NR+NS-N)-by-(NR+NS-N) block contains
C             the reduced part of A corresponding to ALPHA-stable
C             eigenvalues.
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
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of HSV contain the
C             Hankel singular values of the phase system corresponding
C             to the ALPHA-stable part of the original system.
C             The Hankel singular values are ordered decreasingly.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value of TOL1 lies
C             in the interval [0.00001,0.001].
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = NS*EPS, where NS is the number of
C             ALPHA-stable eigenvalues of A and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C             TOL1 < 1.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the phase system (see METHOD) corresponding
C             to the ALPHA-stable part of the given system.
C             The recommended value is TOL2 = NS*EPS.
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C             TOL2 < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(1,2*N)
C             On exit with INFO = 0, IWORK(1) contains the order of the
C             minimal realization of the system.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) contains RCOND, the reciprocal
C             condition number of the U11 matrix from the expression
C             used to compute the solution X = U21*inv(U11) of the
C             Riccati equation for spectral factorization.
C             A small value RCOND indicates possible ill-conditioning
C             of the respective Riccati equation.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 2*N*N + MB*(N+P) + MAX( 2, N*(MAX(N,MB,P)+5),
C                                    2*N*P+MAX(P*(MB+2),10*N*(N+1) ) ),
C             where MB = M if BETA = 0 and MB = M+P if BETA > 0.
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension 2*N
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than NSMIN, the sum of the order of the
C                   ALPHA-unstable part and the order of a minimal
C                   realization of the ALPHA-stable part of the given
C                   system; in this case, the resulting NR is set equal
C                   to NSMIN;
C             = 2:  with ORDSEL = 'F', the selected order NR corresponds
C                   to repeated singular values for the ALPHA-stable
C                   part, which are neither all included nor all
C                   excluded from the reduced model; in this case, the
C                   resulting NR is automatically decreased to exclude
C                   all repeated singular values;
C             = 3:  with ORDSEL = 'F', the selected order NR is less
C                   than the order of the ALPHA-unstable part of the
C                   given system; in this case NR is set equal to the
C                   order of the ALPHA-unstable part.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the computation of the ordered real Schur form of A
C                   failed;
C             = 2:  the reduction of the Hamiltonian matrix to real
C                   Schur form failed;
C             = 3:  the reordering of the real Schur form of the
C                   Hamiltonian matrix failed;
C             = 4:  the Hamiltonian matrix has less than N stable
C                   eigenvalues;
C             = 5:  the coefficient matrix U11 in the linear system
C                   X*U11 = U21 to determine X is singular to working
C                   precision;
C             = 6:  BETA = 0 and D has not a maximal row rank;
C             = 7:  the computation of Hankel singular values failed;
C             = 8:  the separation of the ALPHA-stable/unstable diagonal
C                   blocks failed because of very close eigenvalues;
C             = 9:  the resulting order of reduced stable part is less
C                   than the number of unstable zeros of the stable
C                   part.
C     METHOD
C
C     Let be the following linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t),                      (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09HD determines for
C     the given system (1), the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                  (2)
C
C     such that
C
C          INFNORM[inv(conj(W))*(G-Gr)] <=
C                       (1+HSV(NR+NS-N+1)) / (1-HSV(NR+NS-N+1)) + ...
C                       + (1+HSV(NS)) / (1-HSV(NS)) - 1,
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, W is the right, minimum
C     phase spectral factor satisfying
C
C         G1*conj(G1) = conj(W)* W,                      (3)
C
C     G1 is the NS-order ALPHA-stable part of G, and INFNORM(G) is the
C     infinity-norm of G. HSV(1), ... , HSV(NS) are the Hankel-singular
C     values of the stable part of the phase system (Ap,Bp,Cp)
C     with the transfer-function matrix
C
C          P = inv(conj(W))*G1.
C
C     If BETA > 0, then the model reduction is performed on [G BETA*I]
C     instead of G. This is the recommended approach to be used when D
C     has not a maximal row rank or when a certain balance between
C     relative and absolute approximation errors is desired. For
C     increasingly large values of BETA, the obtained reduced system
C     assymptotically approaches that computed by using the
C     Balance & Truncate or Singular Perturbation Approximation methods.
C
C     Note: conj(G)  denotes either G'(-s) for a continuous-time system
C           or G'(1/z) for a discrete-time system.
C           inv(G) is the inverse of G.
C
C     The following procedure is used to reduce a given G:
C
C     1) Decompose additively G as
C
C          G = G1 + G2,
C
C        such that G1 = (As,Bs,Cs,D) has only ALPHA-stable poles and
C        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles.
C
C     2) Determine G1r, a reduced order approximation of the
C        ALPHA-stable part G1 using the balancing stochastic method
C        in conjunction with either the B&T [1,2] or SPA methods [3].
C
C     3) Assemble the reduced model Gr as
C
C           Gr = G1r + G2.
C
C     Note: The employed stochastic truncation algorithm [2,3] has the
C     property that right half plane zeros of G1 remain as right half
C     plane zeros of G1r. Thus, the order can not be chosen smaller than
C     the sum of the number of unstable poles of G and the number of
C     unstable zeros of G1.
C
C     The reduction of the ALPHA-stable part G1 is done as follows.
C
C     If JOB = 'B', the square-root stochastic Balance & Truncate
C     method of [1] is used.
C     For an ALPHA-stable continuous-time system (DICO = 'C'),
C     the resulting reduced model is stochastically balanced.
C
C     If JOB = 'F', the balancing-free square-root version of the
C     stochastic Balance & Truncate method [1] is used to reduce
C     the ALPHA-stable part G1.
C
C     If JOB = 'S', the stochastic balancing method is used to reduce
C     the ALPHA-stable part G1, in conjunction with the square-root
C     version of the Singular Perturbation Approximation method [3,4].
C
C     If JOB = 'P', the stochastic balancing method is used to reduce
C     the ALPHA-stable part G1, in conjunction with the balancing-free
C     square-root version of the Singular Perturbation Approximation
C     method [3,4].
C
C     REFERENCES
C
C     [1] Varga A. and Fasol K.H.
C         A new square-root balancing-free stochastic truncation model
C         reduction algorithm.
C         Proc. 12th IFAC World Congress, Sydney, 1993.
C
C     [2] Safonov M. G. and Chiang R. Y.
C         Model reduction for robust control: a Schur relative error
C         method.
C         Int. J. Adapt. Contr. Sign. Proc., vol. 2, pp. 259-272, 1988.
C
C     [3] Green M. and Anderson B. D. O.
C         Generalized balanced stochastic truncation.
C         Proc. 29-th CDC, Honolulu, Hawaii, pp. 476-481, 1990.
C
C     [4] Varga A.
C         Balancing-free square-root algorithm for computing
C         singular perturbation approximations.
C         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
C         Vol. 2, pp. 1062-1065.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root or
C     balancing-free square-root techniques. The effectiveness of the
C     accuracy enhancing technique depends on the accuracy of the
C     solution of a Riccati equation. An ill-conditioned Riccati
C     solution typically results when [D BETA*I] is nearly
C     rank deficient.
C                                      3
C     The algorithm requires about 100N  floating point operations.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000.
C     D. Sima, University of Bucharest, May 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2000.
C     Partly based on the RASP routine SRBFS, by A. Varga, 1992.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000.
C              Oct. 2001.
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
      DOUBLE PRECISION  ZERO, ONE, TWO, TWOBY3, C100
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                    TWOBY3 = TWO/3.0D0, C100 = 100.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK,
     $                  M, N, NR, NS, P
      DOUBLE PRECISION  ALPHA, BETA, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*)
      LOGICAL           BWORK(*)
C     .. Local Scalars ..
      LOGICAL           BTA, DISCR, FIXORD, LEQUIL, SPA
      INTEGER           IERR, IWARNL, KB, KD, KT, KTI, KU, KW, KWI, KWR,
     $                  LW, LWR, MB, N2, NMR, NN, NRA, NU, NU1, WRKOPT
      DOUBLE PRECISION  EPSM, MAXRED, RICOND, SCALEC, SCALEO
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB04MD, AB09HY, AB09IX, DLACPY, DLASET, TB01ID,
     $                  TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LEQUIL = LSAME( EQUIL,  'S' )
      BTA    = LSAME( JOB,    'B' ) .OR. LSAME( JOB, 'F' )
      SPA    = LSAME( JOB,    'S' ) .OR. LSAME( JOB, 'P' )
      MB = M
      IF( BETA.GT.ZERO ) MB = M + P
      LW = 2*N*N + MB*(N+P) + MAX( 2, N*(MAX( N, MB, P )+5),
     $                 2*N*P+MAX( P*(MB+2), 10*N*(N+1) ) )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LEQUIL .OR. LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 .OR. ( BETA.EQ.ZERO .AND. P.GT.M ) ) THEN
         INFO = -7
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -8
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -9
      ELSE IF( BETA.LT.ZERO ) THEN
         INFO = -10
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -18
      ELSE IF( TOL1.GE.ONE ) THEN
         INFO = -21
      ELSE IF( ( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 )
     $         .OR. TOL2.GE.ONE ) THEN
         INFO = -22
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -25
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09HD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 .OR.
     $   ( BTA .AND. FIXORD .AND. NR.EQ.0 ) ) THEN
         NR = 0
         NS = 0
         IWORK(1) = 0
         DWORK(1) = TWO
         DWORK(2) = ONE
         RETURN
      END IF
C
      IF( LEQUIL ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C        Workspace: N.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
      END IF
C
C     Allocate working storage.
C
      NN  = N*N
      KU  = 1
      KWR = KU + NN
      KWI = KWR + N
      KW  = KWI + N
      LWR = LDWORK - KW + 1
C
C     Reduce A to a block-diagonal real Schur form, with the
C     ALPHA-unstable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation A <- inv(T)*A*T and
C     apply the transformation to B and C: B <- inv(T)*B and C <- C*T.
C
C     Workspace needed:      N*(N+2);
C     Additional workspace:  need   3*N;
C                            prefer larger.
C
      CALL TB01KD( DICO, 'Unstable', 'General', N, M, P, ALPHA, A, LDA,
     $             B, LDB, C, LDC, NU, DWORK(KU), N, DWORK(KWR),
     $             DWORK(KWI), DWORK(KW), LWR, IERR )
C
      IF( IERR.NE.0 ) THEN
         IF( IERR.NE.3 ) THEN
            INFO = 1
         ELSE
            INFO = 8
         END IF
         RETURN
      END IF
C
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
      IWARNL = 0
      NS = N - NU
      IF( FIXORD ) THEN
         NRA = MAX( 0, NR-NU )
         IF( NR.LT.NU )
     $      IWARNL = 3
      ELSE
         NRA = 0
      END IF
C
C     Finish if the system is completely unstable.
C
      IF( NS.EQ.0 ) THEN
         NR = NU
         IWORK(1) = NS
         DWORK(1) = WRKOPT
         DWORK(2) = ONE
         RETURN
      END IF
C
      NU1 = NU + 1
C
C     Allocate working storage.
C
      N2  = N + N
      KB  = 1
      KD  = KB  + N*MB
      KT  = KD  + P*MB
      KTI = KT  + N*N
      KW  = KTI + N*N
C
C     Form [B 0] and [D BETA*I].
C
      CALL DLACPY( 'F', NS, M, B(NU1,1), LDB, DWORK(KB), N )
      CALL DLACPY( 'F', P, M, D, LDD, DWORK(KD), P )
      IF( BETA.GT.ZERO ) THEN
         CALL DLASET( 'F', NS, P, ZERO, ZERO, DWORK(KB+N*M), N )
         CALL DLASET( 'F', P,  P, ZERO, BETA, DWORK(KD+P*M), P )
      END IF
C
C     For discrete-time case, apply the discrete-to-continuous bilinear
C     transformation to the stable part.
C
      IF( DISCR ) THEN
C
C        Real workspace:    need  N, prefer larger;
C        Integer workspace: need  N.
C
         CALL AB04MD( 'Discrete', NS, MB, P, ONE, ONE, A(NU1,NU1), LDA,
     $                DWORK(KB), N, C(1,NU1), LDC, DWORK(KD), P,
     $                IWORK, DWORK(KT), LDWORK-KT+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(KT) ) + KT - 1 )
      END IF
C
C     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R
C     of the controllability and observability Grammians, respectively.
C     Real workspace:    need  2*N*N + MB*(N+P)+
C                              MAX( 2, N*(MAX(N,MB,P)+5),
C                                   2*N*P+MAX(P*(MB+2), 10*N*(N+1) ) );
C                        prefer larger.
C     Integer workspace: need  2*N.
C
      CALL AB09HY( NS, MB, P, A(NU1,NU1), LDA, DWORK(KB), N,
     $             C(1,NU1), LDC, DWORK(KD), P, SCALEC, SCALEO,
     $             DWORK(KTI), N, DWORK(KT), N, IWORK, DWORK(KW),
     $             LDWORK-KW+1, BWORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      RICOND = DWORK(KW+1)
C
C     Compute a BTA or SPA of the stable part.
C     Real workspace:  need  2*N*N + MB*(N+P)+
C                            MAX( 1, 2*N*N+5*N, N*MAX(MB,P) ).
C
      EPSM = DLAMCH( 'Epsilon' )
      CALL AB09IX( 'C', JOB, 'Schur', ORDSEL, NS, MB, P, NRA, SCALEC,
     $             SCALEO, A(NU1,NU1), LDA, DWORK(KB), N, C(1,NU1), LDC,
     $             DWORK(KD), P, DWORK(KTI), N, DWORK(KT), N, NMR, HSV,
     $             MAX( TOL1, N*EPSM ), TOL2, IWORK, DWORK(KW),
     $             LDWORK-KW+1, IWARN, IERR )
      IWARN = MAX( IWARN, IWARNL )
      IF( IERR.NE.0 ) THEN
         INFO = 7
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Check if the resulting order is greater than the number of
C     unstable zeros (this check is implicit by looking at Hankel
C     singular values equal to 1).
C
      IF( NRA.LT.NS .AND. HSV(NRA+1).GE.ONE-EPSM**TWOBY3 ) THEN
         INFO = 9
         RETURN
      END IF
C
C     For discrete-time case, apply the continuous-to-discrete
C     bilinear transformation.
C
      IF( DISCR ) THEN
         CALL AB04MD( 'Continuous', NRA, MB, P, ONE, ONE,
     $                A(NU1,NU1), LDA, DWORK(KB), N, C(1,NU1), LDC,
     $                DWORK(KD), P, IWORK, DWORK, LDWORK, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      END IF
C
      CALL DLACPY( 'F', NRA, M, DWORK(KB), N, B(NU1,1), LDB )
      CALL DLACPY( 'F', P, M, DWORK(KD), P, D, LDD )
C
      NR = NRA + NU
C
      IWORK(1) = NMR
      DWORK(1) = WRKOPT
      DWORK(2) = RICOND
C
      RETURN
C *** Last line of AB09HD ***
      END
