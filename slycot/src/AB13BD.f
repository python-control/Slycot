      DOUBLE PRECISION FUNCTION AB13BD( DICO, JOBN, N, M, P, A, LDA,
     $                                  B, LDB, C, LDC, D, LDD, NQ, TOL,
     $                                  DWORK, LDWORK, IWARN, INFO)
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
C     To compute the H2 or L2 norm of the transfer-function matrix G
C     of the system (A,B,C,D). G must not have poles on the imaginary
C     axis, for a continuous-time system, or on the unit circle, for
C     a discrete-time system. If the H2-norm is computed, the system
C     must be stable.
C
C     FUNCTION VALUE
C
C     AB13BD   DOUBLE PRECISION
C              The H2-norm of G, if JOBN = 'H', or the L2-norm of G,
C              if JOBN = 'L' (if INFO = 0).
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
C     JOBN    CHARACTER*1
C             Specifies the norm to be computed as follows:
C             = 'H':  the H2-norm;
C             = 'L':  the L2-norm.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A, the number of rows of the
C             matrix B, and the number of columns of the matrix C.
C             N represents the dimension of the state vector.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrices B and D.
C             M represents the dimension of input vector.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of the matrices C and D.
C             P represents the dimension of output vector.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix of the system.
C             On exit, the leading NQ-by-NQ part of this array contains
C             the state dynamics matrix (in a real Schur form) of the
C             numerator factor Q of the right coprime factorization with
C             inner denominator of G (see METHOD).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix of the system.
C             On exit, the leading NQ-by-M part of this array contains
C             the input/state matrix of the numerator factor Q of the
C             right coprime factorization with inner denominator of G
C             (see METHOD).
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix of the system.
C             On exit, the leading P-by-NQ part of this array contains
C             the state/output matrix of the numerator factor Q of the
C             right coprime factorization with inner denominator of G
C             (see METHOD).
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the input/output matrix of the system.
C             If DICO = 'C', D must be a null matrix.
C             On exit, the leading P-by-M part of this array contains
C             the input/output matrix of the numerator factor Q of
C             the right coprime factorization with inner denominator
C             of G (see METHOD).
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     NQ      (output) INTEGER
C             The order of the resulting numerator Q of the right
C             coprime factorization with inner denominator of G (see
C             METHOD).
C             Generally, NQ = N - NS, where NS is the number of
C             uncontrollable unstable eigenvalues.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of
C             B are considered zero (used for controllability tests).
C             If the user sets TOL <= 0, then an implicitly computed,
C             default tolerance, defined by  TOLDEF = N*EPS*NORM(B),
C             is used instead, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH) and NORM(B) denotes
C             the 1-norm of B.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.
C             LDWORK >= MAX( 1, M*(N+M) + MAX( N*(N+5), M*(M+2), 4*P ),
C                               N*( MAX( N, P ) + 4 ) + MIN( N, P ) ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = K:  K violations of the numerical stability condition
C                   occured during the assignment of eigenvalues in
C                   computing the right coprime factorization with inner
C                   denominator of G (see the SLICOT subroutine SB08DD).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A to a real Schur form failed;
C             = 2:  a failure was detected during the reordering of the
C                   real Schur form of A, or in the iterative process
C                   for reordering the eigenvalues of Z'*(A + B*F)*Z
C                   along the diagonal (see SLICOT routine SB08DD);
C             = 3:  if DICO = 'C' and the matrix A has a controllable
C                   eigenvalue on the imaginary axis, or DICO = 'D'
C                   and A has a controllable eigenvalue on the unit
C                   circle;
C             = 4:  the solution of Lyapunov equation failed because
C                   the equation is singular;
C             = 5:  if DICO = 'C' and D is a nonzero matrix;
C             = 6:  if JOBN = 'H' and the system is unstable.
C
C     METHOD
C
C     The subroutine is based on the algorithms proposed in [1] and [2].
C
C     If the given transfer-function matrix G is unstable, then a right
C     coprime factorization with inner denominator of G is first
C     computed
C               -1
C        G = Q*R  ,
C
C     where Q and R are stable transfer-function matrices and R is
C     inner. If G is stable, then Q = G and R = I.
C     Let (AQ,BQ,CQ,DQ) be the state-space representation of Q.
C
C     If DICO = 'C', then the L2-norm of G is computed as
C
C        NORM2(G) = NORM2(Q) = SQRT(TRACE(BQ'*X*BQ)),
C
C     where X satisfies the continuous-time Lyapunov equation
C
C        AQ'*X + X*AQ + CQ'*CQ = 0.
C
C     If DICO = 'D', then the l2-norm of G is computed as
C
C        NORM2(G) = NORM2(Q) = SQRT(TRACE(BQ'*X*BQ+DQ'*DQ)),
C
C     where X satisfies the discrete-time Lyapunov equation
C
C        AQ'*X*AQ - X + CQ'*CQ = 0.
C
C     REFERENCES
C
C     [1] Varga A.
C         On computing 2-norms of transfer-function matrices.
C         Proc. 1992 ACC, Chicago, June 1992.
C
C     [2] Varga A.
C         A Schur method for computing coprime factorizations with
C         inner denominators and applications in model reduction.
C         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993.
C
C     NUMERICAL ASPECTS
C                                            3
C     The algorithm requires no more than 14N  floating point
C     operations.
C
C     CONTRIBUTOR
C
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine SL2NRM.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest.
C     Jan. 2003, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     Coprime factorization, Lyapunov equation, multivariable system,
C     state-space model, system norms.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOBN
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK, M,
     $                  N, NQ, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           DISCR
      INTEGER           KCR, KDR, KRW, KTAU, KU, MXNP, NR
      DOUBLE PRECISION  S2NORM, SCALE, WRKOPT
C     .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLANGE, DLAPY2
      EXTERNAL          DLANGE, DLAPY2, LSAME
C     .. External subroutines ..
      EXTERNAL          DLACPY, DTRMM, SB03OU, SB08DD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     .. Executable Statements ..
C
      DISCR = LSAME( DICO, 'D' )
      INFO  = 0
      IWARN = 0
C
C     Check the scalar input parameters.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( JOBN, 'H' ) .OR. LSAME( JOBN, 'L' ) ) )
     $      THEN
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
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.MAX( 1, M*( N + M ) +
     $                             MAX( N*( N + 5 ), M*( M + 2 ), 4*P ),
     $                           N*( MAX( N, P ) + 4 ) + MIN( N, P ) ) )
     $      THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 )THEN
C
C        Error return.
C
         AB13BD   = ZERO
         CALL XERBLA( 'AB13BD', -INFO )
         RETURN
      END IF
C
C     Compute the Frobenius norm of D.
C
      S2NORM = DLANGE( 'Frobenius', P, M, D, LDD, DWORK )
      IF( .NOT.DISCR .AND. S2NORM.NE.ZERO ) THEN
         AB13BD   = ZERO
         INFO = 5
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NQ = 0
         AB13BD   = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
      KCR = 1
      KDR = KCR + M*N
      KRW = KDR + M*M
C
C     Compute the right coprime factorization with inner denominator
C     of G.
C
C     Workspace needed:      M*(N+M);
C     Additional workspace:  need MAX( N*(N+5), M*(M+2), 4*M, 4*P );
C                            prefer larger.
C
      CALL SB08DD( DICO, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, NQ,
     $             NR, DWORK(KCR), M, DWORK(KDR), M, TOL, DWORK(KRW),
     $             LDWORK-KRW+1, IWARN, INFO )
      IF( INFO.NE.0 )
     $   RETURN
C
      WRKOPT = DWORK(KRW) + DBLE( KRW-1 )
C
C     Check stability.
C
      IF( LSAME( JOBN, 'H' ) .AND. NR.GT.0 ) THEN
         INFO = 6
         RETURN
      END IF
C
      IF( NQ.GT.0 ) THEN
         KU   = 1
         MXNP = MAX( NQ, P )
         KTAU = NQ*MXNP + 1
         KRW  = KTAU + MIN( NQ, P )
C
C        Find X, the solution of Lyapunov equation.
C
C        Workspace needed:      N*MAX(N,P) + MIN(N,P);
C        Additional workspace:  4*N;
C                               prefer larger.
C
         CALL DLACPY( 'Full', P, NQ, C, LDC, DWORK(KU), MXNP )
         CALL SB03OU( DISCR, .FALSE., NQ, P, A, LDA, DWORK(KU), MXNP,
     $                DWORK(KTAU), DWORK(KU), NQ, SCALE, DWORK(KRW),
     $                LDWORK-KRW+1, INFO )
         IF( INFO.NE.0 ) THEN
            IF( INFO.EQ.1 ) THEN
               INFO = 4
            ELSE IF( INFO.EQ.2 ) THEN
               INFO = 3
            END IF
            RETURN
         END IF
C
         WRKOPT = MAX( WRKOPT, DWORK(KRW) + DBLE( KRW-1 ) )
C
C        Add the contribution of BQ'*X*BQ.
C
C        Workspace needed:      N*(N+M).
C
         KTAU = NQ*NQ + 1
         CALL DLACPY( 'Full', NQ, M, B, LDB, DWORK(KTAU), NQ )
         CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', NQ, M,
     $               ONE, DWORK(KU), NQ, DWORK(KTAU), NQ )
         IF( NR.GT.0 )
     $      S2NORM = DLANGE( 'Frobenius', P, M, D, LDD, DWORK )
         S2NORM = DLAPY2( S2NORM, DLANGE( 'Frobenius', NQ, M,
     $                                    DWORK(KTAU), NQ, DWORK )
     $                               / SCALE )
      END IF
C
      AB13BD = S2NORM
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB13BD ***
      END
