      SUBROUTINE AB13ED( N, A, LDA, LOW, HIGH, TOL, DWORK, LDWORK,
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
C     To estimate beta(A), the 2-norm distance from a real matrix A to
C     the nearest complex matrix with an eigenvalue on the imaginary
C     axis. The estimate is given as
C
C            LOW <= beta(A) <= HIGH,
C
C     where either
C
C            (1 + TOL) * LOW >= HIGH,
C
C     or
C
C            LOW = 0   and   HIGH = delta,
C
C     and delta is a small number approximately equal to the square root
C     of machine precision times the Frobenius norm (Euclidean norm)
C     of A. If A is stable in the sense that all eigenvalues of A lie
C     in the open left half complex plane, then beta(A) is the distance
C     to the nearest unstable complex matrix, i.e., the complex
C     stability radius.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     LOW     (output) DOUBLE PRECISION
C             A lower bound for beta(A).
C
C     HIGH    (output) DOUBLE PRECISION
C             An upper bound for beta(A).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Specifies the accuracy with which LOW and HIGH approximate
C             beta(A). If the user sets TOL to be less than SQRT(EPS),
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH), then the tolerance is taken to be
C             SQRT(EPS).
C             The recommended value is TOL = 9, which gives an estimate
C             of beta(A) correct to within an order of magnitude.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, 3*N*(N+1) ).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the QR algorithm (LAPACK Library routine DHSEQR)
C                   fails to converge; this error is very rare.
C
C     METHOD
C
C     Let beta(A) be the 2-norm distance from a real matrix A to the
C     nearest complex matrix with an eigenvalue on the imaginary axis.
C     It is known that beta(A) = minimum of the smallest singular
C     value of (A - jwI), where I is the identity matrix and j**2 = -1,
C     and the minimum is taken over all real w.
C     The algorithm computes a lower bound LOW and an upper bound HIGH
C     for beta(A) by a bisection method in the following way. Given a
C     non-negative real number sigma, the Hamiltonian matrix H(sigma)
C     is constructed:
C
C                       |   A      -sigma*I |     | A   G  |
C           H(sigma) =  |                   | :=  |        | .
C                       | sigma*I    -A'    |     | F  -A' |
C
C     It can be shown [1] that H(sigma) has an eigenvalue whose real
C     part is zero if and only if sigma >= beta. Any lower and upper
C     bounds on beta(A) can be improved by choosing a number between
C     them and checking to see if H(sigma) has an eigenvalue with zero
C     real part.  This decision is made by computing the eigenvalues of
C     H(sigma) using the square reduced algorithm of Van Loan [2].
C
C     REFERENCES
C
C     [1] Byers, R.
C         A bisection method for measuring the distance of a stable
C         matrix to the unstable matrices.
C         SIAM J. Sci. Stat. Comput., Vol. 9, No. 5, pp. 875-880, 1988.
C
C     [2] Van Loan, C.F.
C         A symplectic method for approximating all the eigenvalues of a
C         Hamiltonian matrix.
C         Linear Algebra and its Applications, Vol 61, 233-251, 1984.
C
C     NUMERICAL ASPECTS
C
C     Due to rounding errors the computed values of LOW and HIGH can be
C     proven to satisfy
C
C            LOW - p(n) * sqrt(e) * norm(A) <= beta(A)
C     and
C            beta(A) <= HIGH + p(n) * sqrt(e) * norm(A),
C
C     where p(n) is a modest polynomial of degree 3, e is the machine
C     precision and norm(A) is the Frobenius norm of A, see [1].
C     The recommended value for TOL is 9 which gives an estimate of
C     beta(A) correct to within an order of magnitude.
C     AB13ED requires approximately 38*N**3 flops for TOL = 9.
C
C     CONTRIBUTOR
C
C     R. Byers, the routines BISEC and BISEC0 (January, 1995).
C
C     REVISIONS
C
C     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999.
C     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2003.
C
C     KEYWORDS
C
C     Distances, eigenvalue, eigenvalue perturbation, norms, stability
C     radius.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HIGH, LOW, TOL
      INTEGER           INFO, LDA, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*)
C     .. Local Scalars ..
      INTEGER           I, IA2, IAA, IGF, IHI, ILO, IWI, IWK, IWR,
     $                  JWORK, MINWRK, N2
      DOUBLE PRECISION  ANRM, SEPS, SFMN, SIGMA, TAU, TEMP, TOL1, TOL2
      LOGICAL           RNEG, SUFWRK
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1), DUMMY2(1,1)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBAL, DGEMM, DHSEQR, DLACPY, DSYMM,
     $                  DSYMV, MA02ED, MB04ZD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO   = 0
      MINWRK = 3*N*( N + 1 )
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LDWORK.LT.MAX( 1, MINWRK ) ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB13ED', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      LOW = ZERO
      IF ( N.EQ.0 ) THEN
         HIGH = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Indices for splitting the work array.
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      N2  = N*N
      IGF = 1
      IA2 = IGF + N2 + N
      IAA = IA2 + N2
      IWK = IAA + N2
      IWR = IAA
      IWI = IWR + N
C
      SUFWRK = LDWORK-IWK.GE.N2
C
C     Computation of the tolerances and the treshold for termination of
C     the bisection method. SEPS is the square root of the machine
C     precision.
C
      SFMN = DLAMCH( 'Safe minimum' )
      SEPS = SQRT( DLAMCH( 'Epsilon' ) )
      TAU  = ONE + MAX( TOL, SEPS )
      ANRM = DLANGE( 'Frobenius', N, N, A, LDA, DWORK )
      TOL1 = SEPS * ANRM
      TOL2 = TOL1 * DBLE( 2*N )
C
C     Initialization of the bisection method.
C
      HIGH = ANRM
C
C     WHILE ( HIGH > TAU*MAX( TOL1, LOW ) ) DO
   10 IF ( HIGH.GT.( TAU*MAX( TOL1, LOW ) ) ) THEN
         SIGMA = SQRT( HIGH ) * SQRT( MAX( TOL1, LOW ) )
C
C        Set up H(sigma).
C        Workspace: N*(N+1)+2*N*N.
C
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(IAA), N )
         DWORK(IGF)   =  SIGMA
         DWORK(IGF+N) = -SIGMA
         DUMMY(1) = ZERO
         CALL DCOPY( N-1, DUMMY, 0, DWORK(IGF+1), 1 )
C
         DO 20 I = IGF, IA2 - N - 2, N + 1
            CALL DCOPY( N+1, DWORK(I), 1, DWORK(I+N+1), 1 )
   20    CONTINUE
C
C        Computation of the eigenvalues by the square reduced algorithm.
C        Workspace: N*(N+1)+2*N*N+2*N.
C
         CALL MB04ZD( 'No vectors', N, DWORK(IAA), N, DWORK(IGF), N,
     $                DUMMY2, 1, DWORK(IWK), INFO )
C
C        Form the matrix A*A + F*G.
C        Workspace: need   N*(N+1)+2*N*N+N;
C                   prefer N*(N+1)+3*N*N.
C
         JWORK = IA2
         IF ( SUFWRK )
     $      JWORK = IWK
C
         CALL DLACPY( 'Lower', N, N, DWORK(IGF), N, DWORK(JWORK), N )
         CALL MA02ED( 'Lower', N, DWORK(JWORK), N )
C
         IF ( SUFWRK ) THEN
C
C           Use BLAS 3 calculation.
C
            CALL DSYMM( 'Left', 'Upper', N, N, ONE, DWORK(IGF+N), N,
     $                  DWORK(JWORK), N, ZERO, DWORK(IA2), N )
         ELSE
C
C           Use BLAS 2 calculation.
C
            DO 30 I = 1, N
               CALL DSYMV( 'Upper', N, ONE, DWORK(IGF+N), N,
     $                     DWORK(IA2+N*(I-1)), 1, ZERO, DWORK(IWK), 1 )
               CALL DCOPY( N, DWORK(IWK), 1, DWORK(IA2+N*(I-1)), 1 )
   30       CONTINUE
C
         END IF
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, N, ONE,
     $               DWORK(IAA), N, DWORK(IAA), N, ONE, DWORK(IA2), N )
C
C        Find the eigenvalues of A*A + F*G.
C        Workspace: N*(N+1)+N*N+3*N.
C
         JWORK = IWI + N
         CALL DGEBAL( 'Scale', N, DWORK(IA2), N, ILO, IHI, DWORK(JWORK),
     $                I )
         CALL DHSEQR( 'Eigenvalues', 'NoSchurVectors', N, ILO, IHI,
     $                DWORK(IA2), N, DWORK(IWR), DWORK(IWI), DUMMY2, 1,
     $                DWORK(JWORK), N, INFO )
C
         IF ( INFO.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
C
C        (DWORK(IWR+i), DWORK(IWI+i)), i = 0,...,N-1, contain the
C        squares of the eigenvalues of H(sigma).
C
         I = 0
         RNEG = .FALSE.
C        WHILE ( ( DWORK(IWR+i),DWORK(IWI+i) ) not real positive
C                .AND. I < N ) DO
   40    IF ( .NOT.RNEG .AND. I.LT.N ) THEN
            TEMP = ABS( DWORK(IWI+I) )
            IF ( TOL1.GT.SFMN ) TEMP = TEMP / TOL1
            RNEG = ( ( DWORK(IWR+I).LT.ZERO ) .AND. ( TEMP.LE.TOL2 ) )
            I = I + 1
            GO TO 40
C           END WHILE 40
         END IF

         IF ( RNEG ) THEN
            HIGH = SIGMA
         ELSE
            LOW = SIGMA
         END IF
         GO TO 10
C        END WHILE 10
      END IF
C
C     Set optimal workspace dimension.
C
      DWORK(1) = DBLE( MAX( 4*N2 + N, MINWRK ) )
C
C *** Last line of AB13ED ***
      END
