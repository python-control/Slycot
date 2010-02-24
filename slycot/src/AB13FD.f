      SUBROUTINE AB13FD( N, A, LDA, BETA, OMEGA, TOL, DWORK, LDWORK,
     $                   CWORK, LCWORK, INFO )
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
C     To compute beta(A), the 2-norm distance from a real matrix A to
C     the nearest complex matrix with an eigenvalue on the imaginary
C     axis. If A is stable in the sense that all eigenvalues of A lie
C     in the open left half complex plane, then beta(A) is the complex
C     stability radius, i.e., the distance to the nearest unstable
C     complex matrix. The value of beta(A) is the minimum of the
C     smallest singular value of (A - jwI), taken over all real w.
C     The value of w corresponding to the minimum is also computed.
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
C     BETA    (output) DOUBLE PRECISION
C             The computed value of beta(A), which actually is an upper
C             bound.
C
C     OMEGA   (output) DOUBLE PRECISION
C             The value of w such that the smallest singular value of
C             (A - jwI) equals beta(A).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Specifies the accuracy with which beta(A) is to be
C             calculated. (See the Numerical Aspects section below.)
C             If the user sets TOL to be less than EPS, where EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH),
C             then the tolerance is taken to be EPS.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C             If DWORK(1) is not needed, the first 2*N*N entries of
C             DWORK may overlay CWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, 3*N*(N+2) ).
C             For optimum performance LDWORK should be larger.
C
C     CWORK   COMPLEX*16 array, dimension (LCWORK)
C             On exit, if INFO = 0, CWORK(1) returns the optimal value
C             of LCWORK.
C             If CWORK(1) is not needed, the first N*N entries of
C             CWORK may overlay DWORK.
C
C     LCWORK  INTEGER
C             The length of the array CWORK.
C             LCWORK >= MAX( 1, N*(N+3) ).
C             For optimum performance LCWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the routine fails to compute beta(A) within the
C                   specified tolerance. Nevertheless, the returned
C                   value is an upper bound on beta(A);
C             = 2:  either the QR or SVD algorithm (LAPACK Library
C                   routines DHSEQR, DGESVD or ZGESVD) fails to
C                   converge; this error is very rare.
C
C     METHOD
C
C     AB13FD combines the methods of [1] and [2] into a provably
C     reliable, quadratically convergent algorithm. It uses the simple
C     bisection strategy of [1] to find an interval which contains
C     beta(A), and then switches to the modified bisection strategy of
C     [2] which converges quadratically to a minimizer. Note that the
C     efficiency of the strategy degrades if there are several local
C     minima that are near or equal the global minimum.
C
C     REFERENCES
C
C     [1] Byers, R.
C         A bisection method for measuring the distance of a stable
C         matrix to the unstable matrices.
C         SIAM J. Sci. Stat. Comput., Vol. 9, No. 5, pp. 875-880, 1988.
C
C     [2] Boyd, S. and Balakrishnan, K.
C         A regularity result for the singular values of a transfer
C         matrix and a quadratically convergent algorithm for computing
C         its L-infinity norm.
C         Systems and Control Letters, Vol. 15, pp. 1-7, 1990.
C
C     NUMERICAL ASPECTS
C
C     In the presence of rounding errors, the computed function value
C     BETA  satisfies
C
C           beta(A) <= BETA + epsilon,
C
C           BETA/(1+TOL) - delta <= MAX(beta(A), SQRT(2*N*EPS)*norm(A)),
C
C     where norm(A) is the Frobenius norm of A,
C
C           epsilon = p(N) * EPS * norm(A),
C     and
C           delta   = p(N) * SQRT(EPS) * norm(A),
C
C     and p(N) is a low degree polynomial. It is recommended to choose
C     TOL greater than SQRT(EPS). Although rounding errors can cause
C     AB13FD to fail for smaller values of TOL, nevertheless, it usually
C     succeeds. Regardless of success or failure, the first inequality
C     holds.
C
C     CONTRIBUTORS
C
C     R. Byers, the routines QSEC and QSEC0 (January, 1995).
C
C     REVISIONS
C
C     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999.
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2002,
C     Jan. 2003.
C
C     KEYWORDS
C
C     complex stability radius, distances, eigenvalue, eigenvalue
C     perturbation, norms.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MAXIT
      PARAMETER         ( MAXIT = 50 )
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      COMPLEX*16        CONE
      PARAMETER         ( CONE = ( 1.0D0, 0.0D0 ) )
C     .. Scalar Arguments ..
      INTEGER           INFO, LCWORK, LDA, LDWORK, N
      DOUBLE PRECISION  BETA, OMEGA, TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*)
      COMPLEX*16        CWORK(*)
C     .. Local Scalars ..
      INTEGER           I, IA2, IAA, IGF, IHI, ILO, ITNUM, IWI, IWK,
     $                  IWR, JWORK, KOM, LBEST, MINWRK, N2
      DOUBLE PRECISION  EPS, LOW, OM, OM1, OM2, SFMN, SIGMA, SV, TAU,
     $                  TEMP, TOL1
      LOGICAL           SUFWRK
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1), DUMMY2(1,1)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE, MB03NY
      EXTERNAL          DLAMCH, DLANGE, MB03NY
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBAL, DGEMM, DHSEQR, DLACPY, DSYMM,
     $                  DSYMV, MA02ED, MB04ZD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO   = 0
      MINWRK = 3*N*( N + 2 )
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LDWORK.LT.MAX( 1, MINWRK ) ) THEN
         INFO = -8
      ELSE IF( LCWORK.LT.MAX( 1, N*( N + 3 ) ) ) THEN
         INFO = -10
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB13FD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      OMEGA = ZERO
      IF ( N.EQ.0 ) THEN
         BETA = ZERO
         DWORK(1) = ONE
         CWORK(1) = CONE
         RETURN
      END IF
C
C     Indices for splitting the work array.
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of workspace needed at that point in the code,
C     as well as the preferred amount for good performance.)
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
C     Computation of the tolerances. EPS is the machine precision.
C
      SFMN = DLAMCH( 'Safe minimum' )
      EPS  = DLAMCH( 'Epsilon' )
      TOL1 = SQRT( EPS * DBLE( 2*N ) ) *
     $       DLANGE( 'Frobenius', N, N, A, LDA, DWORK )
      TAU  = ONE + MAX( TOL, EPS )
C
C     Initialization, upper bound at known critical point.
C     Workspace: need N*(N+1)+5*N; prefer larger.
C
      KOM = 2
      LOW = ZERO
      CALL DLACPY( 'All', N, N, A, LDA, DWORK(IGF), N )
      BETA = MB03NY( N, OMEGA, DWORK(IGF), N, DWORK(IGF+N2),
     $               DWORK(IA2), LDWORK-IA2, CWORK, LCWORK, INFO )
      IF ( INFO.NE.0 )
     $   RETURN
      LBEST = MAX( MINWRK, INT( DWORK(IA2) ) - IA2 + 1, 4*N2 + N )
C
      ITNUM = 1
C     WHILE ( ITNUM <= MAXIT and BETA > TAU*MAX( TOL1, LOW ) ) DO
   10 IF ( ( ITNUM.LE.MAXIT ) .AND.
     $     ( BETA.GT.TAU*MAX( TOL1, LOW ) ) ) THEN
         IF ( KOM.EQ.2 ) THEN
            SIGMA = BETA/TAU
         ELSE
            SIGMA = SQRT( BETA ) * SQRT( MAX( TOL1, LOW ) )
         END IF
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
            INFO = 2
            RETURN
         END IF
C
C        Count negative real axis squared eigenvalues. If there are two,
C        then the valley is isolated, and next approximate minimizer is
C        mean of the square roots.
C
         KOM = 0
         DO 40 I = 0, N - 1
            TEMP = ABS( DWORK(IWI+I) )
            IF ( TOL1.GT.SFMN ) TEMP = TEMP / TOL1
            IF ( ( DWORK(IWR+I).LT.ZERO ) .AND. ( TEMP.LE.TOL1 ) ) THEN
               KOM = KOM + 1
               OM = SQRT( -DWORK(IWR+I) )
               IF ( KOM.EQ.1 ) OM1 = OM
               IF ( KOM.EQ.2 ) OM2 = OM
            END IF
   40    CONTINUE
C
         IF ( KOM.EQ.0 ) THEN
            LOW = SIGMA
         ELSE
C
C           In exact arithmetic KOM = 1 is impossible, but if tau is
C           close enough to one, MB04ZD may miss the initial near zero
C           eigenvalue.
C           Workspace, real:    need   3*N*(N+2);  prefer larger;
C                      complex: need     N*(N+3);  prefer larger.
C
            IF ( KOM.EQ.2 ) THEN
               OM = OM1 + ( OM2 - OM1 ) / TWO
            ELSE IF ( KOM.EQ.1 .AND. ITNUM.EQ.1 ) THEN
               OM  = OM1 / TWO
               KOM = 2
            END IF
C
            CALL DLACPY( 'All', N, N, A, LDA, DWORK(IGF), N )
            SV = MB03NY( N, OM, DWORK(IGF), N, DWORK(IGF+N2),
     $                   DWORK(IA2), LDWORK-IA2, CWORK, LCWORK, INFO )
            IF ( INFO.NE.0 )
     $         RETURN
            IF ( BETA.GT.SV ) THEN
               BETA  = SV
               OMEGA = OM
            ELSE
               INFO = 1
               RETURN
            END IF
         END IF
         ITNUM = ITNUM + 1
         GO TO 10
C        END WHILE 10
      END IF
C
      IF ( BETA .GT. TAU*MAX( TOL1, LOW ) ) THEN
C
C        Failed to meet bounds within MAXIT iterations.
C
         INFO = 1
         RETURN
      END IF
C
C     Set optimal real workspace dimension (complex workspace is already
C     set by MB03NY).
C
      DWORK(1) = LBEST
C
      RETURN
C *** Last line of AB13FD ***
      END
