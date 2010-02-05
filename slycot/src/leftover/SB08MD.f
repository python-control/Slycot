      SUBROUTINE SB08MD( ACONA, DA, A, RES, E, DWORK, LDWORK, INFO )
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
C     To compute a real polynomial E(s) such that
C
C        (a)  E(-s) * E(s) = A(-s) * A(s) and
C        (b)  E(s) is stable - that is, all the zeros of E(s) have
C             non-positive real parts,
C
C     which corresponds to computing the spectral factorization of the
C     real polynomial A(s) arising from continuous optimality problems.
C
C     The input polynomial may be supplied either in the form
C
C        A(s) = a(0) + a(1) * s + ... + a(DA) * s**DA
C
C     or as
C
C        B(s) = A(-s) * A(s)
C             = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA)        (1)
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ACONA   CHARACTER*1
C             Indicates whether the coefficients of A(s) or B(s) =
C             A(-s) * A(s) are to be supplied as follows:
C             = 'A':  The coefficients of A(s) are to be supplied;
C             = 'B':  The coefficients of B(s) are to be supplied.
C
C     Input/Output Parameters
C
C     DA      (input) INTEGER
C             The degree of the polynomials A(s) and E(s).  DA >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (DA+1)
C             On entry, this array must contain either the coefficients
C             of the polynomial A(s) in increasing powers of s if
C             ACONA = 'A', or the coefficients of the polynomial B(s) in
C             increasing powers of s**2 (see equation (1)) if ACONA =
C             'B'.
C             On exit, this array contains the coefficients of the
C             polynomial B(s) in increasing powers of s**2.
C
C     RES     (output) DOUBLE PRECISION
C             An estimate of the accuracy with which the coefficients of
C             the polynomial E(s) have been computed (see also METHOD
C             and NUMERICAL ASPECTS).
C
C     E       (output) DOUBLE PRECISION array, dimension (DA+1)
C             The coefficients of the spectral factor E(s) in increasing
C             powers of s.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 5*DA+5.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if on entry, A(I) = 0.0, for I = 1,2,...,DA+1.
C             = 2:  if on entry, ACONA = 'B' but the supplied
C                   coefficients of the polynomial B(s) are not the
C                   coefficients of A(-s) * A(s) for some real A(s);
C                   in this case, RES and E are unassigned;
C             = 3:  if the iterative process (see METHOD) has failed to
C                   converge in 30 iterations;
C             = 4:  if the last computed iterate (see METHOD) is
C                   unstable. If ACONA = 'B', then the supplied
C                   coefficients of the polynomial B(s) may not be the
C                   coefficients of A(-s) * A(s) for some real A(s).
C
C     METHOD
C         _                                               _
C     Let A(s) be the conjugate polynomial of A(s), i.e., A(s) = A(-s).
C
C     The method used by the routine is based on applying the
C     Newton-Raphson iteration to the function
C               _       _
C        F(e) = A * A - e * e,
C
C     which leads to the iteration formulae (see [1]):
C
C        _(i)   (i)  _(i)   (i)     _      )
C        q   * x   + x   * q    = 2 A * A  )
C                                          )   for i = 0, 1, 2,...
C         (i+1)    (i)   (i)               )
C        q     = (q   + x   )/2            )
C
C                    (0)         DA
C     Starting from q   = (1 + s)   (which has no zeros in the closed
C                                                  (1)   (2)   (3)
C     right half-plane), the sequence of iterates q   , q   , q   ,...
C     converges to a solution of F(e) = 0 which has no zeros in the
C     open right half-plane.
C
C     The iterates satisfy the following conditions:
C
C              (i)
C        (a)  q   is a stable polynomial (no zeros in the closed right
C             half-plane) and
C
C              (i)        (i-1)
C        (b)  q   (1) <= q     (1).
C
C                                       (i-1)                       (i)
C     The iterative process stops with q     , (where i <= 30)  if q
C     violates either (a) or (b), or if the condition
C                       _(i) (i)  _
C        (c)  RES  = ||(q   q   - A A)|| < tol,
C
C     is satisfied, where || . || denotes the largest coefficient of
C                     _(i) (i)  _
C     the polynomial (q   q   - A A) and tol is an estimate of the
C                                                    _(i)  (i)
C     rounding error in the computed coefficients of q    q   . If there
C     is no convergence after 30 iterations then the routine returns
C     with the Error Indicator (INFO) set to 3, and the value of RES may
C     indicate whether or not the last computed iterate is close to the
C     solution.
C
C     If ACONA = 'B', then it is possible that the equation e(-s) *
C     e(s) = B(s) has no real solution, which will be the case if A(1)
C     < 0 or if ( -1)**DA * A(DA+1) < 0.
C
C     REFERENCES
C
C     [1] Vostry, Z.
C         New Algorithm for Polynomial Spectral Factorization with
C         Quadratic Convergence II.
C         Kybernetika, 12, pp. 248-259, 1976.
C
C     NUMERICAL ASPECTS
C
C     The conditioning of the problem depends upon the distance of the
C     zeros of A(s) from the imaginary axis and on their multiplicity.
C     For a well-conditioned problem the accuracy of the computed
C     coefficients of E(s) is of the order of RES. However, for problems
C     with zeros near the imaginary axis or with multiple zeros, the
C     value of RES may be an overestimate of the true accuracy.
C
C     FURTHER COMMENTS
C
C     In order for the problem e(-s) * e(s) = B(s) to have a real
C     solution e(s), it is necessary and sufficient that B(j*omega)
C     >= 0 for any purely imaginary argument j*omega (see [1]).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB08AD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Factorization, Laplace transform, optimal control, optimal
C     filtering, polynomial operations, spectral factorization, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         ACONA
      INTEGER           DA, INFO, LDWORK
      DOUBLE PRECISION  RES
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), DWORK(*), E(*)
C     .. Local Scalars ..
      LOGICAL           CONV, LACONA, STABLE
      INTEGER           BINC, DA1, I, I0, J, K, LAMBDA, LAY, LAYEND,
     $                  LDIF, LPHEND, LPHI, LQ, M, NC
      DOUBLE PRECISION  A0, EPS, MU, MUJ, SI, SIGNI, SIGNI0, SIGNJ,
     $                  SIMIN1, SQRTA0, SQRTMJ, SQRTMU, TOLPHI, W, XDA
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           IDAMAX
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, SB08MY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MOD, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      LACONA = LSAME( ACONA, 'A' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LACONA .AND. .NOT.LSAME( ACONA, 'B' ) ) THEN
         INFO = -1
      ELSE IF( DA.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDWORK.LT.5*DA + 5 ) THEN
         INFO = -7
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB08MD', -INFO )
         RETURN
      END IF
C
      IF ( .NOT.LACONA ) THEN
         CALL DCOPY( DA+1, A, 1, E, 1 )
      ELSE
         W = ZERO
         CALL SB08MY( DA, A, E, W )
      END IF
C
C     Reduce E such that the first and the last element are non-zero.
C
      DA1 = DA + 1
C
C     WHILE ( DA1 >= 1 and E(DA1) = 0 ) DO
   20 IF ( DA1.GE.1 ) THEN
         IF ( E(DA1).EQ.ZERO ) THEN
            DA1 = DA1 - 1
            GO TO 20
         END IF
      END IF
C     END WHILE 20
C
      DA1 = DA1 - 1
      IF ( DA1.LT.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
      I0 = 1
C
C     WHILE ( E(I0) = 0 ) DO
   40 IF ( E(I0).EQ.ZERO ) THEN
         I0 = I0 + 1
         GO TO 40
      END IF
C     END WHILE 40
C
      I0 = I0 - 1
      IF ( I0.NE.0 ) THEN
         IF ( MOD( I0, 2 ).EQ.0 ) THEN
            SIGNI0 = ONE
         ELSE
            SIGNI0 = -ONE
         END IF
C
         DO 60 I = 1, DA1 - I0 + 1
            E(I) = SIGNI0*E(I+I0)
   60    CONTINUE
C
         DA1 = DA1 - I0
      END IF
      IF ( MOD( DA1, 2 ).EQ.0 ) THEN
         SIGNI = ONE
      ELSE
         SIGNI = -ONE
      END IF
      NC = DA1 + 1
      IF ( ( E(1).LT.ZERO ) .OR. ( ( E(NC)*SIGNI ).LT.ZERO ) ) THEN
         INFO = 2
         RETURN
      END IF
C
C     Initialization.
C
      EPS = DLAMCH( 'Epsilon' )
      SI = ONE/DLAMCH( 'Safe minimum' )
      LQ = 1
      LAY = LQ + NC
      LAMBDA = LAY + NC
      LPHI = LAMBDA + NC
      LDIF = LPHI + NC
C
      A0 = E(1)
      BINC = 1
C
C     Computation of the starting polynomial and scaling of the input
C     polynomial.
C
      MU = ( A0/ABS( E(NC) ) )**( ONE/DBLE( DA1 ) )
      MUJ = ONE
C
      DO 80 J = 1, NC
         W = E(J)*MUJ/A0
         A(J) = W
         E(J) = BINC
         DWORK(LQ+J-1) = BINC
         MUJ = MUJ*MU
         BINC = BINC*( NC - J )/J
   80 CONTINUE
C
      CONV = .FALSE.
      STABLE = .TRUE.
C
C     The contents of the arrays is, cf [1],
C
C     E : the last computed stable polynomial q   ;
C                                              i-1
C     DWORK(LAY+1,...,LAY+DA1-1)  : a'(1), ..., a'(DA1-1), these values
C                                   are changed during the computation
C                                   into y;
C          (LAMBDA+1,...,LAMBDA+DA1-2) : lambda(1), ..., lambda(DA1-2),
C                                        the factors of the Routh
C                                        stability test, (lambda(i) is
C                                        P(i) in [1]);
C          (LPHI+1,...,LPHI+DA1-1) : phi(1), ..., phi(DA1-1), the values
C                                    phi(i,j), see [1], scheme (11);
C          (LDIF,...,LDIF+DA1) : the coeffs of q (-s) * q (s) - b(s).
C                                               i        i
C     DWORK(LQ,...,LQ+DA1) : the last computed polynomial q .
C                                                          i
      I = 0
C
C     WHILE ( I < 30 and CONV = FALSE and STABLE = TRUE ) DO
  100 IF ( I.LT.30 .AND. .NOT.CONV .AND. STABLE ) THEN
         I = I + 1
         CALL DCOPY( NC, A, 1, DWORK(LAY), 1 )
         CALL DCOPY( NC, DWORK(LQ), 1, DWORK(LPHI), 1 )
         M = DA1/2
         LAYEND = LAY + DA1
         LPHEND = LPHI + DA1
         XDA = A(NC)/DWORK(LQ+DA1)
C
         DO 120 K = 1, M
            DWORK(LAY+K) = DWORK(LAY+K) - DWORK(LPHI+2*K)
            DWORK(LAYEND-K) = DWORK(LAYEND-K) - DWORK(LPHEND-2*K)*XDA
  120    CONTINUE
C
C        Computation of lambda(k) and y(k).
C
         K = 1
C
C        WHILE ( K <= DA1 - 2 and STABLE = TRUE ) DO
  140    IF ( ( K.LE.( DA1 - 2 ) ) .AND. STABLE ) THEN
            IF ( DWORK(LPHI+K).LE.ZERO ) STABLE = .FALSE.
            IF ( STABLE ) THEN
               W = DWORK(LPHI+K-1)/DWORK(LPHI+K)
               DWORK(LAMBDA+K) = W
               CALL DAXPY( ( DA1 - K )/2, -W, DWORK(LPHI+K+2), 2,
     $                                        DWORK(LPHI+K+1), 2 )
               W = DWORK(LAY+K)/DWORK(LPHI+K)
               DWORK(LAY+K) = W
               CALL DAXPY( ( DA1 - K )/2, -W, DWORK(LPHI+K+2), 2,
     $                                        DWORK(LAY+K+1), 1 )
               K = K + 1
            END IF
            GO TO 140
         END IF
C        END WHILE 140
C
         IF ( DWORK(LPHI+DA1-1).LE.ZERO ) THEN
            STABLE = .FALSE.
         ELSE
            DWORK(LAY+DA1-1) = DWORK(LAY+DA1-1)/DWORK(LPHI+DA1-1)
         END IF
C
C        STABLE = The polynomial q    is stable.
C                                 i-1
         IF ( STABLE ) THEN
C
C           Computation of x  and q .
C                           i      i
C
            DO 160 K = DA1 - 2, 1, -1
               W = DWORK(LAMBDA+K)
               CALL DAXPY( ( DA1 - K )/2, -W, DWORK(LAY+K+1), 2,
     $                                        DWORK(LAY+K), 2 )
  160       CONTINUE
C
            DWORK(LAY+DA1) = XDA
C
            CALL DCOPY( NC, DWORK(LQ), 1, E, 1 )
            SIMIN1 = SI
            SI = DWORK(LQ)
            SIGNJ = -ONE
C
            DO 180 J = 1, DA1
               W = HALF*( DWORK(LQ+J) + SIGNJ*DWORK(LAY+J) )
               DWORK(LQ+J) = W
               SI = SI + W
               SIGNJ = -SIGNJ
  180       CONTINUE
C
            TOLPHI = EPS
            CALL SB08MY( DA1, E, DWORK(LDIF), TOLPHI )
            CALL DAXPY( NC, -ONE, A, 1, DWORK(LDIF), 1 )
            RES = ABS( DWORK( IDAMAX( NC, DWORK(LDIF), 1 ) + LDIF-1 ) )
C
C           Convergency test.
C
            IF ( ( SI.GT.SIMIN1 ) .OR. ( RES.LT.TOLPHI ) ) THEN
               CONV = .TRUE.
            END IF
            GO TO 100
         END IF
      END IF
C     END WHILE 100
C
C     Backscaling.
C
      MU = ONE/MU
      SQRTA0 = SQRT( A0 )
      SQRTMU = SQRT( MU )
      MUJ = ONE
      SQRTMJ = ONE
C
      DO 200 J = 1, NC
         E(J) = E(J)*SQRTA0*SQRTMJ
         A(J) = A(J)*A0*MUJ
         MUJ = MUJ*MU
         SQRTMJ = SQRTMJ*SQRTMU
  200 CONTINUE
C
      IF ( I0.NE.0 ) THEN
C
         DO 220 J = NC, 1, -1
            E(I0+J) = E(J)
            A(I0+J) = SIGNI0*A(J)
  220    CONTINUE
C
         DO 240 J = 1, I0
            E(J) = ZERO
            A(J) = ZERO
  240    CONTINUE
C
      END IF
C
      IF ( .NOT.CONV ) THEN
         IF ( STABLE ) THEN
            INFO = 3
         ELSE
            INFO = 4
         END IF
      END IF
C
      RETURN
C *** Last line of SB08MD ***
      END
