      SUBROUTINE MB03MD( N, L, THETA, Q, E, Q2, E2, PIVMIN, TOL, RELTOL,
     $                   IWARN, INFO )
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
C     To compute an upper bound THETA using a bisection method such that
C     the bidiagonal matrix
C
C              |q(1) e(1)  0    ...   0   |
C              | 0   q(2) e(2)        .   |
C          J = | .                    .   |
C              | .                  e(N-1)|
C              | 0   ...        ...  q(N) |
C
C     has precisely L singular values less than or equal to THETA plus
C     a given tolerance TOL.
C
C     This routine is mainly intended to be called only by other SLICOT
C     routines.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the bidiagonal matrix J.  N >= 0.
C
C     L       (input/output) INTEGER
C             On entry, L must contain the number of singular values
C             of J which must be less than or equal to the upper bound
C             computed by the routine.  0 <= L <= N.
C             On exit, L may be increased if the L-th smallest singular
C             value of J has multiplicity greater than 1. In this case,
C             L is increased by the number of singular values of J which
C             are larger than its L-th smallest one and approach the
C             L-th smallest singular value of J within a distance less
C             than TOL.
C             If L has been increased, then the routine returns with
C             IWARN set to 1.
C
C     THETA   (input/output) DOUBLE PRECISION
C             On entry, THETA must contain an initial estimate for the
C             upper bound to be computed. If THETA < 0.0 on entry, then
C             one of the following default values is used.
C             If L = 0, THETA is set to 0.0 irrespective of the input
C             value of THETA; if L = 1, then THETA is taken as
C             MIN(ABS(Q(i))), for i = 1,2,...,N; otherwise, THETA is
C             taken as ABS(Q(N-L+1)).
C             On exit, THETA contains the computed upper bound such that
C             the bidiagonal matrix J has precisely L singular values
C             less than or equal to THETA + TOL.
C
C     Q       (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the diagonal elements q(1),
C             q(2),...,q(N) of the bidiagonal matrix J. That is,
C             Q(i) = J(i,i) for i = 1,2,...,N.
C
C     E       (input) DOUBLE PRECISION array, dimension (N-1)
C             This array must contain the superdiagonal elements
C             e(1),e(2),...,e(N-1) of the bidiagonal matrix J. That is,
C             E(k) = J(k,k+1) for k = 1,2,...,N-1.
C
C     Q2      (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the squares of the diagonal
C             elements q(1),q(2),...,q(N) of the bidiagonal matrix J.
C             That is, Q2(i) = J(i,i)**2 for i = 1,2,...,N.
C
C     E2      (input) DOUBLE PRECISION array, dimension (N-1)
C             This array must contain the squares of the superdiagonal
C             elements e(1),e(2),...,e(N-1) of the bidiagonal matrix J.
C             That is, E2(k) = J(k,k+1)**2 for k = 1,2,...,N-1.
C
C     PIVMIN  (input) DOUBLE PRECISION
C             The minimum absolute value of a "pivot" in the Sturm
C             sequence loop.
C             PIVMIN >= max( max( |q(i)|, |e(k)| )**2*sf_min, sf_min ),
C             where i = 1,2,...,N, k = 1,2,...,N-1, and sf_min is at
C             least the smallest number that can divide one without
C             overflow (see LAPACK Library routine DLAMCH).
C             Note that this condition is not checked by the routine.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             This parameter defines the multiplicity of singular values
C             by considering all singular values within an interval of
C             length TOL as coinciding. TOL is used in checking how many
C             singular values are less than or equal to THETA. Also in
C             computing an appropriate upper bound THETA by a bisection
C             method, TOL is used as a stopping criterion defining the
C             minimum (absolute) subinterval width.  TOL >= 0.
C
C     RELTOL  DOUBLE PRECISION
C             This parameter specifies the minimum relative width of an
C             interval. When an interval is narrower than TOL, or than
C             RELTOL times the larger (in magnitude) endpoint, then it
C             is considered to be sufficiently small and bisection has
C             converged.
C             RELTOL >= BASE * EPS, where BASE is machine radix and EPS
C             is machine precision (see LAPACK Library routine DLAMCH).
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warnings;
C             = 1:  if the value of L has been increased as the L-th
C                   smallest singular value of J coincides with the
C                   (L+1)-th smallest one.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Let s(i), i = 1,2,...,N, be the N non-negative singular values of
C     the bidiagonal matrix J arranged so that s(1) >= ... >= s(N) >= 0.
C     The routine then computes an upper bound T such that s(N-L) > T >=
C     s(N-L+1) as follows (see [2]).
C     First, if the initial estimate of THETA is not specified by the
C     user then the routine initialises THETA to be an estimate which
C     is close to the requested value of THETA if s(N-L) >> s(N-L+1).
C     Second, a bisection method (see [1, 8.5]) is used which generates
C     a sequence of shrinking intervals [Y,Z] such that either THETA in
C     [Y,Z] was found (so that J has L singular values less than or
C     equal to THETA), or
C
C        (number of s(i) <= Y) < L < (number of s(i) <= Z).
C
C     This bisection method is applied to an associated 2N-by-2N
C     symmetric tridiagonal matrix T" whose eigenvalues (see [1]) are
C     given by s(1),s(2),...,s(N),-s(1),-s(2),...,-s(N). One of the
C     starting values for the bisection method is the initial value of
C     THETA. If this value is an upper bound, then the initial lower
C     bound is set to zero, else the initial upper bound is computed
C     from the Gershgorin Circle Theorem [1, Theorem 7.2-1], applied to
C     T". The computation of the "number of s(i) <= Y (or Z)" is
C     achieved by calling SLICOT Library routine MB03ND, which applies
C     Sylvester's Law of Inertia or equivalently Sturm sequences
C     [1, 8.5] to the associated matrix T". If
C
C        Z - Y <= MAX( TOL, PIVMIN, RELTOL*MAX( ABS( Y ), ABS( Z ) ) )
C
C     at some stage of the bisection method, then at least two singular
C     values of J lie in the interval [Y,Z] within a distance less than
C     TOL from each other. In this case, s(N-L) and s(N-L+1) are assumed
C     to coincide, the upper bound T is set to the value of Z, the value
C     of L is increased and IWARN is set to 1.
C
C     REFERENCES
C
C     [1] Golub, G.H. and Van Loan, C.F.
C         Matrix Computations.
C         The Johns Hopkins University Press, Baltimore, Maryland, 1983.
C
C     [2] Van Huffel, S. and Vandewalle, J.
C         The Partial Total Least Squares Algorithm.
C         J. Comput. and Appl. Math., 21, pp. 333-341, 1988.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB03AD by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     June 16, 1997, Oct. 26, 2003.
C
C     KEYWORDS
C
C     Bidiagonal matrix, singular values.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, TWO
      PARAMETER         ( ZERO = 0.0D0, TWO = 2.0D0 )
      DOUBLE PRECISION  FUDGE
      PARAMETER         ( FUDGE = TWO )
C     .. Scalar Arguments ..
      INTEGER           INFO, IWARN, L, N
      DOUBLE PRECISION  PIVMIN, RELTOL, THETA, TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), E2(*), Q(*), Q2(*)
C     .. Local Scalars ..
      INTEGER           I, NUM, NUMZ
      DOUBLE PRECISION  H, TH, Y, Z
C     .. External Functions ..
      INTEGER           MB03ND
      DOUBLE PRECISION  DLAMCH, MB03MY
      EXTERNAL          DLAMCH, MB03MY, MB03ND
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX
C     .. Executable Statements ..
C
C     Test some input scalar arguments.
C
      IWARN = 0
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( L.LT.0 .OR. L.GT.N ) THEN
         INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $    RETURN
C
C     Step 1: initialisation of THETA.
C             -----------------------
      IF ( L.EQ.0 ) THETA = ZERO
      IF ( THETA.LT.ZERO ) THEN
         IF ( L.EQ.1 ) THEN
C
C           An upper bound which is close if S(N-1) >> S(N):
C
            THETA = MB03MY( N, Q, 1 )
            IF ( N.EQ.1 )
     $         RETURN
         ELSE
C
C           An experimentally established estimate which is good if
C           S(N-L) >> S(N-L+1):
C
            THETA = ABS( Q(N-L+1) )
         END IF
      END IF
C
C     Step 2: Check quality of initial estimate THETA.
C             ---------------------------------------
      NUM = MB03ND( N, THETA, Q2, E2, PIVMIN, INFO )
      IF ( NUM.EQ.L )
     $   RETURN
C
C     Step 3: initialisation starting values for bisection method.
C             ---------------------------------------------------
C     Let S(i), i=1,...,N, be the singular values of J in decreasing
C     order. Then, the computed Y and Z will be such that
C     (number of S(i) <= Y) < L < (number of S(i) <= Z).
C
      IF ( NUM.LT.L ) THEN
         TH = ABS( Q(1) )
         Z = ZERO
         Y = THETA
         NUMZ = N
C
         DO 20 I = 1, N - 1
            H = ABS( Q(I+1) )
            Z  = MAX( MAX( TH, H ) + ABS( E(I) ), Z )
            TH = H
   20    CONTINUE
C
C        Widen the Gershgorin interval a bit for machines with sloppy
C        arithmetic.
C
         Z = Z + FUDGE*ABS( Z )*DLAMCH( 'Epsilon' )*DBLE( N )
     $         + FUDGE*PIVMIN
      ELSE
         Z = THETA
         Y = ZERO
         NUMZ = NUM
      END IF
C
C     Step 4: Bisection method for finding the upper bound on the L
C             smallest singular values of the bidiagonal.
C             ------------------------------------------
C     A sequence of subintervals [Y,Z] is produced such that
C         (number of S(i) <= Y) < L < (number of S(i) <= Z).
C     NUM : number of S(i) <= TH,
C     NUMZ: number of S(i) <= Z.
C
C     WHILE ( ( NUM .NE. L ) .AND.
C             ( ( Z-Y ) .GT. MAX( TOL, PIVMIN, RELTOL*ABS( Z ) ) ) ) DO
   40 IF ( ( NUM.NE.L ) .AND.
     $     ( ABS( Z-Y ).GT.MAX( TOL, PIVMIN,
     $                          RELTOL*MAX( ABS( Y ), ABS( Z ) ) ) ) )
     $      THEN
         TH = ( Y + Z )/TWO
         NUM = MB03ND( N, TH, Q2, E2, PIVMIN, INFO )
         IF ( NUM.LT.L ) THEN
            Y = TH
         ELSE
            Z = TH
            NUMZ = NUM
         END IF
         GO TO 40
      END IF
C     END WHILE 40
C
C     If NUM <> L and ( Z - Y ) <= TOL, then at least two singular
C     values of J lie in the interval [Y,Z] within a distance less than
C     TOL from each other. S(N-L) and S(N-L+1) are then assumed to
C     coincide. L is increased, and a warning is given.
C
      IF ( NUM.NE.L ) THEN
         L = NUMZ
         THETA = Z
         IWARN = 1
      ELSE
         THETA = TH
      END IF
C
      RETURN
C *** Last line of MB03MD ***
      END
