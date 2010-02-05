      SUBROUTINE MC01WD( DP, P, U1, U2, Q, INFO )
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
C     To compute, for a given real polynomial P(x) and a quadratic
C     polynomial B(x), the quotient polynomial Q(x) and the linear
C     remainder polynomial R(x) such that
C
C        P(x) = B(x) * Q(x) + R(x),
C
C                                 2
C     where B(x) = u1 + u2 * x + x , R(x) = q(1) + q(2) * (u2 + x)
C     and u1, u2, q(1) and q(2) are real scalars.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DP      (input) INTEGER
C             The degree of the polynomial P(x).  DP >= 0.
C
C     P       (input) DOUBLE PRECISION array, dimension (DP+1)
C             This array must contain the coefficients of P(x) in
C             increasing powers of x.
C
C     U1      (input) DOUBLE PRECISION
C             The value of the constant term of the quadratic
C             polynomial B(x).
C
C     U2      (input) DOUBLE PRECISION
C             The value of the coefficient of x of the quadratic
C             polynomial B(x).
C
C     Q       (output) DOUBLE PRECISION array, dimension (DP+1)
C             If DP >= 1 on entry, then elements Q(1) and Q(2) contain
C             the coefficients q(1) and q(2), respectively, of the
C             remainder polynomial R(x), and the next (DP-1) elements
C             of this array contain the coefficients of the quotient
C             polynomial Q(x) in increasing powers of x.
C             If DP = 0 on entry, then element Q(1) contains the
C             coefficient q(1) of the remainder polynomial R(x) = q(1);
C             Q(x) is the zero polynomial.
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
C     Given the real polynomials
C
C                DP           i                           2
C        P(x) = SUM p(i+1) * x  and B(x) = u1 + u2 * x + x
C               i=0
C
C     the routine uses the recurrence relationships
C
C        q(DP+1) = p(DP+1),
C
C        q(DP) = p(DP) - u2 * q(DP+1) and
C
C        q(i)  = p(i) - u2 * q(i+1) - u1 * q(i+2) for i = DP-1, ..., 1
C
C     to determine the coefficients of the quotient polynomial
C
C               DP-2          i
C        Q(x) = SUM q(i+3) * x
C               i=0
C
C     and the remainder polynomial
C
C        R(x) = q(1) + q(2) * (u2 + x).
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01KD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, polynomial operations,
C     quadratic polynomial.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           DP, INFO
      DOUBLE PRECISION  U1, U2
C     .. Array Arguments ..
      DOUBLE PRECISION  P(*), Q(*)
C     .. Local Scalars ..
      INTEGER           I, N
      DOUBLE PRECISION  A, B, C
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      IF ( DP.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'MC01WD', -INFO )
         RETURN
      END IF
C
      INFO = 0
      N = DP + 1
      Q(N) = P(N)
      IF ( N.GT.1 ) THEN
         B = Q(N)
         Q(N-1) = P(N-1) - U2*B
         IF ( N.GT.2 ) THEN
            A = Q(N-1)
C
            DO 20 I = N - 2, 1, -1
               C = P(I) - U2*A - U1*B
               Q(I) = C
               B = A
               A = C
   20       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of MC01WD ***
      END
