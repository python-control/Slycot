      SUBROUTINE MC01MD( DP, ALPHA, K, P, Q, INFO )
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
C     To calculate, for a given real polynomial P(x) and a real scalar
C     alpha, the leading K coefficients of the shifted polynomial
C                                                               K-1
C        P(x) = q(1) + q(2) * (x-alpha) + ... + q(K) * (x-alpha)   + ...
C
C     using Horner's algorithm.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DP      (input) INTEGER
C             The degree of the polynomial P(x).  DP >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar value alpha of the problem.
C
C     K       (input) INTEGER
C             The number of coefficients of the shifted polynomial to be
C             computed.  1 <= K <= DP+1.
C
C     P       (input) DOUBLE PRECISION array, dimension (DP+1)
C             This array must contain the coefficients of P(x) in
C             increasing powers of x.
C
C     Q       (output) DOUBLE PRECISION array, dimension (DP+1)
C             The leading K elements of this array contain the first
C             K coefficients of the shifted polynomial in increasing
C             powers of (x - alpha), and the next (DP-K+1) elements
C             are used as internal workspace.
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
C     Given the real polynomial
C                                         2                    DP
C        P(x) = p(1) + p(2) * x + p(3) * x  + ... + p(DP+1) * x  ,
C
C     the routine computes the leading K coefficients of the shifted
C     polynomial
C                                                                   K-1
C        P(x) = q(1) + q(2) * (x - alpha) + ... + q(K) * (x - alpha)
C
C     as follows.
C
C     Applying Horner's algorithm (see [1]) to P(x), i.e. dividing P(x)
C     by (x-alpha), yields
C
C        P(x) = q(1) + (x-alpha) * D(x),
C
C     where q(1) is the value of the constant term of the shifted
C     polynomial and D(x) is the quotient polynomial of degree (DP-1)
C     given by
C                                         2                     DP-1
C        D(x) = d(2) + d(3) * x + d(4) * x  + ... +  d(DP+1) * x    .
C
C     Applying Horner's algorithm to D(x) and subsequent quotient
C     polynomials yields q(2) and q(3), q(4), ..., q(K) respectively.
C
C     It follows immediately that q(1) = P(alpha), and in general
C                (i-1)
C        q(i) = P     (alpha) / (i - 1)! for i = 1, 2, ..., K.
C
C     REFERENCES
C
C     [1] STOER, J. and BULIRSCH, R.
C         Introduction to Numerical Analysis.
C         Springer-Verlag. 1980.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01AD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, polynomial operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DP, INFO, K
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  P(*), Q(*)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          DCOPY, XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( DP.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LE.0 .OR. K.GT.DP+1 ) THEN
         INFO = -3
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC01MD', -INFO )
         RETURN
      END IF
C
      CALL DCOPY( DP+1, P, 1, Q, 1 )
      IF ( DP.EQ.0 .OR. ALPHA.EQ.ZERO )
     $   RETURN
C
      DO 40 J = 1, K
C
         DO 20 I = DP, J, -1
            Q(I) = Q(I) + ALPHA*Q(I+1)
   20    CONTINUE
C
   40 CONTINUE
C
      RETURN
C *** Last line of MC01MD ***
      END
