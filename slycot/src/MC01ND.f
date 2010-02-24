      SUBROUTINE MC01ND( DP, XR, XI, P, VR, VI, INFO )
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
C     To compute the value of the real polynomial P(x) at a given
C     complex point x = x0 using Horner's algorithm.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DP      (input) INTEGER
C             The degree of the polynomial P(x).  DP >= 0.
C
C     XR      (input) DOUBLE PRECISION
C     XI      (input) DOUBLE PRECISION
C             The real and imaginary parts, respectively, of x0.
C
C     P       (input) DOUBLE PRECISION array, dimension (DP+1)
C             This array must contain the coefficients of the polynomial
C             P(x) in increasing powers of x.
C
C     VR      (output) DOUBLE PRECISION
C     VI      (output) DOUBLE PRECISION
C             The real and imaginary parts, respectively, of P(x0).
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
C                                         2                   DP
C        P(x) = p(1) + p(2) * x + p(3) * x + ... + p(DP+1) * x  ,
C
C     the routine computes the value of P(x0) using the recursion
C
C        q(DP+1) = p(DP+1),
C        q(i) = x0*q(i+1) + p(i) for i = DP, DP-1, ..., 1,
C
C     which is known as Horner's algorithm (see [1]). Then q(1) = P(x0).
C
C     REFERENCES
C
C     [1] STOER, J and BULIRSCH, R.
C         Introduction to Numerical Analysis.
C         Springer-Verlag. 1980.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires DP operations for real arguments and 4*DP
C     for complex arguments.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01BD by Serge Steer.
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
      INTEGER           DP, INFO
      DOUBLE PRECISION  VI, VR, XI, XR
C     .. Array Arguments ..
      DOUBLE PRECISION  P(*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  T
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      IF( DP.LT.0 ) THEN
         INFO = -1
C
C        Error return.
C
         CALL XERBLA( 'MC01ND', -INFO )
         RETURN
      END IF
C
      INFO = 0
      VR = P(DP+1)
      VI = ZERO
C
      IF ( DP.EQ.0 )
     $   RETURN
C
      IF ( XI.EQ.ZERO ) THEN
C
C        X real.
C
         DO 20 I = DP, 1, -1
            VR = VR*XR + P(I)
   20    CONTINUE
C
      ELSE
C
C        X complex.
C
         DO 40 I = DP, 1, -1
            T  = VR*XR - VI*XI + P(I)
            VI = VI*XR + VR*XI
            VR = T
   40    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MC01ND ***
      END
