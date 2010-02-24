      SUBROUTINE SB08NY( DA, A, B, EPSB )
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
C     To compute the coefficients of B(z) = A(1/z) * A(z) and a norm for
C     the accuracy of the computed coefficients.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DA      (input) INTEGER
C             The degree of the polynomials A(z) and B(z).  DA >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (DA+1)
C             This array must contain the coefficients of the polynomial
C             A(z) in increasing powers of z.
C
C     B       (output) DOUBLE PRECISION array, dimension (DA+1)
C             This array contains the coefficients of the polynomial
C             B(z).
C
C     EPSB    (output) DOUBLE PRECISION
C             A value used for checking the accuracy of the computed
C             coefficients.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB08BZ by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Laplace transform, polynomial operations, spectral factorization.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  THREE
      PARAMETER         ( THREE = 3.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DA
      DOUBLE PRECISION  EPSB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH
      EXTERNAL          DDOT, DLAMCH
C     .. Executable Statements ..
C
      DO 20 I = 1, DA + 1
         B(I) = DDOT( DA-I+2, A(1), 1, A(I), 1 )
   20 CONTINUE
C
      EPSB = THREE*DLAMCH( 'Epsilon' )*B(1)
C
      RETURN
C *** Last line of SB08NY ***
      END
