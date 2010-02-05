      SUBROUTINE SB08MY( DA, A, B, EPSB )
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
C     To compute the coefficients of B(s) = A(s) * A(-s) and a norm
C     for the accuracy of the computed coefficients.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DA      (input) INTEGER
C             The degree of the polynomials A(s) and B(s).  DA >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (DA+1)
C             This array must contain the coefficients of the polynomial
C             A(s) in increasing powers of s.
C
C     B       (output) DOUBLE PRECISION array, dimension (DA+1)
C             This array contains the coefficients of the polynomial
C             B(s) in increasing powers of s**2.
C
C     EPSB    (input/output) DOUBLE PRECISION
C             On entry, EPSB must contain the machine precision (see
C             LAPACK Library routine DLAMCH).
C             On exit, EPSB contains an updated value, using a norm
C             for the accuracy of the computed coefficients.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB08AZ by A.J. Geurts.
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
      DOUBLE PRECISION  ZERO, ONE, TWO, THREE
      PARAMETER         ( ZERO  = 0.0D0, ONE = 1.0D0, TWO=2.0D0,
     $                    THREE = 3.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DA
      DOUBLE PRECISION  EPSB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*)
C     .. Local Scalars ..
      INTEGER           I, K
      DOUBLE PRECISION  MAXSA, SA, SABS, SIGNI, SIGNK, TERM
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
      SIGNI = ONE
      MAXSA = ZERO
C
      DO 40 I = 0, DA
         SABS = A(I+1)**2
         SA = SIGNI*SABS
         SIGNK = -TWO*SIGNI
C
         DO 20 K = 1, MIN( I, DA - I )
            TERM = SIGNK*A(I-K+1)*A(I+K+1)
            SA = SA + TERM
            SABS = SABS + ABS( TERM )
            SIGNK = -SIGNK
   20    CONTINUE
C
         B(I+1) = SA
         MAXSA = MAX( MAXSA, SABS )
         SIGNI = -SIGNI
   40 CONTINUE
C
      EPSB = THREE*MAXSA*EPSB
C
      RETURN
C *** Last line of SB08MY ***
      END
