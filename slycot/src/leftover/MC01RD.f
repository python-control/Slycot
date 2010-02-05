      SUBROUTINE MC01RD( DP1, DP2, DP3, ALPHA, P1, P2, P3, INFO )
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
C     To compute the coefficients of the polynomial
C
C        P(x) = P1(x) * P2(x) + alpha * P3(x),
C
C     where P1(x), P2(x) and P3(x) are given real polynomials and alpha
C     is a real scalar.
C
C     Each of the polynomials P1(x), P2(x) and P3(x) may be the zero
C     polynomial.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DP1     (input) INTEGER
C             The degree of the polynomial P1(x).  DP1 >= -1.
C
C     DP2     (input) INTEGER
C             The degree of the polynomial P2(x).  DP2 >= -1.
C
C     DP3     (input/output) INTEGER
C             On entry, the degree of the polynomial P3(x).  DP3 >= -1.
C             On exit, the degree of the polynomial P(x).
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar value alpha of the problem.
C
C     P1      (input) DOUBLE PRECISION array, dimension (lenp1)
C             where lenp1 = DP1 + 1 if DP1 >= 0 and 1 otherwise.
C             If DP1 >= 0, then this array must contain the
C             coefficients of P1(x) in increasing powers of x.
C             If DP1 = -1, then P1(x) is taken to be the zero
C             polynomial, P1 is not referenced and can be supplied
C             as a dummy array.
C
C     P2      (input) DOUBLE PRECISION array, dimension (lenp2)
C             where lenp2 = DP2 + 1 if DP2 >= 0 and 1 otherwise.
C             If DP2 >= 0, then this array must contain the
C             coefficients of P2(x) in increasing powers of x.
C             If DP2 = -1, then P2(x) is taken to be the zero
C             polynomial, P2 is not referenced and can be supplied
C             as a dummy array.
C
C     P3      (input/output) DOUBLE PRECISION array, dimension (lenp3)
C             where lenp3 = MAX(DP1+DP2,DP3,0) + 1.
C             On entry, if DP3 >= 0, then this array must contain the
C             coefficients of P3(x) in increasing powers of x.
C             On entry, if DP3 = -1, then P3(x) is taken to be the zero
C             polynomial.
C             On exit, the leading (DP3+1) elements of this array
C             contain the coefficients of P(x) in increasing powers of x
C             unless DP3 = -1 on exit, in which case the coefficients of
C             P(x) (the zero polynomial) are not stored in the array.
C             This is the case, for instance, when ALPHA = 0.0 and
C             P1(x) or P2(x) is the zero polynomial.
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
C     Given real polynomials
C
C                DP1           i           DP2           i
C        P1(x) = SUM a(i+1) * x ,  P2(x) = SUM b(i+1) * x  and
C                i=0                       i=0
C
C                DP3           i
C        P3(x) = SUM c(i+1) * x ,
C                i=0
C
C     the routine computes the coefficents of P(x) = P1(x) * P2(x) +
C                     DP3            i
C     alpha * P3(x) = SUM  d(i+1) * x  as follows.
C                     i=0
C
C     Let e(i) = c(i) for 1 <= i <= DP3+1 and e(i) = 0 for i > DP3+1.
C     Then if DP1 >= DP2,
C
C                i
C        d(i) = SUM a(k) * b(i-k+1) + f(i), for i = 1, ..., DP2+1,
C               k=1
C
C                 i
C        d(i)  = SUM a(k) * b(i-k+1) + f(i), for i = DP2+2, ..., DP1+1
C               k=i-DP2
C
C     and
C                DP1+1
C        d(i)  = SUM a(k) * b(i-k+1) + f(i) for i = DP1+2,...,DP1+DP2+1,
C               k=i-DP2
C
C     where f(i) = alpha * e(i).
C
C     Similar formulas hold for the case DP1 < DP2.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01FD by C. Klimann and
C     A.J. Geurts.
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
      INTEGER           DP1, DP2, DP3, INFO
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  P1(*), P2(*), P3(*)
C     .. Local Scalars ..
      INTEGER           D1, D2, D3, DMAX, DMIN, DSUM, E3, I, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( DP1.LT.-1 ) THEN
         INFO = -1
      ELSE IF( DP2.LT.-1 ) THEN
         INFO = -2
      ELSE IF( DP3.LT.-1 ) THEN
         INFO = -3
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC01RD', -INFO )
         RETURN
      END IF
C
C     Computation of the exact degree of the polynomials, i.e., Di such
C     that either Di = -1 or Pi(Di+1) is non-zero.
C
      D1 = DP1
C     WHILE ( D1 >= 0 and P1(D1+1) = 0 ) DO
   20 IF ( D1.GE.0 ) THEN
         IF ( P1(D1+1).EQ.ZERO ) THEN
            D1 = D1 - 1
            GO TO 20
         END IF
      END IF
C     END WHILE 20
      D2 = DP2
C     WHILE ( D2 >= 0 and P2(D2+1) = 0 ) DO
   40 IF ( D2.GE.0 ) THEN
         IF ( P2(D2+1).EQ.ZERO ) THEN
            D2 = D2 - 1
            GO TO 40
         END IF
      END IF
C     END WHILE 40
      IF ( ALPHA.EQ.ZERO ) THEN
         D3 = -1
      ELSE
         D3 = DP3
      END IF
C     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO
   60 IF ( D3.GE.0 ) THEN
         IF ( P3(D3+1).EQ.ZERO ) THEN
            D3 = D3 - 1
            GO TO 60
         END IF
      END IF
C     END WHILE 60
C
C     Computation of P3(x) := ALPHA * P3(x).
C
      CALL DSCAL( D3+1, ALPHA, P3, 1 )
C
      IF ( ( D1.EQ.-1 ) .OR. ( D2.EQ.-1 ) ) THEN
         DP3 = D3
         RETURN
      END IF
C
C     P1(x) and P2(x) are non-zero polynomials.
C
      DSUM = D1 + D2
      DMAX = MAX( D1, D2 )
      DMIN = DSUM - DMAX
C
      IF ( D3.LT.DSUM ) THEN
         P3(D3+2) = ZERO
         CALL DCOPY( DSUM-D3-1, P3(D3+2), 0, P3(D3+3), 1 )
         D3 = DSUM
      END IF
C
      IF ( ( D1.EQ.0 ) .OR. ( D2.EQ.0 ) ) THEN
C
C        D1 or D2 is zero.
C
         IF ( D1.NE.0 ) THEN
            CALL DAXPY( D1+1, P2(1), P1, 1, P3, 1 )
         ELSE
            CALL DAXPY( D2+1, P1(1), P2, 1, P3, 1 )
         END IF
      ELSE
C
C        D1 and D2 are both nonzero.
C
C        First part of the computation.
C
         DO 80 I = 1,  DMIN + 1
            P3(I) = P3(I) + DDOT( I, P1, 1, P2, -1 )
   80    CONTINUE
C
C        Second part of the computation.
C
         DO 100 I = DMIN + 2, DMAX + 1
            IF ( D1.GT.D2 ) THEN
               K = I - D2
               P3(I) = P3(I) + DDOT( DMIN+1, P1(K), 1, P2, -1 )
            ELSE
               K = I - D1
               P3(I) = P3(I) + DDOT( DMIN+1, P2(K), -1, P1, 1 )
            END IF
  100    CONTINUE
C
C        Third part of the computation.
C
         E3 = DSUM + 2
C
         DO 120 I = DMAX + 2, DSUM + 1
            J = E3 - I
            K = I - DMIN
            L = I - DMAX
            IF ( D1.GT.D2 ) THEN
               P3(I) = P3(I) + DDOT( J, P1(K), 1, P2(L), -1 )
            ELSE
               P3(I) = P3(I) + DDOT( J, P1(L), -1, P2(K), 1 )
            END IF
  120    CONTINUE
C
      END IF
C
C     Computation of the exact degree of P3(x).
C
C     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO
  140 IF ( D3.GE.0 ) THEN
         IF ( P3(D3+1).EQ.ZERO ) THEN
            D3 = D3 - 1
            GO TO 140
         END IF
      END IF
C     END WHILE 140
      DP3 = D3
C
      RETURN
C *** Last line of MC01RD ***
      END
