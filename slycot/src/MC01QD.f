      SUBROUTINE MC01QD( DA, DB, A, B, RQ, IWARN, INFO )
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
C     To compute, for two given real polynomials A(x) and B(x), the
C     quotient polynomial Q(x) and the remainder polynomial R(x) of
C     A(x) divided by B(x).
C
C     The polynomials Q(x) and R(x) satisfy the relationship
C
C        A(x) = B(x) * Q(x) + R(x),
C
C     where the degree of R(x) is less than the degree of B(x).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DA      (input) INTEGER
C             The degree of the numerator polynomial A(x).  DA >= -1.
C
C     DB      (input/output) INTEGER
C             On entry, the degree of the denominator polynomial B(x).
C             DB >= 0.
C             On exit, if B(DB+1) = 0.0 on entry, then DB contains the
C             index of the highest power of x for which B(DB+1) <> 0.0.
C
C     A       (input) DOUBLE PRECISION array, dimension (DA+1)
C             This array must contain the coefficients of the
C             numerator polynomial A(x) in increasing powers of x
C             unless DA = -1 on entry, in which case A(x) is taken
C             to be the zero polynomial.
C
C     B       (input) DOUBLE PRECISION array, dimension (DB+1)
C             This array must contain the coefficients of the
C             denominator polynomial B(x) in increasing powers of x.
C
C     RQ      (output) DOUBLE PRECISION array, dimension (DA+1)
C             If DA < DB on exit, then this array contains the
C             coefficients of the remainder polynomial R(x) in
C             increasing powers of x; Q(x) is the zero polynomial.
C             Otherwise, the leading DB elements of this array contain
C             the coefficients of R(x) in increasing powers of x, and
C             the next (DA-DB+1) elements contain the coefficients of
C             Q(x) in increasing powers of x.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = k:  if the degree of the denominator polynomial B(x) has
C                   been reduced to (DB - k) because B(DB+1-j) = 0.0 on
C                   entry for j = 0, 1, ..., k-1 and B(DB+1-k) <> 0.0.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if on entry, DB >= 0 and B(i) = 0.0, where
C                   i = 1, 2, ..., DB+1.
C
C     METHOD
C
C     Given real polynomials
C                                                  DA
C        A(x) = a(1) + a(2) * x + ... + a(DA+1) * x
C
C     and
C                                                  DB
C        B(x) = b(1) + b(2) * x + ... + b(DB+1) * x
C
C     where b(DB+1) is non-zero, the routine computes the coeffcients of
C     the quotient polynomial
C                                                     DA-DB
C        Q(x) = q(1) + q(2) * x + ... + q(DA-DB+1) * x
C
C     and the remainder polynomial
C                                                DB-1
C        R(x) = r(1) + r(2) * x + ... + r(DB) * x
C
C     such that A(x) = B(x) * Q(x) + R(x).
C
C     The algorithm used is synthetic division of polynomials (see [1]),
C     which involves the following steps:
C
C        (a) compute q(k+1) = a(DB+k+1) / b(DB+1)
C
C     and
C
C        (b) set a(j) = a(j) - q(k+1) * b(j-k) for j = k+1, ..., DB+k.
C
C     Steps (a) and (b) are performed for k = DA-DB, DA-DB-1, ..., 0 and
C     the algorithm terminates with r(i) = a(i) for i = 1, 2, ..., DB.
C
C     REFERENCES
C
C     [1] Knuth, D.E.
C         The Art of Computer Programming, (Vol. 2, Seminumerical
C         Algorithms).
C         Addison-Wesley, Reading, Massachusetts (2nd Edition), 1981.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01ED by A.J. Geurts.
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
      INTEGER           DA, DB, INFO, IWARN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*), RQ(*)
C     .. Local Scalars ..
      INTEGER           N
      DOUBLE PRECISION  Q
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      IWARN = 0
      INFO = 0
      IF( DA.LT.-1 ) THEN
         INFO = -1
      ELSE IF( DB.LT.0 ) THEN
         INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC01QD', -INFO )
         RETURN
      END IF
C
C     WHILE ( DB >= 0 and B(DB+1) = 0 ) DO
   20 IF ( DB.GE.0 ) THEN
         IF ( B(DB+1).EQ.ZERO ) THEN
            DB = DB - 1
            IWARN = IWARN + 1
            GO TO 20
         END IF
      END IF
C     END WHILE 20
      IF ( DB.EQ.-1 ) THEN
         INFO = 1
         RETURN
      END IF
C
C     B(x) is non-zero.
C
      IF ( DA.GE.0 ) THEN
         N = DA
         CALL DCOPY( N+1, A, 1, RQ, 1 )
C        WHILE ( N >= DB ) DO
   40    IF ( N.GE.DB ) THEN
            IF ( RQ(N+1).NE.ZERO ) THEN
               Q = RQ(N+1)/B(DB+1)
               CALL DAXPY( DB, -Q, B, 1, RQ(N-DB+1), 1 )
               RQ(N+1) = Q
            END IF
            N = N - 1
            GO TO 40
         END IF
C        END WHILE 40
      END IF
C
      RETURN
C *** Last line of MC01QD ***
      END
