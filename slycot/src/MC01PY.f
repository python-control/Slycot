      SUBROUTINE MC01PY( K, REZ, IMZ, P, DWORK, INFO )
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
C     To compute the coefficients of a real polynomial P(x) from its
C     zeros. The coefficients are stored in decreasing order of the
C     powers of x.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     K       (input) INTEGER
C             The number of zeros (and hence the degree) of P(x).
C             K >= 0.
C
C     REZ     (input) DOUBLE PRECISION array, dimension (K)
C     IMZ     (input) DOUBLE PRECISION array, dimension (K)
C             The real and imaginary parts of the i-th zero of P(x)
C             must be stored in REZ(i) and IMZ(i), respectively, where
C             i = 1, 2, ..., K. The zeros may be supplied in any order,
C             except that complex conjugate zeros must appear
C             consecutively.
C
C     P       (output) DOUBLE PRECISION array, dimension (K+1)
C             This array contains the coefficients of P(x) in decreasing
C             powers of x.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (K)
C             If K = 0, this array is not referenced.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, (REZ(i),IMZ(i)) is a complex zero but
C                   (REZ(i-1),IMZ(i-1)) is not its conjugate.
C
C     METHOD
C
C     The routine computes the coefficients of the real K-th degree
C     polynomial P(x) as
C
C        P(x) = (x - r(1)) * (x - r(2)) * ... * (x - r(K))
C
C     where r(i) = (REZ(i),IMZ(i)).
C
C     Note that REZ(i) = REZ(j) and IMZ(i) = -IMZ(j) if r(i) and r(j)
C     form a complex conjugate pair (where i <> j), and that IMZ(i) = 0
C     if r(i) is real.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2002.
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, K
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), IMZ(*), P(*), REZ(*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  U, V
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      IF( K.LT.0 ) THEN
         INFO = -1
C
C        Error return.
C
         CALL XERBLA( 'MC01PY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      INFO = 0
      P(1) = ONE
      IF ( K.EQ.0 )
     $   RETURN
C
      I = 1
C     WHILE ( I <= K ) DO
   20 IF ( I.LE.K ) THEN
         U = REZ(I)
         V = IMZ(I)
         DWORK(I) = ZERO
C
         IF ( V.EQ.ZERO ) THEN
            CALL DAXPY( I, -U, P, 1, DWORK, 1 )
C
         ELSE
            IF ( I.EQ.K ) THEN
               INFO = K
               RETURN
            ELSE IF ( ( U.NE.REZ(I+1) ) .OR. ( V.NE.-IMZ(I+1) ) ) THEN
               INFO = I + 1
               RETURN
            END IF
C
            DWORK(I+1) = ZERO
            CALL DAXPY( I, -(U + U),  P, 1, DWORK, 1 )
            CALL DAXPY( I, U**2+V**2, P, 1, DWORK(2), 1 )
            I = I + 1
         END IF
C
         CALL DCOPY( I, DWORK, 1, P(2), 1 )
         I = I + 1
         GO TO 20
      END IF
C     END WHILE 20
C
      RETURN
C *** Last line of MC01PY ***
      END
