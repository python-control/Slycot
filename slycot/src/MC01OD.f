      SUBROUTINE MC01OD( K, REZ, IMZ, REP, IMP, DWORK, INFO )
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
C     To compute the coefficients of a complex polynomial P(x) from its
C     zeros.
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
C             i = 1, 2, ..., K. The zeros may be supplied in any order.
C
C     REP     (output) DOUBLE PRECISION array, dimension (K+1)
C     IMP     (output) DOUBLE PRECISION array, dimension (K+1)
C             These arrays contain the real and imaginary parts,
C             respectively, of the coefficients of P(x) in increasing
C             powers of x. If K = 0, then REP(1) is set to one and
C             IMP(1) is set to zero.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (2*K+2)
C             If K = 0, this array is not referenced.
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
C     The routine computes the coefficients of the complex K-th degree
C     polynomial P(x) as
C
C        P(x) = (x - r(1)) * (x - r(2)) * ... * (x - r(K))
C
C     where r(i) = (REZ(i),IMZ(i)), using real arithmetic.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01CD by Alan Brown and
C     A.J. Geurts.
C
C     REVISIONS
C
C     V. Sima, May 2002.
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
      DOUBLE PRECISION  DWORK(*), IMP(*), IMZ(*), REP(*), REZ(*)
C     .. Local Scalars ..
      INTEGER           I, K2
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
         CALL XERBLA( 'MC01OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      INFO   = 0
      REP(1) = ONE
      IMP(1) = ZERO
      IF ( K.EQ.0 )
     $   RETURN
C
      K2 = K + 2
C
      DO 20 I = 1, K
         U = REZ(I)
         V = IMZ(I)
         DWORK(1)  = ZERO
         DWORK(K2) = ZERO
         CALL DCOPY( I, REP, 1, DWORK(2), 1 )
         CALL DCOPY( I, IMP, 1, DWORK(K2+1), 1 )
C
         IF ( U.NE.ZERO ) THEN
            CALL DAXPY( I, -U, REP, 1, DWORK, 1 )
            CALL DAXPY( I, -U, IMP, 1, DWORK(K2), 1 )
         END IF
C
         IF ( V.NE.ZERO ) THEN
            CALL DAXPY( I,  V, IMP, 1, DWORK, 1 )
            CALL DAXPY( I, -V, REP, 1, DWORK(K2), 1 )
         END IF
C
         CALL DCOPY( I+1, DWORK, 1, REP, 1 )
         CALL DCOPY( I+1, DWORK(K2), 1, IMP, 1 )
   20 CONTINUE
C
      RETURN
C *** Last line of MC01OD ***
      END
