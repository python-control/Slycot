      SUBROUTINE MB02SZ( N, H, LDH, IPIV, INFO )
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
C     To compute an LU factorization of a complex n-by-n upper
C     Hessenberg matrix H using partial pivoting with row interchanges.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix H.  N >= 0.
C
C     H       (input/output) COMPLEX*16 array, dimension (LDH,N)
C             On entry, the n-by-n upper Hessenberg matrix to be
C             factored.
C             On exit, the factors L and U from the factorization
C             H = P*L*U; the unit diagonal elements of L are not stored,
C             and L is lower bidiagonal.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= max(1,N).
C
C     IPIV    (output) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= i <= N, row i of the matrix
C             was interchanged with row IPIV(i).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, U(i,i) is exactly zero. The
C                   factorization has been completed, but the factor U
C                   is exactly singular, and division by zero will occur
C                   if it is used to solve a system of equations.
C
C     METHOD
C
C     The factorization has the form
C        H = P * L * U
C     where P is a permutation matrix, L is lower triangular with unit
C     diagonal elements (and one nonzero subdiagonal), and U is upper
C     triangular.
C
C     This is the right-looking Level 2 BLAS version of the algorithm
C     (adapted after ZGETF2).
C
C     REFERENCES
C
C     -
C
C     NUMERICAL ASPECTS
C                                2
C     The algorithm requires 0( N ) complex operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TB01FX by A.J. Laub, University of
C     Southern California, United States of America, May 1980.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000,
C     Jan. 2005.
C
C     KEYWORDS
C
C     Frequency response, Hessenberg form, matrix algebra.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, N
C     .. Array Arguments ..
      INTEGER           IPIV(*)
      COMPLEX*16        H(LDH,*)
C     .. Local Scalars ..
      INTEGER           J, JP
C     .. External Functions ..
      DOUBLE PRECISION  DCABS1
      EXTERNAL          DCABS1
C     .. External Subroutines ..
      EXTERNAL          XERBLA, ZAXPY, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02SZ', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      DO 10 J = 1, N
C
C        Find pivot and test for singularity.
C
         JP = J
         IF ( J.LT.N ) THEN
            IF ( DCABS1( H( J+1, J ) ).GT.DCABS1( H( J, J ) ) )
     $         JP = J + 1
         END IF
         IPIV( J ) = JP
         IF( H( JP, J ).NE.ZERO ) THEN
C
C           Apply the interchange to columns J:N.
C
            IF( JP.NE.J )
     $         CALL ZSWAP( N-J+1, H( J, J ), LDH, H( JP, J ), LDH )
C
C           Compute element J+1 of J-th column.
C
            IF( J.LT.N )
     $         H( J+1, J ) = H( J+1, J )/H( J, J )
C
         ELSE IF( INFO.EQ.0 ) THEN
C
            INFO = J
         END IF
C
         IF( J.LT.N ) THEN
C
C           Update trailing submatrix.
C
            CALL ZAXPY( N-J, -H( J+1, J ), H( J, J+1 ), LDH,
     $                  H( J+1, J+1 ), LDH )
         END IF
   10 CONTINUE
      RETURN
C *** Last line of MB02SZ ***
      END
