      SUBROUTINE MA02GD( N, A, LDA, K1, K2, IPIV, INCX )
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
C     To perform a series of column interchanges on the matrix A.
C     One column interchange is initiated for each of columns K1 through
C     K2 of A. This is useful for solving linear systems X*A = B, when
C     the matrix A has already been factored by LAPACK Library routine
C     DGETRF.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of rows of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,*)
C             On entry, the leading N-by-M part of this array must
C             contain the matrix A to which the column interchanges will
C             be applied, where M is the largest element of IPIV(K), for
C             K = K1, ..., K2.
C             On exit, the leading N-by-M part of this array contains
C             the permuted matrix.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     K1      (input) INTEGER
C             The first element of IPIV for which a column interchange
C             will be done.
C
C     K2      (input) INTEGER
C             The last element of IPIV for which a column interchange
C             will be done.
C
C     IPIV    (input) INTEGER array, dimension (K1+(K2-K1)*abs(INCX))
C             The vector of interchanging (pivot) indices.  Only the
C             elements in positions K1 through K2 of IPIV are accessed.
C             IPIV(K) = L implies columns K and L are to be
C             interchanged.
C
C     INCX    (input) INTEGER
C             The increment between successive values of IPIV.
C             If INCX is negative, the interchanges are applied in
C             reverse order.
C
C     METHOD
C
C     The columns IPIV(K) and K are swapped for K = K1, ..., K2, for
C     INCX = 1 (and similarly, for INCX <> 1).
C
C     FURTHER COMMENTS
C
C     This routine is the column-oriented counterpart of the LAPACK
C     Library routine DLASWP. The LAPACK Library routine DLAPMT cannot
C     be used in this context. To solve the system X*A = B, where A and
C     B are N-by-N and M-by-N, respectively, the following statements
C     can be used:
C
C         CALL DGETRF( N, N, A, LDA, IPIV, INFO )
C         CALL DTRSM( 'R', 'U', 'N', 'N', M, N, ONE, A, LDA, B, LDB )
C         CALL DTRSM( 'R', 'L', 'N', 'U', M, N, ONE, A, LDA, B, LDB )
C         CALL MA02GD( M, B, LDB, 1, N, IPIV, -1 )
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2008.
C
C     KEYWORDS
C
C     Elementary matrix operations, linear algebra.
C
C    ******************************************************************
C
C      .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
C      ..
C      .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
C     ..
C     .. Local Scalars ..
      INTEGER            J, JP, JX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSWAP
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF( INCX.EQ.0 .OR. N.EQ.0 )
     $   RETURN
C
C     Interchange column J with column IPIV(J) for each of columns K1
C     through K2.
C
      IF( INCX.GT.0 ) THEN
         JX = K1
      ELSE
         JX = 1 + ( 1-K2 )*INCX
      END IF
C
      IF( INCX.EQ.1 ) THEN
C
         DO 10 J = K1, K2
            JP = IPIV( J )
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   10    CONTINUE
C
      ELSE IF( INCX.GT.1 ) THEN
C
         DO 20 J = K1, K2
            JP = IPIV( JX )
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
            JX = JX + INCX
   20    CONTINUE
C
      ELSE IF( INCX.LT.0 ) THEN
C
         DO 30 J = K2, K1, -1
            JP = IPIV( JX )
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
            JX = JX + INCX
   30    CONTINUE
C
      END IF
C
      RETURN
C
C *** Last line of MA02GD ***
      END
