      SUBROUTINE MB01XY( UPLO, N, A, LDA, INFO )
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
C     To compute the matrix product U' * U or L * L', where U and L are
C     upper and lower triangular matrices, respectively, stored in the
C     corresponding upper or lower triangular part of the array A.
C
C     If UPLO = 'U' then the upper triangle of the result is stored,
C     overwriting the matrix U in A.
C     If UPLO = 'L' then the lower triangle of the result is stored,
C     overwriting the matrix L in A.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies which triangle (U or L) is given in the array A,
C             as follows:
C             = 'U':  the upper triangular part U is given;
C             = 'L':  the lower triangular part L is given.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the triangular matrices U or L.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, if UPLO = 'U', the leading N-by-N upper
C             triangular part of this array must contain the upper
C             triangular matrix U.
C             On entry, if UPLO = 'L', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangular matrix L.
C             On exit, if UPLO = 'U', the leading N-by-N upper
C             triangular part of this array contains the upper
C             triangular part of the product U' * U. The strictly lower
C             triangular part is not referenced.
C             On exit, if UPLO = 'L', the leading N-by-N lower
C             triangular part of this array contains the lower
C             triangular part of the product L * L'. The strictly upper
C             triangular part is not referenced.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,N).
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
C     The matrix product U' * U or L * L' is computed using BLAS 2 and
C     BLAS 1 operations (an unblocked algorithm).
C
C     FURTHER COMMENTS
C
C     This routine is a counterpart of LAPACK Library routine DLAUU2,
C     which computes the matrix product U * U' or L' * L.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   AII
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO  = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01XY', -INFO )
         RETURN
      END IF
C
C     Quick return, if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( UPPER ) THEN
C
C        Compute the product U' * U.
C
         A( N, N ) = DDOT( N, A( 1, N ), 1, A( 1, N ), 1 )
C
         DO 10 I = N-1, 2, -1
            AII = A( I, I )
            A( I, I ) = DDOT( I, A( 1, I ), 1, A( 1, I ), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA,
     $                  A( 1, I ), 1, AII, A( I, I+1 ), LDA )
   10    CONTINUE
C
         IF( N.GT.1 ) THEN
            AII = A( 1, 1 )
            CALL DSCAL( N, AII, A( 1, 1 ), LDA )
         END IF
C
      ELSE
C
C        Compute the product L * L'.
C
         A( N, N ) = DDOT( N, A( N, 1 ), LDA, A( N, 1 ), LDA )
C
         DO 20 I = N-1, 2, -1
            AII = A( I, I )
            A( I, I ) = DDOT( I, A( I, 1 ), LDA, A( I, 1 ), LDA )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, A( I+1, 1 ),
     $                  LDA, A( I, 1 ), LDA, AII, A( I+1, I ), 1 )
   20    CONTINUE
C
         IF( N.GT.1 ) THEN
            AII = A( 1, 1 )
            CALL DSCAL( N, AII, A( 1, 1 ), 1 )
         END IF
      END IF
C
      RETURN
C
C *** Last line of MB01XY ***
      END
