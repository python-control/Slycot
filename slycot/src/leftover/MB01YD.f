      SUBROUTINE MB01YD( UPLO, TRANS, N, K, L, ALPHA, BETA, A, LDA, C,
     $                   LDC, INFO )
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
C     To perform the symmetric rank k operations
C
C        C := alpha*op( A )*op( A )' + beta*C,
C
C     where alpha and beta are scalars, C is an n-by-n symmetric matrix,
C     op( A ) is an n-by-k matrix, and op( A ) is one of
C
C        op( A ) = A   or   op( A ) = A'.
C
C     The matrix A has l nonzero codiagonals, either upper or lower.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the symmetric matrix C
C             is given and computed, as follows:
C             = 'U':  the upper triangular part is given/computed;
C             = 'L':  the lower triangular part is given/computed.
C             UPLO also defines the pattern of the matrix A (see below).
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( A ) to be used, as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix C.  N >= 0.
C
C     K       (input) INTEGER
C             The number of columns of the matrix op( A ).  K >= 0.
C
C     L       (input) INTEGER
C             If UPLO = 'U', matrix A has L nonzero subdiagonals.
C             If UPLO = 'L', matrix A has L nonzero superdiagonals.
C             MAX(0,NR-1) >= L >= 0, if UPLO = 'U',
C             MAX(0,NC-1) >= L >= 0, if UPLO = 'L',
C             where NR and NC are the numbers of rows and columns of the
C             matrix A, respectively.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then the array A is
C             not referenced.
C
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. When beta is zero then the array C need
C             not be set before entry.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,NC), where
C             NC is K when TRANS = 'N', and is N otherwise.
C             If TRANS = 'N', the leading N-by-K part of this array must
C             contain the matrix A, otherwise the leading K-by-N part of
C             this array must contain the matrix A.
C             If UPLO = 'U', only the upper triangular part and the
C             first L subdiagonals are referenced, and the remaining
C             subdiagonals are assumed to be zero.
C             If UPLO = 'L', only the lower triangular part and the
C             first L superdiagonals are referenced, and the remaining
C             superdiagonals are assumed to be zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,NR),
C             where NR = N, if TRANS = 'N', and NR = K, otherwise.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry with UPLO = 'U', the leading N-by-N upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix C.
C             On entry with UPLO = 'L', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix C.
C             On exit, the leading N-by-N upper triangular part (if
C             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of
C             this array contains the corresponding triangular part of
C             the updated matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,N).
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
C     The calculations are efficiently performed taking the symmetry
C     and structure into account.
C
C     FURTHER COMMENTS
C
C     The matrix A may have the following patterns, when n = 7, k = 5,
C     and l = 2 are used for illustration:
C
C     UPLO = 'U', TRANS = 'N'         UPLO = 'L', TRANS = 'N'
C
C            [ x x x x x ]                   [ x x x 0 0 ]
C            [ x x x x x ]                   [ x x x x 0 ]
C            [ x x x x x ]                   [ x x x x x ]
C        A = [ 0 x x x x ],              A = [ x x x x x ],
C            [ 0 0 x x x ]                   [ x x x x x ]
C            [ 0 0 0 x x ]                   [ x x x x x ]
C            [ 0 0 0 0 x ]                   [ x x x x x ]
C
C     UPLO = 'U', TRANS = 'T'         UPLO = 'L', TRANS = 'T'
C
C            [ x x x x x x x ]               [ x x x 0 0 0 0 ]
C            [ x x x x x x x ]               [ x x x x 0 0 0 ]
C        A = [ x x x x x x x ],          A = [ x x x x x 0 0 ].
C            [ 0 x x x x x x ]               [ x x x x x x 0 ]
C            [ 0 0 x x x x x ]               [ x x x x x x x ]
C
C     If N = K, the matrix A is upper or lower triangular, for L = 0,
C     and upper or lower Hessenberg, for L = 1.
C
C     This routine is a specialization of the BLAS 3 routine DSYRK.
C     BLAS 1 calls are used when appropriate, instead of in-line code,
C     in order to increase the efficiency. If the matrix A is full, or
C     its zero triangle has small order, an optimized DSYRK code could
C     be faster than MB01YD.
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          TRANS, UPLO
      INTEGER            INFO, LDA, LDC, K, L, N
      DOUBLE PRECISION   ALPHA, BETA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            TRANSP, UPPER
      INTEGER            I, J, M, NCOLA, NROWA
      DOUBLE PRECISION   TEMP
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DLASCL, DLASET, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO   = 0
      UPPER  = LSAME( UPLO,  'U' )
      TRANSP = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      IF( TRANSP )THEN
         NROWA = K
         NCOLA = N
      ELSE
         NROWA = N
         NCOLA = K
      END IF
C
      IF( UPPER )THEN
         M = NROWA
      ELSE
         M = NCOLA
      END IF
C
      IF(      .NOT.( UPPER  .OR. LSAME( UPLO,  'L' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( TRANSP .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( L.LT.0 .OR. L.GT.MAX( 0, M-1 ) ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01YD', -INFO )
         RETURN
      END IF
C
C     Quick return, if possible.
C
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
      IF ( ALPHA.EQ.ZERO ) THEN
         IF ( BETA.EQ.ZERO ) THEN
C
C           Special case when both alpha = 0 and beta = 0.
C
            CALL DLASET( UPLO, N, N, ZERO, ZERO, C, LDC )
         ELSE
C
C           Special case alpha = 0.
C
            CALL DLASCL( UPLO, 0, 0, ONE, BETA, N, N, C, LDC, INFO )
         END IF
         RETURN
      END IF
C
C     General case: alpha <> 0.
C
      IF ( .NOT.TRANSP ) THEN
C
C        Form  C := alpha*A*A' + beta*C.
C
         IF ( UPPER ) THEN
C
            DO 30 J = 1, N
               IF ( BETA.EQ.ZERO ) THEN
C
                  DO 10 I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
C
               ELSE IF ( BETA.NE.ONE ) THEN
                  CALL DSCAL ( J, BETA, C( 1, J ), 1 )
               END IF
C
               DO 20 M = MAX( 1, J-L ), K
                  CALL DAXPY ( MIN( J, L+M ), ALPHA*A( J, M ),
     $                         A( 1, M ), 1, C( 1, J ), 1 )
   20          CONTINUE
C
   30       CONTINUE
C
         ELSE
C
            DO 60 J = 1, N
               IF ( BETA.EQ.ZERO ) THEN
C
                  DO 40 I = J, N
                     C( I, J ) = ZERO
   40             CONTINUE
C
               ELSE IF ( BETA.NE.ONE ) THEN
                  CALL DSCAL ( N-J+1, BETA, C( J, J ), 1 )
               END IF
C
               DO 50 M = 1, MIN( J+L, K )
                  CALL DAXPY ( N-J+1, ALPHA*A( J, M ), A( J, M ), 1,
     $                         C( J, J ), 1 )
   50          CONTINUE
C
   60       CONTINUE
C
         END IF
C
      ELSE
C
C        Form  C := alpha*A'*A + beta*C.
C
         IF ( UPPER ) THEN
C
            DO 80 J = 1, N
C
               DO 70 I = 1, J
                  TEMP = ALPHA*DDOT ( MIN( J+L, K ), A( 1, I ), 1,
     $                                A( 1, J ), 1 )
                  IF ( BETA.EQ.ZERO ) THEN
                     C( I, J ) = TEMP
                  ELSE
                     C( I, J ) = TEMP + BETA*C( I, J )
                  END IF
   70          CONTINUE
C
   80       CONTINUE
C
         ELSE
C
            DO 100 J = 1, N
C
               DO 90 I = J, N
                  M = MAX( 1, I-L )
                  TEMP = ALPHA*DDOT ( K-M+1, A( M, I ), 1, A( M, J ),
     $                                1 )
                  IF ( BETA.EQ.ZERO ) THEN
                     C( I, J ) = TEMP
                  ELSE
                     C( I, J ) = TEMP + BETA*C( I, J )
                  END IF
   90          CONTINUE
C
  100       CONTINUE
C
         END IF
C
      END IF
C
      RETURN
C
C *** Last line of MB01YD ***
      END
