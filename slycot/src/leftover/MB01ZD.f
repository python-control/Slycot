      SUBROUTINE MB01ZD( SIDE, UPLO, TRANST, DIAG, M, N, L, ALPHA, T,
     $                   LDT, H, LDH, INFO )
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
C     To compute the matrix product
C
C        H := alpha*op( T )*H,   or   H := alpha*H*op( T ),
C
C     where alpha is a scalar, H is an m-by-n upper or lower
C     Hessenberg-like matrix (with l nonzero subdiagonals or
C     superdiagonals, respectively), T is a unit, or non-unit,
C     upper or lower triangular matrix, and op( T ) is one of
C
C        op( T ) = T   or   op( T ) = T'.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies whether the triangular matrix T appears on the
C             left or right in the matrix product, as follows:
C             = 'L':  the product alpha*op( T )*H is computed;
C             = 'R':  the product alpha*H*op( T ) is computed.
C
C     UPLO    CHARACTER*1
C             Specifies the form of the matrices T and H, as follows:
C             = 'U':  the matrix T is upper triangular and the matrix H
C                     is upper Hessenberg-like;
C             = 'L':  the matrix T is lower triangular and the matrix H
C                     is lower Hessenberg-like.
C
C     TRANST  CHARACTER*1
C             Specifies the form of op( T ) to be used, as follows:
C             = 'N':  op( T ) = T;
C             = 'T':  op( T ) = T';
C             = 'C':  op( T ) = T'.
C
C     DIAG    CHARACTER*1.
C             Specifies whether or not T is unit triangular, as follows:
C             = 'U':  the matrix T is assumed to be unit triangular;
C             = 'N':  the matrix T is not assumed to be unit triangular.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of H.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of H.  N >= 0.
C
C     L       (input) INTEGER
C             If UPLO = 'U', matrix H has L nonzero subdiagonals.
C             If UPLO = 'L', matrix H has L nonzero superdiagonals.
C             MAX(0,M-1) >= L >= 0, if UPLO = 'U';
C             MAX(0,N-1) >= L >= 0, if UPLO = 'L'.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then T is not
C             referenced and H need not be set before entry.
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT,k), where
C             k is m when SIDE = 'L' and is n when SIDE = 'R'.
C             If UPLO = 'U', the leading k-by-k upper triangular part
C             of this array must contain the upper triangular matrix T
C             and the strictly lower triangular part is not referenced.
C             If UPLO = 'L', the leading k-by-k lower triangular part
C             of this array must contain the lower triangular matrix T
C             and the strictly upper triangular part is not referenced.
C             Note that when DIAG = 'U', the diagonal elements of T are
C             not referenced either, but are assumed to be unity.
C
C     LDT     INTEGER
C             The leading dimension of array T.
C             LDT >= MAX(1,M), if SIDE = 'L';
C             LDT >= MAX(1,N), if SIDE = 'R'.
C
C     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C             On entry, if UPLO = 'U', the leading M-by-N upper
C             Hessenberg part of this array must contain the upper
C             Hessenberg-like matrix H.
C             On entry, if UPLO = 'L', the leading M-by-N lower
C             Hessenberg part of this array must contain the lower
C             Hessenberg-like matrix H.
C             On exit, the leading M-by-N part of this array contains
C             the matrix product alpha*op( T )*H, if SIDE = 'L',
C             or alpha*H*op( T ), if SIDE = 'R'. If TRANST = 'N', this
C             product has the same pattern as the given matrix H;
C             the elements below the L-th subdiagonal (if UPLO = 'U'),
C             or above the L-th superdiagonal (if UPLO = 'L'), are not
C             referenced in this case. If TRANST = 'T', the elements
C             below the (N+L)-th row (if UPLO = 'U', SIDE = 'R', and
C             M > N+L), or at the right of the (M+L)-th column
C             (if UPLO = 'L', SIDE = 'L', and N > M+L), are not set to
C             zero nor referenced.
C
C     LDH     INTEGER
C             The leading dimension of array H.  LDH >= max(1,M).
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
C     The calculations are efficiently performed taking the problem
C     structure into account.
C
C     FURTHER COMMENTS
C
C     The matrix H may have the following patterns, when m = 7, n = 6,
C     and l = 2 are used for illustration:
C
C               UPLO = 'U'                    UPLO = 'L'
C
C            [ x x x x x x ]               [ x x x 0 0 0 ]
C            [ x x x x x x ]               [ x x x x 0 0 ]
C            [ x x x x x x ]               [ x x x x x 0 ]
C        H = [ 0 x x x x x ],          H = [ x x x x x x ].
C            [ 0 0 x x x x ]               [ x x x x x x ]
C            [ 0 0 0 x x x ]               [ x x x x x x ]
C            [ 0 0 0 0 x x ]               [ x x x x x x ]
C
C     The products T*H or H*T have the same pattern as H, but the
C     products T'*H or H*T' may be full matrices.
C
C     If m = n, the matrix H is upper or lower triangular, for l = 0,
C     and upper or lower Hessenberg, for l = 1.
C
C     This routine is a specialization of the BLAS 3 routine DTRMM.
C     BLAS 1 calls are used when appropriate, instead of in-line code,
C     in order to increase the efficiency. If the matrix H is full, or
C     its zero triangle has small order, an optimized DTRMM code could
C     be faster than MB01ZD.
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          DIAG, SIDE, TRANST, UPLO
      INTEGER            INFO, L, LDH, LDT, M, N
      DOUBLE PRECISION   ALPHA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), T( LDT, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, TRANS, UPPER
      INTEGER            I, I1, I2, J, K, M2, NROWT
      DOUBLE PRECISION   TEMP
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      LSIDE  = LSAME( SIDE,   'L' )
      UPPER  = LSAME( UPLO,   'U' )
      TRANS  = LSAME( TRANST, 'T' ) .OR. LSAME( TRANST, 'C' )
      NOUNIT = LSAME( DIAG,   'N' )
      IF( LSIDE )THEN
         NROWT = M
      ELSE
         NROWT = N
      END IF
C
      IF( UPPER )THEN
         M2 = M
      ELSE
         M2 = N
      END IF
C
      INFO   = 0
      IF(      .NOT.( LSIDE  .OR. LSAME( SIDE,   'R' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER  .OR. LSAME( UPLO,   'L' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( TRANS  .OR. LSAME( TRANST, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( NOUNIT .OR. LSAME( DIAG,   'U' ) ) ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( L.LT.0 .OR. L.GT.MAX( 0, M2-1 ) ) THEN
         INFO = -7
      ELSE IF( LDT.LT.MAX( 1, NROWT ) ) THEN
         INFO = -10
      ELSE IF( LDH.LT.MAX( 1, M ) )THEN
         INFO = -12
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01ZD', -INFO )
         RETURN
      END IF
C
C     Quick return, if possible.
C
      IF( MIN( M, N ).EQ.0 )
     $   RETURN
C
C     Also, when alpha = 0.
C
      IF( ALPHA.EQ.ZERO ) THEN
C
         DO 20, J = 1, N
            IF( UPPER ) THEN
               I1 = 1
               I2 = MIN( J+L, M )
            ELSE
               I1 = MAX( 1, J-L )
               I2 = M
            END IF
C
            DO 10, I = I1, I2
               H( I, J ) = ZERO
   10       CONTINUE
C
   20    CONTINUE
C
         RETURN
      END IF
C
C     Start the operations.
C
      IF( LSIDE )THEN
         IF( .NOT.TRANS ) THEN
C
C           Form  H := alpha*T*H.
C
            IF( UPPER ) THEN
C
               DO 40, J = 1, N
C
                  DO 30, K = 1, MIN( J+L, M )
                     IF( H( K, J ).NE.ZERO ) THEN
                        TEMP = ALPHA*H( K, J )
                        CALL DAXPY ( K-1, TEMP, T( 1, K ), 1, H( 1, J ),
     $                               1 )
                        IF( NOUNIT )
     $                     TEMP = TEMP*T( K, K )
                        H( K, J ) = TEMP
                     END IF
   30             CONTINUE
C
   40          CONTINUE
C
            ELSE
C
               DO 60, J = 1, N
C
                  DO 50 K = M, MAX( 1, J-L ), -1
                     IF( H( K, J ).NE.ZERO ) THEN
                        TEMP      = ALPHA*H( K, J )
                        H( K, J ) = TEMP
                        IF( NOUNIT )
     $                     H( K, J ) = H( K, J )*T( K, K )
                        CALL DAXPY ( M-K, TEMP, T( K+1, K ), 1,
     $                                          H( K+1, J ), 1 )
                     END IF
   50             CONTINUE
C
   60          CONTINUE
C
            END IF
C
         ELSE
C
C           Form  H := alpha*T'*H.
C
            IF( UPPER ) THEN
C
               DO 80, J = 1, N
                  I1 = J + L
C
                  DO 70, I = M, 1, -1
                     IF( I.GT.I1 ) THEN
                        TEMP = DDOT( I1, T( 1, I ), 1, H( 1, J ), 1 )
                     ELSE
                        TEMP = H( I, J )
                        IF( NOUNIT )
     $                     TEMP = TEMP*T( I, I )
                        TEMP = TEMP + DDOT( I-1, T( 1, I ), 1,
     $                                           H( 1, J ), 1 )
                     END IF
                     H( I, J ) = ALPHA*TEMP
   70             CONTINUE
C
   80          CONTINUE
C
            ELSE
C
               DO 100, J = 1, MIN( M+L, N )
                  I1 = J - L
C
                  DO 90, I = 1, M
                     IF( I.LT.I1 ) THEN
                        TEMP = DDOT( M-I1+1, T( I1, I ), 1, H( I1, J ),
     $                               1 )
                     ELSE
                        TEMP = H( I, J )
                        IF( NOUNIT )
     $                     TEMP = TEMP*T( I, I )
                        TEMP = TEMP + DDOT( M-I, T( I+1, I ), 1,
     $                                           H( I+1, J ), 1 )
                     END IF
                     H( I, J ) = ALPHA*TEMP
   90             CONTINUE
C
  100          CONTINUE
C
            END IF
C
         END IF
C
      ELSE
C
         IF( .NOT.TRANS ) THEN
C
C           Form  H := alpha*H*T.
C
            IF( UPPER ) THEN
C
               DO 120, J = N, 1, -1
                  I2   = MIN( J+L, M )
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*T( J, J )
                  CALL DSCAL ( I2, TEMP, H( 1, J ), 1 )
C
                  DO 110, K = 1, J - 1
                     CALL DAXPY ( I2, ALPHA*T( K, J ), H( 1, K ), 1,
     $                                                 H( 1, J ), 1 )
  110             CONTINUE
C
  120          CONTINUE
C
            ELSE
C
               DO 140, J = 1, N
                  I1   = MAX( 1, J-L )
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*T( J, J )
                  CALL DSCAL ( M-I1+1, TEMP, H( I1, J ), 1 )
C
                  DO 130, K = J + 1, N
                     CALL DAXPY ( M-I1+1, ALPHA*T( K, J ), H( I1, K ),
     $                            1, H( I1, J ), 1 )
  130             CONTINUE
C
  140          CONTINUE
C
            END IF
C
         ELSE
C
C           Form  H := alpha*H*T'.
C
            IF( UPPER ) THEN
               M2 = MIN( N+L, M )
C
               DO 170, K = 1, N
                  I1 = MIN( K+L, M )
                  I2 = MIN( K+L, M2 )
C
                  DO 160, J = 1, K - 1
                     IF( T( J, K ).NE.ZERO ) THEN
                        TEMP = ALPHA*T( J, K )
                        CALL DAXPY ( I1, TEMP, H( 1, K ), 1, H( 1, J ),
     $                               1 )
C
                        DO 150, I = I1 + 1, I2
                           H( I, J ) = TEMP*H( I, K )
  150                   CONTINUE
C
                     END IF
  160             CONTINUE
C
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*T( K, K )
                  IF( TEMP.NE.ONE )
     $               CALL DSCAL( I2, TEMP, H( 1, K ), 1 )
  170          CONTINUE
C
            ELSE
C
               DO 200, K = N, 1, -1
                  I1 = MAX( 1, K-L )
                  I2 = MAX( 1, K-L+1 )
                  M2 = MIN( M, I2-1 )
C
                  DO 190, J = K + 1, N
                     IF( T( J, K ).NE.ZERO ) THEN
                        TEMP = ALPHA*T( J, K )
                        CALL DAXPY ( M-I2+1, TEMP, H( I2, K ), 1,
     $                               H( I2, J ), 1 )
C
                        DO 180, I = I1, M2
                           H( I, J ) = TEMP*H( I, K )
  180                   CONTINUE
C
                     END IF
  190             CONTINUE
C
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*T( K, K )
                  IF( TEMP.NE.ONE )
     $               CALL DSCAL( M-I1+1, TEMP, H( I1, K ), 1 )
  200          CONTINUE
C
            END IF
C
         END IF
C
      END IF
C
      RETURN
C
C *** Last line of MB01ZD ***
      END
