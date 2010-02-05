      SUBROUTINE MB01MD( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                   INCY )
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
C     To perform the matrix-vector operation
C
C        y := alpha*A*x + beta*y,
C
C     where alpha and beta are scalars, x and y are vectors of length
C     n and A is an n-by-n skew-symmetric matrix.
C
C     This is a modified version of the vanilla implemented BLAS
C     routine DSYMV written by Jack Dongarra, Jeremy Du Croz,
C     Sven Hammarling, and Richard Hanson.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies whether the upper or lower triangular part of
C             the array A is to be referenced as follows:
C             = 'U':  only the strictly upper triangular part of A is to
C                     be referenced;
C             = 'L':  only the strictly lower triangular part of A is to
C                     be referenced.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. If alpha is zero the array A is not
C             referenced.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry with UPLO = 'U', the leading N-by-N part of this
C             array must contain the strictly upper triangular part of
C             the matrix A. The lower triangular part of this array is
C             not referenced.
C             On entry with UPLO = 'L', the leading N-by-N part of this
C             array must contain the strictly lower triangular part of
C             the matrix A. The upper triangular part of this array is
C             not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N)
C
C     X       (input) DOUBLE PRECISION array, dimension
C             ( 1 + ( N - 1 )*abs( INCX ) ).
C             On entry, elements 1, INCX+1, .., ( N - 1 )*INCX + 1 of
C             this array must contain the elements of the vector X.
C
C     INCX    (input) INTEGER
C             The increment for the elements of X. IF INCX < 0 then the
C             elements of X are accessed in reversed order.  INCX <> 0.
C
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. If beta is zero then Y need not be set on
C             input.
C
C     Y       (input/output) DOUBLE PRECISION array, dimension
C             ( 1 + ( N - 1 )*abs( INCY ) ).
C             On entry, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of
C             this array must contain the elements of the vector Y.
C             On exit, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of
C             this array contain the updated elements of the vector Y.
C
C     INCY    (input) INTEGER
C             The increment for the elements of Y. IF INCY < 0 then the
C             elements of Y are accessed in reversed order.  INCY <> 0.
C
C     NUMERICAL ASPECTS
C
C     Though being almost identical with the vanilla implementation
C     of the BLAS routine DSYMV the performance of this routine could
C     be significantly lower in the case of vendor supplied, highly
C     optimized BLAS.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DSKMV).
C
C     KEYWORDS
C
C     Elementary matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER          UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF ( N.LT.0 )THEN
         INFO = 2
      ELSE IF ( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF ( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF ( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF ( INFO.NE.0 )THEN
         CALL XERBLA( 'MB01MD', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF ( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF ( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      IF ( BETA.NE.ONE )THEN
         IF ( INCY.EQ.1 )THEN
            IF ( BETA.EQ.ZERO )THEN
               DO 10 I = 1, N
                  Y(I) = ZERO
   10          CONTINUE
            ELSE
               DO 20 I = 1, N
                  Y(I) = BETA*Y(I)
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF ( BETA.EQ.ZERO )THEN
               DO 30 I = 1, N
                  Y(IY) = ZERO
                  IY = IY + INCY
   30          CONTINUE
            ELSE
               DO 40 I = 1, N
                  Y(IY) = BETA*Y(IY)
                  IY = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
C
C     Quick return if possible.
C
      IF ( ALPHA.EQ.ZERO )
     $   RETURN
      IF ( LSAME( UPLO, 'U' ) )THEN
C
C        Form y when A is stored in upper triangle.
C
         IF ( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60 J = 2, N
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y(I)  = Y(I)  + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(I)
   50          CONTINUE
               Y(J) = Y(J) - ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX + INCX
            JY = KY + INCY
            DO 80 J = 2, N
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70 I = 1, J - 1
                  Y(IY) = Y(IY) + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(IX)
                  IX    = IX + INCX
                  IY    = IY + INCY
   70          CONTINUE
               Y(JY) = Y(JY) - ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form y when A is stored in lower triangle.
C
         IF ( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) )THEN
            DO 100 J = 1, N - 1
               TEMP1 = ALPHA*X(J)
               TEMP2 = ZERO
               DO 90 I = J + 1, N
                  Y(I)  = Y(I)  + TEMP1*A(I,J)
                  TEMP2 = TEMP2 + A(I,J)*X(I)
   90          CONTINUE
               Y(J) = Y(J) - ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120 J = 1, N - 1
               TEMP1 = ALPHA*X(JX)
               TEMP2 = ZERO
               IX = JX
               IY = JY
               DO 110 I = J + 1, N
                  IX = IX + INCX
                  IY = IY + INCY
                  Y(IY ) = Y(IY) + TEMP1*A(I,J)
                  TEMP2  = TEMP2 + A(I,J)*X(IX)
  110          CONTINUE
               Y(JY) = Y(JY) - ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  120       CONTINUE
         END IF
      END IF
C *** Last line of MB01MD ***
      END
