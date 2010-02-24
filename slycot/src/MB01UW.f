      SUBROUTINE MB01UW( SIDE, TRANS, M, N, ALPHA, H, LDH, A, LDA,
     $                   DWORK, LDWORK, INFO )
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
C     To compute one of the matrix products
C
C        A : = alpha*op( H ) * A, or A : = alpha*A * op( H ),
C
C     where alpha is a scalar, A is an m-by-n matrix, H is an upper
C     Hessenberg matrix, and op( H ) is one of
C
C        op( H ) = H   or   op( H ) = H',  the transpose of H.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies whether the Hessenberg matrix H appears on the
C             left or right in the matrix product as follows:
C             = 'L':  A := alpha*op( H ) * A;
C             = 'R':  A := alpha*A * op( H ).
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( H ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( H ) = H;
C             = 'T':  op( H ) = H';
C             = 'C':  op( H ) = H'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then H is not
C             referenced and A need not be set before entry.
C
C     H       (input) DOUBLE PRECISION array, dimension (LDH,k)
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C             On entry with SIDE = 'L', the leading M-by-M upper
C             Hessenberg part of this array must contain the upper
C             Hessenberg matrix H.
C             On entry with SIDE = 'R', the leading N-by-N upper
C             Hessenberg part of this array must contain the upper
C             Hessenberg matrix H.
C             The elements below the subdiagonal are not referenced,
C             except possibly for those in the first column, which
C             could be overwritten, but are restored on exit.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= max(1,k),
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading M-by-N part of this array contains
C             the computed product.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, alpha <> 0, and LDWORK >= M*N > 0,
C             DWORK contains a copy of the matrix A, having the leading
C             dimension M.
C             This array is not referenced when alpha = 0.
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= 0,   if  alpha =  0 or MIN(M,N) = 0;
C             LDWORK >= M-1, if  SIDE  = 'L';
C             LDWORK >= N-1, if  SIDE  = 'R'.
C             For maximal efficiency LDWORK should be at least M*N.
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
C     The required matrix product is computed in two steps. In the first
C     step, the upper triangle of H is used; in the second step, the
C     contribution of the subdiagonal is added. If the workspace can
C     accomodate a copy of A, a fast BLAS 3 DTRMM operation is used in
C     the first step.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, January 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
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
C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS
      INTEGER           INFO, LDA, LDH, LDWORK, M, N
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), H(LDH,*)
C     .. Local Scalars ..
      LOGICAL           LSIDE, LTRANS
      INTEGER           I, J, JW
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLACPY, DLASCL, DLASET, DSCAL, DSWAP,
     $                  DTRMM, DTRMV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO   = 0
      LSIDE  = LSAME( SIDE,  'L' )
      LTRANS = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      IF(      ( .NOT.LSIDE  ).AND.( .NOT.LSAME( SIDE,  'R' ) ) )THEN
         INFO = -1
      ELSE IF( ( .NOT.LTRANS ).AND.( .NOT.LSAME( TRANS, 'N' ) ) )THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDH.LT.1 .OR. ( LSIDE .AND. LDH.LT.M ) .OR.
     $                  ( .NOT.LSIDE .AND. LDH.LT.N ) ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( LDWORK.LT.0 .OR.
     $       ( ALPHA.NE.ZERO .AND. MIN( M, N ).GT.0 .AND.
     $            ( ( LSIDE .AND. LDWORK.LT.M-1 ) .OR.
     $         ( .NOT.LSIDE .AND. LDWORK.LT.N-1 ) ) ) ) THEN
         INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01UW', -INFO )
         RETURN
      END IF
C
C     Quick return, if possible.
C
      IF ( MIN( M, N ).EQ.0 ) THEN
         RETURN
      ELSE IF ( LSIDE ) THEN
         IF ( M.EQ.1 ) THEN
            CALL DSCAL( N, ALPHA*H(1,1), A, LDA )
            RETURN
         END IF
      ELSE
         IF ( N.EQ.1 ) THEN
            CALL DSCAL( M, ALPHA*H(1,1), A, 1 )
            RETURN
         END IF
      END IF
C
      IF( ALPHA.EQ.ZERO ) THEN
C
C        Set A to zero and return.
C
         CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
         RETURN
      END IF
C
      IF( LDWORK.GE.M*N ) THEN
C
C        Enough workspace for a fast BLAS 3 calculation.
C        Save A in the workspace and compute one of the matrix products
C          A : = alpha*op( triu( H ) ) * A, or
C          A : = alpha*A * op( triu( H ) ),
C        involving the upper triangle of H.
C
         CALL DLACPY( 'Full', M, N, A, LDA, DWORK, M )
         CALL DTRMM( SIDE, 'Upper', TRANS, 'Non-unit', M, N, ALPHA, H,
     $               LDH, A, LDA )
C
C        Add the contribution of the subdiagonal of H.
C        If SIDE = 'L', the subdiagonal of H is swapped with the
C        corresponding elements in the first column of H, and the
C        calculations are organized for column operations.
C
         IF( LSIDE ) THEN
            IF( M.GT.2 )
     $         CALL DSWAP( M-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
            IF( LTRANS ) THEN
               JW = 1
               DO 20 J = 1, N
                  JW = JW + 1
                  DO 10 I = 1, M - 1
                     A( I, J ) = A( I, J ) +
     $                           ALPHA*H( I+1, 1 )*DWORK( JW )
                     JW = JW + 1
   10             CONTINUE
   20          CONTINUE
            ELSE
               JW = 0
               DO 40 J = 1, N
                  JW = JW + 1
                  DO 30 I = 2, M
                     A( I, J ) = A( I, J ) +
     $                           ALPHA*H( I, 1 )*DWORK( JW )
                     JW = JW + 1
   30             CONTINUE
   40          CONTINUE
            END IF
            IF( M.GT.2 )
     $         CALL DSWAP( M-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
C
         ELSE
C
            IF( LTRANS ) THEN
               JW = 1
               DO 50 J = 1, N - 1
                  IF ( H( J+1, J ).NE.ZERO )
     $               CALL DAXPY( M, ALPHA*H( J+1, J ), DWORK( JW ), 1,
     $                           A( 1, J+1 ), 1 )
                  JW = JW + M
   50          CONTINUE
            ELSE
               JW = M + 1
               DO 60 J = 1, N - 1
                  IF ( H( J+1, J ).NE.ZERO )
     $               CALL DAXPY( M, ALPHA*H( J+1, J ), DWORK( JW ), 1,
     $                           A( 1, J ), 1 )
                  JW = JW + M
   60          CONTINUE
            END IF
         END IF
C
      ELSE
C
C        Use a BLAS 2 calculation.
C
         IF( LSIDE ) THEN
            IF( M.GT.2 )
     $         CALL DSWAP( M-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
            IF( LTRANS ) THEN
               DO 80 J = 1, N
C
C                 Compute the contribution of the subdiagonal of H to
C                 the j-th column of the product.
C
                  DO 70 I = 1, M - 1
                     DWORK( I ) = H( I+1, 1 )*A( I+1, J )
   70             CONTINUE
C
C                 Multiply the upper triangle of H by the j-th column
C                 of A, and add to the above result.
C
                  CALL DTRMV( 'Upper', TRANS, 'Non-unit', M, H, LDH,
     $                        A( 1, J ), 1 )
                  CALL DAXPY( M-1, ONE, DWORK, 1, A( 1, J ), 1 )
   80          CONTINUE
C
            ELSE
               DO 100 J = 1, N
C
C                 Compute the contribution of the subdiagonal of H to
C                 the j-th column of the product.
C
                  DO 90 I = 1, M - 1
                     DWORK( I ) = H( I+1, 1 )*A( I, J )
   90             CONTINUE
C
C                 Multiply the upper triangle of H by the j-th column
C                 of A, and add to the above result.
C
                  CALL DTRMV( 'Upper', TRANS, 'Non-unit', M, H, LDH,
     $                        A( 1, J  ), 1 )
                  CALL DAXPY( M-1, ONE, DWORK, 1, A( 2, J ), 1 )
  100          CONTINUE
            END IF
            IF( M.GT.2 )
     $         CALL DSWAP( M-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
C
         ELSE
C
C           Below, row-wise calculations are used for A.
C
            IF( N.GT.2 )
     $         CALL DSWAP( N-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
            IF( LTRANS ) THEN
               DO 120 I = 1, M
C
C                 Compute the contribution of the subdiagonal of H to
C                 the i-th row of the product.
C
                  DO 110 J = 1, N - 1
                     DWORK( J ) = A( I, J )*H( J+1, 1 )
  110             CONTINUE
C
C                 Multiply the i-th row of A by the upper triangle of H,
C                 and add to the above result.
C
                  CALL DTRMV( 'Upper', 'NoTranspose', 'Non-unit', N, H,
     $                        LDH, A( I, 1 ), LDA )
                  CALL DAXPY( N-1, ONE, DWORK, 1, A( I, 2 ), LDA )
  120          CONTINUE
C
            ELSE
               DO 140 I = 1, M
C
C                 Compute the contribution of the subdiagonal of H to
C                 the i-th row of the product.
C
                  DO 130 J = 1, N - 1
                     DWORK( J ) = A( I, J+1 )*H( J+1, 1 )
  130             CONTINUE
C
C                 Multiply the i-th row of A by the upper triangle of H,
C                 and add to the above result.
C
                  CALL DTRMV( 'Upper', 'Transpose', 'Non-unit', N, H,
     $                        LDH, A( I, 1 ), LDA )
                  CALL DAXPY( N-1, ONE, DWORK, 1, A( I, 1 ), LDA )
  140          CONTINUE
            END IF
            IF( N.GT.2 )
     $         CALL DSWAP( N-2, H( 3, 2 ), LDH+1, H( 3, 1 ), 1 )
C
         END IF
C
C        Scale the result by alpha.
C
         IF ( ALPHA.NE.ONE )
     $      CALL DLASCL( 'General', 0, 0, ONE, ALPHA, M, N, A, LDA,
     $                   INFO )
      END IF
      RETURN
C *** Last line of MB01UW ***
      END
