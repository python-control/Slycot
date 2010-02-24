      SUBROUTINE MB01WD( DICO, UPLO, TRANS, HESS, N, ALPHA, BETA, R,
     $                   LDR, A, LDA, T, LDT, INFO )
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
C     To compute the matrix formula
C     _
C     R = alpha*( op( A )'*op( T )'*op( T ) + op( T )'*op( T )*op( A ) )
C         + beta*R,                                                  (1)
C
C     if DICO = 'C', or
C     _
C     R = alpha*( op( A )'*op( T )'*op( T )*op( A ) -  op( T )'*op( T ))
C         + beta*R,                                                  (2)
C                                                             _
C     if DICO = 'D', where alpha and beta are scalars, R, and R are
C     symmetric matrices, T is a triangular matrix, A is a general or
C     Hessenberg matrix, and op( M ) is one of
C
C        op( M ) = M   or   op( M ) = M'.
C
C     The result is overwritten on R.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the formula to be evaluated, as follows:
C             = 'C':  formula (1), "continuous-time" case;
C             = 'D':  formula (2), "discrete-time" case.
C
C     UPLO    CHARACTER*1
C             Specifies which triangles of the symmetric matrix R and
C             triangular matrix T are given, as follows:
C             = 'U':  the upper triangular parts of R and T are given;
C             = 'L':  the lower triangular parts of R and T are given;
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( M ) to be used, as follows:
C             = 'N':  op( M ) = M;
C             = 'T':  op( M ) = M';
C             = 'C':  op( M ) = M'.
C
C     HESS    CHARACTER*1
C             Specifies the form of the matrix A, as follows:
C             = 'F':  matrix A is full;
C             = 'H':  matrix A is Hessenberg (or Schur), either upper
C                     (if UPLO = 'U'), or lower (if UPLO = 'L').
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices R, A, and T.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then the arrays A
C             and T are not referenced.
C
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. When beta is zero then the array R need
C             not be set before entry.
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N)
C             On entry with UPLO = 'U', the leading N-by-N upper
C             triangular part of this array must contain the upper
C             triangular part of the symmetric matrix R.
C             On entry with UPLO = 'L', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangular part of the symmetric matrix R.
C             On exit, the leading N-by-N upper triangular part (if
C             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of
C             this array contains the corresponding triangular part of
C                                 _
C             the computed matrix R.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A. If HESS = 'H' the elements below the
C             first subdiagonal, if UPLO = 'U', or above the first
C             superdiagonal, if UPLO = 'L', need not be set to zero,
C             and are not referenced if DICO = 'D'.
C             On exit, the leading N-by-N part of this array contains
C             the following matrix product
C                alpha*T'*T*A, if TRANS = 'N', or
C                alpha*A*T*T', otherwise,
C             if DICO = 'C', or
C                T*A, if TRANS = 'N', or
C                A*T, otherwise,
C             if DICO = 'D' (and in this case, these products have a
C             Hessenberg form, if HESS = 'H').
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT,N)
C             If UPLO = 'U', the leading N-by-N upper triangular part of
C             this array must contain the upper triangular matrix T and
C             the strictly lower triangular part need not be set to zero
C             (and it is not referenced).
C             If UPLO = 'L', the leading N-by-N lower triangular part of
C             this array must contain the lower triangular matrix T and
C             the strictly upper triangular part need not be set to zero
C             (and it is not referenced).
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= MAX(1,N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -k, the k-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The matrix expression (1) or (2) is efficiently evaluated taking
C     the structure into account. BLAS 3 operations (DTRMM, DSYRK and
C     their specializations) are used throughout.
C
C     NUMERICAL ASPECTS
C
C     If A is a full matrix, the algorithm requires approximately
C      3
C     N  operations, if DICO = 'C';
C            3
C     7/6 x N  operations, if DICO = 'D'.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, HESS, TRANS, UPLO
      INTEGER           INFO, LDA, LDR, LDT, N
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), R(LDR,*), T(LDT,*)
C     .. Local Scalars ..
      LOGICAL           DISCR, REDUC, TRANSP, UPPER
      CHARACTER         NEGTRA, SIDE
      INTEGER           I, INFO2, J
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLASCL, DLASET, DSYRK, DTRMM, MB01YD, MB01ZD,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      DISCR  = LSAME( DICO,  'D' )
      UPPER  = LSAME( UPLO,  'U' )
      TRANSP = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      REDUC  = LSAME( HESS,  'H' )
C
      IF(      .NOT.( DISCR  .OR. LSAME( DICO,  'C' ) ) )THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER  .OR. LSAME( UPLO,  'L' ) ) )THEN
         INFO = -2
      ELSE IF( .NOT.( TRANSP .OR. LSAME( TRANS, 'N' ) ) )THEN
         INFO = -3
      ELSE IF( .NOT.( REDUC  .OR. LSAME( HESS,  'F' ) ) )THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -13
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01WD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      IF ( ALPHA.EQ.ZERO ) THEN
         IF ( BETA.EQ.ZERO ) THEN
C
C           Special case when both alpha = 0 and beta = 0.
C
            CALL DLASET( UPLO, N, N, ZERO, ZERO, R, LDR )
         ELSE
C
C           Special case alpha = 0.
C
            IF ( BETA.NE.ONE )
     $         CALL DLASCL( UPLO, 0, 0, ONE, BETA, N, N, R, LDR, INFO2 )
         END IF
         RETURN
      END IF
C
C     General case: alpha <> 0.
C
C     Compute (in A) T*A, if TRANS = 'N', or
C                    A*T, otherwise.
C
      IF ( TRANSP ) THEN
         SIDE   = 'R'
         NEGTRA = 'N'
      ELSE
         SIDE   = 'L'
         NEGTRA = 'T'
      END IF
C
      IF ( REDUC .AND. N.GT.2 ) THEN
         CALL MB01ZD( SIDE, UPLO, 'NoTranspose', 'Non-unit', N, N, 1,
     $                ONE, T, LDT, A, LDA, INFO2 )
      ELSE
         CALL DTRMM( SIDE, UPLO, 'NoTranspose', 'Non-unit', N, N, ONE,
     $               T, LDT, A, LDA )
      END IF
C
      IF( .NOT.DISCR ) THEN
C
C        Compute (in A) alpha*T'*T*A, if TRANS = 'N', or
C                       alpha*A*T*T', otherwise.
C
         IF ( REDUC .AND. N.GT.2 ) THEN
            CALL MB01ZD( SIDE, UPLO, 'Transpose', 'Non-unit', N, N, 1,
     $                   ALPHA, T, LDT, A, LDA, INFO2 )
         ELSE
            CALL DTRMM( SIDE, UPLO, 'Transpose', 'Non-unit', N, N,
     $                  ALPHA, T, LDT, A, LDA )
         END IF
C
C        Compute the required triangle of the result, using symmetry.
C
         IF ( UPPER ) THEN
            IF ( BETA.EQ.ZERO ) THEN
C
               DO 20 J = 1, N
                  DO 10 I = 1, J
                     R( I, J ) = A( I, J ) + A( J, I )
   10             CONTINUE
   20          CONTINUE
C
            ELSE
C
               DO 40 J = 1, N
                  DO 30 I = 1, J
                     R( I, J ) = A( I, J ) + A( J, I ) + BETA*R( I, J )
   30             CONTINUE
   40          CONTINUE
C
            END IF
C
         ELSE
C
            IF ( BETA.EQ.ZERO ) THEN
C
               DO 60 J = 1, N
                  DO 50 I = J, N
                     R( I, J ) = A( I, J ) + A( J, I )
   50             CONTINUE
   60          CONTINUE
C
            ELSE
C
               DO 80 J = 1, N
                  DO 70 I = J, N
                     R( I, J ) = A( I, J ) + A( J, I ) + BETA*R( I, J )
   70             CONTINUE
   80          CONTINUE
C
            END IF
C
         END IF
C
      ELSE
C
C        Compute (in R) alpha*A'*T'*T*A + beta*R, if TRANS = 'N', or
C                       alpha*A*T*T'*A' + beta*R, otherwise.
C
         IF ( REDUC .AND. N.GT.2 ) THEN
            CALL MB01YD( UPLO, NEGTRA, N, N, 1, ALPHA, BETA, A, LDA, R,
     $                   LDR, INFO2 )
         ELSE
            CALL DSYRK( UPLO, NEGTRA, N, N, ALPHA, A, LDA, BETA, R,
     $                  LDR )
         END IF
C
C        Compute (in R) -alpha*T'*T + R, if TRANS = 'N', or
C                       -alpha*T*T' + R, otherwise.
C
         CALL MB01YD( UPLO, NEGTRA, N, N, 0, -ALPHA, ONE, T, LDT, R,
     $                LDR, INFO2 )
C
      END IF
C
      RETURN
C *** Last line of MB01WD ***
      END
