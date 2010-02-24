      SUBROUTINE MB01RW( UPLO, TRANS, M, N, A, LDA, Z, LDZ, DWORK,
     $                   INFO )
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
C     To compute the transformation of the symmetric matrix A by the
C     matrix Z in the form
C
C        A := op(Z)*A*op(Z)',
C
C     where op(Z) is either Z or its transpose, Z'.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Specifies whether the upper or lower triangle of A
C             is stored:
C             = 'U':  Upper triangle of A is stored;
C             = 'L':  Lower triangle of A is stored.
C
C     TRANS   CHARACTER*1
C             Specifies whether op(Z) is Z or its transpose Z':
C             = 'N':  op(Z) = Z;
C             = 'T':  op(Z) = Z'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of the resulting symmetric matrix op(Z)*A*op(Z)'
C             and the number of rows of the matrix Z, if TRANS = 'N',
C             or the number of columns of the matrix Z, if TRANS = 'T'.
C             M >= 0.
C
C     N       (input) INTEGER
C             The order of the symmetric matrix A and the number of
C             columns of the matrix Z, if TRANS = 'N', or the number of
C             rows of the matrix Z, if TRANS = 'T'.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C             (LDA,MAX(M,N))
C             On entry, the leading N-by-N upper or lower triangular
C             part of this array must contain the upper (UPLO = 'U')
C             or lower (UPLO = 'L') triangular part of the symmetric
C             matrix A.
C             On exit, the leading M-by-M upper or lower triangular
C             part of this array contains the upper (UPLO = 'U') or
C             lower (UPLO = 'L') triangular part of the symmetric
C             matrix op(Z)*A*op(Z)'.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,M,N).
C
C     Z       (input) DOUBLE PRECISION array, dimension (LDQ,K)
C             where K = N if TRANS = 'N' and K = M if TRANS = 'T'.
C             The leading M-by-N part, if TRANS = 'N', or N-by-M part,
C             if TRANS = 'T', of this array contains the matrix Z.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.
C             LDZ >= MAX(1,M) if TRANS = 'N' and
C             LDZ >= MAX(1,N) if TRANS = 'T'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     FURTHER COMMENTS
C
C     This is a simpler, BLAS 2 version for MB01RD.
C
C     CONTRIBUTOR
C
C     A. Varga, DLR, Feb. 1995.
C
C     REVISIONS
C
C     April 1998 (T. Penzl).
C     Sep. 1998 (V. Sima).
C
C    ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, LDA, LDZ, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           NOTTRA, UPPER
      INTEGER           I, J
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     .. Executable Statements
C
      NOTTRA = LSAME( TRANS, 'N' )
      UPPER  = LSAME( UPLO,  'U' )
C
      INFO = 0
      IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L') ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( NOTTRA .OR. LSAME( TRANS, 'T') ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M, N ) ) THEN
         INFO = -6
      ELSE IF( (      NOTTRA .AND. LDZ.LT.MAX( 1, M ) )   .OR.
     $         ( .NOT.NOTTRA .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB01RW', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
C
      IF ( NOTTRA ) THEN
C
C        Compute Z*A*Z'.
C
         IF ( UPPER ) THEN
C
C           Compute Z*A in A (M-by-N).
C
            DO 10 J = 1, N
               CALL DCOPY( J-1, A(1,J), 1, DWORK, 1 )
               CALL DCOPY( N-J+1, A(J,J), LDA, DWORK(J), 1 )
               CALL DGEMV( TRANS, M, N, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(1,J), 1 )
   10       CONTINUE
C
C           Compute A*Z' in the upper triangular part of A.
C
            DO 20 I = 1, M
               CALL DCOPY( N, A(I,1), LDA, DWORK, 1 )
               CALL DGEMV( TRANS, M-I+1, N, ONE, Z(I,1), LDZ, DWORK, 1,
     $                     ZERO, A(I,I), LDA )
   20       CONTINUE
C
         ELSE
C
C           Compute A*Z' in A (N-by-M).
C
            DO 30 I = 1, N
               CALL DCOPY( I-1, A(I,1), LDA, DWORK, 1 )
               CALL DCOPY( N-I+1, A(I,I), 1, DWORK(I), 1 )
               CALL DGEMV( TRANS, M, N, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(I,1), LDA )
   30       CONTINUE
C
C           Compute Z*A in the lower triangular part of A.
C
            DO 40 J = 1, M
               CALL DCOPY( N, A(1,J), 1, DWORK, 1 )
               CALL DGEMV( TRANS, M-J+1, N, ONE, Z(J,1), LDZ, DWORK, 1,
     $                     ZERO, A(J,J), 1 )
   40       CONTINUE
C
         END IF
      ELSE
C
C        Compute Z'*A*Z.
C
         IF ( UPPER ) THEN
C
C           Compute Z'*A in A (M-by-N).
C
            DO 50 J = 1, N
               CALL DCOPY( J-1, A(1,J), 1, DWORK, 1 )
               CALL DCOPY( N-J+1, A(J,J), LDA, DWORK(J), 1 )
               CALL DGEMV( TRANS, N, M, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(1,J), 1 )
   50       CONTINUE
C
C           Compute A*Z in the upper triangular part of A.
C
            DO 60 I = 1, M
               CALL DCOPY( N, A(I,1), LDA, DWORK, 1 )
               CALL DGEMV( TRANS, N, M-I+1, ONE, Z(1,I), LDZ, DWORK, 1,
     $                     ZERO, A(I,I), LDA )
   60       CONTINUE
C
         ELSE
C
C           Compute A*Z in A (N-by-M).
C
            DO 70 I = 1, N
               CALL DCOPY( I-1, A(I,1), LDA, DWORK, 1 )
               CALL DCOPY( N-I+1, A(I,I), 1, DWORK(I), 1 )
               CALL DGEMV( TRANS, N, M, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                     A(I,1), LDA )
   70       CONTINUE
C
C           Compute Z'*A in the lower triangular part of A.
C
            DO 80 J = 1, M
               CALL DCOPY( N, A(1,J), 1, DWORK, 1 )
               CALL DGEMV( TRANS, N, M-J+1, ONE, Z(1,J), LDZ, DWORK, 1,
     $                     ZERO, A(J,J), 1 )
   80       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of MB01RW ***
      END
