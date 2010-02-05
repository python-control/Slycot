      SUBROUTINE MB02OD( SIDE, UPLO, TRANS, DIAG, NORM, M, N, ALPHA, A,
     $                   LDA, B, LDB, RCOND, TOL, IWORK, DWORK, INFO )
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
C     To solve (if well-conditioned) one of the matrix equations
C
C        op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C
C     where alpha is a scalar, X and B are m-by-n matrices, A is a unit,
C     or non-unit, upper or lower triangular matrix and op( A ) is one
C     of
C
C        op( A ) = A   or   op( A ) = A'.
C
C     An estimate of the reciprocal of the condition number of the
C     triangular matrix A, in either the 1-norm or the infinity-norm, is
C     also computed as
C
C        RCOND = 1 / ( norm(A) * norm(inv(A)) ).
C
C     and the specified matrix equation is solved only if RCOND is
C     larger than a given tolerance TOL.  In that case, the matrix X is
C     overwritten on B.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies whether op( A ) appears on the left or right
C             of X as follows:
C             = 'L':  op( A )*X = alpha*B;
C             = 'R':  X*op( A ) = alpha*B.
C
C     UPLO    CHARACTER*1
C             Specifies whether the matrix A is an upper or lower
C             triangular matrix as follows:
C             = 'U':  A is an upper triangular matrix;
C             = 'L':  A is a lower triangular matrix.
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( A ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     DIAG    CHARACTER*1
C             Specifies whether or not A is unit triangular as follows:
C             = 'U':  A is assumed to be unit triangular;
C             = 'N':  A is not assumed to be unit triangular.
C
C     NORM    CHARACTER*1
C             Specifies whether the 1-norm condition number or the
C             infinity-norm condition number is required:
C             = '1' or 'O':  1-norm;
C             = 'I':         Infinity-norm.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of B.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of B.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar  alpha. When alpha is zero then A is not
C             referenced and B need not be set before entry.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,k),
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C             On entry with UPLO = 'U', the leading k-by-k upper
C             triangular part of this array must contain the upper
C             triangular matrix and the strictly lower triangular part
C             of A is not referenced.
C             On entry with UPLO = 'L', the leading k-by-k lower
C             triangular part of this array must contain the lower
C             triangular matrix and the strictly upper triangular part
C             of A is not referenced.
C             Note that when DIAG = 'U', the diagonal elements of A are
C             not referenced either, but are assumed to be unity.
C
C     LDA     INTEGER
C             The leading dimension of array A.
C             LDA >= max(1,M) when SIDE = 'L';
C             LDA >= max(1,N) when SIDE = 'R'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand side matrix B.
C             On exit, if INFO = 0, the leading M-by-N part of this
C             array contains the solution matrix X.
C             Otherwise, this array is not modified by the routine.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= max(1,M).
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal of the condition number of the matrix A,
C             computed as RCOND = 1/(norm(A) * norm(inv(A))).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the matrix A. If the user sets TOL > 0, then the given
C             value of TOL is used as a lower bound for the reciprocal
C             condition number of that matrix; a matrix whose estimated
C             condition number is less than 1/TOL is considered to be
C             nonsingular. If the user sets TOL <= 0, then an implicitly
C             computed, default tolerance, defined by TOLDEF = k*k*EPS,
C             is used instead, where EPS is the machine precision (see
C             LAPACK Library routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (k)
C
C     DWORK   DOUBLE PRECISION array, dimension (3*k)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the matrix A is numerically singular, i.e. the
C                   condition number estimate of A (in the specified
C                   norm) exceeds 1/TOL.
C
C     METHOD
C
C     An estimate of the reciprocal of the condition number of the
C     triangular matrix A (in the specified norm) is computed, and if
C     this estimate is larger then the given (or default) tolerance,
C     the specified matrix equation is solved using Level 3 BLAS
C     routine DTRSM.
C
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C                             2
C     The algorithm requires k N/2 operations.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     February 20, 1998.
C
C     KEYWORDS
C
C     Condition number, matrix algebra, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, M, N
      DOUBLE PRECISION   ALPHA, RCOND, TOL
C     .. Array Arguments ..
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL            LSIDE, ONENRM
      INTEGER            NROWA
      DOUBLE PRECISION   TOLDEF
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL           DTRCON, DTRSM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
C     .. Executable Statements ..
C
      LSIDE  = LSAME( SIDE, 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.LSAME( SIDE, 'R' ) ) )THEN
         INFO = -1
      ELSE IF( ( .NOT.LSAME( UPLO,  'U' ) ).AND.
     $         ( .NOT.LSAME( UPLO,  'L' ) )      )THEN
         INFO = -2
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = -3
      ELSE IF( ( .NOT.LSAME( DIAG,  'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG,  'N' ) )      )THEN
         INFO = -4
      ELSE IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 )THEN
         INFO = -6
      ELSE IF( N.LT.0 )THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, M ) )THEN
         INFO = -12
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( NROWA.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
C
      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO )
     $   TOLDEF = DBLE( NROWA*NROWA )*DLAMCH( 'Epsilon' )
C
      CALL DTRCON( NORM, UPLO, DIAG, NROWA, A, LDA, RCOND, DWORK,
     $             IWORK, INFO )
C
      IF ( RCOND.GT.TOLDEF ) THEN
         CALL DTRSM( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B,
     $               LDB )
      ELSE
         INFO = 1
      END IF
C *** Last line of MB02OD ***
      END
