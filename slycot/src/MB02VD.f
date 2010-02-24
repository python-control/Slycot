      SUBROUTINE MB02VD( TRANS, M, N, A, LDA, IPIV, B, LDB, INFO )
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
C     To compute the solution to a real system of linear equations
C        X * op(A) = B,
C     where op(A) is either A or its transpose, A is an N-by-N matrix,
C     and X and B are M-by-N matrices.
C     The LU decomposition with partial pivoting and row interchanges,
C     A = P * L * U, is used, where P is a permutation matrix, L is unit
C     lower triangular, and U is upper triangular.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANS   CHARACTER*1
C             Specifies the form of op(A) to be used as follows:
C             = 'N':  op(A) = A;
C             = 'T':  op(A) = A';
C             = 'C':  op(A) = A'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix B.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix B, and the order of
C             the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the coefficient matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the factors L and U from the factorization A = P*L*U;
C             the unit diagonal elements of L are not stored.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     IPIV    (output) INTEGER array, dimension (N)
C             The pivot indices that define the permutation matrix P;
C             row i of the matrix was interchanged with row IPIV(i).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right hand side matrix B.
C             On exit, if INFO = 0, the leading M-by-N part of this
C             array contains the solution matrix X.
C
C     LDB     (input) INTEGER
C             The leading dimension of the array B.  LDB >= max(1,M).
C
C     INFO    (output) INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, U(i,i) is exactly zero.  The
C                   factorization has been completed, but the factor U
C                   is exactly singular, so the solution could not be
C                   computed.
C
C     METHOD
C
C     The LU decomposition with partial pivoting and row interchanges is
C     used to factor A as
C        A = P * L * U,
C     where P is a permutation matrix, L is unit lower triangular, and
C     U is upper triangular.  The factored form of A is then used to
C     solve the system of equations X * A = B or X * A' = B.
C
C     FURTHER COMMENTS
C
C     This routine enables to solve the system X * A = B or X * A' = B
C     as easily and efficiently as possible; it is similar to the LAPACK
C     Library routine DGESV, which solves A * X = B.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, linear algebra.
C
C    ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, M, N
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            TRAN
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGETRF, DTRSM, MA02GD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input parameters.
C
      INFO = 0
      TRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      IF( .NOT.TRAN .AND. .NOT.LSAME( TRANS, 'N' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -8
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02VD', -INFO )
         RETURN
      END IF
C
C     Compute the LU factorization of A.
C
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
C
      IF( INFO.EQ.0 ) THEN
         IF( TRAN ) THEN
C
C           Compute X = B * A**(-T).
C
            CALL MA02GD( M, B, LDB, 1, N, IPIV, 1 )
            CALL DTRSM(  'Right', 'Lower', 'Transpose', 'Unit', M, N,
     $                   ONE, A, LDA, B, LDB )
            CALL DTRSM(  'Right', 'Upper', 'Transpose', 'NonUnit', M,
     $                   N, ONE, A, LDA, B, LDB )
         ELSE
C
C           Compute X = B * A**(-1).
C
            CALL DTRSM(  'Right', 'Upper', 'NoTranspose', 'NonUnit', M,
     $                   N, ONE, A, LDA, B, LDB )
            CALL DTRSM(  'Right', 'Lower', 'NoTranspose', 'Unit', M, N,
     $                   ONE, A, LDA, B, LDB )
            CALL MA02GD( M, B, LDB, 1, N, IPIV, -1 )
         END IF
      END IF
      RETURN
C
C *** Last line of MB02VD ***
      END
