      SUBROUTINE MB03PD( JOBRQ, M, N, A, LDA, JPVT, RCOND, SVLMAX, TAU,
     $                   RANK, SVAL, DWORK, INFO )
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
C     To compute (optionally) a rank-revealing RQ factorization of a
C     real general M-by-N matrix  A,  which may be rank-deficient,
C     and estimate its effective rank using incremental condition
C     estimation.
C
C     The routine uses an RQ factorization with row pivoting:
C        P * A = R * Q,  where  R = [ R11 R12 ],
C                                   [  0  R22 ]
C     with R22 defined as the largest trailing submatrix whose estimated
C     condition number is less than 1/RCOND.  The order of R22, RANK,
C     is the effective rank of A.
C
C     MB03PD  does not perform any scaling of the matrix A.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBRQ   CHARACTER*1
C             = 'R':  Perform an RQ factorization with row pivoting;
C             = 'N':  Do not perform the RQ factorization (but assume
C                     that it has been done outside).
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C             ( LDA, N )
C             On entry with JOBRQ = 'R', the leading M-by-N part of this
C             array must contain the given matrix A.
C             On exit with JOBRQ = 'R',
C             if M <= N, the upper triangle of the subarray
C             A(1:M,N-M+1:N) contains the M-by-M upper triangular
C             matrix R;
C             if M >= N, the elements on and above the (M-N)-th
C             subdiagonal contain the M-by-N upper trapezoidal matrix R;
C             the remaining elements, with the array TAU, represent the
C             orthogonal matrix Q as a product of min(M,N) elementary
C             reflectors (see METHOD).
C             On entry and on exit with JOBRQ = 'N',
C             if M <= N, the upper triangle of the subarray
C             A(1:M,N-M+1:N) must contain the M-by-M upper triangular
C             matrix R;
C             if M >= N, the elements on and above the (M-N)-th
C             subdiagonal must contain the M-by-N upper trapezoidal
C             matrix R;
C             the remaining elements are not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     JPVT    (input/output) INTEGER array, dimension ( M )
C             On entry with JOBRQ = 'R', if JPVT(i) <> 0, the i-th row
C             of A is a final row, otherwise it is a free row. Before
C             the RQ factorization of A, all final rows are permuted
C             to the trailing positions; only the remaining free rows
C             are moved as a result of row pivoting during the
C             factorization.  For rank determination it is preferable
C             that all rows be free.
C             On exit with JOBRQ = 'R', if JPVT(i) = k, then the i-th
C             row of P*A was the k-th row of A.
C             Array JPVT is not referenced when JOBRQ = 'N'.
C
C     RCOND   (input) DOUBLE PRECISION
C             RCOND is used to determine the effective rank of A, which
C             is defined as the order of the largest trailing triangular
C             submatrix R22 in the RQ factorization with pivoting of A,
C             whose estimated condition number is less than 1/RCOND.
C             RCOND >= 0.
C             NOTE that when SVLMAX > 0, the estimated rank could be
C             less than that defined above (see SVLMAX).
C
C     SVLMAX  (input) DOUBLE PRECISION
C             If A is a submatrix of another matrix B, and the rank
C             decision should be related to that matrix, then SVLMAX
C             should be an estimate of the largest singular value of B
C             (for instance, the Frobenius norm of B).  If this is not
C             the case, the input value SVLMAX = 0 should work.
C             SVLMAX >= 0.
C
C     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) )
C             On exit with JOBRQ = 'R', the leading min(M,N) elements of
C             TAU contain the scalar factors of the elementary
C             reflectors.
C             Array TAU is not referenced when JOBRQ = 'N'.
C
C     RANK    (output) INTEGER
C             The effective (estimated) rank of A, i.e. the order of
C             the submatrix R22.
C
C     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
C             The estimates of some of the singular values of the
C             triangular factor R:
C             SVAL(1): largest singular value of
C                      R(M-RANK+1:M,N-RANK+1:N);
C             SVAL(2): smallest singular value of
C                      R(M-RANK+1:M,N-RANK+1:N);
C             SVAL(3): smallest singular value of R(M-RANK:M,N-RANK:N),
C                      if RANK < MIN( M, N ), or of
C                      R(M-RANK+1:M,N-RANK+1:N), otherwise.
C             If the triangular factorization is a rank-revealing one
C             (which will be the case if the trailing rows were well-
C             conditioned), then SVAL(1) will also be an estimate for
C             the largest singular value of A, and SVAL(2) and SVAL(3)
C             will be estimates for the RANK-th and (RANK+1)-st singular
C             values of A, respectively.
C             By examining these values, one can confirm that the rank
C             is well defined with respect to the chosen value of RCOND.
C             The ratio SVAL(1)/SVAL(2) is an estimate of the condition
C             number of R(M-RANK+1:M,N-RANK+1:N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension ( LDWORK )
C             where LDWORK = max( 1, 3*M ),           if JOBRQ = 'R';
C                   LDWORK = max( 1, 3*min( M, N ) ), if JOBRQ = 'N'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The routine computes or uses an RQ factorization with row
C     pivoting of A,  P * A = R * Q,  with  R  defined above, and then
C     finds the largest trailing submatrix whose estimated condition
C     number is less than 1/RCOND, taking the possible positive value of
C     SVLMAX into account.  This is performed using an adaptation of the
C     LAPACK incremental condition estimation scheme and a slightly
C     modified rank decision test.
C
C     The matrix Q is represented as a product of elementary reflectors
C
C        Q = H(1) H(2) . . . H(k), where k = min(m,n).
C
C     Each H(i) has the form
C
C        H = I - tau * v * v'
C
C     where tau is a real scalar, and v is a real vector with
C     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit
C     in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
C
C     The matrix P is represented in jpvt as follows: If
C        jpvt(j) = i
C     then the jth row of P is the ith canonical unit vector.
C
C     REFERENCES
C
C     [1] Bischof, C.H. and P. Tang.
C         Generalizing Incremental Condition Estimation.
C         LAPACK Working Notes 32, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-132,
C         May 1991.
C
C     [2] Bischof, C.H. and P. Tang.
C         Robust Incremental Condition Estimation.
C         LAPACK Working Notes 33, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-133,
C         May 1991.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C
C     REVISIONS
C
C     Nov. 1997
C
C     KEYWORDS
C
C     Eigenvalue problem, matrix operations, orthogonal transformation,
C     singular values.
C
C    ******************************************************************
C
C     .. Parameters ..
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          JOBRQ
      INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
C     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), SVAL( 3 ), TAU( * ), DWORK( * )
C     .. Local Scalars ..
      LOGICAL            LJOBRQ
      INTEGER            I, ISMAX, ISMIN, JWORK, MN
      DOUBLE PRECISION   C1, C2, S1, S2, SMAX, SMAXPR, SMIN, SMINPR
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAIC1, MB04GD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      LJOBRQ = LSAME( JOBRQ, 'R' )
      MN = MIN( M, N )
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( .NOT.LJOBRQ .AND. .NOT.LSAME( JOBRQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( RCOND.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -8
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MN.EQ.0 ) THEN
         RANK = 0
         SVAL( 1 ) = ZERO
         SVAL( 2 ) = ZERO
         SVAL( 3 ) = ZERO
         RETURN
      END IF
C
      IF ( LJOBRQ ) THEN
C
C        Compute RQ factorization with row pivoting of A:
C           P * A = R * Q
C        Workspace 3*M. Details of Householder rotations stored in TAU.
C
         CALL MB04GD( M, N, A, LDA, JPVT, TAU, DWORK( 1 ), INFO )
      END IF
C
C     Determine RANK using incremental condition estimation.
C        Workspace 3*min(M,N).
C
      SMAX = ABS( A( M, N ) )
      IF( SMAX.EQ.ZERO .OR. SVLMAX*RCOND.GT.SMAX ) THEN
         RANK = 0
         SVAL( 1 ) = SMAX
         SVAL( 2 ) = ZERO
         SVAL( 3 ) = ZERO
      ELSE
         ISMIN = MN
         ISMAX = 2*MN
         JWORK = ISMAX + 1
         DWORK( ISMIN ) = ONE
         DWORK( ISMAX ) = ONE
         RANK = 1
         SMIN = SMAX
         SMINPR = SMIN
C
   10    CONTINUE
         IF( RANK.LT.MN ) THEN
            CALL DCOPY ( RANK, A( M-RANK, N-RANK+1 ), LDA,
     $                   DWORK( JWORK ), 1 )
            CALL DLAIC1( IMIN, RANK, DWORK( ISMIN ), SMIN,
     $                   DWORK( JWORK ), A( M-RANK, N-RANK ), SMINPR,
     $                   S1, C1 )
            CALL DLAIC1( IMAX, RANK, DWORK( ISMAX ), SMAX,
     $                   DWORK( JWORK ), A( M-RANK, N-RANK ), SMAXPR,
     $                   S2, C2 )
C
            IF( SVLMAX*RCOND.LE.SMAXPR ) THEN
               IF( SVLMAX*RCOND.LE.SMINPR ) THEN
                  IF( SMAXPR*RCOND.LE.SMINPR ) THEN
                     DO 20 I = 1, RANK
                        DWORK( ISMIN+I-1 ) = S1*DWORK( ISMIN+I-1 )
                        DWORK( ISMAX+I-1 ) = S2*DWORK( ISMAX+I-1 )
   20                CONTINUE
                     ISMIN = ISMIN - 1
                     ISMAX = ISMAX - 1
                     DWORK( ISMIN ) = C1
                     DWORK( ISMAX ) = C2
                     SMIN = SMINPR
                     SMAX = SMAXPR
                     RANK = RANK + 1
                     GO TO 10
                  END IF
               END IF
            END IF
         END IF
         SVAL( 1 ) = SMAX
         SVAL( 2 ) = SMIN
         SVAL( 3 ) = SMINPR
      END IF
C
      RETURN
C *** Last line of MB03PD ***
      END
