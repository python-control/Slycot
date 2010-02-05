      SUBROUTINE MB03PY( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
     $                   TAU, DWORK, INFO )
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
C     To compute a rank-revealing RQ factorization of a real general
C     M-by-N matrix  A,  which may be rank-deficient, and estimate its
C     effective rank using incremental condition estimation.
C
C     The routine uses a truncated RQ factorization with row pivoting:
C                                   [ R11 R12 ]
C        P * A = R * Q,  where  R = [         ],
C                                   [  0  R22 ]
C     with R22 defined as the largest trailing upper triangular
C     submatrix whose estimated condition number is less than 1/RCOND.
C     The order of R22, RANK, is the effective rank of A.  Condition
C     estimation is performed during the RQ factorization process.
C     Matrix R11 is full (but of small norm), or empty.
C
C     MB03PY  does not perform any scaling of the matrix A.
C
C     ARGUMENTS
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
C             On entry, the leading M-by-N part of this array must
C             contain the given matrix A.
C             On exit, the upper triangle of the subarray
C             A(M-RANK+1:M,N-RANK+1:N) contains the RANK-by-RANK upper
C             triangular matrix R22;  the remaining elements in the last
C             RANK  rows, with the array TAU, represent the orthogonal
C             matrix Q as a product of  RANK  elementary reflectors
C             (see METHOD).  The first  M-RANK  rows contain the result
C             of the RQ factorization process used.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     RCOND   (input) DOUBLE PRECISION
C             RCOND is used to determine the effective rank of A, which
C             is defined as the order of the largest trailing triangular
C             submatrix R22 in the RQ factorization with pivoting of A,
C             whose estimated condition number is less than 1/RCOND.
C             0 <= RCOND <= 1.
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
C     RANK    (output) INTEGER
C             The effective (estimated) rank of A, i.e., the order of
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
C     JPVT    (output) INTEGER array, dimension ( M )
C             If JPVT(i) = k, then the i-th row of P*A was the k-th row
C             of A.
C
C     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) )
C             The trailing  RANK  elements of TAU contain the scalar
C             factors of the elementary reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension ( 3*M-1 )
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
C     The routine computes a truncated RQ factorization with row
C     pivoting of A,  P * A = R * Q,  with  R  defined above, and,
C     during this process, finds the largest trailing submatrix whose
C     estimated condition number is less than 1/RCOND, taking the
C     possible positive value of SVLMAX into account.  This is performed
C     using an adaptation of the LAPACK incremental condition estimation
C     scheme and a slightly modified rank decision test.  The
C     factorization process stops when  RANK  has been determined.
C
C     The matrix Q is represented as a product of elementary reflectors
C
C        Q = H(k-rank+1) H(k-rank+2) . . . H(k), where k = min(m,n).
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
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001,
C     Jan. 2009.
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
      DOUBLE PRECISION   ZERO, ONE, P05
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, P05 = 0.05D+0 )
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
C     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), DWORK( * ), SVAL( 3 ), TAU( * )
C     .. Local Scalars ..
      INTEGER            I, ISMAX, ISMIN, ITEMP, J, JWORK, K, MKI, NKI,
     $                   PVT
      DOUBLE PRECISION   AII, C1, C2, S1, S2, SMAX, SMAXPR, SMIN,
     $                   SMINPR, TEMP, TEMP2
C     ..
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DNRM2
      EXTERNAL           DNRM2, IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAIC1, DLARF, DLARFG, DSCAL, DSWAP,
     $                   XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( RCOND.LT.ZERO .OR. RCOND.GT.ONE ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03PY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         RANK = 0
         SVAL( 1 ) = ZERO
         SVAL( 2 ) = ZERO
         SVAL( 3 ) = ZERO
         RETURN
      END IF
C
      ISMIN = M
      ISMAX = ISMIN + M
      JWORK = ISMAX + 1
C
C     Initialize partial row norms and pivoting vector. The first m
C     elements of DWORK store the exact row norms. The already used
C     trailing part is then overwritten by the condition estimator.
C
      DO 10 I = 1, M
         DWORK( I ) = DNRM2( N, A( I, 1 ), LDA )
         DWORK( M+I ) = DWORK( I )
         JPVT( I ) = I
   10 CONTINUE
C
C     Compute factorization and determine RANK using incremental
C     condition estimation.
C
      RANK = 0
C
   20 CONTINUE
      IF( RANK.LT.K ) THEN
         I = K - RANK
C
C        Determine ith pivot row and swap if necessary.
C
         MKI = M - RANK
         NKI = N - RANK
         PVT = IDAMAX( MKI, DWORK, 1 )
C
         IF( PVT.NE.MKI ) THEN
            CALL DSWAP( N, A( PVT, 1 ), LDA, A( MKI, 1 ), LDA )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( MKI )
            JPVT( MKI ) = ITEMP
            DWORK( PVT )   = DWORK( MKI )
            DWORK( M+PVT ) = DWORK( M+MKI )
         END IF
C
         IF( NKI.GT.1 ) THEN
C
C           Save A(m-k+i,n-k+i) and generate elementary reflector H(i)
C           to annihilate A(m-k+i,1:n-k+i-1), k = min(m,n).
C
            AII = A( MKI, NKI )
            CALL DLARFG( NKI, A( MKI, NKI ), A( MKI, 1 ), LDA, TAU( I )
     $                 )
         END IF
C
         IF( RANK.EQ.0 ) THEN
C
C           Initialize; exit if matrix is zero (RANK = 0).
C
            SMAX = ABS( A( M, N ) )
            IF ( SMAX.EQ.ZERO ) THEN
               SVAL( 1 ) = ZERO
               SVAL( 2 ) = ZERO
               SVAL( 3 ) = ZERO
               RETURN
            END IF
            SMIN = SMAX
            SMAXPR = SMAX
            SMINPR = SMIN
            C1 = ONE
            C2 = ONE
         ELSE
C
C           One step of incremental condition estimation.
C
            CALL DCOPY ( RANK, A( MKI, NKI+1 ), LDA, DWORK( JWORK ), 1 )
            CALL DLAIC1( IMIN, RANK, DWORK( ISMIN ), SMIN,
     $                   DWORK( JWORK ), A( MKI, NKI ), SMINPR, S1, C1 )
            CALL DLAIC1( IMAX, RANK, DWORK( ISMAX ), SMAX,
     $                   DWORK( JWORK ), A( MKI, NKI ), SMAXPR, S2, C2 )
         END IF
C
         IF( SVLMAX*RCOND.LE.SMAXPR ) THEN
            IF( SVLMAX*RCOND.LE.SMINPR ) THEN
               IF( SMAXPR*RCOND.LE.SMINPR ) THEN
C
                  IF( MKI.GT.1 ) THEN
C
C                    Continue factorization, as rank is at least RANK.
C                    Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right.
C
                     AII = A( MKI, NKI )
                     A( MKI, NKI ) = ONE
                     CALL DLARF( 'Right', MKI-1, NKI, A( MKI, 1 ), LDA,
     $                           TAU( I ), A, LDA, DWORK( JWORK ) )
                     A( MKI, NKI ) = AII
C
C                    Update partial row norms.
C
                     DO 30 J = 1, MKI - 1
                        IF( DWORK( J ).NE.ZERO ) THEN
                           TEMP  = ONE -
     $                             ( ABS( A( J, NKI ) )/DWORK( J ) )**2
                           TEMP  = MAX( TEMP, ZERO )
                           TEMP2 = ONE + P05*TEMP*
     $                                 ( DWORK( J ) / DWORK( M+J ) )**2
                           IF( TEMP2.EQ.ONE ) THEN
                              DWORK( J )   = DNRM2( NKI-1, A( J, 1 ),
     $                                              LDA )
                              DWORK( M+J ) = DWORK( J )
                           ELSE
                              DWORK( J ) = DWORK( J )*SQRT( TEMP )
                           END IF
                        END IF
   30                CONTINUE
C
                  END IF
C
                  DO 40 I = 1, RANK
                     DWORK( ISMIN+I-1 ) = S1*DWORK( ISMIN+I-1 )
                     DWORK( ISMAX+I-1 ) = S2*DWORK( ISMAX+I-1 )
   40             CONTINUE
C
                  IF( RANK.GT.0 ) THEN
                     ISMIN = ISMIN - 1
                     ISMAX = ISMAX - 1
                  END IF
                  DWORK( ISMIN ) = C1
                  DWORK( ISMAX ) = C2
                  SMIN = SMINPR
                  SMAX = SMAXPR
                  RANK = RANK + 1
                  GO TO 20
               END IF
            END IF
         END IF
      END IF
C
C     Restore the changed part of the (M-RANK)-th row and set SVAL.
C
      IF ( RANK.LT.K .AND. NKI.GT.1 ) THEN
         CALL DSCAL( NKI-1, -A( MKI, NKI )*TAU( I ), A( MKI, 1 ), LDA )
         A( MKI, NKI ) = AII
      END IF
      SVAL( 1 ) = SMAX
      SVAL( 2 ) = SMIN
      SVAL( 3 ) = SMINPR
C
      RETURN
C *** Last line of MB03PY ***
      END
