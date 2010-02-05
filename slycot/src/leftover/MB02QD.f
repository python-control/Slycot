      SUBROUTINE MB02QD( JOB, INIPER, M, N, NRHS, RCOND, SVLMAX, A, LDA,
     $                   B, LDB, Y, JPVT, RANK, SVAL, DWORK, LDWORK,
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
C     To compute a solution, optionally corresponding to specified free
C     elements, to a real linear least squares problem:
C
C         minimize || A * X - B ||
C
C     using a complete orthogonal factorization of the M-by-N matrix A,
C     which may be rank-deficient.
C
C     Several right hand side vectors b and solution vectors x can be
C     handled in a single call; they are stored as the columns of the
C     M-by-NRHS right hand side matrix B and the N-by-NRHS solution
C     matrix X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies whether or not a standard least squares solution
C             must be computed, as follows:
C             = 'L':  Compute a standard least squares solution (Y = 0);
C             = 'F':  Compute a solution with specified free elements
C                     (given in Y).
C
C     INIPER  CHARACTER*1
C             Specifies whether an initial column permutation, defined
C             by JPVT, must be performed, as follows:
C             = 'P':  Perform an initial column permutation;
C             = 'N':  Do not perform an initial column permutation.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     NRHS    (input) INTEGER
C             The number of right hand sides, i.e., the number of
C             columns of the matrices B and X.  NRHS >= 0.
C
C     RCOND   (input) DOUBLE PRECISION
C             RCOND is used to determine the effective rank of A, which
C             is defined as the order of the largest leading triangular
C             submatrix R11 in the QR factorization with pivoting of A,
C             whose estimated condition number is less than 1/RCOND.
C             0 <= RCOND <= 1.
C
C     SVLMAX  (input) DOUBLE PRECISION
C             If A is a submatrix of another matrix C, and the rank
C             decision should be related to that matrix, then SVLMAX
C             should be an estimate of the largest singular value of C
C             (for instance, the Frobenius norm of C).  If this is not
C             the case, the input value SVLMAX = 0 should work.
C             SVLMAX >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the given matrix A.
C             On exit, the leading M-by-N part of this array contains
C             details of its complete orthogonal factorization:
C             the leading RANK-by-RANK upper triangular part contains
C             the upper triangular factor T11 (see METHOD);
C             the elements below the diagonal, with the entries 2 to
C             min(M,N)+1 of the array DWORK, represent the orthogonal
C             matrix Q as a product of min(M,N) elementary reflectors
C             (see METHOD);
C             the elements of the subarray A(1:RANK,RANK+1:N), with the
C             next RANK entries of the array DWORK, represent the
C             orthogonal matrix Z as a product of RANK elementary
C             reflectors (see METHOD).
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,NRHS)
C             On entry, the leading M-by-NRHS part of this array must
C             contain the right hand side matrix B.
C             On exit, the leading N-by-NRHS part of this array contains
C             the solution matrix X.
C             If M >= N and RANK = N, the residual sum-of-squares for
C             the solution in the i-th column is given by the sum of
C             squares of elements N+1:M in that column.
C             If NRHS = 0, this array is not referenced, and the routine
C             returns the effective rank of A, and its QR factorization.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,M,N).
C
C     Y       (input) DOUBLE PRECISION array, dimension ( N*NRHS )
C             If JOB = 'F', the elements Y(1:(N-RANK)*NRHS) are used as
C             free elements in computing the solution (see METHOD).
C             The remaining elements are not referenced.
C             If JOB = 'L', or NRHS = 0, this array is not referenced.
C
C     JPVT    (input/output) INTEGER array, dimension (N)
C             On entry with INIPER = 'P', if JPVT(i) <> 0, the i-th
C             column of A is an initial column, otherwise it is a free
C             column.  Before the QR factorization of A, all initial
C             columns are permuted to the leading positions; only the
C             remaining free columns are moved as a result of column
C             pivoting during the factorization.
C             If INIPER = 'N', JPVT need not be set on entry.
C             On exit, if JPVT(i) = k, then the i-th column of A*P
C             was the k-th column of A.
C
C     RANK    (output) INTEGER
C             The effective rank of A, i.e., the order of the submatrix
C             R11.  This is the same as the order of the submatrix T11
C             in the complete orthogonal factorization of A.
C
C     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
C             The estimates of some of the singular values of the
C             triangular factor R11:
C             SVAL(1): largest singular value of  R(1:RANK,1:RANK);
C             SVAL(2): smallest singular value of R(1:RANK,1:RANK);
C             SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1),
C                      if RANK < MIN( M, N ), or of R(1:RANK,1:RANK),
C                      otherwise.
C             If the triangular factorization is a rank-revealing one
C             (which will be the case if the leading columns were well-
C             conditioned), then SVAL(1) will also be an estimate for
C             the largest singular value of A, and SVAL(2) and SVAL(3)
C             will be estimates for the RANK-th and (RANK+1)-st singular
C             values of A, respectively.
C             By examining these values, one can confirm that the rank
C             is well defined with respect to the chosen value of RCOND.
C             The ratio SVAL(1)/SVAL(2) is an estimate of the condition
C             number of R(1:RANK,1:RANK).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension LDWORK
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, and the entries 2 to min(M,N) + RANK + 1
C             contain the scalar factors of the elementary reflectors
C             used in the complete orthogonal factorization of A.
C             Among the entries 2 to min(M,N) + 1, only the first RANK
C             elements are useful, if INIPER = 'N'.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= max( min(M,N)+3*N+1, 2*min(M,N)+NRHS )
C             For optimum performance LDWORK should be larger.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     If INIPER = 'P', the routine first computes a QR factorization
C     with column pivoting:
C         A * P = Q * [ R11 R12 ]
C                     [  0  R22 ]
C     with R11 defined as the largest leading submatrix whose estimated
C     condition number is less than 1/RCOND.  The order of R11, RANK,
C     is the effective rank of A.
C     If INIPER = 'N', the effective rank is estimated during a
C     truncated QR factorization (with column pivoting) process, and
C     the submatrix R22 is not upper triangular, but full and of small
C     norm. (See SLICOT Library routines MB03OD or MB03OY, respectively,
C     for further details.)
C
C     Then, R22 is considered to be negligible, and R12 is annihilated
C     by orthogonal transformations from the right, arriving at the
C     complete orthogonal factorization:
C        A * P = Q * [ T11 0 ] * Z
C                    [  0  0 ]
C     The solution is then
C        X = P * Z' [ inv(T11)*Q1'*B ]
C                   [        Y       ]
C     where Q1 consists of the first RANK columns of Q, and Y contains
C     free elements (if JOB = 'F'), or is zero (if JOB = 'L').
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     FURTHER COMMENTS
C
C     Significant gain in efficiency is possible for small-rank problems
C     using truncated QR factorization (option INIPER = 'N').
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, Technical University of Sofia, Oct. 1998,
C     modification of the LAPACK routine DGELSX.
C     V. Sima, Katholieke Universiteit Leuven, Jan. 1999, SLICOT Library
C     version.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005.
C
C     KEYWORDS
C
C     Least squares problems, QR factorization.
C
C     ******************************************************************
C
      DOUBLE PRECISION   ZERO, ONE, DONE, NTDONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, DONE = ZERO,
     $                     NTDONE = ONE )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          INIPER, JOB
      INTEGER            INFO, LDA, LDB, LDWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
C     ..
C     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), DWORK( * ),
     $                   SVAL( 3 ), Y ( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            LEASTS, PERMUT
      INTEGER            I, IASCL, IBSCL, J, K, MAXWRK, MINWRK, MN
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM, T1, T2
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLABAD, DLACPY, DLASCL, DLASET, DORMQR, DORMRZ,
     $                   DTRSM, DTZRZF, MB03OD, MB03OY, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      MN = MIN( M, N )
      LEASTS = LSAME( JOB, 'L' )
      PERMUT = LSAME( INIPER, 'P' )
C
C     Test the input scalar arguments.
C
      INFO   = 0
      MINWRK = MAX( MN + 3*N + 1, 2*MN + NRHS )
      IF( .NOT. ( LEASTS .OR. LSAME( JOB, 'F' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( PERMUT .OR. LSAME( INIPER, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( RCOND.LT.ZERO .OR. RCOND.GT.ONE ) THEN
         INFO = -6
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -11
      ELSE IF( LDWORK.LT.MINWRK ) THEN
         INFO = -17
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02QD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MN.EQ.0 ) THEN
         RANK = 0
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Get machine parameters.
C
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
C     Scale A, B if max entries outside range [SMLNUM,BIGNUM].
C
      ANRM = DLANGE( 'M', M, N, A, LDA, DWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
C
C        Scale matrix norm up to SMLNUM.
C
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
C
C        Scale matrix norm down to BIGNUM.
C
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
C
C        Matrix all zero. Return zero solution.
C
         IF( NRHS.GT.0 )
     $      CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RANK = 0
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
      IF( NRHS.GT.0 ) THEN
         BNRM = DLANGE( 'M', M, NRHS, B, LDB, DWORK )
         IBSCL = 0
         IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
C
C           Scale matrix norm up to SMLNUM.
C
            CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB,
     $                   INFO )
            IBSCL = 1
         ELSE IF( BNRM.GT.BIGNUM ) THEN
C
C           Scale matrix norm down to BIGNUM.
C
            CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB,
     $                   INFO )
            IBSCL = 2
         END IF
      END IF
C
C     Compute a rank-revealing QR factorization of A and estimate its
C     effective rank using incremental condition estimation:
C        A * P = Q * R.
C     Workspace need   min(M,N)+3*N+1;
C               prefer min(M,N)+2*N+N*NB.
C     Details of Householder transformations stored in DWORK(1:MN).
C     (Note: Comments in the code beginning "Workspace:" describe the
C      minimal amount of workspace needed at that point in the code,
C      as well as the preferred amount for good performance.
C      NB refers to the optimal block size for the immediately
C      following subroutine, as returned by ILAENV.)
C
      MAXWRK = MINWRK
      IF( PERMUT ) THEN
         CALL MB03OD( 'Q', M, N, A, LDA, JPVT, RCOND, SVLMAX,
     $                DWORK( 1 ), RANK, SVAL, DWORK( MN+1 ), LDWORK-MN,
     $                INFO )
         MAXWRK = MAX( MAXWRK, INT( DWORK( MN+1 ) ) + MN )
      ELSE
         CALL MB03OY( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
     $                DWORK( 1 ), DWORK( MN+1 ), INFO )
      END IF
C
C     Logically partition R = [ R11 R12 ]
C                             [  0  R22 ],
C     where R11 = R(1:RANK,1:RANK).
C
C     [R11,R12] = [ T11, 0 ] * Z.
C
C     Details of Householder transformations stored in DWORK(MN+1:2*MN).
C     Workspace need   3*min(M,N);
C               prefer 2*min(M,N)+min(M,N)*NB.
C
      IF( RANK.LT.N ) THEN
         CALL DTZRZF( RANK, N, A, LDA, DWORK( MN+1 ), DWORK( 2*MN+1 ),
     $                LDWORK-2*MN, INFO )
         MAXWRK = MAX( MAXWRK, INT( DWORK( 2*MN+1 ) ) + 2*MN )
      END IF
C
      IF( NRHS.GT.0 ) THEN
C
C        B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS).
C
C        Workspace: need   2*min(M,N)+NRHS;
C                   prefer   min(M,N)+NRHS*NB.
C
         CALL DORMQR( 'Left', 'Transpose', M, NRHS, MN, A, LDA,
     $                DWORK( 1 ), B, LDB, DWORK( 2*MN+1 ), LDWORK-2*MN,
     $                INFO )
         MAXWRK = MAX( MAXWRK, INT( DWORK( 2*MN+1 ) ) + 2*MN )
C
C        B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS).
C
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', RANK,
     $               NRHS, ONE, A, LDA, B, LDB )
C
         IF( RANK.LT.N ) THEN
C
C           Set B(RANK+1:N,1:NRHS).
C
            IF( LEASTS ) THEN
               CALL DLASET( 'Full', N-RANK, NRHS, ZERO, ZERO,
     $                      B(RANK+1,1), LDB )
            ELSE
               CALL DLACPY( 'Full', N-RANK, NRHS, Y, N-RANK,
     $                      B(RANK+1,1), LDB )
            END IF
C
C           B(1:N,1:NRHS) := Z' * B(1:N,1:NRHS).
C
C           Workspace need   2*min(M,N)+NRHS;
C                     prefer 2*min(M,N)+NRHS*NB.
C
            CALL DORMRZ( 'Left', 'Transpose', N, NRHS, RANK, N-RANK, A,
     $                   LDA, DWORK( MN+1 ), B, LDB, DWORK( 2*MN+1 ),
     $                   LDWORK-2*MN, INFO )
            MAXWRK = MAX( MAXWRK, INT( DWORK( 2*MN+1 ) ) + 2*MN )
         END IF
C
C        Additional workspace: NRHS.
C
C        B(1:N,1:NRHS) := P * B(1:N,1:NRHS).
C
         DO 50 J = 1, NRHS
            DO 20 I = 1, N
               DWORK( 2*MN+I ) = NTDONE
   20       CONTINUE
            DO 40 I = 1, N
               IF( DWORK( 2*MN+I ).EQ.NTDONE ) THEN
                  IF( JPVT( I ).NE.I ) THEN
                     K  = I
                     T1 = B( K, J )
                     T2 = B( JPVT( K ), J )
   30                CONTINUE
                     B( JPVT( K ), J ) = T1
                     DWORK( 2*MN+K ) = DONE
                     T1 = T2
                     K  = JPVT( K )
                     T2 = B( JPVT( K ), J )
                     IF( JPVT( K ).NE.I )
     $                  GO TO 30
                     B( I, J ) = T1
                     DWORK( 2*MN+K ) = DONE
                  END IF
               END IF
   40       CONTINUE
   50    CONTINUE
C
C        Undo scaling for B.
C
         IF( IBSCL.EQ.1 ) THEN
            CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB,
     $                   INFO )
         ELSE IF( IBSCL.EQ.2 ) THEN
            CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB,
     $                   INFO )
         END IF
      END IF
C
C     Undo scaling for A.
C
      IF( IASCL.EQ.1 ) THEN
         IF( NRHS.GT.0 )
     $      CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB,
     $                   INFO )
         CALL DLASCL( 'U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         IF( NRHS.GT.0 )
     $      CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB,
     $                   INFO )
         CALL DLASCL( 'U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      END IF
C
      DO 60 I = MN + RANK, 1, -1
         DWORK( I+1 ) = DWORK( I )
   60 CONTINUE
C
      DWORK( 1 ) = MAXWRK
      RETURN
C *** Last line of MB02QD ***
      END
