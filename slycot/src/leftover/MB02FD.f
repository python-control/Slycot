      SUBROUTINE MB02FD( TYPET, K, N, P, S, T, LDT, R, LDR, DWORK,
     $                   LDWORK, INFO )
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
C     To compute the incomplete Cholesky (ICC) factor of a symmetric
C     positive definite (s.p.d.) block Toeplitz matrix T, defined by
C     either its first block row, or its first block column, depending
C     on the routine parameter TYPET.
C
C     By subsequent calls of this routine, further rows / columns of
C     the Cholesky factor can be added.
C     Furthermore, the generator of the Schur complement of the leading
C     (P+S)*K-by-(P+S)*K block in T is available, which can be used,
C     e.g., for measuring the quality of the ICC factorization.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPET   CHARACTER*1
C             Specifies the type of T, as follows:
C             = 'R':  T contains the first block row of an s.p.d. block
C                     Toeplitz matrix; the ICC factor R is upper
C                     trapezoidal;
C             = 'C':  T contains the first block column of an s.p.d.
C                     block Toeplitz matrix; the ICC factor R is lower
C                     trapezoidal; this choice leads to better
C                     localized memory references and hence a faster
C                     algorithm.
C             Note:   in the sequel, the notation x / y means that
C                     x corresponds to TYPET = 'R' and y corresponds to
C                     TYPET = 'C'.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of rows / columns in T, which should be equal
C             to the blocksize.  K >= 0.
C
C     N       (input)  INTEGER
C             The number of blocks in T.  N >= 0.
C
C     P       (input)  INTEGER
C             The number of previously computed block rows / columns
C             of R.  0 <= P <= N.
C
C     S       (input)  INTEGER
C             The number of block rows / columns of R to compute.
C             0 <= S <= N-P.
C
C     T       (input/output)  DOUBLE PRECISION array, dimension
C             (LDT,(N-P)*K) / (LDT,K)
C             On entry, if P = 0, then the leading K-by-N*K / N*K-by-K
C             part of this array must contain the first block row /
C             column of an s.p.d. block Toeplitz matrix.
C             If P > 0, the leading K-by-(N-P)*K / (N-P)*K-by-K must
C             contain the negative generator of the Schur complement of
C             the leading P*K-by-P*K part in T, computed from previous
C             calls of this routine.
C             On exit, if INFO = 0, then the leading K-by-(N-P)*K /
C             (N-P)*K-by-K part of this array contains, in the first
C             K-by-K block, the upper / lower Cholesky factor of
C             T(1:K,1:K), in the following S-1 K-by-K blocks, the
C             Householder transformations applied during the process,
C             and in the remaining part, the negative generator of the
C             Schur complement of the leading (P+S)*K-by(P+S)*K part
C             in T.
C
C     LDT     INTEGER
C             The leading dimension of the array T.
C             LDT >= MAX(1,K),        if TYPET = 'R';
C             LDT >= MAX(1,(N-P)*K),  if TYPET = 'C'.
C
C     R       (input/output)  DOUBLE PRECISION array, dimension
C             (LDR, N*K)       / (LDR, S*K )     if P = 0;
C             (LDR, (N-P+1)*K) / (LDR, (S+1)*K ) if P > 0.
C             On entry, if P > 0, then the leading K-by-(N-P+1)*K /
C             (N-P+1)*K-by-K part of this array must contain the
C             nonzero blocks of the last block row / column in the
C             ICC factor from a previous call of this routine. Note that
C             this part is identical with the positive generator of
C             the Schur complement of the leading P*K-by-P*K part in T.
C             If P = 0, then R is only an output parameter.
C             On exit, if INFO = 0 and P = 0, then the leading
C             S*K-by-N*K / N*K-by-S*K part of this array contains the
C             upper / lower trapezoidal ICC factor.
C             On exit, if INFO = 0 and P > 0, then the leading
C             (S+1)*K-by-(N-P+1)*K / (N-P+1)*K-by-(S+1)*K part of this
C             array contains the upper / lower trapezoidal part of the
C             P-th to (P+S)-th block rows / columns of the ICC factor.
C             The elements in the strictly lower / upper trapezoidal
C             part are not referenced.
C
C     LDR     INTEGER
C             The leading dimension of the array R.
C             LDR >= MAX(1, S*K ),        if TYPET = 'R' and P = 0;
C             LDR >= MAX(1, (S+1)*K ),    if TYPET = 'R' and P > 0;
C             LDR >= MAX(1, N*K ),        if TYPET = 'C' and P = 0;
C             LDR >= MAX(1, (N-P+1)*K ),  if TYPET = 'C' and P > 0.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -11,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,(N+1)*K,4*K),   if P = 0;
C             LDWORK >= MAX(1,(N-P+2)*K,4*K), if P > 0.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction algorithm failed; the Toeplitz matrix
C                   associated with T is not (numerically) positive
C                   definite in its leading (P+S)*K-by-(P+S)*K part.
C
C     METHOD
C
C     Householder transformations and modified hyperbolic rotations
C     are used in the Schur algorithm [1], [2].
C
C     REFERENCES
C
C     [1] Kailath, T. and Sayed, A.
C         Fast Reliable Algorithms for Matrices with Structure.
C         SIAM Publications, Philadelphia, 1999.
C
C     [2] Kressner, D. and Van Dooren, P.
C         Factorizations and linear system solvers for matrices with
C         Toeplitz structure.
C         SLICOT Working Note 2000-2, 2000.
C
C     NUMERICAL ASPECTS
C
C     The implemented method is numerically stable.
C                               3
C     The algorithm requires 0(K S (N-P)) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, April 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001,
C     Mar. 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, Householder transformation, matrix
C     operations, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TYPET
      INTEGER           INFO, K, LDR, LDT, LDWORK, N, P, S
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), R(LDR,*), T(LDT,*)
C     .. Local Scalars ..
      INTEGER           COUNTR, I, IERR, MAXWRK, ST, STARTR
      LOGICAL           ISROW
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DPOTRF, DTRSM, MB02CX, MB02CY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      ISROW = LSAME( TYPET, 'R' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( ISROW .OR. LSAME( TYPET, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF ( K.LT.0 ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE IF ( P.LT.0 .OR. P.GT.N ) THEN
         INFO = -4
      ELSE IF ( S.LT.0 .OR. S.GT.( N-P ) ) THEN
         INFO = -5
      ELSE IF ( LDT.LT.1 .OR. ( ISROW .AND. LDT.LT.K ) .OR.
     $                   ( .NOT.ISROW .AND. LDT.LT.( N-P )*K ) ) THEN
         INFO = -7
      ELSE IF ( LDR.LT.1 .OR.
     $        ( ISROW .AND. P.EQ.0 .AND. ( LDR.LT.S*K ) ) .OR.
     $        ( ISROW .AND. P.GT.0 .AND. ( LDR.LT.( S+1 )*K ) ) .OR.
     $   ( .NOT.ISROW .AND. P.EQ.0 .AND. ( LDR.LT.N*K ) ) .OR.
     $   ( .NOT.ISROW .AND. P.GT.0 .AND. ( LDR.LT.( N-P+1 )*K ) ) ) THEN
        INFO = -9
      ELSE
         IF ( P.EQ.0 ) THEN
            COUNTR = ( N + 1 )*K
         ELSE
            COUNTR = ( N - P + 2 )*K
         END IF
         COUNTR = MAX( COUNTR, 4*K )
         IF ( LDWORK.LT.MAX( 1, COUNTR ) ) THEN
            DWORK(1) = MAX( 1, COUNTR )
            INFO = -11
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02FD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( K, N, S ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      MAXWRK = 1
C
      IF ( ISROW ) THEN
C
         IF ( P.EQ.0 ) THEN
C
C           T is the first block row of a block Toeplitz matrix.
C           Bring T to proper form by triangularizing its first block.
C
            CALL DPOTRF( 'Upper', K, T, LDT, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            IF ( N.GT.1 )
     $         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'NonUnit', K,
     $                     (N-1)*K, ONE, T, LDT, T(1,K+1), LDT )
            CALL DLACPY( 'Upper', K, N*K, T, LDT, R, LDR )
C
            IF ( S.EQ.1 ) THEN
               DWORK(1) = ONE
               RETURN
            END IF
C
            ST = 2
            COUNTR = ( N - 1 )*K
         ELSE
            ST = 1
            COUNTR = ( N - P )*K
         END IF
C
         STARTR = 1
C
         DO 10 I = ST, S
            CALL DLACPY( 'Upper', K, COUNTR, R(STARTR,STARTR), LDR,
     $                   R(STARTR+K,STARTR+K), LDR )
            STARTR = STARTR + K
            COUNTR = COUNTR - K
            CALL MB02CX( 'Row', K, K, K, R(STARTR,STARTR), LDR,
     $                   T(1,STARTR), LDT, DWORK, 3*K, DWORK(3*K+1),
     $                   LDWORK-3*K, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            MAXWRK = MAX( MAXWRK, INT( DWORK(3*K+1) ) + 3*K )
            CALL MB02CY( 'Row', 'NoStructure', K, K, COUNTR, K,
     $                   R(STARTR,STARTR+K), LDR, T(1,STARTR+K), LDT,
     $                   T(1,STARTR), LDT, DWORK, 3*K, DWORK(3*K+1),
     $                   LDWORK-3*K, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(3*K+1) ) + 3*K )
  10     CONTINUE
C
      ELSE
C
         IF ( P.EQ.0 ) THEN
C
C           T is the first block column of a block Toeplitz matrix.
C           Bring T to proper form by triangularizing its first block.
C
            CALL DPOTRF( 'Lower', K, T, LDT, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            IF ( N.GT.1 )
     $         CALL DTRSM( 'Right', 'Lower', 'Transpose', 'NonUnit',
     $                     (N-1)*K, K, ONE, T, LDT, T(K+1,1), LDT )
            CALL DLACPY( 'Lower', N*K, K, T, LDT, R, LDR )
C
            IF ( S.EQ.1 ) THEN
               DWORK(1) = ONE
               RETURN
            END IF
C
            ST = 2
            COUNTR = ( N - 1 )*K
         ELSE
            ST = 1
            COUNTR = ( N - P )*K
         END IF
C
         STARTR = 1
C
         DO 20 I = ST, S
            CALL DLACPY( 'Lower', COUNTR, K, R(STARTR,STARTR), LDR,
     $                   R(STARTR+K,STARTR+K), LDR )
            STARTR = STARTR + K
            COUNTR = COUNTR - K
            CALL MB02CX( 'Column', K, K, K, R(STARTR,STARTR), LDR,
     $                   T(STARTR,1), LDT, DWORK, 3*K, DWORK(3*K+1),
     $                   LDWORK-3*K, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            MAXWRK = MAX( MAXWRK, INT( DWORK(3*K+1) ) + 3*K )
            CALL MB02CY( 'Column', 'NoStructure', K, K, COUNTR, K,
     $                   R(STARTR+K,STARTR), LDR, T(STARTR+K,1), LDT,
     $                   T(STARTR,1), LDT, DWORK, 3*K, DWORK(3*K+1),
     $                   LDWORK-3*K, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(3*K+1) ) + 3*K )
  20     CONTINUE
C
      END IF
C
      DWORK(1) = MAXWRK
C
      RETURN
C
C *** Last line of MB02FD ***
      END
