      SUBROUTINE MB02ED( TYPET, K, N, NRHS, T, LDT, B, LDB, DWORK,
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
C     To solve a system of linear equations  T*X = B  or  X*T = B  with
C     a symmetric positive definite (s.p.d.) block Toeplitz matrix T.
C     T is defined either by its first block row or its first block
C     column, depending on the parameter TYPET.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPET   CHARACTER*1
C             Specifies the type of T, as follows:
C             = 'R':  T contains the first block row of an s.p.d. block
C                     Toeplitz matrix, and the system X*T = B is solved;
C             = 'C':  T contains the first block column of an s.p.d.
C                     block Toeplitz matrix, and the system T*X = B is
C                     solved.
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
C     NRHS    (input)  INTEGER
C             The number of right hand sides.  NRHS >= 0.
C
C     T       (input/output)  DOUBLE PRECISION array, dimension
C             (LDT,N*K) / (LDT,K)
C             On entry, the leading K-by-N*K / N*K-by-K part of this
C             array must contain the first block row / column of an
C             s.p.d. block Toeplitz matrix.
C             On exit, if  INFO = 0  and  NRHS > 0,  then the leading
C             K-by-N*K / N*K-by-K part of this array contains the last
C             row / column of the Cholesky factor of inv(T).
C
C     LDT     INTEGER
C             The leading dimension of the array T.
C             LDT >= MAX(1,K),    if TYPET = 'R';
C             LDT >= MAX(1,N*K),  if TYPET = 'C'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,N*K) / (LDB,NRHS)
C             On entry, the leading NRHS-by-N*K / N*K-by-NRHS part of
C             this array must contain the right hand side matrix B.
C             On exit, the leading NRHS-by-N*K / N*K-by-NRHS part of
C             this array contains the solution matrix X.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,NRHS),  if TYPET = 'R';
C             LDB >= MAX(1,N*K),   if TYPET = 'C'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -10,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*K*K+(N+2)*K).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction algorithm failed. The Toeplitz matrix
C                   associated with T is not (numerically) positive
C                   definite.
C
C     METHOD
C
C     Householder transformations, modified hyperbolic rotations and
C     block Gaussian eliminations are used in the Schur algorithm [1],
C     [2].
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
C     The implemented method is numerically equivalent with forming
C     the Cholesky factor R and the inverse Cholesky factor of T, using
C     the generalized Schur algorithm, and solving the systems of
C     equations  R*X = L*B  or  X*R = B*L by a blocked backward
C     substitution algorithm.
C                               3 2    2 2
C     The algorithm requires 0(K N  + K N NRHS) floating point
C     operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Chemnitz, Germany, December 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C     February 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, Householder transformation, matrix
C     operations, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TYPET
      INTEGER           INFO, K, LDB, LDT, LDWORK, N, NRHS
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), DWORK(*), T(LDT,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, MAXWRK, STARTH, STARTI, STARTN,
     $                  STARTR, STARTT
      LOGICAL           ISROW
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLACPY, DLASET, DPOTRF, DTRMM, DTRSM,
     $                  MB02CX, MB02CY, XERBLA
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
      ELSE IF ( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF ( LDT.LT.1 .OR. ( ISROW .AND. LDT.LT.K ) .OR.
     $                   ( .NOT.ISROW .AND. LDT.LT.N*K ) ) THEN
         INFO = -6
      ELSE IF ( LDB.LT.1 .OR. ( ISROW .AND. LDB.LT.NRHS ) .OR.
     $                   ( .NOT.ISROW .AND. LDB.LT.N*K ) ) THEN
         INFO = -8
      ELSE IF ( LDWORK.LT.MAX( 1, N*K*K + ( N + 2 )*K ) ) THEN
         DWORK(1) = MAX( 1, N*K*K + ( N + 2 )*K )
         INFO = -10
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02ED', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( K, N, NRHS ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      MAXWRK = 0
      STARTN = 1
      STARTT = N*K*K  + 1
      STARTH = STARTT + 3*K
C
      IF ( ISROW ) THEN
C
C        T is the first block row of a block Toeplitz matrix.
C        Bring T to proper form by triangularizing its first block.
C
         CALL DPOTRF( 'Upper', K, T, LDT, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The matrix is not positive definite.
C
            INFO = 1
            RETURN
         END IF
C
         IF ( N.GT.1 )
     $      CALL DTRSM( 'Left', 'Upper', 'Transpose', 'NonUnit', K,
     $                  (N-1)*K, ONE, T, LDT, T(1,K+1), LDT )
C
C        Initialize the generator, do the first Schur step and set
C        B = -B.
C        T contains the nonzero blocks of the positive parts in the
C        generator and the inverse generator.
C        DWORK(STARTN) contains the nonzero blocks of the negative parts
C        in the generator and the inverse generator.
C
         CALL DTRSM( 'Right', 'Upper', 'NonTranspose', 'NonUnit', NRHS,
     $               K, ONE, T, LDT, B, LDB )
         IF ( N.GT.1 )
     $      CALL DGEMM( 'NonTranspose', 'NonTranspose', NRHS, (N-1)*K,
     $                  K, ONE, B, LDB, T(1,K+1), LDT, -ONE, B(1,K+1),
     $                  LDB )
C
         CALL DLASET( 'All', K, K, ZERO, ONE, DWORK(STARTN), K )
         CALL DTRSM(  'Left', 'Upper', 'Transpose', 'NonUnit', K, K,
     $                ONE, T, LDT, DWORK(STARTN), K )
         IF ( N.GT.1 )
     $      CALL DLACPY( 'All', K, (N-1)*K, T(1,K+1), LDT,
     $                   DWORK(STARTN+K*K), K )
         CALL DLACPY( 'All', K, K, DWORK(STARTN), K, T(1,(N-1)*K+1),
     $                LDT )
C
         CALL DTRMM ( 'Right', 'Lower', 'NonTranspose', 'NonUnit', NRHS,
     $                K, ONE, T(1,(N-1)*K+1), LDT, B, LDB )
C
C        Processing the generator.
C
         DO 10  I = 2, N
            STARTR = ( I - 1 )*K + 1
            STARTI = ( N - I )*K + 1
C
C           Transform the generator of T to proper form.
C
            CALL MB02CX( 'Row', K, K, K, T, LDT,
     $                   DWORK(STARTN+(I-1)*K*K), K, DWORK(STARTT), 3*K,
     $                   DWORK(STARTH), LDWORK-STARTH+1, IERR )
C
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
            CALL MB02CY( 'Row', 'NoStructure', K, K, (N-I)*K, K,
     $                    T(1,K+1), LDT, DWORK(STARTN+I*K*K), K,
     $                    DWORK(STARTN+(I-1)*K*K), K, DWORK(STARTT),
     $                    3*K, DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Block Gaussian eliminates the i-th block in B.
C
            CALL DTRSM( 'Right', 'Upper', 'NonTranspose', 'NonUnit',
     $                  NRHS, K, -ONE, T, LDT, B(1,STARTR), LDB )
            IF ( N.GT.I )
     $         CALL DGEMM( 'NonTranspose', 'NonTranspose', NRHS,
     $                     (N-I)*K, K, ONE, B(1,STARTR), LDB, T(1,K+1),
     $                     LDT, ONE, B(1,STARTR+K), LDB )
C
C           Apply hyperbolic transformations on the negative generator.
C
            CALL DLASET( 'All', K, K, ZERO, ZERO, T(1,STARTI), LDT )
            CALL MB02CY( 'Row', 'NoStructure', K, K, (I-1)*K, K,
     $                   T(1,STARTI), LDT, DWORK(STARTN), K,
     $                   DWORK(STARTN+(I-1)*K*K), K, DWORK(STARTT), 3*K,
     $                   DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Note that  DWORK(STARTN+(I-1)*K*K)  serves simultaneously
C           as the transformation container as well as the new block in
C           the negative generator.
C
            CALL MB02CY( 'Row', 'Triangular', K, K, K, K,
     $                   T(1,(N-1)*K+1), LDT, DWORK(STARTN+(I-1)*K*K),
     $                   K, DWORK(STARTN+(I-1)*K*K), K, DWORK(STARTT),
     $                   3*K, DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Finally the Gaussian elimination is applied on the inverse
C           generator.
C
            CALL DGEMM( 'NonTranspose', 'NonTranspose', NRHS, (I-1)*K,
     $                  K, ONE, B(1,STARTR), LDB, T(1,STARTI), LDT, ONE,
     $                  B, LDB )
            CALL DTRMM( 'Right', 'Lower', 'NonTranspose', 'NonUnit',
     $                  NRHS, K, ONE, T(1,(N-1)*K+1), LDT, B(1,STARTR),
     $                  LDB )
   10    CONTINUE
C
      ELSE
C
C        T is the first block column of a block Toeplitz matrix.
C        Bring T to proper form by triangularizing its first block.
C
         CALL DPOTRF( 'Lower', K, T, LDT, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The matrix is not positive definite.
C
            INFO = 1
            RETURN
         END IF
C
         IF ( N.GT.1 )
     $      CALL DTRSM( 'Right', 'Lower', 'Transpose', 'NonUnit',
     $                  (N-1)*K, K, ONE, T, LDT, T(K+1,1), LDT )
C
C        Initialize the generator, do the first Schur step and set
C        B = -B.
C        T contains the nonzero blocks of the positive parts in the
C        generator and the inverse generator.
C        DWORK(STARTN) contains the nonzero blocks of the negative parts
C        in the generator and the inverse generator.
C
         CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit', K,
     $               NRHS, ONE, T, LDT, B, LDB )
         IF ( N.GT.1 )
     $      CALL DGEMM( 'NonTranspose', 'NonTranspose', (N-1)*K, NRHS,
     $                  K, ONE, T(K+1,1), LDT, B, LDB, -ONE, B(K+1,1),
     $                  LDB )
C
         CALL DLASET( 'All', K, K, ZERO, ONE, DWORK(STARTN), N*K )
         CALL DTRSM(  'Right', 'Lower', 'Transpose', 'NonUnit', K, K,
     $                ONE, T, LDT, DWORK(STARTN), N*K )
         IF ( N.GT.1 )
     $      CALL DLACPY( 'All', (N-1)*K, K, T(K+1,1), LDT,
     $                   DWORK(STARTN+K), N*K )
         CALL DLACPY( 'All', K, K, DWORK(STARTN), N*K, T((N-1)*K+1,1),
     $                LDT )
C
         CALL DTRMM ( 'Left', 'Upper', 'NonTranspose', 'NonUnit', K,
     $                NRHS, ONE, T((N-1)*K+1,1), LDT, B, LDB )
C
C        Processing the generator.
C
         DO 20  I = 2, N
            STARTR = ( I - 1 )*K + 1
            STARTI = ( N - I )*K + 1
C
C           Transform the generator of T to proper form.
C
            CALL MB02CX( 'Column', K, K, K, T, LDT,
     $                   DWORK(STARTN+(I-1)*K), N*K, DWORK(STARTT), 3*K,
     $                   DWORK(STARTH), LDWORK-STARTH+1, IERR )
C
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
            CALL MB02CY( 'Column', 'NoStructure', K, K, (N-I)*K, K,
     $                   T(K+1,1), LDT, DWORK(STARTN+I*K), N*K,
     $                   DWORK(STARTN+(I-1)*K), N*K, DWORK(STARTT),
     $                   3*K, DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Block Gaussian eliminates the i-th block in B.
C
            CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit', K,
     $                  NRHS, -ONE, T, LDT, B(STARTR,1), LDB )
            IF ( N.GT.I )
     $         CALL DGEMM( 'NonTranspose', 'NonTranspose', (N-I)*K,
     $                     NRHS, K, ONE, T(K+1,1), LDT, B(STARTR,1),
     $                     LDB, ONE, B(STARTR+K,1), LDB )
C
C           Apply hyperbolic transformations on the negative generator.
C
            CALL DLASET( 'All', K, K, ZERO, ZERO, T(STARTI,1), LDT )
            CALL MB02CY( 'Column', 'NoStructure', K, K, (I-1)*K, K,
     $                   T(STARTI,1), LDT, DWORK(STARTN), N*K,
     $                   DWORK(STARTN+(I-1)*K), N*K, DWORK(STARTT), 3*K,
     $                   DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Note that  DWORK(STARTN+(I-1)*K)  serves simultaneously
C           as the transformation container as well as the new block in
C           the negative generator.
C
            CALL MB02CY( 'Column', 'Triangular', K, K, K, K,
     $                   T((N-1)*K+1,1), LDT, DWORK(STARTN+(I-1)*K),
     $                   N*K, DWORK(STARTN+(I-1)*K), N*K, DWORK(STARTT),
     $                   3*K, DWORK(STARTH), LDWORK-STARTH+1, IERR )
            MAXWRK = MAX( MAXWRK, INT( DWORK(STARTH) ) )
C
C           Finally the Gaussian elimination is applied on the inverse
C           generator.
C
            CALL DGEMM( 'NonTranspose', 'NonTranspose', (I-1)*K, NRHS,
     $                  K, ONE, T(STARTI,1), LDT, B(STARTR,1), LDB, ONE,
     $                  B, LDB )
            CALL DTRMM( 'Left', 'Upper', 'NonTranspose', 'NonUnit',
     $                  K, NRHS, ONE, T((N-1)*K+1,1), LDT, B(STARTR,1),
     $                  LDB )
C
   20    CONTINUE
C
      END IF
C
      DWORK(1) = MAX( 1, STARTH - 1 + MAXWRK )
C
      RETURN
C
C *** Last line of MB02ED ***
      END
