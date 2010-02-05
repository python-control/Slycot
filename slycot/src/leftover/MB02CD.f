      SUBROUTINE MB02CD( JOB, TYPET, K, N, T, LDT, G, LDG, R, LDR, L,
     $                   LDL, CS, LCS, DWORK, LDWORK, INFO )
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
C     To compute the Cholesky factor and the generator and/or the
C     Cholesky factor of the inverse of a symmetric positive definite
C     (s.p.d.) block Toeplitz matrix T, defined by either its first
C     block row, or its first block column, depending on the routine
C     parameter TYPET. Transformation information is stored.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the output of the routine, as follows:
C             = 'G':  only computes the generator G of the inverse;
C             = 'R':  computes the generator G of the inverse and the
C                     Cholesky factor R of T, i.e., if TYPET = 'R',
C                     then R'*R = T, and if TYPET = 'C', then R*R' = T;
C             = 'L':  computes the generator G and the Cholesky factor L
C                     of the inverse, i.e., if TYPET = 'R', then
C                     L'*L = inv(T), and if TYPET = 'C', then
C                     L*L' = inv(T);
C             = 'A':  computes the generator G, the Cholesky factor L
C                     of the inverse and the Cholesky factor R of T;
C             = 'O':  only computes the Cholesky factor R of T.
C
C     TYPET   CHARACTER*1
C             Specifies the type of T, as follows:
C             = 'R':  T contains the first block row of an s.p.d. block
C                     Toeplitz matrix; if demanded, the Cholesky factors
C                     R and L are upper and lower triangular,
C                     respectively, and G contains the transposed
C                     generator of the inverse;
C             = 'C':  T contains the first block column of an s.p.d.
C                     block Toeplitz matrix; if demanded, the Cholesky
C                     factors R and L are lower and upper triangular,
C                     respectively, and G contains the generator of the
C                     inverse. This choice results in a column oriented
C                     algorithm which is usually faster.
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
C     T       (input/output)  DOUBLE PRECISION array, dimension
C             (LDT,N*K) / (LDT,K)
C             On entry, the leading K-by-N*K / N*K-by-K part of this
C             array must contain the first block row / column of an
C             s.p.d. block Toeplitz matrix.
C             On exit, if INFO = 0, then the leading K-by-N*K / N*K-by-K
C             part of this array contains, in the first K-by-K block,
C             the upper / lower Cholesky factor of T(1:K,1:K), and in
C             the remaining part, the Householder transformations
C             applied during the process.
C
C     LDT     INTEGER
C             The leading dimension of the array T.
C             LDT >= MAX(1,K),    if TYPET = 'R';
C             LDT >= MAX(1,N*K),  if TYPET = 'C'.
C
C     G       (output)  DOUBLE PRECISION array, dimension
C             (LDG,N*K) / (LDG,2*K)
C             If INFO = 0 and JOB = 'G', 'R', 'L', or 'A', the leading
C             2*K-by-N*K / N*K-by-2*K part of this array contains, in
C             the first K-by-K block of the second block row / column,
C             the lower right block of L (necessary for updating
C             factorizations in SLICOT Library routine MB02DD), and
C             in the remaining part, the generator of the inverse of T.
C             Actually, to obtain a generator one has to set
C                 G(K+1:2*K, 1:K) = 0,    if TYPET = 'R';
C                 G(1:K, K+1:2*K) = 0,    if TYPET = 'C'.
C
C     LDG     INTEGER
C             The leading dimension of the array G.
C             LDG >= MAX(1,2*K),  if TYPET = 'R' and
C                                    JOB = 'G', 'R', 'L', or 'A';
C             LDG >= MAX(1,N*K),  if TYPET = 'C' and
C                                    JOB = 'G', 'R', 'L', or 'A';
C             LDG >= 1,           if JOB = 'O'.
C
C     R       (output)  DOUBLE PRECISION array, dimension (LDR,N*K)
C             If INFO = 0 and JOB = 'R', 'A', or 'O', then the leading
C             N*K-by-N*K part of this array contains the upper / lower
C             Cholesky factor of T.
C             The elements in the strictly lower / upper triangular part
C             are not referenced.
C
C     LDR     INTEGER
C             The leading dimension of the array R.
C             LDR >= MAX(1,N*K),  if JOB = 'R', 'A', or 'O';
C             LDR >= 1,           if JOB = 'G', or 'L'.
C
C     L       (output)  DOUBLE PRECISION array, dimension (LDL,N*K)
C             If INFO = 0 and JOB = 'L', or 'A', then the leading
C             N*K-by-N*K part of this array contains the lower / upper
C             Cholesky factor of the inverse of T.
C             The elements in the strictly upper / lower triangular part
C             are not referenced.
C
C     LDL     INTEGER
C             The leading dimension of the array L.
C             LDL >= MAX(1,N*K),  if JOB = 'L', or 'A';
C             LDL >= 1,           if JOB = 'G', 'R', or 'O'.
C
C     CS      (output)  DOUBLE PRECISION array, dimension (LCS)
C             If INFO = 0, then the leading 3*(N-1)*K part of this
C             array contains information about the hyperbolic rotations
C             and Householder transformations applied during the
C             process. This information is needed for updating the
C             factorizations in SLICOT Library routine MB02DD.
C
C     LCS     INTEGER
C             The length of the array CS.  LCS >= 3*(N-1)*K.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -16,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,(N-1)*K).
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
C                               3 2
C     The algorithm requires 0(K N ) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2000,
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
      CHARACTER         JOB, TYPET
      INTEGER           INFO, K, LCS, LDG, LDL, LDR, LDT, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), G(LDG, *), L(LDL,*), R(LDR,*),
     $                  T(LDT,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, MAXWRK, STARTI, STARTR, STARTT
      LOGICAL           COMPG, COMPL, COMPR, ISROW
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DLASET, DPOTRF, DTRSM, MB02CX, MB02CY,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      COMPL = LSAME( JOB, 'L' ) .OR. LSAME( JOB, 'A' )
      COMPG = LSAME( JOB, 'G' ) .OR. LSAME( JOB, 'R' ) .OR. COMPL
      COMPR = LSAME( JOB, 'R' ) .OR. LSAME( JOB, 'A' ) .OR.
     $        LSAME( JOB, 'O' )
      ISROW = LSAME( TYPET, 'R' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( COMPG .OR. COMPR ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( ISROW .OR. LSAME( TYPET, 'C' ) ) ) THEN
         INFO = -2
      ELSE IF ( K.LT.0 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( LDT.LT.1 .OR. ( ISROW .AND. LDT.LT.K ) .OR.
     $                   ( .NOT.ISROW .AND. LDT.LT.N*K ) ) THEN
         INFO = -6
      ELSE IF ( LDG.LT.1 .OR.
     $           ( COMPG .AND. ( ( ISROW .AND. LDG.LT.2*K )
     $                 .OR. ( .NOT.ISROW .AND. LDG.LT.N*K ) ) ) ) THEN
         INFO = -8
      ELSE IF ( LDR.LT.1 .OR. ( COMPR .AND. ( LDR.LT.N*K ) ) ) THEN
         INFO = -10
      ELSE IF ( LDL.LT.1 .OR. ( COMPL .AND. ( LDL.LT.N*K ) ) ) THEN
         INFO = -12
      ELSE IF ( LCS.LT.3*( N - 1 )*K ) THEN
         INFO = -14
      ELSE IF ( LDWORK.LT.MAX( 1, ( N - 1 )*K ) ) THEN
         DWORK(1) = MAX( 1, ( N - 1 )*K )
         INFO = -16
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02CD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( K, N ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      MAXWRK = 1
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
C        Initialize the output matrices.
C
         IF ( COMPG ) THEN
            CALL DLASET( 'All', 2*K, N*K, ZERO, ZERO, G, LDG )
            CALL DLASET( 'All', 1, K, ONE, ONE, G(K+1,1), LDG+1 )
            CALL DTRSM(  'Left', 'Upper', 'Transpose', 'NonUnit', K, K,
     $                   ONE, T, LDT, G(K+1,1), LDG )
            IF ( N.GT.1 )
     $         CALL DLACPY( 'Upper', K, (N-1)*K, T, LDT, G(K+1,K+1),
     $                      LDG )
            CALL DLACPY( 'Lower', K, K, G(K+1,1), LDG, G, LDG )
         END IF
C
         IF ( COMPL ) THEN
            CALL DLACPY( 'Lower', K, K, G(K+1,1), LDG, L, LDL )
         END IF
C
         IF ( COMPR ) THEN
            CALL DLACPY( 'Upper', K, N*K, T, LDT, R, LDR )
         END IF
C
C        Processing the generator.
C
         IF ( COMPG ) THEN
C
C           Here we use G as working array for holding the generator.
C           T contains the second row of the generator.
C           G contains in its first block row the second row of the
C           inverse generator.
C           The second block row of G is partitioned as follows:
C
C           [ First block of the inverse generator, ...
C             First row of the generator, ...
C             The rest of the blocks of the inverse generator ]
C
C           The reason for the odd partitioning is that the first block
C           of the inverse generator will be thrown out at the end and
C           we want to avoid reordering.
C
C           (N-1)*K locations of DWORK are used by SLICOT Library
C           routine MB02CY.
C
            DO 10  I = 2, N
               STARTR = ( I - 1 )*K + 1
               STARTI = ( N - I + 1 )*K + 1
               STARTT = 3*( I - 2 )*K + 1
C
C              Transformations acting on the generator:
C
               CALL MB02CX( 'Row', K, K, K, G(K+1,K+1), LDG,
     $                      T(1,STARTR), LDT, CS(STARTT), 3*K, DWORK,
     $                      LDWORK, IERR )
C
               IF ( IERR.NE.0 )  THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               IF ( N.GT.I )  THEN
                  CALL MB02CY( 'Row', 'NoStructure', K, K, (N-I)*K, K,
     $                         G(K+1,2*K+1), LDG, T(1,STARTR+K), LDT,
     $                         T(1,STARTR), LDT, CS(STARTT), 3*K, DWORK,
     $                         LDWORK, IERR )
                  MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               END IF
C
               IF ( COMPR ) THEN
                  CALL DLACPY( 'Upper', K, (N-I+1)*K, G(K+1,K+1), LDG,
     $                         R(STARTR,STARTR), LDR)
               END IF
C
C              Transformations acting on the inverse generator:
C
               CALL DLASET( 'All', K, K, ZERO, ZERO, G(K+1,STARTI),
     $                      LDG )
               CALL MB02CY( 'Row', 'Triangular', K, K, K, K, G(K+1,1),
     $                      LDG, G(1,STARTR), LDG, T(1,STARTR), LDT,
     $                      CS(STARTT), 3*K, DWORK, LDWORK, IERR )
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
               CALL MB02CY( 'Row', 'NoStructure', K, K, (I-1)*K, K,
     $                      G(K+1,STARTI), LDG, G, LDG, T(1,STARTR),
     $                      LDT, CS(STARTT), 3*K, DWORK, LDWORK, IERR )
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
               IF ( COMPL ) THEN
                  CALL DLACPY( 'All', K, (I-1)*K, G(K+1,STARTI), LDG,
     $                         L(STARTR,1), LDL )
                  CALL DLACPY( 'Lower', K, K, G(K+1,1), LDG,
     $                         L(STARTR,(I-1)*K+1), LDL )
               END IF
   10       CONTINUE
C
         ELSE
C
C           Here R is used as working array for holding the generator.
C           Again, T contains the second row of the generator.
C           The current row of R contains the first row of the
C           generator.
C
            IF ( N.GT.1 )
     $         CALL DLACPY( 'Upper', K, (N-1)*K, T, LDT, R(K+1,K+1),
     $                      LDR )
C
            DO 20  I = 2, N
               STARTR = ( I - 1 )*K + 1
               STARTT = 3*( I - 2 )*K + 1
               CALL MB02CX( 'Row', K, K, K, R(STARTR,STARTR), LDR,
     $                      T(1,STARTR), LDT, CS(STARTT), 3*K, DWORK,
     $                      LDWORK, IERR )
               IF ( IERR.NE.0 )  THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               IF ( N.GT.I ) THEN
                  CALL MB02CY( 'Row', 'NoStructure', K, K, (N-I)*K, K,
     $                         R(STARTR,STARTR+K), LDR, T(1,STARTR+K),
     $                         LDT, T(1,STARTR), LDT, CS(STARTT), 3*K,
     $                         DWORK, LDWORK, IERR )
                  MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
                  CALL DLACPY( 'Upper', K, (N-I)*K, R(STARTR,STARTR),
     $                          LDR, R(STARTR+K,STARTR+K), LDR )
               END IF
   20       CONTINUE
C
         END IF
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
C        Initialize the output matrices.
C
         IF ( COMPG ) THEN
            CALL DLASET( 'All', N*K, 2*K, ZERO, ZERO, G, LDG )
            CALL DLASET( 'All', 1, K, ONE, ONE, G(1,K+1), LDG+1 )
            CALL DTRSM(  'Right', 'Lower', 'Transpose', 'NonUnit', K, K,
     $                   ONE, T, LDT, G(1,K+1), LDG )
            IF ( N.GT.1 )
     $         CALL DLACPY( 'Lower', (N-1)*K, K, T, LDT, G(K+1,K+1),
     $                      LDG )
            CALL DLACPY( 'Upper', K, K, G(1,K+1), LDG, G, LDG )
         END IF
C
         IF ( COMPL ) THEN
            CALL DLACPY( 'Upper', K, K, G(1,K+1), LDG, L, LDL )
         END IF
C
         IF ( COMPR ) THEN
            CALL DLACPY( 'Lower', N*K, K, T, LDT, R, LDR )
         END IF
C
C        Processing the generator.
C
         IF ( COMPG ) THEN
C
C           Here we use G as working array for holding the generator.
C           T contains the second column of the generator.
C           G contains in its first block column the second column of
C           the inverse generator.
C           The second block column of G is partitioned as follows:
C
C           [ First block of the inverse generator; ...
C             First column of the generator; ...
C             The rest of the blocks of the inverse generator ]
C
C           The reason for the odd partitioning is that the first block
C           of the inverse generator will be thrown out at the end and
C           we want to avoid reordering.
C
C           (N-1)*K locations of DWORK are used by SLICOT Library
C           routine MB02CY.
C
            DO 30  I = 2, N	
               STARTR = ( I - 1 )*K + 1
               STARTI = ( N - I + 1 )*K + 1
               STARTT = 3*( I - 2 )*K + 1
C
C              Transformations acting on the generator:
C
               CALL MB02CX( 'Column', K, K, K, G(K+1,K+1), LDG,
     $                      T(STARTR,1), LDT, CS(STARTT), 3*K, DWORK,
     $                      LDWORK, IERR )
C
               IF ( IERR.NE.0 )  THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               IF ( N.GT.I ) THEN
                  CALL MB02CY( 'Column', 'NoStructure', K, K, (N-I)*K,
     $                         K, G(2*K+1,K+1), LDG, T(STARTR+K,1), LDT,
     $                         T(STARTR,1), LDT, CS(STARTT), 3*K, DWORK,
     $                         LDWORK, IERR )
                  MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               END IF
C
               IF ( COMPR ) THEN
                  CALL DLACPY( 'Lower', (N-I+1)*K, K, G(K+1,K+1), LDG,
     $                         R(STARTR,STARTR), LDR)
               END IF
C
C              Transformations acting on the inverse generator:
C
               CALL DLASET( 'All', K, K, ZERO, ZERO, G(STARTI,K+1),
     $                      LDG )
               CALL MB02CY( 'Column', 'Triangular', K, K, K, K,
     $                      G(1,K+1), LDG, G(STARTR,1), LDG,
     $                      T(STARTR,1), LDT, CS(STARTT), 3*K, DWORK,
     $                      LDWORK, IERR )
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
               CALL MB02CY( 'Column', 'NoStructure', K, K, (I-1)*K, K,
     $                      G(STARTI,K+1), LDG, G, LDG, T(STARTR,1),
     $                      LDT, CS(STARTT), 3*K, DWORK, LDWORK, IERR )
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
               IF ( COMPL ) THEN
                  CALL DLACPY( 'All', (I-1)*K, K, G(STARTI,K+1), LDG,
     $                         L(1,STARTR), LDL )
                  CALL DLACPY( 'Upper', K, K, G(1,K+1), LDG,
     $                         L((I-1)*K+1,STARTR), LDL )
               END IF
   30       CONTINUE
C
         ELSE
C
C           Here R is used as working array for holding the generator.
C           Again, T contains the second column of the generator.
C           The current column of R contains the first column of the
C           generator.
C
            IF ( N.GT.1 )
     $         CALL DLACPY( 'Lower', (N-1)*K, K, T, LDT, R(K+1,K+1),
     $                      LDR )
C
            DO 40  I = 2, N
               STARTR = ( I - 1 )*K + 1
               STARTT = 3*( I - 2 )*K + 1
               CALL MB02CX( 'Column', K, K, K, R(STARTR,STARTR), LDR,
     $                      T(STARTR,1), LDT, CS(STARTT), 3*K, DWORK,
     $                      LDWORK, IERR )
               IF ( IERR.NE.0 )  THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
               MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
               IF ( N.GT.I ) THEN
                  CALL MB02CY( 'Column', 'NoStructure', K, K, (N-I)*K,
     $                         K, R(STARTR+K,STARTR), LDR,
     $                         T(STARTR+K,1), LDT, T(STARTR,1), LDT,
     $                         CS(STARTT), 3*K, DWORK, LDWORK, IERR )
                  MAXWRK = MAX( MAXWRK, INT( DWORK(1) ) )
C
                  CALL DLACPY( 'Lower', (N-I)*K, K, R(STARTR,STARTR),
     $                         LDR, R(STARTR+K,STARTR+K), LDR )
               END IF
   40       CONTINUE
C
         END IF
      END IF
C
      DWORK(1) = MAXWRK
C
      RETURN
C
C *** Last line of MB02CD ***
      END
