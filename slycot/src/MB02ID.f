      SUBROUTINE MB02ID( JOB, K, L, M, N, RB, RC, TC, LDTC, TR, LDTR, B,
     $                   LDB, C, LDC, DWORK, LDWORK, INFO )
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
C     To solve the overdetermined or underdetermined real linear systems
C     involving an M*K-by-N*L block Toeplitz matrix T that is specified
C     by its first block column and row. It is assumed that T has full
C     rank.
C     The following options are provided:
C
C     1. If JOB = 'O' or JOB = 'A' :  find the least squares solution of
C        an overdetermined system, i.e., solve the least squares problem
C
C                  minimize || B - T*X ||.                           (1)
C
C     2. If JOB = 'U' or JOB = 'A' :  find the minimum norm solution of
C        the undetermined system
C                   T
C                  T * X = C.                                        (2)
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the problem to be solved as follows
C             = 'O':  solve the overdetermined system (1);
C             = 'U':  solve the underdetermined system (2);
C             = 'A':  solve (1) and (2).
C
C     Input/Output Parameters
C
C     K       (input) INTEGER
C             The number of rows in the blocks of T.  K >= 0.
C
C     L       (input) INTEGER
C             The number of columns in the blocks of T.  L >= 0.
C
C     M       (input) INTEGER
C             The number of blocks in the first block column of T.
C             M >= 0.
C
C     N       (input) INTEGER
C             The number of blocks in the first block row of T.
C             0 <= N <= M*K / L.
C
C     RB      (input) INTEGER
C             If JOB = 'O' or 'A', the number of columns in B.  RB >= 0.
C
C     RC      (input) INTEGER
C             If JOB = 'U' or 'A', the number of columns in C.  RC >= 0.
C
C     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L)
C             On entry, the leading M*K-by-L part of this array must
C             contain the first block column of T.
C
C     LDTC    INTEGER
C             The leading dimension of the array TC.  LDTC >= MAX(1,M*K)
C
C     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L)
C             On entry, the leading K-by-(N-1)*L part of this array must
C             contain the 2nd to the N-th blocks of the first block row
C             of T.
C
C     LDTR    INTEGER
C             The leading dimension of the array TR.  LDTR >= MAX(1,K).
C
C     B       (input/output)  DOUBLE PRECISION array, dimension (LDB,RB)
C             On entry, if JOB = 'O' or JOB = 'A', the leading M*K-by-RB
C             part of this array must contain the right hand side
C             matrix B of the overdetermined system (1).
C             On exit, if JOB = 'O' or JOB = 'A', the leading N*L-by-RB
C             part of this array contains the solution of the
C             overdetermined system (1).
C             This array is not referenced if JOB = 'U'.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,M*K),  if JOB = 'O'  or  JOB = 'A';
C             LDB >= 1,           if JOB = 'U'.
C
C     C       (input)  DOUBLE PRECISION array, dimension (LDC,RC)
C             On entry, if JOB = 'U' or JOB = 'A', the leading N*L-by-RC
C             part of this array must contain the right hand side
C             matrix C of the underdetermined system (2).
C             On exit, if JOB = 'U' or JOB = 'A', the leading M*K-by-RC
C             part of this array contains the solution of the
C             underdetermined system (2).
C             This array is not referenced if JOB = 'O'.
C
C     LDC     INTEGER
C             The leading dimension of the array C.
C             LDB >= 1,           if JOB = 'O';
C             LDB >= MAX(1,M*K),  if JOB = 'U'  or  JOB = 'A'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -17,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             Let x = MAX( 2*N*L*(L+K) + (6+N)*L,(N*L+M*K+1)*L + M*K )
C             and y = N*M*K*L + N*L, then
C             if MIN( M,N ) = 1 and JOB = 'O',
C                         LDWORK >= MAX( y + MAX( M*K,RB ),1 );
C             if MIN( M,N ) = 1 and JOB = 'U',
C                         LDWORK >= MAX( y + MAX( M*K,RC ),1 );
C             if MIN( M,N ) = 1 and JOB = 'A',
C                         LDWORK >= MAX( y +MAX( M*K,MAX( RB,RC ),1 );
C             if MIN( M,N ) > 1 and JOB = 'O',
C                         LDWORK >= MAX( x,N*L*RB + 1 );
C             if MIN( M,N ) > 1 and JOB = 'U',
C                         LDWORK >= MAX( x,N*L*RC + 1 );
C             if MIN( M,N ) > 1 and JOB = 'A',
C                         LDWORK >= MAX( x,N*L*MAX( RB,RC ) + 1 ).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction algorithm failed. The Toeplitz matrix
C                   associated with T is (numerically) not of full rank.
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
C     The algorithm requires O( L*L*K*(N+M)*log(N+M) + N*N*L*L*(L+K) )
C     and additionally
C
C     if JOB = 'O' or JOB = 'A',
C                  O( (K*L+RB*L+K*RB)*(N+M)*log(N+M) + N*N*L*L*RB );
C     if JOB = 'U' or JOB = 'A',
C                  O( (K*L+RC*L+K*RC)*(N+M)*log(N+M) + N*N*L*L*RC );
C
C     floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001.
C     D. Kressner, Technical Univ. Berlin, Germany, July 2002.
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004.
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
      CHARACTER         JOB
      INTEGER           INFO, K, L, LDB, LDC, LDTC, LDTR, LDWORK, M, N,
     $                  RB, RC
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), C(LDC,*), DWORK(LDWORK), TC(LDTC,*),
     $                  TR(LDTR,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, KK, LEN, NB, NBMIN, PDI, PDW, PNI, PNR,
     $                  PPI, PPR, PT, RNK, WRKMIN, WRKOPT, X, Y
      LOGICAL           COMPO, COMPU
C     .. Local Arrays ..
      INTEGER           IPVT(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGELS, DGEMM, DGEQRF, DLACPY, DLASET, DORGQR,
     $                  DTRMM, DTRSM, DTRTRI, MA02AD, MB02CU, MB02CV,
     $                  MB02KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      COMPO = LSAME( JOB, 'O' ) .OR. LSAME( JOB, 'A' )
      COMPU = LSAME( JOB, 'U' ) .OR. LSAME( JOB, 'A' )
      X = MAX( 2*N*L*( L + K ) + ( 6 + N )*L,
     $         ( N*L + M*K + 1 )*L + M*K )
      Y = N*M*K*L + N*L
      IF ( MIN( M, N ).EQ.1 ) THEN
         WRKMIN = MAX( M*K, 1 )
         IF ( COMPO )  WRKMIN = MAX( WRKMIN, RB )
         IF ( COMPU )  WRKMIN = MAX( WRKMIN, RC )
         WRKMIN = MAX( Y + WRKMIN, 1 )
      ELSE
         WRKMIN = X
         IF ( COMPO )  WRKMIN = MAX( WRKMIN, N*L*RB + 1 )
         IF ( COMPU )  WRKMIN = MAX( WRKMIN, N*L*RC + 1 )
      END IF
      WRKOPT = 1
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( COMPO .OR. COMPU ) ) THEN
         INFO = -1
      ELSE IF ( K.LT.0 ) THEN
         INFO = -2
      ELSE IF ( L.LT.0 ) THEN
         INFO = -3
      ELSE IF ( M.LT.0 ) THEN
         INFO = -4
      ELSE IF ( N.LT.0 .OR. ( N*L ).GT.( M*K ) ) THEN
         INFO = -5
      ELSE IF ( COMPO .AND. RB.LT.0 ) THEN
         INFO = -6
      ELSE IF ( COMPU .AND. RC.LT.0  ) THEN
         INFO = -7
      ELSE IF ( LDTC.LT.MAX( 1, M*K ) ) THEN
         INFO = -9
      ELSE IF ( LDTR.LT.MAX( 1, K ) ) THEN
         INFO = -11
      ELSE IF ( LDB.LT.1 .OR. ( COMPO .AND. LDB.LT.M*K ) ) THEN
         INFO = -13
      ELSE IF ( LDC.LT.1 .OR. ( COMPU .AND. LDC.LT.M*K ) ) THEN
         INFO = -15
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -17
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02ID', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( COMPO .AND. MIN( N*L, RB ).EQ.0 ) THEN
         COMPO = .FALSE.
      END IF
      IF( COMPU .AND. MIN( N*L, RC ).EQ.0 ) THEN
         CALL DLASET( 'Full', M*K, RC, ZERO, ZERO, C, LDC )
         COMPU = .FALSE.
      END IF
      IF ( .NOT.( COMPO .OR. COMPU ) ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Check cases M = 1 or N = 1.
C
      IF ( MIN( M, N ).EQ.1 ) THEN
         PDW = K*L*M*N
         IF ( COMPO ) THEN
            CALL DLACPY( 'All', M*K, L, TC, LDTC, DWORK, M*K )
            CALL DLACPY( 'All', K, (N-1)*L, TR, LDTR, DWORK(K*L+1),
     $                   M*K )
            CALL DGELS( 'NonTranspose', M*K, N*L, RB, DWORK, M*K, B,
     $                  LDB, DWORK(PDW+1), LDWORK-PDW, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+1) ) + PDW )
         END IF
         IF ( COMPU ) THEN
            CALL DLACPY( 'All', M*K, L, TC, LDTC, DWORK, M*K )
            CALL DLACPY( 'All', K, (N-1)*L, TR, LDTR, DWORK(K*L+1),
     $                   M*K )
            CALL DGELS( 'Transpose', M*K, N*L, RC, DWORK, M*K, C, LDC,
     $                  DWORK(PDW+1), LDWORK-PDW, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+1) ) + PDW )
         END IF
         DWORK(1) = DBLE( WRKOPT )
         RETURN
      END IF
C
C     Step 1:  Compute the generator.
C
      IF ( COMPO ) THEN
         CALL MB02KD( 'Column', 'Transpose', K, L, M, N, RB, ONE, ZERO,
     $                TC, LDTC, TR, LDTR, B, LDB, DWORK, N*L,
     $                DWORK(N*L*RB+1), LDWORK-N*L*RB, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(N*L*RB+1) ) + N*L*RB )
         CALL DLACPY( 'All', N*L, RB, DWORK, N*L, B, LDB )
      END IF
C
      PDW = N*L*L + 1
      CALL DLACPY( 'All', M*K, L, TC, LDTC, DWORK(PDW), M*K )
      CALL DGEQRF( M*K, L, DWORK(PDW), M*K, DWORK(PDW+M*K*L),
     $             DWORK(PDW+(M*K+1)*L), LDWORK-PDW-(M*K+1)*L-1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+(M*K+1)*L) ) +
     $                      PDW + (M*K+1)*L - 1 )
C
      DO 10  I = PDW, PDW + M*K*L - 1, M*K + 1
         IF ( DWORK(I).EQ.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
   10 CONTINUE
C
      CALL MA02AD( 'Upper', L, L, DWORK(PDW), M*K, DWORK, N*L )
      CALL DORGQR( M*K, L, L, DWORK(PDW), M*K, DWORK(PDW+M*K*L),
     $             DWORK(PDW+(M*K+1)*L), LDWORK-PDW-(M*K+1)*L-1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+(M*K+1)*L) ) +
     $                      PDW + (M*K+1)*L - 1 )
      CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, L, ONE, ZERO,
     $              TC, LDTC, TR, LDTR, DWORK(PDW), M*K, DWORK(L+1),
     C              N*L, DWORK(PDW+M*K*L), LDWORK-PDW-M*K*L+1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+M*K*L) ) + PDW + M*K*L - 1 )
      PPR = N*L*L + 1
      PNR = N*L*( L + K ) + 1
      CALL MA02AD( 'All', K, (N-1)*L, TR, LDTR, DWORK(PPR+L), N*L )
      CALL DLACPY( 'All', (N-1)*L, L, DWORK(L+1), N*L, DWORK(PNR+L),
     $             N*L )
      PT  = ( M - 1 )*K + 1
      PDW = PNR + N*L*L + L
C
      DO 30  I = 1, MIN( M, N-1 )
         CALL MA02AD( 'All', K, L, TC(PT,1), LDTC, DWORK(PDW), N*L )
         PT  = PT  - K
         PDW = PDW + L
   30 CONTINUE
C
      PT = 1
C
      DO 40  I = M + 1, N - 1
         CALL MA02AD( 'All', K, L, TR(1,PT), LDTR, DWORK(PDW), N*L )
         PT  = PT  + L
         PDW = PDW + L
   40 CONTINUE
C
      IF ( COMPO ) THEN
C
C        Apply the first reduction step to T'*B.
C
         CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit',
     $               L, RB, ONE, DWORK, N*L, B, LDB )
         CALL DGEMM( 'NoTranspose', 'NoTranspose', (N-1)*L, RB, L, ONE,
     $               DWORK(L+1), N*L, B, LDB, -ONE, B(L+1,1), LDB )
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'NonUnit', L,
     $               RB, ONE, DWORK, N*L, B, LDB )
      END IF
C
      IF ( COMPU ) THEN
C
C        Apply the first reduction step to C.
C
         CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit',
     $               L, RC, ONE, DWORK, N*L, C, LDC )
         CALL DGEMM( 'NoTranspose', 'NoTranspose', (N-1)*L, RC, L, ONE,
     $               DWORK(L+1), N*L, C, LDC, -ONE, C(L+1,1), LDC )
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'NonUnit', L,
     $               RC, ONE, DWORK, N*L, C, LDC )
      END IF
C
      PDI = ( N - 1 )*L + 1
      CALL DLACPY( 'Lower', L, L, DWORK, N*L, DWORK(PDI), N*L )
      CALL DTRTRI( 'Lower', 'NonUnit', L, DWORK(PDI), N*L, IERR )
      CALL MA02AD( 'Lower', L-1, L, DWORK(PDI+1), N*L,
     $             DWORK((2*N-1)*L+1), N*L )
      CALL DLASET( 'Lower', L-1, L, ZERO, ZERO, DWORK(PDI+1), N*L )
      CALL DLACPY( 'Upper', L, L, DWORK(PDI), N*L, DWORK(PNR), N*L )
      CALL DLASET( 'Lower', L-1, L, ZERO, ZERO, DWORK(PNR+1), N*L )
      CALL DLASET( 'All', L, K, ZERO, ZERO, DWORK(PPR), N*L )
      CALL DLASET( 'All', L, K, ZERO, ZERO, DWORK(PNR+N*L*L), N*L )
C
      PPI = PPR
      PPR = PPR + L
      PNI = PNR
      PNR = PNR + L
      PDW = 2*N*L*( L + K ) + 1
      LEN = ( N - 1 )*L
C
C     Determine block size for the involved block Householder
C     transformations.
C
      NB = MIN( ILAENV( 1, 'DGELQF', ' ', N*L, L, -1, -1 ), L )
      KK = PDW + 6*L - 1
      WRKOPT = MAX( WRKOPT, KK + N*L*NB )
      KK = LDWORK - KK
      IF ( KK.LT.N*L*NB )  NB = KK / ( N*L )
      NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', N*L, L, -1, -1 ) )
      IF ( NB.LT.NBMIN )  NB = 0
C
      DO 50  I = L + 1, N*L, L
         CALL MB02CU( 'Column', L, L+K, L+K, NB, DWORK, N*L, DWORK(PPR),
     $                N*L, DWORK(PNR), N*L, RNK, IPVT, DWORK(PDW), ZERO,
     $                DWORK(PDW+6*L), LDWORK-PDW-6*L+1, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The rank condition is (numerically) not
C                          satisfied.
C
            INFO = 1
            RETURN
         END IF
         CALL MB02CV( 'Column', 'NoStructure', L, LEN-L, L+K, L+K, NB,
     $                -1, DWORK, N*L, DWORK(PPR), N*L, DWORK(PNR), N*L,
     $                DWORK(L+1), N*L, DWORK(PPR+L), N*L, DWORK(PNR+L),
     $                N*L, DWORK(PDW), DWORK(PDW+6*L), LDWORK-PDW-6*L+1,
     $                IERR )
         PDI = PDI - L
         IF ( COMPO ) THEN
C
C           Block Gaussian elimination to B.
C
            CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit',
     $                  L, RB, -ONE, DWORK, N*L, B(I,1), LDB )
            IF ( LEN.GT.L ) THEN
               CALL DGEMM( 'NonTranspose', 'NonTranspose', LEN-L, RB, L,
     $                     ONE, DWORK(L+1), N*L, B(I,1), LDB, ONE,
     $                     B(I+L,1), LDB )
            END IF
         END IF
         IF ( COMPU ) THEN
C
C           Block Gaussian elimination to C.
C
            CALL DTRSM( 'Left', 'Lower', 'NonTranspose', 'NonUnit',
     $                  L, RC, -ONE, DWORK, N*L, C(I,1), LDC )
            IF ( LEN.GT.L ) THEN
               CALL DGEMM( 'NonTranspose', 'NonTranspose', LEN-L, RC, L,
     $                     ONE, DWORK(L+1), N*L, C(I,1), LDC, ONE,
     $                     C(I+L,1), LDC )
            END IF
         END IF
         CALL DLASET( 'All', L, L, ZERO, ZERO, DWORK(PDI), N*L )
         CALL MB02CV( 'Column', 'Triangular', L, I+L-1, L+K, L+K, NB,
     $                -1, DWORK, N*L, DWORK(PPR), N*L, DWORK(PNR), N*L,
     $                DWORK(PDI), N*L, DWORK(PPI), N*L, DWORK(PNI), N*L,
     $                DWORK(PDW), DWORK(PDW+6*L), LDWORK-PDW-6*L+1,
     $                IERR )
         IF ( COMPO ) THEN
C
C           Apply block Gaussian elimination to B.
C
            CALL DGEMM( 'NoTranspose', 'NoTranspose', I-1, RB, L, ONE,
     $                  DWORK(PDI), N*L, B(I,1), LDB, ONE, B, LDB )
            CALL DTRMM( 'Left', 'Upper', 'NonTranspose', 'NonUnit', L,
     $                  RB, ONE, DWORK((N-1)*L+1), N*L, B(I,1), LDB )
         END IF
         IF ( COMPU ) THEN
C
C           Apply block Gaussian elimination to C.
C
            CALL DGEMM( 'NonTranspose', 'NonTranspose', I-1, RC, L, ONE,
     $                  DWORK(PDI), N*L, C(I,1), LDC, ONE, C, LDC )
            CALL DTRMM( 'Left', 'Upper', 'NonTranspose', 'NonUnit', L,
     $                  RC, ONE, DWORK((N-1)*L+1), N*L, C(I,1), LDC )
         END IF
         LEN = LEN - L
         PNR = PNR + L
         PPR = PPR + L
   50 CONTINUE
C
      IF ( COMPU ) THEN
         CALL MB02KD( 'Column', 'NonTranspose', K, L, M, N, RC, ONE,
     $                ZERO, TC, LDTC, TR, LDTR, C, LDC, DWORK, M*K,
     $                DWORK(M*K*RC+1), LDWORK-M*K*RC, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(M*K*RC+1) ) + M*K*RC )
         CALL DLACPY( 'All', M*K, RC, DWORK, M*K, C, LDC )
      END IF
      DWORK(1) = DBLE( WRKOPT )
      RETURN
C
C *** Last line of MB02ID ***
      END
