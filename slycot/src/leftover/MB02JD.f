      SUBROUTINE MB02JD( JOB, K, L, M, N, P, S, TC, LDTC, TR, LDTR, Q,
     $                   LDQ, R, LDR, DWORK, LDWORK, INFO )
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
C     To compute a lower triangular matrix R and a matrix Q with
C     Q^T Q = I such that
C                                    T
C                           T  =  Q R ,
C
C     where T is a K*M-by-L*N block Toeplitz matrix with blocks of size
C     (K,L). The first column of T will be denoted by TC and the first
C     row by TR. It is assumed that the first MIN(M*K, N*L) columns of T
C     have full rank.
C
C     By subsequent calls of this routine the factors Q and R can be
C     computed block column by block column.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the output of the routine as follows:
C             = 'Q':  computes Q and R;
C             = 'R':  only computes R.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of rows in one block of T.  K >= 0.
C
C     L       (input)  INTEGER
C             The number of columns in one block of T.  L >= 0.
C
C     M       (input)  INTEGER
C             The number of blocks in one block column of T.  M >= 0.
C
C     N       (input)  INTEGER
C             The number of blocks in one block row of T.  N >= 0.
C
C     P       (input)  INTEGER
C             The number of previously computed block columns of R.
C             P*L < MIN( M*K,N*L ) + L and P >= 0.
C
C     S       (input)  INTEGER
C             The number of block columns of R to compute.
C             (P+S)*L < MIN( M*K,N*L ) + L and S >= 0.
C
C     TC      (input) DOUBLE PRECISION array, dimension (LDTC, L)
C             On entry, if P = 0, the leading M*K-by-L part of this
C             array must contain the first block column of T.
C
C     LDTC    INTEGER
C             The leading dimension of the array TC.
C             LDTC >= MAX(1,M*K).
C
C     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L)
C             On entry, if P = 0, the leading K-by-(N-1)*L part of this
C             array must contain the first block row of T without the
C             leading K-by-L block.
C
C     LDTR    INTEGER
C             The leading dimension of the array TR.
C             LDTR >= MAX(1,K).
C
C     Q       (input/output)  DOUBLE PRECISION array, dimension
C                             (LDQ,MIN( S*L, MIN( M*K,N*L )-P*L ))
C             On entry, if JOB = 'Q'  and  P > 0, the leading M*K-by-L
C             part of this array must contain the last block column of Q
C             from a previous call of this routine.
C             On exit, if JOB = 'Q'  and  INFO = 0, the leading
C             M*K-by-MIN( S*L, MIN( M*K,N*L )-P*L ) part of this array
C             contains the P-th to (P+S)-th block columns of the factor
C             Q.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.
C             LDQ >= MAX(1,M*K), if JOB = 'Q';
C             LDQ >= 1,          if JOB = 'R'.
C
C     R       (input/output)  DOUBLE PRECISION array, dimension
C                             (LDR,MIN( S*L, MIN( M*K,N*L )-P*L ))
C             On entry, if P > 0, the leading (N-P+1)*L-by-L
C             part of this array must contain the nozero part of the
C             last block column of R from a previous call of this
C             routine.
C             One exit, if INFO = 0, the leading
C             MIN( N, N-P+1 )*L-by-MIN( S*L, MIN( M*K,N*L )-P*L )
C             part of this array contains the nonzero parts of the P-th
C             to (P+S)-th block columns of the lower triangular
C             factor R.
C             Note that elements in the strictly upper triangular part
C             will not be referenced.
C
C     LDR     INTEGER
C             The leading dimension of the array R.
C             LDR >= MAX( 1, MIN( N, N-P+1 )*L )
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C             On exit, if INFO = -17,  DWORK(1) returns the minimum
C             value of LDWORK.
C             If JOB = 'Q', the first 1 + ( (N-1)*L + M*K )*( 2*K + L )
C             elements of DWORK should be preserved during successive
C             calls of the routine.
C             If JOB = 'R', the first 1 + (N-1)*L*( 2*K + L ) elements
C             of DWORK should be preserved during successive calls of
C             the routine.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             JOB = 'Q':
C                LDWORK >= 1 + ( M*K + ( N - 1 )*L )*( L + 2*K ) + 6*L
C                            + MAX( M*K,( N - MAX( 1,P )*L ) );
C             JOB = 'R':
C                If P = 0,
C                   LDWORK >= MAX( 1 + ( N - 1 )*L*( L + 2*K ) + 6*L
C                                    + (N-1)*L, M*K*( L + 1 ) + L );
C                If P > 0,
C                   LDWORK >= 1 + (N-1)*L*( L + 2*K ) + 6*L + (N-P)*L.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the full rank condition for the first MIN(M*K, N*L)
C                   columns of T is (numerically) violated.
C
C     METHOD
C
C     Block Householder transformations and modified hyperbolic
C     rotations are used in the Schur algorithm [1], [2].
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
C     The implemented method yields a factor R which has comparable
C     accuracy with the Cholesky factor of T^T * T. Q is implicitly
C     computed from the formula Q = T * inv(R^T R) * R, i.e., for ill
C     conditioned problems this factor is of very limited value.
C                                 2
C     The algorithm requires 0(K*L *M*N) floating point operations.
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
      INTEGER           INFO, K, L, LDQ, LDR, LDTC, LDTR, LDWORK,
     $                  M, N, P, S
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), Q(LDQ,*), R(LDR,*), TC(LDTC,*),
     $                  TR(LDTR,*)
C     .. Local Scalars ..
      INTEGER           COLR, I, IERR, KK, LEN, NB, NBMIN, PDQ, PDW,
     $                  PNQ, PNR, PRE, PT, RNK, SHFR, STPS, WRKMIN,
     $                  WRKOPT
      LOGICAL           COMPQ
C     .. Local Arrays ..
      INTEGER           IPVT(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEQRF, DLACPY, DLASET, DORGQR, MA02AD, MB02CU,
     $                  MB02CV, MB02KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      COMPQ = LSAME( JOB, 'Q' )
      IF ( COMPQ ) THEN
         WRKMIN = 1 + ( M*K + ( N - 1 )*L )*( L + 2*K ) + 6*L
     $              + MAX( M*K, ( N - MAX( 1, P ) )*L )
      ELSE
         WRKMIN = 1 + ( N - 1 )*L*( L + 2*K ) + 6*L
     $              + ( N - MAX( P, 1 ) )*L
         IF ( P.EQ.0 ) THEN
            WRKMIN = MAX( WRKMIN, M*K*( L + 1 ) + L )
         END IF
      END IF
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( COMPQ .OR. LSAME( JOB, 'R' ) ) ) THEN
         INFO = -1
      ELSE IF ( K.LT.0 ) THEN
         INFO = -2
      ELSE IF ( L.LT.0 ) THEN
         INFO = -3
      ELSE IF ( M.LT.0 ) THEN
         INFO = -4
      ELSE IF ( N.LT.0 ) THEN
         INFO = -5
      ELSE IF ( P*L.GE.MIN( M*K, N*L ) + L .OR. P.LT.0 ) THEN
         INFO = -6
      ELSE IF ( ( P + S )*L.GE.MIN( M*K, N*L ) + L .OR. S.LT.0  ) THEN
         INFO = -7
      ELSE IF ( LDTC.LT.MAX( 1, M*K ) ) THEN
         INFO = -9
      ELSE IF ( LDTR.LT.MAX( 1, K ) ) THEN
         INFO = -11
      ELSE IF ( LDQ.LT.1 .OR. ( COMPQ .AND. LDQ.LT.M*K ) ) THEN
         INFO = -13
      ELSE IF ( LDR.LT.MAX( 1, MIN( N, N - P + 1 )*L ) ) THEN
         INFO = -15
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -17
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'MB02JD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( M, N, K*L, S ) .EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Catch M*K <= L.
C
      WRKOPT = 1
      IF ( M*K.LE.L ) THEN
         CALL DLACPY( 'All', M*K, L, TC, LDTC, DWORK, M*K )
         PDW = M*K*L + 1
         CALL DGEQRF( M*K, L, DWORK, M*K, DWORK(PDW),
     $                DWORK(PDW+M*K), LDWORK-PDW-M*K+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+M*K) ) + PDW + M*K - 1 )
         CALL MA02AD( 'Upper part', M*K, L, DWORK, M*K, R, LDR )
         CALL DORGQR( M*K, M*K, M*K, DWORK, M*K, DWORK(PDW),
     $                DWORK(PDW+M*K), LDWORK-PDW-M*K+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+M*K) ) + PDW + M*K - 1 )
         IF ( COMPQ ) THEN
            CALL DLACPY( 'All', M*K, M*K, DWORK, M*K, Q, LDQ )
         END IF
         PDW = M*K*M*K + 1
         IF ( N.GT.1 ) THEN
            CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, M*K, ONE,
     $                   ZERO, TC, LDTC, TR, LDTR, DWORK, M*K, R(L+1,1),
     $                   LDR, DWORK(PDW), LDWORK-PDW+1, IERR )
         END IF
         WRKOPT   = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         DWORK(1) = DBLE( WRKOPT )
         RETURN
      END IF
C
C     Compute the generator if P = 0.
C
      IF ( P.EQ.0 ) THEN
C
C        1st column of the generator.
C
         IF ( COMPQ ) THEN
            CALL DLACPY( 'All', M*K, L, TC, LDTC, Q, LDQ )
            CALL DGEQRF( M*K, L, Q, LDQ, DWORK, DWORK(L+1),
     $                   LDWORK-L, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(L+1) ) + L )
            CALL MA02AD( 'Upper part', L, L, Q, LDQ, R, LDR )
            CALL DORGQR( M*K, L, L, Q, LDQ, DWORK, DWORK(L+1), LDWORK-L,
     $                  IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(L+1) ) + L )
            IF ( N.GT.1 ) THEN
               CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, L, ONE,
     $                      ZERO, TC, LDTC, TR, LDTR, Q, LDQ, R(L+1,1),
     $                      LDR, DWORK, LDWORK, IERR )
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
         ELSE
            PDW = M*K*L + 1
            CALL DLACPY( 'All', M*K, L, TC, LDTC, DWORK, M*K )
            CALL DGEQRF( M*K, L, DWORK, M*K, DWORK(PDW), DWORK(PDW+L),
     $                   LDWORK-PDW-L+1, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+L) ) + PDW + L - 1 )
            CALL MA02AD( 'Upper part', L, L, DWORK, M*K, R, LDR )
            CALL DORGQR( M*K, L, L, DWORK, M*K, DWORK(PDW),
     $                   DWORK(PDW+L), LDWORK-PDW-L+1, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+L) ) + PDW + L - 1 )
            IF ( N.GT.1 ) THEN
               CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, L, ONE,
     $                      ZERO, TC, LDTC, TR, LDTR, DWORK, M*K,
     $                      R(L+1,1), LDR, DWORK(PDW), LDWORK-PDW+1,
     $                      IERR )
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         END IF
C
C        Quick return if N = 1.
C
         IF ( N.EQ.1 ) THEN
            DWORK(1) = DBLE( WRKOPT )
            RETURN
         END IF
C
C        2nd column of the generator.
C
         PNR = ( N - 1 )*L*K + 2
         CALL MA02AD( 'All', K, (N-1)*L, TR, LDTR, DWORK(2), (N-1)*L )
C
C        3rd and 4th column of the generator.
C
         CALL DLACPY( 'All', (N-1)*L, L, R(L+1,1), LDR, DWORK(PNR),
     $                (N-1)*L )
         PT  = ( M - 1 )*K + 1
         PDW = PNR + ( N - 1 )*L*L
C
         DO 10  I = 1, MIN( M, N-1 )
            CALL MA02AD( 'All', K, L, TC(PT,1), LDTC, DWORK(PDW),
     $                   (N-1)*L )
            PT  = PT  - K
            PDW = PDW + L
   10    CONTINUE
C
         PT = 1
C
         DO 20  I = M + 1, N - 1
            CALL MA02AD( 'All', K, L, TR(1,PT), LDTR, DWORK(PDW),
     $                   (N-1)*L )
            PT  = PT  + L
            PDW = PDW + L
   20    CONTINUE
C
         IF ( COMPQ ) THEN
            PDQ = ( 2*K + L )*( N - 1 )*L + 2
            PDW = ( 2*K + L )*( ( N - 1 )*L + M*K ) + 2
            PNQ = PDQ + M*K*K
            CALL DLASET( 'All', K, K, ZERO, ONE, DWORK(PDQ), M*K )
            CALL DLASET( 'All', (M-1)*K, K, ZERO, ZERO, DWORK(PDQ+K),
     $                   M*K )
            CALL DLACPY( 'All', M*K, L, Q, LDQ, DWORK(PNQ), M*K )
            CALL DLASET( 'All', M*K, K, ZERO, ZERO, DWORK(PNQ+M*L*K),
     $                   M*K )
         ELSE
            PDW = ( 2*K + L )*( N - 1 )*L + 2
         END IF
         PRE  = 1
         STPS = S - 1
      ELSE
C
C        Set workspace pointers.
C
         PNR = ( N - 1 )*L*K + 2
         IF ( COMPQ ) THEN
            PDQ = ( 2*K + L )*( N - 1 )*L + 2
            PDW = ( 2*K + L )*( ( N - 1 )*L + M*K ) + 2
            PNQ = PDQ + M*K*K
         ELSE
            PDW = ( 2*K + L )*( N - 1 )*L + 2
         END IF
         PRE  = P
         STPS = S
      END IF
C
C     Determine suitable size for the block Housholder reflectors.
C
      IF ( COMPQ ) THEN
         LEN = MAX( L + M*K, ( N - PRE + 1 )*L )
      ELSE
         LEN = ( N - PRE + 1 )*L
      END IF
      NB = MIN( ILAENV( 1, 'DGELQF', ' ', LEN, L, -1, -1 ), L )
      KK = PDW + 6*L - 1
      WRKOPT = MAX( WRKOPT, KK + LEN*NB )
      KK = LDWORK - KK
      IF ( KK.LT.LEN*NB )  NB = KK / LEN
      NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', LEN, L, -1, -1 ) )
      IF ( NB.LT.NBMIN )  NB = 0
      COLR = L + 1
C
C     Generator reduction process.
C
      LEN  = ( N - PRE )*L
      SHFR = ( PRE - 1 )*L
      DO 30  I = PRE, PRE + STPS - 1
C
C        IF M*K < N*L the last block might have less than L columns.
C
         KK = MIN( L, M*K - I*L )
         CALL DLACPY( 'Lower', LEN, KK, R(COLR-L,COLR-L), LDR,
     $                R(COLR,COLR), LDR )
         CALL MB02CU( 'Column', KK, KK+K, L+K, NB, R(COLR,COLR), LDR,
     $                 DWORK(SHFR+2), (N-1)*L, DWORK(PNR+SHFR), (N-1)*L,
     $                 RNK, IPVT, DWORK(PDW), ZERO, DWORK(PDW+6*L),
     $                 LDWORK-PDW-6*L+1, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The rank condition is (numerically) not
C                          satisfied.
C
            INFO = 1
            RETURN
         END IF
         IF ( LEN.GT.KK ) THEN
            CALL MB02CV( 'Column', 'NoStructure', KK, LEN-KK, KK+K, L+K,
     $                   NB, -1, R(COLR,COLR), LDR, DWORK(SHFR+2),
     $                   (N-1)*L, DWORK(PNR+SHFR), (N-1)*L,
     $                   R(COLR+KK,COLR), LDR, DWORK(SHFR+KK+2),
     $                   (N-1)*L, DWORK(PNR+SHFR+KK), (N-1)*L,
     $                   DWORK(PDW), DWORK(PDW+6*L), LDWORK-PDW-6*L+1,
     $                   IERR )
         END IF
         IF ( COMPQ ) THEN
            CALL DLASET( 'All', K, KK, ZERO, ZERO, Q(1,COLR), LDQ )
            IF ( M.GT.1 ) THEN
               CALL DLACPY( 'All', (M-1)*K, KK, Q(1,COLR-L), LDQ,
     $                      Q(K+1,COLR), LDQ )
            END IF
            CALL MB02CV( 'Column', 'NoStructure', KK, M*K, KK+K, L+K,
     $                   NB, -1, R(COLR,COLR), LDR, DWORK(SHFR+2),
     $                   (N-1)*L, DWORK(PNR+SHFR), (N-1)*L, Q(1,COLR),
     $                   LDQ, DWORK(PDQ), M*K, DWORK(PNQ), M*K,
     $                   DWORK(PDW), DWORK(PDW+6*L), LDWORK-PDW-6*L+1,
     $                   IERR )
         END IF
         LEN  = LEN  - L
         COLR = COLR + L
         SHFR = SHFR + L
   30 CONTINUE
C
      DWORK(1) = DBLE( WRKOPT )
      RETURN
C
C *** Last line of MB02JD ***
      END
