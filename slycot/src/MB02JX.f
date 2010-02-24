      SUBROUTINE MB02JX( JOB, K, L, M, N, TC, LDTC, TR, LDTR, RNK, Q,
     $                   LDQ, R, LDR, JPVT, TOL1, TOL2, DWORK, LDWORK,
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
C     To compute a low rank QR factorization with column pivoting of a
C     K*M-by-L*N block Toeplitz matrix T with blocks of size (K,L);
C     specifically,
C                                     T
C                           T P =  Q R ,
C
C     where R is lower trapezoidal, P is a block permutation matrix
C     and Q^T Q = I. The number of columns in R is equivalent to the
C     numerical rank of T with respect to the given tolerance TOL1.
C     Note that the pivoting scheme is local, i.e., only columns
C     belonging to the same block in T are permuted.
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
C     TC      (input) DOUBLE PRECISION array, dimension (LDTC, L)
C             The leading M*K-by-L part of this array must contain
C             the first block column of T.
C
C     LDTC    INTEGER
C             The leading dimension of the array TC.
C             LDTC >= MAX(1,M*K).
C
C     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L)
C             The leading K-by-(N-1)*L part of this array must contain
C             the first block row of T without the leading K-by-L
C             block.
C
C     LDTR    INTEGER
C             The leading dimension of the array TR.  LDTR >= MAX(1,K).
C
C     RNK     (output)  INTEGER
C             The number of columns in R, which is equivalent to the
C             numerical rank of T.
C
C     Q       (output)  DOUBLE PRECISION array, dimension (LDQ,RNK)
C             If JOB = 'Q', then the leading M*K-by-RNK part of this
C             array contains the factor Q.
C             If JOB = 'R', then this array is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.
C             LDQ >= MAX(1,M*K),  if JOB = 'Q';
C             LDQ >= 1,           if JOB = 'R'.
C
C     R       (output)  DOUBLE PRECISION array, dimension (LDR,RNK)
C             The leading N*L-by-RNK part of this array contains the
C             lower trapezoidal factor R.
C
C     LDR     INTEGER
C             The leading dimension of the array R.
C             LDR >= MAX(1,N*L)
C
C     JPVT    (output)  INTEGER array, dimension (MIN(M*K,N*L))
C             This array records the column pivoting performed.
C             If JPVT(j) = k, then the j-th column of T*P was
C             the k-th column of T.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If TOL1 >= 0.0, the user supplied diagonal tolerance;
C             if TOL1 < 0.0, a default diagonal tolerance is used.
C
C     TOL2    DOUBLE PRECISION
C             If TOL2 >= 0.0, the user supplied offdiagonal tolerance;
C             if TOL2 < 0.0, a default offdiagonal tolerance is used.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK;  DWORK(2) and DWORK(3) return the used values
C             for TOL1 and TOL2, respectively.
C             On exit, if INFO = -19,  DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 3, ( M*K + ( N - 1 )*L )*( L + 2*K ) + 9*L
C                                 + MAX(M*K,(N-1)*L) ),    if JOB = 'Q';
C             LDWORK >= MAX( 3, ( N - 1 )*L*( L + 2*K + 1 ) + 9*L,
C                                 M*K*( L + 1 ) + L ),     if JOB = 'R'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  due to perturbations induced by roundoff errors, or
C                   removal of nearly linearly dependent columns of the
C                   generator, the Schur algorithm encountered a
C                   situation where a diagonal element in the negative
C                   generator is larger in magnitude than the
C                   corresponding diagonal element in the positive
C                   generator (modulo TOL1);
C             = 2:  due to perturbations induced by roundoff errors, or
C                   removal of nearly linearly dependent columns of the
C                   generator, the Schur algorithm encountered a
C                   situation where diagonal elements in the positive
C                   and negative generator are equal in magnitude
C                   (modulo TOL1), but the offdiagonal elements suggest
C                   that these columns are not linearly dependent
C                   (modulo TOL2*ABS(diagonal element)).
C
C     METHOD
C
C     Householder transformations and modified hyperbolic rotations
C     are used in the Schur algorithm [1], [2].
C     If, during the process, the hyperbolic norm of a row in the
C     leading part of the generator is found to be less than or equal
C     to TOL1, then this row is not reduced. If the difference of the
C     corresponding columns has a norm less than or equal to TOL2 times
C     the magnitude of the leading element, then this column is removed
C     from the generator, as well as from R. Otherwise, the algorithm
C     breaks down. TOL1 is set to norm(TC)*sqrt(eps) and TOL2 is set
C     to N*L*sqrt(eps) by default.
C     If M*K > L, the columns of T are permuted so that the diagonal
C     elements in one block column of R have decreasing magnitudes.
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
C     The algorithm requires 0(K*RNK*L*M*N) floating point operations.
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
      INTEGER           INFO, K, L, LDQ, LDR, LDTC, LDTR, LDWORK, M, N,
     $                  RNK
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), Q(LDQ,*), R(LDR,*), TC(LDTC,*),
     $                  TR(LDTR,*)
      INTEGER           JPVT(*)
C     .. Local Scalars ..
      LOGICAL           COMPQ, LAST
      INTEGER           CPCOL, GAP, I, IERR, J, JJ, JWORK, KK, LEN, MK,
     $                  NZC, PDP, PDQ, PDW, PNQ, PNR, PP, PPR, PT, RDEF,
     $                  RRDF, RRNK, WRKMIN, WRKOPT
      DOUBLE PRECISION  LTOL1, LTOL2, NRM, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DNRM2
      EXTERNAL          DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEQP3, DGEQRF, DLACPY, DLASET,
     $                  DORGQR, DSCAL, DSWAP, DTRMV, MA02AD, MB02CU,
     $                  MB02CV, MB02KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO   = 0
      WRKOPT = 3
      MK     = M*K
      COMPQ  = LSAME( JOB, 'Q' )
      IF ( COMPQ ) THEN
         WRKMIN = MAX( 3, ( MK + ( N - 1 )*L )*( L + 2*K ) + 9*L +
     $                    MAX( MK, ( N - 1 )*L ) )
      ELSE
         WRKMIN = MAX( 3, MAX ( ( N - 1 )*L*( L + 2*K + 1 ) + 9*L,
     $                    MK*( L + 1 ) + L ) )
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
      ELSE IF ( LDTC.LT.MAX( 1, MK ) ) THEN
         INFO = -7
      ELSE IF ( LDTR.LT.MAX( 1, K ) ) THEN
         INFO = -9
      ELSE IF ( LDQ.LT.1 .OR. ( COMPQ .AND. LDQ.LT.MK ) ) THEN
         INFO = -12
      ELSE IF ( LDR.LT.MAX( 1, N*L ) ) THEN
         INFO = -14
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -19
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02JX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( M, N, K, L ).EQ.0 ) THEN
         RNK      = 0
         DWORK(1) = DBLE( WRKOPT )
         DWORK(2) = ZERO
         DWORK(3) = ZERO
         RETURN
      END IF
C
      WRKOPT = WRKMIN
C
      IF ( MK.LE.L ) THEN
C
C        Catch M*K <= L.
C
         CALL DLACPY( 'All', MK, L, TC, LDTC, DWORK, MK )
         PDW   = MK*L + 1
         JWORK = PDW  + MK
         CALL DGEQRF( MK, L, DWORK, MK, DWORK(PDW), DWORK(JWORK),
     $                LDWORK-JWORK+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         CALL MA02AD( 'Upper part', MK, L, DWORK, MK, R, LDR )
         CALL DORGQR( MK, MK, MK, DWORK, MK, DWORK(PDW),
     $                DWORK(JWORK), LDWORK-JWORK+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IF ( COMPQ )
     $      CALL DLACPY( 'All', MK, MK, DWORK, MK, Q, LDQ )
         PDW = MK*MK + 1
         IF ( N.GT.1 ) THEN
            CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, MK, ONE,
     $                   ZERO, TC, LDTC, TR, LDTR, DWORK, MK, R(L+1,1),
     $                   LDR, DWORK(PDW), LDWORK-PDW+1, IERR )
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
         DO 10  I = 1, MK
            JPVT(I) = I
   10    CONTINUE
C
         RNK      = MK
         DWORK(1) = DBLE( WRKOPT )
         DWORK(2) = ZERO
         DWORK(3) = ZERO
         RETURN
      END IF
C
C     Compute the generator:
C
C     1st column of the generator.
C
      DO 20  I = 1, L
         JPVT(I) = 0
   20 CONTINUE
C
      LTOL1 = TOL1
      LTOL2 = TOL2
C
      IF ( COMPQ ) THEN
         CALL DLACPY( 'All', MK, L, TC, LDTC, Q, LDQ )
         CALL DGEQP3( MK, L, Q, LDQ, JPVT, DWORK, DWORK(L+1),
     $                LDWORK-L, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(L+1) ) + L )
C
         IF ( LTOL1.LT.ZERO ) THEN
C
C           Compute default tolerance LTOL1.
C
C           Estimate the 2-norm of the first block column of the
C           matrix with 5 power iterations.
C
            TEMP = ONE / SQRT( DBLE( L ) )
            CALL DLASET( 'All', L, 1, TEMP, TEMP, DWORK(L+1), 1 )
C
            DO 30  I = 1, 5
               CALL DTRMV( 'Upper', 'NonTranspose', 'NonUnit', L, Q,
     $                     LDQ, DWORK(L+1), 1 )
               CALL DTRMV( 'Upper', 'Transpose', 'NonUnit', L, Q, LDQ,
     $                     DWORK(L+1), 1 )
               NRM = DNRM2( L, DWORK(L+1), 1 )
               CALL DSCAL( L, ONE/NRM, DWORK(L+1), 1 )
   30       CONTINUE
C
            LTOL1 = SQRT( NRM*DLAMCH( 'Epsilon' ) )
         END IF
C
         I = L
C
   40    CONTINUE
         IF ( ABS( Q(I,I) ).LE.LTOL1 ) THEN
            I = I - 1
            IF ( I.GT.0 )  GO TO 40
         END IF
C
         RRNK = I
         RRDF = L - RRNK
         CALL MA02AD( 'Upper', RRNK, L, Q, LDQ, R, LDR )
         IF ( RRNK.GT.1 )
     $      CALL DLASET( 'Upper', L-1, RRNK-1, ZERO, ZERO, R(1,2), LDR )
         CALL DORGQR( MK, L, RRNK, Q, LDQ, DWORK, DWORK(L+1), LDWORK-L,
     $                IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(L+1) ) + L )
         IF ( N.GT.1 ) THEN
            CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, RRNK, ONE,
     $                   ZERO, TC, LDTC, TR, LDTR, Q, LDQ, R(L+1,1),
     $                   LDR, DWORK, LDWORK, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
         END IF
C
      ELSE
C
         PDW   = MK*L + 1
         JWORK = PDW  + L
         CALL DLACPY( 'All', MK, L, TC, LDTC, DWORK, MK )
         CALL DGEQP3( MK, L, DWORK, MK, JPVT, DWORK(PDW),
     $                DWORK(JWORK), LDWORK-JWORK+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
         IF ( LTOL1.LT.ZERO ) THEN
C
C           Compute default tolerance LTOL1.
C
C           Estimate the 2-norm of the first block column of the
C           matrix with 5 power iterations.
C
            TEMP = ONE / SQRT( DBLE( L ) )
            CALL DLASET( 'All', L, 1, TEMP, TEMP, DWORK(JWORK), 1 )
C
            DO 50  I = 1, 5
               CALL DTRMV( 'Upper', 'NonTranspose', 'NonUnit', L, DWORK,
     $                     MK, DWORK(JWORK), 1 )
               CALL DTRMV( 'Upper', 'Transpose', 'NonUnit', L, DWORK,
     $                     MK, DWORK(JWORK), 1 )
               NRM = DNRM2( L, DWORK(JWORK), 1 )
               CALL DSCAL( L, ONE/NRM, DWORK(JWORK), 1 )
   50       CONTINUE
C
            LTOL1 = SQRT( NRM*DLAMCH( 'Epsilon' ) )
         END IF
C
         RRNK = L
         I = ( L - 1 )*MK + L
C
   60    CONTINUE
         IF ( ABS( DWORK(I) ).LE.LTOL1 ) THEN
            RRNK = RRNK - 1
            I = I - MK - 1
            IF ( I.GT.0 )  GO TO 60
         END IF
C
         RRDF = L - RRNK
         CALL MA02AD( 'Upper part', RRNK, L, DWORK, MK, R, LDR )
         IF ( RRNK.GT.1 )
     $      CALL DLASET( 'Upper', L-1, RRNK-1, ZERO, ZERO, R(1,2), LDR )
         CALL DORGQR( MK, L, RRNK, DWORK, MK, DWORK(PDW),
     $                DWORK(JWORK), LDWORK-JWORK+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IF ( N.GT.1 ) THEN
            CALL MB02KD( 'Row', 'Transpose', K, L, M, N-1, RRNK, ONE,
     $                   ZERO, TC, LDTC, TR, LDTR, DWORK, MK, R(L+1,1),
     $                   LDR, DWORK(PDW), LDWORK-PDW+1, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         END IF
      END IF
C
C     Quick return if N = 1.
C
      IF ( N.EQ.1 ) THEN
         RNK      = RRNK
         DWORK(1) = DBLE( WRKOPT )
         DWORK(2) = LTOL1
         DWORK(3) = ZERO
         RETURN
      END IF
C
C     Compute default tolerance LTOL2.
C
      IF ( LTOL2.LT.ZERO )
     $   LTOL2 = DBLE( N*L )*SQRT( DLAMCH( 'Epsilon' ) )
C
      DO 70  J = 1, L
         CALL DCOPY( RRNK, R(J,1), LDR, R(L+JPVT(J),RRNK+1), LDR )
   70 CONTINUE
C
      IF ( N.GT.2 )
     $   CALL DLACPY( 'All', (N-2)*L, RRNK, R(L+1,1), LDR,
     $                R(2*L+1,RRNK+1), LDR )
C
C     2nd column of the generator.
C
      IF ( RRDF.GT.0 )
     $   CALL MA02AD( 'All', MIN( RRDF, K ), (N-1)*L, TR, LDTR,
     $                R(L+1,2*RRNK+1), LDR )
      IF ( K.GT.RRDF )
     $   CALL MA02AD( 'All', K-RRDF, (N-1)*L, TR(RRDF+1,1), LDTR, DWORK,
     $                (N-1)*L )
C
C     3rd column of the generator.
C
      PNR = ( N - 1 )*L*MAX( 0, K-RRDF ) + 1
      CALL DLACPY( 'All', (N-1)*L, RRNK, R(L+1,1), LDR, DWORK(PNR),
     $             (N-1)*L )
C
C     4th column of the generator.
C
      PDW = PNR + ( N - 1 )*L*RRNK
      PT  = ( M - 1 )*K + 1
C
      DO 80  I = 1, MIN( M, N-1 )
         CALL MA02AD( 'All', K, L, TC(PT,1), LDTC, DWORK(PDW), (N-1)*L )
         PT  = PT  - K
         PDW = PDW + L
   80 CONTINUE
C
      PT = 1
C
      DO 90  I = M + 1, N - 1
         CALL MA02AD( 'All', K, L, TR(1,PT), LDTR, DWORK(PDW), (N-1)*L )
         PT  = PT  + L
         PDW = PDW + L
   90 CONTINUE
C
      IF ( COMPQ ) THEN
         PDQ = PNR + ( N - 1 )*L*( RRNK + K )
         PNQ = PDQ + MK*MAX( 0, K-RRDF )
         PDW = PNQ + MK*( RRNK + K )
         CALL DLACPY( 'All', MK, RRNK, Q, LDQ, DWORK(PNQ), MK )
         IF ( M.GT.1 )
     $      CALL DLACPY( 'All', (M-1)*K, RRNK, Q, LDQ, Q(K+1,RRNK+1),
     $                   LDQ )
         CALL DLASET( 'All', K, RRNK, ZERO, ZERO, Q(1,RRNK+1), LDQ )
      IF ( RRDF.GT.0 )
     $      CALL DLASET( 'All', MK, RRDF, ZERO, ONE, Q(1,2*RRNK+1),
     $                   LDQ )
         CALL DLASET( 'All', RRDF, MAX( 0, K-RRDF ), ZERO, ZERO,
     $                DWORK(PDQ), MK )
         CALL DLASET( 'All', M*K-RRDF, MAX( 0, K-RRDF ), ZERO, ONE,
     $                DWORK(PDQ+RRDF), MK )
         CALL DLASET( 'All', MK, K, ZERO, ZERO, DWORK(PNQ+MK*RRNK), MK )
      ELSE
         PDW = PNR + ( N - 1 )*L*( RRNK + K )
      END IF
      PPR  = 1
      RNK  = RRNK
      RDEF = RRDF
      LEN  = N*L
      GAP  = N*L - MIN( N*L, MK )
C
C     KK is the number of columns in the leading part of the
C     generator. After sufficiently many rank drops or if
C     M*K < N*L it may be less than L.
C
      KK = MIN( L+K-RDEF, L )
      KK = MIN( KK, MK-L )
C
C     Generator reduction process.
C
      DO 190  I = L + 1, MIN( MK, N*L ), L
         IF ( I+L.LE.MIN( MK, N*L ) ) THEN
            LAST = .FALSE.
         ELSE
            LAST = .TRUE.
         END IF
         PP  = KK  + MAX( K - RDEF, 0 )
         LEN = LEN - L
         CALL MB02CU( 'Deficient', KK, PP, L+K-RDEF, -1, R(I,RNK+1),
     $                 LDR, DWORK(PPR), (N-1)*L, DWORK(PNR), (N-1)*L,
     $                 RRNK, JPVT(I), DWORK(PDW), LTOL1, DWORK(PDW+5*L),
     $                 LDWORK-PDW-5*L+1, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The current generator is indefinite.
C
            INFO = 1
            RETURN
         END IF
C
C        Apply pivoting to other columns of R.
C
         PDP = PDW + 6*L - I
C
         DO 100  J = I, I + KK - 1
            JPVT(J) = JPVT(J) + I - 1
            DWORK(PDP+JPVT(J)) = DBLE(J)
  100    CONTINUE
C
         DO 120  J = I, I + KK - 1
            TEMP = DBLE(J)
            JJ = J-1
C
  110       CONTINUE
            JJ = JJ + 1
            IF ( DWORK(PDP+JJ).NE.TEMP )  GO TO 110
C
            IF ( JJ.NE.J ) THEN
               DWORK(PDP+JJ) = DWORK(PDP+J)
               CALL DSWAP( RNK, R(J,1), LDR, R(JJ,1), LDR )
            END IF
  120    CONTINUE
C
         DO 130  J = I + KK, I + L - 1
            JPVT(J) = J
  130    CONTINUE
C
C        Apply reduction to other rows of R.
C
         IF ( LEN.GT.KK ) THEN
            CALL MB02CV( 'Deficient', 'NoStructure', KK, LEN-KK, PP,
     $                   L+K-RDEF, -1, RRNK, R(I,RNK+1), LDR,
     $                   DWORK(PPR), (N-1)*L, DWORK(PNR), (N-1)*L,
     $                   R(I+KK,RNK+1), LDR, DWORK(PPR+KK), (N-1)*L,
     $                   DWORK(PNR+KK), (N-1)*L, DWORK(PDW),
     $                   DWORK(PDW+5*L), LDWORK-PDW-5*L+1, IERR )
         END IF
C
C        Apply reduction to Q.
C
         IF ( COMPQ ) THEN
            CALL MB02CV( 'Deficient', 'NoStructure', KK, MK, PP,
     $                   L+K-RDEF, -1, RRNK, R(I,RNK+1), LDR,
     $                   DWORK(PPR), (N-1)*L, DWORK(PNR), (N-1)*L,
     $                   Q(1,RNK+1), LDQ, DWORK(PDQ), MK, DWORK(PNQ),
     $                   MK, DWORK(PDW), DWORK(PDW+5*L),
     $                   LDWORK-PDW-5*L+1, IERR )
         END IF
C
C        Inspection of the rank deficient columns:
C        Look for small diagonal entries.
C
         NZC = 0
C
         DO 140  J = KK, RRNK + 1, -1
            IF ( ABS( R(I+J-1,RNK+J) ).LE.LTOL1 )  NZC = NZC + 1
  140    CONTINUE
C
C        The last NZC columns of the generator cannot be removed.
C        Now, decide whether for the other rank deficient columns
C        it is safe to remove.
C
         PT = PNR
C
         DO 150  J = RRNK + 1, KK - NZC
            TEMP = R(I+J-1,RNK+J)
            CALL DSCAL( LEN-J-GAP, TEMP, R(I+J,RNK+J), 1 )
            CALL DAXPY( LEN-J-GAP, -DWORK(PT+J-1), DWORK(PT+J), 1,
     $                  R(I+J,RNK+J), 1 )
            IF ( DNRM2( LEN-J-GAP, R(I+J,RNK+J), 1 )
     $          .GT.LTOL2*ABS( TEMP ) ) THEN
C
C              Unlucky case:
C              It is neither advisable to remove the whole column nor
C              possible to remove the diagonal entries by Hyperbolic
C              rotations.
C
               INFO = 2
               RETURN
            END IF
            PT = PT + ( N - 1 )*L
  150    CONTINUE
C
C        Annihilate unwanted elements in the factor R.
C
         RRDF = KK - RRNK
         CALL DLASET( 'All', I-1, RRNK, ZERO, ZERO, R(1,RNK+1), LDR )
         CALL DLASET( 'Upper', L-1, RRNK-1, ZERO, ZERO, R(I,RNK+2),
     $                LDR )
C
C        Construct the generator for the next step.
C
         IF ( .NOT.LAST ) THEN
C
C           Compute KK for the next step.
C
            KK = MIN( L+K-RDEF-RRDF+NZC, L )
            KK = MIN( KK, MK-I-L+1 )
C
            IF ( KK.LE.0 ) THEN
               RNK = RNK + RRNK
               GO TO 200
            END IF
C
            CALL DLASET( 'All', L, RRDF, ZERO, ZERO, R(I,RNK+RRNK+1),
     $                   LDR )
C
C           The columns with small diagonal entries form parts of the
C           new positive generator.
C
            IF ( ( RRDF-NZC ).GT.0 .AND. NZC.GT.0 ) THEN
               CPCOL = MIN( NZC, KK )
C
               DO 160  J = RNK + RRNK + 1, RNK + RRNK + CPCOL
                  CALL DCOPY( LEN-L, R(I+L,J+RRDF-NZC), 1,
     $                        R(I+L,J), 1 )
  160          CONTINUE
C
            END IF
C
C           Construct the leading parts of the positive generator.
C
            CPCOL = MIN( RRNK, KK-NZC )
            IF ( CPCOL.GT.0 ) THEN
C
               DO 170  J = I, I + L - 1
                  CALL DCOPY( CPCOL, R(J,RNK+1), LDR,
     $                        R(JPVT(J)+L,RNK+RRNK+NZC+1), LDR )
  170          CONTINUE
C
               IF ( LEN.GT.2*L ) THEN
                  CALL DLACPY( 'All', LEN-2*L, CPCOL, R(I+L,RNK+1), LDR,
     $                         R(I+2*L,RNK+RRNK+NZC+1), LDR )
               END IF
            END IF
            PPR = PPR + L
C
C           Refill the leading parts of the positive generator.
C
            CPCOL = MIN( K-RDEF, KK-RRNK-NZC )
            IF ( CPCOL.GT.0 ) THEN
               CALL DLACPY( 'All', LEN-L, CPCOL, DWORK(PPR), (N-1)*L,
     $                      R(I+L,RNK+2*RRNK+NZC+1), LDR )
               PPR = PPR + CPCOL*( N - 1 )*L
            END IF
            PNR = PNR + ( RRDF - NZC )*( N - 1 )*L + L
C
C           Do the same things for Q.
C
            IF ( COMPQ ) THEN
               IF ( ( RRDF - NZC ).GT.0 .AND. NZC.GT.0 ) THEN
                  CPCOL = MIN( NZC, KK )
C
                  DO 180  J = RNK + RRNK + 1, RNK + RRNK + CPCOL
                     CALL DCOPY( MK, Q(1,J+RRDF-NZC), 1, Q(1,J), 1 )
  180             CONTINUE
C
               END IF
               CPCOL = MIN( RRNK, KK-NZC )
               IF ( CPCOL.GT.0 ) THEN
                  CALL DLASET( 'All', K, CPCOL, ZERO, ZERO,
     $                         Q(1,RNK+RRNK+NZC+1), LDQ )
                  IF ( M.GT.1 )
     $               CALL DLACPY( 'All', (M-1)*K, CPCOL, Q(1,RNK+1),
     $                            LDQ, Q(K+1,RNK+RRNK+NZC+1), LDQ )
               END IF
               CPCOL = MIN( K-RDEF, KK-RRNK-NZC )
               IF ( CPCOL.GT.0 ) THEN
                  CALL DLACPY( 'All', MK, CPCOL, DWORK(PDQ), MK,
     $                         Q(1,RNK+2*RRNK+NZC+1), LDQ )
                  PDQ = PDQ + CPCOL*MK
               END IF
               PNQ = PNQ + ( RRDF - NZC )*MK
            END IF
         END IF
         RNK  = RNK  + RRNK
         RDEF = RDEF + RRDF - NZC
  190 CONTINUE
C
  200 CONTINUE
      DWORK(1) = DBLE( WRKOPT )
      DWORK(2) = LTOL1
      DWORK(3) = LTOL2
C
C *** Last line of MB02JX ***
      END
