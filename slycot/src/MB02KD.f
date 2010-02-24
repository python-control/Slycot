      SUBROUTINE MB02KD( LDBLK, TRANS, K, L, M, N, R, ALPHA, BETA,
     $                   TC, LDTC, TR, LDTR, B, LDB, C, LDC, DWORK,
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
C     To compute the matrix product
C
C               C = alpha*op( T )*B + beta*C,
C
C     where alpha and beta are scalars and T is a block Toeplitz matrix
C     specified by its first block column TC and first block row TR;
C     B and C are general matrices of appropriate dimensions.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LDBLK   CHARACTER*1
C             Specifies where the (1,1)-block of T is stored, as
C             follows:
C             = 'C':  in the first block of TC;
C             = 'R':  in the first block of TR.
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( T ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( T ) = T;
C             = 'T':  op( T ) = T';
C             = 'C':  op( T ) = T'.
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
C             The number of blocks in the first block row of T.  N >= 0.
C
C     R       (input) INTEGER
C             The number of columns in B and C.  R >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then TC, TR and B
C             are not referenced.
C
C     BETA    (input) DOUBLE PRECISION
C             The scalar beta. When beta is zero then C need not be set
C             before entry.
C
C     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L)
C             On entry with LDBLK = 'C', the leading M*K-by-L part of
C             this array must contain the first block column of T.
C             On entry with LDBLK = 'R', the leading (M-1)*K-by-L part
C             of this array must contain the 2nd to the M-th blocks of
C             the first block column of T.
C
C     LDTC    INTEGER
C             The leading dimension of the array TC.
C             LDTC >= MAX(1,M*K),      if LDBLK = 'C';
C             LDTC >= MAX(1,(M-1)*K),  if LDBLK = 'R'.
C
C     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,k)
C             where k is (N-1)*L when LDBLK = 'C' and is N*L when
C             LDBLK = 'R'.
C             On entry with LDBLK = 'C', the leading K-by-(N-1)*L part
C             of this array must contain the 2nd to the N-th blocks of
C             the first block row of T.
C             On entry with LDBLK = 'R', the leading K-by-N*L part of
C             this array must contain the first block row of T.
C
C     LDTR    INTEGER
C             The leading dimension of the array TR.  LDTR >= MAX(1,K).
C
C     B       (input)  DOUBLE PRECISION array, dimension (LDB,R)
C             On entry with TRANS = 'N', the leading N*L-by-R part of
C             this array must contain the matrix B.
C             On entry with TRANS = 'T' or TRANS = 'C', the leading
C             M*K-by-R part of this array must contain the matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,N*L),  if TRANS = 'N';
C             LDB >= MAX(1,M*K),  if TRANS = 'T' or TRANS = 'C'.
C
C     C       (input/output)  DOUBLE PRECISION array, dimension (LDC,R)
C             On entry with TRANS = 'N', the leading M*K-by-R part of
C             this array must contain the matrix C.
C             On entry with TRANS = 'T' or TRANS = 'C', the leading
C             N*L-by-R part of this array must contain the matrix C.
C             On exit with TRANS = 'N', the leading M*K-by-R part of
C             this array contains the updated matrix C.
C             On exit with TRANS = 'T' or TRANS = 'C', the leading
C             N*L-by-R part of this array contains the updated matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.
C             LDC >= MAX(1,M*K),  if TRANS = 'N';
C             LDC >= MAX(1,N*L),  if TRANS = 'T' or TRANS = 'C'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -19,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 1.
C             For optimum performance LDWORK should be larger.
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
C     For point Toeplitz matrices or sufficiently large block Toeplitz
C     matrices, this algorithm uses convolution algorithms based on
C     the fast Hartley transforms [1]. Otherwise, TC is copied in
C     reversed order into the workspace such that C can be computed from
C     barely M matrix-by-matrix multiplications.
C
C     REFERENCES
C
C     [1] Van Loan, Charles.
C         Computational frameworks for the fast Fourier transform.
C         SIAM, 1992.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O( (K*L+R*L+K*R)*(N+M)*log(N+M) + K*L*R )
C     floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C     March 2004.
C
C     KEYWORDS
C
C     Convolution, elementary matrix operations,
C     fast Hartley transform, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, THREE, FOUR, THOM50
      PARAMETER         ( ZERO  = 0.0D0, ONE  = 1.0D0, TWO    = 2.0D0,
     $                    THREE = 3.0D0, FOUR = 4.0D0, THOM50 = .95D3 )
C     .. Scalar Arguments ..
      CHARACTER         LDBLK, TRANS
      INTEGER           INFO, K, L, LDB, LDC, LDTC, LDTR, LDWORK, M, N,
     $                  R
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), C(LDC,*), DWORK(*), TC(LDTC,*),
     $                  TR(LDTR,*)
C     .. Local Scalars ..
      LOGICAL           FULLC, LMULT, LTRAN
      CHARACTER*1       WGHT
      INTEGER           DIMB, DIMC, I, ICP, ICQ, IERR, IR, J, JJ, KK,
     $                  LEN, LL, LN, METH, MK, NL, P, P1, P2, PB, PC,
     $                  PDW, PP, PT, Q1, Q2, R1, R2, S1, S2, SHFT, WPOS,
     $                  WRKOPT
      DOUBLE PRECISION  CF, COEF, PARAM, SCAL, SF, T1, T2, TH
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DG01OD, DGEMM, DLACPY, DLASET,
     $                  DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, COS, DBLE, MAX, MIN, SIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      FULLC = LSAME( LDBLK, 'C' )
      LTRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      LMULT = ALPHA.NE.ZERO
      MK    = M*K
      NL    = N*L
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( FULLC .OR. LSAME( LDBLK, 'R' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( LTRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( K.LT.0 ) THEN
         INFO = -3
      ELSE IF ( L.LT.0 ) THEN
         INFO = -4
      ELSE IF ( M.LT.0 ) THEN
         INFO = -5
      ELSE IF ( N.LT.0 ) THEN
         INFO = -6
      ELSE IF ( R.LT.0 ) THEN
         INFO = -7
      ELSE IF ( LMULT .AND. FULLC .AND. LDTC.LT.MAX( 1, MK ) ) THEN
         INFO = -11
      ELSE IF ( LMULT .AND. .NOT.FULLC .AND.
     $          LDTC.LT.MAX( 1,( M - 1 )*K ) ) THEN
         INFO = -11
      ELSE IF ( LMULT .AND. LDTR.LT.MAX( 1, K ) ) THEN
         INFO = -13
      ELSE IF ( LMULT .AND. .NOT.LTRAN .AND. LDB.LT.MAX( 1, NL ) ) THEN
         INFO = -15
      ELSE IF ( LMULT .AND. LTRAN .AND. LDB.LT.MAX( 1, MK ) ) THEN
         INFO = -15
      ELSE IF ( .NOT.LTRAN .AND. LDC.LT.MAX( 1, MK ) ) THEN
         INFO = -17
      ELSE IF ( LTRAN .AND. LDC.LT.MAX( 1, NL ) ) THEN
         INFO = -17
      ELSE IF ( LDWORK.LT.1 ) THEN
         DWORK(1) = ONE
         INFO = -19
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02KD', -INFO )
         RETURN
      END IF
C
C     Scale C beforehand.
C
      IF ( BETA.EQ.ZERO ) THEN
         IF ( LTRAN ) THEN
            CALL DLASET( 'All', NL, R, ZERO, ZERO, C, LDC )
         ELSE
            CALL DLASET( 'All', MK, R, ZERO, ZERO, C, LDC )
         END IF
      ELSE IF ( BETA.NE.ONE ) THEN
         IF ( LTRAN ) THEN
C
            DO 10  I = 1, R
               CALL DSCAL( NL, BETA, C(1,I), 1 )
   10       CONTINUE
C
         ELSE
C
            DO 20  I = 1, R
               CALL DSCAL( MK, BETA, C(1,I), 1 )
   20       CONTINUE
C
         END IF
      END IF
C
C     Quick return if possible.
C
      IF ( .NOT.LMULT .OR. MIN( MK, NL, R ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     The parameter PARAM is the watershed between conventional
C     multiplication and convolution. This is of course depending
C     on the used computer architecture. The lower this value is set
C     the more likely the routine will use convolution to compute
C     op( T )*B. Note that if there is enough workspace available,
C     convolution is always used for point Toeplitz matrices.
C
      PARAM = THOM50
C
C     Decide which method to choose, based on the block sizes and
C     the available workspace.
C
      LEN = 1
      P   = 0
C
   30 CONTINUE
      IF ( LEN.LT.M+N-1 ) THEN
         LEN = LEN*2
         P = P + 1
         GO TO 30
      END IF
C
      COEF = THREE*DBLE( M*N )*DBLE( K*L )*DBLE( R ) /
     $             DBLE( LEN*( K*L + L*R + K*R ) )
C
      IF ( FULLC ) THEN
         P1 = MK*L
         SHFT = 0
      ELSE
         P1 = ( M - 1 )*K*L
         SHFT = 1
      END IF
      IF ( K*L.EQ.1 .AND. MIN( M, N ).GT.1 ) THEN
         WRKOPT = LEN*( 2 + R ) - P
         METH   = 3
      ELSE IF ( ( LEN.LT.M*N ) .AND. ( COEF.GE.PARAM ) ) THEN
         WRKOPT = LEN*( K*L + K*R + L*R + 1 ) - P
         METH   = 3
      ELSE
         METH   = 2
         WRKOPT = P1
      END IF
C
      IF ( LDWORK.LT.WRKOPT )  METH = METH - 1
      IF ( LDWORK.LT.P1 )      METH = 1
C
C     Start computations.
C
      IF ( METH.EQ.1 .AND. .NOT.LTRAN ) THEN
C
C        Method 1 is the most unlucky way to multiply Toeplitz matrices
C        with vectors. Due to the memory restrictions it is not
C        possible to flip TC.
C
         PC = 1
C
         DO 50 I = 1, M
            PT = ( I - 1 - SHFT )*K + 1
            PB = 1
C
            DO 40 J = SHFT + 1, I
               CALL DGEMM( 'No Transpose', 'No Transpose', K, R, L,
     $                     ALPHA, TC(PT,1), LDTC, B(PB,1), LDB, ONE,
     $                     C(PC,1), LDC )
               PT = PT - K
               PB = PB + L
   40       CONTINUE
C
            IF ( N.GT.I-SHFT ) THEN
               CALL DGEMM( 'No Transpose', 'No Transpose', K, R,
     $                     (N-I+SHFT)*L, ALPHA, TR, LDTR, B(PB,1), LDB,
     $                     ONE, C(PC,1), LDC )
            END IF
            PC = PC + K
   50    CONTINUE
C
      ELSE IF ( METH.EQ.1 .AND. LTRAN ) THEN
C
         PB = 1
C
         DO 70 I = 1, M
            PT = ( I - 1 - SHFT )*K + 1
            PC = 1
C
            DO 60 J = SHFT + 1, I
               CALL DGEMM( 'Transpose', 'No Transpose', L, R, K, ALPHA,
     $                     TC(PT,1), LDTC, B(PB,1), LDB, ONE, C(PC,1),
     $                     LDC )
               PT = PT - K
               PC = PC + L
   60       CONTINUE
C
            IF ( N.GT.I-SHFT ) THEN
               CALL DGEMM( 'Transpose', 'No Transpose', (N-I+SHFT)*L,
     $                     R, K, ALPHA, TR, LDTR, B(PB,1), LDB, ONE,
     $                     C(PC,1), LDC )
            END IF
            PB = PB + K
   70    CONTINUE
C
      ELSE IF ( METH.EQ.2 .AND. .NOT.LTRAN ) THEN
C
C        In method 2 TC is flipped resulting in less calls to the BLAS
C        routine DGEMM. Actually this seems often to be the best way to
C        multiply with Toeplitz matrices except the point Toeplitz
C        case.
C
         PT = ( M - 1 - SHFT )*K + 1
C
         DO 80  I = 1, ( M - SHFT )*K*L, K*L
            CALL DLACPY( 'All', K, L, TC(PT,1), LDTC, DWORK(I), K )
            PT = PT - K
   80    CONTINUE
C
         PT = ( M - 1 )*K*L + 1
         PC = 1
C
         DO 90  I = 1, M
            CALL DGEMM( 'No Transpose', 'No Transpose', K, R,
     $                  MIN( I-SHFT, N )*L, ALPHA, DWORK(PT), K, B, LDB,
     $                  ONE, C(PC,1), LDC )
            IF ( N.GT.I-SHFT ) THEN
               CALL DGEMM( 'No Transpose', 'No Transpose', K, R,
     $                     (N-I+SHFT)*L, ALPHA, TR, LDTR,
     $                     B((I-SHFT)*L+1,1), LDB, ONE, C(PC,1), LDC )
            END IF
            PC = PC + K
            PT = PT - K*L
   90    CONTINUE
C
      ELSE IF ( METH.EQ.2 .AND. LTRAN ) THEN
C
         PT = ( M - 1 - SHFT )*K + 1
C
         DO 100  I = 1, ( M - SHFT )*K*L, K*L
            CALL DLACPY( 'All', K, L, TC(PT,1), LDTC, DWORK(I), K )
            PT = PT - K
  100    CONTINUE
C
         PT = ( M - 1 )*K*L + 1
         PB = 1
C
         DO 110  I = 1, M
            CALL DGEMM( 'Tranpose', 'No Transpose', MIN( I-SHFT, N )*L,
     $                  R, K, ALPHA, DWORK(PT), K, B(PB,1), LDB, ONE,
     $                  C, LDC )
            IF ( N.GT.I-SHFT ) THEN
               CALL DGEMM( 'Transpose', 'No Transpose', (N-I+SHFT)*L, R,
     $                     K, ALPHA, TR, LDTR, B(PB,1), LDB, ONE,
     $                     C((I-SHFT)*L+1,1), LDC )
            END IF
            PB = PB + K
            PT = PT - K*L
  110    CONTINUE
C
      ELSE IF ( METH.EQ.3 ) THEN
C
C        In method 3 the matrix-vector product is computed by a suitable
C        block convolution via fast Hartley transforms similar to the
C        SLICOT routine DE01PD.
C
C        Step 1: Copy input data into the workspace arrays.
C
         PDW = 1
         IF ( LTRAN ) THEN
            DIMB = K
            DIMC = L
         ELSE
            DIMB = L
            DIMC = K
         END IF
         PB = LEN*K*L
         PC = LEN*( K*L + DIMB*R )
         IF ( LTRAN ) THEN
            IF ( FULLC ) THEN
               CALL DLACPY( 'All', K, L, TC, LDTC, DWORK, LEN*K )
            END IF
C
            DO 120  I = 1, N - 1 + SHFT
               CALL DLACPY( 'All', K, L, TR(1,(I-1)*L+1), LDTR,
     $                      DWORK((I-SHFT)*K+1), LEN*K )
  120       CONTINUE
C
            PDW = N*K + 1
            R1  = ( LEN - M - N + 1 )*K
            CALL DLASET( 'All', R1, L, ZERO, ZERO, DWORK(PDW), LEN*K )
            PDW = PDW + R1
C
            DO 130  I = ( M - 1 - SHFT )*K + 1, K - SHFT*K + 1, -K
               CALL DLACPY( 'All', K, L, TC(I,1), LDTC,
     $                      DWORK(PDW), LEN*K )
               PDW = PDW + K
  130       CONTINUE
C
            PDW = PB + 1
            CALL DLACPY( 'All', MK, R, B, LDB, DWORK(PDW), LEN*K )
            PDW = PDW + MK
            CALL DLASET( 'All', (LEN-M)*K, R, ZERO, ZERO, DWORK(PDW),
     $                   LEN*K )
         ELSE
            IF ( .NOT.FULLC ) THEN
               CALL DLACPY( 'All', K, L, TR, LDTR, DWORK, LEN*K )
            END IF
            CALL DLACPY( 'All', (M-SHFT)*K, L, TC, LDTC,
     $                   DWORK(SHFT*K+1), LEN*K )
            PDW = MK + 1
            R1  = ( LEN - M - N + 1 )*K
            CALL DLASET( 'All', R1, L, ZERO, ZERO, DWORK(PDW), LEN*K )
            PDW = PDW + R1
C
            DO 140  I = ( N - 2 + SHFT )*L + 1, SHFT*L + 1, -L
               CALL DLACPY( 'All', K, L, TR(1,I), LDTR, DWORK(PDW),
     $                      LEN*K )
               PDW = PDW + K
  140       CONTINUE
C
            PDW = PB + 1
            CALL DLACPY( 'All', NL, R, B, LDB, DWORK(PDW), LEN*L )
            PDW = PDW + NL
            CALL DLASET( 'All', (LEN-N)*L, R, ZERO, ZERO, DWORK(PDW),
     $                   LEN*L )
         END IF
C
C        Take point Toeplitz matrices into extra consideration.
C
         IF ( K*L.EQ.1 ) THEN
            WGHT = 'N'
            CALL DG01OD( 'OutputScrambled', WGHT, LEN, DWORK,
     $                   DWORK(PC+1), IERR )
C
            DO 170  I = PB, PB + LEN*R - 1, LEN
               CALL DG01OD( 'OutputScrambled', WGHT, LEN, DWORK(I+1),
     $                      DWORK(PC+1), IERR )
               SCAL = ALPHA / DBLE( LEN )
               DWORK(I+1) = SCAL*DWORK(I+1)*DWORK(1)
               DWORK(I+2) = SCAL*DWORK(I+2)*DWORK(2)
               SCAL = SCAL / TWO
C
               LN = 1
C
               DO 160 LL = 1, P - 1
                  LN = 2*LN
                  R1 = 2*LN
C
                  DO 150 P1 = LN + 1, LN + LN/2
                     T1 = DWORK(P1) + DWORK(R1)
                     T2 = DWORK(P1) - DWORK(R1)
                     TH = T2*DWORK(I+P1)
                     DWORK(I+P1) = SCAL*( T1*DWORK(I+P1)
     $                                  + T2*DWORK(I+R1) )
                     DWORK(I+R1) = SCAL*( T1*DWORK(I+R1) - TH )
                     R1 = R1 - 1
  150             CONTINUE
C
  160          CONTINUE
C
               CALL DG01OD( 'InputScrambled', WGHT, LEN, DWORK(I+1),
     $                      DWORK(PC+1), IERR )
  170       CONTINUE
C
            PC = PB
            GOTO 420
         END IF
C
C        Step 2: Compute the weights for the Hartley transforms.
C
         PDW = PC
         R1  = 1
         LN  = 1
         TH  = FOUR*ATAN( ONE ) / DBLE( LEN )
C
         DO 190  LL = 1, P - 2
            LN = 2*LN
            TH = TWO*TH
            CF = COS( TH )
            SF = SIN( TH )
            DWORK(PDW+R1)   = CF
            DWORK(PDW+R1+1) = SF
            R1 = R1 + 2
C
            DO 180  I = 1, LN-2, 2
               DWORK(PDW+R1)   = CF*DWORK(PDW+I) - SF*DWORK(PDW+I+1)
               DWORK(PDW+R1+1) = SF*DWORK(PDW+I) + CF*DWORK(PDW+I+1)
               R1 = R1 + 2
  180       CONTINUE
C
  190    CONTINUE
C
         P1 = 3
         Q1 = R1 - 2
C
         DO 210  LL = P - 2, 1, -1
C
            DO 200  I = P1, Q1, 4
               DWORK(PDW+R1)   = DWORK(PDW+I)
               DWORK(PDW+R1+1) = DWORK(PDW+I+1)
               R1 = R1 + 2
  200       CONTINUE
C
            P1 = Q1 + 4
            Q1 = R1 - 2
  210    CONTINUE
C
C        Step 3: Compute the Hartley transforms with scrambled output.
C
         J  = 0
         KK = K
C
C        WHILE   J < (L*LEN*K + R*LEN*DIMB),
C
  220    CONTINUE
C
            LN   = LEN
            WPOS = PDW+1
C
            DO 270  PP = P - 1, 1, -1
               LN = LN / 2
               P2 = 1
               Q2 = LN*KK + 1
               R2 = ( LN/2 )*KK + 1
               S2 = R2 + Q2 - 1
C
               DO 260  I = 0, LEN/( 2*LN ) - 1
C
                  DO 230  IR = 0, KK - 1
                     T1 = DWORK(Q2+IR+J)
                     DWORK(Q2+IR+J) = DWORK(P2+IR+J) - T1
                     DWORK(P2+IR+J) = DWORK(P2+IR+J) + T1
                     T1 = DWORK(S2+IR+J)
                     DWORK(S2+IR+J) = DWORK(R2+IR+J) - T1
                     DWORK(R2+IR+J) = DWORK(R2+IR+J) + T1
  230             CONTINUE
C
                  P1 = P2 + KK
                  Q1 = P1 + LN*KK
                  R1 = Q1 - 2*KK
                  S1 = R1 + LN*KK
C
                  DO 250  JJ = WPOS, WPOS + LN - 3, 2
                     CF = DWORK(JJ)
                     SF = DWORK(JJ+1)
C
                     DO 240  IR = 0, KK-1
                        T1 = DWORK(P1+IR+J) - DWORK(Q1+IR+J)
                        T2 = DWORK(R1+IR+J) - DWORK(S1+IR+J)
                        DWORK(P1+IR+J) = DWORK(P1+IR+J) +
     $                                   DWORK(Q1+IR+J)
                        DWORK(R1+IR+J) = DWORK(R1+IR+J) +
     $                                   DWORK(S1+IR+J)
                        DWORK(Q1+IR+J) =  CF*T1 + SF*T2
                        DWORK(S1+IR+J) = -CF*T2 + SF*T1
  240                CONTINUE
C
                     P1 = P1 + KK
                     Q1 = Q1 + KK
                     R1 = R1 - KK
                     S1 = S1 - KK
  250             CONTINUE
C
                  P2 = P2 + 2*KK*LN
                  Q2 = Q2 + 2*KK*LN
                  R2 = R2 + 2*KK*LN
                  S2 = S2 + 2*KK*LN
  260          CONTINUE
C
               WPOS = WPOS + LN - 2
  270       CONTINUE
C
            DO 290 ICP = KK + 1, LEN*KK, 2*KK
               ICQ = ICP - KK
C
               DO 280 IR = 0, KK - 1
                  T1 = DWORK(ICP+IR+J)
                  DWORK(ICP+IR+J) = DWORK(ICQ+IR+J) - T1
                  DWORK(ICQ+IR+J) = DWORK(ICQ+IR+J) + T1
  280          CONTINUE
C
  290       CONTINUE
C
            J = J + LEN*KK
            IF ( J.EQ.L*LEN*K ) THEN
               KK = DIMB
            END IF
         IF ( J.LT.PC )   GOTO 220
C        END WHILE 220
C
C        Step 4: Compute a Hadamard like product.
C
         CALL DCOPY( LEN-P, DWORK(PDW+1), 1,DWORK(PDW+1+R*LEN*DIMC), 1 )
         PDW  = PDW + R*LEN*DIMC
         SCAL = ALPHA / DBLE( LEN )
         P1 = 1
         R1 = LEN*K*L + 1
         S1 = R1 + LEN*DIMB*R
         IF ( LTRAN ) THEN
            KK = L
            LL = K
         ELSE
            KK = K
            LL = L
         END IF
         CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, SCAL, DWORK(P1),
     $               LEN*K, DWORK(R1), LEN*DIMB, ZERO, DWORK(S1),
     $               LEN*DIMC )
         P1 = P1 + K
         R1 = R1 + DIMB
         S1 = S1 + DIMC
         CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, SCAL, DWORK(P1),
     $               LEN*K, DWORK(R1), LEN*DIMB, ZERO, DWORK(S1),
     $               LEN*DIMC )
         SCAL = SCAL / TWO
         LN = 1
C
         DO 330 PP = 1, P - 1
            LN = 2*LN
            P2 = ( 2*LN - 1 )*K + 1
            R1 = PB + LN*DIMB + 1
            R2 = PB + ( 2*LN - 1 )*DIMB + 1
            S1 = PC + LN*DIMC + 1
            S2 = PC + ( 2*LN - 1 )*DIMC + 1
C
            DO 320  P1 = LN*K + 1, ( LN + LN/2 )*K, K
C
               DO 310 J = 0, LEN*K*( L - 1 ), LEN*K
C
                  DO 300  I = P1, P1 + K - 1
                     T1 = DWORK(P2)
                     DWORK(P2)  = DWORK(J+I) - T1
                     DWORK(J+I) = DWORK(J+I) + T1
                     P2 = P2 + 1
  300             CONTINUE
C
                  P2 = P2 + ( LEN - 1 )*K
  310          CONTINUE
C
               P2 = P2 - LEN*K*L
               CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, SCAL,
     $                     DWORK(P1), LEN*K, DWORK(R1), LEN*DIMB,
     $                     ZERO, DWORK(S1), LEN*DIMC )
               CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, SCAL,
     $                     DWORK(P2), LEN*K, DWORK(R2), LEN*DIMB, ONE,
     $                     DWORK(S1), LEN*DIMC )
               CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, SCAL,
     $                     DWORK(P1), LEN*K, DWORK(R2), LEN*DIMB, ZERO,
     $                     DWORK(S2), LEN*DIMC )
               CALL DGEMM( TRANS, 'No Transpose', KK, R, LL, -SCAL,
     $                     DWORK(P2), LEN*K, DWORK(R1), LEN*DIMB, ONE,
     $                     DWORK(S2), LEN*DIMC )
               P2 = P2 - K
               R1 = R1 + DIMB
               R2 = R2 - DIMB
               S1 = S1 + DIMC
               S2 = S2 - DIMC
  320       CONTINUE
C
  330    CONTINUE
C
C        Step 5: Hartley transform with scrambled input.
C
         DO 410 J = PC, PC + LEN*DIMC*R, LEN*DIMC
C
           DO 350 ICP = DIMC + 1, LEN*DIMC, 2*DIMC
               ICQ = ICP - DIMC
C
               DO 340 IR = 0, DIMC - 1
                  T1 = DWORK(ICP+IR+J)
                  DWORK(ICP+IR+J) = DWORK(ICQ+IR+J) - T1
                  DWORK(ICQ+IR+J) = DWORK(ICQ+IR+J) + T1
  340          CONTINUE
C
  350       CONTINUE
C
            LN   = 1
            WPOS = PDW + LEN - 2*P + 1
C
            DO 400 PP = 1, P - 1
               LN = 2*LN
               P2 = 1
               Q2 = LN*DIMC + 1
               R2 = ( LN/2 )*DIMC + 1
               S2 = R2 + Q2 - 1
C
               DO 390 I = 0, LEN/( 2*LN ) - 1
C
                  DO 360 IR = 0, DIMC - 1
                     T1 = DWORK(Q2+IR +J)
                     DWORK(Q2+IR+J) = DWORK(P2+IR+J) - T1
                     DWORK(P2+IR+J) = DWORK(P2+IR+J) + T1
                     T1 = DWORK(S2+IR+J)
                     DWORK(S2+IR+J) = DWORK(R2+IR+J) - T1
                     DWORK(R2+IR+J) = DWORK(R2+IR+J) + T1
  360             CONTINUE
C
                  P1 = P2 + DIMC
                  Q1 = P1 + LN*DIMC
                  R1 = Q1 - 2*DIMC
                  S1 = R1 + LN*DIMC
C
                  DO 380 JJ = WPOS, WPOS + LN - 3, 2
                     CF = DWORK(JJ)
                     SF = DWORK(JJ+1)
C
                     DO 370 IR = 0, DIMC - 1
                        T1 =  CF*DWORK(Q1+IR+J) + SF*DWORK(S1+IR+J)
                        T2 = -CF*DWORK(S1+IR+J) + SF*DWORK(Q1+IR+J)
                        DWORK(Q1+IR+J) = DWORK(P1+IR+J) - T1
                        DWORK(P1+IR+J) = DWORK(P1+IR+J) + T1
                        DWORK(S1+IR+J) = DWORK(R1+IR+J) - T2
                        DWORK(R1+IR+J) = DWORK(R1+IR+J) + T2
  370                CONTINUE
C
                     P1 = P1 + DIMC
                     Q1 = Q1 + DIMC
                     R1 = R1 - DIMC
                     S1 = S1 - DIMC
  380             CONTINUE
C
                  P2 = P2 + 2*DIMC*LN
                  Q2 = Q2 + 2*DIMC*LN
                  R2 = R2 + 2*DIMC*LN
                  S2 = S2 + 2*DIMC*LN
  390          CONTINUE
C
               WPOS = WPOS - 2*LN + 2
  400       CONTINUE
C
  410    CONTINUE
C
C        Step 6: Copy data from workspace to output.
C
  420    CONTINUE
C
         IF ( LTRAN ) THEN
            I = NL
         ELSE
            I = MK
         END IF
C
         DO 430  J = 0, R - 1
            CALL DAXPY( I, ONE, DWORK(PC+(J*LEN*DIMC) + 1), 1,
     $                  C(1,J+1), 1 )
  430    CONTINUE
C
      END IF
      DWORK(1) = DBLE( MAX( 1, WRKOPT ) )
      RETURN
C
C *** Last line of MB02KD ***
      END
