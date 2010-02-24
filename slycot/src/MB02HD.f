      SUBROUTINE MB02HD( TRIU, K, L, M, ML, N, NU, P, S, TC, LDTC, TR,
     $                   LDTR, RB, LDRB, DWORK, LDWORK, INFO )
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
C     To compute, for a banded K*M-by-L*N block Toeplitz matrix T with
C     block size (K,L), specified by the nonzero blocks of its first
C     block column TC and row TR, a LOWER triangular matrix R (in band
C     storage scheme) such that
C                          T          T
C                         T  T  =  R R .                             (1)
C
C     It is assumed that the first MIN(M*K, N*L) columns of T are
C     linearly independent.
C
C     By subsequent calls of this routine, the matrix R can be computed
C     block column by block column.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRIU    CHARACTER*1
C             Specifies the structure, if any, of the last blocks in TC
C             and TR, as follows:
C             = 'N':  TC and TR have no special structure;
C             = 'T':  TC and TR are upper and lower triangular,
C                     respectively. Depending on the block sizes, two
C                     different shapes of the last blocks in TC and TR
C                     are possible, as illustrated below:
C
C                     1)    TC       TR     2)   TC         TR
C
C                          x x x    x 0 0      x x x x    x 0 0 0
C                          0 x x    x x 0      0 x x x    x x 0 0
C                          0 0 x    x x x      0 0 x x    x x x 0
C                          0 0 0    x x x
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
C             M >= 1.
C
C     ML      (input) INTEGER
C             The lower block bandwidth, i.e., ML + 1 is the number of
C             nonzero blocks in the first block column of T.
C             0 <= ML < M and (ML + 1)*K >= L and
C             if ( M*K <= N*L ),  ML >= M - INT( ( M*K - 1 )/L ) - 1;
C                                 ML >= M - INT( M*K/L ) or
C                                 MOD( M*K, L ) >= K;
C             if ( M*K >= N*L ),  ML*K >= N*( L - K ).
C
C     N       (input) INTEGER
C             The number of blocks in the first block row of T.
C             N >= 1.
C
C     NU      (input) INTEGER
C             The upper block bandwidth, i.e., NU + 1 is the number of
C             nonzero blocks in the first block row of T.
C             If TRIU = 'N',   0 <= NU < N and
C                              (M + NU)*L >= MIN( M*K, N*L );
C             if TRIU = 'T',   MAX(1-ML,0) <= NU < N and
C                              (M + NU)*L >= MIN( M*K, N*L ).
C
C     P       (input)  INTEGER
C             The number of previously computed block columns of R.
C             P*L < MIN( M*K,N*L ) + L and P >= 0.
C
C     S       (input)  INTEGER
C             The number of block columns of R to compute.
C             (P+S)*L < MIN( M*K,N*L ) + L and S >= 0.
C
C     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L)
C             On entry, if P = 0, the leading (ML+1)*K-by-L part of this
C             array must contain the nonzero blocks in the first block
C             column of T.
C
C     LDTC    INTEGER
C             The leading dimension of the array TC.
C             LDTC >= MAX(1,(ML+1)*K),  if P = 0.
C
C     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,NU*L)
C             On entry, if P = 0, the leading K-by-NU*L part of this
C             array must contain the 2nd to the (NU+1)-st blocks of
C             the first block row of T.
C
C     LDTR    INTEGER
C             The leading dimension of the array TR.
C             LDTR >= MAX(1,K),  if P = 0.
C
C     RB      (output)  DOUBLE PRECISION array, dimension
C             (LDRB,MIN( S*L,MIN( M*K,N*L )-P*L ))
C             On exit, if INFO = 0 and TRIU = 'N', the leading
C             MIN( ML+NU+1,N )*L-by-MIN( S*L,MIN( M*K,N*L )-P*L ) part
C             of this array contains the (P+1)-th to (P+S)-th block
C             column of the lower R factor (1) in band storage format.
C             On exit, if INFO = 0 and TRIU = 'T', the leading
C             MIN( (ML+NU)*L+1,N*L )-by-MIN( S*L,MIN( M*K,N*L )-P*L )
C             part of this array contains the (P+1)-th to (P+S)-th block
C             column of the lower R factor (1) in band storage format.
C             For further details regarding the band storage scheme see
C             the documentation of the LAPACK routine DPBTF2.
C
C     LDRB    INTEGER
C             The leading dimension of the array RB.
C             LDRB >= MAX( MIN( ML+NU+1,N )*L,1 ),      if TRIU = 'N';
C             LDRB >= MAX( MIN( (ML+NU)*L+1,N*L ),1 ),  if TRIU = 'T'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -17,  DWORK(1)  returns the minimum
C             value of LDWORK.
C             The first 1 + 2*MIN( ML+NU+1,N )*L*(K+L) elements of DWORK
C             should be preserved during successive calls of the routine.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             Let x = MIN( ML+NU+1,N ), then
C             LDWORK >= 1 + MAX( x*L*L + (2*NU+1)*L*K,
C                                2*x*L*(K+L) + (6+x)*L ),  if P = 0;
C             LDWORK >= 1 + 2*x*L*(K+L) + (6+x)*L,         if P > 0.
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
C     The implemented method yields a factor R which has comparable
C     accuracy with the Cholesky factor of T^T * T.
C     The algorithm requires
C               2                                  2
C           O( L *K*N*( ML + NU ) + N*( ML + NU )*L *( L + K ) )
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
      CHARACTER         TRIU
      INTEGER           INFO, K, L, LDRB, LDTC, LDTR, LDWORK, M, ML, N,
     $                  NU, P, S
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), RB(LDRB,*), TC(LDTC,*),
     $                  TR(LDTR,*)
C     .. Local Scalars ..
      CHARACTER         STRUCT
      INTEGER           COL2, HEAD, I, IERR, J, KK, LEN, LEN2, LENC,
     $                  LENL, LENR, NB, NBMIN, PDC, PDR, PDW, PFR, PNR,
     $                  POSR, PRE, PT, RNK, SIZR, STPS, WRKMIN, WRKOPT,
     $                  X
      LOGICAL           LTRI
C     .. Local Arrays ..
      INTEGER           IPVT(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEQRF, DLACPY, DLASET, DORGQR,
     $                  MA02AD, MB02CU, MB02CV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, MOD
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO = 0
      LTRI = LSAME( TRIU, 'T' )
      X    = MIN( ML + NU + 1, N )
      LENR = X*L
      IF ( LTRI ) THEN
         SIZR = MIN( ( ML + NU )*L + 1, N*L )
      ELSE
         SIZR = LENR
      END IF
      IF ( P.EQ.0 ) THEN
         WRKMIN = 1 + MAX( LENR*L + ( 2*NU + 1 )*L*K,
     $                     2*LENR*( K + L ) + ( 6 + X )*L )
      ELSE
         WRKMIN = 1 + 2*LENR*( K + L ) + ( 6 + X )*L
      END IF
      POSR = 1
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( LTRI .OR. LSAME( TRIU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF ( K.LT.0 ) THEN
         INFO = -2
      ELSE IF ( L.LT.0 ) THEN
         INFO = -3
      ELSE IF ( M.LT.1 ) THEN
         INFO = -4
      ELSE IF ( ML.GE.M .OR. ( ML + 1 )*K.LT.L .OR. ( M*K.LE.N*L .AND.
     $          ( ( ML.LT.M - INT( ( M*K - 1 )/L ) - 1 ) .OR.
     $            ( ML.LT.M - INT( M*K/L ).AND.MOD( M*K, L ).LT.K ) ) )
     $          .OR. ( M*K.GE.N*L .AND. ML*K.LT.N*( L - K ) ) ) THEN
         INFO = -5
      ELSE IF ( N.LT.1 ) THEN
         INFO = -6
      ELSE IF ( NU.GE.N .OR. NU.LT.0 .OR. ( LTRI .AND. NU.LT.1-ML ) .OR.
     $          (M + NU)*L.LT.MIN( M*K, N*L ) ) THEN
         INFO = -7
      ELSE IF ( P.LT.0 .OR. ( P*L - L ).GE.MIN( M*K, N*L ) ) THEN
         INFO = -8
      ELSE IF ( S.LT.0 .OR. ( P + S - 1 )*L.GE.MIN( M*K, N*L ) ) THEN
         INFO = -9
      ELSE IF ( P.EQ.0 .AND. LDTC.LT.MAX( 1, ( ML + 1 )*K ) ) THEN
         INFO = -11
      ELSE IF ( P.EQ.0 .AND. LDTR.LT.MAX( 1, K ) ) THEN
         INFO = -13
      ELSE IF ( LDRB.LT.MAX( SIZR, 1 ) ) THEN
         INFO = 15
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -17
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02HD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( L*K*S.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      WRKOPT = 1
C
C     Compute the generator if P = 0.
C
      IF ( P.EQ.0 ) THEN
C
C        1st column of the generator.
C
         LENC = ( ML + 1 )*K
         LENL = MAX( ML + 1 + MIN( NU, N - M ), 0 )
         PDC  = LENR*L + 1
         PDW  = PDC + LENC*L
C
C        QR decomposition of the nonzero blocks in TC.
C
         CALL DLACPY( 'All', LENC, L, TC, LDTC, DWORK(PDC+1), LENC )
         CALL DGEQRF( LENC, L, DWORK(PDC+1), LENC, DWORK(PDW+1),
     $                DWORK(PDW+L+1), LDWORK-PDW-L, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+L+1) ) + PDW + L )
C
C        The R factor is the transposed of the first block in the
C        generator.
C
         CALL MA02AD( 'Upper part', L, L, DWORK(PDC+1), LENC, DWORK(2),
     $                LENR )
C
C        Get the first block column of the Q factor.
C
         CALL DORGQR( LENC, L, L, DWORK(PDC+1), LENC, DWORK(PDW+1),
     $                DWORK(PDW+L+1), LDWORK-PDW-L, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+L+1) ) + PDW + L )
C
C        Construct a flipped copy of TC for faster multiplication.
C
         PT = LENC - 2*K + 1
C
         DO 10  I = PDW + 1, PDW + ML*K*L, K*L
            CALL DLACPY( 'All', K, L, TC(PT,1), LDTC, DWORK(I), K )
            PT = PT - K
   10    CONTINUE
C
C        Multiply T^T with the first block column of Q.
C
         PDW = I
         PDR = L + 2
         LEN = NU*L
         CALL DLASET( 'All', LENR-L, L, ZERO, ZERO, DWORK(PDR), LENR )
C
         DO 20  I = 1, ML + 1
            CALL DGEMM( 'Transpose', 'NonTranspose', MIN( I-1, N-1 )*L,
     $                  L, K, ONE, DWORK(PDW), K, DWORK(PDC+1), LENC,
     $                  ONE, DWORK(PDR), LENR )
            IF ( LEN.GT.0 ) THEN
               CALL DGEMM( 'Transpose', 'NonTranspose', LEN, L, K, ONE,
     $                     TR, LDTR, DWORK(PDC+1), LENC, ONE,
     $                     DWORK(PDR+(I-1)*L), LENR )
            END IF
            PDW = PDW - K*L
            PDC = PDC + K
            IF ( I.GE.N-NU )  LEN = LEN - L
   20    CONTINUE
C
C        Copy the first block column to R.
C
         IF ( LTRI ) THEN
C
            DO 30  I = 1, L
               CALL DCOPY( MIN( SIZR, N*L - I + 1 ),
     $                     DWORK(( I - 1 )*LENR + I + 1), 1, RB(1,POSR),
     $                     1 )
               POSR = POSR + 1
   30       CONTINUE
C
         ELSE
C
            DO 40  I = 1, L
               CALL DCOPY(  LENR-I+1, DWORK(( I - 1 )*LENR + I + 1), 1,
     $                      RB(1,POSR), 1 )
               IF ( LENR.LT.N*L .AND. I.GT.1 ) THEN
                  CALL DLASET( 'All', I-1, 1, ZERO, ZERO,
     $                         RB(LENR-I+2,POSR), LDRB )
               END IF
               POSR = POSR + 1
   40       CONTINUE
C
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
         PDR = LENR*L + 1
         CALL MA02AD( 'All', K, NU*L, TR, LDTR, DWORK(PDR+1), LENR )
         CALL DLASET( 'All', LENR-NU*L, K, ZERO, ZERO,
     $                DWORK(PDR+NU*L+1), LENR )
C
C        3rd column of the generator.
C
         PNR = PDR + LENR*K
         CALL DLACPY( 'All', LENR-L, L, DWORK(L+2), LENR, DWORK(PNR+1),
     $                LENR )
         CALL DLASET( 'All', L, L, ZERO, ZERO, DWORK(PNR+LENR-L+1),
     $                LENR )
C
C        4th column of the generator.
C
         PFR = PNR + LENR*L
C
         PDW = PFR + MOD( ( M - ML - 1 )*L, LENR )
         PT  = ML*K + 1
         DO 50  I = 1, MIN( ML + 1, LENL )
            CALL MA02AD( 'All', K, L, TC(PT,1), LDTC, DWORK(PDW+1),
     $                   LENR )
            PT  = PT - K
            PDW = PFR + MOD( PDW + L - PFR, LENR )
   50    CONTINUE
         PT = 1
         DO 60  I = ML + 2, LENL
            CALL MA02AD( 'All', K, L, TR(1,PT), LDTR, DWORK(PDW+1),
     $                   LENR )
            PT  = PT + L
            PDW = PFR + MOD( PDW + L - PFR, LENR )
   60    CONTINUE
         PRE  = 1
         STPS = S - 1
      ELSE
         PDR  = LENR*L + 1
         PNR  = PDR + LENR*K
         PFR  = PNR + LENR*L
         PRE  = P
         STPS = S
      END IF
C
      PDW  = PFR + LENR*K
      HEAD = MOD( ( PRE - 1 )*L, LENR )
C
C     Determine block size for the involved block Householder
C     transformations.
C
      NB = MIN( ILAENV( 1, 'DGELQF', ' ', LENR, L, -1, -1 ), L )
      KK = PDW + 6*L
      WRKOPT = MAX( WRKOPT, KK + LENR*NB )
      KK = LDWORK - KK
      IF ( KK.LT.LENR*NB )  NB = KK / LENR
      NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', LENR, L, -1, -1 ) )
      IF ( NB.LT.NBMIN )  NB = 0
C
C     Generator reduction process.
C
      DO 90  I = PRE, PRE + STPS - 1
C
C        The 4th generator column is not used in the first (M-ML) steps.
C
         IF ( I.LT.M-ML ) THEN
            COL2 = L
         ELSE
            COL2 = K + L
         END IF
C
         KK = MIN( L, M*K - I*L )
         CALL MB02CU( 'Column', KK, KK+K, COL2, NB, DWORK(2), LENR,
     $                DWORK(PDR+HEAD+1), LENR, DWORK(PNR+HEAD+1), LENR,
     $                RNK, IPVT, DWORK(PDW+1), ZERO, DWORK(PDW+6*L+1),
     $                LDWORK-PDW-6*L, IERR )
         IF ( IERR.NE.0 )  THEN
C
C           Error return:  The rank condition is (numerically) not
C                          satisfied.
C
            INFO = 1
            RETURN
         END IF
C
         LEN  = MAX( MIN( ( N - I )*L - KK, LENR - HEAD - KK ), 0 )
         LEN2 = MAX( MIN( ( N - I )*L - LEN - KK, HEAD ), 0 )
         IF ( LEN.EQ.( LENR - KK ) ) THEN
            STRUCT = TRIU
         ELSE
            STRUCT = 'N'
         END IF
         CALL MB02CV( 'Column', STRUCT, KK, LEN, KK+K, COL2, NB, -1,
     $                DWORK(2), LENR, DWORK(PDR+HEAD+1), LENR,
     $                DWORK(PNR+HEAD+1), LENR, DWORK(KK+2), LENR,
     $                DWORK(PDR+HEAD+KK+1), LENR, DWORK(PNR+HEAD+KK+1),
     $                LENR, DWORK(PDW+1), DWORK(PDW+6*L+1),
     $                LDWORK-PDW-6*L, IERR )
C
         IF ( ( N - I )*L.GE.LENR ) THEN
            STRUCT = TRIU
         ELSE
            STRUCT = 'N'
         END IF
C
         CALL MB02CV( 'Column', STRUCT, KK, LEN2, KK+K, COL2, NB, -1,
     $                DWORK(2), LENR, DWORK(PDR+HEAD+1), LENR,
     $                DWORK(PNR+HEAD+1), LENR, DWORK(KK+LEN+2), LENR,
     $                DWORK(PDR+1), LENR, DWORK(PNR+1), LENR,
     $                DWORK(PDW+1), DWORK(PDW+6*L+1),
     $                LDWORK-PDW-6*L, IERR )
C
         CALL DLASET( 'All', L, K+COL2, ZERO, ZERO, DWORK(PDR+HEAD+1),
     $                LENR )
C
C        Copy current block column to R.
C
         IF ( LTRI ) THEN
C
            DO 70  J = 1, KK
               CALL DCOPY( MIN( SIZR, (N-I)*L-J+1 ),
     $                      DWORK(( J - 1 )*LENR + J + 1), 1,
     $                      RB(1,POSR), 1 )
               POSR = POSR + 1
   70       CONTINUE
C
         ELSE
C
            DO 80  J = 1, KK
               CALL DCOPY(  MIN( SIZR-J+1, (N-I)*L-J+1 ),
     $                      DWORK(( J - 1 )*LENR + J + 1), 1,
     $                      RB(1,POSR), 1 )
               IF ( LENR.LT.( N - I )*L .AND. J.GT.1 ) THEN
                  CALL DLASET( 'All', J-1, 1, ZERO, ZERO,
     $                         RB(MIN( SIZR-J+1, (N-I)*L-J+1 )+1,POSR),
     $                         LDRB )
               END IF
               POSR = POSR + 1
   80       CONTINUE
C
         END IF
C
         HEAD = MOD( HEAD + L, LENR )
   90 CONTINUE
C
      DWORK(1) = DBLE( WRKOPT )
      RETURN
C
C *** Last line of MB02HD ***
      END
