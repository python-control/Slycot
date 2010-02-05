      SUBROUTINE MB02GD( TYPET, TRIU, K, N, NL, P, S, T, LDT, RB, LDRB,
     $                   DWORK, LDWORK, INFO )
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
C     To compute the Cholesky factor of a banded symmetric positive
C     definite (s.p.d.) block Toeplitz matrix, defined by either its
C     first block row, or its first block column, depending on the
C     routine parameter TYPET.
C
C     By subsequent calls of this routine the Cholesky factor can be
C     computed block column by block column.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPET   CHARACTER*1
C             Specifies the type of T, as follows:
C             = 'R':  T contains the first block row of an s.p.d. block
C                     Toeplitz matrix; the Cholesky factor is upper
C                     triangular;
C             = 'C':  T contains the first block column of an s.p.d.
C                     block Toeplitz matrix; the Cholesky factor is
C                     lower triangular. This choice results in a column
C                     oriented algorithm which is usually faster.
C             Note:   in the sequel, the notation x / y means that
C                     x corresponds to TYPET = 'R' and y corresponds to
C                     TYPET = 'C'.
C
C     TRIU    CHARACTER*1
C             Specifies the structure of the last block in T, as
C             follows:
C             = 'N':  the last block has no special structure;
C             = 'T':  the last block is lower / upper triangular.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of rows / columns in T, which should be equal
C             to the blocksize.  K >= 0.
C
C     N       (input)  INTEGER
C             The number of blocks in T.  N >= 1.
C             If TRIU = 'N',   N >= 1;
C             if TRIU = 'T',   N >= 2.
C
C     NL      (input)  INTEGER
C             The lower block bandwidth, i.e., NL + 1 is the number of
C             nonzero blocks in the first block column of the block
C             Toeplitz matrix.
C             If TRIU = 'N',   0 <= NL < N;
C             if TRIU = 'T',   1 <= NL < N.
C
C     P       (input)  INTEGER
C             The number of previously computed block rows / columns of
C             the Cholesky factor.  0 <= P <= N.
C
C     S       (input)  INTEGER
C             The number of block rows / columns of the Cholesky factor
C             to compute.  0 <= S <= N - P.
C
C     T       (input/output)  DOUBLE PRECISION array, dimension
C             (LDT,(NL+1)*K) / (LDT,K)
C             On entry, if P = 0, the leading K-by-(NL+1)*K /
C             (NL+1)*K-by-K part of this array must contain the first
C             block row / column of an s.p.d. block Toeplitz matrix.
C             On entry, if P > 0, the leading K-by-(NL+1)*K /
C             (NL+1)*K-by-K part of this array must contain the P-th
C             block row / column of the Cholesky factor.
C             On exit, if INFO = 0, then the leading K-by-(NL+1)*K /
C             (NL+1)*K-by-K part of this array contains the (P+S)-th
C             block row / column of the Cholesky factor.
C
C     LDT     INTEGER
C             The leading dimension of the array T.
C             LDT >= MAX(1,K) / MAX(1,(NL+1)*K).
C
C     RB      (input/output)  DOUBLE PRECISION array, dimension
C             (LDRB,MIN(P+NL+S,N)*K) / (LDRB,MIN(P+S,N)*K)
C             On entry, if TYPET = 'R'  and  TRIU = 'N'  and  P > 0,
C             the leading (NL+1)*K-by-MIN(NL,N-P)*K part of this array
C             must contain the (P*K+1)-st to ((P+NL)*K)-th columns
C             of the upper Cholesky factor in banded format from a
C             previous call of this routine.
C             On entry, if TYPET = 'R'  and  TRIU = 'T'  and  P > 0,
C             the leading (NL*K+1)-by-MIN(NL,N-P)*K part of this array
C             must contain the (P*K+1)-st to (MIN(P+NL,N)*K)-th columns
C             of the upper Cholesky factor in banded format from a
C             previous call of this routine.
C             On exit, if TYPET = 'R'  and  TRIU = 'N', the leading
C             (NL+1)*K-by-MIN(NL+S,N-P)*K part of this array contains
C             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the
C             upper Cholesky factor in banded format.
C             On exit, if TYPET = 'R'  and  TRIU = 'T', the leading
C             (NL*K+1)-by-MIN(NL+S,N-P)*K part of this array contains
C             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the
C             upper Cholesky factor in banded format.
C             On exit, if TYPET = 'C'  and  TRIU = 'N', the leading
C             (NL+1)*K-by-MIN(S,N-P)*K part of this array contains
C             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower
C             Cholesky factor in banded format.
C             On exit, if TYPET = 'C'  and  TRIU = 'T', the leading
C             (NL*K+1)-by-MIN(S,N-P)*K part of this array contains
C             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower
C             Cholesky factor in banded format.
C             For further details regarding the band storage scheme see
C             the documentation of the LAPACK routine DPBTF2.
C
C     LDRB    INTEGER
C             The leading dimension of the array RB.
C             If TRIU = 'N',   LDRB >= MAX( (NL+1)*K,1 );
C             if TRIU = 'T',   LDRB >= NL*K+1.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -13,  DWORK(1)  returns the minimum
C             value of LDWORK.
C             The first 1 + ( NL + 1 )*K*K elements of DWORK should be
C             preserved during successive calls of the routine.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 1 + ( NL + 1 )*K*K + NL*K.
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
C                                3
C     The algorithm requires O( K *N*NL ) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
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
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TRIU, TYPET
      INTEGER           INFO, K, LDRB, LDT, LDWORK, N, NL, P, S
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), RB(LDRB,*), T(LDT,*)
C     .. Local Scalars ..
      CHARACTER         STRUCT
      LOGICAL           ISROW, LTRI
      INTEGER           HEAD, I, IERR, J, JJ, KK, LEN, LEN2, LENR, NB,
     $                  NBMIN, PDW, POSR, PRE, RNK, SIZR, STPS, WRKMIN,
     $                  WRKOPT
C     .. Local Arrays ..
      INTEGER           IPVT(1)
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLACPY, DLASET, DPOTRF, DTRSM, MB02CU,
     $                  MB02CV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, MOD
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO = 0
      LTRI = LSAME( TRIU, 'T' )
      LENR = ( NL + 1 )*K
      IF ( LTRI ) THEN
         SIZR = NL*K + 1
      ELSE
         SIZR = LENR
      END IF
      ISROW  = LSAME( TYPET, 'R' )
      WRKMIN = 1 + ( LENR + NL )*K
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( ISROW .OR. LSAME( TYPET, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( LTRI .OR. LSAME( TRIU, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( K.LT.0 ) THEN
         INFO = -3
      ELSE IF ( (      LTRI .AND. N.LT.2 ) .OR.
     $          ( .NOT.LTRI .AND. N.LT.1 ) ) THEN
         INFO = -4
      ELSE IF ( NL.GE.N .OR. (      LTRI .AND. NL.LT.1 ) .OR.
     $                       ( .NOT.LTRI .AND. NL.LT.0 ) ) THEN
         INFO = -5
      ELSE IF ( P.LT.0 .OR. P.GT.N ) THEN
         INFO = -6
      ELSE IF ( S.LT.0 .OR. S.GT.N-P ) THEN
         INFO = -7
      ELSE IF ( (      ISROW .AND. LDT.LT.MAX( 1, K ) ) .OR.
     $          ( .NOT.ISROW .AND. LDT.LT.MAX( 1, LENR ) ) )
     $      THEN
         INFO = -9
      ELSE IF ( (      LTRI .AND. LDRB.LT.SIZR ) .OR.
     $          ( .NOT.LTRI .AND. LDRB.LT.MAX( 1, LENR ) ) )
     $      THEN
         INFO = -11
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -13
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02GD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( S*K.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Compute the generator if P = 0.
C
      IF ( P.EQ.0 ) THEN
         IF ( ISROW ) THEN
            CALL DPOTRF( 'Upper', K, T, LDT, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
            IF ( NL.GT.0 )
     $         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'NonUnit', K,
     $                     NL*K, ONE, T, LDT, T(1,K+1), LDT )
C
C           Copy the first block row to RB.
C
            IF ( LTRI ) THEN
C
               DO 10  I = 1, LENR - K
                  CALL DCOPY( MIN( I, K ), T(1,I), 1,
     $                        RB( MAX( SIZR-I+1, 1 ),I ), 1 )
   10          CONTINUE
C
               DO 20  I = K, 1, -1
                  CALL DCOPY( I, T(K-I+1,LENR-I+1), 1,
     $                        RB( 1,LENR-I+1 ), 1 )
   20          CONTINUE
C
            ELSE
C
               DO 30  I = 1, LENR
                  CALL DCOPY( MIN( I, K ), T(1,I), 1,
     $                        RB( MAX( SIZR-I+1, 1 ),I ), 1 )
   30          CONTINUE
C
            END IF
C
C           Quick return if N = 1.
C
            IF ( N.EQ.1 ) THEN
               DWORK(1) = ONE
               RETURN
            END IF
C
            CALL DLACPY( 'All', K, NL*K, T(1,K+1), LDT, DWORK(2), K )
            CALL DLASET( 'All', K, K, ZERO, ZERO, DWORK(NL*K*K+2), K )
            POSR = K + 1
         ELSE
            CALL DPOTRF( 'Lower', K, T, LDT, IERR )
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
            IF ( NL.GT.0 )
     $         CALL DTRSM( 'Right', 'Lower', 'Transpose', 'NonUnit',
     $                     NL*K, K, ONE, T, LDT, T(K+1,1), LDT )
C
C           Copy the first block column to RB.
C
            POSR = 1
            IF ( LTRI ) THEN
C
               DO 40  I = 1, K
                  CALL DCOPY( SIZR, T(I,I), 1, RB(1,POSR), 1 )
                  POSR = POSR + 1
   40          CONTINUE
C
            ELSE
C
               DO 50  I = 1, K
                  CALL DCOPY(  LENR-I+1, T(I,I), 1, RB(1,POSR), 1 )
                  IF ( LENR.LT.N*K .AND. I.GT.1 ) THEN
                     CALL DLASET( 'All', I-1, 1, ZERO, ZERO,
     $                            RB(LENR-I+2,POSR), LDRB )
                  END IF
                  POSR = POSR + 1
   50          CONTINUE
C
            END IF
C
C           Quick return if N = 1.
C
            IF ( N.EQ.1 ) THEN
               DWORK(1) = ONE
               RETURN
            END IF
C
            CALL DLACPY( 'All', NL*K, K, T(K+1,1), LDT, DWORK(2), LENR )
            CALL DLASET( 'All', K, K, ZERO, ZERO, DWORK(NL*K+2), LENR )
         END IF
         PRE  = 1
         STPS = S - 1
      ELSE
         PRE  = P
         STPS = S
         POSR = 1
      END IF
C
      PDW  = LENR*K + 1
      HEAD = MOD( ( PRE - 1 )*K, LENR )
C
C     Determine block size for the involved block Householder
C     transformations.
C
      IF ( ISROW ) THEN
         NB = MIN( ILAENV( 1, 'DGEQRF', ' ', K, LENR, -1, -1 ), K )
      ELSE
         NB = MIN( ILAENV( 1, 'DGELQF', ' ', LENR, K, -1, -1 ), K )
      END IF
      KK = PDW + 4*K
      WRKOPT = KK + LENR*NB
      KK = LDWORK - KK
      IF ( KK.LT.LENR*NB )  NB = KK / LENR
      IF ( ISROW ) THEN
         NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', K, LENR, -1, -1 ) )
      ELSE
         NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', LENR, K, -1, -1 ) )
      END IF
      IF ( NB.LT.NBMIN )  NB = 0
C
C     Generator reduction process.
C
      IF ( ISROW ) THEN
C
         DO 90  I = PRE, PRE + STPS - 1
            CALL MB02CU( 'Row', K, K, K, NB, T, LDT, DUM, 1,
     $                   DWORK(HEAD*K+2), K, RNK, IPVT, DWORK(PDW+1),
     $                   ZERO, DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The positive definiteness is (numerically)
C                             not satisfied.
C
               INFO = 1
               RETURN
            END IF
C
            LEN  = MAX( MIN( ( N - I )*K - K, LENR - HEAD - K ), 0 )
            LEN2 = MAX( MIN( ( N - I )*K - LEN - K, HEAD ), 0 )
            IF ( LEN.EQ.( LENR-K ) ) THEN
               STRUCT = TRIU
            ELSE
               STRUCT = 'N'
            END IF
            CALL MB02CV( 'Row', STRUCT, K, LEN, K, K, NB, -1, DUM, 1,
     $                   DUM, 1, DWORK(HEAD*K+2), K, T(1,K+1), LDT,
     $                   DUM, 1, DWORK((HEAD+K)*K+2), K, DWORK(PDW+1),
     $                   DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            IF ( ( N - I )*K.GE.LENR ) THEN
               STRUCT = TRIU
            ELSE
               STRUCT = 'N'
            END IF
            CALL MB02CV( 'Row', STRUCT, K, LEN2, K, K, NB, -1, DUM, 1,
     $                   DUM, 1, DWORK(HEAD*K+2), K, T(1,K+LEN+1), LDT,
     $                   DUM, 1, DWORK(2), K, DWORK(PDW+1),
     $                   DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            CALL DLASET( 'All', K, K, ZERO, ZERO, DWORK(HEAD*K+2), K )
C
C           Copy current block row to RB.
C
            IF ( LTRI ) THEN
C
               DO 60  J = 1, MIN( LEN + LEN2 + K, LENR - K )
                  CALL DCOPY(  MIN( J, K ), T(1,J), 1,
     $                         RB(MAX( SIZR-J+1, 1 ),POSR+J-1 ), 1 )
   60          CONTINUE
C
               IF ( LEN+LEN2+K.GE.LENR ) THEN
C
                  DO 70  JJ = K, 1, -1
                     CALL DCOPY(  JJ, T(K-JJ+1,LENR-JJ+1), 1,
     $                            RB(1,POSR+LENR-JJ), 1 )
   70             CONTINUE
C
               END IF
               POSR = POSR + K
C
            ELSE
C
               DO 80  J = 1, LEN + LEN2 + K
                  CALL DCOPY(  MIN( J, K ), T(1,J), 1,
     $                         RB(MAX( SIZR-J+1, 1 ),POSR+J-1), 1 )
                  IF ( J.GT.LENR-K ) THEN
                     CALL DLASET( 'All', SIZR-J, 1, ZERO, ZERO,
     $                            RB(1,POSR+J-1), 1 )
                  END IF
   80          CONTINUE
C
               POSR = POSR + K
            END IF
            HEAD = MOD( HEAD + K, LENR )
   90    CONTINUE
C
      ELSE
C
         DO 120  I = PRE, PRE + STPS - 1
C
            CALL MB02CU( 'Column', K, K, K, NB, T, LDT, DUM, 1,
     $                   DWORK(HEAD+2), LENR, RNK, IPVT, DWORK(PDW+1),
     $                   ZERO, DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            IF ( IERR.NE.0 )  THEN
C
C              Error return:  The positive definiteness is (numerically)
C                             not satisfied.
C
               INFO = 1
               RETURN
            END IF
C
            LEN  = MAX( MIN( ( N - I )*K - K, LENR - HEAD - K ), 0 )
            LEN2 = MAX( MIN( ( N - I )*K - LEN - K, HEAD ), 0 )
            IF ( LEN.EQ.( LENR-K ) ) THEN
               STRUCT = TRIU
            ELSE
               STRUCT = 'N'
            END IF
            CALL MB02CV( 'Column', STRUCT, K, LEN, K, K, NB, -1, DUM,
     $                   1, DUM, 1, DWORK(HEAD+2), LENR, T(K+1,1), LDT,
     $                   DUM, 1, DWORK(HEAD+K+2), LENR, DWORK(PDW+1),
     $                   DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            IF ( ( N - I )*K.GE.LENR ) THEN
               STRUCT = TRIU
            ELSE
               STRUCT = 'N'
            END IF
            CALL MB02CV( 'Column', STRUCT, K, LEN2, K, K, NB, -1, DUM,
     $                   1, DUM, 1, DWORK(HEAD+2), LENR, T(K+LEN+1,1),
     $                   LDT, DUM, 1, DWORK(2), LENR, DWORK(PDW+1),
     $                   DWORK(PDW+4*K+1), LDWORK-PDW-4*K, IERR )
C
            CALL DLASET( 'All', K, K, ZERO, ZERO, DWORK(HEAD+2), LENR )
C
C           Copy current block column to RB.
C
            IF ( LTRI ) THEN
C
               DO 100  J = 1, K
                  CALL DCOPY( MIN( SIZR, (N-I)*K-J+1 ), T(J,J), 1,
     $                        RB(1,POSR), 1 )
                  POSR = POSR + 1
  100          CONTINUE
C
            ELSE
C
               DO 110  J = 1, K
                  CALL DCOPY( MIN( SIZR-J+1, (N-I)*K-J+1 ), T(J,J), 1,
     $                        RB(1,POSR), 1 )
                  IF ( LENR.LT.(N-I)*K ) THEN
                     CALL DLASET( 'All', J-1, 1, ZERO, ZERO,
     $                             RB(MIN( SIZR-J+1, (N-I)*K-J+1 ) + 1,
     $                                POSR), LDRB )
                  END IF
                  POSR = POSR + 1
  110          CONTINUE
C
            END IF
            HEAD = MOD( HEAD + K, LENR )
  120    CONTINUE
C
      END IF
      DWORK(1) = DBLE( WRKOPT )
      RETURN
C
C *** Last line of MB02GD ***
      END
