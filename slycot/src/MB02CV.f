      SUBROUTINE MB02CV( TYPEG, STRUCG, K, N, P, Q, NB, RNK, A1, LDA1,
     $                   A2, LDA2, B, LDB, F1, LDF1, F2, LDF2, G, LDG,
     $                   CS, DWORK, LDWORK, INFO )
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
C     To apply the transformations created by the SLICOT Library routine
C     MB02CU on other columns / rows of the generator, contained in the
C     arrays F1, F2 and G.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPEG   CHARACTER*1
C             Specifies the type of the generator, as follows:
C             = 'D':  generator is column oriented and rank
C                     deficient;
C             = 'C':  generator is column oriented and not rank
C                     deficient;
C             = 'R':  generator is row oriented and not rank
C                     deficient.
C             Note that this parameter must be equivalent with the
C             used TYPEG in the call of MB02CU.
C
C     STRUCG  CHARACTER*1
C             Information about the structure of the generators,
C             as follows:
C             = 'T':  the trailing block of the positive generator
C                     is upper / lower triangular, and the trailing
C                     block of the negative generator is zero;
C             = 'N':  no special structure to mention.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of rows in A1 to be processed.  K >= 0.
C
C     N       (input)  INTEGER
C             If TYPEG = 'D'  or  TYPEG = 'C', the number of rows in F1;
C             if TYPEG = 'R', the number of columns in F1.  N >= 0.
C
C     P       (input)  INTEGER
C             The number of columns of the positive generator.  P >= K.
C
C     Q       (input)  INTEGER
C             The number of columns in B.
C             If TYPEG = 'D',        Q >= K;
C             If TYPEG = 'C' or 'R', Q >= 0.
C
C     NB      (input)  INTEGER
C             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies
C             the block size to be used in the blocked parts of the
C             algorithm. NB must be equivalent with the used block size
C             in the routine MB02CU.
C
C     RNK     (input)  INTEGER
C             If TYPEG = 'D', the number of linearly independent columns
C             in the generator as returned by MB02CU.  0 <= RNK <= K.
C             If TYPEG = 'C' or 'R', the value of this parameter is
C             irrelevant.
C
C     A1      (input)  DOUBLE PRECISION array, dimension
C             (LDA1, K)
C             On entry, if TYPEG = 'D', the leading K-by-K part of this
C             array must contain the matrix A1 as returned by MB02CU.
C             If TYPEG = 'C' or 'R', this array is not referenced.
C
C     LDA1    INTEGER
C             The leading dimension of the array A1.
C             If TYPEG = 'D',                   LDA1 >= MAX(1,K);
C             if TYPEG = 'C'  or  TYPEG = 'R',  LDA1 >= 1.
C
C     A2      (input)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K);
C             if TYPEG = 'R',                   dimension (LDA2, K).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-(P-K) part of this array must contain the matrix
C             A2 as returned by MB02CU.
C             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of
C             this array must contain the matrix A2 as returned by
C             MB02CU.
C
C     LDA2    INTEGER
C             The leading dimension of the array A2.
C             If P = K,                  LDA2 >= 1;
C             If P > K and (TYPEG = 'D' or TYPEG = 'C'),
C                                        LDA2 >= MAX(1,K);
C             if P > K and TYPEG = 'R',  LDA2 >= P-K.
C
C     B       (input)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q);
C             if TYPEG = 'R',                   dimension (LDB, K).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-Q part of this array must contain the matrix B as
C             returned by MB02CU.
C             On entry, if TYPEG = 'R', the leading Q-by-K part of this
C             array must contain the matrix B as returned by MB02CU.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             If Q = 0,                  LDB >= 1;
C             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'),
C                                        LDB >= MAX(1,K);
C             if Q > 0 and TYPEG = 'R',  LDB >= Q.
C
C     F1      (input/output)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF1, K);
C             if TYPEG = 'R',                   dimension (LDF1, N).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-K part of this array must contain the first part
C             of the positive generator to be processed.
C             On entry, if TYPEG = 'R', the leading K-by-N part of this
C             array must contain the first part of the positive
C             generator to be processed.
C             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-K part of this array contains the first part of the
C             transformed positive generator.
C             On exit, if TYPEG = 'R', the leading K-by-N part of this
C             array contains the first part of the transformed positive
C             generator.
C
C     LDF1    INTEGER
C             The leading dimension of the array F1.
C             If TYPEG = 'D'  or  TYPEG = 'C',   LDF1 >= MAX(1,N);
C             if TYPEG = 'R',                    LDF1 >= MAX(1,K).
C
C     F2      (input/output)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF2, P-K);
C             if TYPEG = 'R',                   dimension (LDF2, N).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-(P-K) part of this array must contain the second part
C             of the positive generator to be processed.
C             On entry, if TYPEG = 'R', the leading (P-K)-by-N part of
C             this array must contain the second part of the positive
C             generator to be processed.
C             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-(P-K) part of this array contains the second part of
C             the transformed positive generator.
C             On exit, if TYPEG = 'R', the leading (P-K)-by-N part of
C             this array contains the second part of the transformed
C             positive generator.
C
C     LDF2    INTEGER
C             The leading dimension of the array F2.
C             If P = K,                  LDF2 >= 1;
C             If P > K and (TYPEG = 'D' or TYPEG = 'C'),
C                                        LDF2 >= MAX(1,N);
C             if P > K and TYPEG = 'R',  LDF2 >= P-K.
C
C     G       (input/output)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDG, Q);
C             if TYPEG = 'R',                   dimension (LDG, N).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-Q part of this array must contain the negative part
C             of the generator to be processed.
C             On entry, if TYPEG = 'R', the leading Q-by-N part of this
C             array must contain the negative part of the generator to
C             be processed.
C             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             N-by-Q part of this array contains the transformed
C             negative generator.
C             On exit, if TYPEG = 'R', the leading Q-by-N part of this
C             array contains the transformed negative generator.
C
C     LDG     INTEGER
C             The leading dimension of the array G.
C             If Q = 0,                  LDG >= 1;
C             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'),
C                                        LDG >= MAX(1,N);
C             if Q > 0 and TYPEG = 'R',  LDG >= Q.
C
C     CS      (input)  DOUBLE PRECISION array, dimension (x)
C             If TYPEG = 'D' and P = K,                   x = 3*K;
C             If TYPEG = 'D' and P > K,                   x = 5*K;
C             If (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K;
C             If (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K.
C             On entry, the first x elements of this array must contain
C             Givens and modified hyperbolic rotation parameters, and
C             scalar factors of the Householder transformations as
C             returned by MB02CU.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = -23,  DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             TYPEG = 'D':               LDWORK >= MAX(1,N);
C             (TYPEG = 'C' or TYPEG = 'R')  and  NB <= 0:
C                                        LDWORK >= MAX(1,N);
C             (TYPEG = 'C' or TYPEG = 'R')  and  NB >= 1:
C                                        LDWORK >= MAX(1,( N + K )*NB).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0(N*K*( P + Q )) floating point operations.
C
C     METHOD
C
C     The Householder transformations and modified hyperbolic rotations
C     computed by SLICOT Library routine MB02CU are applied to the
C     corresponding parts of the matrices F1, F2 and G.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C     March 2004, March 2007.
C
C     KEYWORDS
C
C     Elementary matrix operations, Householder transformation, matrix
C     operations, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER          STRUCG, TYPEG
      INTEGER            INFO, K, LDA1, LDA2, LDB, LDF1, LDF2, LDG,
     $                   LDWORK, N, NB, P, Q, RNK
C     .. Array Arguments ..
      DOUBLE PRECISION   A1(LDA1,*), A2(LDA2,*), B(LDB,*), CS(*),
     $                   DWORK(*), F1(LDF1,*), F2(LDF2,*), G(LDG,*)
C     .. Local Scalars ..
      INTEGER            COL2, I, IB, J, JJ, LEN, NBL, POS, PST2,
     $                   WRKMIN
      DOUBLE PRECISION   ALPHA, BETA, C, S, TAU, TEMP
      LOGICAL            LRDEF, LTRI, LCOL
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARF, DLARFB, DLARFT, DROT, DSCAL,
     $                   XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      COL2  = MAX( 0, P - K )
      LRDEF = LSAME( TYPEG,  'D' )
      LCOL  = LSAME( TYPEG,  'C' )
      LTRI  = LSAME( STRUCG, 'T' )
      IF ( LRDEF ) THEN
         WRKMIN = MAX( 1, N )
      ELSE
         IF ( NB.GE.1 ) THEN
            WRKMIN = MAX( 1, ( N + K )*NB )
         ELSE
            WRKMIN = MAX( 1, N )
         END IF
      END IF
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( LCOL .OR. LRDEF .OR. LSAME( TYPEG, 'R' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( LTRI .OR. LSAME( STRUCG, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( K.LT.0 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( P.LT.K ) THEN
         INFO = -5
      ELSE IF ( Q.LT.0 .OR. ( LRDEF .AND. Q.LT.K ) ) THEN
         INFO = -6
      ELSE IF ( LRDEF .AND. ( RNK.LT.0 .OR. RNK.GT.K ) ) THEN
         INFO = -8
      ELSE IF ( ( LDA1.LT.1 ) .OR. ( LRDEF .AND. LDA1.LT.K ) ) THEN
         INFO = -10
      ELSE IF ( ( ( P.EQ.K ) .AND. LDA2.LT.1 ) .OR.
     $          ( ( P.GT.K ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDA2.LT.MAX( 1, K ) ) ) .OR.
     $          ( ( P.GT.K ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDA2.LT.( P-K ) ) ) ) THEN
         INFO = -12
      ELSE IF ( ( ( Q.EQ.0 ) .AND. LDB.LT.1 ) .OR.
     $          ( ( Q.GT.0 ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDB.LT.MAX( 1, K ) ) ) .OR.
     $          ( ( Q.GT.0 ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDB.LT.Q ) ) ) THEN
         INFO = -14
      ELSE IF ( ( LRDEF .OR. LCOL ) .AND. LDF1.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF ( (.NOT.( LRDEF .OR. LCOL ) ) .AND. LDF1.LT.MAX( 1, K ) )
     $      THEN
         INFO = -16
      ELSE IF ( ( ( P.EQ.K ) .AND. LDF2.LT.1 ) .OR.
     $          ( ( P.GT.K ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDF2.LT.MAX( 1, N ) ) ) .OR.
     $          ( ( P.GT.K ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDF2.LT.( P-K ) ) ) ) THEN
         INFO = -18
      ELSE IF ( ( ( Q.EQ.0 ) .AND. LDG.LT.1 ) .OR.
     $          ( ( Q.GT.0 ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDG.LT.MAX( 1, N ) ) ) .OR.
     $          ( ( Q.GT.0 ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDG.LT.Q ) ) ) THEN
         INFO = -20
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -23
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02CV', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( K, N ).EQ.0 .OR.
     $     ( ( .NOT.LRDEF ) .AND. Q.EQ.0 .AND. P.EQ.K ) ) THEN
         RETURN
      END IF
C
      IF ( LRDEF ) THEN
C
C        Deficient generator.
C
         IF ( COL2.EQ.0 ) THEN
            PST2 = 2*K
         ELSE
            PST2 = 4*K
         END IF
C
         DO 10  I = 1, RNK
C
C           Apply elementary reflectors.
C
            IF ( COL2.GT.1 ) THEN
               TAU     = A2(I,1)
               A2(I,1) = ONE
               CALL DLARF( 'Right', N, COL2, A2(I,1), LDA2, TAU, F2,
     $                     LDF2, DWORK )
               A2(I,1) = TAU
            END IF
C
            IF ( K.GT.I ) THEN
               ALPHA   = A1(I,I)
               A1(I,I) = ONE
               CALL DLARF( 'Right', N, K-I+1, A1(I,I), LDA1, CS(PST2+I),
     $                     F1(1,I), LDF1, DWORK )
               A1(I,I) = ALPHA
            END IF
C
            IF ( COL2.GT.0 ) THEN
               C = CS(2*K+I*2-1)
               S = CS(2*K+I*2)
               CALL DROT( N, F1(1,I), 1, F2, 1, C, S )
            END IF
C
            IF ( Q.GT.1 ) THEN
               TAU    = B(I,1)
               B(I,1) = ONE
               CALL DLARF( 'Right', N, Q, B(I,1), LDB, TAU,
     $                     G, LDG, DWORK )
               B(I,1) = TAU
            END IF
C
C           Apply hyperbolic rotation.
C
            C = CS(I*2-1)
            S = CS(I*2)
            CALL DSCAL( N, ONE/C, F1(1,I), 1 )
            CALL DAXPY( N,  -S/C, G(1,1), 1, F1(1,I), 1 )
            CALL DSCAL( N,     C, G(1,1), 1 )
            CALL DAXPY( N,    -S, F1(1,I), 1, G(1,1), 1 )
   10    CONTINUE
C
         LEN = Q
         POS = 1
C
         DO 20 J = RNK + 1, K
C
C           Apply the reductions working on singular rows.
C
            IF ( COL2.GT.1 ) THEN
               TAU     = A2(J,1)
               A2(J,1) = ONE
               CALL DLARF( 'Right', N, COL2, A2(J,1), LDA2, TAU, F2,
     $                     LDF2, DWORK )
               A2(J,1) = TAU
            END IF
            IF ( K.GT.J ) THEN
               ALPHA   = A1(J,J)
               A1(J,J) = ONE
               CALL DLARF( 'Right', N, K-J+1, A1(J,J), LDA1, CS(PST2+J),
     $                     F1(1,J), LDF1, DWORK )
               A1(J,J) = ALPHA
            END IF
            IF ( COL2.GT.0 ) THEN
               C = CS(2*K+J*2-1)
               S = CS(2*K+J*2)
               CALL DROT( N, F1(1,J), 1, F2, 1, C, S )
            END IF
            IF ( LEN.GT.1 ) THEN
               BETA     = B(J,POS)
               B(J,POS) = ONE
               CALL DLARF( 'Right', N, LEN, B(J,POS), LDB, CS(J*2-1),
     $                     G(1,POS), LDG, DWORK )
               B(J,POS) = BETA
            END IF
            LEN = LEN - 1
            POS = POS + 1
   20    CONTINUE
C
      ELSE IF ( LCOL ) THEN
C
C        Column oriented and not deficient generator.
C
C        Apply an LQ like hyperbolic/orthogonal blocked decomposition.
C
         IF ( LTRI ) THEN
            LEN = MAX( N - K, 0 )
         ELSE
            LEN = N
         END IF
         IF ( COL2.GT.0 ) THEN
C
            NBL = MIN( COL2, NB )
            IF ( NBL.GT.0 ) THEN
C
C              Blocked version.
C
               DO 50  I = 1, K - NBL + 1, NBL
                  IB = MIN( K-I+1, NBL )
                  CALL DLARFT( 'Forward', 'Rowwise', COL2, IB, A2(I,1),
     $                         LDA2, CS(4*K+I), DWORK, N+K )
                  CALL DLARFB( 'Right', 'No Transpose', 'Forward',
     $                         'Rowwise', LEN, COL2, IB, A2(I,1),
     $                         LDA2, DWORK, N+K, F2, LDF2,
     $                         DWORK(IB+1), N+K )
C
                  DO 40  J = I, I + IB - 1
                     TAU     = A2(J,1)
                     A2(J,1) = ONE
                     CALL DLARF( 'Right', LEN, MIN( COL2, J-I+1 ),
     $                           A2(J,1), LDA2, TAU, F2, LDF2, DWORK )
                     A2(J,1) = TAU
                     C = CS(2*K+J*2-1)
                     S = CS(2*K+J*2)
                     CALL DROT( LEN, F1(1,J), 1, F2, 1, C, S )
                     IF ( LTRI ) THEN
                        LEN  = LEN + 1
                        TEMP = F1(LEN,J)
                        F1(LEN,J) =  C*TEMP
                        F2(LEN,1) = -S*TEMP
C
                        DO 30  JJ = 2, COL2
                           F2(LEN,JJ) = ZERO
   30                   CONTINUE
C
                     END IF
   40             CONTINUE
C
   50          CONTINUE
C
            ELSE
               I = 1
            END IF
C
C           Unblocked version for the last or only block.
C
            DO 70  J = I, K
               IF ( COL2.GT.1 ) THEN
                  TAU     = A2(J,1)
                  A2(J,1) = ONE
                  CALL DLARF( 'Right', LEN, COL2, A2(J,1), LDA2, TAU,
     $                        F2, LDF2, DWORK )
                  A2(J,1) = TAU
               END IF
C
               C = CS(2*K+J*2-1)
               S = CS(2*K+J*2)
               CALL DROT( LEN, F1(1,J), 1, F2, 1, C, S )
               IF ( LTRI ) THEN
                  LEN = LEN + 1
                  TEMP = F1(LEN,J)
                  F1(LEN,J) =  C*TEMP
                  F2(LEN,1) = -S*TEMP
C
                  DO 60  JJ = 2, COL2
                     F2(LEN,JJ) = ZERO
   60             CONTINUE
C
               END IF
   70       CONTINUE
C
            PST2 = 5*K
         ELSE
            PST2 = 2*K
         END IF
C
         IF ( LTRI ) THEN
            LEN = N - K
         ELSE
            LEN = N
         END IF
C
         NBL = MIN( Q, NB )
         IF ( NBL.GT.0 ) THEN
C
C           Blocked version.
C
            DO 100  I = 1, K - NBL + 1, NBL
               IB = MIN( K-I+1, NBL )
               CALL DLARFT( 'Forward', 'Rowwise', Q, IB, B(I,1),
     $                      LDB, CS(PST2+I), DWORK, N+K )
               CALL DLARFB( 'Right', 'NonTranspose', 'Forward',
     $                      'Rowwise', LEN, Q, IB, B(I,1),
     $                      LDB, DWORK, N+K, G, LDG,
     $                      DWORK(IB+1), N+K )
C
               DO 90  J = I, I + IB - 1
                  TAU    = B(J,1)
                  B(J,1) = ONE
                  CALL DLARF( 'Right', LEN, J-I+1, B(J,1), LDB,
     $                         TAU, G, LDG, DWORK )
                  B(J,1) = TAU
C
C                 Apply hyperbolic rotation.
C
                  C = CS(J*2-1)
                  S = CS(J*2)
                  CALL DSCAL( LEN, ONE/C, F1(1,J), 1 )
                  CALL DAXPY( LEN,  -S/C, G, 1, F1(1,J), 1 )
                  CALL DSCAL( LEN,     C, G, 1 )
                  CALL DAXPY( LEN,    -S, F1(1,J), 1, G, 1 )
                  IF ( LTRI ) THEN
                     LEN = LEN + 1
                     G(LEN,1)  = -S/C*F1(LEN,J)
                     F1(LEN,J) = F1(LEN,J) / C
C
                     DO 80  JJ = 2, Q
                        G(LEN,JJ) = ZERO
   80                CONTINUE
C
                  END IF
   90          CONTINUE
C
  100       CONTINUE
C
         ELSE
            I = 1
         END IF
C
C        Unblocked version for the last or only block.
C
         DO 120  J = I, K
            IF ( Q.GT.1 ) THEN
               TAU    = B(J,1)
               B(J,1) = ONE
               CALL DLARF( 'Right', LEN, Q, B(J,1), LDB, TAU,
     $                     G, LDG, DWORK )
               B(J,1) = TAU
            END IF
            IF ( Q.GT.0 ) THEN
C
C              Apply hyperbolic rotation.
C
               C = CS(J*2-1)
               S = CS(J*2)
               CALL DSCAL( LEN, ONE/C, F1(1,J), 1 )
               CALL DAXPY( LEN,  -S/C, G, 1, F1(1,J), 1 )
               CALL DSCAL( LEN,     C, G, 1 )
               CALL DAXPY( LEN,    -S, F1(1,J), 1, G, 1 )
               IF ( LTRI ) THEN
                  LEN = LEN + 1
                  G(LEN,1)  = -S/C*F1(LEN,J)
                  F1(LEN,J) = F1(LEN,J) / C
C
                  DO 110  JJ = 2, Q
                     G(LEN,JJ) = ZERO
  110             CONTINUE
C
               END IF
            END IF
  120    CONTINUE
C
      ELSE
C
C        Row oriented and not deficient generator.
C
         IF ( LTRI ) THEN
            LEN = MAX( N - K, 0 )
         ELSE
            LEN = N
         END IF
C
         IF ( COL2.GT.0 ) THEN
            NBL = MIN( NB, COL2 )
            IF ( NBL.GT.0 ) THEN
C
C              Blocked version.
C
               DO 150  I = 1, K - NBL + 1, NBL
                  IB = MIN( K-I+1, NBL )
                  CALL DLARFT( 'Forward', 'Columnwise', COL2, IB,
     $                         A2(1,I), LDA2, CS(4*K+I), DWORK, N+K )
                  CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                         'Columnwise', COL2, LEN, IB, A2(1,I),
     $                         LDA2, DWORK, N+K, F2, LDF2,
     $                         DWORK(IB+1), N+K )
C
                  DO 140  J = I, I + IB - 1
                     TAU     = A2(1,J)
                     A2(1,J) = ONE
                     CALL DLARF( 'Left', MIN( COL2, J-I+1 ), LEN,
     $                           A2(1,J), 1, TAU, F2, LDF2, DWORK )
                     A2(1,J) = TAU
                     C = CS(2*K+J*2-1)
                     S = CS(2*K+J*2)
                     CALL DROT( LEN, F1(J,1), LDF1, F2, LDF2, C, S )
                     IF ( LTRI ) THEN
                        LEN  = LEN + 1
                        TEMP = F1(J,LEN)
                        F1(J,LEN) =  C*TEMP
                        F2(1,LEN) = -S*TEMP
C
                        DO 130  JJ = 2, COL2
                           F2(JJ,LEN) = ZERO
  130                   CONTINUE
C
                     END IF
  140             CONTINUE
C
  150          CONTINUE
C
            ELSE
               I = 1
            END IF
C
C           Unblocked version for the last or only block.
C
            DO 170  J = I, K
               IF ( COL2.GT.1 ) THEN
                  TAU     = A2(1,J)
                  A2(1,J) = ONE
                  CALL DLARF( 'Left', COL2, LEN, A2(1,J), 1, TAU,
     $                        F2, LDF2, DWORK )
                  A2(1,J) = TAU
               END IF
C
               C = CS(2*K+J*2-1)
               S = CS(2*K+J*2)
               CALL DROT( LEN, F1(J,1), LDF1, F2, LDF2, C, S )
               IF ( LTRI ) THEN
                  LEN  = LEN + 1
                  TEMP = F1(J,LEN)
                  F1(J,LEN) =  C*TEMP
                  F2(1,LEN) = -S*TEMP
C
                  DO 160  JJ = 2, COL2
                     F2(JJ,LEN) = ZERO
  160             CONTINUE
C
               END IF
  170       CONTINUE
C
            PST2 = 5*K
         ELSE
            PST2 = 2*K
         END IF
C
         IF ( LTRI ) THEN
            LEN = N - K
         ELSE
            LEN = N
         END IF
C
         NBL = MIN( Q, NB )
         IF ( NBL.GT.0 ) THEN
C
C           Blocked version.
C
            DO 200  I = 1, K - NBL + 1, NBL
               IB = MIN( K-I+1, NBL )
               CALL DLARFT( 'Forward', 'Columnwise', Q, IB, B(1,I),
     $                      LDB, CS(PST2+I), DWORK, N+K )
               CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', Q, LEN, IB, B(1,I),
     $                      LDB, DWORK, N+K, G, LDG,
     $                      DWORK(IB+1), N+K )
C
               DO 190  J = I, I + IB - 1
                  TAU    = B(1,J)
                  B(1,J) = ONE
                  CALL DLARF( 'Left', J-I+1, LEN, B(1,J), 1,
     $                         TAU, G, LDG, DWORK )
                  B(1,J) = TAU
C
C                 Apply hyperbolic rotation.
C
                  C = CS(J*2-1)
                  S = CS(J*2)
                  CALL DSCAL( LEN, ONE/C, F1(J,1), LDF1 )
                  CALL DAXPY( LEN,  -S/C, G, LDG, F1(J,1), LDF1 )
                  CALL DSCAL( LEN,     C, G, LDG )
                  CALL DAXPY( LEN,    -S, F1(J,1), LDF1, G, LDG )
                  IF ( LTRI ) THEN
                     LEN = LEN + 1
                     G(1,LEN)  = -S/C*F1(J,LEN)
                     F1(J,LEN) = F1(J,LEN) / C
C
                     DO 180  JJ = 2, Q
                        G(JJ,LEN) = ZERO
  180                CONTINUE
C
                  END IF
  190          CONTINUE
C
  200       CONTINUE
C
         ELSE
            I = 1
         END IF
C
C        Unblocked version for the last or only block.
C
         DO 220  J = I, K
            IF ( Q.GT.1 ) THEN
               TAU    = B(1,J)
               B(1,J) = ONE
               CALL DLARF( 'Left', Q, LEN, B(1,J), 1, TAU,
     $                     G, LDG, DWORK )
               B(1,J) = TAU
            END IF
            IF ( Q.GT.0 ) THEN
C
C              Apply hyperbolic rotation.
C
               C = CS(J*2-1)
               S = CS(J*2)
               CALL DSCAL( LEN, ONE/C, F1(J,1), LDF1 )
               CALL DAXPY( LEN,  -S/C, G, LDG, F1(J,1), LDF1 )
               CALL DSCAL( LEN,     C, G, LDG )
               CALL DAXPY( LEN,    -S, F1(J,1), LDF1, G, LDG )
               IF ( LTRI ) THEN
                  LEN = LEN + 1
                  G(1,LEN)  = -S/C*F1(J,LEN)
                  F1(J,LEN) = F1(J,LEN) / C
C
                  DO 210  JJ = 2, Q
                     G(JJ,LEN) = ZERO
  210             CONTINUE
C
               END IF
            END IF
  220    CONTINUE
C
      END IF
C
C *** Last line of MB02CV ***
      END
