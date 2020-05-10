      SUBROUTINE MB04WD( TRANQ1, TRANQ2, M, N, K, Q1, LDQ1, Q2, LDQ2,
     $                   CS, TAU, DWORK, LDWORK, INFO )
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
C     To generate a matrix Q with orthogonal columns (spanning an
C     isotropic subspace), which is defined as the first n columns
C     of a product of symplectic reflectors and Givens rotators,
C
C         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )
C             diag( H(2),H(2) ) G(2) diag( F(2),F(2) )
C                               ....
C             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).
C
C     The matrix Q is returned in terms of its first 2*M rows
C
C                      [  op( Q1 )   op( Q2 ) ]
C                  Q = [                      ].
C                      [ -op( Q2 )   op( Q1 ) ]
C
C     Blocked version of the SLICOT Library routine MB04WU.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANQ1  CHARACTER*1
C             Specifies the form of op( Q1 ) as follows:
C             = 'N':  op( Q1 ) = Q1;
C             = 'T':  op( Q1 ) = Q1';
C             = 'C':  op( Q1 ) = Q1'.
C
C     TRANQ2  CHARACTER*1
C             Specifies the form of op( Q2 ) as follows:
C             = 'N':  op( Q2 ) = Q2;
C             = 'T':  op( Q2 ) = Q2';
C             = 'C':  op( Q2 ) = Q2'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrices Q1 and Q2. M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrices Q1 and Q2.
C             M >= N >= 0.
C
C     K       (input) INTEGER
C             The number of symplectic Givens rotators whose product
C             partly defines the matrix Q. N >= K >= 0.
C
C     Q1      (input/output) DOUBLE PRECISION array, dimension
C                     (LDQ1,N) if TRANQ1 = 'N',
C                     (LDQ1,M) if TRANQ1 = 'T' or TRANQ1 = 'C'
C             On entry with TRANQ1 = 'N', the leading M-by-K part of
C             this array must contain in its i-th column the vector
C             which defines the elementary reflector F(i).
C             On entry with TRANQ1 = 'T' or TRANQ1 = 'C', the leading
C             K-by-M part of this array must contain in its i-th row
C             the vector which defines the elementary reflector F(i).
C             On exit with TRANQ1 = 'N', the leading M-by-N part of this
C             array contains the matrix Q1.
C             On exit with TRANQ1 = 'T' or TRANQ1 = 'C', the leading
C             N-by-M part of this array contains the matrix Q1'.
C
C     LDQ1    INTEGER
C             The leading dimension of the array Q1.
C             LDQ1 >= MAX(1,M),  if TRANQ1 = 'N';
C             LDQ1 >= MAX(1,N),  if TRANQ1 = 'T' or TRANQ1 = 'C'.
C
C     Q2      (input/output) DOUBLE PRECISION array, dimension
C                     (LDQ2,N) if TRANQ2 = 'N',
C                     (LDQ2,M) if TRANQ2 = 'T' or TRANQ2 = 'C'
C             On entry with TRANQ2 = 'N', the leading M-by-K part of
C             this array must contain in its i-th column the vector
C             which defines the elementary reflector H(i) and, on the
C             diagonal, the scalar factor of H(i).
C             On entry with TRANQ2 = 'T' or TRANQ2 = 'C', the leading
C             K-by-M part of this array must contain in its i-th row the
C             vector which defines the elementary reflector H(i) and, on
C             the diagonal, the scalar factor of H(i).
C             On exit with TRANQ2 = 'N', the leading M-by-N part of this
C             array contains the matrix Q2.
C             On exit with TRANQ2 = 'T' or TRANQ2 = 'C', the leading
C             N-by-M part of this array contains the matrix Q2'.
C
C     LDQ2    INTEGER
C             The leading dimension of the array Q2.
C             LDQ2 >= MAX(1,M),  if TRANQ2 = 'N';
C             LDQ2 >= MAX(1,N),  if TRANQ2 = 'T' or TRANQ2 = 'C'.
C
C     CS      (input) DOUBLE PRECISION array, dimension (2*K)
C             On entry, the first 2*K elements of this array must
C             contain the cosines and sines of the symplectic Givens
C             rotators G(i).
C
C     TAU     (input) DOUBLE PRECISION array, dimension (K)
C             On entry, the first K elements of this array must
C             contain the scalar factors of the elementary reflectors
C             F(i).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK, MAX(M+N,8*N*NB + 15*NB*NB), where NB is
C             the optimal block size determined by the function UE01MD.
C             On exit, if  INFO = -13,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,M+N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSB).
C
C     KEYWORDS
C
C     Elementary matrix operations, orthogonal symplectic matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANQ1, TRANQ2
      INTEGER           INFO, K, LDQ1, LDQ2, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), Q1(LDQ1,*), Q2(LDQ2,*), TAU(*)
C     .. Local Scalars ..
      LOGICAL           LTRQ1, LTRQ2
      INTEGER           I, IB, IERR, KI, KK, NB, NBMIN, NX, PDRS, PDT,
     $                  PDW, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           UE01MD
      EXTERNAL          LSAME, UE01MD
C     .. External Subroutines ..
      EXTERNAL          MB04QC, MB04QF, MB04WU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      LTRQ1 = LSAME( TRANQ1, 'T' ) .OR. LSAME( TRANQ1,'C' )
      LTRQ2 = LSAME( TRANQ2, 'T' ) .OR. LSAME( TRANQ2,'C' )
      NB = UE01MD( 1, 'MB04WD', TRANQ1 // TRANQ2, M, N, K )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( LTRQ1 .OR. LSAME( TRANQ1, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( LTRQ2 .OR. LSAME( TRANQ2, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( M.LT.0 ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -4
      ELSE IF ( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -5
      ELSE IF ( ( LTRQ1 .AND. LDQ1.LT.MAX( 1, N ) ) .OR.
     $     ( .NOT.LTRQ1 .AND. LDQ1.LT.MAX( 1, M ) ) ) THEN
         INFO = -7
      ELSE IF ( ( LTRQ2 .AND. LDQ2.LT.MAX( 1, N ) ) .OR.
     $     ( .NOT.LTRQ2 .AND. LDQ2.LT.MAX( 1, M ) ) ) THEN
         INFO = -9
      ELSE IF ( LDWORK.LT.MAX( 1, M + N ) ) THEN
         DWORK(1) = DBLE( MAX( 1, M + N ) )
         INFO = -13
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04WD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      WRKOPT = M + N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, UE01MD( 3, 'MB04WD', TRANQ1 // TRANQ2, M, N, K ) )
         IF ( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            WRKOPT = MAX( WRKOPT, 8*N*NB + 15*NB*NB )
            IF( LDWORK.LT.WRKOPT ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = INT( ( SQRT( DBLE( 16*N*N + 15*LDWORK ) )
     $                     - DBLE( 4*N ) ) / 15.0D0 )
               NBMIN = MAX( 2, UE01MD( 2, 'MB04WD', TRANQ1 // TRANQ2, M,
     $                                 N, K ) )
            END IF
         END IF
      END IF
C
      IF ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code after the last block.
C        The first kk columns are handled by the block method.
C
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the last or only block.
C
      IF ( KK.LT.N )
     $   CALL MB04WU( TRANQ1, TRANQ2, M-KK, N-KK, K-KK, Q1(KK+1,KK+1),
     $                LDQ1, Q2(KK+1,KK+1), LDQ2, CS(2*KK+1), TAU(KK+1),
     $                DWORK, LDWORK, IERR )
C
C     Blocked code.
C
      IF ( KK.GT.0 ) THEN
         PDRS = 1
         PDT  = PDRS + 6*NB*NB
         PDW  = PDT  + 9*NB*NB
         IF ( LTRQ1.AND.LTRQ2 ) THEN
            DO 10 I = KI + 1, 1, -NB
               IB = MIN( NB, K-I+1 )
               IF ( I+IB.LE.N ) THEN
C
C                 Form the triangular factors of the symplectic block
C                 reflector SH.
C
                  CALL MB04QF( 'Forward', 'Rowwise', 'Rowwise', M-I+1,
     $                         IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         CS(2*I-1), TAU(I), DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, DWORK(PDW) )
C
C                 Apply SH to Q1(i+ib:n,i:m) and Q2(i+ib:n,i:m) from
C                 the right.
C
                  CALL MB04QC( 'Zero Structure', 'Transpose',
     $                         'Transpose', 'No Transpose', 'Forward',
     $                         'Rowwise', 'Rowwise', M-I+1, N-I-IB+1,
     $                         IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         DWORK(PDRS), NB, DWORK(PDT), NB,
     $                         Q2(I+IB,I), LDQ2, Q1(I+IB,I), LDQ1,
     $                         DWORK(PDW) )
               END IF
C
C              Apply SH to columns i:m of the current block.
C
               CALL MB04WU( 'Transpose', 'Transpose', M-I+1, IB, IB,
     $                      Q1(I,I), LDQ1, Q2(I,I), LDQ2, CS(2*I-1),
     $                      TAU(I), DWORK, LDWORK, IERR )
   10       CONTINUE
C
         ELSE IF ( LTRQ1 ) THEN
            DO 20 I = KI + 1, 1, -NB
               IB = MIN( NB, K-I+1 )
               IF ( I+IB.LE.N ) THEN
C
C                 Form the triangular factors of the symplectic block
C                 reflector SH.
C
                  CALL MB04QF( 'Forward', 'Rowwise', 'Columnwise',
     $                         M-I+1, IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         CS(2*I-1), TAU(I), DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, DWORK(PDW) )
C
C                 Apply SH to Q1(i+ib:n,i:m) from the right and to
C                 Q2(i:m,i+ib:n) from the left.
C
                  CALL MB04QC( 'Zero Structure', 'No Transpose',
     $                         'Transpose', 'No Transpose',
     $                         'Forward', 'Rowwise', 'Columnwise',
     $                         M-I+1, N-I-IB+1, IB, Q1(I,I), LDQ1,
     $                         Q2(I,I), LDQ2, DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, Q2(I,I+IB), LDQ2,
     $                         Q1(I+IB,I), LDQ1, DWORK(PDW) )
               END IF
C
C              Apply SH to columns/rows i:m of the current block.
C
               CALL MB04WU( 'Transpose', 'No Transpose', M-I+1, IB, IB,
     $                      Q1(I,I), LDQ1, Q2(I,I), LDQ2, CS(2*I-1),
     $                      TAU(I), DWORK, LDWORK, IERR )
   20       CONTINUE
C
         ELSE IF ( LTRQ2 ) THEN
            DO 30 I = KI + 1, 1, -NB
               IB = MIN( NB, K-I+1 )
               IF ( I+IB.LE.N ) THEN
C
C                 Form the triangular factors of the symplectic block
C                 reflector SH.
C
                  CALL MB04QF( 'Forward', 'Columnwise', 'Rowwise',
     $                         M-I+1, IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         CS(2*I-1), TAU(I), DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, DWORK(PDW) )
C
C                 Apply SH to Q1(i:m,i+ib:n) from the left and to
C                 Q2(i+ib:n,i:m) from the right.
C
                  CALL MB04QC( 'Zero Structure', 'Transpose',
     $                        'No Transpose', 'No Transpose', 'Forward',
     $                         'Columnwise', 'Rowwise', M-I+1, N-I-IB+1,
     $                         IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         DWORK(PDRS), NB, DWORK(PDT), NB,
     $                         Q2(I+IB,I), LDQ2, Q1(I,I+IB), LDQ1,
     $                         DWORK(PDW) )
               END IF
C
C              Apply SH to columns/rows i:m of the current block.
C
               CALL MB04WU( 'No Transpose', 'Transpose', M-I+1, IB, IB,
     $                      Q1(I,I), LDQ1, Q2(I,I), LDQ2, CS(2*I-1),
     $                      TAU(I), DWORK, LDWORK, IERR )
   30       CONTINUE
C
         ELSE
            DO 40 I = KI + 1, 1, -NB
               IB = MIN( NB, K-I+1 )
               IF ( I+IB.LE.N ) THEN
C
C                 Form the triangular factors of the symplectic block
C                 reflector SH.
C
                  CALL MB04QF( 'Forward', 'Columnwise', 'Columnwise',
     $                         M-I+1, IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2,
     $                         CS(2*I-1), TAU(I), DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, DWORK(PDW) )
C
C                 Apply SH to Q1(i:m,i+ib:n) and Q2(i:m,i+ib:n) from
C                 the left.
C
                  CALL MB04QC( 'Zero Structure', 'No Transpose',
     $                         'No Transpose', 'No Transpose',
     $                         'Forward', 'Columnwise', 'Columnwise',
     $                         M-I+1, N-I-IB+1, IB, Q1(I,I), LDQ1,
     $                         Q2(I,I), LDQ2, DWORK(PDRS), NB,
     $                         DWORK(PDT), NB, Q2(I,I+IB), LDQ2,
     $                         Q1(I,I+IB), LDQ1,  DWORK(PDW) )
               END IF
C
C              Apply SH to rows i:m of the current block.
C
               CALL MB04WU( 'No Transpose', 'No Transpose', M-I+1, IB,
     $                      IB, Q1(I,I), LDQ1, Q2(I,I), LDQ2, CS(2*I-1),
     $                      TAU(I), DWORK, LDWORK, IERR )
   40       CONTINUE
         END IF
      END IF
C
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of MB04WD ***
      END
