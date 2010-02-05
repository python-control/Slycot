      SUBROUTINE MB04QB( TRANC, TRAND, TRANQ, STOREV, STOREW, M, N, K,
     $                   V, LDV, W, LDW, C, LDC, D, LDD, CS, TAU, DWORK,
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
C     To overwrite general real m-by-n matrices C and D, or their
C     transposes, with
C
C               [ op(C) ]
C         Q  *  [       ]   if TRANQ = 'N', or
C               [ op(D) ]
C
C          T    [ op(C) ]
C         Q  *  [       ]   if TRANQ = 'T',
C               [ op(D) ]
C
C     where Q is defined as the product of symplectic reflectors and
C     Givens rotators,
C
C         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )
C             diag( H(2),H(2) ) G(2) diag( F(2),F(2) )
C                               ....
C             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).
C
C     Blocked version.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANC   CHARACTER*1
C             Specifies the form of op( C ) as follows:
C             = 'N':  op( C ) = C;
C             = 'T':  op( C ) = C';
C             = 'C':  op( C ) = C'.
C
C     TRAND   CHARACTER*1
C             Specifies the form of op( D ) as follows:
C             = 'N':  op( D ) = D;
C             = 'T':  op( D ) = D';
C             = 'C':  op( D ) = D'.
C
C     TRANQ   CHARACTER*1
C             = 'N':  apply Q;
C             = 'T':  apply Q'.
C
C     STOREV  CHARACTER*1
C             Specifies how the vectors which define the concatenated
C             Householder reflectors contained in V are stored:
C             = 'C':  columnwise;
C             = 'R':  rowwise.
C
C     STOREW  CHARACTER*1
C             Specifies how the vectors which define the concatenated
C             Householder reflectors contained in W are stored:
C             = 'C':  columnwise;
C             = 'R':  rowwise.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrices op(C) and op(D).
C             M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrices op(C) and op(D).
C             N >= 0.
C
C     K       (input) INTEGER
C             The number of elementary reflectors whose product defines
C             the matrix Q.  M >= K >= 0.
C
C     V       (input) DOUBLE PRECISION array, dimension
C                     (LDV,K) if STOREV = 'C',
C                     (LDV,M) if STOREV = 'R'
C             On entry with STOREV = 'C', the leading M-by-K part of
C             this array must contain in its columns the vectors which
C             define the elementary reflectors F(i).
C             On entry with STOREV = 'R', the leading K-by-M part of
C             this array must contain in its rows the vectors which
C             define the elementary reflectors F(i).
C
C     LDV     INTEGER
C             The leading dimension of the array V.
C             LDV >= MAX(1,M),  if STOREV = 'C';
C             LDV >= MAX(1,K),  if STOREV = 'R'.
C
C     W       (input) DOUBLE PRECISION array, dimension
C                     (LDW,K) if STOREW = 'C',
C                     (LDW,M) if STOREW = 'R'
C             On entry with STOREW = 'C', the leading M-by-K part of
C             this array must contain in its columns the vectors which
C             define the elementary reflectors H(i).
C             On entry with STOREW = 'R', the leading K-by-M part of
C             this array must contain in its rows the vectors which
C             define the elementary reflectors H(i).
C
C     LDW     INTEGER
C             The leading dimension of the array W.
C             LDW >= MAX(1,M),  if STOREW = 'C';
C             LDW >= MAX(1,K),  if STOREW = 'R'.
C
C     C       (input/output) DOUBLE PRECISION array, dimension
C                     (LDC,N) if TRANC = 'N',
C                     (LDC,M) if TRANC = 'T' or TRANC = 'C'
C             On entry with TRANC = 'N', the leading M-by-N part of
C             this array must contain the matrix C.
C             On entry with TRANC = 'C' or TRANC = 'T', the leading
C             N-by-M part of this array must contain the transpose of
C             the matrix C.
C             On exit with TRANC = 'N', the leading M-by-N part of
C             this array contains the updated matrix C.
C             On exit with TRANC = 'C' or TRANC = 'T', the leading
C             N-by-M part of this array contains the transpose of the
C             updated matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.
C             LDC >= MAX(1,M),  if TRANC = 'N';
C             LDC >= MAX(1,N),  if TRANC = 'T' or TRANC = 'C'.
C
C     D       (input/output) DOUBLE PRECISION array, dimension
C                     (LDD,N) if TRAND = 'N',
C                     (LDD,M) if TRAND = 'T' or TRAND = 'C'
C             On entry with TRAND = 'N', the leading M-by-N part of
C             this array must contain the matrix D.
C             On entry with TRAND = 'C' or TRAND = 'T', the leading
C             N-by-M part of this array must contain the transpose of
C             the matrix D.
C             On exit with TRAND = 'N', the leading M-by-N part of
C             this array contains the updated matrix D.
C             On exit with TRAND = 'C' or TRAND = 'T', the leading
C             N-by-M part of this array contains the transpose of the
C             updated matrix D.
C
C     LDD     INTEGER
C             The leading dimension of the array D.
C             LDD >= MAX(1,M),  if TRAND = 'N';
C             LDD >= MAX(1,N),  if TRAND = 'T' or TRAND = 'C'.
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
C             value of LDWORK.
C             On exit, if  INFO = -20,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
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
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSMSB).
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
      CHARACTER         STOREV, STOREW, TRANC, TRAND, TRANQ
      INTEGER           INFO, K, LDC, LDD, LDV, LDW, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,*), CS(*), D(LDD,*), DWORK(*), TAU(*),
     $                  V(LDV,*), W(LDW,*)
C     .. Local Scalars ..
      LOGICAL           LCOLV, LCOLW, LTRC, LTRD, LTRQ
      INTEGER           I, IB, IC, ID, IERR, JC, JD, KI, KK, NB, NBMIN,
     $                  NX, PDRS, PDT, PDW, WRKOPT
C     .. External Functions ..
      INTEGER           UE01MD
      LOGICAL           LSAME
      EXTERNAL          LSAME, UE01MD
C     .. External Subroutines ..
      EXTERNAL          MB04QC, MB04QF, MB04QU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      LCOLV = LSAME( STOREV, 'C' )
      LCOLW = LSAME( STOREW, 'C' )
      LTRC  = LSAME( TRANC, 'T' ) .OR. LSAME( TRANC, 'C' )
      LTRD  = LSAME( TRAND, 'T' ) .OR. LSAME( TRAND, 'C' )
      LTRQ  = LSAME( TRANQ, 'T' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( LTRC .OR. LSAME( TRANC, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( LTRD  .OR. LSAME( TRAND,  'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( .NOT.( LTRQ  .OR. LSAME( TRANQ,  'N' ) ) ) THEN
         INFO = -3
      ELSE IF ( .NOT.( LCOLV .OR. LSAME( STOREV, 'R' ) ) ) THEN
         INFO = -4
      ELSE IF ( .NOT.( LCOLW .OR. LSAME( STOREW, 'R' ) ) ) THEN
         INFO = -5
      ELSE IF ( M.LT.0 ) THEN
         INFO = -6
      ELSE IF ( N.LT.0 ) THEN
         INFO = -7
      ELSE IF ( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -8
      ELSE IF ( ( LCOLV .AND. LDV.LT.MAX( 1, M ) ) .OR.
     $     ( .NOT.LCOLV .AND. LDV.LT.MAX( 1, K ) ) ) THEN
         INFO = -10
      ELSE IF ( ( LCOLW .AND. LDW.LT.MAX( 1, M ) ) .OR.
     $     ( .NOT.LCOLW .AND. LDW.LT.MAX( 1, K ) ) ) THEN
         INFO = -12
      ELSE IF ( ( LTRC  .AND. LDC.LT.MAX( 1, N ) ) .OR.
     $     ( .NOT.LTRC  .AND. LDC.LT.MAX( 1, M ) ) ) THEN
         INFO = -14
      ELSE IF ( ( LTRD  .AND. LDD.LT.MAX( 1, N ) ) .OR.
     $     ( .NOT.LTRD  .AND. LDD.LT.MAX( 1, M ) ) ) THEN
         INFO = -16
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = DBLE( MAX( 1, N ) )
         INFO = -20
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04QB', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( K, M, N ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      WRKOPT = N
      NB = UE01MD( 1, 'MB04QB', TRANC // TRAND // TRANQ, M, N, K )
      IF ( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, UE01MD( 3, 'MB04QB', TRANC // TRAND // TRANQ, M,
     $                        N, K ) )
         IF ( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            WRKOPT = MAX( WRKOPT, 9*N*NB + 15*NB*NB )
            IF ( LDWORK.LT.WRKOPT ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = INT( ( SQRT( DBLE( 81*N*N + 60*LDWORK ) )
     $                     - DBLE( 9*N ) ) / 30.0D0 )
               NBMIN = MAX( 2, UE01MD( 2, 'MB04QB', TRANC // TRAND //
     $                                 TRANQ, M, N, K ) )
            END IF
         END IF
      END IF
C
      PDRS = 1
      PDT  = PDRS + 6*NB*NB
      PDW  = PDT  + 9*NB*NB
      IC = 1
      JC = 1
      ID = 1
      JD = 1
C
      IF ( LTRQ ) THEN
C
C        Use blocked code initially.
C
         IF ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
            DO 10 I = 1, K - NX, NB
               IB = MIN( K-I+1, NB )
C
C              Form the triangular factors of the symplectic block
C              reflector SH.
C
               CALL MB04QF( 'Forward', STOREV, STOREW, M-I+1, IB,
     $                      V(I,I), LDV, W(I,I), LDW, CS(2*I-1), TAU(I),
     $                      DWORK(PDRS), NB, DWORK(PDT), NB,
     $                      DWORK(PDW) )
C
C              Apply SH' to [ op(C)(i:m,:); op(D)(i:m,:) ] from the
C              left.
C
               IF ( LTRC ) THEN
                  JC = I
               ELSE
                  IC = I
               END IF
               IF ( LTRD ) THEN
                  JD = I
               ELSE
                  ID = I
               END IF
               CALL MB04QC( 'No Structure', TRANC, TRAND, TRANQ,
     $                      'Forward', STOREV, STOREW, M-I+1, N, IB,
     $                      V(I,I), LDV, W(I,I), LDW, DWORK(PDRS), NB,
     $                      DWORK(PDT), NB, C(IC,JC), LDC, D(ID,JD),
     $                      LDD, DWORK(PDW) )
   10       CONTINUE
         ELSE
            I = 1
         END IF
C
C        Use unblocked code to update last or only block.
C
         IF ( I.LE.K ) THEN
            IF ( LTRC ) THEN
               JC = I
            ELSE
               IC = I
            END IF
            IF ( LTRD ) THEN
               JD = I
            ELSE
               ID = I
            END IF
            CALL MB04QU( TRANC, TRAND, TRANQ, STOREV, STOREW, M-I+1, N,
     $                   K-I+1, V(I,I), LDV, W(I,I), LDW, C(IC,JC), LDC,
     $                   D(ID,JD), LDD, CS(2*I-1), TAU(I), DWORK,
     $                   LDWORK, IERR )
         END IF
      ELSE
         IF ( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C           Use blocked code after the last block.
C           The first kk columns are handled by the block method.
C
            KI = ( ( K-NX-1 ) / NB )*NB
            KK = MIN( K, KI+NB )
         ELSE
            KK = 0
         END IF
C
C        Use unblocked code for the last or only block.
C
         IF ( KK.LT.K ) THEN
            IF ( LTRC ) THEN
               JC = KK + 1
            ELSE
               IC = KK + 1
            END IF
            IF ( LTRD ) THEN
               JD = KK + 1
            ELSE
               ID = KK + 1
            END IF
            CALL MB04QU( TRANC, TRAND, TRANQ, STOREV, STOREW, M-KK, N,
     $                   K-KK, V(KK+1,KK+1), LDV, W(KK+1,KK+1), LDW,
     $                   C(IC,JC), LDC, D(ID,JD), LDD, CS(2*KK+1),
     $                   TAU(KK+1), DWORK, LDWORK, IERR )
         END IF
C
C        Blocked code.
C
         IF ( KK.GT.0 ) THEN
            DO 20 I = KI + 1, 1, -NB
               IB = MIN( NB, K-I+1 )
C
C              Form the triangular factors of the symplectic block
C              reflector SH.
C
               CALL MB04QF( 'Forward', STOREV, STOREW, M-I+1, IB,
     $                      V(I,I), LDV, W(I,I), LDW, CS(2*I-1), TAU(I),
     $                      DWORK(PDRS), NB, DWORK(PDT), NB,
     $                      DWORK(PDW) )
C
C              Apply SH to [ op(C)(i:m,:); op(D)(i:m,:) ] from
C              the left.
C
               IF ( LTRC ) THEN
                  JC = I
               ELSE
                  IC = I
               END IF
               IF ( LTRD ) THEN
                  JD = I
               ELSE
                  ID = I
               END IF
               CALL MB04QC( 'No Structure', TRANC, TRAND, TRANQ,
     $                      'Forward', STOREV, STOREW, M-I+1, N, IB,
     $                      V(I,I), LDV, W(I,I), LDW, DWORK(PDRS), NB,
     $                      DWORK(PDT), NB, C(IC,JC), LDC, D(ID,JD),
     $                      LDD, DWORK(PDW) )
   20       CONTINUE
         END IF
      END IF
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of MB04QB ***
      END
