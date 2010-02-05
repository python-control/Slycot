      SUBROUTINE MB04QC( STRUCT, TRANA, TRANB, TRANQ, DIRECT, STOREV,
     $                   STOREW, M, N, K, V, LDV, W, LDW, RS, LDRS, T,
     $                   LDT, A, LDA, B, LDB, DWORK )
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
C     To apply the orthogonal symplectic block reflector
C
C              [  I+V*T*V'  V*R*S*V'  ]
C         Q =  [                      ]
C              [ -V*R*S*V'  I+V*T*V'  ]
C
C     or its transpose to a real 2m-by-n matrix [ op(A); op(B) ] from
C     the left.
C     The k-by-k upper triangular blocks of the matrices
C
C                                 [ S1 ]       [ T11 T12 T13 ]
C         R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ],
C                                 [ S3 ]       [ T31 T32 T33 ]
C
C     with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular,
C     are stored rowwise in the arrays RS and T, respectively.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     STRUCT  CHARACTER*1
C             Specifies the structure of the first blocks of A and B:
C             = 'Z':  the leading K-by-N submatrices of op(A) and op(B)
C                     are (implicitly) assumed to be zero;
C             = 'N';  no structure to mention.
C
C     TRANA   CHARACTER*1
C             Specifies the form of op( A ) as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     TRANB   CHARACTER*1
C             Specifies the form of op( B ) as follows:
C             = 'N':  op( B ) = B;
C             = 'T':  op( B ) = B';
C             = 'C':  op( B ) = B'.
C
C     DIRECT  CHARACTER*1
C             This is a dummy argument, which is reserved for future
C             extensions of this subroutine. Not referenced.
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
C             The number of rows of the matrices op(A) and op(B).
C             M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrices op(A) and op(B).
C             N >= 0.
C
C     K       (input) INTEGER
C             The order of the triangular matrices defining R, S and T.
C             M >= K >= 0.
C
C     V       (input) DOUBLE PRECISION array, dimension
C                     (LDV,K) if STOREV = 'C',
C                     (LDV,M) if STOREV = 'R'
C             On entry with STOREV = 'C', the leading M-by-K part of
C             this array must contain in its columns the vectors which
C             define the elementary reflector used to form parts of Q.
C             On entry with STOREV = 'R', the leading K-by-M part of
C             this array must contain in its rows the vectors which
C             define the elementary reflector used to form parts of Q.
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
C             define the elementary reflector used to form parts of Q.
C             On entry with STOREW = 'R', the leading K-by-M part of
C             this array must contain in its rows the vectors which
C             define the elementary reflector used to form parts of Q.
C
C     LDW     INTEGER
C             The leading dimension of the array W.
C             LDW >= MAX(1,M),  if STOREW = 'C';
C             LDW >= MAX(1,K),  if STOREW = 'R'.
C
C     RS      (input) DOUBLE PRECISION array, dimension (K,6*K)
C             On entry, the leading K-by-6*K part of this array must
C             contain the upper triangular matrices defining the factors
C             R and S of the symplectic block reflector Q. The
C             (strictly) lower portions of this array are not
C             referenced.
C
C     LDRS    INTEGER
C             The leading dimension of the array RS.  LDRS >= MAX(1,K).
C
C     T       (input) DOUBLE PRECISION array, dimension (K,9*K)
C             On entry, the leading K-by-9*K part of this array must
C             contain the upper triangular matrices defining the factor
C             T of the symplectic block reflector Q. The (strictly)
C             lower portions of this array are not referenced.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= MAX(1,K).
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C                     (LDA,N) if TRANA = 'N',
C                     (LDA,M) if TRANA = 'C' or TRANA = 'T'
C             On entry with TRANA = 'N', the leading M-by-N part of this
C             array must contain the matrix A.
C             On entry with TRANA = 'T' or TRANA = 'C', the leading
C             N-by-M part of this array must contain the matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,M),  if TRANA = 'N';
C             LDA >= MAX(1,N),  if TRANA = 'C' or TRANA = 'T'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C                     (LDB,N) if TRANB = 'N',
C                     (LDB,M) if TRANB = 'C' or TRANB = 'T'
C             On entry with TRANB = 'N', the leading M-by-N part of this
C             array must contain the matrix B.
C             On entry with TRANB = 'T' or TRANB = 'C', the leading
C             N-by-M part of this array must contain the matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,M),  if TRANB = 'N';
C             LDB >= MAX(1,N),  if TRANB = 'C' or TRANB = 'T'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK), where
C             LDWORK >= 8*N*K,   if STRUCT = 'Z',
C             LDWORK >= 9*N*K,   if STRUCT = 'N'.
C
C     REFERENCES
C
C     [1] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 16*( M - K )*N + ( 26*K - 4 )*K*N floating
C     point operations if STRUCT = 'Z' and additional ( 12*K + 2 )*K*N
C     floating point operations if STRUCT = 'N'.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAESB).
C
C     KEYWORDS
C
C     Elementary matrix operations, orthogonal symplectic matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DIRECT, STOREV, STOREW, STRUCT, TRANA, TRANB,
     $                  TRANQ
      INTEGER           K, LDA, LDB, LDRS, LDT, LDV, LDW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), RS(LDRS,*),
     $                  T(LDT,*), V(LDV,*), W(LDW,*)
C     .. Local Scalars ..
      LOGICAL           LA1B1, LCOLV, LCOLW, LTRA, LTRB, LTRQ
      INTEGER           I, ITEMP, PDW1, PDW2, PDW3, PDW4, PDW5, PDW6,
     $                  PDW7, PDW8, PDW9, PR1, PR2, PR3, PS1, PS2, PS3,
     $                  PT11, PT12, PT13, PT21, PT22, PT23, PT31, PT32,
     $                  PT33
      DOUBLE PRECISION  FACT
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DLASET, DTRMM
C
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF ( M.LE.0 .OR. N.LE.0 )
     $   RETURN
      LA1B1 = LSAME( STRUCT, 'N' )
      LCOLV = LSAME( STOREV, 'C' )
      LCOLW = LSAME( STOREW, 'C' )
      LTRA  = LSAME( TRANA,  'T' ) .OR. LSAME( TRANA, 'C' )
      LTRB  = LSAME( TRANB,  'T' ) .OR. LSAME( TRANB, 'C' )
      LTRQ  = LSAME( TRANQ,  'T' ) .OR. LSAME( TRANQ, 'C' )
C
      PR1  = 1
      PR2  = PR1 + K
      PR3  = PR2 + K
      PS1  = PR3 + K
      PS2  = PS1 + K
      PS3  = PS2 + K
      PT11 = 1
      PT12 = PT11 + K
      PT13 = PT12 + K
      PT21 = PT13 + K
      PT22 = PT21 + K
      PT23 = PT22 + K
      PT31 = PT23 + K
      PT32 = PT31 + K
      PT33 = PT32 + K
      PDW1 = 1
      PDW2 = PDW1 + N*K
      PDW3 = PDW2 + N*K
      PDW4 = PDW3 + N*K
      PDW5 = PDW4 + N*K
      PDW6 = PDW5 + N*K
      PDW7 = PDW6 + N*K
      PDW8 = PDW7 + N*K
      PDW9 = PDW8 + N*K
C
C     Update the matrix A.
C
      IF ( LA1B1 ) THEN
C
C        NZ1) DW7 := A1'
C
         IF ( LTRA ) THEN
            DO 10  I = 1, K
               CALL DCOPY( N, A(1,I), 1, DWORK(PDW7+(I-1)*N), 1 )
   10       CONTINUE
         ELSE
            DO 20  I = 1, N
               CALL DCOPY( K, A(1,I), 1, DWORK(PDW7+I-1), N )
   20       CONTINUE
         END IF
C
C        NZ2) DW1 := DW7*W1
C
         CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW1), 1 )
         IF ( LCOLW ) THEN
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                  K, ONE, W, LDW, DWORK(PDW1), N )
         ELSE
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N,
     $                  K, ONE, W, LDW, DWORK(PDW1), N )
         END IF
C
C        NZ3) DW2 := DW7*V1
C
         CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW2), 1 )
         IF ( LCOLV ) THEN
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                  K, ONE, V, LDV, DWORK(PDW2), N )
         ELSE
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N,
     $                  K, ONE, V, LDV, DWORK(PDW2), N )
         END IF
         FACT = ONE
      ELSE
         FACT = ZERO
      END IF
C
C     1) DW1 := A2'*W2
C
      IF ( M.GT.K ) THEN
         IF ( LTRA.AND.LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  A(1,K+1), LDA, W(K+1,1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE IF ( LTRA ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, K, M-K, ONE,
     $                  A(1,K+1), LDA, W(1,K+1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE IF ( LCOLW ) THEN
            CALL DGEMM( 'Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  A(K+1,1), LDA, W(K+1,1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                  A(K+1,1), LDA, W(1,K+1), LDW, FACT, DWORK(PDW1),
     $                  N )
         END IF
      ELSE IF ( .NOT.LA1B1 ) THEN
         CALL DLASET( 'All', N, K, ZERO, ZERO, DWORK(PDW1), N )
      END IF
C
C     2) DW2 := A2'*V2
C
      IF ( M.GT.K ) THEN
         IF ( LTRA.AND.LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  A(1,K+1), LDA, V(K+1,1), LDV, FACT, DWORK(PDW2),
     $                  N )
         ELSE IF ( LTRA ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, K, M-K, ONE,
     $                  A(1,K+1), LDA, V(1,K+1), LDV, FACT, DWORK(PDW2),
     $                  N )
         ELSE IF ( LCOLV ) THEN
            CALL DGEMM( 'Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  A(K+1,1), LDA, V(K+1,1), LDV, FACT, DWORK(PDW2),
     $                  N )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                  A(K+1,1), LDA, V(1,K+1), LDV, FACT, DWORK(PDW2),
     $                  N )
         END IF
      ELSE IF ( .NOT.LA1B1 ) THEN
         CALL DLASET( 'All', N, K, ZERO, ZERO, DWORK(PDW2), N )
      END IF
C
      IF ( LTRQ ) THEN
C
C        3) DW3 := DW1*T11
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW3), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT11), LDT, DWORK(PDW3), N )
C
C        4) DW4 := DW2*T31
C
         CALL DCOPY( N*(K-1), DWORK(PDW2), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, T(1,PT31+1), LDT, DWORK(PDW4), N )
C
C        5) DW3 := DW3 + DW4
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW4), 1, DWORK(PDW3+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ4) DW8 := DW7*T21
C
            CALL DCOPY( N*(K-1), DWORK(PDW7), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                  K-1, ONE, T(1,PT21+1), LDT, DWORK(PDW8), N )
C
C           NZ5) DW3 := DW3 + DW8
C
            CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW3+N), 1 )
         END IF
C
C        6) DW4 := DW1*T12
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K, ONE, T(1,PT12), LDT, DWORK(PDW4), N )
C
C        7) DW5 := DW2*T32
C
         CALL DCOPY( N*(K-1), DWORK(PDW2), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, T(1,PT32+1), LDT, DWORK(PDW5), N )
C
C        8) DW4 := DW4 + DW5
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW5), 1, DWORK(PDW4+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ6) DW8 := DW7*T22
C
            CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                   K, ONE, T(1,PT22), LDT, DWORK(PDW8), N )
C
C           NZ7) DW4 := DW4 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW4), 1 )
         END IF
C
C        9) DW5 := DW2*T33
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT33), LDT, DWORK(PDW5), N )
C
C        10) DW6 := DW1*T13
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW6), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K, ONE, T(1,PT13), LDT, DWORK(PDW6), N )
C
C        11) DW5 := DW5 + DW6
C
         CALL DAXPY( N*K, ONE, DWORK(PDW6), 1, DWORK(PDW5), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ8) DW8 := DW7*T23
C
            CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                  K, ONE, T(1,PT23), LDT, DWORK(PDW8), N )
C
C           NZ9) DW5 := DW5 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW5), 1 )
         END IF
C
C        12) DW1 := DW1*R1
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PR1), LDRS, DWORK(PDW1), N )
C
C        13) DW2 := DW2*R3
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, RS(1,PR3+1), LDRS, DWORK(PDW2), N )
C
C        14) DW1 := DW1 + DW2
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW2), 1, DWORK(PDW1+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ10) DW7 := DW7*R2
C
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $                  ONE, RS(1,PR2), LDRS, DWORK(PDW7), N )
C
C           NZ11) DW1 := DW1 + DW7
C
            CALL DAXPY( N*K, ONE, DWORK(PDW7), 1, DWORK(PDW1), 1 )
         END IF
C
C        Swap Pointers PDW1 <-> PDW2
C
         ITEMP = PDW2
         PDW2 = PDW1
         PDW1 = ITEMP
      ELSE
C
C        3) DW3 := DW1*T11'
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW3), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT11), LDT, DWORK(PDW3), N )
C
C        4) DW4 := DW2*T13'
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT13), LDT, DWORK(PDW4), N )
C
C        5) DW3 := DW3 + DW4
C
         CALL DAXPY( N*K, ONE, DWORK(PDW4), 1, DWORK(PDW3), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ4) DW8 := DW7*T12'
C
            CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, T(1,PT12), LDT, DWORK(PDW8), N )
C
C           NZ5) DW3 := DW3 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW3), 1 )
         END IF
C
C        6) DW4 := DW2*T23'
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT23), LDT, DWORK(PDW4), N )
C
C        7) DW5 := DW1*T21'
C
         CALL DCOPY( N*(K-1), DWORK(PDW1+N), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, T(1,PT21+1), LDT, DWORK(PDW5), N )
C
C        8) DW4 := DW4 + DW5
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW5), 1, DWORK(PDW4), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ6) DW8 := DW7*T22'
C
            CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, T(1,PT22), LDT, DWORK(PDW8), N )
C
C           NZ7) DW4 := DW4 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW4), 1 )
         END IF
C
C        9) DW5 := DW2*T33'
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT33), LDT, DWORK(PDW5), N )
C
C        10) DW6 := DW1*T31'
C
         CALL DCOPY( N*(K-1), DWORK(PDW1+N), 1, DWORK(PDW6), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, T(1,PT31+1), LDT, DWORK(PDW6), N )
C
C        11) DW5 := DW5 + DW6
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW6), 1, DWORK(PDW5), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ8) DW8 := DW7*T32'
C
            CALL DCOPY( N*(K-1), DWORK(PDW7+N), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N,
     $                  K-1, ONE, T(1,PT32+1), LDT, DWORK(PDW8), N )
C
C           NZ9) DW5 := DW5 + DW8
C
            CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW5), 1 )
         END IF
C
C        12) DW1 := DW1*S1'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, RS(1,PS1+1), LDRS, DWORK(PDW1+N), N )
C
C        13) DW2 := DW2*S3'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PS3), LDRS, DWORK(PDW2), N )
C
C        14) DW2 := DW1 + DW2
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW1+N), 1, DWORK(PDW2), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ10) DW7 := DW7*S2'
C
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, RS(1,PS2), LDRS, DWORK(PDW7), N )
C
C           NZ11) DW2 := DW2 + DW7
C
            CALL DAXPY( N*K, ONE, DWORK(PDW7), 1, DWORK(PDW2), 1 )
         END IF
      END IF
C
      IF ( LA1B1 ) THEN
C
C        NZ12) DW9 := B1'
C
         IF ( LTRB ) THEN
            DO 30  I = 1, K
               CALL DCOPY( N, B(1,I), 1, DWORK(PDW9+(I-1)*N), 1 )
   30       CONTINUE
         ELSE
            DO 40  I = 1, N
               CALL DCOPY( K, B(1,I), 1, DWORK(PDW9+I-1), N )
   40       CONTINUE
         END IF
C
C        NZ13) DW1 := DW9*W1
C
         CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW1), 1 )
         IF ( LCOLW ) THEN
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                  K, ONE, W, LDW, DWORK(PDW1), N )
         ELSE
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N,
     $                  K, ONE, W, LDW, DWORK(PDW1), N )
         END IF
C
C        NZ14) DW6 := DW9*V1
C
         CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW6), 1 )
         IF ( LCOLV ) THEN
            CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                  K, ONE, V, LDV, DWORK(PDW6), N )
         ELSE
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N,
     $                  K, ONE, V, LDV, DWORK(PDW6), N )
         END IF
      END IF
C
C     15) DW1 := B2'*W2
C
      IF ( M.GT.K ) THEN
         IF ( LTRB.AND.LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  B(1,K+1), LDB, W(K+1,1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE IF ( LTRB ) THEN
C
C           Critical Position
C
            CALL DGEMM( 'No Transpose', 'Transpose', N, K, M-K, ONE,
     $                  B(1,K+1), LDB, W(1,K+1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE IF ( LCOLW ) THEN
            CALL DGEMM( 'Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  B(K+1,1), LDB, W(K+1,1), LDW, FACT, DWORK(PDW1),
     $                  N )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                  B(K+1,1), LDB, W(1,K+1), LDW, FACT, DWORK(PDW1),
     $                  N )
         END IF
      ELSE IF ( .NOT.LA1B1 ) THEN
         CALL DLASET( 'All', N, K, ZERO, ZERO, DWORK(PDW1), N )
      END IF
C
C     16) DW6 := B2'*V2
C
      IF ( M.GT.K ) THEN
         IF ( LTRB.AND.LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  B(1,K+1), LDB, V(K+1,1), LDV, FACT, DWORK(PDW6),
     $                  N )
         ELSE IF ( LTRB ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, K, M-K, ONE,
     $                  B(1,K+1), LDB, V(1,K+1), LDV, FACT, DWORK(PDW6),
     $                  N )
         ELSE IF ( LCOLV ) THEN
            CALL DGEMM( 'Transpose', 'No Transpose', N, K, M-K, ONE,
     $                  B(K+1,1), LDB, V(K+1,1), LDV, FACT, DWORK(PDW6),
     $                  N )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                  B(K+1,1), LDB, V(1,K+1), LDV, FACT, DWORK(PDW6),
     $                  N )
         END IF
      ELSE IF ( .NOT.LA1B1 ) THEN
         CALL DLASET( 'All', N, K, ZERO, ZERO, DWORK(PDW6), N )
      END IF
C
      IF ( LTRQ ) THEN
C
C        17) DW7 := DW1*R1
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW7), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PR1), LDRS, DWORK(PDW7), N )
C
C        18) DW8 := DW6*R3
C
         CALL DCOPY( N*(K-1), DWORK(PDW6), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, RS(1,PR3+1), LDRS, DWORK(PDW8), N )
C
C        19) DW7 := DW7 + DW8
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW7+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ15) DW8 := DW9*R2
C
            CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $                  ONE, RS(1,PR2), LDRS, DWORK(PDW8), N )
C
C           NZ16) DW7 := DW7 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW7), 1 )
         END IF
C
C        20) DW8 := DW7*S1
C
         CALL DCOPY( N*(K-1), DWORK(PDW7), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, RS(1,PS1+1), LDRS, DWORK(PDW8), N )
C
C        21) DW3 := DW3 - DW8
C
         CALL DAXPY( N*(K-1), -ONE, DWORK(PDW8), 1, DWORK(PDW3+N), 1 )
C
C        22) DW8 := DW7*S3
C
         CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K, ONE, RS(1,PS3), LDRS, DWORK(PDW8), N )
C
C        23) DW5 := DW5 - DW8
C
         CALL DAXPY( N*K, -ONE, DWORK(PDW8), 1, DWORK(PDW5), 1 )
C
C        24) DW7 := DW7*S2
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               -ONE, RS(1,PS2), LDRS, DWORK(PDW7), N )
      ELSE
C
C        17) DW7 := DW6*S3'
C
         CALL DCOPY( N*K, DWORK(PDW6), 1, DWORK(PDW7), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PS3), LDRS, DWORK(PDW7), N )
C
C        18) DW8 := DW1*S1'
C
         CALL DCOPY( N*(K-1), DWORK(PDW1+N), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, RS(1,PS1+1), LDRS, DWORK(PDW8), N )
C
C        19) DW7 := DW7 + DW8
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW7), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ15) DW8 := DW9*S2'
C
            CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, RS(1,PS2), LDRS, DWORK(PDW8), N )
C
C           NZ16) DW7 := DW7 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW7), 1 )
         END IF
C
C        20) DW8 := DW7*R1'
C
         CALL DCOPY( N*K, DWORK(PDW7), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PR1), LDRS, DWORK(PDW8), N )
C
C        21) DW3 := DW3 + DW8
C
         CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW3), 1 )
C
C        22) DW8 := DW7*R3'
C
         CALL DCOPY( N*(K-1), DWORK(PDW7+N), 1, DWORK(PDW8), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, RS(1,PR3+1), LDRS, DWORK(PDW8), N )
C
C        23) DW5 := DW5 + DW8
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW5), 1 )
C
C        24) DW7 := DW7*R2'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $               ONE, RS(1,PR2), LDRS, DWORK(PDW7), N )
      END IF
C
C     25) A2 := A2 + W2*DW3'
C
      IF ( M.GT.K ) THEN
         IF ( LTRA.AND.LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW3), N, W(K+1,1), LDW, ONE, A(1,K+1),
     $                  LDA )
         ELSE IF ( LTRA ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW3), N, W(1,K+1), LDW, ONE, A(1,K+1),
     $                  LDA )
         ELSE IF ( LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', M-K, N, K, ONE,
     $                  W(K+1,1), LDW, DWORK(PDW3), N, ONE, A(K+1,1),
     $                  LDA )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, ONE,
     $                  W(1,K+1), LDW, DWORK(PDW3), N, ONE, A(K+1,1),
     $                  LDA )
         END IF
      END IF
C
C     26) A2 := A2 + V2*DW5'
C
      IF ( M.GT.K ) THEN
         IF ( LTRA.AND.LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW5), N, V(K+1,1), LDV, ONE, A(1,K+1),
     $                  LDA )
         ELSE IF ( LTRA ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW5), N, V(1,K+1), LDV, ONE, A(1,K+1),
     $                  LDA )
         ELSE IF ( LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', M-K, N, K, ONE,
     $                  V(K+1,1), LDV, DWORK(PDW5), N, ONE, A(K+1,1),
     $                  LDA )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, ONE,
     $                  V(1,K+1), LDV, DWORK(PDW5), N, ONE, A(K+1,1),
     $                  LDA )
         END IF
      END IF
C
C     27) DW4 := DW4 + DW7
C
      CALL DAXPY( N*K, ONE, DWORK(PDW7), 1, DWORK(PDW4), 1 )
C
C     28) DW3 := DW3*W1'
C
      IF ( LCOLW ) THEN
         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, ONE,
     $               W, LDW, DWORK(PDW3), N )
      ELSE
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $               ONE, W, LDW, DWORK(PDW3), N )
      END IF
C
C     29) DW4 := DW4 + DW3
C
      CALL DAXPY( N*K, ONE, DWORK(PDW3), 1, DWORK(PDW4), 1 )
C
C     30) DW5 := DW5*V1'
C
      IF ( LCOLV ) THEN
         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, ONE,
     $               V, LDV, DWORK(PDW5), N )
      ELSE
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $               ONE, V, LDV, DWORK(PDW5), N )
      END IF
C
C     31) DW4 := DW4 + DW5
C
      CALL DAXPY( N*K, ONE, DWORK(PDW5), 1, DWORK(PDW4), 1 )
C
C     32) A1 := A1 + DW4'
C
      IF ( LA1B1 ) THEN
         IF ( LTRA ) THEN
            DO 50  I = 1, K
               CALL DAXPY( N, ONE, DWORK(PDW4+(I-1)*N), 1, A(1,I), 1 )
   50       CONTINUE
         ELSE
            DO 60  I = 1, N
               CALL DAXPY( K, ONE, DWORK(PDW4+I-1), N, A(1,I), 1 )
   60       CONTINUE
         END IF
      ELSE
         IF ( LTRA ) THEN
            DO 70  I = 1, K
               CALL DCOPY( N, DWORK(PDW4+(I-1)*N), 1, A(1,I), 1 )
   70       CONTINUE
         ELSE
            DO 80  I = 1, N
               CALL DCOPY( K, DWORK(PDW4+I-1), N, A(1,I), 1 )
   80       CONTINUE
         END IF
      END IF
C
C     Update the matrix B.
C
      IF ( LTRQ ) THEN
C
C        33) DW3 := DW1*T11
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW3), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT11), LDT, DWORK(PDW3), N )
C
C        34) DW4 := DW6*T31
C
         CALL DCOPY( N*(K-1), DWORK(PDW6), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, T(1,PT31+1), LDT, DWORK(PDW4), N )
C
C        35) DW3 := DW3 + DW4
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW4), 1, DWORK(PDW3+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ17) DW8 := DW9*T21
C
            CALL DCOPY( N*(K-1), DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                  K-1, ONE, T(1,PT21+1), LDT, DWORK(PDW8), N )
C
C           NZ18) DW3 := DW3 + DW8
C
            CALL DAXPY( N*(K-1), ONE, DWORK(PDW8), 1, DWORK(PDW3+N), 1 )
         END IF
C
C        36) DW4 := DW2*S1
C
         CALL DCOPY( N*(K-1), DWORK(PDW2), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, RS(1,PS1+1), LDRS, DWORK(PDW4), N )
C
C        37) DW3 := DW3 + DW4
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW4), 1, DWORK(PDW3+N), 1 )
C
C        38) DW4 := DW1*T12
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT12), LDT, DWORK(PDW4), N )
C
C        38) DW5 := DW6*T32
C
         CALL DCOPY( N*(K-1), DWORK(PDW6), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $               K-1, ONE, T(1,PT32+1), LDT, DWORK(PDW5), N )
C
C        40) DW4 := DW4 + DW5
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW5), 1, DWORK(PDW4+N), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ19) DW8 := DW9*T22
C
            CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                  K, ONE, T(1,PT22), LDT, DWORK(PDW8), N )
C
C           NZ20) DW4 := DW4 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW4), 1 )
         END IF
C
C        41) DW5 := DW2*S2
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PS2), LDRS, DWORK(PDW5), N )
C
C        42) DW4 := DW4 + DW5
C
         CALL DAXPY( N*K, ONE, DWORK(PDW5), 1, DWORK(PDW4), 1 )
C
C        43) DW6 := DW6*T33
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT33), LDT, DWORK(PDW6), N )
C
C        44) DW1 := DW1*T13
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT13), LDT, DWORK(PDW1), N )
C
C        45) DW6 := DW6 + DW1
C
         CALL DAXPY( N*K, ONE, DWORK(PDW1), 1, DWORK(PDW6), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ19) DW9 := DW9*T23
C
            CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N,
     $                  K, ONE, T(1,PT23), LDT, DWORK(PDW9), N )
C
C           NZ20) DW6 := DW6 + DW9
C
            CALL DAXPY( N*K, ONE, DWORK(PDW9), 1, DWORK(PDW6), 1 )
         END IF
C
C        46) DW2 := DW2*S3
C
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PS3), LDRS, DWORK(PDW2), N )
C
C        45) DW6 := DW6 + DW2
C
         CALL DAXPY( N*K, ONE, DWORK(PDW2), 1, DWORK(PDW6), 1 )
      ELSE
C
C        33) DW3 := DW1*T11'
C
         CALL DCOPY( N*K, DWORK(PDW1), 1, DWORK(PDW3), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT11), LDT, DWORK(PDW3), N )
C
C        34) DW4 := DW6*T13'
C
         CALL DCOPY( N*K, DWORK(PDW6), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT13), LDT, DWORK(PDW4), N )
C
C        35) DW3 := DW3 + DW4
C
         CALL DAXPY( N*K, ONE, DWORK(PDW4), 1, DWORK(PDW3), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ17) DW8 := DW9*T12'
C
            CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, T(1,PT12), LDT, DWORK(PDW8), N )
C
C           NZ18) DW3 := DW3 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW3), 1 )
         END IF
C
C        36) DW4 := DW2*R1'
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, RS(1,PR1), LDRS, DWORK(PDW4), N )
C
C        37) DW3 := DW3 - DW4
C
         CALL DAXPY( N*K, -ONE, DWORK(PDW4), 1, DWORK(PDW3), 1 )
C
C        38) DW4 := DW6*T23'
C
         CALL DCOPY( N*K, DWORK(PDW6), 1, DWORK(PDW4), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT23), LDT, DWORK(PDW4), N )
C
C        39) DW5 := DW1*T21'
C
         CALL DCOPY( N*(K-1), DWORK(PDW1+N), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, T(1,PT21+1), LDT, DWORK(PDW5), N )
C
C        40) DW4 := DW4 + DW5
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW5), 1, DWORK(PDW4), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ19) DW8 := DW9*T22'
C
            CALL DCOPY( N*K, DWORK(PDW9), 1, DWORK(PDW8), 1 )
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $                  ONE, T(1,PT22), LDT, DWORK(PDW8), N )
C
C           NZ20) DW4 := DW4 + DW8
C
            CALL DAXPY( N*K, ONE, DWORK(PDW8), 1, DWORK(PDW4), 1 )
         END IF
C
C        41) DW5 := DW2*R2'
C
         CALL DCOPY( N*K, DWORK(PDW2), 1, DWORK(PDW5), 1 )
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $               ONE, RS(1,PR2), LDRS, DWORK(PDW5), N )
C
C        42) DW4 := DW4 - DW5
C
         CALL DAXPY( N*K, -ONE, DWORK(PDW5), 1, DWORK(PDW4), 1 )
C
C        43) DW6 := DW6*T33'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K,
     $               ONE, T(1,PT33), LDT, DWORK(PDW6), N )
C
C        44) DW1 := DW1*T31'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, T(1,PT31+1), LDT, DWORK(PDW1+N), N )
C
C        45) DW6 := DW6 + DW1
C
         CALL DAXPY( N*(K-1), ONE, DWORK(PDW1+N), 1, DWORK(PDW6), 1 )
C
         IF ( LA1B1 ) THEN
C
C           NZ19) DW9 := DW9*T32'
C
            CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N,
     $                  K-1, ONE, T(1,PT32+1), LDT, DWORK(PDW9+N), N )
C
C           NZ20) DW6 := DW6 + DW9
C
            CALL DAXPY( N*(K-1), ONE, DWORK(PDW9+N), 1, DWORK(PDW6), 1 )
         END IF
C
C        46) DW2 := DW2*R3'
C
         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit', N, K-1,
     $               ONE, RS(1,PR3+1), LDRS, DWORK(PDW2+N), N )
C
C        45) DW6 := DW6 - DW2
C
         CALL DAXPY( N*(K-1), -ONE, DWORK(PDW2+N), 1, DWORK(PDW6), 1 )
      END IF
C
C     46) B2 := B2 + W2*DW3'
C
      IF ( M.GT.K ) THEN
         IF ( LTRB.AND.LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW3), N, W(K+1,1), LDW, ONE, B(1,K+1),
     $                  LDB )
         ELSE IF ( LTRB ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW3), N, W(1,K+1), LDW, ONE, B(1,K+1),
     $                  LDB )
         ELSE IF ( LCOLW ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', M-K, N, K, ONE,
     $                  W(K+1,1), LDW, DWORK(PDW3), N, ONE, B(K+1,1),
     $                  LDB )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, ONE,
     $                  W(1,K+1), LDW, DWORK(PDW3), N, ONE, B(K+1,1),
     $                  LDB )
         END IF
      END IF
C
C     47) B2 := B2 + V2*DW6'
C
      IF ( M.GT.K ) THEN
         IF ( LTRB.AND.LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW6), N, V(K+1,1), LDV, ONE, B(1,K+1),
     $                  LDB )
         ELSE IF ( LTRB ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', N, M-K, K, ONE,
     $                  DWORK(PDW6), N, V(1,K+1), LDV, ONE, B(1,K+1),
     $                  LDB )
         ELSE IF ( LCOLV ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', M-K, N, K, ONE,
     $                  V(K+1,1), LDV, DWORK(PDW6), N, ONE, B(K+1,1),
     $                  LDB )
         ELSE
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, ONE,
     $                  V(1,K+1), LDV, DWORK(PDW6), N, ONE, B(K+1,1),
     $                  LDB )
         END IF
      END IF
C
C     48) DW3 := DW3*W1'
C
      IF ( LCOLW ) THEN
         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, ONE,
     $               W, LDW, DWORK(PDW3), N )
      ELSE
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $               ONE, W, LDW, DWORK(PDW3), N )
      END IF
C
C     49) DW4 := DW4 + DW3
C
      CALL DAXPY( N*K, ONE, DWORK(PDW3), 1, DWORK(PDW4), 1 )
C
C     50) DW6 := DW6*V1'
C
      IF ( LCOLV ) THEN
         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, ONE,
     $               V, LDV, DWORK(PDW6), N )
      ELSE
         CALL DTRMM( 'Right', 'Upper', 'No Transpose', 'Unit', N, K,
     $               ONE, V, LDV, DWORK(PDW6), N )
      END IF
C
C     51) DW4 := DW4 + DW6
C
      CALL DAXPY( N*K, ONE, DWORK(PDW6), 1, DWORK(PDW4), 1 )
C
C     52) B1 := B1 + DW4'
C
      IF ( LA1B1 ) THEN
         IF ( LTRB ) THEN
            DO 90  I = 1, K
               CALL DAXPY( N, ONE, DWORK(PDW4+(I-1)*N), 1, B(1,I), 1 )
   90       CONTINUE
         ELSE
            DO 100  I = 1, N
               CALL DAXPY( K, ONE, DWORK(PDW4+I-1), N, B(1,I), 1 )
  100       CONTINUE
         END IF
      ELSE
         IF ( LTRB ) THEN
            DO 110  I = 1, K
               CALL DCOPY( N, DWORK(PDW4+(I-1)*N), 1, B(1,I), 1 )
  110       CONTINUE
         ELSE
            DO 120  I = 1, N
               CALL DCOPY( K, DWORK(PDW4+I-1), N, B(1,I), 1 )
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C *** Last line of MB04QC ***
      END
