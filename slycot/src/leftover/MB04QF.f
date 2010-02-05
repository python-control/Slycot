      SUBROUTINE MB04QF( DIRECT, STOREV, STOREW, N, K, V, LDV, W, LDW,
     $                   CS, TAU, RS, LDRS, T, LDT, DWORK )
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
C     To form the triangular block factors R, S and T of a symplectic
C     block reflector SH, which is defined as a product of 2k
C     concatenated Householder reflectors and k Givens rotators,
C
C         SH = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )
C              diag( H(2),H(2) ) G(2) diag( F(2),F(2) )
C                                ....
C              diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).
C
C     The upper triangular blocks of the matrices
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
C     DIRECT  CHARACTER*1
C             This is a dummy argument, which is reserved for future
C             extensions of this subroutine. Not referenced.
C
C     STOREV  CHARACTER*1
C             Specifies how the vectors which define the concatenated
C             Householder F(i) reflectors are stored:
C             = 'C':  columnwise;
C             = 'R':  rowwise.
C
C     STOREW  CHARACTER*1
C             Specifies how the vectors which define the concatenated
C             Householder H(i) reflectors are stored:
C             = 'C':  columnwise;
C             = 'R':  rowwise.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the Householder reflectors F(i) and H(i).
C             N >= 0.
C
C     K       (input) INTEGER
C             The number of Givens rotators.  K >= 1.
C
C     V       (input) DOUBLE PRECISION array, dimension
C                     (LDV,K) if STOREV = 'C',
C                     (LDV,N) if STOREV = 'R'
C             On entry with STOREV = 'C', the leading N-by-K part of
C             this array must contain in its i-th column the vector
C             which defines the elementary reflector F(i).
C             On entry with STOREV = 'R', the leading K-by-N part of
C             this array must contain in its i-th row the vector
C             which defines the elementary reflector F(i).
C
C     LDV     INTEGER
C             The leading dimension of the array V.
C             LDV >= MAX(1,N),  if STOREV = 'C';
C             LDV >= K,         if STOREV = 'R'.
C
C     W       (input) DOUBLE PRECISION array, dimension
C                     (LDW,K) if STOREW = 'C',
C                     (LDW,N) if STOREW = 'R'
C             On entry with STOREW = 'C', the leading N-by-K part of
C             this array must contain in its i-th column the vector
C             which defines the elementary reflector H(i).
C             On entry with STOREV = 'R', the leading K-by-N part of
C             this array must contain in its i-th row the vector
C             which defines the elementary reflector H(i).
C
C     LDW     INTEGER
C             The leading dimension of the array W.
C             LDW >= MAX(1,N),  if STOREW = 'C';
C             LDW >= K,         if STOREW = 'R'.
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
C     RS      (output) DOUBLE PRECISION array, dimension (K,6*K)
C             On exit, the leading K-by-6*K part of this array contains
C             the upper triangular matrices defining the factors R and
C             S of the symplectic block reflector SH. The (strictly)
C             lower portions of this array are not used.
C
C     LDRS    INTEGER
C             The leading dimension of the array RS.  LDRS >= K.
C
C     T       (output) DOUBLE PRECISION array, dimension (K,9*K)
C             On exit, the leading K-by-9*K part of this array contains
C             the upper triangular matrices defining the factor T of the
C             symplectic block reflector SH. The (strictly) lower
C             portions of this array are not used.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= K.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (3*K)
C
C     REFERENCES
C
C     [1] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires ( 4*K - 2 )*K*N + 19/3*K*K*K + 1/2*K*K
C     + 43/6*K - 4 floating point operations.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAEST).
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
      CHARACTER         DIRECT, STOREV, STOREW
      INTEGER           K, LDRS, LDT, LDV, LDW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), RS(LDRS,*), T(LDT,*),
     $                  TAU(*), V(LDV,*), W(LDW,*)
C     .. Local Scalars ..
      LOGICAL           LCOLV, LCOLW
      INTEGER           I, J, K2, PR1, PR2, PR3, PS1, PS2, PS3, PT11,
     $                  PT12, PT13, PT21, PT22, PT23, PT31, PT32, PT33
      DOUBLE PRECISION  CM1, TAUI, VII, WII
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, DTRMV
C
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      LCOLV = LSAME( STOREV, 'C' )
      LCOLW = LSAME( STOREW, 'C' )
C
      K2  = K + K
      PR1 = 0
      PR2 = PR1 + K
      PR3 = PR2 + K
      PS1 = PR3 + K
      PS2 = PS1 + K
      PS3 = PS2 + K
C
      PT11 = 0
      PT12 = PT11 + K
      PT13 = PT12 + K
      PT21 = PT13 + K
      PT22 = PT21 + K
      PT23 = PT22 + K
      PT31 = PT23 + K
      PT32 = PT31 + K
      PT33 = PT32 + K
C
      DO 90  I = 1, K
         TAUI = TAU(I)
         VII = V(I,I)
         V(I,I) = ONE
         WII = W(I,I)
         W(I,I) = ONE
         IF ( WII.EQ.ZERO ) THEN
            DO 10  J = 1, I
               T(J,PT11+I) = ZERO
   10       CONTINUE
            DO 20  J = 1, I-1
               T(J,PT21+I) = ZERO
   20       CONTINUE
            DO 30  J = 1, I-1
               T(J,PT31+I) = ZERO
   30       CONTINUE
            DO 40  J = 1, I-1
               RS(J,PS1+I) = ZERO
   40       CONTINUE
         ELSE
C
C           Treat first Householder reflection.
C
            IF ( LCOLV.AND.LCOLW ) THEN
C
C              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -WII, W(I,1), LDW,
     $                     W(I,I), 1, ZERO, DWORK, 1 )
C
C              Compute t2 = -wii * V(i:n,1:i-1)' * W(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -WII, V(I,1), LDV,
     $                     W(I,I), 1, ZERO, DWORK(K+1), 1 )
            ELSE IF ( LCOLV ) THEN
C
C              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -WII, W(1,I),
     $                     LDW, W(I,I), LDW, ZERO, DWORK, 1 )
C
C              Compute t2 = -wii * V(i:n,1:i-1)' * W(i,i:n)'.
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -WII, V(I,1), LDV,
     $                     W(I,I), LDW, ZERO, DWORK(K+1), 1 )
            ELSE IF ( LCOLW ) THEN
C
C              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -WII, W(I,1), LDW,
     $                     W(I,I), 1, ZERO, DWORK, 1 )
C
C              Compute t2 = -wii * V(1:i-1,i:n) * W(i:n,i).
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -WII, V(1,I),
     $                     LDV, W(I,I), 1, ZERO, DWORK(K+1), 1 )
            ELSE
C
C              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -WII, W(1,I),
     $                     LDW, W(I,I), LDW, ZERO, DWORK, 1 )
C
C              Compute t2 = -wii * V(1:i-1,i:n) * W(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -WII, V(1,I),
     $                     LDV, W(I,I), LDW, ZERO, DWORK(K+1), 1 )
            END IF
C
C           T11(1:i-1,i) := T11(1:i-1,1:i-1)*t1 + T13(1:i-1,1:i-1)*t2
C
            CALL DCOPY( I-1, DWORK, 1, T(1,PT11+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT11+1), LDT, T(1,PT11+I), 1 )
            CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT13+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT13+1), LDT, T(1,PT13+I), 1 )
            CALL DAXPY( I-1, ONE, T(1,PT13+I), 1, T(1,PT11+I), 1 )
            T(I,PT11+I) = -WII
C
            IF ( I.GT.1 ) THEN
C
C              T21(1:i-1,i) := T21(1:i-1,1:i-1)*t1 + T23(1:i-1,1:i-1)*t2
C
               CALL DCOPY( I-2, DWORK(2), 1, T(1,PT21+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-2,
     $                     T(1,PT21+2), LDT, T(1,PT21+I), 1 )
               T(I-1, PT21+I) = ZERO
               CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT23+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                     T(1,PT23+1), LDT, T(1,PT23+I), 1 )
               CALL DAXPY( I-1, ONE, T(1,PT23+I), 1, T(1,PT21+I), 1 )
C
C              T31(1:i-1,i) := T31(1:i-1,1:i-1)*t1 + T33(1:i-1,1:i-1)*t2
C
               CALL DCOPY( I-2, DWORK(2), 1, T(1,PT31+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-2,
     $                     T(1,PT31+2), LDT, T(1,PT31+I), 1 )
               T(I-1, PT31+I) = ZERO
               CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT33+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                     T(1,PT33+1), LDT, T(1,PT33+I), 1 )
               CALL DAXPY( I-1, ONE, T(1,PT33+I), 1, T(1,PT31+I), 1 )
C
C              S1(1:i-1,i) := S1(1:i-1,1:i-1)*t1 + S3(1:i-1,1:i-1)*t2
C
               CALL DCOPY( I-2, DWORK(2), 1, RS(1,PS1+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-2,
     $                     RS(1,PS1+2), LDRS, RS(1,PS1+I), 1 )
               RS(I-1, PS1+I) = ZERO
               CALL DCOPY( I-1, DWORK(K+1), 1, RS(1,PS3+I), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                     RS(1,PS3+1), LDRS, RS(1,PS3+I), 1 )
               CALL DAXPY( I-1, ONE, RS(1,PS3+I), 1, RS(1,PS1+I), 1 )
            END IF
         END IF
C
C        Treat Givens rotation.
C
         CM1 = CS(2*I-1) - ONE
         IF ( LCOLW ) THEN
            CALL DCOPY( I, W(I,1), LDW, DWORK, 1 )
         ELSE
            CALL DCOPY( I, W(1,I), 1, DWORK, 1 )
         END IF
         IF ( LCOLV ) THEN
            CALL DCOPY( I-1, V(I,1), LDV, DWORK(K+1), 1 )
         ELSE
            CALL DCOPY( I-1, V(1,I), 1, DWORK(K+1), 1 )
         END IF
C
C        R1(1:i,i) = T11(1:i,1:i) * dwork(1:i)
C                    + [ T13(1:i-1,1:i-1) * dwork(k+1:k+i-1); 0 ]
C
         CALL DCOPY( I, DWORK, 1, RS(1,PR1+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I,
     $               T(1,PT11+1), LDT, RS(1,PR1+I), 1 )
         CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT13+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               T(1,PT13+1), LDT, T(1,PT13+I), 1 )
         CALL DAXPY( I-1, ONE, T(1,PT13+I), 1, RS(1,PR1+I), 1 )
C
C        R2(1:i-1,i) = T21(1:i-1,2:i) * W(i,2:i)
C                      + T23(1:i-1,1:i-1) * V(i,1:i-1)
C
         CALL DCOPY( I-1, DWORK(2), 1, RS(1,PR2+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               T(1,PT21+2), LDT, RS(1,PR2+I), 1 )
         CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT23+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               T(1,PT23+1), LDT, T(1,PT23+I), 1 )
         CALL DAXPY( I-1, ONE, T(1,PT23+I), 1, RS(1,PR2+I), 1 )
C
C        R3(1:i-1,i) = T31(1:i-1,2:i) * dwork(2:i)
C                      + T33(1:i-1,1:i-1) * dwork(k+1:k+i-1)
C
         CALL DCOPY( I-1, DWORK(2), 1, RS(1,PR3+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               T(1,PT31+2), LDT, RS(1,PR3+I), 1 )
         CALL DCOPY( I-1, DWORK(K+1), 1, T(1,PT33+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               T(1,PT33+1), LDT, T(1,PT33+I), 1 )
         CALL DAXPY( I-1, ONE, T(1,PT33+I), 1, RS(1,PR3+I), 1 )
C
C        S2(1:i-1,i) = S1(1:i-1,2:i) * dwork(2:i)
C                      + S3(1:i-1,1:i-1) * dwork(k+1:k+i-1)
C
         CALL DCOPY( I-1, DWORK(2), 1, RS(1,PS2+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               RS(1,PS1+2), LDRS, RS(1,PS2+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               RS(1,PS3+1), LDRS, DWORK(K+1), 1 )
         CALL DAXPY( I-1, ONE, DWORK(K+1), 1, RS(1,PS2+I), 1 )
         RS(I,PS2+I) = -CS(2*I)
C
C        T12(1:i,i) = [ R1(1:i-1,1:i-1)*S2(1:i-1,i); 0 ]
C                     + (c-1) * R1(1:i,i)
C
         CALL DCOPY( I-1, RS(1,PS2+I), 1, T(1,PT12+I), 1 )
         CALL DSCAL( I-1, CM1, RS(1,PS2+I), 1)
         CALL DSCAL( I-1, CS(2*I), T(1,PT12+I), 1 )
         CALL DCOPY( I-1, T(1,PT12+I), 1, T(1,PT22+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $               RS(1,PR1+1), LDRS, T(1,PT12+I), 1 )
         T(I,PT12+I) = ZERO
         CALL DAXPY( I, CM1, RS(1,PR1+I), 1, T(1,PT12+I), 1 )
C
C        T22(1:i-1,i) = R2(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R2(1:i-1,i)
C
         IF (I.GT.1)
     $      CALL DCOPY( I-2, T(2,PT22+I), 1, T(1,PT32+I), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Unit diagonal', I-1,
     $               RS(1,PR2+1), LDRS, T(1,PT22+I), 1 )
         CALL DAXPY( I-1, CM1, RS(1,PR2+I), 1, T(1,PT22+I), 1 )
         T(I,PT22+I) = CM1
C
C        T32(1:i-1,i) = R3(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R3(1:i-1,i)
C
         IF ( I.GT.1 ) THEN
            CALL DTRMV( 'Upper', 'No transpose', 'Non-Unit', I-2,
     $                  RS(1,PR3+2), LDRS, T(1,PT32+I), 1 )
            T(I-1,PT32+I) = ZERO
            CALL DAXPY( I-1, CM1, RS(1,PR3+I), 1, T(1,PT32+I), 1 )
         END IF
C
         IF ( TAUI.EQ.ZERO ) THEN
            DO 50  J = 1, I
               T(J,PT13+I) = ZERO
   50       CONTINUE
            DO 60  J = 1, I
               T(J,PT23+I) = ZERO
   60       CONTINUE
            DO 70  J = 1, I
               T(J,PT33+I) = ZERO
   70       CONTINUE
            DO 80  J = 1, I
               RS(J,PS3+I) = ZERO
   80       CONTINUE
         ELSE
C
C           Treat second Householder reflection.
C
            IF ( LCOLV.AND.LCOLW ) THEN
C
C              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I, -TAUI, W(I,1),
     $                     LDW, V(I,I), 1, ZERO, DWORK, 1 )
C
C              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -TAUI, V(I,1),
     $                     LDV, V(I,I), 1, ZERO, DWORK(K2+1), 1 )
            ELSE IF ( LCOLV ) THEN
C
C              Compute t1 = -tau(i) * W(1:i,i:n) * V(i:n,i).
C
               CALL DGEMV( 'No Transpose', I, N-I+1, -TAUI, W(1,I),
     $                     LDW, V(I,I), 1, ZERO, DWORK, 1 )
C
C              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i).
C
               CALL DGEMV( 'Transpose', N-I+1, I-1, -TAUI, V(I,1),
     $                     LDV, V(I,I), 1, ZERO, DWORK(K2+1), 1 )
            ELSE IF ( LCOLW ) THEN
C
C              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i,i:n)'.
C
               CALL DGEMV( 'Transpose', N-I+1, I, -TAUI, W(I,1),
     $                     LDW, V(I,I), LDV, ZERO, DWORK, 1 )
C
C              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -TAUI, V(1,I),
     $                     LDV, V(I,I), LDV, ZERO, DWORK(K2+1), 1 )
            ELSE
C
C              Compute t1 = -tau(i) * W(1:i,i:n) * V(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I, N-I+1, -TAUI, W(1,I),
     $                     LDW, V(I,I), LDV, ZERO, DWORK, 1 )
C
C              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'.
C
               CALL DGEMV( 'No Transpose', I-1, N-I+1, -TAUI, V(1,I),
     $                     LDV, V(I,I), LDV, ZERO, DWORK(K2+1), 1 )
            END IF
C
C           T13(1:i,i) := T11(1:i,1:i)*t1 - tau(i)*T12(1:i,i)
C                                         + [T13(1:i-1,1:i-1)*t2;0]
C
            CALL DCOPY( I-1, DWORK(K2+1), 1, T(1,PT13+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT13+1), LDT, T(1,PT13+I), 1 )
            T(I,PT13+I) = ZERO
            CALL DCOPY( I, DWORK, 1, DWORK(K+1), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I,
     $                  T(1,PT11+1), LDT, DWORK(K+1), 1 )
            CALL DAXPY( I, ONE, DWORK(K+1), 1, T(1,PT13+I), 1 )
            CALL DAXPY( I, -TAUI, T(1,PT12+I), 1, T(1,PT13+I), 1 )
C
C           T23(1:i,i) := T21(1:i,1:i)*t1 - tau(i)*T22(1:i,i)
C                                         + [T23(1:i-1,1:i-1)*t2;0]
C
            CALL DCOPY( I-1, DWORK(K2+1), 1, T(1,PT23+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT23+1), LDT, T(1,PT23+I), 1 )
            T(I,PT23+I) = ZERO
            CALL DCOPY( I-1, DWORK(2), 1, DWORK(K+1), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT21+2), LDT, DWORK(K+1), 1 )
            CALL DAXPY( I-1, ONE, DWORK(K+1), 1, T(1,PT23+I), 1 )
            CALL DAXPY( I, -TAUI, T(1,PT22+I), 1, T(1,PT23+I), 1 )
C
C           T33(1:i,i) := T31(1:i,1:i)*t1 - tau(i)*T32(1:i,i)
C                                         + [T33(1:i-1,1:i-1)*t2;0]
C
            CALL DCOPY( I-1, DWORK(K2+1), 1, T(1,PT33+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT33+1), LDT, T(1,PT33+I), 1 )
            CALL DCOPY( I-1, DWORK(2), 1, DWORK(K+1), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  T(1,PT31+2), LDT, DWORK(K+1), 1 )
            CALL DAXPY( I-1, ONE, DWORK(K+1), 1, T(1,PT33+I), 1 )
            CALL DAXPY( I-1, -TAUI, T(1,PT32+I), 1, T(1,PT33+I), 1 )
            T(I,PT33+I) = -TAUI
C
C           S3(1:i,i) := S1(1:i,1:i)*t1 - tau(i)*S2(1:i,i)
C                                       + [S3(1:i-1,1:i-1)*t2;0]
C
            CALL DCOPY( I-1, DWORK(K2+1), 1, RS(1,PS3+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  RS(1,PS3+1), LDRS, RS(1,PS3+I), 1 )
            CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
     $                  RS(1,PS1+2), LDRS, DWORK(2), 1 )
            CALL DAXPY( I-1, ONE, DWORK(2), 1, RS(1,PS3+I), 1 )
            RS(I,PS3+I) = ZERO
            CALL DAXPY( I, -TAUI, RS(1,PS2+I), 1, RS(1,PS3+I), 1 )
         END IF
         W(I,I) = WII
         V(I,I) = VII
   90 CONTINUE
C
      RETURN
C *** Last line of MB04QF ***
      END
