      SUBROUTINE MB04WU( TRANQ1, TRANQ2, M, N, K, Q1, LDQ1, Q2, LDQ2,
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
C             value of LDWORK.
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
C     [1] Bunse-Gerstner, A.
C         Matrix factorizations for symplectic QR-like methods.
C         Linear Algebra Appl., 83, pp. 49-77, 1986.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSQ).
C
C     KEYWORDS
C
C     Elementary matrix operations, orthogonal symplectic matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANQ1, TRANQ2
      INTEGER           INFO, K, LDQ1, LDQ2, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), Q1(LDQ1,*), Q2(LDQ2,*), TAU(*)
C     .. Local Scalars ..
      LOGICAL           LTRQ1, LTRQ2
      INTEGER           I, J
      DOUBLE PRECISION  NU
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLARF, DLASET, DROT, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      LTRQ1 = LSAME( TRANQ1,'T' ) .OR. LSAME( TRANQ1,'C' )
      LTRQ2 = LSAME( TRANQ2,'T' ) .OR. LSAME( TRANQ2,'C' )
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
      ELSE IF ( LDWORK.LT.MAX( 1,M + N ) ) THEN
         DWORK(1) = DBLE( MAX( 1,M + N ) )
         INFO = -13
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04WU', -INFO )
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
C     Initialize columns K+1:N to columns of the unit matrix.
C
      DO 20 J = K + 1, N
         DO 10 I = 1, M
            Q1(I,J) = ZERO
   10    CONTINUE
         Q1(J,J) = ONE
   20 CONTINUE
      CALL DLASET( 'All', M, N-K, ZERO, ZERO, Q2(1,K+1), LDQ2 )
C
      IF ( LTRQ1.AND.LTRQ2 ) THEN
         DO 50 I = K, 1, -1
C
C           Apply F(I) to Q1(I+1:N,I:M) and Q2(I+1:N,I:M) from the
C           right.
C
            CALL DCOPY( M-I+1, Q2(I,I), LDQ2, DWORK, 1 )
            IF ( I.LT.N ) THEN
               Q1(I,I) = ONE
               CALL DLARF( 'Right', N-I, M-I+1, Q1(I,I), LDQ1, TAU(I),
     $                     Q1(I+1,I), LDQ1, DWORK(M+1) )
               CALL DLARF( 'Right', N-I, M-I+1, Q1(I,I), LDQ1, TAU(I),
     $                     Q2(I+1,I), LDQ2, DWORK(M+1) )
            END IF
            IF ( I.LT.M )
     $         CALL DSCAL( M-I, -TAU(I), Q1(I,I+1), LDQ1 )
            Q1(I,I) = ONE - TAU(I)
C
C           Set Q1(I,1:I-1) and Q2(I,1:M) to zero.
C
            DO 30 J = 1, I - 1
               Q1(I,J) = ZERO
   30       CONTINUE
            DO 40 J = 1, M
               Q2(I,J) = ZERO
   40       CONTINUE
C
C           Apply G(I) to Q1(I:N,I) and Q2(I:N,I) from the right.
C
            CALL DROT( N-I+1, Q1(I,I), 1, Q2(I,I), 1, CS(2*I-1),
     $                 CS(2*I) )
C
C           Apply H(I) to Q1(I:N,I:M) and Q2(I:N,I:M) from the right.
C
            NU = DWORK(1)
            DWORK(1) = ONE
            CALL DLARF( 'Right', N-I+1, M-I+1, DWORK, 1, NU, Q1(I,I),
     $                  LDQ1, DWORK(M+1) )
            CALL DLARF( 'Right', N-I+1, M-I+1, DWORK, 1, NU, Q2(I,I),
     $                  LDQ2, DWORK(M+1) )
   50    CONTINUE
      ELSE IF ( LTRQ1 ) THEN
         DO 80 I = K, 1, -1
C
C           Apply F(I) to Q1(I+1:N,I:M) from the right and to
C           Q2(I:M,I+1:N) from the left.
C
            CALL DCOPY( M-I+1, Q2(I,I), 1, DWORK, 1 )
            IF ( I.LT.N ) THEN
               Q1(I,I) = ONE
               CALL DLARF( 'Right', N-I, M-I+1, Q1(I,I), LDQ1, TAU(I),
     $                     Q1(I+1,I), LDQ1, DWORK(M+1) )
               CALL DLARF( 'Left', M-I+1, N-I, Q1(I,I), LDQ1, TAU(I),
     $                     Q2(I,I+1), LDQ2, DWORK(M+1) )
            END IF
            IF ( I.LT.M )
     $         CALL DSCAL( M-I, -TAU(I), Q1(I,I+1), LDQ1 )
            Q1(I,I) = ONE - TAU(I)
C
C           Set Q1(I,1:I-1) and Q2(1:M,I) to zero.
C
            DO 60 J = 1, I - 1
               Q1(I,J) = ZERO
   60       CONTINUE
            DO 70 J = 1, M
               Q2(J,I) = ZERO
   70       CONTINUE
C
C           Apply G(I) to Q1(I:N,I) from the right and to Q2(I,I:N)
C           from the left.
C
            CALL DROT( N-I+1, Q1(I,I), 1, Q2(I,I), LDQ2, CS(2*I-1),
     $                 CS(2*I) )
C
C           Apply H(I) to Q1(I:N,I:M) from the right and to Q2(I:M,I:N)
C           from the left.
C
            NU = DWORK(1)
            DWORK(1) = ONE
            CALL DLARF( 'Right', N-I+1, M-I+1, DWORK, 1, NU, Q1(I,I),
     $                  LDQ1, DWORK(M+1) )
            CALL DLARF( 'Left', M-I+1, N-I+1, DWORK, 1, NU, Q2(I,I),
     $                  LDQ2, DWORK(M+1) )
   80    CONTINUE
      ELSE IF ( LTRQ2 ) THEN
         DO 110 I = K, 1, -1
C
C           Apply F(I) to Q1(I:M,I+1:N) from the left and to
C           Q2(I+1:N,I:M) from the right.
C
            CALL DCOPY( M-I+1, Q2(I,I), LDQ2, DWORK, 1 )
            IF ( I.LT.N ) THEN
               Q1(I,I) = ONE
               CALL DLARF( 'Left', M-I+1, N-I, Q1(I,I), 1, TAU(I),
     $                     Q1(I,I+1), LDQ1, DWORK(M+1) )
               CALL DLARF( 'Right', N-I, M-I+1, Q1(I,I), 1, TAU(I),
     $                     Q2(I+1,I), LDQ2, DWORK(M+1) )
            END IF
            IF ( I.LT.M )
     $         CALL DSCAL( M-I, -TAU(I), Q1(I+1,I), 1 )
            Q1(I,I) = ONE - TAU(I)
C
C           Set Q1(1:I-1,I) and Q2(I,1:M) to zero.
C
            DO 90 J = 1, I - 1
               Q1(J,I) = ZERO
   90       CONTINUE
            DO 100 J = 1, M
               Q2(I,J) = ZERO
  100       CONTINUE
C
C           Apply G(I) to Q1(I,I:N) from the left and to Q2(I:N,I)
C           from the right.
C
            CALL DROT( N-I+1, Q1(I,I), LDQ1, Q2(I,I), 1, CS(2*I-1),
     $                 CS(2*I) )
C
C           Apply H(I) to Q1(I:M,I:N) from the left and to Q2(I:N,I:M)
C           from the left.
C
            NU = DWORK(1)
            DWORK(1) = ONE
            CALL DLARF( 'Left', M-I+1, N-I+1, DWORK, 1, NU, Q1(I,I),
     $                  LDQ1, DWORK(M+1) )
            CALL DLARF( 'Right', N-I+1, M-I+1, DWORK, 1, NU, Q2(I,I),
     $                  LDQ2, DWORK(M+1) )
  110    CONTINUE
      ELSE
         DO 140 I = K, 1, -1
C
C           Apply F(I) to Q1(I:M,I+1:N) and Q2(I:M,I+1:N) from the left.
C
            CALL DCOPY( M-I+1, Q2(I,I), 1, DWORK, 1 )
            IF ( I.LT.N ) THEN
               Q1(I,I) = ONE
               CALL DLARF( 'Left', M-I+1, N-I, Q1(I,I), 1, TAU(I),
     $                     Q1(I,I+1), LDQ1, DWORK(M+1) )
               CALL DLARF( 'Left', M-I+1, N-I, Q1(I,I), 1, TAU(I),
     $                     Q2(I,I+1), LDQ2, DWORK(M+1) )
            END IF
            IF ( I.LT.M )
     $         CALL DSCAL( M-I, -TAU(I), Q1(I+1,I), 1 )
            Q1(I,I) = ONE - TAU(I)
C
C           Set Q1(1:I-1,I) and Q2(1:M,I) to zero.
C
            DO 120 J = 1, I - 1
               Q1(J,I) = ZERO
  120       CONTINUE
            DO 130 J = 1, M
               Q2(J,I) = ZERO
  130       CONTINUE
C
C           Apply G(I) to Q1(I,I:N) and Q2(I,I:N) from the left.
C
            CALL DROT( N-I+1, Q1(I,I), LDQ1, Q2(I,I), LDQ2, CS(2*I-1),
     $                 CS(2*I) )
C
C           Apply H(I) to Q1(I:M,I:N) and Q2(I:M,I:N) from the left.
C
            NU = DWORK(1)
            DWORK(1) = ONE
            CALL DLARF( 'Left', M-I+1, N-I+1, DWORK, 1, NU, Q1(I,I),
     $                  LDQ1, DWORK(M+1) )
            CALL DLARF( 'Left', M-I+1, N-I+1, DWORK, 1, NU, Q2(I,I),
     $                  LDQ2, DWORK(M+1) )
  140    CONTINUE
      END IF
      DWORK(1) = DBLE( MAX( 1, M+N ) )
C *** Last line of MB04WU ***
      END
