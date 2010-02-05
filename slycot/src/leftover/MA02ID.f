      DOUBLE PRECISION FUNCTION MA02ID( TYP, NORM, N, A, LDA, QG,
     $                                  LDQG, DWORK )
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
C     To compute the value of the one norm, or the Frobenius norm, or
C     the infinity norm, or the element of largest absolute value
C     of a real skew-Hamiltonian matrix
C
C                   [  A   G  ]          T         T
C             X  =  [       T ],   G = -G,   Q = -Q,
C                   [  Q   A  ]
C
C     or of a real Hamiltonian matrix
C
C                   [  A   G  ]          T         T
C             X  =  [       T ],   G =  G,   Q =  Q,
C                   [  Q  -A  ]
C
C     where A, G and Q are real n-by-n matrices.
C
C     Note that for this kind of matrices the infinity norm is equal
C     to the one norm.
C
C     FUNCTION VALUE
C
C     MA02ID  DOUBLE PRECISION
C             The computed norm.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYP     CHARACTER*1
C             Specifies the type of the input matrix X:
C             = 'S':         X is skew-Hamiltonian;
C             = 'H':         X is Hamiltonian.
C
C     NORM    CHARACTER*1
C             Specifies the value to be returned in MA02ID:
C             = '1' or 'O':  one norm of X;
C             = 'F' or 'E':  Frobenius norm of X;
C             = 'I':         infinity norm of X;
C             = 'M':         max(abs(X(i,j)).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1)
C             On entry, the leading N-by-N+1 part of this array must
C             contain in columns 1:N the lower triangular part of the
C             matrix Q and in columns 2:N+1 the upper triangular part
C             of the matrix G. If TYP = 'S', the parts containing the
C             diagonal and the first supdiagonal of this array are not
C             referenced.
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= MAX(1,N).
C
C     Workspace
C
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             where LDWORK >= 2*N when NORM = '1', NORM = 'I' or
C             NORM = 'O'; otherwise, DWORK is not referenced.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLANHA).
C
C     KEYWORDS
C
C     Elementary matrix operations, Hamiltonian matrix, skew-Hamiltonian
C     matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER          NORM, TYP
      INTEGER            LDA, LDQG, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), DWORK(*), QG(LDQG,*)
C     .. Local Scalars ..
      LOGICAL            LSH
      INTEGER            I, J
      DOUBLE PRECISION   DSCL, DSUM, SCALE, SUM, TEMP, VALUE
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLANGE, DLAPY2
      EXTERNAL           DLANGE, DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL           DLASSQ
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
C
C     .. Executable Statements ..
C
      LSH = LSAME( TYP, 'S' )
C
      IF ( N.EQ.0 ) THEN
         VALUE = ZERO
C
      ELSE IF ( LSAME( NORM, 'M' ) .AND. LSH ) THEN
C
C        Find max(abs(A(i,j))).
C
         VALUE = DLANGE( 'MaxElement', N, N, A, LDA, DWORK )
         IF ( N.GT.1 ) THEN
            DO 30  J = 1, N+1
               DO 10  I = 1, J-2
                  VALUE = MAX( VALUE, ABS( QG(I,J) ) )
   10          CONTINUE
               DO 20  I = J+1, N
                  VALUE = MAX( VALUE, ABS( QG(I,J) ) )
   20          CONTINUE
   30       CONTINUE
         END IF
C
      ELSE IF ( LSAME( NORM, 'M' ) ) THEN
C
C        Find max( abs( A(i,j) ), abs( QG(i,j) ) ).
C
         VALUE = MAX( DLANGE( 'MaxElement', N, N, A, LDA, DWORK ),
     $                DLANGE( 'MaxElement', N, N+1, QG, LDQG,
     $                        DWORK ) )
C
      ELSE IF ( ( LSAME( NORM, 'O' ) .OR. ( NORM.EQ.'1' ) .OR.
     $            LSAME( NORM, 'I' ) ) .AND. LSH ) THEN
C
C        Find the column and row sums of A (in one pass).
C
         VALUE = ZERO
         DO 40 I = 1, N
            DWORK(I) = ZERO
   40    CONTINUE
C
         DO 60 J = 1, N
            SUM = ZERO
            DO 50 I = 1, N
               TEMP = ABS( A(I,J) )
               SUM  = SUM + TEMP
               DWORK(I) = DWORK(I) + TEMP
   50       CONTINUE
            DWORK(N+J) = SUM
   60    CONTINUE
C
C        Compute the maximal absolute column sum.
C
         DO 90 J = 1, N+1
            DO 70  I = 1, J-2
               TEMP = ABS( QG(I,J) )
               DWORK(I) = DWORK(I) + TEMP
               DWORK(J-1) = DWORK(J-1) + TEMP
   70       CONTINUE
            IF ( J.LT.N+1 ) THEN
               SUM = DWORK(N+J)
               DO 80  I = J+1, N
                  TEMP = ABS( QG(I,J) )
                  SUM  = SUM + TEMP
                  DWORK(N+I) = DWORK(N+I) + TEMP
   80          CONTINUE
               VALUE = MAX( VALUE, SUM )
            END IF
   90    CONTINUE
         DO 100 I = 1, N
            VALUE = MAX( VALUE, DWORK(I) )
  100    CONTINUE
C
      ELSE IF ( LSAME( NORM, 'O' ) .OR. ( NORM.EQ.'1' ) .OR.
     $          LSAME( NORM, 'I' ) ) THEN
C
C        Find the column and row sums of A (in one pass).
C
         VALUE = ZERO
         DO 110 I = 1, N
            DWORK(I) = ZERO
  110   CONTINUE
C
         DO 130 J = 1, N
            SUM = ZERO
            DO 120 I = 1, N
               TEMP = ABS( A(I,J) )
               SUM  = SUM + TEMP
               DWORK(I) = DWORK(I) + TEMP
  120       CONTINUE
            DWORK(N+J) = SUM
  130    CONTINUE
C
C        Compute the maximal absolute column sum.
C
         DO 160 J = 1, N+1
            DO 140  I = 1, J-2
               TEMP = ABS( QG(I,J) )
               DWORK(I) = DWORK(I) + TEMP
               DWORK(J-1) = DWORK(J-1) + TEMP
  140       CONTINUE
            IF ( J.GT.1 )
     $         DWORK(J-1) = DWORK(J-1) + ABS( QG(J-1,J) )
            IF ( J.LT.N+1 ) THEN
               SUM = DWORK(N+J) + ABS( QG(J,J) )
               DO 150 I = J+1, N
                  TEMP = ABS( QG(I,J) )
                  SUM  = SUM + TEMP
                  DWORK(N+I) = DWORK(N+I) + TEMP
  150          CONTINUE
               VALUE = MAX( VALUE, SUM )
            END IF
  160    CONTINUE
         DO 170 I = 1, N
            VALUE = MAX( VALUE, DWORK(I) )
  170    CONTINUE
C
      ELSE IF ( ( LSAME( NORM, 'F' ) .OR.
     $            LSAME( NORM, 'E' ) ) .AND. LSH ) THEN
C
C        Find normF(A).
C
         SCALE = ZERO
         SUM = ONE
         DO 180 J = 1, N
            CALL DLASSQ( N, A(1,J), 1, SCALE, SUM )
  180    CONTINUE
C
C        Add normF(G) and normF(Q).
C
         DO 190 J = 1, N+1
            IF ( J.GT.2 )
     $         CALL DLASSQ( J-2, QG(1,J), 1, SCALE, SUM )
            IF ( J.LT.N )
     $         CALL DLASSQ( N-J, QG(J+1,J), 1, SCALE, SUM )
  190    CONTINUE
         VALUE = SQRT( TWO )*SCALE*SQRT( SUM )
      ELSE IF ( LSAME( NORM, 'F' ) .OR. LSAME( NORM, 'E' ) ) THEN
         SCALE = ZERO
         SUM = ONE
         DO 200 J = 1, N
            CALL DLASSQ( N, A(1,J), 1, SCALE, SUM )
  200    CONTINUE
         DSCL = ZERO
         DSUM = ONE
         DO 210 J = 1, N+1
            IF ( J.GT.1 ) THEN
               CALL DLASSQ( J-2, QG(1,J), 1, SCALE, SUM )
               CALL DLASSQ( 1, QG(J-1,J), 1, DSCL, DSUM )
            END IF
            IF ( J.LT.N+1 ) THEN
               CALL DLASSQ( 1, QG(J,J), 1, DSCL, DSUM )
               CALL DLASSQ( N-J, QG(J+1,J), 1, SCALE, SUM )
            END IF
  210    CONTINUE
         VALUE = DLAPY2( SQRT( TWO )*SCALE*SQRT( SUM ),
     $                   DSCL*SQRT( DSUM ) )
      END IF
C
      MA02ID = VALUE
      RETURN
C *** Last line of MA02ID ***
      END
