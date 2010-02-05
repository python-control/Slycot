      SUBROUTINE AB05PD( OVER, N1, M, P, N2, ALPHA, A1, LDA1, B1, LDB1,
     $                   C1, LDC1, D1, LDD1, A2, LDA2, B2, LDB2, C2,
     $                   LDC2, D2, LDD2, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, INFO)
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
C     To compute the state-space model G = (A,B,C,D) corresponding to
C     the sum G = G1 + alpha*G2, where G1 = (A1,B1,C1,D1) and
C     G2 = (A2,B2,C2,D2).  G, G1, and G2 are the transfer-function
C     matrices of the corresponding state-space models.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     OVER    CHARACTER*1
C             Indicates whether the user wishes to overlap pairs of
C             arrays, as follows:
C             = 'N':  Do not overlap;
C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
C                     C1 and C, and D1 and D, i.e. the same name is
C                     effectively used for each pair (for all pairs)
C                     in the routine call.  In this case, setting
C                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD
C                     will give maximum efficiency.
C
C     Input/Output Parameters
C
C     N1      (input) INTEGER
C             The number of state variables in the first system, i.e.
C             the order of the matrix A1, the number of rows of B1 and
C             the number of columns of C1.  N1 >= 0.
C
C     M       (input) INTEGER
C             The number of input variables of the two systems, i.e. the
C             number of columns of matrices B1, D1, B2 and D2.  M >= 0.
C
C     P       (input) INTEGER
C             The number of output variables of the two systems, i.e.
C             the number of rows of matrices C1, D1, C2 and D2.  P >= 0.
C
C     N2      (input) INTEGER
C             The number of state variables in the second system, i.e.
C             the order of the matrix A2, the number of rows of B2 and
C             the number of columns of C2.  N2 >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The coefficient multiplying G2.
C
C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
C             The leading N1-by-N1 part of this array must contain the
C             state transition matrix A1 for the first system.
C
C     LDA1    INTEGER
C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
C
C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M)
C             The leading N1-by-M part of this array must contain the
C             input/state matrix B1 for the first system.
C
C     LDB1    INTEGER
C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
C
C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
C             The leading P-by-N1 part of this array must contain the
C             state/output matrix C1 for the first system.
C
C     LDC1    INTEGER
C             The leading dimension of array C1.
C             LDC1 >= MAX(1,P) if N1 > 0.
C             LDC1 >= 1 if N1 = 0.
C
C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M)
C             The leading P-by-M part of this array must contain the
C             input/output matrix D1 for the first system.
C
C     LDD1    INTEGER
C             The leading dimension of array D1.  LDD1 >= MAX(1,P).
C
C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
C             The leading N2-by-N2 part of this array must contain the
C             state transition matrix A2 for the second system.
C
C     LDA2    INTEGER
C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
C
C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M)
C             The leading N2-by-M part of this array must contain the
C             input/state matrix B2 for the second system.
C
C     LDB2    INTEGER
C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
C
C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
C             The leading P-by-N2 part of this array must contain the
C             state/output matrix C2 for the second system.
C
C     LDC2    INTEGER
C             The leading dimension of array C2.
C             LDC2 >= MAX(1,P) if N2 > 0.
C             LDC2 >= 1 if N2 = 0.
C
C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M)
C             The leading P-by-M part of this array must contain the
C             input/output matrix D2 for the second system.
C
C     LDD2    INTEGER
C             The leading dimension of array D2.  LDD2 >= MAX(1,P).
C
C     N       (output) INTEGER
C             The number of state variables (N1 + N2) in the resulting
C             system, i.e. the order of the matrix A, the number of rows
C             of B and the number of columns of C.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
C             The leading N-by-N part of this array contains the state
C             transition matrix A for the resulting system.
C             The array A can overlap A1 if OVER = 'O'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array contains the
C             input/state matrix B for the resulting system.
C             The array B can overlap B1 if OVER = 'O'.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
C             The leading P-by-N part of this array contains the
C             state/output matrix C for the resulting system.
C             The array C can overlap C1 if OVER = 'O'.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P) if N1+N2 > 0.
C             LDC >= 1 if N1+N2 = 0.
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array contains the
C             input/output matrix D for the resulting system.
C             The array D can overlap D1 if OVER = 'O'.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
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
C     The matrices of the resulting systems are determined as:
C
C           ( A1   0  )             ( B1 )
C       A = (         ) ,       B = (    ) ,
C           ( 0    A2 )             ( B2 )
C
C       C = ( C1  alpha*C2 ) ,  D = D1 + alpha*D2 .
C
C     REFERENCES
C
C     None
C
C     NUMERICAL ASPECTS
C
C     None
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Research Establishment,
C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
C     Belgium, Nov. 1996.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
C     Feb. 2004.
C
C     KEYWORDS
C
C     Multivariable system, state-space model, state-space
C     representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO=0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         OVER
      INTEGER           INFO, LDA, LDA1, LDA2, LDB, LDB1, LDB2, LDC,
     $                  LDC1, LDC2, LDD, LDD1, LDD2, M, N, N1, N2, P
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), A1(LDA1,*), A2(LDA2,*), B(LDB,*),
     $                  B1(LDB1,*), B2(LDB2,*), C(LDC,*), C1(LDC1,*),
     $                  C2(LDC2,*), D(LDD,*), D1(LDD1,*), D2(LDD2,*)
C     .. Local Scalars ..
      LOGICAL           LOVER
      INTEGER           I, J, N1P1
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLACPY, DLASCL, DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      LOVER = LSAME( OVER, 'O' )
      N = N1 + N2
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( .NOT.LOVER .AND. .NOT.LSAME( OVER, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N1.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( N2.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA1.LT.MAX( 1, N1 ) ) THEN
         INFO = -8
      ELSE IF( LDB1.LT.MAX( 1, N1 ) ) THEN
         INFO = -10
      ELSE IF( ( N1.GT.0 .AND. LDC1.LT.MAX( 1, P ) ) .OR.
     $         ( N1.EQ.0 .AND. LDC1.LT.1 ) ) THEN
         INFO = -12
      ELSE IF( LDD1.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDA2.LT.MAX( 1, N2 ) ) THEN
         INFO = -16
      ELSE IF( LDB2.LT.MAX( 1, N2 ) ) THEN
         INFO = -18
      ELSE IF( ( N2.GT.0 .AND. LDC2.LT.MAX( 1, P ) ) .OR.
     $         ( N2.EQ.0 .AND. LDC2.LT.1 ) ) THEN
         INFO = -20
      ELSE IF( LDD2.LT.MAX( 1, P ) ) THEN
         INFO = -22
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -25
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -27
      ELSE IF( ( N.GT.0 .AND. LDC.LT.MAX( 1, P ) ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -29
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -31
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB05PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, MIN( M, P ) ).EQ.0 )
     $   RETURN
C
      N1P1 = N1 + 1
C
C                       ( A1   0  )
C     Construct     A = (         ) .
C                       ( 0    A2 )
C
      IF ( LOVER .AND. LDA1.LE.LDA ) THEN
         IF ( LDA1.LT.LDA ) THEN
C
            DO 20 J = N1, 1, -1
               DO 10 I = N1, 1, -1
                  A(I,J) = A1(I,J)
   10          CONTINUE
   20       CONTINUE
C
         END IF
      ELSE
         CALL DLACPY( 'F', N1, N1, A1, LDA1, A, LDA )
      END IF
C
      IF ( N2.GT.0 ) THEN
         CALL DLASET( 'F', N1, N2, ZERO, ZERO, A(1,N1P1), LDA )
         CALL DLASET( 'F', N2, N1, ZERO, ZERO, A(N1P1,1), LDA )
         CALL DLACPY( 'F', N2, N2, A2, LDA2, A(N1P1,N1P1), LDA )
      END IF
C
C                        ( B1 )
C     Construct      B = (    ) .
C                        ( B2 )
C
      IF ( LOVER .AND. LDB1.LE.LDB ) THEN
         IF ( LDB1.LT.LDB ) THEN
C
            DO 40 J = M, 1, -1
               DO 30 I = N1, 1, -1
                  B(I,J) = B1(I,J)
   30          CONTINUE
   40       CONTINUE
C
         END IF
      ELSE
         CALL DLACPY( 'F', N1, M, B1, LDB1, B, LDB )
      END IF
C
      IF ( N2.GT.0 )
     $   CALL DLACPY( 'F', N2, M, B2, LDB2, B(N1P1,1), LDB )
C
C     Construct      C = ( C1 alpha*C2 ) .
C
      IF ( LOVER .AND. LDC1.LE.LDC ) THEN
         IF ( LDC1.LT.LDC ) THEN
C
            DO 60 J = N1, 1, -1
               DO 50 I = P, 1, -1
                  C(I,J) = C1(I,J)
   50          CONTINUE
   60       CONTINUE
C
         END IF
      ELSE
         CALL DLACPY( 'F', P, N1, C1, LDC1, C, LDC )
      END IF
C
      IF ( N2.GT.0 ) THEN
         CALL DLACPY( 'F', P, N2, C2, LDC2, C(1,N1P1), LDC )
         IF ( ALPHA.NE.ONE )
     $      CALL DLASCL( 'G', 0, 0, ONE, ALPHA, P, N2, C(1,N1P1), LDC,
     $                   INFO )
      END IF
C
C     Construct       D = D1 + alpha*D2 .
C
      IF ( LOVER .AND. LDD1.LE.LDD ) THEN
         IF ( LDD1.LT.LDD ) THEN
C
            DO 80 J = M, 1, -1
               DO 70 I = P, 1, -1
                  D(I,J) = D1(I,J)
   70          CONTINUE
   80       CONTINUE
C
         END IF
      ELSE
         CALL DLACPY( 'F', P, M, D1, LDD1, D, LDD )
      END IF
C
      DO 90 J = 1, M
         CALL DAXPY( P, ALPHA, D2(1,J), 1, D(1,J), 1 )
   90 CONTINUE
C
      RETURN
C *** Last line of AB05PD ***
      END
