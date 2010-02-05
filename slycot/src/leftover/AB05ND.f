      SUBROUTINE AB05ND( OVER, N1, M1, P1, N2, ALPHA, A1, LDA1, B1,
     $                   LDB1, C1, LDC1, D1, LDD1, A2, LDA2, B2, LDB2,
     $                   C2, LDC2, D2, LDD2, N, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, IWORK, DWORK, LDWORK, INFO )
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
C     To obtain the state-space model (A,B,C,D) for the feedback
C     inter-connection of two systems, each given in state-space form.
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
C             the order of the matrix A1.  N1 >= 0.
C
C     M1      (input) INTEGER
C             The number of input variables for the first system and the
C             number of output variables from the second system.
C             M1 >= 0.
C
C     P1      (input) INTEGER
C             The number of output variables from the first system and
C             the number of input variables for the second system.
C             P1 >= 0.
C
C     N2      (input) INTEGER
C             The number of state variables in the second system, i.e.
C             the order of the matrix A2.  N2 >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             A coefficient multiplying the transfer-function matrix
C             (or the output equation) of the second system.
C             ALPHA = +1 corresponds to positive feedback, and
C             ALPHA = -1 corresponds to negative feedback.
C
C     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1)
C             The leading N1-by-N1 part of this array must contain the
C             state transition matrix A1 for the first system.
C
C     LDA1    INTEGER
C             The leading dimension of array A1.  LDA1 >= MAX(1,N1).
C
C     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1)
C             The leading N1-by-M1 part of this array must contain the
C             input/state matrix B1 for the first system.
C
C     LDB1    INTEGER
C             The leading dimension of array B1.  LDB1 >= MAX(1,N1).
C
C     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1)
C             The leading P1-by-N1 part of this array must contain the
C             state/output matrix C1 for the first system.
C
C     LDC1    INTEGER
C             The leading dimension of array C1.
C             LDC1 >= MAX(1,P1) if N1 > 0.
C             LDC1 >= 1 if N1 = 0.
C
C     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1)
C             The leading P1-by-M1 part of this array must contain the
C             input/output matrix D1 for the first system.
C
C     LDD1    INTEGER
C             The leading dimension of array D1.  LDD1 >= MAX(1,P1).
C
C     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2)
C             The leading N2-by-N2 part of this array must contain the
C             state transition matrix A2 for the second system.
C
C     LDA2    INTEGER
C             The leading dimension of array A2.  LDA2 >= MAX(1,N2).
C
C     B2      (input) DOUBLE PRECISION array, dimension (LDB2,P1)
C             The leading N2-by-P1 part of this array must contain the
C             input/state matrix B2 for the second system.
C
C     LDB2    INTEGER
C             The leading dimension of array B2.  LDB2 >= MAX(1,N2).
C
C     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2)
C             The leading M1-by-N2 part of this array must contain the
C             state/output matrix C2 for the second system.
C
C     LDC2    INTEGER
C             The leading dimension of array C2.
C             LDC2 >= MAX(1,M1) if N2 > 0.
C             LDC2 >= 1 if N2 = 0.
C
C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1)
C             The leading M1-by-P1 part of this array must contain the
C             input/output matrix D2 for the second system.
C
C     LDD2    INTEGER
C             The leading dimension of array D2.  LDD2 >= MAX(1,M1).
C
C     N       (output) INTEGER
C             The number of state variables (N1 + N2) in the connected
C             system, i.e. the order of the matrix A, the number of rows
C             of B and the number of columns of C.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
C             The leading N-by-N part of this array contains the state
C             transition matrix A for the connected system.
C             The array A can overlap A1 if OVER = 'O'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1)
C             The leading N-by-M1 part of this array contains the
C             input/state matrix B for the connected system.
C             The array B can overlap B1 if OVER = 'O'.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
C             The leading P1-by-N part of this array contains the
C             state/output matrix C for the connected system.
C             The array C can overlap C1 if OVER = 'O'.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P1) if N1+N2 > 0.
C             LDC >= 1 if N1+N2 = 0.
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1)
C             The leading P1-by-M1 part of this array contains the
C             input/output matrix D for the connected system.
C             The array D can overlap D1 if OVER = 'O'.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P1).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (P1)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.             If OVER = 'N',
C             LDWORK >= MAX(1, P1*P1, M1*M1, N1*P1), and if OVER = 'O',
C             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*M1, N1*P1) ),
C                                                        if M1 <= N*N2;
C             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*(M1+1), N1*P1) ),
C                                                        if M1 >  N*N2.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C             > 0:  if INFO = i, 1 <= i <= P1, the system is not
C                   completely controllable. That is, the matrix
C                   (I + ALPHA*D1*D2) is exactly singular (the element
C                   U(i,i) of the upper triangular factor of LU
C                   factorization is exactly zero), possibly due to
C                   rounding errors.
C
C     METHOD
C
C     After feedback inter-connection of the two systems,
C
C     X1'     = A1*X1 + B1*U1
C     Y1      = C1*X1 + D1*U1
C
C     X2'     = A2*X2 + B2*U2
C     Y2      = C2*X2 + D2*U2
C
C     (where  '  denotes differentiation with respect to time)
C
C     the following state-space model will be obtained:
C
C     X'      = A*X  +  B*U
C     Y       = C*X  +  D*U
C
C     where       U = U1 + alpha*Y2,    X  =  ( X1 ),
C                 Y = Y1 = U2,                ( X2 )
C
C     matrix  A  has the form
C
C     ( A1  -  alpha*B1*E12*D2*C1       -  alpha*B1*E12*C2    ),
C     (        B2*E21*C1            A2  -  alpha*B2*E21*D1*C2 )
C
C     matrix  B  has the form
C
C     (  B1*E12    ),
C     (  B2*E21*D1 )
C
C     matrix  C  has the form
C
C     (  E21*C1     -  alpha*E21*D1*C2 ),
C
C     matrix D  has the form
C
C     (  E21*D1 ),
C
C     E21  =  ( I + alpha*D1*D2 )-INVERSE and
C     E12  =  ( I + alpha*D2*D1 )-INVERSE = I - alpha*D2*E21*D1.
C
C     Taking N1 = 0 and/or N2 = 0 on the routine call will solve the
C     constant plant and/or constant feedback cases.
C
C     REFERENCES
C
C     None
C
C     NUMERICAL ASPECTS
C
C     None
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB05BD by C.J.Benson, Kingston
C     Polytechnic, United Kingdom, January 1982.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
C     Feb. 2004.
C
C     KEYWORDS
C
C     Continuous-time system, multivariable system, state-space model,
C     state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO=0.0D0, ONE=1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         OVER
      INTEGER           INFO, LDA, LDA1, LDA2, LDB, LDB1, LDB2, LDC,
     $                  LDC1, LDC2, LDD, LDD1, LDD2, LDWORK, M1, N, N1,
     $                  N2, P1
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), A1(LDA1,*), A2(LDA2,*), B(LDB,*),
     $                  B1(LDB1,*), B2(LDB2,*), C(LDC,*), C1(LDC1,*),
     $                  C2(LDC2,*), D(LDD,*), D1(LDD1,*), D2(LDD2,*),
     $                  DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LOVER
      INTEGER           I, J, LDW, LDWM1
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DGETRF, DGETRS, DLACPY,
     $                  DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      LOVER = LSAME( OVER, 'O' )
      LDWM1 = MAX( 1, M1 )
      N = N1 + N2
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( .NOT.LOVER .AND. .NOT.LSAME( OVER, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N1.LT.0 ) THEN
         INFO = -2
      ELSE IF( M1.LT.0 ) THEN
         INFO = -3
      ELSE IF( P1.LT.0 ) THEN
         INFO = -4
      ELSE IF( N2.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA1.LT.MAX( 1, N1 ) ) THEN
         INFO = -8
      ELSE IF( LDB1.LT.MAX( 1, N1 ) ) THEN
         INFO = -10
      ELSE IF( ( N1.GT.0 .AND. LDC1.LT.MAX( 1, P1 ) ) .OR.
     $         ( N1.EQ.0 .AND. LDC1.LT.1 ) ) THEN
         INFO = -12
      ELSE IF( LDD1.LT.MAX( 1, P1 ) ) THEN
         INFO = -14
      ELSE IF( LDA2.LT.MAX( 1, N2 ) ) THEN
         INFO = -16
      ELSE IF( LDB2.LT.MAX( 1, N2 ) ) THEN
         INFO = -18
      ELSE IF( ( N2.GT.0 .AND. LDC2.LT.LDWM1 ) .OR.
     $         ( N2.EQ.0 .AND. LDC2.LT.1 ) ) THEN
         INFO = -20
      ELSE IF( LDD2.LT.LDWM1 ) THEN
         INFO = -22
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -25
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -27
      ELSE IF( ( N.GT.0 .AND. LDC.LT.MAX( 1, P1 ) ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -29
      ELSE IF( LDD.LT.MAX( 1, P1 ) ) THEN
         INFO = -31
      ELSE
         LDW = MAX( P1*P1, M1*M1, N1*P1 )
         IF( LOVER ) THEN
            IF( M1.GT.N*N2 )
     $         LDW = MAX( LDW, M1*( M1 + 1 ) )
            LDW = N1*P1 + LDW
         END IF
         IF( LDWORK.LT.MAX( 1, LDW ) )
     $      INFO = -34
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB05ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, MIN( M1, P1 ) ).EQ.0 )
     $   RETURN
C
      IF ( P1.GT.0 ) THEN
C
C        Form  ( I  +  alpha * D1 * D2 ).
C
         CALL DLASET( 'F', P1, P1, ZERO, ONE, DWORK, P1 )
         CALL DGEMM ( 'No transpose', 'No transpose', P1, P1, M1, ALPHA,
     $                D1, LDD1, D2, LDD2, ONE, DWORK, P1 )
C
C        Factorize this matrix.
C
         CALL DGETRF( P1, P1, DWORK, P1, IWORK, INFO )
C
         IF ( INFO.NE.0 )
     $      RETURN
C
C        Form  E21 * D1.
C
         IF ( LOVER .AND. LDD1.LE.LDD ) THEN
            IF ( LDD1.LT.LDD ) THEN
C
               DO 20 J = M1, 1, -1
                  DO 10 I = P1, 1, -1
                     D(I,J) = D1(I,J)
   10             CONTINUE
   20          CONTINUE
C
            END IF
         ELSE
            CALL DLACPY( 'F', P1, M1, D1, LDD1, D, LDD )
         END IF
C
         CALL DGETRS( 'No transpose', P1, M1, DWORK, P1, IWORK, D, LDD,
     $                INFO )
         IF ( N1.GT.0 ) THEN
C
C           Form  E21 * C1.
C
            IF ( LOVER ) THEN
C
C              First save  C1.
C
               LDW = LDW - P1*N1 + 1
               CALL DLACPY( 'F', P1, N1, C1, LDC1, DWORK(LDW), P1 )
C
               IF ( LDC1.NE.LDC )
     $            CALL DLACPY( 'F', P1, N1, DWORK(LDW), P1, C, LDC )
            ELSE
               CALL DLACPY( 'F', P1, N1, C1, LDC1, C, LDC )
            END IF
C
            CALL DGETRS( 'No transpose', P1, N1, DWORK, P1, IWORK,
     $                   C, LDC, INFO )
         END IF
C
C        Form  E12  =  I  -  alpha * D2 * ( E21 * D1 ).
C
         CALL DLASET( 'F', M1, M1, ZERO, ONE, DWORK, LDWM1 )
         CALL DGEMM ( 'No transpose', 'No transpose', M1, M1, P1,
     $                -ALPHA, D2, LDD2, D, LDD, ONE, DWORK, LDWM1 )
C
      ELSE
         CALL DLASET( 'F', M1, M1, ZERO, ONE, DWORK, LDWM1 )
      END IF
C
      IF ( LOVER .AND. LDA1.LE.LDA ) THEN
         IF ( LDA1.LT.LDA ) THEN
C
            DO 40 J = N1, 1, -1
               DO 30 I = N1, 1, -1
                  A(I,J) = A1(I,J)
   30          CONTINUE
   40       CONTINUE
C
         END IF
      ELSE
         CALL DLACPY( 'F', N1, N1, A1, LDA1, A, LDA )
      END IF
C
      IF ( N1.GT.0 .AND. M1.GT.0 ) THEN
C
C        Form  B1 * E12.
C
         IF ( LOVER ) THEN
C
C           Use the blocks (1,2) and (2,2) of A as workspace.
C
            IF ( N1*M1.LE.N*N2 ) THEN
C
C              Use BLAS 3 code.
C
               CALL DLACPY( 'F', N1, M1, B1, LDB1, A(1,N1+1), N1 )
               CALL DGEMM ( 'No transpose', 'No transpose', N1, M1, M1,
     $                      ONE, A(1,N1+1), N1, DWORK, LDWM1, ZERO, B,
     $                      LDB )
            ELSE IF ( LDB1.LT.LDB ) THEN
C
               DO 60 J = M1, 1, -1
                  DO 50 I = N1, 1, -1
                     B(I,J) = B1(I,J)
   50             CONTINUE
   60          CONTINUE
C
               IF ( M1.LE.N*N2 ) THEN
C
C                 Use BLAS 2 code.
C
                  DO 70 J = 1, N1
                     CALL DCOPY( M1, B(J,1), LDB, A(1,N1+1), 1 )
                     CALL DGEMV( 'Transpose', M1, M1, ONE, DWORK, LDWM1,
     $                           A(1,N1+1), 1, ZERO, B(J,1), LDB )
   70             CONTINUE
C
               ELSE
C
C                 Use additional workspace.
C
                  DO 80 J = 1, N1
                     CALL DCOPY( M1, B(J,1), LDB, DWORK(M1*M1+1), 1 )
                     CALL DGEMV( 'Transpose', M1, M1, ONE, DWORK, LDWM1,
     $                           DWORK(M1*M1+1), 1, ZERO, B(J,1), LDB )
   80             CONTINUE
C
               END IF
C
            ELSE IF ( M1.LE.N*N2 ) THEN
C
C              Use BLAS 2 code.
C
               DO 90 J = 1, N1
                  CALL DCOPY( M1, B1(J,1), LDB1, A(1,N1+1), 1 )
                  CALL DGEMV( 'Transpose', M1, M1, ONE, DWORK, LDWM1,
     $                        A(1,N1+1), 1, ZERO, B(J,1), LDB )
   90          CONTINUE
C
            ELSE
C
C              Use additional workspace.
C
               DO 100 J = 1, N1
                  CALL DCOPY( M1, B1(J,1), LDB1, DWORK(M1*M1+1), 1 )
                  CALL DGEMV( 'Transpose', M1, M1, ONE, DWORK, LDWM1,
     $                        DWORK(M1*M1+1), 1, ZERO, B(J,1), LDB )
  100          CONTINUE
C
            END IF
         ELSE
            CALL DGEMM ( 'No transpose', 'No transpose', N1, M1, M1,
     $                   ONE, B1, LDB1, DWORK, LDWM1, ZERO, B, LDB )
         END IF
      END IF
C
      IF ( N2.GT.0 ) THEN
C
C        Complete matrices  B  and  C.
C
         IF ( P1.GT.0 ) THEN
            CALL DGEMM ( 'No transpose', 'No transpose', N2, M1, P1,
     $                   ONE, B2, LDB2, D, LDD, ZERO, B(N1+1,1), LDB )
            CALL DGEMM ( 'No transpose', 'No transpose', P1, N2, M1,
     $                   -ALPHA, D, LDD, C2, LDC2, ZERO, C(1,N1+1), LDC
     $                 )
         ELSE IF ( M1.GT.0 ) THEN
            CALL DLASET( 'F', N2, M1, ZERO, ZERO, B(N1+1,1), LDB )
         END IF
      END IF
C
      IF ( N1.GT.0 .AND. P1.GT.0 ) THEN
C
C        Form upper left quadrant of  A.
C
         CALL DGEMM ( 'No transpose', 'No transpose', N1, P1, M1,
     $                -ALPHA, B, LDB, D2, LDD2, ZERO, DWORK, N1 )
C
         IF ( LOVER ) THEN
            CALL DGEMM ( 'No transpose', 'No transpose', N1, N1, P1,
     $                   ONE, DWORK, N1, DWORK(LDW), P1, ONE, A, LDA )
         ELSE
            CALL DGEMM ( 'No transpose', 'No transpose', N1, N1, P1,
     $                   ONE, DWORK, N1, C1, LDC1, ONE, A, LDA )
         END IF
      END IF
C
      IF ( N2.GT.0 ) THEN
C
C        Form lower right quadrant of  A.
C
         CALL DLACPY( 'F', N2, N2, A2, LDA2, A(N1+1,N1+1), LDA )
         IF ( M1.GT.0 )
     $      CALL DGEMM ( 'No transpose', 'No transpose', N2, N2, M1,
     $                   -ALPHA, B(N1+1,1), LDB, C2, LDC2, ONE,
     $                   A(N1+1,N1+1), LDA )
C
C        Complete the matrix  A.
C
         CALL DGEMM ( 'No transpose', 'No transpose', N2, N1, P1,
     $                ONE, B2, LDB2, C, LDC, ZERO, A(N1+1,1), LDA )
         CALL DGEMM ( 'No transpose', 'No transpose', N1, N2, M1,
     $                -ALPHA, B, LDB, C2, LDC2, ZERO, A(1,N1+1), LDA )
      END IF
C
      RETURN
C *** Last line of AB05ND ***
      END
