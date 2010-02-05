      SUBROUTINE AB05MD( UPLO, OVER, N1, M1, P1, N2, P2, A1, LDA1, B1,
     $                   LDB1, C1, LDC1, D1, LDD1, A2, LDA2, B2, LDB2,
     $                   C2, LDC2, D2, LDD2, N, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, DWORK, LDWORK, INFO )
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
C     To obtain the state-space model (A,B,C,D) for the cascaded
C     inter-connection of two systems, each given in state-space form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Indicates whether the user wishes to obtain the matrix A
C             in the upper or lower block diagonal form, as follows:
C             = 'U':  Obtain A in the upper block diagonal form;
C             = 'L':  Obtain A in the lower block diagonal form.
C
C     OVER    CHARACTER*1
C             Indicates whether the user wishes to overlap pairs of
C             arrays, as follows:
C             = 'N':  Do not overlap;
C             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B,
C                     C1 and C, and D1 and D (for UPLO = 'L'), or A2
C                     and A, B2 and B, C2 and C, and D2 and D (for
C                     UPLO = 'U'), i.e. the same name is effectively
C                     used for each pair (for all pairs) in the routine
C                     call.  In this case, setting LDA1 = LDA,
C                     LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD, or
C                     LDA2 = LDA, LDB2 = LDB, LDC2 = LDC, and LDD2 = LDD
C                     will give maximum efficiency.
C
C     Input/Output Parameters
C
C     N1      (input) INTEGER
C             The number of state variables in the first system, i.e.
C             the order of the matrix A1.  N1 >= 0.
C
C     M1      (input) INTEGER
C             The number of input variables for the first system.
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
C     P2      (input) INTEGER
C             The number of output variables from the second system.
C             P2 >= 0.
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
C             The leading P2-by-N2 part of this array must contain the
C             state/output matrix C2 for the second system.
C
C     LDC2    INTEGER
C             The leading dimension of array C2.
C             LDC2 >= MAX(1,P2) if N2 > 0.
C             LDC2 >= 1 if N2 = 0.
C
C     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1)
C             The leading P2-by-P1 part of this array must contain the
C             input/output matrix D2 for the second system.
C
C     LDD2    INTEGER
C             The leading dimension of array D2.  LDD2 >= MAX(1,P2).
C
C     N       (output) INTEGER
C             The number of state variables (N1 + N2) in the resulting
C             system, i.e. the order of the matrix A, the number of rows
C             of B and the number of columns of C.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2)
C             The leading N-by-N part of this array contains the state
C             transition matrix A for the cascaded system.
C             If OVER = 'O', the array A can overlap A1, if UPLO = 'L',
C             or A2, if UPLO = 'U'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N1+N2).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M1)
C             The leading N-by-M1 part of this array contains the
C             input/state matrix B for the cascaded system.
C             If OVER = 'O', the array B can overlap B1, if UPLO = 'L',
C             or B2, if UPLO = 'U'.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N1+N2).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2)
C             The leading P2-by-N part of this array contains the
C             state/output matrix C for the cascaded system.
C             If OVER = 'O', the array C can overlap C1, if UPLO = 'L',
C             or C2, if UPLO = 'U'.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P2) if N1+N2 > 0.
C             LDC >= 1 if N1+N2 = 0.
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M1)
C             The leading P2-by-M1 part of this array contains the
C             input/output matrix D for the cascaded system.
C             If OVER = 'O', the array D can overlap D1, if UPLO = 'L',
C             or D2, if UPLO = 'U'.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P2).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             The array DWORK is not referenced if OVER = 'N'.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, P1*MAX(N1, M1, N2, P2) ) if OVER = 'O'.
C             LDWORK >= 1 if OVER = 'N'.
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
C     After cascaded inter-connection of the two systems
C
C     X1'     = A1*X1 + B1*U
C     V       = C1*X1 + D1*U
C
C     X2'     = A2*X2 + B2*V
C     Y       = C2*X2 + D2*V
C
C     (where  '  denotes differentiation with respect to time)
C
C     the following state-space model will be obtained:
C
C     X'      = A*X + B*U
C     Y       = C*X + D*U
C
C     where matrix  A  has the form   ( A1     0 ),
C                                     ( B2*C1  A2)
C
C           matrix  B  has the form  (  B1   ),
C                                    ( B2*D1 )
C
C           matrix  C  has the form  ( D2*C1  C2 ) and
C
C           matrix  D  has the form  ( D2*D1 ).
C
C     This form is returned by the routine when UPLO = 'L'.  Note that
C     when A1 and A2 are block lower triangular, the resulting state
C     matrix is also block lower triangular.
C
C     By applying a similarity transformation to the system above,
C     using the matrix  ( 0  I ),  where  I  is the identity matrix of
C                       ( J  0 )
C     order  N2,  and  J  is the identity matrix of order  N1,  the
C     system matrices become
C
C           A = ( A2  B2*C1 ),
C               ( 0     A1  )
C
C           B = ( B2*D1 ),
C               (  B1   )
C
C           C = ( C2  D2*C1 ) and
C
C           D = ( D2*D1 ).
C
C     This form is returned by the routine when UPLO = 'U'.  Note that
C     when A1 and A2 are block upper triangular (for instance, in the
C     real Schur form), the resulting state matrix is also block upper
C     triangular.
C
C     REFERENCES
C
C     None
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires P1*(N1+M1)*(N2+P2) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and
C                  A. Varga, German Aerospace Research Establishment,
C                  Oberpfaffenhofen, Germany, Nov. 1996.
C     Supersedes Release 2.0 routine AB05AD by C.J.Benson, Kingston
C     Polytechnic, United Kingdom, January 1982.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
C     Feb. 2004.
C
C     KEYWORDS
C
C     Cascade control, continuous-time system, multivariable
C     system, state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         OVER, UPLO
      INTEGER           INFO, LDA, LDA1, LDA2, LDB, LDB1, LDB2, LDC,
     $                  LDC1, LDC2, LDD, LDD1, LDD2, LDWORK, M1, N, N1,
     $                  N2, P1, P2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), A1(LDA1,*), A2(LDA2,*), B(LDB,*),
     $                  B1(LDB1,*), B2(LDB2,*), C(LDC,*), C1(LDC1,*),
     $                  C2(LDC2,*), D(LDD,*), D1(LDD1,*), D2(LDD2,*),
     $                  DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LOVER, LUPLO
      INTEGER           I, I1, I2, J, LDWN2, LDWP1, LDWP2
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLACPY, DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      LOVER = LSAME( OVER, 'O' )
      LUPLO = LSAME( UPLO, 'L' )
      N = N1 + N2
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( .NOT.LUPLO .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOVER .AND. .NOT.LSAME( OVER, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N1.LT.0 ) THEN
         INFO = -3
      ELSE IF( M1.LT.0 ) THEN
         INFO = -4
      ELSE IF( P1.LT.0 ) THEN
         INFO = -5
      ELSE IF( N2.LT.0 ) THEN
         INFO = -6
      ELSE IF( P2.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDA1.LT.MAX( 1, N1 ) ) THEN
         INFO = -9
      ELSE IF( LDB1.LT.MAX( 1, N1 ) ) THEN
         INFO = -11
      ELSE IF( ( N1.GT.0 .AND. LDC1.LT.MAX( 1, P1 ) ) .OR.
     $         ( N1.EQ.0 .AND. LDC1.LT.1 ) ) THEN
         INFO = -13
      ELSE IF( LDD1.LT.MAX( 1, P1 ) ) THEN
         INFO = -15
      ELSE IF( LDA2.LT.MAX( 1, N2 ) ) THEN
         INFO = -17
      ELSE IF( LDB2.LT.MAX( 1, N2 ) ) THEN
         INFO = -19
      ELSE IF( ( N2.GT.0 .AND. LDC2.LT.MAX( 1, P2 ) ) .OR.
     $         ( N2.EQ.0 .AND. LDC2.LT.1 ) ) THEN
         INFO = -21
      ELSE IF( LDD2.LT.MAX( 1, P2 ) ) THEN
         INFO = -23
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -26
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -28
      ELSE IF( ( N.GT.0 .AND. LDC.LT.MAX( 1, P2 ) ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -30
      ELSE IF( LDD.LT.MAX( 1, P2 ) ) THEN
         INFO = -32
      ELSE IF( ( LOVER.AND.LDWORK.LT.MAX( 1, P1*MAX( N1, M1, N2, P2 )) )
     $.OR.( .NOT.LOVER.AND.LDWORK.LT.1 ) ) THEN
         INFO = -34
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB05MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, MIN( M1, P2 ) ).EQ.0 )
     $   RETURN
C
C     Set row/column indices for storing the results.
C
      IF ( LUPLO ) THEN
         I1 = 1
         I2 = MIN( N1 + 1, N )
      ELSE
         I1 = MIN( N2 + 1, N )
         I2 = 1
      END IF
C
      LDWN2 = MAX( 1, N2 )
      LDWP1 = MAX( 1, P1 )
      LDWP2 = MAX( 1, P2 )
C
C     Construct the cascaded system matrices, taking the desired block
C     structure and possible overwriting into account.
C
C     Form the diagonal blocks of matrix  A.
C
      IF ( LUPLO ) THEN
C
C        Lower block diagonal structure.
C
         IF ( LOVER .AND. LDA1.LE.LDA ) THEN
            IF ( LDA1.LT.LDA ) THEN
C
               DO 20 J = N1, 1, -1
                  DO 10 I = N1, 1, -1
                     A(I,J) = A1(I,J)
   10             CONTINUE
   20          CONTINUE
C
            END IF
         ELSE
            CALL DLACPY( 'F', N1, N1, A1, LDA1, A, LDA )
         END IF
         IF ( N2.GT.0 )
     $      CALL DLACPY( 'F', N2, N2, A2, LDA2, A(I2,I2), LDA )
      ELSE
C
C        Upper block diagonal structure.
C
         IF ( LOVER .AND. LDA2.LE.LDA ) THEN
            IF ( LDA2.LT.LDA ) THEN
C
               DO 40 J = N2, 1, -1
                  DO 30 I = N2, 1, -1
                     A(I,J) = A2(I,J)
   30             CONTINUE
   40          CONTINUE
C
            END IF
         ELSE
            CALL DLACPY( 'F', N2, N2, A2, LDA2, A, LDA )
         END IF
         IF ( N1.GT.0 )
     $      CALL DLACPY( 'F', N1, N1, A1, LDA1, A(I1,I1), LDA )
      END IF
C
C     Form the off-diagonal blocks of matrix  A.
C
      IF ( MIN( N1, N2 ).GT.0 ) THEN
         CALL DLASET( 'F', N1, N2, ZERO, ZERO, A(I1,I2), LDA )
         CALL DGEMM ( 'No transpose', 'No transpose', N2, N1, P1, ONE,
     $                B2, LDB2, C1, LDC1, ZERO, A(I2,I1), LDA )
      END IF
C
      IF ( LUPLO ) THEN
C
C        Form the matrix  B.
C
         IF ( LOVER .AND. LDB1.LE.LDB ) THEN
            IF ( LDB1.LT.LDB ) THEN
C
               DO 60 J = M1, 1, -1
                  DO 50 I = N1, 1, -1
                     B(I,J) = B1(I,J)
   50             CONTINUE
   60          CONTINUE
C
            END IF
         ELSE
            CALL DLACPY( 'F', N1, M1, B1, LDB1, B, LDB )
         END IF
C
         IF ( MIN( N2, M1 ).GT.0 )
     $      CALL DGEMM ( 'No transpose', 'No transpose', N2, M1, P1,
     $                   ONE, B2, LDB2, D1, LDD1, ZERO, B(I2,1), LDB )
C
C        Form the matrix  C.
C
         IF ( N1.GT.0 ) THEN
            IF ( LOVER ) THEN
C
C              Workspace:  P1*N1.
C
               CALL DLACPY( 'F', P1, N1, C1, LDC1, DWORK, LDWP1 )
               CALL DGEMM ( 'No transpose', 'No transpose', P2, N1, P1,
     $                      ONE, D2, LDD2, DWORK, LDWP1, ZERO, C, LDC )
            ELSE
               CALL DGEMM ( 'No transpose', 'No transpose', P2, N1, P1,
     $                      ONE, D2, LDD2, C1, LDC1, ZERO, C, LDC )
            END IF
         END IF
C
         IF ( MIN( P2, N2 ).GT.0 )
     $      CALL DLACPY( 'F', P2, N2, C2, LDC2, C(1,I2), LDC )
C
C        Now form the matrix  D.
C
         IF ( LOVER ) THEN
C
C           Workspace:  P1*M1.
C
            CALL DLACPY( 'F', P1, M1, D1, LDD1, DWORK, LDWP1 )
            CALL DGEMM ( 'No transpose', 'No transpose', P2, M1, P1,
     $                   ONE, D2, LDD2, DWORK, LDWP1, ZERO, D, LDD )
         ELSE
            CALL DGEMM ( 'No transpose', 'No transpose', P2, M1, P1,
     $                    ONE, D2, LDD2, D1, LDD1, ZERO, D, LDD )
         END IF
C
      ELSE
C
C        Form the matrix  B.
C
         IF ( LOVER ) THEN
C
C           Workspace:  N2*P1.
C
            CALL DLACPY( 'F', N2, P1, B2, LDB2, DWORK, LDWN2 )
            IF ( MIN( N2, M1 ).GT.0 )
     $         CALL DGEMM ( 'No transpose', 'No transpose', N2, M1, P1,
     $                      ONE, DWORK, LDWN2, D1, LDD1, ZERO, B(I2,1),
     $                      LDB )
         ELSE
            CALL DGEMM ( 'No transpose', 'No transpose', N2, M1, P1,
     $                   ONE, B2, LDB2, D1, LDD1, ZERO, B, LDB )
         END IF
C
         IF ( MIN( N1, M1 ).GT.0 )
     $      CALL DLACPY( 'F', N1, M1, B1, LDB1, B(I1,1), LDB )
C
C        Form the matrix  C.
C
         IF ( LOVER .AND. LDC2.LE.LDC ) THEN
            IF ( LDC2.LT.LDC ) THEN
C
               DO 80 J = N2, 1, -1
                  DO 70 I = P2, 1, -1
                     C(I,J) = C2(I,J)
   70             CONTINUE
   80          CONTINUE
C
            END IF
         ELSE
            CALL DLACPY( 'F', P2, N2, C2, LDC2, C, LDC )
         END IF
C
         IF ( MIN( P2, N1 ).GT.0 )
     $      CALL DGEMM ( 'No transpose', 'No transpose', P2, N1, P1,
     $                   ONE, D2, LDD2, C1, LDC1, ZERO, C(1,I1), LDC )
C
C        Now form the matrix  D.
C
         IF ( LOVER ) THEN
C
C           Workspace:  P2*P1.
C
            CALL DLACPY( 'F', P2, P1, D2, LDD2, DWORK, LDWP2 )
            CALL DGEMM ( 'No transpose', 'No transpose', P2, M1, P1,
     $                   ONE, DWORK, LDWP2, D1, LDD1, ZERO, D, LDD )
         ELSE
            CALL DGEMM ( 'No transpose', 'No transpose', P2, M1, P1,
     $                    ONE, D2, LDD2, D1, LDD1, ZERO, D, LDD )
         END IF
      END IF
C
      RETURN
C *** Last line of AB05MD ***
      END
