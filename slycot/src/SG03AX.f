      SUBROUTINE SG03AX( TRANS, N, A, LDA, E, LDE, X, LDX, SCALE, INFO )
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
C     To solve for X either the reduced generalized discrete-time
C     Lyapunov equation
C
C         T            T
C        A  * X * A - E  * X * E  =  SCALE * Y                       (1)
C
C     or
C
C                 T            T
C        A * X * A  - E * X * E   =  SCALE * Y                       (2)
C
C     where the right hand side Y is symmetric. A, E, Y, and the
C     solution X are N-by-N matrices. The pencil A - lambda * E must be
C     in generalized Schur form (A upper quasitriangular, E upper
C     triangular). SCALE is an output scale factor, set to avoid
C     overflow in X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANS   CHARACTER*1
C             Specifies whether the transposed equation is to be solved
C             or not:
C             = 'N':  Solve equation (1);
C             = 'T':  Solve equation (2).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N upper Hessenberg part of this array
C             must contain the quasitriangular matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             The leading N-by-N upper triangular part of this array
C             must contain the matrix E.
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, the leading N-by-N part of this array must
C             contain the right hand side matrix Y of the equation. Only
C             the upper triangular part of this matrix need be given.
C             On exit, the leading N-by-N part of this array contains
C             the solution matrix X of the equation.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in X.
C             (0 < SCALE <= 1)
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  equation is (almost) singular to working precision;
C                   perturbed values were used to solve the equation
C                   (but the matrices A and E are unchanged).
C
C     METHOD
C
C     The solution X of (1) or (2) is computed via block back
C     substitution or block forward substitution, respectively. (See
C     [1] and [2] for details.)
C
C     REFERENCES
C
C     [1] Bartels, R.H., Stewart, G.W.
C         Solution of the equation A X + X B = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Penzl, T.
C         Numerical solution of generalized Lyapunov equations.
C         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
C
C     NUMERICAL ASPECTS
C
C     8/3 * N**3 flops are required by the routine. Note that we count a
C     single floating point arithmetic operation as one flop.
C
C     The algorithm is backward stable if the eigenvalues of the pencil
C     A - lambda * E are real. Otherwise, linear systems of order at
C     most 4 are involved into the computation. These systems are solved
C     by Gauss elimination with complete pivoting. The loss of stability
C     of the Gauss elimination with complete pivoting is rarely
C     encountered in practice.
C
C     CONTRIBUTOR
C
C     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998.
C
C     REVISIONS
C
C     Sep. 1998 (V. Sima).
C     Dec. 1998 (V. Sima).
C
C     KEYWORDS
C
C     Lyapunov equation
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  MONE, ONE, ZERO
      PARAMETER         ( MONE = -1.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANS
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, LDA, LDE, LDX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), X(LDX,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AK11, AK12, AK21, AK22, AL11, AL12, AL21, AL22,
     $                  EK11, EK12, EK22, EL11, EL12, EL22, SCALE1
      INTEGER           DIMMAT, I, INFO1, KB, KH, KL, LB, LH, LL
      LOGICAL           NOTRNS
C     .. Local Arrays ..
      DOUBLE PRECISION  MAT(4,4), RHS(4), TM(2,2)
      INTEGER           PIV1(4), PIV2(4)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DGEMV, DSCAL, MB02UU,
     $                  MB02UV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Decode input parameter.
C
      NOTRNS = LSAME( TRANS, 'N' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( NOTRNS .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSEIF ( N .LT. 0 ) THEN
         INFO = -2
      ELSEIF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -4
      ELSEIF ( LDE .LT. MAX( 1, N ) ) THEN
         INFO = -6
      ELSEIF ( LDX .LT. MAX( 1, N ) ) THEN
         INFO = -8
      ELSE
         INFO = 0
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SG03AX', -INFO )
         RETURN
      END IF
C
      SCALE = ONE
C
C     Quick return if possible.
C
      IF ( N .EQ. 0 ) RETURN
C
      IF ( NOTRNS ) THEN
C
C        Solve equation (1).
C
C        Outer Loop. Compute block row X(KL:KH,:). KB denotes the number
C        of rows in this block row.
C
         KL = 0
         KB = 1
C        WHILE ( KL+KB .LE. N ) DO
   20    IF ( KL+KB .LE. N ) THEN
            KL = KL + KB
            IF ( KL .EQ. N ) THEN
               KB = 1
            ELSE
               IF ( A(KL+1,KL) .NE. ZERO ) THEN
                  KB = 2
               ELSE
                  KB = 1
               END IF
            END IF
            KH = KL + KB - 1
C
C           Copy elements of solution already known by symmetry.
C
C              X(KL:KH,1:KL-1) = X(1:KL-1,KL:KH)'
C
            IF ( KL .GT. 1 ) THEN
               DO 40 I = KL, KH
                  CALL DCOPY( KL-1, X(1,I), 1, X(I,1), LDX )
   40          CONTINUE
            END IF
C
C           Inner Loop. Compute block X(KL:KH,LL:LH). LB denotes the
C           number of columns in this block.
C
            LL = KL - 1
            LB = 1
C           WHILE ( LL+LB .LE. N ) DO
   60       IF ( LL+LB .LE. N ) THEN
               LL = LL + LB
               IF ( LL .EQ. N ) THEN
                  LB = 1
               ELSE
                  IF ( A(LL+1,LL) .NE. ZERO ) THEN
                     LB = 2
                  ELSE
                     LB = 1
                  END IF
               END IF
               LH = LL + LB - 1
C
C              Update right hand sides (I).
C
C                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) -
C                    A(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*A(1:LL-1,LL:LH))
C
C                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) +
C                    E(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*E(1:LL-1,LL:LH))
C
               IF ( LL .GT. 1 ) THEN
                  CALL DGEMM( 'N', 'N', KB, LB, LL-1, ONE, X(KL,1), LDX,
     $                        A(1,LL), LDA, ZERO, TM, 2 )
                  CALL DGEMM( 'T', 'N', LH-KL+1, LB, KB, MONE, A(KL,KL),
     $                        LDA, TM, 2, ONE, X(KL,LL), LDX )
                  CALL DGEMM( 'N', 'N', KB, LB, LL-1, ONE, X(KL,1),
     $                        LDX, E(1,LL), LDE, ZERO, TM, 2 )
                  CALL DGEMM( 'T', 'N', LH-KH+1, LB, KB, ONE, E(KL,KH),
     $                        LDE, TM, 2, ONE, X(KH,LL), LDX )
                  IF ( KB .EQ. 2 ) CALL DAXPY( LB, E(KL,KL), TM, 2,
     $                                         X(KL,LL), LDX )
               END IF
C
C              Solve small Sylvester equations of order at most (2,2).
C
               IF ( KB.EQ.1 .AND. LB.EQ.1 ) THEN
C
                  DIMMAT = 1
C
                  MAT(1,1) = A(LL,LL)*A(KL,KL) - E(LL,LL)*E(KL,KL)
C
                  RHS(1) = X(KL,LL)
C
               ELSEIF ( KB.EQ.2 .AND. LB.EQ.1 ) THEN
C
                  DIMMAT = 2
C
                  AK11 = A(KL,KL)
                  AK12 = A(KL,KH)
                  AK21 = A(KH,KL)
                  AK22 = A(KH,KH)
C
                  AL11 = A(LL,LL)
C
                  EK11 = E(KL,KL)
                  EK12 = E(KL,KH)
                  EK22 = E(KH,KH)
C
                  EL11 = E(LL,LL)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL11*AK21
                  MAT(2,1) = AL11*AK12 - EL11*EK12
                  MAT(2,2) = AL11*AK22 - EL11*EK22
C
                  RHS(1) = X(KL,LL)
                  RHS(2) = X(KH,LL)
C
               ELSEIF ( KB.EQ.1 .AND. LB.EQ.2 ) THEN
C
                  DIMMAT = 2
C
                  AK11 = A(KL,KL)
C
                  AL11 = A(LL,LL)
                  AL12 = A(LL,LH)
                  AL21 = A(LH,LL)
                  AL22 = A(LH,LH)
C
                  EK11 = E(KL,KL)
C
                  EL11 = E(LL,LL)
                  EL12 = E(LL,LH)
                  EL22 = E(LH,LH)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL21*AK11
                  MAT(2,1) = AL12*AK11 - EL12*EK11
                  MAT(2,2) = AL22*AK11 - EL22*EK11
C
                  RHS(1) = X(KL,LL)
                  RHS(2) = X(KL,LH)
C
               ELSE
C
                  DIMMAT = 4
C
                  AK11 = A(KL,KL)
                  AK12 = A(KL,KH)
                  AK21 = A(KH,KL)
                  AK22 = A(KH,KH)
C
                  AL11 = A(LL,LL)
                  AL12 = A(LL,LH)
                  AL21 = A(LH,LL)
                  AL22 = A(LH,LH)
C
                  EK11 = E(KL,KL)
                  EK12 = E(KL,KH)
                  EK22 = E(KH,KH)
C
                  EL11 = E(LL,LL)
                  EL12 = E(LL,LH)
                  EL22 = E(LH,LH)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL11*AK21
                  MAT(1,3) = AL21*AK11
                  MAT(1,4) = AL21*AK21
C
                  MAT(2,1) = AL11*AK12 - EL11*EK12
                  MAT(2,2) = AL11*AK22 - EL11*EK22
                  MAT(2,3) = AL21*AK12
                  MAT(2,4) = AL21*AK22
C
                  MAT(3,1) = AL12*AK11 - EL12*EK11
                  MAT(3,2) = AL12*AK21
                  MAT(3,3) = AL22*AK11 - EL22*EK11
                  MAT(3,4) = AL22*AK21
C
                  MAT(4,1) = AL12*AK12 - EL12*EK12
                  MAT(4,2) = AL12*AK22 - EL12*EK22
                  MAT(4,3) = AL22*AK12 - EL22*EK12
                  MAT(4,4) = AL22*AK22 - EL22*EK22
C
                  RHS(1) = X(KL,LL)
                  IF ( KL .EQ. LL ) THEN
                     RHS(2) = X(KL,KH)
                  ELSE
                     RHS(2) = X(KH,LL)
                  END IF
                  RHS(3) = X(KL,LH)
                  RHS(4) = X(KH,LH)
C
               END IF
C
               CALL MB02UV( DIMMAT, MAT, 4, PIV1, PIV2, INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 1
               CALL MB02UU( DIMMAT, MAT, 4, RHS, PIV1, PIV2, SCALE1 )
C
C              Scaling.
C
               IF ( SCALE1 .NE. ONE ) THEN
                  DO 80 I = 1, N
                     CALL DSCAL( N, SCALE1, X(1,I), 1 )
   80             CONTINUE
                  SCALE = SCALE*SCALE1
               END IF
C
               IF ( LB.EQ.1 .AND. KB.EQ.1 ) THEN
                  X(KL,LL) = RHS(1)
               ELSEIF ( LB.EQ.1 .AND. KB.EQ.2 ) THEN
                  X(KL,LL) = RHS(1)
                  X(KH,LL) = RHS(2)
               ELSEIF ( LB.EQ.2 .AND. KB.EQ.1 ) THEN
                  X(KL,LL) = RHS(1)
                  X(KL,LH) = RHS(2)
               ELSE
                  X(KL,LL) = RHS(1)
                  X(KH,LL) = RHS(2)
                  X(KL,LH) = RHS(3)
                  X(KH,LH) = RHS(4)
               END IF
C
C              Update right hand sides (II).
C
C              X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) -
C                 A(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*A(LL:LH,LL:LH))
C
C              X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) +
C                 E(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*E(LL:LH,LL:LH))
C
               IF ( KL .LT. LL ) THEN
                  CALL DGEMM( 'N', 'N', KB, LB, LB, ONE, X(KL,LL), LDX,
     $                        A(LL,LL), LDA, ZERO, TM, 2 )
                  CALL DGEMM( 'T', 'N', LH-KH, LB, KB, MONE, A(KL,KH+1),
     $                        LDA, TM, 2, ONE, X(KH+1,LL), LDX )
                  IF ( LB .EQ. 2 ) THEN
                     CALL DCOPY( KB, X(KL,LL), 1, TM, 1 )
                     CALL DSCAL( KB, E(LL,LL), TM, 1 )
                  END IF
                  CALL DGEMV( 'N', KB, LB, ONE, X(KL,LL), LDX, E(LL,LH),
     $                        1, ZERO, TM(1,LB), 1 )
                  CALL DGEMM( 'T', 'N', LH-KH, LB, KB, ONE, E(KL,KH+1),
     $                        LDE, TM, 2, ONE, X(KH+1,LL), LDX )
               END IF
C
            GOTO 60
            END IF
C           END WHILE 60
C
         GOTO 20
         END IF
C        END WHILE 20
C
      ELSE
C
C        Solve equation (2).
C
C        Outer Loop. Compute block column X(:,LL:LH). LB denotes the
C        number of columns in this block column.
C
         LL = N + 1
C        WHILE ( LL .GT. 1 ) DO
  100    IF ( LL .GT. 1 ) THEN
            LH = LL - 1
            IF ( LH .EQ. 1 ) THEN
               LB = 1
            ELSE
               IF ( A(LL-1,LL-2) .NE. ZERO ) THEN
                  LB = 2
               ELSE
                  LB = 1
               END IF
            END IF
            LL = LL - LB
C
C           Copy elements of solution already known by symmetry.
C
C              X(LH+1:N,LL:LH) = X(LL:LH,LH+1:N)'
C
            IF ( LH .LT. N ) THEN
               DO 120 I = LL, LH
                  CALL DCOPY( N-LH, X(I,LH+1), LDX, X(LH+1,I), 1 )
  120          CONTINUE
            END IF
C
C           Inner Loop. Compute block X(KL:KH,LL:LH). KB denotes the
C           number of rows in this block.
C
            KL = LH + 1
C           WHILE ( KL .GT. 1 ) DO
  140       IF ( KL .GT. 1 ) THEN
               KH = KL - 1
               IF ( KH .EQ. 1 ) THEN
                  KB = 1
               ELSE
                  IF ( A(KL-1,KL-2) .NE. ZERO ) THEN
                     KB =2
                  ELSE
                     KB = 1
                  END IF
               END IF
               KL = KL - KB
C
C              Update right hand sides (I).
C
C                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) -
C                    (A(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*A(KL:LH,LL:LH)'
C
C                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) +
C                    (E(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*E(KL:LH,LL:LH)'
C
               IF ( KH .LT. N ) THEN
                  CALL DGEMM( 'N', 'N', KB, LB, N-KH, ONE, A(KL,KH+1),
     $                        LDA, X(KH+1,LL), LDX, ZERO, TM, 2 )
                  CALL DGEMM( 'N', 'T', KB, LH-KL+1, LB, MONE, TM, 2,
     $                        A(KL,LL), LDA, ONE, X(KL,KL), LDX )
                  CALL DGEMM( 'N', 'N', KB, LB, N-KH, ONE, E(KL,KH+1),
     $                        LDE, X(KH+1,LL), LDX, ZERO, TM, 2 )
                  CALL DGEMM( 'N', 'T', KB, LL-KL+1, LB, ONE, TM, 2,
     $                        E(KL,LL), LDE, ONE, X(KL,KL), LDX )
                  IF ( LB .EQ. 2 ) CALL DAXPY( KB, E(LH,LH), TM(1,2), 1,
     $                                         X(KL,LH), 1 )
               END IF
C
C              Solve small Sylvester equations of order at most (2,2).
C
               IF ( KB.EQ.1 .AND. LB.EQ.1 ) THEN
C
                  DIMMAT = 1
C
                  MAT(1,1) = A(LL,LL)*A(KL,KL) - E(LL,LL)*E(KL,KL)
C
                  RHS(1) = X(KL,LL)
C
               ELSEIF ( KB.EQ.2 .AND. LB.EQ.1 ) THEN
C
                  DIMMAT = 2
C
                  AK11 = A(KL,KL)
                  AK12 = A(KL,KH)
                  AK21 = A(KH,KL)
                  AK22 = A(KH,KH)
C
                  AL11 = A(LL,LL)
C
                  EK11 = E(KL,KL)
                  EK12 = E(KL,KH)
                  EK22 = E(KH,KH)
C
                  EL11 = E(LL,LL)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL11*AK12 - EL11*EK12
                  MAT(2,1) = AL11*AK21
                  MAT(2,2) = AL11*AK22 - EL11*EK22
C
                  RHS(1) = X(KL,LL)
                  RHS(2) = X(KH,LL)
C
               ELSEIF ( KB.EQ.1 .AND. LB.EQ.2 ) THEN
C
                  DIMMAT = 2
C
                  AK11 = A(KL,KL)
C
                  AL11 = A(LL,LL)
                  AL12 = A(LL,LH)
                  AL21 = A(LH,LL)
                  AL22 = A(LH,LH)
C
                  EK11 = E(KL,KL)
C
                  EL11 = E(LL,LL)
                  EL12 = E(LL,LH)
                  EL22 = E(LH,LH)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL12*AK11 - EL12*EK11
                  MAT(2,1) = AL21*AK11
                  MAT(2,2) = AL22*AK11 - EL22*EK11
C
                  RHS(1) = X(KL,LL)
                  RHS(2) = X(KL,LH)
C
               ELSE
C
                  DIMMAT = 4
C
                  AK11 = A(KL,KL)
                  AK12 = A(KL,KH)
                  AK21 = A(KH,KL)
                  AK22 = A(KH,KH)
C
                  AL11 = A(LL,LL)
                  AL12 = A(LL,LH)
                  AL21 = A(LH,LL)
                  AL22 = A(LH,LH)
C
                  EK11 = E(KL,KL)
                  EK12 = E(KL,KH)
                  EK22 = E(KH,KH)
C
                  EL11 = E(LL,LL)
                  EL12 = E(LL,LH)
                  EL22 = E(LH,LH)
C
                  MAT(1,1) = AL11*AK11 - EL11*EK11
                  MAT(1,2) = AL11*AK12 - EL11*EK12
                  MAT(1,3) = AL12*AK11 - EL12*EK11
                  MAT(1,4) = AL12*AK12 - EL12*EK12
C
                  MAT(2,1) = AL11*AK21
                  MAT(2,2) = AL11*AK22 - EL11*EK22
                  MAT(2,3) = AL12*AK21
                  MAT(2,4) = AL12*AK22 - EL12*EK22
C
                  MAT(3,1) = AL21*AK11
                  MAT(3,2) = AL21*AK12
                  MAT(3,3) = AL22*AK11 - EL22*EK11
                  MAT(3,4) = AL22*AK12 - EL22*EK12
C
                  MAT(4,1) = AL21*AK21
                  MAT(4,2) = AL21*AK22
                  MAT(4,3) = AL22*AK21
                  MAT(4,4) = AL22*AK22 - EL22*EK22
C
                  RHS(1) = X(KL,LL)
                  IF ( KL .EQ. LL ) THEN
                     RHS(2) = X(KL,KH)
                  ELSE
                     RHS(2) = X(KH,LL)
                  END IF
                  RHS(3) = X(KL,LH)
                  RHS(4) = X(KH,LH)
C
               END IF
C
               CALL MB02UV( DIMMAT, MAT, 4, PIV1, PIV2, INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 1
               CALL MB02UU( DIMMAT, MAT, 4, RHS, PIV1, PIV2, SCALE1 )
C
C              Scaling.
C
               IF ( SCALE1 .NE. ONE ) THEN
                  DO 160 I = 1, N
                     CALL DSCAL( N, SCALE1, X(1,I), 1 )
  160             CONTINUE
                  SCALE = SCALE*SCALE1
               END IF
C
               IF ( LB.EQ.1 .AND. KB.EQ.1 ) THEN
                  X(KL,LL) = RHS(1)
               ELSEIF ( LB.EQ.1 .AND. KB.EQ.2 ) THEN
                  X(KL,LL) = RHS(1)
                  X(KH,LL) = RHS(2)
               ELSEIF ( LB.EQ.2 .AND. KB.EQ.1 ) THEN
                  X(KL,LL) = RHS(1)
                  X(KL,LH) = RHS(2)
               ELSE
                  X(KL,LL) = RHS(1)
                  X(KH,LL) = RHS(2)
                  X(KL,LH) = RHS(3)
                  X(KH,LH) = RHS(4)
               END IF
C
C              Update right hand sides (II).
C
C                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) -
C                    (A(KL:KH,KL:KH)*X(KL:KH,LL:LH))*A(KL:LL-1,LL:LH)'
C
C                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) +
C                    (E(KL:KH,KL:KH)*X(KL:KH,LL:LH))*E(KL:LL-1,LL:LH)'
C
               IF ( KL .LT. LL ) THEN
                  CALL DGEMM( 'N', 'N', KB, LB, KB, ONE, A(KL,KL), LDA,
     $                        X(KL,LL), LDX, ZERO, TM, 2 )
                  CALL DGEMM( 'N', 'T', KB, LL-KL, LB, MONE, TM, 2,
     $                        A(KL,LL), LDA, ONE, X(KL,KL), LDX )
                  CALL DGEMV( 'T', KB, LB, ONE, X(KL,LL), LDX, E(KL,KL),
     $                        LDE, ZERO, TM, 2 )
                  IF ( KB .EQ. 2 ) THEN
                     CALL DCOPY( LB, X(KH,LL), LDX, TM(2,1), 2 )
                     CALL DSCAL( LB, E(KH,KH), TM(2,1), 2 )
                  END IF
                  CALL DGEMM( 'N', 'T', KB, LL-KL, LB, ONE, TM, 2,
     $                        E(KL,LL), LDE, ONE, X(KL,KL), LDX )
               END IF
C
            GOTO 140
            END IF
C           END WHILE 140
C
         GOTO 100
         END IF
C        END WHILE 100
C
      END IF
C
      RETURN
C *** Last line of SG03AX ***
      END
