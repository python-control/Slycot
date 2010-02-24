      SUBROUTINE SG03BW( TRANS, M, N, A, LDA, C, LDC, E, LDE, D, LDD, X,
     $                   LDX, SCALE, INFO )
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
C     To solve for X the generalized Sylvester equation
C
C         T            T
C        A  * X * C + E  * X * D  =  SCALE * Y,                      (1)
C
C     or the transposed equation
C
C                 T            T
C        A * X * C  + E * X * D   =  SCALE * Y,                      (2)
C
C     where A and E are real M-by-M matrices, C and D are real N-by-N
C     matrices, X and Y are real M-by-N matrices. N is either 1 or 2.
C     The pencil A - lambda * E must be in generalized real Schur form
C     (A upper quasitriangular, E upper triangular). SCALE is an output
C     scale factor, set to avoid overflow in X.
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
C     M       (input) INTEGER
C             The order of the matrices A and E.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrices C and D.  N = 1 or N = 2.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,M)
C             The leading M-by-M part of this array must contain the
C             upper quasitriangular matrix A. The elements below the
C             upper Hessenberg part are not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,M).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading N-by-N part of this array must contain the
C             matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,M)
C             The leading M-by-M part of this array must contain the
C             upper triangular matrix E. The elements below the main
C             diagonal are not referenced.
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= MAX(1,M).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,N)
C             The leading N-by-N part of this array must contain the
C             matrix D.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right hand side matrix Y.
C             On exit, the leading M-by-N part of this array contains
C             the solution matrix X.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= MAX(1,M).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in X.
C             0 < SCALE <= 1.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the generalized Sylvester equation is (nearly)
C                   singular to working precision;  perturbed values
C                   were used to solve the equation (but the matrices
C                   A, C, D, and E are unchanged).
C
C     METHOD
C
C     The method used by the routine is based on a generalization of the
C     algorithm due to Bartels and Stewart [1]. See also [2] and [3] for
C     details.
C
C     REFERENCES
C
C     [1] Bartels, R.H., Stewart, G.W.
C         Solution of the equation A X + X B = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Gardiner, J.D., Laub, A.J., Amato, J.J., Moler, C.B.
C         Solution of the Sylvester Matrix Equation
C         A X B**T + C X D**T = E.
C         A.C.M. Trans. Math. Soft., vol. 18, no. 2, pp. 223-231, 1992.
C
C     [3] Penzl, T.
C         Numerical solution of generalized Lyapunov equations.
C         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
C
C     NUMERICAL ASPECTS
C
C     The routine requires about 2 * N * M**2 flops. Note that we count
C     a single floating point arithmetic operation as one flop.
C
C     The algorithm is backward stable if the eigenvalues of the pencil
C     A - lambda * E are real. Otherwise, linear systems of order at
C     most 4 are involved into the computation. These systems are solved
C     by Gauss elimination with complete pivoting. The loss of stability
C     of the Gauss elimination with complete pivoting is rarely
C     encountered in practice.
C
C     FURTHER COMMENTS
C
C     When near singularity is detected, perturbed values are used
C     to solve the equation (but the given matrices are unchanged).
C
C     CONTRIBUTOR
C
C     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998.
C
C     REVISIONS
C
C     Sep. 1998 (V. Sima).
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
      INTEGER           INFO, LDA, LDC, LDD, LDE, LDX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*), D(LDD,*), E(LDE,*), X(LDX,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  SCALE1
      INTEGER           DIMMAT, I, INFO1, J, MA, MAI, MAJ, MB, ME
      LOGICAL           NOTRNS
C     .. Local Arrays ..
      DOUBLE PRECISION  MAT(4,4), RHS(4), TM(2,2)
      INTEGER           PIV1(4), PIV2(4)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DSCAL, MB02UU, MB02UV, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     Decode input parameters.
C
      NOTRNS = LSAME( TRANS, 'N' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( NOTRNS .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSEIF ( M .LT. 0 ) THEN
         INFO = -2
      ELSEIF ( N .NE. 1 .AND. N .NE. 2 ) THEN
         INFO = -3
      ELSEIF ( LDA .LT. MAX( 1, M ) ) THEN
         INFO = -5
      ELSEIF ( LDC .LT. MAX( 1, N ) ) THEN
         INFO = -7
      ELSEIF ( LDE .LT. MAX( 1, M ) ) THEN
         INFO = -9
      ELSEIF ( LDD .LT. MAX( 1, N ) ) THEN
         INFO = -11
      ELSEIF ( LDX .LT. MAX( 1, M ) ) THEN
         INFO = -13
      ELSE
         INFO = 0
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SG03BW', -INFO )
         RETURN
      END IF
C
      SCALE = ONE
C
C     Quick return if possible.
C
      IF ( M .EQ. 0 )
     $    RETURN
C
      IF ( NOTRNS ) THEN
C
C        Solve equation (1).
C
C        Compute block row X(MA:ME,:). MB denotes the number of rows in
C        this block row.
C
         ME = 0
C        WHILE ( ME .NE. M ) DO
   20    IF ( ME .NE. M ) THEN
            MA = ME + 1
            IF ( MA .EQ. M ) THEN
               ME = M
               MB = 1
            ELSE
               IF ( A(MA+1,MA) .EQ. ZERO ) THEN
                  ME = MA
                  MB = 1
               ELSE
                  ME = MA + 1
                  MB = 2
               END IF
            END IF
C
C           Assemble Kronecker product system of linear equations with
C           matrix
C
C              MAT = kron(C',A(MA:ME,MA:ME)') + kron(D',E(MA:ME,MA:ME)')
C
C           and right hand side
C
C              RHS = vec(X(MA:ME,:)).
C
            IF ( N .EQ. 1 ) THEN
               DIMMAT = MB
               DO 60 I = 1, MB
                  MAI = MA + I - 1
                  DO 40 J = 1, MB
                     MAJ = MA + J - 1
                     MAT(I,J) = C(1,1)*A(MAJ,MAI)
                     IF ( MAJ .LE. MAI )
     $                  MAT(I,J) = MAT(I,J) + D(1,1)*E(MAJ,MAI)
   40             CONTINUE
                  RHS(I) = X(MAI,1)
   60          CONTINUE
            ELSE
               DIMMAT = 2*MB
               DO 100 I = 1, MB
                  MAI = MA + I - 1
                  DO 80 J = 1, MB
                     MAJ = MA + J - 1
                     MAT(I,J) = C(1,1)*A(MAJ,MAI)
                     MAT(MB+I,J) = C(1,2)*A(MAJ,MAI)
                     MAT(I,MB+J) = C(2,1)*A(MAJ,MAI)
                     MAT(MB+I,MB+J) = C(2,2)*A(MAJ,MAI)
                     IF ( MAJ .LE. MAI ) THEN
                        MAT(I,J) = MAT(I,J) + D(1,1)*E(MAJ,MAI)
                        MAT(MB+I,J) = MAT(MB+I,J) + D(1,2)*E(MAJ,MAI)
                        MAT(I,MB+J) = MAT(I,MB+J) + D(2,1)*E(MAJ,MAI)
                        MAT(MB+I,MB+J) = MAT(MB+I,MB+J) +
     $                                   D(2,2)*E(MAJ,MAI)
                     END IF
   80             CONTINUE
                  RHS(I) = X(MAI,1)
                  RHS(MB+I) = X(MAI,2)
  100          CONTINUE
            END IF
C
C           Solve the system of linear equations.
C
            CALL MB02UV( DIMMAT, MAT, 4, PIV1, PIV2, INFO1 )
            IF ( INFO1 .NE. 0 )
     $         INFO = 1
            CALL MB02UU( DIMMAT, MAT, 4, RHS, PIV1, PIV2, SCALE1 )
            IF ( SCALE1 .NE. ONE ) THEN
               SCALE = SCALE1*SCALE
               DO 120 I = 1, N
                  CALL DSCAL( M, SCALE1, X(1,I), 1 )
  120          CONTINUE
            END IF
C
            IF ( N .EQ. 1 ) THEN
               DO 140 I = 1, MB
                  MAI = MA + I - 1
                  X(MAI,1) = RHS(I)
  140          CONTINUE
            ELSE
               DO 160 I = 1, MB
                  MAI = MA + I - 1
                  X(MAI,1) = RHS(I)
                  X(MAI,2) = RHS(MB+I)
  160          CONTINUE
            END IF
C
C           Update right hand sides.
C
C           X(ME+1:M,:) = X(ME+1:M,:) - A(MA:ME,ME+1:M)'*X(MA:ME,:)*C
C
C           X(ME+1:M,:) = X(ME+1:M,:) - E(MA:ME,ME+1:M)'*X(MA:ME,:)*D
C
            IF ( ME .LT. M ) THEN
               CALL DGEMM( 'N', 'N', MB, N, N, ONE, X(MA,1), LDX, C,
     $                     LDC, ZERO, TM, 2 )
               CALL DGEMM( 'T', 'N', M-ME, N, MB, MONE, A(MA,ME+1),
     $                     LDA, TM, 2, ONE, X(ME+1,1), LDX )
               CALL DGEMM( 'N', 'N', MB, N, N, ONE, X(MA,1), LDX, D,
     $                     LDD, ZERO, TM, 2 )
               CALL DGEMM( 'T', 'N', M-ME, N, MB, MONE, E(MA,ME+1), LDE,
     $                     TM, 2, ONE, X(ME+1,1), LDX )
            END IF
C
         GOTO 20
         END IF
C        END WHILE 20
C
      ELSE
C
C        Solve equation (2).
C
C        Compute block row X(MA:ME,:). MB denotes the number of rows in
C        this block row.
C
         MA = M + 1
C        WHILE ( MA .NE. 1 ) DO
  180    IF ( MA .NE. 1 ) THEN
            ME = MA - 1
            IF ( ME .EQ. 1 ) THEN
               MA = 1
               MB = 1
            ELSE
               IF ( A(ME,ME-1) .EQ. ZERO ) THEN
                  MA = ME
                  MB = 1
               ELSE
                  MA = ME - 1
                  MB = 2
               END IF
            END IF
C
C           Assemble Kronecker product system of linear equations with
C           matrix
C
C              MAT = kron(C,A(MA:ME,MA:ME)) + kron(D,E(MA:ME,MA:ME))
C
C           and right hand side
C
C              RHS = vec(X(MA:ME,:)).
C
            IF ( N .EQ. 1 ) THEN
               DIMMAT = MB
               DO 220 I = 1, MB
                  MAI = MA + I - 1
                  DO 200 J = 1, MB
                     MAJ = MA + J - 1
                     MAT(I,J) = C(1,1)*A(MAI,MAJ)
                     IF ( MAJ .GE. MAI )
     $                  MAT(I,J) = MAT(I,J) + D(1,1)*E(MAI,MAJ)
  200             CONTINUE
                  RHS(I) = X(MAI,1)
  220          CONTINUE
            ELSE
               DIMMAT = 2*MB
               DO 260 I = 1, MB
                  MAI = MA + I - 1
                  DO 240 J = 1, MB
                     MAJ = MA + J - 1
                     MAT(I,J) = C(1,1)*A(MAI,MAJ)
                     MAT(MB+I,J) = C(2,1)*A(MAI,MAJ)
                     MAT(I,MB+J) = C(1,2)*A(MAI,MAJ)
                     MAT(MB+I,MB+J) = C(2,2)*A(MAI,MAJ)
                     IF ( MAJ .GE. MAI ) THEN
                        MAT(I,J) = MAT(I,J) + D(1,1)*E(MAI,MAJ)
                        MAT(MB+I,J) = MAT(MB+I,J) + D(2,1)*E(MAI,MAJ)
                        MAT(I,MB+J) = MAT(I,MB+J) + D(1,2)*E(MAI,MAJ)
                        MAT(MB+I,MB+J) = MAT(MB+I,MB+J) +
     $                                   D(2,2)*E(MAI,MAJ)
                     END IF
  240             CONTINUE
                  RHS(I) = X(MAI,1)
                  RHS(MB+I) = X(MAI,2)
  260          CONTINUE
            END IF
C
C           Solve the system of linear equations.
C
            CALL MB02UV( DIMMAT, MAT, 4, PIV1, PIV2, INFO1 )
            IF ( INFO1 .NE. 0 )
     $         INFO = 1
            CALL MB02UU( DIMMAT, MAT, 4, RHS, PIV1, PIV2, SCALE1 )
            IF ( SCALE1 .NE. ONE ) THEN
               SCALE = SCALE1*SCALE
               DO 280 I = 1, N
                  CALL DSCAL( M, SCALE1, X(1,I), 1 )
  280          CONTINUE
            END IF
C
            IF ( N .EQ. 1 ) THEN
               DO 300 I = 1, MB
                  MAI = MA + I - 1
                  X(MAI,1) = RHS(I)
  300          CONTINUE
            ELSE
               DO 320 I = 1, MB
                  MAI = MA + I - 1
                  X(MAI,1) = RHS(I)
                  X(MAI,2) = RHS(MB+I)
  320          CONTINUE
            END IF
C
C           Update right hand sides.
C
C              X(1:MA-1,:) = X(1:MA-1,:) - A(1:MA-1,MA:ME)*X(MA:ME,:)*C'
C
C              X(1:MA-1,:) = X(1:MA-1,:) - E(1:MA-1,MA:ME)*X(MA:ME,:)*D'
C
            IF ( MA .GT. 1 ) THEN
               CALL DGEMM( 'N', 'T', MB, N, N, ONE, X(MA,1), LDX, C,
     $                     LDC, ZERO, TM, 2 )
               CALL DGEMM( 'N', 'N', MA-1, N, MB, MONE, A(1,MA), LDA,
     $                     TM, 2, ONE, X, LDX )
               CALL DGEMM( 'N', 'T', MB, N, N, ONE, X(MA,1), LDX, D,
     $                     LDD, ZERO, TM, 2 )
               CALL DGEMM( 'N', 'N', MA-1, N, MB, MONE, E(1,MA), LDE,
     $                     TM, 2, ONE, X, LDX )
            END IF
C
         GOTO 180
         END IF
C        END WHILE 180
C
      END IF
C
      RETURN
C *** Last line of SG03BW ***
      END
