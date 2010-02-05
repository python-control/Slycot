      SUBROUTINE SG03BX( DICO, TRANS, A, LDA, E, LDE, B, LDB, U, LDU,
     $                   SCALE, M1, LDM1, M2, LDM2, INFO )
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
C     To solve for X = op(U)**T * op(U) either the generalized c-stable
C     continuous-time Lyapunov equation
C
C             T                    T
C        op(A)  * X * op(E) + op(E)  * X * op(A)
C
C                 2        T
C        = - SCALE  * op(B)  * op(B),                                (1)
C
C     or the generalized d-stable discrete-time Lyapunov equation
C
C             T                    T
C        op(A)  * X * op(A) - op(E)  * X * op(E)
C
C                 2        T
C        = - SCALE  * op(B)  * op(B),                                (2)
C
C     where op(K) is either K or K**T for K = A, B, E, U. The Cholesky
C     factor U of the solution is computed without first finding X.
C
C     Furthermore, the auxiliary matrices
C
C                                   -1        -1
C        M1 := op(U) * op(A) * op(E)   * op(U)
C
C                           -1        -1
C        M2 := op(B) * op(E)   * op(U)
C
C     are computed in a numerically reliable way.
C
C     The matrices A, B, E, M1, M2, and U are real 2-by-2 matrices. The
C     pencil A - lambda * E must have a pair of complex conjugate
C     eigenvalues. The eigenvalues must be in the open right half plane
C     (in the continuous-time case) or inside the unit circle (in the
C     discrete-time case).
C
C     The resulting matrix U is upper triangular. The entries on its
C     main diagonal are non-negative. SCALE is an output scale factor
C     set to avoid overflow in U.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies whether the continuous-time or the discrete-time
C             equation is to be solved:
C             = 'C':  Solve continuous-time equation (1);
C             = 'D':  Solve discrete-time equation (2).
C
C     TRANS   CHARACTER*1
C             Specifies whether the transposed equation is to be solved
C             or not:
C             = 'N':  op(K) = K,     K = A, B, E, U;
C             = 'T':  op(K) = K**T,  K = A, B, E, U.
C
C     Input/Output Parameters
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,2)
C             The leading 2-by-2 part of this array must contain the
C             matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= 2.
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,2)
C             The leading 2-by-2 upper triangular part of this array
C             must contain the matrix E.
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= 2.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,2)
C             The leading 2-by-2 upper triangular part of this array
C             must contain the matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= 2.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,2)
C             The leading 2-by-2 part of this array contains the upper
C             triangular matrix U.
C
C     LDU     INTEGER
C             The leading dimension of the array U.  LDU >= 2.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in U.
C             0 < SCALE <= 1.
C
C     M1      (output) DOUBLE PRECISION array, dimension (LDM1,2)
C             The leading 2-by-2 part of this array contains the
C             matrix M1.
C
C     LDM1    INTEGER
C             The leading dimension of the array M1.  LDM1 >= 2.
C
C     M2      (output) DOUBLE PRECISION array, dimension (LDM2,2)
C             The leading 2-by-2 part of this array contains the
C             matrix M2.
C
C     LDM2    INTEGER
C             The leading dimension of the array M2.  LDM2 >= 2.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 2:  the eigenvalues of the pencil A - lambda * E are not
C                   a pair of complex conjugate numbers;
C             = 3:  the eigenvalues of the pencil A - lambda * E are
C                   not in the open right half plane (in the continuous-
C                   time case) or inside the unit circle (in the
C                   discrete-time case).
C
C     METHOD
C
C     The method used by the routine is based on a generalization of the
C     method due to Hammarling ([1], section 6) for Lyapunov equations
C     of order 2. A more detailed description is given in [2].
C
C     REFERENCES
C
C     [1] Hammarling, S.J.
C         Numerical solution of the stable, non-negative definite
C         Lyapunov equation.
C         IMA J. Num. Anal., 2, pp. 303-323, 1982.
C
C     [2] Penzl, T.
C         Numerical solution of generalized Lyapunov equations.
C         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
C
C     FURTHER COMMENTS
C
C     If the solution matrix U is singular, the matrices M1 and M2 are
C     properly set (see [1], equation (6.21)).
C
C     CONTRIBUTOR
C
C     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998.
C
C     REVISIONS
C
C     Sep. 1998 (V. Sima).
C     Dec. 1998 (V. Sima).
C     July 2003 (V. Sima; suggested by Klaus Schnepper).
C     Oct. 2003 (A. Varga).
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  MONE, ONE, TWO, ZERO
      PARAMETER         ( MONE = -1.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                    ZERO = 0.0D+0)
C     .. Scalar Arguments ..
      CHARACTER         DICO, TRANS
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, LDA, LDB, LDE, LDM1, LDM2, LDU
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), E(LDE,*), M1(LDM1,*),
     $                  M2(LDM2,*), U(LDU,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, B11, B12I, B12R, B22, BETAI, BETAR,
     $                  BIGNUM, CI, CR, EPS, L, LAMI, LAMR, SCALE1,
     $                  SCALE2, SI, SMLNUM, SR, T, V, W, XR, XI, YR, YI
      LOGICAL           ISCONT, ISTRNS
C     .. Local Arrays ..
      DOUBLE PRECISION  AA(2,2), AI(2,2), AR(2,2), BB(2,2), BI(2,2),
     $                  BR(2,2), EE(2,2), EI(2,2), ER(2,2), M1I(2,2),
     $                  M1R(2,2), M2I(2,2), M2R(2,2), QBI(2,2),
     $                  QBR(2,2), QI(2,2), QR(2,2), QUI(2,2), QUR(2,2),
     $                  TI(2,2), TR(2,2), UI(2,2), UR(2,2), ZI(2,2),
     $                  ZR(2,2)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLAPY2
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLABAD, DLADIV, DLAG2,
     $                  SG03BY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C
C     Decode input parameters.
C
      ISTRNS = LSAME( TRANS, 'T' )
      ISCONT = LSAME( DICO,  'C' )
C
C     Do not check input parameters for errors.
C
C     Set constants to control overflow.
C
      EPS    = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )/EPS
      BIGNUM = ONE/SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
      INFO  = 0
      SCALE = ONE
C
C     Make copies of A, E, and B.
C
      AA(1,1) = A(1,1)
      AA(2,1) = A(2,1)
      AA(1,2) = A(1,2)
      AA(2,2) = A(2,2)
      EE(1,1) = E(1,1)
      EE(2,1) = ZERO
      EE(1,2) = E(1,2)
      EE(2,2) = E(2,2)
      BB(1,1) = B(1,1)
      BB(2,1) = ZERO
      BB(1,2) = B(1,2)
      BB(2,2) = B(2,2)
C
C     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be
C     solved, transpose the matrices A, E, B with respect to the
C     anti-diagonal. This results in a non-transposed equation.
C
      IF ( ISTRNS ) THEN
         V = AA(1,1)
         AA(1,1) = AA(2,2)
         AA(2,2) = V
         V = EE(1,1)
         EE(1,1) = EE(2,2)
         EE(2,2) = V
         V = BB(1,1)
         BB(1,1) = BB(2,2)
         BB(2,2) = V
      END IF
C
C     Perform QZ-step to transform the pencil A - lambda * E to
C     generalized Schur form. The main diagonal of the Schur factor of E
C     is real and positive.
C
C     Compute eigenvalues (LAMR + LAMI * I, LAMR - LAMI * I).
C
      T = MAX( EPS*MAX( ABS( EE(1,1) ), ABS( EE(1,2) ),
     $                  ABS( EE(2,2) ) ), SMLNUM )
      IF ( MIN( ABS( EE(1,1) ), ABS( EE(2,2) ) ) .LT. T ) THEN
         INFO = 3
         RETURN
      END IF
      CALL DLAG2( AA, 2, EE, 2, SMLNUM*EPS, SCALE1, SCALE2, LAMR,
     $            W, LAMI )
      IF (LAMI .LE. ZERO) THEN
         INFO = 2
         RETURN
      END IF
C
C     Compute right orthogonal transformation matrix Q.
C
      CALL SG03BY( SCALE1*AA(1,1) - EE(1,1)*LAMR, -EE(1,1)*LAMI,
     $             SCALE1*AA(2,1), ZERO, CR, CI, SR, SI, L )
      QR(1,1) =  CR
      QR(1,2) =  SR
      QR(2,1) = -SR
      QR(2,2) =  CR
      QI(1,1) = -CI
      QI(1,2) = -SI
      QI(2,1) = -SI
      QI(2,2) =  CI
C
C     A := Q * A
C
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, QR, 2, AA, 2, ZERO, AR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, QI, 2, AA, 2, ZERO, AI, 2 )
C
C     E := Q * E
C
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, QR, 2, EE, 2, ZERO, ER, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, QI, 2, EE, 2, ZERO, EI, 2 )
C
C     Compute left orthogonal transformation matrix Z.
C
      CALL SG03BY( ER(2,2), EI(2,2), ER(2,1), EI(2,1), CR, CI, SR, SI,
     $             L )
      ZR(1,1) =  CR
      ZR(1,2) =  SR
      ZR(2,1) = -SR
      ZR(2,2) =  CR
      ZI(1,1) =  CI
      ZI(1,2) = -SI
      ZI(2,1) = -SI
      ZI(2,2) = -CI
C
C     E := E * Z
C
      CALL DGEMV( 'T', 2, 2,  ONE, ZR, 2, ER, 2, ZERO, TR, 2 )
      CALL DGEMV( 'T', 2, 2, MONE, ZI, 2, EI, 2,  ONE, TR, 2 )
      CALL DGEMV( 'T', 2, 2,  ONE, ZI, 2, ER, 2, ZERO, TI, 2 )
      CALL DGEMV( 'T', 2, 2,  ONE, ZR, 2, EI, 2,  ONE, TI, 2 )
      CALL DCOPY( 2, TR, 2, ER, 2 )
      CALL DCOPY( 2, TI, 2, EI, 2 )
      ER(2,1) = ZERO
      ER(2,2) = L
      EI(2,1) = ZERO
      EI(2,2) = ZERO
C
C     Make main diagonal entries of E real and positive.
C     (Note:  Z and E are altered.)
C
      V = DLAPY2( ER(1,1), EI(1,1) )
      CALL DLADIV( V, ZERO, ER(1,1), EI(1,1), XR, XI )
      ER(1,1) = V
      EI(1,1) = ZERO
      YR = ZR(1,1)
      YI = ZI(1,1)
      ZR(1,1) = XR*YR - XI*YI
      ZI(1,1) = XR*YI + XI*YR
      YR = ZR(2,1)
      YI = ZI(2,1)
      ZR(2,1) = XR*YR - XI*YI
      ZI(2,1) = XR*YI + XI*YR
C
C     A := A * Z
C
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, AR, 2, ZR, 2, ZERO, TR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, MONE, AI, 2, ZI, 2,  ONE, TR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, AR, 2, ZI, 2, ZERO, TI, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, AI, 2, ZR, 2,  ONE, TI, 2 )
      CALL DCOPY( 4, TR, 1, AR, 1 )
      CALL DCOPY( 4, TI, 1, AI, 1 )
C
C     End of QZ-step.
C
C     B := B * Z
C
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, BB, 2, ZR, 2, ZERO, BR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, ONE, BB, 2, ZI, 2, ZERO, BI, 2 )
C
C     Overwrite B with the upper triangular matrix of its
C     QR-factorization. The elements on the main diagonal are real
C     and non-negative.
C
      CALL SG03BY( BR(1,1), BI(1,1), BR(2,1), BI(2,1), CR, CI, SR, SI,
     $             L )
      QBR(1,1) =  CR
      QBR(1,2) =  SR
      QBR(2,1) = -SR
      QBR(2,2) =  CR
      QBI(1,1) = -CI
      QBI(1,2) = -SI
      QBI(2,1) = -SI
      QBI(2,2) =  CI
      CALL DGEMV( 'N', 2, 2,  ONE, QBR, 2, BR(1,2), 1, ZERO, TR, 1 )
      CALL DGEMV( 'N', 2, 2, MONE, QBI, 2, BI(1,2), 1,  ONE, TR, 1 )
      CALL DGEMV( 'N', 2, 2,  ONE, QBI, 2, BR(1,2), 1, ZERO, TI, 1 )
      CALL DGEMV( 'N', 2, 2,  ONE, QBR, 2, BI(1,2), 1,  ONE, TI, 1 )
      CALL DCOPY( 2, TR, 1, BR(1,2), 1 )
      CALL DCOPY( 2, TI, 1, BI(1,2), 1 )
      BR(1,1) = L
      BR(2,1) = ZERO
      BI(1,1) = ZERO
      BI(2,1) = ZERO
      V = DLAPY2( BR(2,2), BI(2,2) )
      IF ( V .GE. MAX( EPS*MAX( BR(1,1), DLAPY2( BR(1,2), BI(1,2) ) ),
     $                 SMLNUM ) ) THEN
         CALL DLADIV( V, ZERO, BR(2,2), BI(2,2), XR, XI )
         BR(2,2) = V
         YR = QBR(2,1)
         YI = QBI(2,1)
         QBR(2,1) = XR*YR - XI*YI
         QBI(2,1) = XR*YI + XI*YR
         YR = QBR(2,2)
         YI = QBI(2,2)
         QBR(2,2) = XR*YR - XI*YI
         QBI(2,2) = XR*YI + XI*YR
      ELSE
         BR(2,2) = ZERO
      END IF
      BI(2,2) = ZERO
C
C     Compute the Cholesky factor of the solution of the reduced
C     equation. The solution may be scaled to avoid overflow.
C
      IF ( ISCONT ) THEN
C
C        Continuous-time equation.
C
C        Step I:  Compute U(1,1). Set U(2,1) = 0.
C
         V = -TWO*( AR(1,1)*ER(1,1) + AI(1,1)*EI(1,1) )
         IF ( V .LE. ZERO ) THEN
            INFO = 3
            RETURN
         END IF
         V = SQRT( V )
         T = TWO*ABS( BR(1,1) )*SMLNUM
         IF ( T .GT. V ) THEN
            SCALE1  = V/T
            SCALE   = SCALE1*SCALE
            BR(1,1) = SCALE1*BR(1,1)
            BR(1,2) = SCALE1*BR(1,2)
            BI(1,2) = SCALE1*BI(1,2)
            BR(2,2) = SCALE1*BR(2,2)
         END IF
         UR(1,1) = BR(1,1)/V
         UI(1,1) = ZERO
         UR(2,1) = ZERO
         UI(2,1) = ZERO
C
C        Step II:  Compute U(1,2).
C
         T = MAX( EPS*MAX( BR(2,2), DLAPY2( BR(1,2), BI(1,2) ) ),
     $            SMLNUM )
         IF ( ABS( BR(1,1) ) .LT. T ) THEN
            UR(1,2) = ZERO
            UI(1,2) = ZERO
         ELSE
            XR = AR(1,1)*ER(1,2) + AI(1,1)*EI(1,2)
            XI = AI(1,1)*ER(1,2) - AR(1,1)*EI(1,2)
            XR = XR + AR(1,2)*ER(1,1) + AI(1,2)*EI(1,1)
            XI = XI - AI(1,2)*ER(1,1) + AR(1,2)*EI(1,1)
            XR = -BR(1,2)*V - XR*UR(1,1)
            XI =  BI(1,2)*V - XI*UR(1,1)
            YR =  AR(2,2)*ER(1,1) + AI(2,2)*EI(1,1)
            YI = -AI(2,2)*ER(1,1) + AR(2,2)*EI(1,1)
            YR = YR + ER(2,2)*AR(1,1) + EI(2,2)*AI(1,1)
            YI = YI - EI(2,2)*AR(1,1) + ER(2,2)*AI(1,1)
            T  = TWO*DLAPY2( XR, XI )*SMLNUM
            IF ( T .GT. DLAPY2( YR, YI ) ) THEN
               SCALE1  = DLAPY2( YR, YI )/T
               SCALE   = SCALE1*SCALE
               BR(1,1) = SCALE1*BR(1,1)
               BR(1,2) = SCALE1*BR(1,2)
               BI(1,2) = SCALE1*BI(1,2)
               BR(2,2) = SCALE1*BR(2,2)
               UR(1,1) = SCALE1*UR(1,1)
               XR = SCALE1*XR
               XI = SCALE1*XI
            END IF
            CALL DLADIV( XR, XI, YR, YI, UR(1,2), UI(1,2) )
            UI(1,2) = -UI(1,2)
         END IF
C
C        Step III:  Compute U(2,2).
C
         XR = ( ER(1,2)*UR(1,1) + ER(2,2)*UR(1,2) - EI(2,2)*UI(1,2) )*V
         XI = (-EI(1,2)*UR(1,1) - ER(2,2)*UI(1,2) - EI(2,2)*UR(1,2) )*V
         T  = TWO*DLAPY2( XR, XI )*SMLNUM
         IF ( T .GT. DLAPY2( ER(1,1), EI(1,1) ) ) THEN
            SCALE1  = DLAPY2( ER(1,1), EI(1,1) )/T
            SCALE   = SCALE1*SCALE
            UR(1,1) = SCALE1*UR(1,1)
            UR(1,2) = SCALE1*UR(1,2)
            UI(1,2) = SCALE1*UI(1,2)
            BR(1,1) = SCALE1*BR(1,1)
            BR(1,2) = SCALE1*BR(1,2)
            BI(1,2) = SCALE1*BI(1,2)
            BR(2,2) = SCALE1*BR(2,2)
            XR = SCALE1*XR
            XI = SCALE1*XI
         END IF
         CALL DLADIV( XR, XI, ER(1,1), -EI(1,1), YR, YI )
         YR =  BR(1,2) - YR
         YI = -BI(1,2) - YI
         V  = -TWO*( AR(2,2)*ER(2,2) + AI(2,2)*EI(2,2) )
         IF ( V .LE. ZERO ) THEN
            INFO = 3
            RETURN
         END IF
         V = SQRT( V )
         W = DLAPY2( DLAPY2( BR(2,2), BI(2,2) ), DLAPY2( YR, YI ) )
         T = TWO*W*SMLNUM
         IF ( T .GT. V ) THEN
            SCALE1  = V/T
            SCALE   = SCALE1*SCALE
            UR(1,1) = SCALE1*UR(1,1)
            UR(1,2) = SCALE1*UR(1,2)
            UI(1,2) = SCALE1*UI(1,2)
            BR(1,1) = SCALE1*BR(1,1)
            BR(1,2) = SCALE1*BR(1,2)
            BI(1,2) = SCALE1*BI(1,2)
            BR(2,2) = SCALE1*BR(2,2)
            W = SCALE1*W
         END IF
         UR(2,2) = W/V
         UI(2,2) = ZERO
C
C        Compute matrices M1 and M2 for the reduced equation.
C
         M1R(2,1) = ZERO
         M1I(2,1) = ZERO
         M2R(2,1) = ZERO
         M2I(2,1) = ZERO
         CALL DLADIV( AR(1,1), AI(1,1), ER(1,1), EI(1,1), BETAR, BETAI )
         M1R(1,1) =  BETAR
         M1I(1,1) =  BETAI
         M1R(2,2) =  BETAR
         M1I(2,2) = -BETAI
         ALPHA = SQRT( -TWO*BETAR )
         M2R(1,1) = ALPHA
         M2I(1,1) = ZERO
         V  = ER(1,1)*ER(2,2)
         XR = ( -BR(1,1)*ER(1,2) + ER(1,1)*BR(1,2) )/V
         XI = ( -BR(1,1)*EI(1,2) + ER(1,1)*BI(1,2) )/V
         YR =  XR - ALPHA*UR(1,2)
         YI = -XI + ALPHA*UI(1,2)
         IF ( ( YR.NE.ZERO ) .OR. ( YI.NE.ZERO ) ) THEN
            M2R(1,2) =  YR/UR(2,2)
            M2I(1,2) = -YI/UR(2,2)
            M2R(2,2) =  BR(2,2)/( ER(2,2)*UR(2,2) )
            M2I(2,2) =  ZERO
            M1R(1,2) = -ALPHA*M2R(1,2)
            M1I(1,2) = -ALPHA*M2I(1,2)
         ELSE
            M2R(1,2) = ZERO
            M2I(1,2) = ZERO
            M2R(2,2) = ALPHA
            M2I(2,2) = ZERO
            M1R(1,2) = ZERO
            M1I(1,2) = ZERO
         END IF
      ELSE
C
C        Discrete-time equation.
C
C        Step I:  Compute U(1,1). Set U(2,1) = 0.
C
         V = ER(1,1)**2 + EI(1,1)**2 - AR(1,1)**2 - AI(1,1)**2
         IF ( V .LE. ZERO ) THEN
            INFO = 3
            RETURN
         END IF
         V = SQRT( V )
         T = TWO*ABS( BR(1,1) )*SMLNUM
         IF ( T .GT. V ) THEN
            SCALE1  = V/T
            SCALE   = SCALE1*SCALE
            BR(1,1) = SCALE1*BR(1,1)
            BR(1,2) = SCALE1*BR(1,2)
            BI(1,2) = SCALE1*BI(1,2)
            BR(2,2) = SCALE1*BR(2,2)
         END IF
         UR(1,1) = BR(1,1)/V
         UI(1,1) = ZERO
         UR(2,1) = ZERO
         UI(2,1) = ZERO
C
C        Step II:  Compute U(1,2).
C
         T = MAX( EPS*MAX( BR(2,2), DLAPY2( BR(1,2), BI(1,2) ) ),
     $            SMLNUM )
         IF ( ABS( BR(1,1) ) .LT. T ) THEN
            UR(1,2) = ZERO
            UI(1,2) = ZERO
         ELSE
            XR =  AR(1,1)*AR(1,2) + AI(1,1)*AI(1,2)
            XI =  AI(1,1)*AR(1,2) - AR(1,1)*AI(1,2)
            XR =  XR - ER(1,2)*ER(1,1) - EI(1,2)*EI(1,1)
            XI =  XI + EI(1,2)*ER(1,1) - ER(1,2)*EI(1,1)
            XR = -BR(1,2)*V - XR*UR(1,1)
            XI =  BI(1,2)*V - XI*UR(1,1)
            YR =  AR(2,2)*AR(1,1) + AI(2,2)*AI(1,1)
            YI = -AI(2,2)*AR(1,1) + AR(2,2)*AI(1,1)
            YR = YR - ER(2,2)*ER(1,1) - EI(2,2)*EI(1,1)
            YI = YI + EI(2,2)*ER(1,1) - ER(2,2)*EI(1,1)
            T  = TWO*DLAPY2( XR, XI )*SMLNUM
            IF ( T .GT. DLAPY2( YR, YI ) ) THEN
               SCALE1  = DLAPY2( YR, YI )/T
               SCALE   = SCALE1*SCALE
               BR(1,1) = SCALE1*BR(1,1)
               BR(1,2) = SCALE1*BR(1,2)
               BI(1,2) = SCALE1*BI(1,2)
               BR(2,2) = SCALE1*BR(2,2)
               UR(1,1) = SCALE1*UR(1,1)
               XR = SCALE1*XR
               XI = SCALE1*XI
            END IF
            CALL DLADIV( XR, XI, YR, YI, UR(1,2), UI(1,2) )
            UI(1,2) = -UI(1,2)
         END IF
C
C        Step III:  Compute U(2,2).
C
         XR =  ER(1,2)*UR(1,1) + ER(2,2)*UR(1,2) - EI(2,2)*UI(1,2)
         XI = -EI(1,2)*UR(1,1) - ER(2,2)*UI(1,2) - EI(2,2)*UR(1,2)
         YR =  AR(1,2)*UR(1,1) + AR(2,2)*UR(1,2) - AI(2,2)*UI(1,2)
         YI = -AI(1,2)*UR(1,1) - AR(2,2)*UI(1,2) - AI(2,2)*UR(1,2)
         V  = ER(2,2)**2 + EI(2,2)**2 - AR(2,2)**2 - AI(2,2)**2
         IF ( V .LE. ZERO ) THEN
            INFO = 3
            RETURN
         END IF
         V = SQRT( V )
         T = MAX( ABS( BR(2,2) ), ABS( BR(1,2) ), ABS( BI(1,2) ),
     $            ABS( XR ), ABS( XI ), ABS( YR ), ABS( YI) )
         IF ( T .LE. SMLNUM ) T = ONE
         W = ( BR(2,2)/T )**2 + ( BR(1,2)/T )**2 + ( BI(1,2)/T )**2 -
     $       ( XR/T )**2 - ( XI/T )**2 + ( YR/T )**2 + ( YI/T )**2
         IF ( W .LT. ZERO ) THEN
            INFO = 3
            RETURN
         END IF
         W = T*SQRT( W )
         T = TWO*W*SMLNUM
         IF ( T .GT. V ) THEN
            SCALE1  = V/T
            SCALE   = SCALE1*SCALE
            UR(1,1) = SCALE1*UR(1,1)
            UR(1,2) = SCALE1*UR(1,2)
            UI(1,2) = SCALE1*UI(1,2)
            BR(1,1) = SCALE1*BR(1,1)
            BR(1,2) = SCALE1*BR(1,2)
            BI(1,2) = SCALE1*BI(1,2)
            BR(2,2) = SCALE1*BR(2,2)
            W = SCALE1*W
         END IF
         UR(2,2) = W/V
         UI(2,2) = ZERO
C
C        Compute matrices M1 and M2 for the reduced equation.
C
         B11  = BR(1,1)/ER(1,1)
         T    = ER(1,1)*ER(2,2)
         B12R = ( ER(1,1)*BR(1,2) - BR(1,1)*ER(1,2) )/T
         B12I = ( ER(1,1)*BI(1,2) - BR(1,1)*EI(1,2) )/T
         B22  = BR(2,2)/ER(2,2)
         M1R(2,1) = ZERO
         M1I(2,1) = ZERO
         M2R(2,1) = ZERO
         M2I(2,1) = ZERO
         CALL DLADIV( AR(1,1), AI(1,1), ER(1,1), EI(1,1), BETAR, BETAI )
         M1R(1,1) =  BETAR
         M1I(1,1) =  BETAI
         M1R(2,2) =  BETAR
         M1I(2,2) = -BETAI
         V = DLAPY2( BETAR, BETAI )
         ALPHA = SQRT( ( ONE - V )*( ONE + V ) )
         M2R(1,1) = ALPHA
         M2I(1,1) = ZERO
         XR = ( AI(1,1)*EI(1,2) - AR(1,1)*ER(1,2) )/T + AR(1,2)/ER(2,2)
         XI = ( AR(1,1)*EI(1,2) + AI(1,1)*ER(1,2) )/T - AI(1,2)/ER(2,2)
         XR = -TWO*BETAI*B12I - B11*XR
         XI = -TWO*BETAI*B12R - B11*XI
         V  =  ONE + ( BETAI - BETAR )*( BETAI + BETAR )
         W  = -TWO*BETAI*BETAR
         CALL DLADIV( XR, XI, V, W, YR, YI )
         IF ( ( YR.NE.ZERO ) .OR. ( YI.NE.ZERO ) ) THEN
            M2R(1,2) =  ( YR*BETAR - YI*BETAI )/UR(2,2)
            M2I(1,2) = -( YI*BETAR + YR*BETAI )/UR(2,2)
            M2R(2,2) =  B22/UR(2,2)
            M2I(2,2) =  ZERO
            M1R(1,2) = -ALPHA*YR/UR(2,2)
            M1I(1,2) =  ALPHA*YI/UR(2,2)
         ELSE
            M2R(1,2) = ZERO
            M2I(1,2) = ZERO
            M2R(2,2) = ALPHA
            M2I(2,2) = ZERO
            M1R(1,2) = ZERO
            M1I(1,2) = ZERO
         END IF
      END IF
C
C     Transform U back:  U := U * Q.
C     (Note:  Z is used as workspace.)
C
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, UR, 2, QR, 2, ZERO, ZR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, MONE, UI, 2, QI, 2,  ONE, ZR, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, UR, 2, QI, 2, ZERO, ZI, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, UI, 2, QR, 2,  ONE, ZI, 2 )
C
C     Overwrite U with the upper triangular matrix of its
C     QR-factorization. The elements on the main diagonal are real
C     and non-negative.
C
      CALL SG03BY( ZR(1,1), ZI(1,1), ZR(2,1), ZI(2,1), CR, CI, SR, SI,
     $             L )
      QUR(1,1) =  CR
      QUR(1,2) =  SR
      QUR(2,1) = -SR
      QUR(2,2) =  CR
      QUI(1,1) = -CI
      QUI(1,2) = -SI
      QUI(2,1) = -SI
      QUI(2,2) =  CI
      CALL DGEMV( 'N', 2, 2,  ONE, QUR, 2, ZR(1,2), 1, ZERO,  U(1,2), 1)
      CALL DGEMV( 'N', 2, 2, MONE, QUI, 2, ZI(1,2), 1,  ONE,  U(1,2), 1)
      CALL DGEMV( 'N', 2, 2,  ONE, QUI, 2, ZR(1,2), 1, ZERO, UI(1,2), 1)
      CALL DGEMV( 'N', 2, 2,  ONE, QUR, 2, ZI(1,2), 1,  ONE, UI(1,2), 1)
      U(1,1) = L
      U(2,1) = ZERO
      V = DLAPY2( U(2,2), UI(2,2) )
      IF ( V .NE. ZERO ) THEN
         CALL DLADIV( V, ZERO, U(2,2), UI(2,2), XR, XI )
         YR = QUR(2,1)
         YI = QUI(2,1)
         QUR(2,1) = XR*YR - XI*YI
         QUI(2,1) = XR*YI + XI*YR
         YR = QUR(2,2)
         YI = QUI(2,2)
         QUR(2,2) = XR*YR - XI*YI
         QUI(2,2) = XR*YI + XI*YR
      END IF
      U(2,2) = V
C
C     Transform the matrices M1 and M2 back.
C
C        M1 := QU * M1 * QU**H
C        M2 := QB**H * M2 * QU**H
C
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M1R, 2, QUR, 2, ZERO, TR, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M1I, 2, QUI, 2,  ONE, TR, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2, MONE, M1R, 2, QUI, 2, ZERO, TI, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M1I, 2, QUR, 2,  ONE, TI, 2 )
      CALL DGEMM( 'N', 'N', 2, 2, 2,  ONE, QUR, 2, TR,  2, ZERO, M1,
     $            LDM1 )
      CALL DGEMM( 'N', 'N', 2, 2, 2, MONE, QUI, 2, TI,  2,  ONE, M1,
     $            LDM1 )
C
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M2R, 2, QUR, 2, ZERO, TR, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M2I, 2, QUI, 2,  ONE, TR, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2, MONE, M2R, 2, QUI, 2, ZERO, TI, 2 )
      CALL DGEMM( 'N', 'T', 2, 2, 2,  ONE, M2I, 2, QUR, 2,  ONE, TI, 2 )
      CALL DGEMM( 'T', 'N', 2, 2, 2,  ONE, QBR, 2,  TR, 2, ZERO, M2,
     $            LDM2 )
      CALL DGEMM( 'T', 'N', 2, 2, 2,  ONE, QBI, 2,  TI, 2,  ONE, M2,
     $            LDM2 )
C
C     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be
C     solved, transpose the matrix U with respect to the
C     anti-diagonal and the matrices M1, M2 with respect to the diagonal
C     and the anti-diagonal.
C
      IF ( ISTRNS ) THEN
         V = U(1,1)
         U(1,1) = U(2,2)
         U(2,2) = V
         V = M1(1,1)
         M1(1,1) = M1(2,2)
         M1(2,2) = V
         V = M2(1,1)
         M2(1,1) = M2(2,2)
         M2(2,2) = V
      END IF
C
      RETURN
C *** Last line of SG03BX ***
      END
