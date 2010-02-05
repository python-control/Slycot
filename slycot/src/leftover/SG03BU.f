      SUBROUTINE SG03BU( TRANS, N, A, LDA, E, LDE, B, LDB, SCALE,
     $                   DWORK, INFO )
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
C     To compute the Cholesky factor U of the matrix X, X = U**T * U or
C     X = U * U**T, which is the solution of the generalized d-stable
C     discrete-time Lyapunov equation
C
C         T            T                  2    T
C        A  * X * A - E  * X * E = - SCALE  * B  * B,                (1)
C
C     or the transposed equation
C
C                 T            T          2        T
C        A * X * A  - E * X * E  = - SCALE  * B * B ,                (2)
C
C     respectively, where A, E, B, and U are real N-by-N matrices. The
C     Cholesky factor U of the solution is computed without first
C     finding X. The pencil A - lambda * E must be in generalized Schur
C     form ( A upper quasitriangular, E upper triangular ). Moreover, it
C     must be d-stable, i.e. the moduli of its eigenvalues must be less
C     than one. B must be an upper triangular matrix with non-negative
C     entries on its main diagonal.
C
C     The resulting matrix U is upper triangular. The entries on its
C     main diagonal are non-negative. SCALE is an output scale factor
C     set to avoid overflow in U.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANS   CHARACTER*1
C             Specifies whether equation (1) or equation (2) is to be
C             solved:
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
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N upper triangular part of this
C             array must contain the matrix B.
C             On exit, the leading N-by-N upper triangular part of this
C             array contains the solution matrix U.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in U.
C             0 < SCALE <= 1.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (6*N-6)
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the generalized Sylvester equation to be solved in
C                   step II (see METHOD) is (nearly) singular to working
C                   precision;  perturbed values were used to solve the
C                   equation (but the matrices A and E are unchanged);
C             = 2:  the generalized Schur form of the pencil
C                   A - lambda * E contains a 2-by-2 main diagonal block
C                   whose eigenvalues are not a pair of conjugate
C                   complex numbers;
C             = 3:  the pencil A - lambda * E is not d-stable, i.e.
C                   there are eigenvalues outside the open unit circle;
C             = 4:  the LAPACK routine DSYEVX utilized to factorize M3
C                   failed to converge. This error is unlikely to occur.
C
C     METHOD
C
C     The method [2] used by the routine is an extension of Hammarling's
C     algorithm [1] to generalized Lyapunov equations.
C
C     We present the method for solving equation (1). Equation (2) can
C     be treated in a similar fashion. For simplicity, assume SCALE = 1.
C
C     The matrix A is an upper quasitriangular matrix, i.e. it is a
C     block triangular matrix with square blocks on the main diagonal
C     and the block order at most 2. We use the following partitioning
C     for the matrices A, E, B and the solution matrix U
C
C               ( A11   A12 )        ( E11   E12 )
C           A = (           ),   E = (           ),
C               (   0   A22 )        (   0   E22 )
C
C               ( B11   B12 )        ( U11   U12 )
C           B = (           ),   U = (           ).                  (3)
C               (   0   B22 )        (   0   U22 )
C
C     The size of the (1,1)-blocks is 1-by-1 (iff A(2,1) = 0.0) or
C     2-by-2.
C
C     We compute U11 and U12**T in three steps.
C
C     Step I:
C
C        From (1) and (3) we get the 1-by-1 or 2-by-2 equation
C
C                T      T                   T      T
C             A11  * U11  * U11 * A11  - E11  * U11  * U11 * E11
C
C                    T
C             = - B11  * B11.
C
C        For brevity, details are omitted here. The technique for
C        computing U11 is similar to those applied to standard Lyapunov
C        equations in Hammarling's algorithm ([1], section 6).
C
C        Furthermore, the auxiliary matrices M1 and M2 defined as
C        follows
C
C                               -1      -1
C           M1 = U11 * A11 * E11   * U11
C
C                         -1      -1
C           M2 = B11 * E11   * U11
C
C        are computed in a numerically reliable way.
C
C     Step II:
C
C        We solve for U12**T the generalized Sylvester equation
C
C              T      T           T      T
C           A22  * U12  * M1 - E22  * U12
C
C                  T           T      T      T      T
C           = - B12  * M2 + E12  * U11  - A12  * U11  * M1.
C
C     Step III:
C
C        One can show that
C
C              T      T                  T      T
C           A22  * U22  * U22 * A22 - E22  * U22  * U22 * E22  =
C
C                T              T
C           - B22  * B22 - y * y                                     (4)
C
C        holds, where y is defined as follows
C
C                  T      T      T      T
C           w = A12  * U11  + A22  * U12
C
C                    T
C           y = ( B12   w ) * M3EV,
C
C        where M3EV is a matrix which fulfils
C
C                ( I-M2*M2**T   -M2*M1**T )              T
C           M3 = (                        ) = M3EV * M3EV .
C                (  -M1*M2**T  I-M1*M1**T )
C
C        M3 is positive semidefinite and its rank is equal to the size
C        of U11. Therefore, a matrix M3EV can be found by solving the
C        symmetric eigenvalue problem for M3 such that y consists of
C        either 1 or 2 rows.
C
C        If B22_tilde is the square triangular matrix arising from the
C        QR-factorization
C
C               ( B22_tilde )     ( B22  )
C           Q * (           )  =  (      ),
C               (     0     )     ( y**T )
C
C        then
C
C                T              T                T
C           - B22  * B22 - y * y   =  - B22_tilde  * B22_tilde.
C
C        Replacing the right hand side in (4) by the term
C        - B22_tilde**T * B22_tilde leads to a generalized Lyapunov
C        equation of lower dimension compared to (1).
C
C     The solution U of the equation (1) can be obtained by recursive
C     application of the steps I to III.
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
C     NUMERICAL ASPECTS
C
C     The routine requires 2*N**3 flops. Note that we count a single
C     floating point arithmetic operation as one flop.
C
C     FURTHER COMMENTS
C
C     The Lyapunov equation may be very ill-conditioned. In particular,
C     if the pencil A - lambda * E has a pair of almost reciprocal
C     eigenvalues, then the Lyapunov equation will be ill-conditioned.
C     Perturbed values were used to solve the equation.
C     A condition estimate can be obtained from the routine SG03AD.
C     When setting the error indicator INFO, the routine does not test
C     for near instability in the equation but only for exact
C     instability.
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
      DOUBLE PRECISION  HALF, MONE, ONE, TWO, ZERO
      PARAMETER         ( HALF = 0.5D+0, MONE = -1.0D0, ONE = 1.0D+0,
     $                    TWO = 2.0D+0,  ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANS
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, LDA, LDB, LDE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), E(LDE,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGNUM, C, DELTA1, EPS, S, SCALE1, SMLNUM, UFLT,
     $                  X, Z
      INTEGER           I, INFO1, J, KB, KH, KL, LDWS, M, UIIPT, WPT,
     $                  YPT
      LOGICAL           NOTRNS
C     .. Local Arrays ..
      DOUBLE PRECISION  M1(2,2), M2(2,2), M3(4,4), M3C(4,4), M3EW(4),
     $                  RW(32), TM(2,2), UI(2,2)
      INTEGER           IW(24)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLABAD, DLACPY, DLASET,
     $                  DROT, DROTG, DSCAL, DSYEVX, DSYRK, SG03BW,
     $                  SG03BX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
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
      ELSEIF ( LDB .LT. MAX( 1, N ) ) THEN
         INFO = -8
      ELSE
         INFO = 0
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SG03BU', -INFO )
         RETURN
      END IF
C
      SCALE = ONE
C
C     Quick return if possible.
C
      IF ( N .EQ. 0 )
     $    RETURN
C
C     Set constants to control overflow.
C
      EPS  = DLAMCH( 'P' )
      UFLT = DLAMCH( 'S' )
      SMLNUM = UFLT/EPS
      BIGNUM = ONE/SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
C     Set work space pointers and leading dimension of matrices in
C     work space.
C
      UIIPT = 1
      WPT = 2*N-1
      YPT = 4*N-3
      LDWS = N-1
C
      IF ( NOTRNS ) THEN
C
C        Solve equation (1).
C
C        Main Loop. Compute block row U(KL:KH,KL:N). KB denotes the
C        number of rows in this block row.
C
         KH = 0
C        WHILE ( KH .LT. N ) DO
   20    IF ( KH .LT. N ) THEN
            KL = KH + 1
            IF ( KL .EQ. N ) THEN
               KH = N
               KB = 1
            ELSE
               IF ( A(KL+1,KL) .EQ. ZERO ) THEN
                  KH = KL
                  KB = 1
               ELSE
                  KH = KL + 1
                  KB = 2
               END IF
            END IF
C
C           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary
C                   matrices M1 and M2. (For the moment the result
C                   U(KL:KH,KL:KH) is stored in UI).
C
            IF ( KB .EQ. 1 ) THEN
               DELTA1 = E(KL,KL)**2 - A(KL,KL)**2
               IF ( DELTA1 .LE. ZERO ) THEN
                  INFO = 3
                  RETURN
               END IF
               DELTA1 = SQRT( DELTA1 )
               Z = TWO*ABS( B(KL,KL) )*SMLNUM
               IF ( Z .GT. DELTA1 ) THEN
                  SCALE1 = DELTA1/Z
                  SCALE = SCALE1*SCALE
                  DO 40 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
   40             CONTINUE
               END IF
               UI(1,1) = B(KL,KL)/DELTA1
               M1(1,1) = A(KL,KL)/E(KL,KL)
               M2(1,1) = DELTA1/E(KL,KL)
            ELSE
C
C              If a pair of complex conjugate eigenvalues occurs, apply
C              (complex) Hammarling algorithm for the 2-by-2 problem.
C
               CALL SG03BX( 'D', 'N', A(KL,KL), LDA, E(KL,KL), LDE,
     $                      B(KL,KL), LDB, UI, 2, SCALE1, M1, 2, M2, 2,
     $                      INFO1 )
               IF ( INFO1 .NE. 0 ) THEN
                  INFO = INFO1
                  RETURN
               END IF
               IF ( SCALE1 .NE. ONE ) THEN
                  SCALE = SCALE1*SCALE
                  DO 60 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
   60             CONTINUE
               END IF
            END IF
C
            IF ( KH .LT. N ) THEN
C
C              STEP II: Compute U(KL:KH,KH+1:N) by solving a generalized
C                       Sylvester equation. (For the moment the result
C                       U(KL:KH,KH+1:N) is stored in the workspace.)
C
C              Form right hand side of the Sylvester equation.
C
               CALL DGEMM( 'T', 'N', N-KH, KB, KB, MONE, B(KL,KH+1),
     $                     LDB, M2, 2, ZERO, DWORK(UIIPT), LDWS )
               CALL DGEMM( 'T', 'T', N-KH, KB, KB, ONE, E(KL,KH+1),
     $                     LDE, UI, 2, ONE, DWORK(UIIPT), LDWS )
               CALL DGEMM( 'T', 'N', KB, KB, KB, ONE, UI, 2, M1, 2,
     $                     ZERO, TM, 2 )
               CALL DGEMM( 'T', 'N', N-KH, KB, KB, MONE, A(KL,KH+1),
     $                     LDA, TM, 2, ONE, DWORK(UIIPT), LDWS )
C
C              Solve generalized Sylvester equation.
C
               CALL DLASET( 'A', KB, KB, ZERO, MONE, TM, 2 )
               CALL SG03BW( 'N', N-KH, KB, A(KH+1,KH+1), LDA, M1, 2,
     $                      E(KH+1,KH+1), LDE, TM, 2, DWORK(UIIPT),
     $                      LDWS, SCALE1, INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 1
               IF ( SCALE1 .NE. ONE ) THEN
                  SCALE = SCALE1*SCALE
                  DO 80 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
   80             CONTINUE
                  CALL DSCAL( 4, SCALE1, UI(1,1), 1 )
               END IF
C
C              STEP III: Form the right hand side matrix
C                        B(KH+1:N,KH+1:N) of the (smaller) Lyapunov
C                        equation to be solved during the next pass of
C                        the main loop.
C
C              Compute auxiliary matrices M3 and Y. The factorization
C              M3 = M3C * M3C**T is found by solving the symmetric
C              eigenvalue problem.
C
               CALL DLASET( 'U', 2*KB, 2*KB, ZERO, ONE, M3, 4 )
               CALL DSYRK(  'U', 'N', KB, KB, MONE, M2, 2, ONE, M3, 4 )
               CALL DGEMM(  'N', 'T', KB, KB, KB, MONE, M2, 2, M1, 2,
     $                      ZERO, M3(1,KB+1), 4 )
               CALL DSYRK(  'U', 'N', KB, KB, MONE, M1, 2, ONE,
     $                      M3(KB+1,KB+1), 4 )
               CALL DSYEVX( 'V', 'V', 'U', 2*KB, M3, 4, HALF, TWO, 1, 4,
     $                      TWO*UFLT, M, M3EW, M3C, 4, RW, 32, IW(5),
     $                      IW, INFO1 )
               IF ( INFO1 .NE. 0 ) THEN
                  INFO = 4
                  RETURN
               END IF
               CALL DGEMM( 'T', 'N', N-KH, KB, KB, ONE, B(KL,KH+1), LDB,
     $                     M3C, 4, ZERO, DWORK(YPT), LDWS )
               CALL DGEMM( 'T', 'T', N-KH, KB, KB, ONE, A(KL,KH+1), LDA,
     $                     UI, 2, ZERO, DWORK(WPT), LDWS )
               DO 100 I = 1, N-KH
                  CALL DGEMV( 'T', MIN( I+1, N-KH ), KB, ONE,
     $                        DWORK(UIIPT), LDWS, A(KH+1,KH+I), 1, ONE,
     $                        DWORK(WPT+I-1), LDWS )
  100          CONTINUE
               CALL DGEMM( 'N', 'N', N-KH, KB, KB, ONE, DWORK(WPT),
     $                     LDWS, M3C(KB+1,1), 4, ONE, DWORK(YPT), LDWS )
C
C              Overwrite B(KH+1:N,KH+1:N) with the triangular matrix
C              from the QR-factorization of the (N-KH+KB)-by-(N-KH)
C              matrix
C
C                          (  B(KH+1:N,KH+1:N)  )
C                          (                    )
C                          (       Y**T         ) .
C
               DO 140 J = 1, KB
                  DO 120 I = 1, N-KH
                     X = B(KH+I,KH+I)
                     Z = DWORK(YPT+I-1+(J-1)*LDWS)
                     CALL DROTG( X, Z, C, S )
                     CALL DROT(  N-KH-I+1, B(KH+I,KH+I), LDB,
     $                           DWORK(YPT+I-1+(J-1)*LDWS), 1, C, S )
  120             CONTINUE
  140          CONTINUE
C
C              Make main diagonal elements of B(KH+1:N,KH+1:N) positive.
C
               DO 160 I = KH+1, N
                  IF ( B(I,I) .LT. ZERO )
     $               CALL DSCAL( N-I+1, MONE, B(I,I), LDB )
  160          CONTINUE
C
C              Overwrite right hand side with the part of the solution
C              computed in step II.
C
               DO 180 J = KL, KH
                  CALL DCOPY( N-KH, DWORK(UIIPT+(J-KL)*LDWS), 1,
     $                        B(J,KH+1), LDB )
  180          CONTINUE
            END IF
C
C           Overwrite right hand side with the part of the solution
C           computed in step I.
C
            CALL DLACPY( 'U', KB, KB, UI, 2, B(KL,KL), LDB )
C
         GOTO 20
         END IF
C        END WHILE 20
C
      ELSE
C
C        Solve equation (2).
C
C        Main Loop. Compute block column U(1:KH,KL:KH). KB denotes the
C        number of columns in this block column.
C
         KL = N + 1
C        WHILE ( KL .GT. 1 ) DO
  200    IF ( KL .GT. 1 ) THEN
            KH = KL - 1
            IF ( KH .EQ. 1 ) THEN
               KL = 1
               KB = 1
            ELSE
               IF ( A(KH,KH-1) .EQ. ZERO ) THEN
                  KL = KH
                  KB = 1
               ELSE
                  KL = KH - 1
                  KB = 2
               END IF
            END IF
C
C           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary
C                   matrices M1 and M2. (For the moment the result
C                   U(KL:KH,KL:KH) is stored in UI).
C
            IF ( KB .EQ. 1 ) THEN
               DELTA1 = E(KL,KL)**2 - A(KL,KL)**2
               IF ( DELTA1 .LE. ZERO ) THEN
                  INFO = 3
                  RETURN
               END IF
               DELTA1 = SQRT( DELTA1 )
               Z = TWO*ABS( B(KL,KL) )*SMLNUM
               IF ( Z .GT. DELTA1 ) THEN
                  SCALE1 = DELTA1/Z
                  SCALE = SCALE1*SCALE
                  DO 220 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
  220             CONTINUE
               END IF
               UI(1,1) = B(KL,KL)/DELTA1
               M1(1,1) = A(KL,KL)/E(KL,KL)
               M2(1,1) = DELTA1/E(KL,KL)
            ELSE
C
C              If a pair of complex conjugate eigenvalues occurs, apply
C              (complex) Hammarling algorithm for the 2-by-2 problem.
C
               CALL SG03BX( 'D', 'T', A(KL,KL), LDA, E(KL,KL), LDE,
     $                      B(KL,KL), LDB, UI, 2, SCALE1, M1, 2, M2, 2,
     $                      INFO1 )
               IF ( INFO1 .NE. 0 ) THEN
                  INFO = INFO1
                  RETURN
               END IF
               IF ( SCALE1 .NE. ONE ) THEN
                  SCALE = SCALE1*SCALE
                  DO 240 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
  240             CONTINUE
               END IF
            END IF
C
            IF ( KL .GT. 1 ) THEN
C
C              STEP II: Compute U(1:KL-1,KL:KH) by solving a generalized
C                       Sylvester equation. (For the moment the result
C                       U(1:KL-1,KL:KH) is stored in the workspace.)
C
C              Form right hand side of the Sylvester equation.
C
               CALL DGEMM( 'N', 'T', KL-1, KB, KB, MONE, B(1,KL), LDB,
     $                     M2, 2, ZERO, DWORK(UIIPT), LDWS )
               CALL DGEMM( 'N', 'N', KL-1, KB, KB, ONE, E(1,KL), LDE,
     $                     UI, 2, ONE, DWORK(UIIPT), LDWS )
               CALL DGEMM( 'N', 'T', KB, KB, KB, ONE, UI, 2, M1, 2,
     $                     ZERO, TM, 2 )
               CALL DGEMM( 'N', 'N', KL-1, KB, KB, MONE, A(1,KL), LDA,
     $                     TM, 2, ONE, DWORK(UIIPT), LDWS )
C
C              Solve generalized Sylvester equation.
C
               CALL DLASET( 'A', KB, KB, ZERO, MONE, TM, 2 )
               CALL SG03BW( 'T', KL-1, KB, A, LDA, M1, 2, E, LDE, TM, 2,
     $                      DWORK(UIIPT), LDWS, SCALE1, INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 1
               IF ( SCALE1 .NE. ONE ) THEN
                  SCALE = SCALE1*SCALE
                  DO 260 I = 1, N
                     CALL DSCAL( I, SCALE1, B(1,I), 1 )
  260             CONTINUE
                  CALL DSCAL( 4, SCALE1, UI(1,1), 1 )
               END IF
C
C              STEP III: Form the right hand side matrix
C                        B(1:KL-1,1:KL-1) of the (smaller) Lyapunov
C                        equation to be solved during the next pass of
C                        the main loop.
C
C              Compute auxiliary matrices M3 and Y. The factorization
C              M3 = M3C * M3C**T is found by solving the symmetric
C              eigenvalue problem.
C
               CALL DLASET( 'U', 2*KB, 2*KB, ZERO, ONE, M3, 4 )
               CALL DSYRK(  'U', 'T', KB, KB, MONE, M2, 2, ONE, M3, 4 )
               CALL DGEMM(  'T', 'N', KB, KB, KB, MONE, M2, 2, M1, 2,
     $                      ZERO, M3(1,KB+1), 4 )
               CALL DSYRK(  'U', 'T', KB, KB, MONE, M1, 2, ONE,
     $                      M3(KB+1,KB+1), 4 )
               CALL DSYEVX( 'V', 'V', 'U', 2*KB, M3, 4, HALF, TWO, 1, 4,
     $                      TWO*UFLT, M, M3EW, M3C, 4, RW, 32, IW(5),
     $                      IW, INFO1 )
               IF ( INFO1 .NE. 0 ) THEN
                  INFO = 4
                  RETURN
               END IF
               CALL DGEMM( 'N', 'N', KL-1, KB, KB, ONE, B(1,KL), LDB,
     $                     M3C, 4, ZERO, DWORK(YPT), LDWS )
               CALL DGEMM( 'N', 'N', KL-1, KB, KB, ONE, A(1,KL), LDA,
     $                     UI, 2, ZERO, DWORK(WPT), LDWS )
               DO 280 I = 1, KL-1
                  CALL DGEMV( 'T', MIN( KL-I+1, KL-1 ), KB, ONE,
     $                        DWORK(MAX( UIIPT, UIIPT+I-2 )), LDWS,
     $                        A(I,MAX( I-1, 1 )), LDA, ONE,
     $                        DWORK(WPT+I-1), LDWS )
  280          CONTINUE
               CALL DGEMM( 'N', 'N', KL-1, KB, KB, ONE, DWORK(WPT),
     $                     LDWS, M3C(KB+1,1), 4, ONE, DWORK(YPT), LDWS )
C
C              Overwrite B(1:KL-1,1:KL-1) with the triangular matrix
C              from the RQ-factorization of the (KL-1)-by-KH matrix
C
C                          (                        )
C                          (  B(1:KL-1,1:KL-1)   Y  )
C                          (                        ).
C
               DO 320 J = 1, KB
                  DO 300 I = KL-1, 1, -1
                     X = B(I,I)
                     Z = DWORK(YPT+I-1+(J-1)*LDWS)
                     CALL DROTG( X, Z, C, S )
                     CALL DROT(  I, B(1,I), 1, DWORK(YPT+(J-1)*LDWS), 1,
     $                           C, S )
  300             CONTINUE
  320          CONTINUE
C
C              Make main diagonal elements of B(1:KL-1,1:KL-1) positive.
C
               DO 340 I = 1, KL-1
                  IF ( B(I,I) .LT. ZERO )
     $               CALL DSCAL( I, MONE, B(1,I), 1 )
  340          CONTINUE
C
C              Overwrite right hand side with the part of the solution
C              computed in step II.
C
               CALL DLACPY( 'A', KL-1, KB, DWORK(UIIPT), LDWS, B(1,KL),
     $                      LDB )
C
            END IF
C
C           Overwrite right hand side with the part of the solution
C           computed in step I.
C
            CALL DLACPY( 'U', KB, KB, UI, 2, B(KL,KL), LDB )
C
         GOTO 200
         END IF
C        END WHILE 200
C
      END IF
C
      RETURN
C *** Last line of SG03BU ***
      END
