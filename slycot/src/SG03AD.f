      SUBROUTINE SG03AD( DICO, JOB, FACT, TRANS, UPLO, N, A, LDA, E,
     $                   LDE, Q, LDQ, Z, LDZ, X, LDX, SCALE, SEP, FERR,
     $                   ALPHAR, ALPHAI, BETA, IWORK, DWORK, LDWORK,
     $                   INFO )
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
C     To solve for X either the generalized continuous-time Lyapunov
C     equation
C
C             T                T
C        op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1)
C
C     or the generalized discrete-time Lyapunov equation
C
C             T                T
C        op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2)
C
C     where op(M) is either M or M**T for M = A, E and the right hand
C     side Y is symmetric. A, E, Y, and the solution X are N-by-N
C     matrices. SCALE is an output scale factor, set to avoid overflow
C     in X.
C
C     Estimates of the separation and the relative forward error norm
C     are provided.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies which type of the equation is considered:
C             = 'C':  Continuous-time equation (1);
C             = 'D':  Discrete-time equation (2).
C
C     JOB     CHARACTER*1
C             Specifies if the solution is to be computed and if the
C             separation is to be estimated:
C             = 'X':  Compute the solution only;
C             = 'S':  Estimate the separation only;
C             = 'B':  Compute the solution and estimate the separation.
C
C     FACT    CHARACTER*1
C             Specifies whether the generalized real Schur
C             factorization of the pencil A - lambda * E is supplied
C             on entry or not:
C             = 'N':  Factorization is not supplied;
C             = 'F':  Factorization is supplied.
C
C     TRANS   CHARACTER*1
C             Specifies whether the transposed equation is to be solved
C             or not:
C             = 'N':  op(A) = A,    op(E) = E;
C             = 'T':  op(A) = A**T, op(E) = E**T.
C
C     UPLO    CHARACTER*1
C             Specifies whether the lower or the upper triangle of the
C             array X is needed on input:
C             = 'L':  Only the lower triangle is needed on input;
C             = 'U':  Only the upper triangle is needed on input.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, if FACT = 'F', then the leading N-by-N upper
C             Hessenberg part of this array must contain the
C             generalized Schur factor A_s of the matrix A (see
C             definition (3) in section METHOD). A_s must be an upper
C             quasitriangular matrix. The elements below the upper
C             Hessenberg part of the array A are not referenced.
C             If FACT = 'N', then the leading N-by-N part of this
C             array must contain the matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the generalized Schur factor A_s of the matrix A. (A_s is
C             an upper quasitriangular matrix.)
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, if FACT = 'F', then the leading N-by-N upper
C             triangular part of this array must contain the
C             generalized Schur factor E_s of the matrix E (see
C             definition (4) in section METHOD). The elements below the
C             upper triangular part of the array E are not referenced.
C             If FACT = 'N', then the leading N-by-N part of this
C             array must contain the coefficient matrix E of the
C             equation.
C             On exit, the leading N-by-N part of this array contains
C             the generalized Schur factor E_s of the matrix E. (E_s is
C             an upper triangular matrix.)
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, if FACT = 'F', then the leading N-by-N part of
C             this array must contain the orthogonal matrix Q from
C             the generalized Schur factorization (see definitions (3)
C             and (4) in section METHOD).
C             If FACT = 'N', Q need not be set on entry.
C             On exit, the leading N-by-N part of this array contains
C             the orthogonal matrix Q from the generalized Schur
C             factorization.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, if FACT = 'F', then the leading N-by-N part of
C             this array must contain the orthogonal matrix Z from
C             the generalized Schur factorization (see definitions (3)
C             and (4) in section METHOD).
C             If FACT = 'N', Z need not be set on entry.
C             On exit, the leading N-by-N part of this array contains
C             the orthogonal matrix Z from the generalized Schur
C             factorization.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, if JOB = 'B' or 'X', then the leading N-by-N
C             part of this array must contain the right hand side matrix
C             Y of the equation. Either the lower or the upper
C             triangular part of this array is needed (see mode
C             parameter UPLO).
C             If JOB = 'S', X is not referenced.
C             On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then
C             the leading N-by-N part of this array contains the
C             solution matrix X of the equation.
C             If JOB = 'S', X is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= MAX(1,N).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor set to avoid overflow in X.
C             (0 < SCALE <= 1)
C
C     SEP     (output) DOUBLE PRECISION
C             If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then
C             SEP contains an estimate of the separation of the
C             Lyapunov operator.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an
C             estimated forward error bound for the solution X. If XTRUE
C             is the true solution, FERR estimates the relative error
C             in the computed solution, measured in the Frobenius norm:
C             norm(X - XTRUE) / norm(XTRUE)
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C             If FACT = 'N' and INFO = 0, 3, or 4, then
C             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the
C             eigenvalues of the matrix pencil A - lambda * E.
C             If FACT = 'F', ALPHAR, ALPHAI, and BETA are not
C             referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N**2)
C             IWORK is not referenced if JOB = 'X'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK. The following table
C             contains the minimal work space requirements depending
C             on the choice of JOB and FACT.
C
C                    JOB        FACT    |  LDWORK
C                    -------------------+-------------------
C                    'X'        'F'     |  MAX(1,N)
C                    'X'        'N'     |  MAX(1,4*N,8*N+16)
C                    'B', 'S'   'F'     |  MAX(1,2*N**2)
C                    'B', 'S'   'N'     |  MAX(1,2*N**2,4*N,8*N+16)
C
C             For optimum performance, LDWORK should be larger.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  FACT = 'F' and the matrix contained in the upper
C                   Hessenberg part of the array A is not in upper
C                   quasitriangular form;
C             = 2:  FACT = 'N' and the pencil A - lambda * E cannot be
C                   reduced to generalized Schur form: LAPACK routine
C                   DGGES has failed to converge;
C             = 3:  DICO = 'D' and the pencil A - lambda * E has a
C                   pair of reciprocal eigenvalues. That is, lambda_i =
C                   1/lambda_j for some i and j, where lambda_i and
C                   lambda_j are eigenvalues of A - lambda * E. Hence,
C                   equation (2) is singular;  perturbed values were
C                   used to solve the equation (but the matrices A and
C                   E are unchanged);
C             = 4:  DICO = 'C' and the pencil A - lambda * E has a
C                   degenerate pair of eigenvalues. That is, lambda_i =
C                   -lambda_j for some i and j, where lambda_i and
C                   lambda_j are eigenvalues of A - lambda * E. Hence,
C                   equation (1) is singular;  perturbed values were
C                   used to solve the equation (but the matrices A and
C                   E are unchanged).
C
C     METHOD
C
C     A straightforward generalization [3] of the method proposed by
C     Bartels and Stewart [1] is utilized to solve (1) or (2).
C
C     First the pencil A - lambda * E is reduced to real generalized
C     Schur form A_s - lambda * E_s by means of orthogonal
C     transformations (QZ-algorithm):
C
C        A_s = Q**T * A * Z   (upper quasitriangular)                (3)
C
C        E_s = Q**T * E * Z   (upper triangular).                    (4)
C
C     If FACT = 'F', this step is omitted. Assuming SCALE = 1 and
C     defining
C
C              ( Z**T * Y * Z   :   TRANS = 'N'
C        Y_s = <
C              ( Q**T * Y * Q   :   TRANS = 'T'
C
C
C              ( Q**T * X * Q    if TRANS = 'N'
C        X_s = <                                                     (5)
C              ( Z**T * X * Z    if TRANS = 'T'
C
C     leads to the reduced Lyapunov equation
C
C               T                      T
C        op(A_s)  X_s op(E_s) + op(E_s)  X_s op(A_s) = Y_s,          (6)
C
C     or
C               T                      T
C        op(A_s)  X_s op(A_s) - op(E_s)  X_s op(E_s) = Y_s,          (7)
C
C     which are equivalent to (1) or (2), respectively. The solution X_s
C     of (6) or (7) is computed via block back substitution (if TRANS =
C     'N') or block forward substitution (if TRANS = 'T'), where the
C     block order is at most 2. (See [1] and [3] for details.)
C     Equation (5) yields the solution matrix X.
C
C     For fast computation the estimates of the separation and the
C     forward error are gained from (6) or (7) rather than (1) or
C     (2), respectively. We consider (6) and (7) as special cases of the
C     generalized Sylvester equation
C
C        R * X * S + U * X * V = Y,                                  (8)
C
C     whose separation is defined as follows
C
C        sep = sep(R,S,U,V) =   min   || R * X * S + U * X * V || .
C                            ||X|| = 1                           F
C                                 F
C
C     Equation (8) is equivalent to the system of linear equations
C
C        K * vec(X) = (kron(S**T,R) + kron(V**T,U)) * vec(X) = vec(Y),
C
C     where kron is the Kronecker product of two matrices and vec
C     is the mapping that stacks the columns of a matrix. If K is
C     nonsingular then
C
C        sep = 1 / ||K**(-1)|| .
C                             2
C
C     We estimate ||K**(-1)|| by a method devised by Higham [2]. Note
C     that this method yields an estimation for the 1-norm but we use it
C     as an approximation for the 2-norm. Estimates for the forward
C     error norm are provided by
C
C        FERR = 2 * EPS * ||A_s||  * ||E_s||  / sep
C                                F          F
C
C     in the continuous-time case (1) and
C
C        FERR = EPS * ( ||A_s|| **2 + ||E_s|| **2 ) / sep
C                              F             F
C
C     in the discrete-time case (2).
C     The reciprocal condition number, RCOND, of the Lyapunov equation
C     can be estimated by FERR/EPS.
C
C     REFERENCES
C
C     [1] Bartels, R.H., Stewart, G.W.
C         Solution of the equation A X + X B = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Higham, N.J.
C         FORTRAN codes for estimating the one-norm of a real or complex
C         matrix, with applications to condition estimation.
C         A.C.M. Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, 1988.
C
C     [3] Penzl, T.
C         Numerical solution of generalized Lyapunov equations.
C         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
C
C     NUMERICAL ASPECTS
C
C     The number of flops required by the routine is given by the
C     following table. Note that we count a single floating point
C     arithmetic operation as one flop. c is an integer number of modest
C     size (say 4 or 5).
C
C                   |  FACT = 'F'            FACT = 'N'
C        -----------+------------------------------------------
C        JOB = 'B'  |  (26+8*c)/3 * N**3     (224+8*c)/3 * N**3
C        JOB = 'S'  |  8*c/3 * N**3          (198+8*c)/3 * N**3
C        JOB = 'X'  |  26/3 * N**3           224/3 * N**3
C
C     The algorithm is backward stable if the eigenvalues of the pencil
C     A - lambda * E are real. Otherwise, linear systems of order at
C     most 4 are involved into the computation. These systems are solved
C     by Gauss elimination with complete pivoting. The loss of stability
C     of the Gauss elimination with complete pivoting is rarely
C     encountered in practice.
C
C     The Lyapunov equation may be very ill-conditioned. In particular,
C     if DICO = 'D' and the pencil A - lambda * E has a pair of almost
C     reciprocal eigenvalues, or DICO = 'C' and the pencil has an almost
C     degenerate pair of eigenvalues, then the Lyapunov equation will be
C     ill-conditioned. Perturbed values were used to solve the equation.
C     Ill-conditioning can be detected by a very small value of the
C     reciprocal condition number RCOND.
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
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         ( ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOB, TRANS, UPLO
      DOUBLE PRECISION  FERR, SCALE, SEP
      INTEGER           INFO, LDA, LDE, LDQ, LDWORK, LDX, LDZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ALPHAI(*), ALPHAR(*), BETA(*),
     $                  DWORK(*), E(LDE,*), Q(LDQ,*), X(LDX,*),
     $                  Z(LDZ,*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      CHARACTER         ETRANS
      DOUBLE PRECISION  EST, EPS, NORMA, NORME, SCALE1
      INTEGER           I, INFO1, KASE, MINWRK, OPTWRK
      LOGICAL           ISDISC, ISFACT, ISTRAN, ISUPPR, WANTBH, WANTSP,
     $                  WANTX
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGGES, DLACON, MB01RD, MB01RW, SG03AX,
     $                  SG03AY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C     .. Executable Statements ..
C
C     Decode input parameters.
C
      ISDISC = LSAME( DICO,  'D' )
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      ISFACT = LSAME( FACT,  'F' )
      ISTRAN = LSAME( TRANS, 'T' )
      ISUPPR = LSAME( UPLO,  'U' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( ISDISC .OR. LSAME( DICO,  'C' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( WANTX .OR. WANTSP .OR. WANTBH ) ) THEN
         INFO = -2
      ELSEIF ( .NOT.( ISFACT .OR. LSAME( FACT,  'N' ) ) ) THEN
         INFO = -3
      ELSEIF ( .NOT.( ISTRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -4
      ELSEIF ( .NOT.( ISUPPR .OR. LSAME( UPLO,  'L' ) ) ) THEN
         INFO = -5
      ELSEIF ( N .LT. 0 ) THEN
         INFO = -6
      ELSEIF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -8
      ELSEIF ( LDE .LT. MAX( 1, N ) ) THEN
         INFO = -10
      ELSEIF ( LDQ .LT. MAX( 1, N ) ) THEN
         INFO = -12
      ELSEIF ( LDZ .LT. MAX( 1, N ) ) THEN
         INFO = -14
      ELSEIF ( LDX .LT. MAX( 1, N ) ) THEN
         INFO = -16
      ELSE
         INFO = 0
      END IF
      IF ( INFO .EQ. 0 ) THEN
C
C        Compute minimal workspace.
C
         IF ( WANTX ) THEN
            IF ( ISFACT ) THEN
               MINWRK = MAX( N, 1 )
            ELSE
               MINWRK = MAX( 8*N+16, 1 )
            END IF
         ELSE
            IF ( ISFACT ) THEN
               MINWRK = MAX( 2*N*N, 1 )
            ELSE
               MINWRK = MAX( 2*N*N, 8*N+16, 1 )
            END IF
         END IF
         IF ( MINWRK .GT. LDWORK ) THEN
            INFO = -25
         END IF
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SG03AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N .EQ. 0 ) THEN
         SCALE = ONE
         IF ( .NOT.WANTX ) SEP = ZERO
         IF ( WANTBH ) FERR = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( ISFACT ) THEN
C
C        Make sure the upper Hessenberg part of A is quasitriangular.
C
         DO 20 I = 1, N-2
            IF ( A(I+1,I).NE.ZERO .AND. A(I+2,I+1).NE.ZERO ) THEN
               INFO = 1
               RETURN
            END IF
   20    CONTINUE
      END IF
C
      IF ( .NOT.ISFACT ) THEN
C
C        Reduce A - lambda * E to generalized Schur form.
C
C           A := Q**T * A * Z   (upper quasitriangular)
C           E := Q**T * E * Z   (upper triangular)
C
C        ( Workspace: >= MAX(1,4*N) )
C
         CALL DGGES( 'Vectors', 'Vectors', 'N', 0, N, A, LDA, E, LDE,
     $               SDIM, ALPHAR,
     $               ALPHAI, BETA, Q, LDQ, Z, LDZ, DWORK, LDWORK,
     $               0, INFO1 )
         IF ( INFO1 .NE. 0 ) THEN
            INFO = 2
            RETURN
         END IF
         OPTWRK = INT( DWORK(1) )
      ELSE
         OPTWRK = MINWRK
      END IF
C
      IF ( WANTBH .OR. WANTX ) THEN
C
C        Transform right hand side.
C
C           X := Z**T * X * Z  or  X := Q**T * X * Q
C
C        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2.
C
C        ( Workspace: >= N )
C
         IF ( LDWORK .LT. N*N ) THEN
            IF ( ISTRAN ) THEN
               CALL MB01RW( UPLO, 'Transpose', N, N, X, LDX, Q, LDQ,
     $                      DWORK, INFO1 )
            ELSE
               CALL MB01RW( UPLO, 'Transpose', N, N, X, LDX, Z, LDZ,
     $                      DWORK, INFO1 )
            END IF
         ELSE
            IF ( ISTRAN ) THEN
               CALL MB01RD( UPLO, 'Transpose', N, N, ZERO, ONE, X, LDX,
     $                      Q, LDQ, X, LDX, DWORK, LDWORK, INFO )
            ELSE
               CALL MB01RD( UPLO, 'Transpose', N, N, ZERO, ONE, X, LDX,
     $                      Z, LDZ, X, LDX, DWORK, LDWORK, INFO )
            END IF
         END IF
         IF ( .NOT.ISUPPR ) THEN
            DO 40 I = 1, N-1
               CALL DCOPY( N-I, X(I+1,I), 1, X(I,I+1), LDX )
   40       CONTINUE
         END IF
         OPTWRK = MAX( OPTWRK, N*N )
C
C        Solve reduced generalized Lyapunov equation.
C
         IF ( ISDISC ) THEN
            CALL SG03AX( TRANS, N, A, LDA, E, LDE, X, LDX, SCALE, INFO1)
            IF ( INFO1 .NE. 0 )
     $         INFO = 3
         ELSE
            CALL SG03AY( TRANS, N, A, LDA, E, LDE, X, LDX, SCALE, INFO1)
            IF ( INFO1 .NE. 0 )
     $         INFO = 4
         END IF
C
C        Transform the solution matrix back.
C
C           X := Q * X * Q**T  or  X := Z * X * Z**T.
C
C        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2.
C
C        ( Workspace: >= N )
C
         IF ( LDWORK .LT. N*N ) THEN
            IF ( ISTRAN ) THEN
               CALL MB01RW( 'Upper', 'NoTranspose', N, N, X, LDX, Z,
     $                      LDZ, DWORK, INFO1 )
            ELSE
               CALL MB01RW( 'Upper', 'NoTranspose', N, N, X, LDX, Q,
     $                      LDQ, DWORK, INFO1 )
            END IF
         ELSE
            IF ( ISTRAN ) THEN
               CALL MB01RD( 'Upper', 'NoTranspose', N, N, ZERO, ONE, X,
     $                      LDX, Z, LDZ, X, LDX, DWORK, LDWORK, INFO )
            ELSE
               CALL MB01RD( 'Upper', 'NoTranspose', N, N, ZERO, ONE, X,
     $                      LDX, Q, LDQ, X, LDX, DWORK, LDWORK, INFO )
            END IF
         END IF
         DO 60 I = 1, N-1
            CALL DCOPY( N-I, X(I,I+1), LDX, X(I+1,I), 1 )
   60    CONTINUE
      END IF
C
      IF ( WANTBH .OR. WANTSP ) THEN
C
C        Estimate the 1-norm of the inverse Kronecker product matrix
C        belonging to the reduced generalized Lyapunov equation.
C
C        ( Workspace: 2*N*N )
C
         EST = ZERO
         KASE = 0
   80    CONTINUE
         CALL DLACON( N*N, DWORK(N*N+1), DWORK, IWORK, EST, KASE )
         IF ( KASE .NE. 0 ) THEN
            IF ( ( KASE.EQ.1 .AND. .NOT.ISTRAN ) .OR.
     $           ( KASE.NE.1 .AND. ISTRAN ) ) THEN
               ETRANS = 'N'
            ELSE
               ETRANS = 'T'
            END IF
            IF ( ISDISC ) THEN
               CALL SG03AX( ETRANS, N, A, LDA, E, LDE, DWORK, N, SCALE1,
     $                      INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 3
            ELSE
               CALL SG03AY( ETRANS, N, A, LDA, E, LDE, DWORK, N, SCALE1,
     $                      INFO1 )
               IF ( INFO1 .NE. 0 )
     $            INFO = 4
            END IF
         GOTO 80
         END IF
         SEP = SCALE1/EST
      END IF
C
C     Estimate the relative forward error.
C
C     ( Workspace: 2*N )
C
      IF ( WANTBH ) THEN
         EPS = DLAMCH( 'Precision' )
         DO 100 I = 1, N
            DWORK(I) = DNRM2( MIN( I+1, N ), A(1,I), 1 )
            DWORK(N+I) = DNRM2( I, E(1,I), 1 )
  100    CONTINUE
         NORMA = DNRM2( N, DWORK, 1 )
         NORME = DNRM2( N, DWORK(N+1), 1 )
         IF ( ISDISC ) THEN
            FERR = ( NORMA**2 + NORME**2 )*EPS/SEP
         ELSE
            FERR = TWO*NORMA*NORME*EPS/SEP
         END IF
      END IF
C
      DWORK(1) = DBLE( MAX( OPTWRK, MINWRK ) )
      RETURN
C *** Last line of SG03AD ***
      END
