      SUBROUTINE SB03QD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A,
     $                   LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEP,
     $                   RCOND, FERR, IWORK, DWORK, LDWORK, INFO )
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
C     To estimate the conditioning and compute an error bound on the
C     solution of the real continuous-time Lyapunov matrix equation
C
C         op(A)'*X + X*op(A) = scale*C
C
C     where op(A) = A or A' (A**T) and C is symmetric (C = C**T). The
C     matrix A is N-by-N, the right hand side C and the solution X are
C     N-by-N symmetric matrices, and scale is a given scale factor.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'C':  Compute the reciprocal condition number only;
C             = 'E':  Compute the error bound only;
C             = 'B':  Compute both the reciprocal condition number and
C                     the error bound.
C
C     FACT    CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix A is supplied on entry, as follows:
C             = 'F':  On entry, T and U (if LYAPUN = 'O') contain the
C                     factors from the real Schur factorization of the
C                     matrix A;
C             = 'N':  The Schur factorization of A will be computed
C                     and the factors will be stored in T and U (if
C                     LYAPUN = 'O').
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     UPLO    CHARACTER*1
C             Specifies which part of the symmetric matrix C is to be
C             used, as follows:
C             = 'U':  Upper triangular part;
C             = 'L':  Lower triangular part.
C
C     LYAPUN  CHARACTER*1
C             Specifies whether or not the original Lyapunov equations
C             should be solved in the iterative estimation process,
C             as follows:
C             = 'O':  Solve the original Lyapunov equations, updating
C                     the right-hand sides and solutions with the
C                     matrix U, e.g., X <-- U'*X*U;
C             = 'R':  Solve reduced Lyapunov equations only, without
C                     updating the right-hand sides and solutions.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X and C.  N >= 0.
C
C     SCALE   (input) DOUBLE PRECISION
C             The scale factor, scale, set by a Lyapunov solver.
C             0 <= SCALE <= 1.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             If FACT = 'N' or LYAPUN = 'O', the leading N-by-N part of
C             this array must contain the original matrix A.
C             If FACT = 'F' and LYAPUN = 'R', A is not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,N), if FACT = 'N' or  LYAPUN = 'O';
C             LDA >= 1,        if FACT = 'F' and LYAPUN = 'R'.
C
C     T       (input/output) DOUBLE PRECISION array, dimension
C             (LDT,N)
C             If FACT = 'F', then on entry the leading N-by-N upper
C             Hessenberg part of this array must contain the upper
C             quasi-triangular matrix T in Schur canonical form from a
C             Schur factorization of A.
C             If FACT = 'N', then this array need not be set on input.
C             On exit, (if INFO = 0 or INFO = N+1, for FACT = 'N') the
C             leading N-by-N upper Hessenberg part of this array
C             contains the upper quasi-triangular matrix T in Schur
C             canonical form from a Schur factorization of A.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= MAX(1,N).
C
C     U       (input or output) DOUBLE PRECISION array, dimension
C             (LDU,N)
C             If LYAPUN = 'O' and FACT = 'F', then U is an input
C             argument and on entry, the leading N-by-N part of this
C             array must contain the orthogonal matrix U from a real
C             Schur factorization of A.
C             If LYAPUN = 'O' and FACT = 'N', then U is an output
C             argument and on exit, if INFO = 0 or INFO = N+1, it
C             contains the orthogonal N-by-N matrix from a real Schur
C             factorization of A.
C             If LYAPUN = 'R', the array U is not referenced.
C
C     LDU     INTEGER
C             The leading dimension of the array U.
C             LDU >= 1,        if LYAPUN = 'R';
C             LDU >= MAX(1,N), if LYAPUN = 'O'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             If UPLO = 'U', the leading N-by-N upper triangular part of
C             this array must contain the upper triangular part of the
C             matrix C of the original Lyapunov equation (with
C             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov
C             equation (with matrix T), if LYAPUN = 'R'.
C             If UPLO = 'L', the leading N-by-N lower triangular part of
C             this array must contain the lower triangular part of the
C             matrix C of the original Lyapunov equation (with
C             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov
C             equation (with matrix T), if LYAPUN = 'R'.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,N).
C
C     X       (input) DOUBLE PRECISION array, dimension (LDX,N)
C             The leading N-by-N part of this array must contain the
C             symmetric solution matrix X of the original Lyapunov
C             equation (with matrix A), if LYAPUN = 'O', or of the
C             reduced Lyapunov equation (with matrix T), if
C             LYAPUN = 'R'.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= MAX(1,N).
C
C     SEP     (output) DOUBLE PRECISION
C             If JOB = 'C' or JOB = 'B', the estimated quantity
C             sep(op(A),-op(A)').
C             If N = 0, or X = 0, or JOB = 'E', SEP is not referenced.
C
C     RCOND   (output) DOUBLE PRECISION
C             If JOB = 'C' or JOB = 'B', an estimate of the reciprocal
C             condition number of the continuous-time Lyapunov equation.
C             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively.
C             If JOB = 'E', RCOND is not referenced.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'E' or JOB = 'B', an estimated forward error
C             bound for the solution X. If XTRUE is the true solution,
C             FERR bounds the magnitude of the largest entry in
C             (X - XTRUE) divided by the magnitude of the largest entry
C             in X.
C             If N = 0 or X = 0, FERR is set to 0.
C             If JOB = 'C', FERR is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N*N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the
C             optimal value of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             If JOB = 'C', then
C             LDWORK >= MAX(1,2*N*N),         if FACT = 'F';
C             LDWORK >= MAX(1,2*N*N,5*N),     if FACT = 'N'.
C             If JOB = 'E', or JOB = 'B', and LYAPUN  = 'O', then
C             LDWORK >= MAX(1,3*N*N),         if FACT = 'F';
C             LDWORK >= MAX(1,3*N*N,5*N),     if FACT = 'N'.
C             If JOB = 'E', or JOB = 'B', and LYAPUN  = 'R', then
C             LDWORK >= MAX(1,3*N*N+N-1),     if FACT = 'F';
C             LDWORK >= MAX(1,3*N*N+N-1,5*N), if FACT = 'N'.
C             For optimum performance LDWORK should sometimes be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, i <= N, the QR algorithm failed to
C                   complete the reduction to Schur canonical form (see
C                   LAPACK Library routine DGEES); on exit, the matrix
C                   T(i+1:N,i+1:N) contains the partially converged
C                   Schur form, and DWORK(i+1:N) and DWORK(N+i+1:2*N)
C                   contain the real and imaginary parts, respectively,
C                   of the converged eigenvalues; this error is unlikely
C                   to appear;
C             = N+1:  if the matrices T and -T' have common or very
C                   close eigenvalues; perturbed values were used to
C                   solve Lyapunov equations, but the matrix T, if given
C                   (for FACT = 'F'), is unchanged.
C
C     METHOD
C
C     The condition number of the continuous-time Lyapunov equation is
C     estimated as
C
C     cond = (norm(Theta)*norm(A) + norm(inv(Omega))*norm(C))/norm(X),
C
C     where Omega and Theta are linear operators defined by
C
C     Omega(W) = op(A)'*W + W*op(A),
C     Theta(W) = inv(Omega(op(W)'*X + X*op(W))).
C
C     The routine estimates the quantities
C
C     sep(op(A),-op(A)') = 1 / norm(inv(Omega))
C
C     and norm(Theta) using 1-norm condition estimators.
C
C     The forward error bound is estimated using a practical error bound
C     similar to the one proposed in [1].
C
C     REFERENCES
C
C     [1] Higham, N.J.
C         Perturbation theory and backward error for AX-XB=C.
C         BIT, vol. 33, pp. 124-136, 1993.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C     The accuracy of the estimates obtained depends on the solution
C     accuracy and on the properties of the 1-norm estimator.
C
C     FURTHER COMMENTS
C
C     The option LYAPUN = 'R' may occasionally produce slightly worse
C     or better estimates, and it is much faster than the option 'O'.
C     When SEP is computed and it is zero, the routine returns
C     immediately, with RCOND and FERR (if requested) set to 0 and 1,
C     respectively. In this case, the equation is singular.
C
C     CONTRIBUTORS
C
C     P. Petkov, Tech. University of Sofia, December 1998.
C     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, March 2003.
C
C     KEYWORDS
C
C     Lyapunov equation, orthogonal transformation, real Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D0,
     $                     THREE = 3.0D0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEP
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                   T( LDT, * ), U( LDU, * ), X( LDX, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            JOBB, JOBC, JOBE, LOWER, NOFACT, NOTRNA,
     $                   UPDATE
      CHARACTER          SJOB, TRANAT
      INTEGER            I, IABS, IRES, IWRK, IXBS, J, JJ, JX, LDW, NN,
     $                   SDIM, WRKOPT
      DOUBLE PRECISION   ANORM, CNORM, DENOM, EPS, EPSN, TEMP, THNORM,
     $                   TMAX, XANORM, XNORM
C     ..
C     .. Local Arrays ..
      LOGICAL            BWORK( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, SELECT
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANHS, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANHS, DLANSY, LSAME, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEES, DLACPY, DLASET, DSYR2K, MB01UD,
     $                   MB01UW, SB03QX, SB03QY, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      JOBC   = LSAME( JOB,    'C' )
      JOBE   = LSAME( JOB,    'E' )
      JOBB   = LSAME( JOB,    'B' )
      NOFACT = LSAME( FACT,   'N' )
      NOTRNA = LSAME( TRANA,  'N' )
      LOWER  = LSAME( UPLO,   'L' )
      UPDATE = LSAME( LYAPUN, 'O' )
C
      NN = N*N
      IF( JOBC ) THEN
         LDW = 2*NN
      ELSE
         LDW = 3*NN
      END IF
      IF( .NOT.( JOBC .OR. UPDATE ) )
     $   LDW = LDW + N - 1
C
      INFO = 0
      IF( .NOT.( JOBB .OR. JOBC .OR. JOBE ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( NOFACT .OR. LSAME( FACT, 'F' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( NOTRNA .OR. LSAME( TRANA,  'T' ) .OR.
     $                            LSAME( TRANA,  'C' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( LOWER  .OR. LSAME( UPLO,   'U' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT.( UPDATE .OR. LSAME( LYAPUN, 'R' ) ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( SCALE.LT.ZERO .OR. SCALE.GT.ONE ) THEN
         INFO = -7
      ELSE IF( LDA.LT.1 .OR.
     $       ( LDA.LT.N .AND. ( UPDATE .OR. NOFACT ) ) ) THEN
         INFO = -9
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDU.LT.1 .OR. ( LDU.LT.N .AND. UPDATE ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.1 .OR.
     $       ( LDWORK.LT.LDW .AND. .NOT.NOFACT ) .OR.
     $       ( LDWORK.LT.MAX( LDW, 5*N ) .AND. NOFACT ) ) THEN
         INFO = -23
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03QD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         IF( .NOT.JOBE )
     $      RCOND = ONE
         IF( .NOT.JOBC )
     $      FERR  = ZERO
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Compute the 1-norm of the matrix X.
C
      XNORM = DLANSY( '1-norm', UPLO, N, X, LDX, DWORK )
      IF( XNORM.EQ.ZERO ) THEN
C
C        The solution is zero.
C
         IF( .NOT.JOBE )
     $      RCOND = ZERO
         IF( .NOT.JOBC )
     $      FERR  = ZERO
         DWORK( 1 ) = DBLE( N )
         RETURN
      END IF
C
C     Compute the 1-norm of A or T.
C
      IF( NOFACT .OR. UPDATE ) THEN
         ANORM  = DLANGE( '1-norm', N, N, A, LDA, DWORK )
      ELSE
         ANORM  = DLANHS( '1-norm', N, T, LDT, DWORK )
      END IF
C
C     For the special case A = 0, set SEP and RCOND to 0.
C     For the special case A = I, set SEP to 2 and RCOND to 1.
C     A quick test is used in general.
C
      IF( ANORM.EQ.ONE ) THEN
         IF( NOFACT .OR. UPDATE ) THEN
            CALL DLACPY( 'Full', N, N, A, LDA, DWORK, N )
         ELSE
            CALL DLACPY( 'Full', N, N, T, LDT, DWORK, N )
            IF( N.GT.2 )
     $         CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, DWORK( 3 ),
     $                      N )
         END IF
         DWORK( NN+1 ) = ONE
         CALL DAXPY( N, -ONE, DWORK( NN+1 ), 0, DWORK, N+1 )
         IF( DLANGE( 'Max', N, N, DWORK, N, DWORK ).EQ.ZERO ) THEN
            IF( .NOT.JOBE ) THEN
               SEP   = TWO
               RCOND = ONE
            END IF
            IF( JOBC ) THEN
               DWORK( 1 ) = DBLE( NN + 1 )
               RETURN
            ELSE
C
C              Set FERR for the special case A = I.
C
               CALL DLACPY( UPLO, N, N, X, LDX, DWORK, N )
C
               IF( LOWER ) THEN
                  DO 10 J = 1, N
                     CALL DAXPY( N-J+1, -SCALE/TWO, C( J, J ), 1,
     $                           DWORK( (J-1)*N+J ), 1 )
   10             CONTINUE
               ELSE
                  DO 20 J = 1, N
                     CALL DAXPY( J, -SCALE/TWO, C( 1, J ), 1,
     $                           DWORK( (J-1)*N+1 ), 1 )
   20             CONTINUE
               END IF
C
               FERR = MIN( ONE, DLANSY( '1-norm', UPLO, N, DWORK, N,
     $                                  DWORK( NN+1 ) ) / XNORM )
               DWORK( 1 ) = DBLE( NN + N )
               RETURN
            END IF
         END IF
C
      ELSE IF( ANORM.EQ.ZERO ) THEN
         IF( .NOT.JOBE ) THEN
            SEP   = ZERO
            RCOND = ZERO
         END IF
         IF( .NOT.JOBC )
     $      FERR = ONE
         DWORK( 1 ) = DBLE( N )
         RETURN
      END IF
C
C     General case.
C
      CNORM = DLANSY( '1-norm', UPLO, N, C, LDC, DWORK )
C
C     Workspace usage.
C
      IABS = 0
      IXBS = IABS + NN
      IRES = IXBS + NN
      IWRK = IRES + NN
      WRKOPT = 0
C
      IF( NOFACT ) THEN
C
C        Compute the Schur factorization of A, A = U*T*U'.
C        Workspace:  need   5*N;
C                    prefer larger.
C        (Note: Comments in the code beginning "Workspace:" describe the
C        minimal amount of real workspace needed at that point in the
C        code, as well as the preferred amount for good performance.)
C
         CALL DLACPY( 'Full', N, N, A, LDA, T, LDT )
         IF( UPDATE ) THEN
            SJOB = 'V'
         ELSE
            SJOB = 'N'
         END IF
         CALL DGEES( SJOB, 'Not ordered', SELECT, N, T, LDT, SDIM,
     $               DWORK( 1 ), DWORK( N+1 ), U, LDU, DWORK( 2*N+1 ),
     $               LDWORK-2*N, BWORK, INFO )
         IF( INFO.GT.0 )
     $      RETURN
         WRKOPT = INT( DWORK( 2*N+1 ) ) + 2*N
      END IF
C
      IF( .NOT.JOBE ) THEN
C
C        Estimate sep(op(A),-op(A)') = sep(op(T),-op(T)') and
C        norm(Theta).
C        Workspace 2*N*N.
C
         CALL SB03QY( 'Both', TRANA, LYAPUN, N, T, LDT, U, LDU, X, LDX,
     $                SEP, THNORM, IWORK, DWORK, LDWORK, INFO )
C
         WRKOPT = MAX( WRKOPT, 2*NN )
C
C        Return if the equation is singular.
C
         IF( SEP.EQ.ZERO ) THEN
            RCOND = ZERO
            IF( JOBB )
     $         FERR  = ONE
            DWORK( 1 ) = DBLE( WRKOPT )
            RETURN
         END IF
C
C        Estimate the reciprocal condition number.
C
         TMAX = MAX( SEP, XNORM, ANORM )
         IF( TMAX.LE.ONE ) THEN
            TEMP  =     SEP*XNORM
            DENOM = ( SCALE*CNORM ) + ( SEP*ANORM )*THNORM
         ELSE
            TEMP  =   (   SEP / TMAX )*( XNORM / TMAX )
            DENOM = ( ( SCALE / TMAX )*( CNORM / TMAX ) ) +
     $              ( (   SEP / TMAX )*( ANORM / TMAX ) )*THNORM
         END IF
         IF( TEMP.GE.DENOM ) THEN
            RCOND = ONE
         ELSE
            RCOND = TEMP / DENOM
         END IF
      END IF
C
      IF( .NOT.JOBC ) THEN
C
C        Form a triangle of the residual matrix
C        R = op(A)'*X + X*op(A) - scale*C, or
C        R = op(T)'*X + X*op(T) - scale*C,
C        exploiting the symmetry.
C        Workspace 3*N*N.
C
         IF( NOTRNA ) THEN
            TRANAT = 'T'
         ELSE
            TRANAT = 'N'
         END IF
C
         IF( UPDATE ) THEN
C
            CALL DLACPY( UPLO, N, N, C, LDC, DWORK( IRES+1 ), N )
            CALL DSYR2K( UPLO, TRANAT, N, N, ONE, A, LDA, X, LDX,
     $                   -SCALE, DWORK( IRES+1 ), N )
         ELSE
            CALL MB01UD( 'Right', TRANA, N, N, ONE, T, LDT, X, LDX,
     $                   DWORK( IRES+1 ), N, INFO )
            JJ = IRES + 1
            IF( LOWER ) THEN
               DO 30 J = 1, N
                  CALL DAXPY( N-J+1, ONE, DWORK( JJ ), N, DWORK( JJ ),
     $                        1 )
                  CALL DAXPY( N-J+1, -SCALE, C( J, J ), 1, DWORK( JJ ),
     $                        1 )
                  JJ = JJ + N + 1
   30          CONTINUE
            ELSE
               DO 40 J = 1, N
                  CALL DAXPY( J, ONE, DWORK( IRES+J ), N, DWORK( JJ ),
     $                        1 )
                  CALL DAXPY( J, -SCALE, C( 1, J ), 1, DWORK( JJ ), 1 )
                  JJ = JJ + N
   40          CONTINUE
            END IF
         END IF
C
         WRKOPT = MAX( WRKOPT, 3*NN )
C
C        Get the machine precision.
C
         EPS  = DLAMCH( 'Epsilon' )
         EPSN = EPS*DBLE( N + 3 )
         TEMP = EPS*THREE*SCALE
C
C        Add to abs(R) a term that takes account of rounding errors in
C        forming R:
C          abs(R) := abs(R) + EPS*(3*scale*abs(C) +
C                    (n+3)*(abs(op(A))'*abs(X) + abs(X)*abs(op(A)))), or
C          abs(R) := abs(R) + EPS*(3*scale*abs(C) +
C                    (n+3)*(abs(op(T))'*abs(X) + abs(X)*abs(op(T)))),
C        where EPS is the machine precision.
C
         DO 60 J = 1, N
            DO 50 I = 1, N
               DWORK( IXBS+(J-1)*N+I ) = ABS( X( I, J ) )
   50       CONTINUE
   60    CONTINUE
C
         IF( LOWER ) THEN
            DO 80 J = 1, N
               DO 70 I = J, N
                  DWORK( IRES+(J-1)*N+I ) = TEMP*ABS( C( I, J ) ) +
     $                   ABS( DWORK( IRES+(J-1)*N+I ) )
   70          CONTINUE
   80       CONTINUE
         ELSE
            DO 100 J = 1, N
               DO 90 I = 1, J
                  DWORK( IRES+(J-1)*N+I ) = TEMP*ABS( C( I, J ) ) +
     $                   ABS( DWORK( IRES+(J-1)*N+I ) )
   90          CONTINUE
  100       CONTINUE
         END IF
C
         IF( UPDATE ) THEN
C
C           Workspace 3*N*N.
C
            DO 120 J = 1, N
               DO 110 I = 1, N
                  DWORK( IABS+(J-1)*N+I ) = ABS( A( I, J ) )
  110          CONTINUE
  120       CONTINUE
C
            CALL DSYR2K( UPLO, TRANAT, N, N, EPSN, DWORK( IABS+1 ), N,
     $                   DWORK( IXBS+1 ), N, ONE,  DWORK( IRES+1 ), N )
         ELSE
C
C           Workspace 3*N*N + N - 1.
C
            DO 140 J = 1, N
               DO 130 I = 1, MIN( J+1, N )
                  DWORK( IABS+(J-1)*N+I ) = ABS( T( I, J ) )
  130          CONTINUE
  140       CONTINUE
C
            CALL MB01UW( 'Left', TRANAT, N, N, EPSN, DWORK( IABS+1 ),
     $                   N, DWORK( IXBS+1), N, DWORK( IWRK+1 ),
     $                   LDWORK-IWRK, INFO )
            JJ = IRES + 1
            JX = IXBS + 1
            IF( LOWER ) THEN
               DO 150 J = 1, N
                  CALL DAXPY( N-J+1, ONE, DWORK( JX ), N, DWORK( JX ),
     $                        1 )
                  CALL DAXPY( N-J+1, ONE, DWORK( JX ), 1, DWORK( JJ ),
     $                        1 )
                  JJ = JJ + N + 1
                  JX = JX + N + 1
  150          CONTINUE
            ELSE
               DO 160 J = 1, N
                  CALL DAXPY( J, ONE, DWORK( IXBS+J ), N, DWORK( JX ),
     $                        1 )
                  CALL DAXPY( J, ONE, DWORK( JX ), 1, DWORK( JJ ), 1 )
                  JJ = JJ + N
                  JX = JX + N
  160          CONTINUE
            END IF
C
            WRKOPT = MAX( WRKOPT, 3*NN + N - 1 )
         END IF
C
C        Compute forward error bound, using matrix norm estimator.
C        Workspace 3*N*N.
C
         XANORM = DLANSY( 'Max', UPLO, N, X, LDX, DWORK )
C
         CALL SB03QX( TRANA, UPLO, LYAPUN, N, XANORM, T, LDT, U, LDU,
     $                DWORK( IRES+1 ), N, FERR, IWORK, DWORK, IRES,
     $                INFO )
      END IF
C
      DWORK( 1 ) = DBLE( WRKOPT )
      RETURN
C
C *** Last line of SB03QD ***
      END
