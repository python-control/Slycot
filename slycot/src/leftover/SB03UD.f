      SUBROUTINE SB03UD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A,
     $                   LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEPD,
     $                   RCOND, FERR, WR, WI, IWORK, DWORK, LDWORK,
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
C     To solve the real discrete-time Lyapunov matrix equation
C
C            op(A)'*X*op(A) - X = scale*C,
C
C     estimate the conditioning, and compute an error bound on the
C     solution X, where op(A) = A or A' (A**T), the matrix A is N-by-N,
C     the right hand side C and the solution X are N-by-N symmetric
C     matrices (C = C', X = X'), and scale is an output scale factor,
C     set less than or equal to 1 to avoid overflow in X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'X':  Compute the solution only;
C             = 'S':  Compute the separation only;
C             = 'C':  Compute the reciprocal condition number only;
C             = 'E':  Compute the error bound only;
C             = 'A':  Compute all: the solution, separation, reciprocal
C                     condition number, and the error bound.
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
C             Specifies whether or not the original or "reduced"
C             Lyapunov equations should be solved, as follows:
C             = 'O':  Solve the original Lyapunov equations, updating
C                     the right-hand sides and solutions with the
C                     matrix U, e.g., X <-- U'*X*U;
C             = 'R':  Solve reduced Lyapunov equations only, without
C                     updating the right-hand sides and solutions.
C                     This means that a real Schur form T of A appears
C                     in the equation, instead of A.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X, and C.  N >= 0.
C
C     SCALE   (input or output) DOUBLE PRECISION
C             If JOB = 'C' or JOB = 'E', SCALE is an input argument:
C             the scale factor, set by a Lyapunov solver.
C             0 <= SCALE <= 1.
C             If JOB = 'X' or JOB = 'A', SCALE is an output argument:
C             the scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C             If JOB = 'S', this argument is not used.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             If FACT = 'N' or (LYAPUN = 'O' and JOB <> 'X'), the
C             leading N-by-N part of this array must contain the
C             original matrix A.
C             If FACT = 'F' and (LYAPUN = 'R' or JOB = 'X'), A is
C             not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,N), if FACT = 'N' or LYAPUN = 'O' and
C                                               JOB <> 'X';
C             LDA >= 1,        otherwise.
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
C             The contents of array T is not modified if FACT = 'F'.
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
C             If JOB <> 'S' and UPLO = 'U', the leading N-by-N upper
C             triangular part of this array must contain the upper
C             triangular part of the matrix C of the original Lyapunov
C             equation (with matrix A), if LYAPUN = 'O', or of the
C             reduced Lyapunov equation (with matrix T), if
C             LYAPUN = 'R'.
C             If JOB <> 'S' and UPLO = 'L', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangular part of the matrix C of the original Lyapunov
C             equation (with matrix A), if LYAPUN = 'O', or of the
C             reduced Lyapunov equation (with matrix T), if
C             LYAPUN = 'R'.
C             The remaining strictly triangular part of this array is
C             used as workspace.
C             If JOB = 'X', then this array may be identified with X
C             in the call of this routine.
C             If JOB = 'S', the array C is not referenced.
C
C     LDC     INTEGER
C             The leading dimension of the array C.
C             LDC >= 1,        if JOB = 'S';
C             LDC >= MAX(1,N), otherwise.
C
C     X       (input or output) DOUBLE PRECISION array, dimension
C             (LDX,N)
C             If JOB = 'C' or 'E', then X is an input argument and on
C             entry, the leading N-by-N part of this array must contain
C             the symmetric solution matrix X of the original Lyapunov
C             equation (with matrix A), if LYAPUN = 'O', or of the
C             reduced Lyapunov equation (with matrix T), if
C             LYAPUN = 'R'.
C             If JOB = 'X' or 'A', then X is an output argument and on
C             exit, if INFO = 0 or INFO = N+1, the leading N-by-N part
C             of this array contains the symmetric solution matrix X of
C             of the original Lyapunov equation (with matrix A), if
C             LYAPUN = 'O', or of the reduced Lyapunov equation (with
C             matrix T), if LYAPUN = 'R'.
C             If JOB = 'S', the array X is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of the array X.
C             LDX >= 1,        if JOB = 'S';
C             LDX >= MAX(1,N), otherwise.
C
C     SEPD    (output) DOUBLE PRECISION
C             If JOB = 'S' or JOB = 'C' or JOB = 'A', and INFO = 0 or
C             INFO = N+1, SEPD contains the estimated separation of the
C             matrices op(A) and op(A)', sepd(op(A),op(A)').
C             If N = 0, or X = 0, or JOB = 'X' or JOB = 'E', SEPD is not
C             referenced.
C
C     RCOND   (output) DOUBLE PRECISION
C             If JOB = 'C' or JOB = 'A', an estimate of the reciprocal
C             condition number of the continuous-time Lyapunov equation.
C             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively.
C             If JOB = 'X' or JOB = 'S' or JOB = 'E', RCOND is not
C             referenced.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'E' or JOB = 'A', and INFO = 0 or INFO = N+1,
C             FERR contains an estimated forward error bound for the
C             solution X. If XTRUE is the true solution, FERR bounds the
C             relative error in the computed solution, measured in the
C             Frobenius norm:  norm(X - XTRUE)/norm(XTRUE).
C             If N = 0 or X = 0, FERR is set to 0.
C             If JOB = 'X' or JOB = 'S' or JOB = 'C', FERR is not
C             referenced.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             If FACT = 'N', and INFO = 0 or INFO = N+1, WR and WI
C             contain the real and imaginary parts, respectively, of the
C             eigenvalues of A.
C             If FACT = 'F', WR and WI are not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N*N)
C             This array is not referenced if JOB = 'X'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the
C             optimal value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If JOB = 'X', then
C             LDWORK >= MAX(1,N*N,2*N),       if FACT = 'F';
C             LDWORK >= MAX(1,N*N,3*N),       if FACT = 'N'.
C             If JOB = 'S', then
C             LDWORK >= MAX(3,2*N*N).
C             If JOB = 'C', then
C             LDWORK >= MAX(3,2*N*N) + N*N.
C             If JOB = 'E', or JOB = 'A', then
C             LDWORK >= MAX(3,2*N*N) + N*N + 2*N.
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
C                   Schur form, and the elements i+1:n of WR and WI
C                   contain the real and imaginary parts, respectively,
C                   of the converged eigenvalues; this error is unlikely
C                   to appear;
C             = N+1:  if the matrix T has almost reciprocal eigenvalues;
C                   perturbed values were used to solve Lyapunov
C                   equations, but the matrix T, if given (for
C                   FACT = 'F'), is unchanged.
C
C     METHOD
C
C     After reducing matrix A to real Schur canonical form (if needed),
C     a discrete-time version of the Bartels-Stewart algorithm is used.
C     A set of equivalent linear algebraic systems of equations of order
C     at most four are formed and solved using Gaussian elimination with
C     complete pivoting.
C
C     The condition number of the discrete-time Lyapunov equation is
C     estimated as
C
C     cond = (norm(Theta)*norm(A) + norm(inv(Omega))*norm(C))/norm(X),
C
C     where Omega and Theta are linear operators defined by
C
C     Omega(W) = op(A)'*W*op(A) - W,
C     Theta(W) = inv(Omega(op(W)'*X*op(A) + op(A)'*X*op(W))).
C
C     The routine estimates the quantities
C
C     sepd(op(A),op(A)') = 1 / norm(inv(Omega))
C
C     and norm(Theta) using 1-norm condition estimators.
C
C     The forward error bound is estimated using a practical error bound
C     similar to the one proposed in [3].
C
C     REFERENCES
C
C     [1] Barraud, A.Y.                   T
C         A numerical algorithm to solve A XA - X = Q.
C         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977.
C
C     [2] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [3] Higham, N.J.
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
C     The "separation" sepd of op(A) and op(A)' can also be defined as
C
C            sepd( op(A), op(A)' ) = sigma_min( T ),
C
C     where sigma_min(T) is the smallest singular value of the
C     N*N-by-N*N matrix
C
C        T = kprod( op(A)', op(A)' ) - I(N**2).
C
C     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the
C     Kronecker product. The routine estimates sigma_min(T) by the
C     reciprocal of an estimate of the 1-norm of inverse(T). The true
C     reciprocal 1-norm of inverse(T) cannot differ from sigma_min(T) by
C     more than a factor of N.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999.
C     This is an extended and improved version of Release 3.0 routine
C     SB03PD.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004.
C
C     KEYWORDS
C
C     Lyapunov equation, orthogonal transformation, real Schur form,
C     Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEPD
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                   T( LDT, * ), U( LDU, * ), WI( * ), WR( * ),
     $                   X( LDX, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            JOBA, JOBC, JOBE, JOBS, JOBX, LOWER, NOFACT,
     $                   NOTRNA, UPDATE
      CHARACTER          CFACT, JOBL, SJOB
      INTEGER            LDW, NN, SDIM
      DOUBLE PRECISION   THNORM
C     ..
C     .. Local Arrays ..
      LOGICAL            BWORK( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, SELECT
      EXTERNAL           LSAME, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEES, DLACPY, DSCAL, MA02ED, MB01RU, SB03MX,
     $                   SB03SD, SB03SY, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX
C     ..
C     .. Executable Statements ..
C
C     Decode option parameters.
C
      JOBX   = LSAME( JOB,    'X' )
      JOBS   = LSAME( JOB,    'S' )
      JOBC   = LSAME( JOB,    'C' )
      JOBE   = LSAME( JOB,    'E' )
      JOBA   = LSAME( JOB,    'A' )
      NOFACT = LSAME( FACT,   'N' )
      NOTRNA = LSAME( TRANA,  'N' )
      LOWER  = LSAME( UPLO,   'L' )
      UPDATE = LSAME( LYAPUN, 'O' )
C
C     Compute workspace.
C
      NN = N*N
      IF( JOBX ) THEN
         IF( NOFACT ) THEN
            LDW = MAX( 1, NN, 3*N )
         ELSE
            LDW = MAX( 1, NN, 2*N )
         END IF
      ELSE IF( JOBS ) THEN
         LDW = MAX( 3, 2*NN )
      ELSE IF( JOBC ) THEN
         LDW = MAX( 3, 2*NN ) + NN
      ELSE
         LDW = MAX( 3, 2*NN ) + NN + 2*N
      END IF
C
C     Test the scalar input parameters.
C
      INFO = 0
      IF( .NOT.( JOBX .OR. JOBS .OR. JOBC .OR. JOBE .OR. JOBA ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( NOFACT .OR. LSAME( FACT,   'F' ) ) ) THEN
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
      ELSE IF( ( JOBC .OR. JOBE ) .AND.
     $         ( SCALE.LT.ZERO .OR. SCALE.GT.ONE ) )THEN
         INFO = -7
      ELSE IF( LDA.LT.1 .OR.
     $       ( LDA.LT.N .AND. ( ( UPDATE .AND. .NOT.JOBX ) .OR.
     $                            NOFACT ) ) ) THEN
         INFO = -9
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDU.LT.1 .OR. ( LDU.LT.N .AND. UPDATE ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.1 .OR. ( .NOT.JOBS .AND. LDC.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDX.LT.1 .OR. ( .NOT.JOBS .AND. LDX.LT.N ) ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.LDW ) THEN
         INFO = -25
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03UD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         IF( JOBX .OR. JOBA )
     $      SCALE = ONE
         IF( JOBC .OR. JOBA )
     $      RCOND = ONE
         IF( JOBE .OR. JOBA )
     $      FERR  = ZERO
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
      IF( NOFACT ) THEN
C
C        Compute the Schur factorization of A.
C        Workspace:  need   3*N;
C                    prefer larger.
C
         CALL DLACPY( 'Full', N, N, A, LDA, T, LDT )
         IF( UPDATE ) THEN
            SJOB = 'V'
         ELSE
            SJOB = 'N'
         END IF
         CALL DGEES( SJOB, 'Not ordered', SELECT, N, T, LDT, SDIM, WR,
     $               WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         IF( INFO.GT.0 )
     $      RETURN
         LDW = MAX( LDW, INT( DWORK( 1 ) ) )
         CFACT = 'F'
      ELSE
         CFACT = FACT
      END IF
C
      IF( JOBX .OR. JOBA ) THEN
C
C        Copy the right-hand side in X.
C
         CALL DLACPY( UPLO, N, N, C, LDC, X, LDX )
C
         IF( UPDATE ) THEN
C
C           Transform the right-hand side.
C           Workspace:  need   N*N.
C
            CALL MB01RU( UPLO, 'Transpose', N, N, ZERO, ONE, X, LDX, U,
     $                   LDU, X, LDX, DWORK, LDWORK, INFO )
            CALL DSCAL( N, HALF, X, LDX+1 )
         END IF
C
C        Fill in the remaining triangle of X.
C
         CALL MA02ED( UPLO, N, X, LDX )
C
C        Solve the transformed equation.
C        Workspace:  2*N.
C
         CALL SB03MX( TRANA, N, T, LDT, X, LDX, SCALE, DWORK, INFO )
         IF( INFO.GT.0 )
     $      INFO = N + 1
C
         IF( UPDATE ) THEN
C
C           Transform back the solution.
C
            CALL MB01RU( UPLO, 'No transpose', N, N, ZERO, ONE, X, LDX,
     $                   U, LDU, X, LDX, DWORK, LDWORK, INFO )
            CALL DSCAL( N, HALF, X, LDX+1 )
C
C           Fill in the remaining triangle of X.
C
            CALL MA02ED( UPLO, N, X, LDX )
         END IF
      END IF
C
      IF( JOBS ) THEN
C
C        Estimate sepd(op(A),op(A)').
C        Workspace:  MAX(3,2*N*N).
C
         CALL SB03SY( 'Separation', TRANA, LYAPUN, N, T, LDT, U, LDU,
     $                DWORK, 1, SEPD, THNORM, IWORK, DWORK, LDWORK,
     $                INFO )
C
      ELSE IF( .NOT.JOBX ) THEN
C
C        Estimate the reciprocal condition and/or the error bound.
C        Workspace:  MAX(3,2*N*N) + N*N + a*N, where:
C                    a = 2, if JOB = 'E' or JOB = 'A';
C                    a = 0, otherwise.
C
         IF( JOBA ) THEN
            JOBL = 'B'
         ELSE
            JOBL = JOB
         END IF
         CALL SB03SD( JOBL, CFACT, TRANA, UPLO, LYAPUN, N, SCALE, A,
     $                LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEPD, RCOND,
     $                FERR, IWORK, DWORK, LDWORK, INFO )
         LDW = MAX( LDW, INT( DWORK( 1 ) ) )
      END IF
C
      DWORK( 1 ) = DBLE( LDW )
C
      RETURN
C *** Last line of SB03UD ***
      END
