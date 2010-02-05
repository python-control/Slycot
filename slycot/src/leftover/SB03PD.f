      SUBROUTINE SB03PD( JOB, FACT, TRANA, N, A, LDA, U, LDU, C, LDC,
     $                   SCALE, SEPD, FERR, WR, WI, IWORK, DWORK,
     $                   LDWORK, INFO )
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
C     To solve the real discrete Lyapunov matrix equation
C
C            op(A)'*X*op(A) - X = scale*C
C
C     and/or estimate the quantity, called separation,
C
C         sepd(op(A),op(A)') = min norm(op(A)'*X*op(A) - X)/norm(X)
C
C     where op(A) = A or A' (A**T) and C is symmetric (C = C').
C     (A' denotes the transpose of the matrix A.) A is N-by-N, the right
C     hand side C and the solution X are N-by-N, and scale is an output
C     scale factor, set less than or equal to 1 to avoid overflow in X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'X':  Compute the solution only;
C             = 'S':  Compute the separation only;
C             = 'B':  Compute both the solution and the separation.
C
C     FACT    CHARACTER*1
C             Specifies whether or not the real Schur factorization
C             of the matrix A is supplied on entry, as follows:
C             = 'F':  On entry, A and U contain the factors from the
C                     real Schur factorization of the matrix A;
C             = 'N':  The Schur factorization of A will be computed
C                     and the factors will be stored in A and U.
C
C     TRANA   CHARACTER*1
C             Specifies the form of op(A) to be used, as follows:
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, X, and C.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A. If FACT = 'F', then A contains
C             an upper quasi-triangular matrix in Schur canonical form.
C             On exit, if INFO = 0 or INFO = N+1, the leading N-by-N
C             part of this array contains the upper quasi-triangular
C             matrix in Schur canonical form from the Shur factorization
C             of A. The contents of array A is not modified if
C             FACT = 'F'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     U       (input or output) DOUBLE PRECISION array, dimension
C             (LDU,N)
C             If FACT = 'F', then U is an input argument and on entry
C             it must contain the orthogonal matrix U from the real
C             Schur factorization of A.
C             If FACT = 'N', then U is an output argument and on exit,
C             if INFO = 0 or INFO = N+1, it contains the orthogonal
C             N-by-N matrix from the real Schur factorization of A.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry with JOB = 'X' or 'B', the leading N-by-N part of
C             this array must contain the symmetric matrix C.
C             On exit with JOB = 'X' or 'B', if INFO = 0 or INFO = N+1,
C             the leading N-by-N part of C has been overwritten by the
C             symmetric solution matrix X.
C             If JOB = 'S', C is not referenced.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= 1,        if JOB = 'S';
C             LDC >= MAX(1,N), otherwise.
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, scale, set less than or equal to 1 to
C             prevent the solution overflowing.
C
C     SEPD    (output) DOUBLE PRECISION
C             If JOB = 'S' or JOB = 'B', and INFO = 0 or INFO = N+1,
C             SEPD contains the estimate in the 1-norm of
C             sepd(op(A),op(A)').
C             If JOB = 'X' or N = 0, SEPD is not referenced.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'B', and INFO = 0 or INFO = N+1, FERR contains
C             an estimated forward error bound for the solution X.
C             If XTRUE is the true solution, FERR bounds the relative
C             error in the computed solution, measured in the Frobenius
C             norm:  norm(X - XTRUE)/norm(XTRUE).
C             If JOB = 'X' or JOB = 'S', FERR is not referenced.
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
C             The length of the array DWORK.  LDWORK >= 1 and
C             If JOB = 'X' then
C                If FACT = 'F', LDWORK >= MAX(N*N,2*N);
C                If FACT = 'N', LDWORK >= MAX(N*N,3*N).
C             If JOB = 'S' or JOB = 'B' then
C                LDWORK >= 2*N*N + 2*N.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, the QR algorithm failed to compute all
C                   the eigenvalues (see LAPACK Library routine DGEES);
C                   elements i+1:n of WR and WI contain eigenvalues
C                   which have converged, and A contains the partially
C                   converged Schur form;
C             = N+1:  if matrix A has almost reciprocal eigenvalues;
C                   perturbed values were used to solve the equation
C                   (but the matrix A is unchanged).
C
C     METHOD
C
C     After reducing matrix A to real Schur canonical form (if needed),
C     a discrete-time version of the Bartels-Stewart algorithm is used.
C     A set of equivalent linear algebraic systems of equations of order
C     at most four are formed and solved using Gaussian elimination with
C     complete pivoting.
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
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     FURTHER COMMENTS
C
C     SEPD is defined as
C
C            sepd( op(A), op(A)' ) = sigma_min( T )
C
C     where sigma_min(T) is the smallest singular value of the
C     N*N-by-N*N matrix
C
C        T = kprod( op(A)', op(A)' ) - I(N**2).
C
C     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the
C     Kronecker product. The program estimates sigma_min(T) by the
C     reciprocal of an estimate of the 1-norm of inverse(T). The true
C     reciprocal 1-norm of inverse(T) cannot differ from sigma_min(T) by
C     more than a factor of N.
C
C     When SEPD is small, small changes in A, C can cause large changes
C     in the solution of the equation. An approximate bound on the
C     maximum relative error in the computed solution is
C
C                            EPS * norm(A)**2 / SEPD
C
C     where EPS is the machine precision.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Supersedes Release 2.0 routine MB03AD by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, October 1982.
C     Based on DGELPD by P. Petkov, Tech. University of Sofia, September
C     1993.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999.
C
C     KEYWORDS
C
C     Lyapunov equation, orthogonal transformation, real Schur form,
C     Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, TRANA
      INTEGER            INFO, LDA, LDC, LDU, LDWORK, N
      DOUBLE PRECISION   FERR, SCALE, SEPD
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                   U( LDU, * ), WI( * ), WR( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOFACT, NOTA, WANTBH, WANTSP, WANTX
      CHARACTER          NOTRA, UPLO
      INTEGER            I, IERR, KASE, LWA, MINWRK, SDIM
      DOUBLE PRECISION   EST, SCALEF
C     ..
C     .. Local Arrays ..
      LOGICAL            BWORK( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, SELECT
      DOUBLE PRECISION   DLAMCH, DLANHS
      EXTERNAL           DLAMCH, DLANHS, LSAME, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEES, DLACON, MB01RD, SB03MX, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      NOFACT = LSAME( FACT,  'N' )
      NOTA   = LSAME( TRANA, 'N' )
C
      INFO = 0
      IF( .NOT.WANTBH .AND. .NOT.WANTSP .AND. .NOT.WANTX ) THEN
         INFO = -1
      ELSE IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTA .AND. .NOT.LSAME( TRANA, 'T' ) .AND.
     $                         .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( WANTSP .AND. LDC.LT.1 .OR.
     $    .NOT.WANTSP .AND. LDC.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
C
C     Compute workspace.
C
      IF( WANTX ) THEN
         IF( NOFACT ) THEN
            MINWRK = MAX( N*N, 3*N )
         ELSE
            MINWRK = MAX( N*N, 2*N )
         END IF
      ELSE
         MINWRK = 2*N*N + 2*N
      END IF
      IF( LDWORK.LT.MAX( 1, MINWRK ) ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB03PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         SCALE = ONE
         IF( WANTBH )
     $      FERR  = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
      LWA = 0
C
      IF( NOFACT ) THEN
C
C        Compute the Schur factorization of A.
C        Workspace:  need   3*N;
C                    prefer larger.
C
         CALL DGEES( 'Vectors', 'Not ordered', SELECT, N, A, LDA, SDIM,
     $               WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         IF( INFO.GT.0 )
     $      RETURN
         LWA = INT( DWORK( 1 ) )
      END IF
C
      IF( .NOT.WANTSP ) THEN
C
C        Transform the right-hand side.
C        Workspace:  need   N*N.
C
         UPLO   = 'U'
         CALL MB01RD( UPLO, 'Transpose', N, N, ZERO, ONE, C, LDC, U,
     $                LDU, C, LDC, DWORK, LDWORK, INFO )
C
         DO 10 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   10    CONTINUE
C
C        Solve the transformed equation.
C        Workspace:  2*N.
C
         CALL SB03MX( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
         IF( INFO.GT.0 )
     $      INFO = N + 1
C
C        Transform back the solution.
C
         CALL MB01RD( UPLO, 'No transpose', N, N, ZERO, ONE, C, LDC, U,
     $                LDU, C, LDC, DWORK, LDWORK, INFO )
C
         DO 20 I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
   20    CONTINUE
C
      END IF
C
      IF( .NOT.WANTX ) THEN
C
C        Estimate sepd(op(A),op(A)').
C        Workspace:  2*N*N + 2*N.
C
         IF( NOTA ) THEN
            NOTRA = 'T'
         ELSE
            NOTRA = 'N'
         END IF
C
         EST = ZERO
         KASE = 0
C        REPEAT
   30    CONTINUE
         CALL DLACON( N*N, DWORK( N*N+1 ), DWORK, IWORK, EST, KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
               CALL SB03MX( TRANA, N, A, LDA, DWORK, N, SCALEF,
     $                      DWORK( 2*N*N + 1 ), IERR )
            ELSE
               CALL SB03MX( NOTRA, N, A, LDA, DWORK, N, SCALEF,
     $                      DWORK( 2*N*N + 1 ), IERR )
            END IF
            GO TO 30
         END IF
C        UNTIL KASE = 0
C
         SEPD = SCALEF / EST
C
         IF( WANTBH ) THEN
C
C           Compute the estimate of the relative error.
C
            FERR = DLAMCH( 'Precision' )*
     $             DLANHS( 'Frobenius', N, A, LDA, DWORK )**2 / SEPD
         END IF
      END IF
C
      DWORK( 1 ) = DBLE( MAX( LWA, MINWRK ) )
C
      RETURN
C *** Last line of SB03PD ***
      END
