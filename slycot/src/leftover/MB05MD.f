      SUBROUTINE MB05MD( BALANC, N, DELTA, A, LDA, V, LDV, Y, LDY, VALR,
     $                   VALI, IWORK, DWORK, LDWORK, INFO )
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
C     To compute exp(A*delta) where A is a real N-by-N non-defective
C     matrix with real or complex eigenvalues and delta is a scalar
C     value. The routine also returns the eigenvalues and eigenvectors
C     of A as well as (if all eigenvalues are real) the matrix product
C     exp(Lambda*delta) times the inverse of the eigenvector matrix
C     of A, where Lambda is the diagonal matrix of eigenvalues.
C     Optionally, the routine computes a balancing transformation to
C     improve the conditioning of the eigenvalues and eigenvectors.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALANC  CHARACTER*1
C             Indicates how the input matrix should be diagonally scaled
C             to improve the conditioning of its eigenvalues as follows:
C             = 'N':  Do not diagonally scale;
C             = 'S':  Diagonally scale the matrix, i.e. replace A by
C                     D*A*D**(-1), where D is a diagonal matrix chosen
C                     to make the rows and columns of A more equal in
C                     norm. Do not permute.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     DELTA   (input) DOUBLE PRECISION
C             The scalar value delta of the problem.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A of the problem.
C             On exit, the leading N-by-N part of this array contains
C             the solution matrix exp(A*delta).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,N).
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,N)
C             The leading N-by-N part of this array contains the
C             eigenvector matrix for A.
C             If the k-th eigenvalue is real the k-th column of the
C             eigenvector matrix holds the eigenvector corresponding
C             to the k-th eigenvalue.
C             Otherwise, the k-th and (k+1)-th eigenvalues form a
C             complex conjugate pair and the k-th and (k+1)-th columns
C             of the eigenvector matrix hold the real and imaginary
C             parts of the eigenvectors corresponding to these
C             eigenvalues as follows.
C             If p and q denote the k-th and (k+1)-th columns of the
C             eigenvector matrix, respectively, then the eigenvector
C             corresponding to the complex eigenvalue with positive
C             (negative) imaginary value is given by
C                                       2
C             p + q*j (p - q*j), where j  = -1.
C
C     LDV     INTEGER
C             The leading dimension of array V.  LDV >= max(1,N).
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY,N)
C             The leading N-by-N part of this array contains an
C             intermediate result for computing the matrix exponential.
C             Specifically, exp(A*delta) is obtained as the product V*Y,
C             where V is the matrix stored in the leading N-by-N part of
C             the array V. If all eigenvalues of A are real, then the
C             leading N-by-N part of this array contains the matrix
C             product exp(Lambda*delta) times the inverse of the (right)
C             eigenvector matrix of A, where Lambda is the diagonal
C             matrix of eigenvalues.
C
C     LDY     INTEGER
C             The leading dimension of array Y.  LDY >= max(1,N).
C
C     VALR    (output) DOUBLE PRECISION array, dimension (N)
C     VALI    (output) DOUBLE PRECISION array, dimension (N)
C             These arrays contain the real and imaginary parts,
C             respectively, of the eigenvalues of the matrix A. The
C             eigenvalues are unordered except that complex conjugate
C             pairs of values appear consecutively with the eigenvalue
C             having positive imaginary part first.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, and if N > 0, DWORK(2) returns the reciprocal
C             condition number of the triangular matrix used to obtain
C             the inverse of the eigenvector matrix.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= max(1,4*N).
C             For good performance, LDWORK must generally be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = i:  if INFO = i, the QR algorithm failed to compute all
C                   the eigenvalues; no eigenvectors have been computed;
C                   elements i+1:N of VALR and VALI contain eigenvalues
C                   which have converged;
C             = N+1:  if the inverse of the eigenvector matrix could not
C                   be formed due to an attempt to divide by zero, i.e.,
C                   the eigenvector matrix is singular;
C             = N+2:  if the matrix A is defective, possibly due to
C                   rounding errors.
C
C     METHOD
C
C     This routine is an implementation of "Method 15" of the set of
C     methods described in reference [1], which uses an eigenvalue/
C     eigenvector decomposition technique. A modification of LAPACK
C     Library routine DGEEV is used for obtaining the right eigenvector
C     matrix. A condition estimate is then employed to determine if the
C     matrix A is near defective and hence the exponential solution is
C     inaccurate. In this case the routine returns with the Error
C     Indicator (INFO) set to N+2, and SLICOT Library routines MB05ND or
C     MB05OD are the preferred alternative routines to be used.
C
C     REFERENCES
C
C     [1] Moler, C.B. and Van Loan, C.F.
C         Nineteen dubious ways to compute the exponential of a matrix.
C         SIAM Review, 20, pp. 801-836, 1978.
C
C     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05AD by M.J. Denham, Kingston
C     Polytechnic, March 1981.
C
C     REVISIONS
C
C     V. Sima, June 13, 1997, April 25, 2003, Feb. 15, 2004.
C
C     KEYWORDS
C
C     Eigenvalue, eigenvector decomposition, matrix exponential.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         BALANC
      INTEGER           INFO, LDA, LDV, LDWORK, LDY, N
      DOUBLE PRECISION  DELTA
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), V(LDV,*), VALI(*), VALR(*),
     $                  Y(LDY,*)
C     .. Local Scalars ..
      LOGICAL           SCALE
      INTEGER           I
      DOUBLE PRECISION  RCOND, TEMPI, TEMPR, WRKOPT
C     .. Local Arrays ..
      DOUBLE PRECISION  TMP(2,2)
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEBAK, DGEMM, DLACPY, DSCAL, DSWAP, DTRCON,
     $                  DTRMM, DTRSM, MB05MY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         COS, EXP, MAX, SIN
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO  = 0
      SCALE = LSAME( BALANC, 'S' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. SCALE ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDY.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDWORK.LT.MAX( 1, 4*N ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB05MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of workspace needed at that point in the code,
C     as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
C     Compute the eigenvalues and right eigenvectors of the real
C     nonsymmetric matrix A; optionally, compute a balancing
C     transformation.
C     Workspace:  need: 4*N.
C
      CALL MB05MY( BALANC, N, A, LDA, VALR, VALI, V, LDV, Y, LDY,
     $             DWORK, LDWORK, INFO )
C
      IF ( INFO.GT.0 )
     $   RETURN
      WRKOPT = DWORK(1)
      IF ( SCALE ) THEN
         DO 10 I = 1, N
            DWORK(I) = DWORK(I+1)
   10    CONTINUE
      END IF
C
C     Exit with INFO = N + 1 if V is exactly singular.
C
      DO 20 I = 1, N
         IF ( V(I,I).EQ.ZERO ) THEN
            INFO = N + 1
            RETURN
         END IF
   20 CONTINUE
C
C     Compute the reciprocal condition number of the triangular matrix.
C
      CALL DTRCON( '1-norm', 'Upper', 'Non unit', N, V, LDV, RCOND,
     $             DWORK(N+1), IWORK, INFO )
C
C     Return if the matrix is singular to working precision.
C
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
         DWORK(2) = RCOND
         INFO = N + 2
         RETURN
      END IF
C
C     Compute the right eigenvector matrix (temporarily) in A.
C
      CALL DLACPY( 'Full', N, N, Y, LDY, A, LDA )
      CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Non unit', N, N,
     $            ONE, V, LDV, A, LDA )
      IF ( SCALE )
     $   CALL DGEBAK( BALANC, 'Right', N, 1, N, DWORK, N, A, LDA, INFO )
C
C     Compute the inverse of the right eigenvector matrix, by solving
C     a set of linear systems, V * X = Y' (if BALANC = 'N').
C
      DO 40 I = 2, N
         CALL DSWAP( I-1, Y(I,1), LDY, Y(1,I), 1 )
   40 CONTINUE
C
      CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non unit', N, N,
     $            ONE, V, LDV, Y, LDY )
      IF( SCALE ) THEN
C
         DO 60 I = 1, N
            TEMPR = ONE / DWORK(I)
            CALL DSCAL( N, TEMPR, Y(1,I), 1 )
   60    CONTINUE
C
      END IF
C
C     Save the right eigenvector matrix in V.
C
      CALL DLACPY( 'Full', N, N, A, LDA, V, LDV )
C
C     Premultiply the inverse eigenvector matrix by the exponential of
C     quasi-diagonal matrix Lambda * DELTA, where Lambda is the matrix
C     of eigenvalues.
C     Note that only real arithmetic is used, taking the special storing
C     of eigenvalues/eigenvectors into account.
C
      I = 0
C     REPEAT
   80 CONTINUE
      I = I + 1
      IF ( VALI(I).EQ.ZERO ) THEN
         TEMPR = EXP( VALR(I)*DELTA )
         CALL DSCAL( N, TEMPR, Y(I,1), LDY )
      ELSE
         TEMPR = VALR(I)*DELTA
         TEMPI = VALI(I)*DELTA
         TMP(1,1) = COS( TEMPI )*EXP( TEMPR )
         TMP(1,2) = SIN( TEMPI )*EXP( TEMPR )
         TMP(2,1) = -TMP(1,2)
         TMP(2,2) =  TMP(1,1)
         CALL DLACPY( 'Full', 2, N, Y(I,1), LDY, DWORK, 2 )
         CALL DGEMM( 'No transpose', 'No transpose', 2, N, 2, ONE,
     $               TMP, 2, DWORK, 2, ZERO, Y(I,1), LDY )
         I = I + 1
      END IF
      IF ( I.LT.N ) GO TO 80
C     UNTIL I = N.
C
C     Compute the matrix exponential as the product V * Y.
C
      CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE, V, LDV,
     $            Y, LDY, ZERO, A, LDA )
C
C     Set optimal workspace dimension and reciprocal condition number.
C
      DWORK(1) = WRKOPT
      DWORK(2) = RCOND
C
      RETURN
C *** Last line of MB05MD ***
      END
