      SUBROUTINE MB02XD( FORM, STOR, UPLO, F, M, N, NRHS, IPAR, LIPAR,
     $                   DPAR, LDPAR, A, LDA, B, LDB, ATA, LDATA, DWORK,
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
C     To solve a set of systems of linear equations, A'*A*X = B, or,
C     in the implicit form, f(A)*X = B, with A'*A or f(A) positive
C     definite, using symmetric Gaussian elimination.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     FORM    CHARACTER*1
C             Specifies the form in which the matrix A is provided, as
C             follows:
C             = 'S' :  standard form, the matrix A is given;
C             = 'F' :  the implicit, function form f(A) is provided.
C             If FORM = 'F', then the routine F is called to compute the
C             matrix A'*A.
C
C     STOR    CHARACTER*1
C             Specifies the storage scheme for the symmetric
C             matrix A'*A, as follows:
C             = 'F' :  full storage is used;
C             = 'P' :  packed storage is used.
C
C     UPLO    CHARACTER*1
C             Specifies which part of the matrix A'*A is stored, as
C             follows:
C             = 'U' :  the upper triagular part is stored;
C             = 'L' :  the lower triagular part is stored.
C
C     Function Parameters
C
C     F       EXTERNAL
C             If FORM = 'F', then F is a subroutine which calculates the
C             value of f(A) = A'*A, for given A.
C             If FORM = 'S', then F is not called.
C
C             F must have the following interface:
C
C             SUBROUTINE F( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, A,
C            $              LDA, ATA, LDATA, DWORK, LDWORK, INFO )
C
C             where
C
C             STOR    (input) CHARACTER*1
C                     Specifies the storage scheme for the symmetric
C                     matrix A'*A, as follows:
C                     = 'F' :  full storage is used;
C                     = 'P' :  packed storage is used.
C
C             UPLO    (input) CHARACTER*1
C                     Specifies which part of the matrix A'*A is stored,
C                     as follows:
C                     = 'U' :  the upper triagular part is stored;
C                     = 'L' :  the lower triagular part is stored.
C
C             N       (input) INTEGER
C                     The order of the matrix A'*A.  N >= 0.
C
C             IPAR    (input) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the matrix A.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 0.
C
C             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C                     The real parameters needed for solving the
C                     problem.
C
C             LDPAR   (input) INTEGER
C                     The length of the array DPAR.  LDPAR >= 0.
C
C             A       (input) DOUBLE PRECISION array, dimension
C                     (LDA, NC), where NC is the number of columns.
C                     The leading NR-by-NC part of this array must
C                     contain the (compressed) representation of the
C                     matrix A, where NR is the number of rows of A
C                     (function of IPAR entries).
C
C             LDA     (input) INTEGER
C                     The leading dimension of the array A.
C                     LDA >= MAX(1,NR).
C
C             ATA     (output) DOUBLE PRECISION array,
C                              dimension (LDATA,N),    if STOR = 'F',
C                              dimension (N*(N+1)/2),  if STOR = 'P'.
C                     The leading N-by-N (if STOR = 'F'), or N*(N+1)/2
C                     (if STOR = 'P') part of this array contains the
C                     upper or lower triangle of the matrix A'*A,
C                     depending on UPLO = 'U', or UPLO = 'L',
C                     respectively, stored either as a two-dimensional,
C                     or one-dimensional array, depending on STOR.
C
C             LDATA   (input) INTEGER
C                     The leading dimension of the array ATA.
C                     LDATA >= MAX(1,N), if STOR = 'F'.
C                     LDATA >= 1,        if STOR = 'P'.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine F.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine F).
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input scalar argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine F. The LAPACK Library routine XERBLA
C                     should be used in conjunction with negative INFO.
C                     INFO must be zero if the subroutine finished
C                     successfully.
C
C             Parameters marked with "(input)" must not be changed.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix A'*A, the number of columns of the
C             matrix A, and the number of rows of the matrix X.  N >= 0.
C
C     NRHS    (input) INTEGER
C             The number of columns of the matrices B and X.  NRHS >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             If FORM = 'F', the integer parameters describing the
C             structure of the matrix A.
C             This parameter is ignored if FORM = 'S'.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 0.
C
C     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C             If FORM = 'F', the real parameters needed for solving
C             the problem.
C             This parameter is ignored if FORM = 'S'.
C
C     LDPAR   (input) INTEGER
C             The length of the array DPAR.  LDPAR >= 0.
C
C     A       (input) DOUBLE PRECISION array,
C                     dimension (LDA, N),  if FORM = 'S',
C                     dimension (LDA, NC), if FORM = 'F', where NC is
C                     the number of columns.
C             If FORM = 'S', the leading M-by-N part of this array
C             must contain the matrix A.
C             If FORM = 'F', the leading NR-by-NC part of this array
C             must contain an appropriate representation of matrix A,
C             where NR is the number of rows.
C             If FORM = 'F', this array is not referenced by this
C             routine itself, except in the call to the routine F.
C
C     LDA     INTEGER
C             The leading dimension of array A.
C             LDA >= MAX(1,M),  if FORM = 'S';
C             LDA >= MAX(1,NR), if FORM = 'F'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB, NRHS)
C             On entry, the leading N-by-NRHS part of this array must
C             contain the right hand side matrix B.
C             On exit, if INFO = 0 and M (or NR) is nonzero, the leading
C             N-by-NRHS part of this array contains the solution X of
C             the set of systems of linear equations A'*A*X = B or
C             f(A)*X = B. If M (or NR) is zero, then B is unchanged.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     ATA     (output) DOUBLE PRECISION array,
C                      dimension (LDATA,N),    if STOR = 'F',
C                      dimension (N*(N+1)/2),  if STOR = 'P'.
C             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if
C             STOR = 'P') part of this array contains the upper or lower
C             triangular Cholesky factor of the matrix A'*A, depending
C             on UPLO = 'U', or UPLO = 'L', respectively, stored either
C             as a two-dimensional, or one-dimensional array, depending
C             on STOR.
C
C     LDATA   INTEGER
C             The leading dimension of the array ATA.
C             LDATA >= MAX(1,N), if STOR = 'F'.
C             LDATA >= 1,        if STOR = 'P'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, then the (i,i) element of the
C                   triangular factor of the matrix A'*A is exactly
C                   zero (the matrix A'*A is exactly singular);
C                   if INFO = j > n, then F returned with INFO = j-n.
C
C     METHOD
C
C     The matrix A'*A is built either directly (if FORM = 'S'), or
C     implicitly, by calling the routine F. Then, A'*A is Cholesky
C     factored and its factor is used to solve the set of systems of
C     linear equations, A'*A*X = B.
C
C     REFERENCES
C
C     [1] Golub, G.H. and van Loan, C.F.
C         Matrix Computations. Third Edition.
C         M. D. Johns Hopkins University Press, Baltimore, 1996.
C
C     [2] Anderson, E., Bai, Z., Bischof, C., Blackford, Demmel, J.,
C         Dongarra, J., Du Croz, J., Greenbaum, A., Hammarling, S.,
C         McKenney, A., Sorensen, D.
C         LAPACK Users' Guide: Third Edition, SIAM, Philadelphia, 1999.
C
C     NUMERICAL ASPECTS
C
C     For speed, this routine does not check for near singularity of the
C     matrix A'*A. If the matrix A is nearly rank deficient, then the
C     computed X could be inaccurate. Estimates of the reciprocal
C     condition numbers of the matrices A and A'*A can be obtained
C     using LAPACK routines DGECON and DPOCON (DPPCON), respectively.
C
C     The approximate number of floating point operations is
C        (M+3)*N**2/2 + N**3/6 + NRHS*N**2, if FORM = 'S',
C        f + N**3/6 + NRHS*N**2,            if FORM = 'F',
C     where M is the number of rows of A, and f is the number of
C     floating point operations required by the subroutine F.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     REVISIONS
C
C     V. Sima, Mar. 2002.
C
C     KEYWORDS
C
C     Linear system of equations, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         FORM, STOR, UPLO
      INTEGER           INFO, LDA, LDATA, LDB, LDPAR, LDWORK, LIPAR, M,
     $                  N, NRHS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ATA(*), B(LDB,*), DPAR(*), DWORK(*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      INTEGER           IERR, J, J1
      LOGICAL           FULL, MAT, UPPER
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DPOTRF, DPOTRS, DPPTRF, DPPTRS, DSYRK, F,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      MAT   = LSAME( FORM, 'S' )
      FULL  = LSAME( STOR, 'F' )
      UPPER = LSAME( UPLO, 'U' )
C
C     Check the scalar input parameters.
C
      INFO = 0
      IF( .NOT.( MAT .OR. LSAME( FORM, 'F' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( FULL .OR. LSAME( STOR, 'P' ) ) ) THEN
         INFO = -2
      ELSEIF ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSEIF ( M.LT.0 ) THEN
         INFO = -5
      ELSEIF ( N.LT.0 ) THEN
         INFO = -6
      ELSEIF ( NRHS.LT.0 ) THEN
         INFO = -7
      ELSEIF ( .NOT. MAT .AND. LIPAR.LT.0 ) THEN
         INFO = -9
      ELSEIF ( .NOT. MAT .AND. LDPAR.LT.0 ) THEN
         INFO = -11
      ELSEIF ( LDA.LT.1 .OR. ( MAT .AND. LDA.LT.M ) ) THEN
         INFO = -13
      ELSEIF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSEIF ( LDATA.LT.1 .OR. ( FULL .AND. LDATA.LT.N ) ) THEN
         INFO = -17
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02XD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. ( MAT .AND. M.EQ.0 ) )
     $   RETURN
C
C     Build a triangle of the matrix A'*A.
C
      IF ( MAT ) THEN
C
C        Matrix A is given in the usual form.
C
         IF ( FULL ) THEN
            CALL DSYRK( UPLO, 'Transpose', N, M, ONE, A, LDA, ZERO,
     $                  ATA, LDATA )
         ELSEIF ( UPPER ) THEN
            J1 = 1
C
            DO 10 J = 1, N
               CALL DGEMV( 'Transpose', M, J, ONE, A, LDA, A(1,J), 1,
     $                     ZERO, ATA(J1), 1 )
               J1 = J1 + J
   10       CONTINUE
C
         ELSE
            J1 = 1
C
            DO 20 J = 1, N
               CALL DGEMV( 'Transpose', M, N-J+1, ONE, A(1,J), LDA,
     $                     A(1,J), 1, ZERO, ATA(J1), 1 )
               J1 = J1 + N - J + 1
   20       CONTINUE
C
         ENDIF
C
      ELSE
C
C        Implicit form, A'*A = f(A).
C
         CALL F( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, A, LDA, ATA,
     $           LDATA, DWORK, LDWORK, IERR )
         IF ( IERR.NE.0 ) THEN
            INFO = N + IERR
            RETURN
         ENDIF
C
      ENDIF
C
C     Factor the matrix A'*A.
C
      IF ( FULL ) THEN
         CALL DPOTRF( UPLO, N, ATA, LDATA, IERR )
      ELSE
         CALL DPPTRF( UPLO, N, ATA, IERR )
      ENDIF
C
      IF ( IERR.NE.0 ) THEN
         INFO = IERR
         RETURN
      ENDIF
C
C     Solve the set of linear systems.
C
      IF ( FULL ) THEN
         CALL DPOTRS( UPLO, N, NRHS, ATA, LDATA, B, LDB, IERR )
      ELSE
         CALL DPPTRS( UPLO, N, NRHS, ATA, B, LDB, IERR )
      ENDIF
C
C *** Last line of MB02XD ***
      END
