/* MB02XD.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb02xd_(char *form, char *stor, char *uplo, S_fp f, 
	integer *m, integer *n, integer *nrhs, integer *ipar, integer *lipar, 
	doublereal *dpar, integer *ldpar, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *ata, integer *ldata, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen form_len, 
	ftnlen stor_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer j, j1;
    static logical mat;
    static integer ierr;
    static logical full;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    dpotrf_(char *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dpptrf_(char *, integer *, doublereal *, integer *, 
	    ftnlen), dpotrs_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), dpptrs_(
	    char *, integer *, integer *, doublereal *, doublereal *, integer 
	    *, integer *, ftnlen);


/*     SLICOT RELEASE 5.0. */

/*     Copyright (c) 2002-2009 NICONET e.V. */

/*     This program is free software: you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation, either version 2 of */
/*     the License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see */
/*     <http://www.gnu.org/licenses/>. */

/*     PURPOSE */

/*     To solve a set of systems of linear equations, A'*A*X = B, or, */
/*     in the implicit form, f(A)*X = B, with A'*A or f(A) positive */
/*     definite, using symmetric Gaussian elimination. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FORM    CHARACTER*1 */
/*             Specifies the form in which the matrix A is provided, as */
/*             follows: */
/*             = 'S' :  standard form, the matrix A is given; */
/*             = 'F' :  the implicit, function form f(A) is provided. */
/*             If FORM = 'F', then the routine F is called to compute the */
/*             matrix A'*A. */

/*     STOR    CHARACTER*1 */
/*             Specifies the storage scheme for the symmetric */
/*             matrix A'*A, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the matrix A'*A is stored, as */
/*             follows: */
/*             = 'U' :  the upper triagular part is stored; */
/*             = 'L' :  the lower triagular part is stored. */

/*     Function Parameters */

/*     F       EXTERNAL */
/*             If FORM = 'F', then F is a subroutine which calculates the */
/*             value of f(A) = A'*A, for given A. */
/*             If FORM = 'S', then F is not called. */

/*             F must have the following interface: */

/*             SUBROUTINE F( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, A, */
/*            $              LDA, ATA, LDATA, DWORK, LDWORK, INFO ) */

/*             where */

/*             STOR    (input) CHARACTER*1 */
/*                     Specifies the storage scheme for the symmetric */
/*                     matrix A'*A, as follows: */
/*                     = 'F' :  full storage is used; */
/*                     = 'P' :  packed storage is used. */

/*             UPLO    (input) CHARACTER*1 */
/*                     Specifies which part of the matrix A'*A is stored, */
/*                     as follows: */
/*                     = 'U' :  the upper triagular part is stored; */
/*                     = 'L' :  the lower triagular part is stored. */

/*             N       (input) INTEGER */
/*                     The order of the matrix A'*A.  N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the matrix A. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*                     The real parameters needed for solving the */
/*                     problem. */

/*             LDPAR   (input) INTEGER */
/*                     The length of the array DPAR.  LDPAR >= 0. */

/*             A       (input) DOUBLE PRECISION array, dimension */
/*                     (LDA, NC), where NC is the number of columns. */
/*                     The leading NR-by-NC part of this array must */
/*                     contain the (compressed) representation of the */
/*                     matrix A, where NR is the number of rows of A */
/*                     (function of IPAR entries). */

/*             LDA     (input) INTEGER */
/*                     The leading dimension of the array A. */
/*                     LDA >= MAX(1,NR). */

/*             ATA     (output) DOUBLE PRECISION array, */
/*                              dimension (LDATA,N),    if STOR = 'F', */
/*                              dimension (N*(N+1)/2),  if STOR = 'P'. */
/*                     The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 */
/*                     (if STOR = 'P') part of this array contains the */
/*                     upper or lower triangle of the matrix A'*A, */
/*                     depending on UPLO = 'U', or UPLO = 'L', */
/*                     respectively, stored either as a two-dimensional, */
/*                     or one-dimensional array, depending on STOR. */

/*             LDATA   (input) INTEGER */
/*                     The leading dimension of the array ATA. */
/*                     LDATA >= MAX(1,N), if STOR = 'F'. */
/*                     LDATA >= 1,        if STOR = 'P'. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine F. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine F). */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input scalar argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine F. The LAPACK Library routine XERBLA */
/*                     should be used in conjunction with negative INFO. */
/*                     INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix A'*A, the number of columns of the */
/*             matrix A, and the number of rows of the matrix X.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of columns of the matrices B and X.  NRHS >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             If FORM = 'F', the integer parameters describing the */
/*             structure of the matrix A. */
/*             This parameter is ignored if FORM = 'S'. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 0. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             If FORM = 'F', the real parameters needed for solving */
/*             the problem. */
/*             This parameter is ignored if FORM = 'S'. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 0. */

/*     A       (input) DOUBLE PRECISION array, */
/*                     dimension (LDA, N),  if FORM = 'S', */
/*                     dimension (LDA, NC), if FORM = 'F', where NC is */
/*                     the number of columns. */
/*             If FORM = 'S', the leading M-by-N part of this array */
/*             must contain the matrix A. */
/*             If FORM = 'F', the leading NR-by-NC part of this array */
/*             must contain an appropriate representation of matrix A, */
/*             where NR is the number of rows. */
/*             If FORM = 'F', this array is not referenced by this */
/*             routine itself, except in the call to the routine F. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. */
/*             LDA >= MAX(1,M),  if FORM = 'S'; */
/*             LDA >= MAX(1,NR), if FORM = 'F'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB, NRHS) */
/*             On entry, the leading N-by-NRHS part of this array must */
/*             contain the right hand side matrix B. */
/*             On exit, if INFO = 0 and M (or NR) is nonzero, the leading */
/*             N-by-NRHS part of this array contains the solution X of */
/*             the set of systems of linear equations A'*A*X = B or */
/*             f(A)*X = B. If M (or NR) is zero, then B is unchanged. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     ATA     (output) DOUBLE PRECISION array, */
/*                      dimension (LDATA,N),    if STOR = 'F', */
/*                      dimension (N*(N+1)/2),  if STOR = 'P'. */
/*             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if */
/*             STOR = 'P') part of this array contains the upper or lower */
/*             triangular Cholesky factor of the matrix A'*A, depending */
/*             on UPLO = 'U', or UPLO = 'L', respectively, stored either */
/*             as a two-dimensional, or one-dimensional array, depending */
/*             on STOR. */

/*     LDATA   INTEGER */
/*             The leading dimension of the array ATA. */
/*             LDATA >= MAX(1,N), if STOR = 'F'. */
/*             LDATA >= 1,        if STOR = 'P'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then the (i,i) element of the */
/*                   triangular factor of the matrix A'*A is exactly */
/*                   zero (the matrix A'*A is exactly singular); */
/*                   if INFO = j > n, then F returned with INFO = j-n. */

/*     METHOD */

/*     The matrix A'*A is built either directly (if FORM = 'S'), or */
/*     implicitly, by calling the routine F. Then, A'*A is Cholesky */
/*     factored and its factor is used to solve the set of systems of */
/*     linear equations, A'*A*X = B. */

/*     REFERENCES */

/*     [1] Golub, G.H. and van Loan, C.F. */
/*         Matrix Computations. Third Edition. */
/*         M. D. Johns Hopkins University Press, Baltimore, 1996. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Blackford, Demmel, J., */
/*         Dongarra, J., Du Croz, J., Greenbaum, A., Hammarling, S., */
/*         McKenney, A., Sorensen, D. */
/*         LAPACK Users' Guide: Third Edition, SIAM, Philadelphia, 1999. */

/*     NUMERICAL ASPECTS */

/*     For speed, this routine does not check for near singularity of the */
/*     matrix A'*A. If the matrix A is nearly rank deficient, then the */
/*     computed X could be inaccurate. Estimates of the reciprocal */
/*     condition numbers of the matrices A and A'*A can be obtained */
/*     using LAPACK routines DGECON and DPOCON (DPPCON), respectively. */

/*     The approximate number of floating point operations is */
/*        (M+3)*N**2/2 + N**3/6 + NRHS*N**2, if FORM = 'S', */
/*        f + N**3/6 + NRHS*N**2,            if FORM = 'F', */
/*     where M is the number of rows of A, and f is the number of */
/*     floating point operations required by the subroutine F. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     V. Sima, Mar. 2002. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    --ipar;
    --dpar;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --ata;
    --dwork;

    /* Function Body */
    mat = lsame_(form, "S", (ftnlen)1, (ftnlen)1);
    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    *info = 0;
    if (! (mat || lsame_(form, "F", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*nrhs < 0) {
	*info = -7;
    } else if (! mat && *lipar < 0) {
	*info = -9;
    } else if (! mat && *ldpar < 0) {
	*info = -11;
    } else if (*lda < 1 || mat && *lda < *m) {
	*info = -13;
    } else if (*ldb < max(1,*n)) {
	*info = -15;
    } else if (*ldata < 1 || full && *ldata < *n) {
	*info = -17;
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02XD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || mat && *m == 0) {
	return 0;
    }

/*     Build a triangle of the matrix A'*A. */

    if (mat) {

/*        Matrix A is given in the usual form. */

	if (full) {
	    dsyrk_(uplo, "Transpose", n, m, &c_b10, &a[a_offset], lda, &c_b11,
		     &ata[1], ldata, (ftnlen)1, (ftnlen)9);
	} else if (upper) {
	    j1 = 1;

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dgemv_("Transpose", m, &j, &c_b10, &a[a_offset], lda, &a[j * 
			a_dim1 + 1], &c__1, &c_b11, &ata[j1], &c__1, (ftnlen)
			9);
		j1 += j;
/* L10: */
	    }

	} else {
	    j1 = 1;

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		dgemv_("Transpose", m, &i__2, &c_b10, &a[j * a_dim1 + 1], lda,
			 &a[j * a_dim1 + 1], &c__1, &c_b11, &ata[j1], &c__1, (
			ftnlen)9);
		j1 = j1 + *n - j + 1;
/* L20: */
	    }

	}

    } else {

/*        Implicit form, A'*A = f(A). */

	(*f)(stor, uplo, n, &ipar[1], lipar, &dpar[1], ldpar, &a[a_offset], 
		lda, &ata[1], ldata, &dwork[1], ldwork, &ierr, (ftnlen)1, (
		ftnlen)1);
	if (ierr != 0) {
	    *info = *n + ierr;
	    return 0;
	}

    }

/*     Factor the matrix A'*A. */

    if (full) {
	dpotrf_(uplo, n, &ata[1], ldata, &ierr, (ftnlen)1);
    } else {
	dpptrf_(uplo, n, &ata[1], &ierr, (ftnlen)1);
    }

    if (ierr != 0) {
	*info = ierr;
	return 0;
    }

/*     Solve the set of linear systems. */

    if (full) {
	dpotrs_(uplo, n, nrhs, &ata[1], ldata, &b[b_offset], ldb, &ierr, (
		ftnlen)1);
    } else {
	dpptrs_(uplo, n, nrhs, &ata[1], &b[b_offset], ldb, &ierr, (ftnlen)1);
    }

/* *** Last line of MB02XD *** */
    return 0;
} /* mb02xd_ */

