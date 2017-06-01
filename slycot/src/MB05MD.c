/* MB05MD.f -- translated by f2c (version 20100827).
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

static doublereal c_b16 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b40 = 0.;

/* Subroutine */ int mb05md_(char *balanc, integer *n, doublereal *delta, 
	doublereal *a, integer *lda, doublereal *v, integer *ldv, doublereal *
	y, integer *ldy, doublereal *valr, doublereal *vali, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen balanc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, y_dim1, y_offset, i__1, i__2;

    /* Builtin functions */
    double exp(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal tmp[4]	/* was [2][2] */;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb05my_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    static doublereal tempi;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal tempr;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgebak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrcon_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal wrkopt;


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

/*     To compute exp(A*delta) where A is a real N-by-N non-defective */
/*     matrix with real or complex eigenvalues and delta is a scalar */
/*     value. The routine also returns the eigenvalues and eigenvectors */
/*     of A as well as (if all eigenvalues are real) the matrix product */
/*     exp(Lambda*delta) times the inverse of the eigenvector matrix */
/*     of A, where Lambda is the diagonal matrix of eigenvalues. */
/*     Optionally, the routine computes a balancing transformation to */
/*     improve the conditioning of the eigenvalues and eigenvectors. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Indicates how the input matrix should be diagonally scaled */
/*             to improve the conditioning of its eigenvalues as follows: */
/*             = 'N':  Do not diagonally scale; */
/*             = 'S':  Diagonally scale the matrix, i.e. replace A by */
/*                     D*A*D**(-1), where D is a diagonal matrix chosen */
/*                     to make the rows and columns of A more equal in */
/*                     norm. Do not permute. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             The scalar value delta of the problem. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A of the problem. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the solution matrix exp(A*delta). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,N). */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,N) */
/*             The leading N-by-N part of this array contains the */
/*             eigenvector matrix for A. */
/*             If the k-th eigenvalue is real the k-th column of the */
/*             eigenvector matrix holds the eigenvector corresponding */
/*             to the k-th eigenvalue. */
/*             Otherwise, the k-th and (k+1)-th eigenvalues form a */
/*             complex conjugate pair and the k-th and (k+1)-th columns */
/*             of the eigenvector matrix hold the real and imaginary */
/*             parts of the eigenvectors corresponding to these */
/*             eigenvalues as follows. */
/*             If p and q denote the k-th and (k+1)-th columns of the */
/*             eigenvector matrix, respectively, then the eigenvector */
/*             corresponding to the complex eigenvalue with positive */
/*             (negative) imaginary value is given by */
/*                                       2 */
/*             p + q*j (p - q*j), where j  = -1. */

/*     LDV     INTEGER */
/*             The leading dimension of array V.  LDV >= max(1,N). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains an */
/*             intermediate result for computing the matrix exponential. */
/*             Specifically, exp(A*delta) is obtained as the product V*Y, */
/*             where V is the matrix stored in the leading N-by-N part of */
/*             the array V. If all eigenvalues of A are real, then the */
/*             leading N-by-N part of this array contains the matrix */
/*             product exp(Lambda*delta) times the inverse of the (right) */
/*             eigenvector matrix of A, where Lambda is the diagonal */
/*             matrix of eigenvalues. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= max(1,N). */

/*     VALR    (output) DOUBLE PRECISION array, dimension (N) */
/*     VALI    (output) DOUBLE PRECISION array, dimension (N) */
/*             These arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the matrix A. The */
/*             eigenvalues are unordered except that complex conjugate */
/*             pairs of values appear consecutively with the eigenvalue */
/*             having positive imaginary part first. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and if N > 0, DWORK(2) returns the reciprocal */
/*             condition number of the triangular matrix used to obtain */
/*             the inverse of the eigenvector matrix. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= max(1,4*N). */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if INFO = i, the QR algorithm failed to compute all */
/*                   the eigenvalues; no eigenvectors have been computed; */
/*                   elements i+1:N of VALR and VALI contain eigenvalues */
/*                   which have converged; */
/*             = N+1:  if the inverse of the eigenvector matrix could not */
/*                   be formed due to an attempt to divide by zero, i.e., */
/*                   the eigenvector matrix is singular; */
/*             = N+2:  if the matrix A is defective, possibly due to */
/*                   rounding errors. */

/*     METHOD */

/*     This routine is an implementation of "Method 15" of the set of */
/*     methods described in reference [1], which uses an eigenvalue/ */
/*     eigenvector decomposition technique. A modification of LAPACK */
/*     Library routine DGEEV is used for obtaining the right eigenvector */
/*     matrix. A condition estimate is then employed to determine if the */
/*     matrix A is near defective and hence the exponential solution is */
/*     inaccurate. In this case the routine returns with the Error */
/*     Indicator (INFO) set to N+2, and SLICOT Library routines MB05ND or */
/*     MB05OD are the preferred alternative routines to be used. */

/*     REFERENCES */

/*     [1] Moler, C.B. and Van Loan, C.F. */
/*         Nineteen dubious ways to compute the exponential of a matrix. */
/*         SIAM Review, 20, pp. 801-836, 1978. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05AD by M.J. Denham, Kingston */
/*     Polytechnic, March 1981. */

/*     REVISIONS */

/*     V. Sima, June 13, 1997, April 25, 2003, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvector decomposition, matrix exponential. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --valr;
    --vali;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    scale = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1);
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || scale)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldv < max(1,*n)) {
	*info = -7;
    } else if (*ldy < max(1,*n)) {
	*info = -9;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n << 2;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -14;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB05MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/*     Compute the eigenvalues and right eigenvectors of the real */
/*     nonsymmetric matrix A; optionally, compute a balancing */
/*     transformation. */
/*     Workspace:  need: 4*N. */

    mb05my_(balanc, n, &a[a_offset], lda, &valr[1], &vali[1], &v[v_offset], 
	    ldv, &y[y_offset], ldy, &dwork[1], ldwork, info, (ftnlen)1);

    if (*info > 0) {
	return 0;
    }
    wrkopt = dwork[1];
    if (scale) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[i__] = dwork[i__ + 1];
/* L10: */
	}
    }

/*     Exit with INFO = N + 1 if V is exactly singular. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (v[i__ + i__ * v_dim1] == 0.) {
	    *info = *n + 1;
	    return 0;
	}
/* L20: */
    }

/*     Compute the reciprocal condition number of the triangular matrix. */

    dtrcon_("1-norm", "Upper", "Non unit", n, &v[v_offset], ldv, &rcond, &
	    dwork[*n + 1], &iwork[1], info, (ftnlen)6, (ftnlen)5, (ftnlen)8);

/*     Return if the matrix is singular to working precision. */

    if (rcond < dlamch_("Epsilon", (ftnlen)7)) {
	dwork[2] = rcond;
	*info = *n + 2;
	return 0;
    }

/*     Compute the right eigenvector matrix (temporarily) in A. */

    dlacpy_("Full", n, n, &y[y_offset], ldy, &a[a_offset], lda, (ftnlen)4);
    dtrmm_("Right", "Upper", "No transpose", "Non unit", n, n, &c_b16, &v[
	    v_offset], ldv, &a[a_offset], lda, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
    if (scale) {
	dgebak_(balanc, "Right", n, &c__1, n, &dwork[1], n, &a[a_offset], lda,
		 info, (ftnlen)1, (ftnlen)5);
    }

/*     Compute the inverse of the right eigenvector matrix, by solving */
/*     a set of linear systems, V * X = Y' (if BALANC = 'N'). */

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	dswap_(&i__2, &y[i__ + y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1);
/* L40: */
    }

    dtrsm_("Left", "Upper", "No transpose", "Non unit", n, n, &c_b16, &v[
	    v_offset], ldv, &y[y_offset], ldy, (ftnlen)4, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
    if (scale) {

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tempr = 1. / dwork[i__];
	    dscal_(n, &tempr, &y[i__ * y_dim1 + 1], &c__1);
/* L60: */
	}

    }

/*     Save the right eigenvector matrix in V. */

    dlacpy_("Full", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)4);

/*     Premultiply the inverse eigenvector matrix by the exponential of */
/*     quasi-diagonal matrix Lambda * DELTA, where Lambda is the matrix */
/*     of eigenvalues. */
/*     Note that only real arithmetic is used, taking the special storing */
/*     of eigenvalues/eigenvectors into account. */

    i__ = 0;
/*     REPEAT */
L80:
    ++i__;
    if (vali[i__] == 0.) {
	tempr = exp(valr[i__] * *delta);
	dscal_(n, &tempr, &y[i__ + y_dim1], ldy);
    } else {
	tempr = valr[i__] * *delta;
	tempi = vali[i__] * *delta;
	tmp[0] = cos(tempi) * exp(tempr);
	tmp[2] = sin(tempi) * exp(tempr);
	tmp[1] = -tmp[2];
	tmp[3] = tmp[0];
	dlacpy_("Full", &c__2, n, &y[i__ + y_dim1], ldy, &dwork[1], &c__2, (
		ftnlen)4);
	dgemm_("No transpose", "No transpose", &c__2, n, &c__2, &c_b16, tmp, &
		c__2, &dwork[1], &c__2, &c_b40, &y[i__ + y_dim1], ldy, (
		ftnlen)12, (ftnlen)12);
	++i__;
    }
    if (i__ < *n) {
	goto L80;
    }
/*     UNTIL I = N. */

/*     Compute the matrix exponential as the product V * Y. */

    dgemm_("No transpose", "No transpose", n, n, n, &c_b16, &v[v_offset], ldv,
	     &y[y_offset], ldy, &c_b40, &a[a_offset], lda, (ftnlen)12, (
	    ftnlen)12);

/*     Set optimal workspace dimension and reciprocal condition number. */

    dwork[1] = wrkopt;
    dwork[2] = rcond;

    return 0;
/* *** Last line of MB05MD *** */
} /* mb05md_ */

