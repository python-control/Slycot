/* MB01YD.f -- translated by f2c (version 20100827).
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

static doublereal c_b8 = 0.;
static integer c__0 = 0;
static doublereal c_b12 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01yd_(char *uplo, char *trans, integer *n, integer *k, 
	integer *l, doublereal *alpha, doublereal *beta, doublereal *a, 
	integer *lda, doublereal *c__, integer *ldc, integer *info, ftnlen 
	uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, m;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical transp;


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

/*     To perform the symmetric rank k operations */

/*        C := alpha*op( A )*op( A )' + beta*C, */

/*     where alpha and beta are scalars, C is an n-by-n symmetric matrix, */
/*     op( A ) is an n-by-k matrix, and op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A'. */

/*     The matrix A has l nonzero codiagonals, either upper or lower. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the symmetric matrix C */
/*             is given and computed, as follows: */
/*             = 'U':  the upper triangular part is given/computed; */
/*             = 'L':  the lower triangular part is given/computed. */
/*             UPLO also defines the pattern of the matrix A (see below). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( A ) to be used, as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix C.  N >= 0. */

/*     K       (input) INTEGER */
/*             The number of columns of the matrix op( A ).  K >= 0. */

/*     L       (input) INTEGER */
/*             If UPLO = 'U', matrix A has L nonzero subdiagonals. */
/*             If UPLO = 'L', matrix A has L nonzero superdiagonals. */
/*             MAX(0,NR-1) >= L >= 0, if UPLO = 'U', */
/*             MAX(0,NC-1) >= L >= 0, if UPLO = 'L', */
/*             where NR and NC are the numbers of rows and columns of the */
/*             matrix A, respectively. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then the array A is */
/*             not referenced. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then the array C need */
/*             not be set before entry. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,NC), where */
/*             NC is K when TRANS = 'N', and is N otherwise. */
/*             If TRANS = 'N', the leading N-by-K part of this array must */
/*             contain the matrix A, otherwise the leading K-by-N part of */
/*             this array must contain the matrix A. */
/*             If UPLO = 'U', only the upper triangular part and the */
/*             first L subdiagonals are referenced, and the remaining */
/*             subdiagonals are assumed to be zero. */
/*             If UPLO = 'L', only the lower triangular part and the */
/*             first L superdiagonals are referenced, and the remaining */
/*             superdiagonals are assumed to be zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,NR), */
/*             where NR = N, if TRANS = 'N', and NR = K, otherwise. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry with UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix C. */
/*             On entry with UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix C. */
/*             On exit, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of */
/*             this array contains the corresponding triangular part of */
/*             the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The calculations are efficiently performed taking the symmetry */
/*     and structure into account. */

/*     FURTHER COMMENTS */

/*     The matrix A may have the following patterns, when n = 7, k = 5, */
/*     and l = 2 are used for illustration: */

/*     UPLO = 'U', TRANS = 'N'         UPLO = 'L', TRANS = 'N' */

/*            [ x x x x x ]                   [ x x x 0 0 ] */
/*            [ x x x x x ]                   [ x x x x 0 ] */
/*            [ x x x x x ]                   [ x x x x x ] */
/*        A = [ 0 x x x x ],              A = [ x x x x x ], */
/*            [ 0 0 x x x ]                   [ x x x x x ] */
/*            [ 0 0 0 x x ]                   [ x x x x x ] */
/*            [ 0 0 0 0 x ]                   [ x x x x x ] */

/*     UPLO = 'U', TRANS = 'T'         UPLO = 'L', TRANS = 'T' */

/*            [ x x x x x x x ]               [ x x x 0 0 0 0 ] */
/*            [ x x x x x x x ]               [ x x x x 0 0 0 ] */
/*        A = [ x x x x x x x ],          A = [ x x x x x 0 0 ]. */
/*            [ 0 x x x x x x ]               [ x x x x x x 0 ] */
/*            [ 0 0 x x x x x ]               [ x x x x x x x ] */

/*     If N = K, the matrix A is upper or lower triangular, for L = 0, */
/*     and upper or lower Hessenberg, for L = 1. */

/*     This routine is a specialization of the BLAS 3 routine DSYRK. */
/*     BLAS 1 calls are used when appropriate, instead of in-line code, */
/*     in order to increase the efficiency. If the matrix A is full, or */
/*     its zero triangle has small order, an optimized DSYRK code could */
/*     be faster than MB01YD. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    transp = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

    if (transp) {
	nrowa = *k;
	ncola = *n;
    } else {
	nrowa = *n;
	ncola = *k;
    }

    if (upper) {
	m = nrowa;
    } else {
	m = ncola;
    }

    if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (transp || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*k < 0) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 0, i__2 = m - 1;
	if (*l < 0 || *l > max(i__1,i__2)) {
	    *info = -5;
	} else if (*lda < max(1,nrowa)) {
	    *info = -9;
	} else if (*ldc < max(1,*n)) {
	    *info = -11;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01YD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return, if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

    if (*alpha == 0.) {
	if (*beta == 0.) {

/*           Special case when both alpha = 0 and beta = 0. */

	    dlaset_(uplo, n, n, &c_b8, &c_b8, &c__[c_offset], ldc, (ftnlen)1);
	} else {

/*           Special case alpha = 0. */

	    dlascl_(uplo, &c__0, &c__0, &c_b12, beta, n, n, &c__[c_offset], 
		    ldc, info, (ftnlen)1);
	}
	return 0;
    }

/*     General case: alpha <> 0. */

    if (! transp) {

/*        Form  C := alpha*A*A' + beta*C. */

	if (upper) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {

		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L10: */
		    }

		} else if (*beta != 1.) {
		    dscal_(&j, beta, &c__[j * c_dim1 + 1], &c__1);
		}

/* Computing MAX */
		i__2 = 1, i__3 = j - *l;
		i__4 = *k;
		for (m = max(i__2,i__3); m <= i__4; ++m) {
/* Computing MIN */
		    i__3 = j, i__5 = *l + m;
		    i__2 = min(i__3,i__5);
		    d__1 = *alpha * a[j + m * a_dim1];
		    daxpy_(&i__2, &d__1, &a[m * a_dim1 + 1], &c__1, &c__[j * 
			    c_dim1 + 1], &c__1);
/* L20: */
		}

/* L30: */
	    }

	} else {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {

		    i__4 = *n;
		    for (i__ = j; i__ <= i__4; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L40: */
		    }

		} else if (*beta != 1.) {
		    i__4 = *n - j + 1;
		    dscal_(&i__4, beta, &c__[j + j * c_dim1], &c__1);
		}

/* Computing MIN */
		i__2 = j + *l;
		i__4 = min(i__2,*k);
		for (m = 1; m <= i__4; ++m) {
		    i__2 = *n - j + 1;
		    d__1 = *alpha * a[j + m * a_dim1];
		    daxpy_(&i__2, &d__1, &a[j + m * a_dim1], &c__1, &c__[j + 
			    j * c_dim1], &c__1);
/* L50: */
		}

/* L60: */
	    }

	}

    } else {

/*        Form  C := alpha*A'*A + beta*C. */

	if (upper) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__4 = j;
		for (i__ = 1; i__ <= i__4; ++i__) {
/* Computing MIN */
		    i__3 = j + *l;
		    i__2 = min(i__3,*k);
		    temp = *alpha * ddot_(&i__2, &a[i__ * a_dim1 + 1], &c__1, 
			    &a[j * a_dim1 + 1], &c__1);
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = temp;
		    } else {
			c__[i__ + j * c_dim1] = temp + *beta * c__[i__ + j * 
				c_dim1];
		    }
/* L70: */
		}

/* L80: */
	    }

	} else {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__4 = *n;
		for (i__ = j; i__ <= i__4; ++i__) {
/* Computing MAX */
		    i__2 = 1, i__3 = i__ - *l;
		    m = max(i__2,i__3);
		    i__2 = *k - m + 1;
		    temp = *alpha * ddot_(&i__2, &a[m + i__ * a_dim1], &c__1, 
			    &a[m + j * a_dim1], &c__1);
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = temp;
		    } else {
			c__[i__ + j * c_dim1] = temp + *beta * c__[i__ + j * 
				c_dim1];
		    }
/* L90: */
		}

/* L100: */
	    }

	}

    }

    return 0;

/* *** Last line of MB01YD *** */
} /* mb01yd_ */

