/* MB01ND.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb01nd_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To perform the skew-symmetric rank 2 operation */

/*          A := alpha*x*y' - alpha*y*x' + A, */

/*     where alpha is a scalar, x and y are vectors of length n and A is */
/*     an n-by-n skew-symmetric matrix. */

/*     This is a modified version of the vanilla implemented BLAS */
/*     routine DSYR2 written by Jack Dongarra, Jeremy Du Croz, */
/*     Sven Hammarling, and Richard Hanson. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the upper or lower triangular part of */
/*             the array A is to be referenced as follows: */
/*             = 'U':  only the strictly upper triangular part of A is to */
/*                     be referenced; */
/*             = 'L':  only the strictly lower triangular part of A is to */
/*                     be referenced. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. If alpha is zero X and Y are not */
/*             referenced. */

/*     X       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCX ) ). */
/*             On entry, elements 1, INCX+1, .., ( N - 1 )*INCX + 1 of */
/*             this array must contain the elements of the vector X. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X. IF INCX < 0 then the */
/*             elements of X are accessed in reversed order.  INCX <> 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCY ) ). */
/*             On entry, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array must contain the elements of the vector Y. */

/*     INCY    (input) INTEGER */
/*             The increment for the elements of Y. IF INCY < 0 then the */
/*             elements of Y are accessed in reversed order.  INCY <> 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry with UPLO = 'U', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the matrix A. The lower triangular part of this array is */
/*             not referenced. */
/*             On entry with UPLO = 'L', the leading N-by-N part of this */
/*             array must contain the strictly lower triangular part of */
/*             the matrix A. The upper triangular part of this array is */
/*             not referenced. */
/*             On exit with UPLO = 'U', the leading N-by-N part of this */
/*             array contains the strictly upper triangular part of the */
/*             updated matrix A. */
/*             On exit with UPLO = 'L', the leading N-by-N part of this */
/*             array contains the strictly lower triangular part of the */
/*             updated matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N) */

/*     NUMERICAL ASPECTS */

/*     Though being almost identical with the vanilla implementation */
/*     of the BLAS routine DSYR2 the performance of this routine could */
/*     be significantly lower in the case of vendor supplied, highly */
/*     optimized BLAS. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DSKR2). */

/*     KEYWORDS */

/*     Elementary matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }

    if (info != 0) {
	xerbla_("MB01ND", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form A when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 - y[i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 - y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form A when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 - y[i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 - y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }
    return 0;
/* *** Last line of MB01ND *** */
} /* mb01nd_ */

