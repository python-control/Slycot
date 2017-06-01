/* MB01MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb01md_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal 
	*beta, doublereal *y, integer *incy, ftnlen uplo_len)
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

/*     To perform the matrix-vector operation */

/*        y := alpha*A*x + beta*y, */

/*     where alpha and beta are scalars, x and y are vectors of length */
/*     n and A is an n-by-n skew-symmetric matrix. */

/*     This is a modified version of the vanilla implemented BLAS */
/*     routine DSYMV written by Jack Dongarra, Jeremy Du Croz, */
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
/*             The scalar alpha. If alpha is zero the array A is not */
/*             referenced. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry with UPLO = 'U', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the matrix A. The lower triangular part of this array is */
/*             not referenced. */
/*             On entry with UPLO = 'L', the leading N-by-N part of this */
/*             array must contain the strictly lower triangular part of */
/*             the matrix A. The upper triangular part of this array is */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N) */

/*     X       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCX ) ). */
/*             On entry, elements 1, INCX+1, .., ( N - 1 )*INCX + 1 of */
/*             this array must contain the elements of the vector X. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X. IF INCX < 0 then the */
/*             elements of X are accessed in reversed order.  INCX <> 0. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. If beta is zero then Y need not be set on */
/*             input. */

/*     Y       (input/output) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCY ) ). */
/*             On entry, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array must contain the elements of the vector Y. */
/*             On exit, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array contain the updated elements of the vector Y. */

/*     INCY    (input) INTEGER */
/*             The increment for the elements of Y. IF INCY < 0 then the */
/*             elements of Y are accessed in reversed order.  INCY <> 0. */

/*     NUMERICAL ASPECTS */

/*     Though being almost identical with the vanilla implementation */
/*     of the BLAS routine DSYMV the performance of this routine could */
/*     be significantly lower in the case of vendor supplied, highly */
/*     optimized BLAS. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DSKMV). */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("MB01MD", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

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

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }

/*     Quick return if possible. */

    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form y when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] -= *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx + *incx;
	    jy = ky + *incy;
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] -= *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form y when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] -= *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] -= *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }
/* *** Last line of MB01MD *** */
    return 0;
} /* mb01md_ */

