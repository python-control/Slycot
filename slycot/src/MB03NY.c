/* MB03NY.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;

doublereal mb03ny_(integer *n, doublereal *omega, doublereal *a, integer *lda,
	 doublereal *s, doublereal *dwork, integer *ldwork, doublecomplex *
	cwork, integer *lcwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, ic;
    static doublereal dummy[1]	/* was [1][1] */;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), zgesvd_(char 
	    *, char *, integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublecomplex zdummy[1]	/* was [1][1] */;


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

/*     To compute the smallest singular value of A - jwI. */

/*     FUNCTION VALUE */

/*     MB03NY  DOUBLE PRECISION */
/*             The smallest singular value of A - jwI (if INFO = 0). */
/*             If N = 0, the function value is set to zero. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the the matrix A.  N >= 0. */

/*     OMEGA   (input) DOUBLE PRECISION */
/*             The constant factor of A - jwI. */

/*     A       (input/workspace) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, if OMEGA = 0, the contents of this array are */
/*             destroyed. Otherwise, this array is unchanged. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     S       (output) DOUBLE PRECISION array, dimension (N) */
/*             The singular values of A - jwI in decreasing order. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX( 1, 5*N ). */
/*             For optimum performance LDWORK should be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0 and OMEGA <> 0, CWORK(1) returns the */
/*             optimal value of LCWORK. */
/*             If OMEGA is zero, this array is not referenced. */

/*     LCWORK  INTEGER */
/*             The length of the array CWORK. */
/*             LCWORK >= 1,                 if OMEGA =  0; */
/*             LCWORK >= MAX( 1, N*N+3*N ), if OMEGA <> 0. */
/*             For optimum performance LCWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  The SVD algorithm (in either LAPACK Library routine */
/*                   DGESVD or ZGESVD) fails to converge; this error is */
/*                   very rare. */

/*     METHOD */

/*     This procedure simply constructs the matrix A - jwI, and calls */
/*     ZGESVD if w is not zero, or DGESVD if w = 0. */

/*     FURTHER COMMENTS */

/*     This routine is not very efficient because it computes all */
/*     singular values, but it is very accurate. The routine is intended */
/*     to be called only from the SLICOT Library routine AB13FD. */

/*     CONTRIBUTOR */

/*     R. Byers, the routine SIGMIN (January, 1995). */

/*     REVISIONS */

/*     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999. */

/*     REVISIONS */

/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Apr. 2002, V. Sima. */

/*     KEYWORDS */

/*     singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --dwork;
    --cwork;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 5;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -7;
	} else if (*lcwork < 1 || *omega != 0. && *lcwork < *n * *n + *n * 3) 
		{
	    *info = -9;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB03NY", &i__1, (ftnlen)6);
	return ret_val;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	ret_val = 0.;
	dwork[1] = 1.;
	if (*omega != 0.) {
	    cwork[1].r = 1., cwork[1].i = 0.;
	}
	return ret_val;
    }

    if (*omega == 0.) {

/*        OMEGA = 0 allows real SVD. */

	dgesvd_("No vectors", "No vectors", n, n, &a[a_offset], n, &s[1], 
		dummy, &c__1, dummy, &c__1, &dwork[1], ldwork, info, (ftnlen)
		10, (ftnlen)10);
	if (*info != 0) {
	    *info = 2;
	    return ret_val;
	}
    } else {

/*        General case, that is complex SVD. */

	ic = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = ic;
		i__4 = i__ + j * a_dim1;
		cwork[i__3].r = a[i__4], cwork[i__3].i = 0.;
		++ic;
/* L10: */
	    }
	    i__2 = (j - 1) * *n + j;
	    i__3 = (j - 1) * *n + j;
	    z__2.r = *omega * 0., z__2.i = *omega * 1.;
	    z__1.r = cwork[i__3].r - z__2.r, z__1.i = cwork[i__3].i - z__2.i;
	    cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
/* L20: */
	}
	i__1 = *lcwork - *n * *n;
	zgesvd_("No vectors", "No vectors", n, n, &cwork[1], n, &s[1], zdummy,
		 &c__1, zdummy, &c__1, &cwork[*n * *n + 1], &i__1, &dwork[1], 
		info, (ftnlen)10, (ftnlen)10);
	if (*info != 0) {
	    *info = 2;
	    return ret_val;
	}
	i__1 = *n * *n + 1;
	d__1 = (doublereal) (*n * *n);
	z__2.r = d__1 * 1., z__2.i = d__1 * 0.;
	z__1.r = cwork[i__1].r + z__2.r, z__1.i = cwork[i__1].i + z__2.i;
	cwork[1].r = z__1.r, cwork[1].i = z__1.i;
	dwork[1] = (doublereal) (*n * 5);
    }

    ret_val = s[*n];

/* *** Last line of MB03NY *** */
    return ret_val;
} /* mb03ny_ */

