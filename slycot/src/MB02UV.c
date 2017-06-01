/* MB02UV.f -- translated by f2c (version 20100827).
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
static doublereal c_b9 = -1.;

/* Subroutine */ int mb02uv_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, integer *jpiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, ip, jp;
    static doublereal eps;
    static integer ipv, jpv;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal smin, xmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bignum, smlnum;


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

/*     To compute an LU factorization, using complete pivoting, of the */
/*     N-by-N matrix A. The factorization has the form A = P * L * U * Q, */
/*     where P and Q are permutation matrices, L is lower triangular with */
/*     unit diagonal elements and U is upper triangular. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A to be factored. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the factors L and U from the factorization A = P*L*U*Q; */
/*             the unit diagonal elements of L are not stored. If U(k, k) */
/*             appears to be less than SMIN, U(k, k) is given the value */
/*             of SMIN, giving a nonsingular perturbed system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1, N). */

/*     IPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the */
/*             matrix has been interchanged with row IPIV(i). */

/*     JPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= j <= N, column j of the */
/*             matrix has been interchanged with column JPIV(j). */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = k:  U(k, k) is likely to produce owerflow if one tries */
/*                   to solve for x in Ax = b. So U is perturbed to get */
/*                   a nonsingular system. This is a warning. */

/*     FURTHER COMMENTS */

/*     In the interests of speed, this routine does not check the input */
/*     for errors. It should only be used to factorize matrices A of */
/*     very small order. */

/*     CONTRIBUTOR */

/*     Bo Kagstrom and Peter Poromaa, Univ. of Umea, Sweden, Nov. 1993. */

/*     REVISIONS */

/*     April 1998 (T. Penzl). */
/*     Sep. 1998 (V. Sima). */
/*     March 1999 (V. Sima). */
/*     March 2004 (V. Sima). */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Set constants to control owerflow. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --jpiv;

    /* Function Body */
    *info = 0;
    eps = dlamch_("Precision", (ftnlen)9);
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

/*     Find max element in matrix A. */

    ipv = 1;
    jpv = 1;
    xmax = 0.;
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	i__2 = *n;
	for (ip = 1; ip <= i__2; ++ip) {
	    if ((d__1 = a[ip + jp * a_dim1], abs(d__1)) > xmax) {
		xmax = (d__1 = a[ip + jp * a_dim1], abs(d__1));
		ipv = ip;
		jpv = jp;
	    }
/* L20: */
	}
/* L40: */
    }
/* Computing MAX */
    d__1 = eps * xmax;
    smin = max(d__1,smlnum);

/*     Swap rows. */

    if (ipv != 1) {
	dswap_(n, &a[ipv + a_dim1], lda, &a[a_dim1 + 1], lda);
    }
    ipiv[1] = ipv;

/*     Swap columns. */

    if (jpv != 1) {
	dswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[a_dim1 + 1], &c__1);
    }
    jpiv[1] = jpv;

/*     Check for singularity. */

    if ((d__1 = a[a_dim1 + 1], abs(d__1)) < smin) {
	*info = 1;
	a[a_dim1 + 1] = smin;
    }
    if (*n > 1) {
	i__1 = *n - 1;
	d__1 = 1. / a[a_dim1 + 1];
	dscal_(&i__1, &d__1, &a[a_dim1 + 2], &c__1);
	i__1 = *n - 1;
	i__2 = *n - 1;
	dger_(&i__1, &i__2, &c_b9, &a[a_dim1 + 2], &c__1, &a[(a_dim1 << 1) + 
		1], lda, &a[(a_dim1 << 1) + 2], lda);
    }

/*     Factorize the rest of A with complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        Find max element in remaining matrix. */

	ipv = i__;
	jpv = i__;
	xmax = 0.;
	i__2 = *n;
	for (jp = i__; jp <= i__2; ++jp) {
	    i__3 = *n;
	    for (ip = i__; ip <= i__3; ++ip) {
		if ((d__1 = a[ip + jp * a_dim1], abs(d__1)) > xmax) {
		    xmax = (d__1 = a[ip + jp * a_dim1], abs(d__1));
		    ipv = ip;
		    jpv = jp;
		}
/* L60: */
	    }
/* L80: */
	}

/*        Swap rows. */

	if (ipv != i__) {
	    dswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
	}
	ipiv[i__] = ipv;

/*        Swap columns. */

	if (jpv != i__) {
	    dswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
	}
	jpiv[i__] = jpv;

/*        Check for almost singularity. */

	if ((d__1 = a[i__ + i__ * a_dim1], abs(d__1)) < smin) {
	    *info = i__;
	    a[i__ + i__ * a_dim1] = smin;
	}
	i__2 = *n - i__;
	d__1 = 1. / a[i__ + i__ * a_dim1];
	dscal_(&i__2, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
	i__2 = *n - i__;
	i__3 = *n - i__;
	dger_(&i__2, &i__3, &c_b9, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[i__ 
		+ (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * a_dim1], 
		lda);
/* L100: */
    }
    if ((d__1 = a[*n + *n * a_dim1], abs(d__1)) < smin) {
	*info = *n;
	a[*n + *n * a_dim1] = smin;
    }

    return 0;
/* *** Last line of MB02UV *** */
} /* mb02uv_ */

