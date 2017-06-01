/* MB02UU.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02uu_(integer *n, doublereal *a, integer *lda, 
	doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal eps, temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal factor, bignum, smlnum;


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

/*     To solve for x in A * x = scale * RHS, using the LU factorization */
/*     of the N-by-N matrix A computed by SLICOT Library routine MB02UV. */
/*     The factorization has the form A = P * L * U * Q, where P and Q */
/*     are permutation matrices, L is unit lower triangular and U is */
/*     upper triangular. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA, N) */
/*             The leading N-by-N part of this array must contain */
/*             the LU part of the factorization of the matrix A computed */
/*             by SLICOT Library routine MB02UV:  A = P * L * U * Q. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1, N). */

/*     RHS     (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the right hand side */
/*             of the system. */
/*             On exit, this array contains the solution of the system. */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the */
/*             matrix has been interchanged with row IPIV(i). */

/*     JPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= j <= N, column j of the */
/*             matrix has been interchanged with column JPIV(j). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, chosen 0 < SCALE <= 1 to prevent */
/*             overflow in the solution. */

/*     FURTHER COMMENTS */

/*     In the interest of speed, this routine does not check the input */
/*     for errors. It should only be used if the order of the matrix A */
/*     is very small. */

/*     CONTRIBUTOR */

/*     Bo Kagstrom and P. Poromaa, Univ. of Umea, Sweden, Nov. 1993. */

/*     REVISIONS */

/*     April 1998 (T. Penzl). */
/*     Sep. 1998 (V. Sima). */

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
    --rhs;
    --ipiv;
    --jpiv;

    /* Function Body */
    eps = dlamch_("Precision", (ftnlen)9);
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

/*     Apply permutations IPIV to RHS. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip = ipiv[i__];
	if (ip != i__) {
	    temp = rhs[i__];
	    rhs[i__] = rhs[ip];
	    rhs[ip] = temp;
	}
/* L20: */
    }

/*     Solve for L part. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__;
	d__1 = -rhs[i__];
	daxpy_(&i__2, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &rhs[i__ + 1],
		 &c__1);
/* L40: */
    }

/*     Solve for U part. */

/*     Check for scaling. */

    factor = (doublereal) (*n) * 2.;
    i__ = 1;
L60:
    if (factor * smlnum * (d__1 = rhs[i__], abs(d__1)) <= (d__2 = a[i__ + i__ 
	    * a_dim1], abs(d__2))) {
	++i__;
	if (i__ <= *n) {
	    goto L60;
	}
	*scale = 1.;
    } else {
	*scale = 1. / factor / (d__1 = rhs[idamax_(n, &rhs[1], &c__1)], abs(
		d__1));
	dscal_(n, scale, &rhs[1], &c__1);
    }

    for (i__ = *n; i__ >= 1; --i__) {
	temp = 1. / a[i__ + i__ * a_dim1];
	rhs[i__] *= temp;
	i__1 = *n;
	for (j = i__ + 1; j <= i__1; ++j) {
	    rhs[i__] -= rhs[j] * (a[i__ + j * a_dim1] * temp);
/* L80: */
	}
/* L100: */
    }

/*     Apply permutations JPIV to the solution (RHS). */

    for (i__ = *n - 1; i__ >= 1; --i__) {
	ip = jpiv[i__];
	if (ip != i__) {
	    temp = rhs[i__];
	    rhs[i__] = rhs[ip];
	    rhs[ip] = temp;
	}
/* L120: */
    }

    return 0;
/* *** Last line of MB02UU *** */
} /* mb02uu_ */

