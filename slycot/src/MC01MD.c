/* MC01MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mc01md_(integer *dp, doublereal *alpha, integer *k, 
	doublereal *p, doublereal *q, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To calculate, for a given real polynomial P(x) and a real scalar */
/*     alpha, the leading K coefficients of the shifted polynomial */
/*                                                               K-1 */
/*        P(x) = q(1) + q(2) * (x-alpha) + ... + q(K) * (x-alpha)   + ... */

/*     using Horner's algorithm. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial P(x).  DP >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar value alpha of the problem. */

/*     K       (input) INTEGER */
/*             The number of coefficients of the shifted polynomial to be */
/*             computed.  1 <= K <= DP+1. */

/*     P       (input) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array must contain the coefficients of P(x) in */
/*             increasing powers of x. */

/*     Q       (output) DOUBLE PRECISION array, dimension (DP+1) */
/*             The leading K elements of this array contain the first */
/*             K coefficients of the shifted polynomial in increasing */
/*             powers of (x - alpha), and the next (DP-K+1) elements */
/*             are used as internal workspace. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given the real polynomial */
/*                                         2                    DP */
/*        P(x) = p(1) + p(2) * x + p(3) * x  + ... + p(DP+1) * x  , */

/*     the routine computes the leading K coefficients of the shifted */
/*     polynomial */
/*                                                                   K-1 */
/*        P(x) = q(1) + q(2) * (x - alpha) + ... + q(K) * (x - alpha) */

/*     as follows. */

/*     Applying Horner's algorithm (see [1]) to P(x), i.e. dividing P(x) */
/*     by (x-alpha), yields */

/*        P(x) = q(1) + (x-alpha) * D(x), */

/*     where q(1) is the value of the constant term of the shifted */
/*     polynomial and D(x) is the quotient polynomial of degree (DP-1) */
/*     given by */
/*                                         2                     DP-1 */
/*        D(x) = d(2) + d(3) * x + d(4) * x  + ... +  d(DP+1) * x    . */

/*     Applying Horner's algorithm to D(x) and subsequent quotient */
/*     polynomials yields q(2) and q(3), q(4), ..., q(K) respectively. */

/*     It follows immediately that q(1) = P(alpha), and in general */
/*                (i-1) */
/*        q(i) = P     (alpha) / (i - 1)! for i = 1, 2, ..., K. */

/*     REFERENCES */

/*     [1] STOER, J. and BULIRSCH, R. */
/*         Introduction to Numerical Analysis. */
/*         Springer-Verlag. 1980. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01AD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    --q;
    --p;

    /* Function Body */
    *info = 0;
    if (*dp < 0) {
	*info = -1;
    } else if (*k <= 0 || *k > *dp + 1) {
	*info = -3;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC01MD", &i__1, (ftnlen)6);
	return 0;
    }

    i__1 = *dp + 1;
    dcopy_(&i__1, &p[1], &c__1, &q[1], &c__1);
    if (*dp == 0 || *alpha == 0.) {
	return 0;
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {

	i__2 = j;
	for (i__ = *dp; i__ >= i__2; --i__) {
	    q[i__] += *alpha * q[i__ + 1];
/* L20: */
	}

/* L40: */
    }

    return 0;
/* *** Last line of MC01MD *** */
} /* mc01md_ */

