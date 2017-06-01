/* MC01RD.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;
static integer c_n1 = -1;

/* Subroutine */ int mc01rd_(integer *dp1, integer *dp2, integer *dp3, 
	doublereal *alpha, doublereal *p1, doublereal *p2, doublereal *p3, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, d1, d2, d3, e3, dmin__, dmax__;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer dsum;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute the coefficients of the polynomial */

/*        P(x) = P1(x) * P2(x) + alpha * P3(x), */

/*     where P1(x), P2(x) and P3(x) are given real polynomials and alpha */
/*     is a real scalar. */

/*     Each of the polynomials P1(x), P2(x) and P3(x) may be the zero */
/*     polynomial. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP1     (input) INTEGER */
/*             The degree of the polynomial P1(x).  DP1 >= -1. */

/*     DP2     (input) INTEGER */
/*             The degree of the polynomial P2(x).  DP2 >= -1. */

/*     DP3     (input/output) INTEGER */
/*             On entry, the degree of the polynomial P3(x).  DP3 >= -1. */
/*             On exit, the degree of the polynomial P(x). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar value alpha of the problem. */

/*     P1      (input) DOUBLE PRECISION array, dimension (lenp1) */
/*             where lenp1 = DP1 + 1 if DP1 >= 0 and 1 otherwise. */
/*             If DP1 >= 0, then this array must contain the */
/*             coefficients of P1(x) in increasing powers of x. */
/*             If DP1 = -1, then P1(x) is taken to be the zero */
/*             polynomial, P1 is not referenced and can be supplied */
/*             as a dummy array. */

/*     P2      (input) DOUBLE PRECISION array, dimension (lenp2) */
/*             where lenp2 = DP2 + 1 if DP2 >= 0 and 1 otherwise. */
/*             If DP2 >= 0, then this array must contain the */
/*             coefficients of P2(x) in increasing powers of x. */
/*             If DP2 = -1, then P2(x) is taken to be the zero */
/*             polynomial, P2 is not referenced and can be supplied */
/*             as a dummy array. */

/*     P3      (input/output) DOUBLE PRECISION array, dimension (lenp3) */
/*             where lenp3 = MAX(DP1+DP2,DP3,0) + 1. */
/*             On entry, if DP3 >= 0, then this array must contain the */
/*             coefficients of P3(x) in increasing powers of x. */
/*             On entry, if DP3 = -1, then P3(x) is taken to be the zero */
/*             polynomial. */
/*             On exit, the leading (DP3+1) elements of this array */
/*             contain the coefficients of P(x) in increasing powers of x */
/*             unless DP3 = -1 on exit, in which case the coefficients of */
/*             P(x) (the zero polynomial) are not stored in the array. */
/*             This is the case, for instance, when ALPHA = 0.0 and */
/*             P1(x) or P2(x) is the zero polynomial. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given real polynomials */

/*                DP1           i           DP2           i */
/*        P1(x) = SUM a(i+1) * x ,  P2(x) = SUM b(i+1) * x  and */
/*                i=0                       i=0 */

/*                DP3           i */
/*        P3(x) = SUM c(i+1) * x , */
/*                i=0 */

/*     the routine computes the coefficents of P(x) = P1(x) * P2(x) + */
/*                     DP3            i */
/*     alpha * P3(x) = SUM  d(i+1) * x  as follows. */
/*                     i=0 */

/*     Let e(i) = c(i) for 1 <= i <= DP3+1 and e(i) = 0 for i > DP3+1. */
/*     Then if DP1 >= DP2, */

/*                i */
/*        d(i) = SUM a(k) * b(i-k+1) + f(i), for i = 1, ..., DP2+1, */
/*               k=1 */

/*                 i */
/*        d(i)  = SUM a(k) * b(i-k+1) + f(i), for i = DP2+2, ..., DP1+1 */
/*               k=i-DP2 */

/*     and */
/*                DP1+1 */
/*        d(i)  = SUM a(k) * b(i-k+1) + f(i) for i = DP1+2,...,DP1+DP2+1, */
/*               k=i-DP2 */

/*     where f(i) = alpha * e(i). */

/*     Similar formulas hold for the case DP1 < DP2. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01FD by C. Klimann and */
/*     A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    --p3;
    --p2;
    --p1;

    /* Function Body */
    *info = 0;
    if (*dp1 < -1) {
	*info = -1;
    } else if (*dp2 < -1) {
	*info = -2;
    } else if (*dp3 < -1) {
	*info = -3;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC01RD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Computation of the exact degree of the polynomials, i.e., Di such */
/*     that either Di = -1 or Pi(Di+1) is non-zero. */

    d1 = *dp1;
/*     WHILE ( D1 >= 0 and P1(D1+1) = 0 ) DO */
L20:
    if (d1 >= 0) {
	if (p1[d1 + 1] == 0.) {
	    --d1;
	    goto L20;
	}
    }
/*     END WHILE 20 */
    d2 = *dp2;
/*     WHILE ( D2 >= 0 and P2(D2+1) = 0 ) DO */
L40:
    if (d2 >= 0) {
	if (p2[d2 + 1] == 0.) {
	    --d2;
	    goto L40;
	}
    }
/*     END WHILE 40 */
    if (*alpha == 0.) {
	d3 = -1;
    } else {
	d3 = *dp3;
    }
/*     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO */
L60:
    if (d3 >= 0) {
	if (p3[d3 + 1] == 0.) {
	    --d3;
	    goto L60;
	}
    }
/*     END WHILE 60 */

/*     Computation of P3(x) := ALPHA * P3(x). */

    i__1 = d3 + 1;
    dscal_(&i__1, alpha, &p3[1], &c__1);

    if (d1 == -1 || d2 == -1) {
	*dp3 = d3;
	return 0;
    }

/*     P1(x) and P2(x) are non-zero polynomials. */

    dsum = d1 + d2;
    dmax__ = max(d1,d2);
    dmin__ = dsum - dmax__;

    if (d3 < dsum) {
	p3[d3 + 2] = 0.;
	i__1 = dsum - d3 - 1;
	dcopy_(&i__1, &p3[d3 + 2], &c__0, &p3[d3 + 3], &c__1);
	d3 = dsum;
    }

    if (d1 == 0 || d2 == 0) {

/*        D1 or D2 is zero. */

	if (d1 != 0) {
	    i__1 = d1 + 1;
	    daxpy_(&i__1, &p2[1], &p1[1], &c__1, &p3[1], &c__1);
	} else {
	    i__1 = d2 + 1;
	    daxpy_(&i__1, &p1[1], &p2[1], &c__1, &p3[1], &c__1);
	}
    } else {

/*        D1 and D2 are both nonzero. */

/*        First part of the computation. */

	i__1 = dmin__ + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p3[i__] += ddot_(&i__, &p1[1], &c__1, &p2[1], &c_n1);
/* L80: */
	}

/*        Second part of the computation. */

	i__1 = dmax__ + 1;
	for (i__ = dmin__ + 2; i__ <= i__1; ++i__) {
	    if (d1 > d2) {
		k = i__ - d2;
		i__2 = dmin__ + 1;
		p3[i__] += ddot_(&i__2, &p1[k], &c__1, &p2[1], &c_n1);
	    } else {
		k = i__ - d1;
		i__2 = dmin__ + 1;
		p3[i__] += ddot_(&i__2, &p2[k], &c_n1, &p1[1], &c__1);
	    }
/* L100: */
	}

/*        Third part of the computation. */

	e3 = dsum + 2;

	i__1 = dsum + 1;
	for (i__ = dmax__ + 2; i__ <= i__1; ++i__) {
	    j = e3 - i__;
	    k = i__ - dmin__;
	    l = i__ - dmax__;
	    if (d1 > d2) {
		p3[i__] += ddot_(&j, &p1[k], &c__1, &p2[l], &c_n1);
	    } else {
		p3[i__] += ddot_(&j, &p1[l], &c_n1, &p2[k], &c__1);
	    }
/* L120: */
	}

    }

/*     Computation of the exact degree of P3(x). */

/*     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO */
L140:
    if (d3 >= 0) {
	if (p3[d3 + 1] == 0.) {
	    --d3;
	    goto L140;
	}
    }
/*     END WHILE 140 */
    *dp3 = d3;

    return 0;
/* *** Last line of MC01RD *** */
} /* mc01rd_ */

