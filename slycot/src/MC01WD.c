/* MC01WD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mc01wd_(integer *dp, doublereal *p, doublereal *u1, 
	doublereal *u2, doublereal *q, integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal a, b, c__;
    static integer i__, n;
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

/*     To compute, for a given real polynomial P(x) and a quadratic */
/*     polynomial B(x), the quotient polynomial Q(x) and the linear */
/*     remainder polynomial R(x) such that */

/*        P(x) = B(x) * Q(x) + R(x), */

/*                                 2 */
/*     where B(x) = u1 + u2 * x + x , R(x) = q(1) + q(2) * (u2 + x) */
/*     and u1, u2, q(1) and q(2) are real scalars. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial P(x).  DP >= 0. */

/*     P       (input) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array must contain the coefficients of P(x) in */
/*             increasing powers of x. */

/*     U1      (input) DOUBLE PRECISION */
/*             The value of the constant term of the quadratic */
/*             polynomial B(x). */

/*     U2      (input) DOUBLE PRECISION */
/*             The value of the coefficient of x of the quadratic */
/*             polynomial B(x). */

/*     Q       (output) DOUBLE PRECISION array, dimension (DP+1) */
/*             If DP >= 1 on entry, then elements Q(1) and Q(2) contain */
/*             the coefficients q(1) and q(2), respectively, of the */
/*             remainder polynomial R(x), and the next (DP-1) elements */
/*             of this array contain the coefficients of the quotient */
/*             polynomial Q(x) in increasing powers of x. */
/*             If DP = 0 on entry, then element Q(1) contains the */
/*             coefficient q(1) of the remainder polynomial R(x) = q(1); */
/*             Q(x) is the zero polynomial. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given the real polynomials */

/*                DP           i                           2 */
/*        P(x) = SUM p(i+1) * x  and B(x) = u1 + u2 * x + x */
/*               i=0 */

/*     the routine uses the recurrence relationships */

/*        q(DP+1) = p(DP+1), */

/*        q(DP) = p(DP) - u2 * q(DP+1) and */

/*        q(i)  = p(i) - u2 * q(i+1) - u1 * q(i+2) for i = DP-1, ..., 1 */

/*     to determine the coefficients of the quotient polynomial */

/*               DP-2          i */
/*        Q(x) = SUM q(i+3) * x */
/*               i=0 */

/*     and the remainder polynomial */

/*        R(x) = q(1) + q(2) * (u2 + x). */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01KD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations, */
/*     quadratic polynomial. */

/*     ****************************************************************** */

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
    if (*dp < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("MC01WD", &i__1, (ftnlen)6);
	return 0;
    }

    *info = 0;
    n = *dp + 1;
    q[n] = p[n];
    if (n > 1) {
	b = q[n];
	q[n - 1] = p[n - 1] - *u2 * b;
	if (n > 2) {
	    a = q[n - 1];

	    for (i__ = n - 2; i__ >= 1; --i__) {
		c__ = p[i__] - *u2 * a - *u1 * b;
		q[i__] = c__;
		b = a;
		a = c__;
/* L20: */
	    }

	}
    }

    return 0;
/* *** Last line of MC01WD *** */
} /* mc01wd_ */

