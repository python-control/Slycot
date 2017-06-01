/* SB08MY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb08my_(integer *da, doublereal *a, doublereal *b, 
	doublereal *epsb)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, k;
    static doublereal sa, sabs, term, maxsa, signi, signk;


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

/*     To compute the coefficients of B(s) = A(s) * A(-s) and a norm */
/*     for the accuracy of the computed coefficients. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(s) and B(s).  DA >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (DA+1) */
/*             This array must contain the coefficients of the polynomial */
/*             A(s) in increasing powers of s. */

/*     B       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             This array contains the coefficients of the polynomial */
/*             B(s) in increasing powers of s**2. */

/*     EPSB    (input/output) DOUBLE PRECISION */
/*             On entry, EPSB must contain the machine precision (see */
/*             LAPACK Library routine DLAMCH). */
/*             On exit, EPSB contains an updated value, using a norm */
/*             for the accuracy of the computed coefficients. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08AZ by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Laplace transform, polynomial operations, spectral factorization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    signi = 1.;
    maxsa = 0.;

    i__1 = *da;
    for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = a[i__ + 1];
	sabs = d__1 * d__1;
	sa = signi * sabs;
	signk = signi * -2.;

/* Computing MIN */
	i__3 = i__, i__4 = *da - i__;
	i__2 = min(i__3,i__4);
	for (k = 1; k <= i__2; ++k) {
	    term = signk * a[i__ - k + 1] * a[i__ + k + 1];
	    sa += term;
	    sabs += abs(term);
	    signk = -signk;
/* L20: */
	}

	b[i__ + 1] = sa;
	maxsa = max(maxsa,sabs);
	signi = -signi;
/* L40: */
    }

    *epsb = maxsa * 3. * *epsb;

    return 0;
/* *** Last line of SB08MY *** */
} /* sb08my_ */

