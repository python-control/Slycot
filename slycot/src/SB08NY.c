/* SB08NY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb08ny_(integer *da, doublereal *a, doublereal *b, 
	doublereal *epsb)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dlamch_(char *, ftnlen);


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

/*     To compute the coefficients of B(z) = A(1/z) * A(z) and a norm for */
/*     the accuracy of the computed coefficients. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(z) and B(z).  DA >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (DA+1) */
/*             This array must contain the coefficients of the polynomial */
/*             A(z) in increasing powers of z. */

/*     B       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             This array contains the coefficients of the polynomial */
/*             B(z). */

/*     EPSB    (output) DOUBLE PRECISION */
/*             A value used for checking the accuracy of the computed */
/*             coefficients. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08BZ by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Laplace transform, polynomial operations, spectral factorization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    i__1 = *da + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *da - i__ + 2;
	b[i__] = ddot_(&i__2, &a[1], &c__1, &a[i__], &c__1);
/* L20: */
    }

    *epsb = dlamch_("Epsilon", (ftnlen)7) * 3. * b[1];

    return 0;
/* *** Last line of SB08NY *** */
} /* sb08ny_ */

