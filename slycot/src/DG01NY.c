/* DG01NY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int dg01ny_(char *indi, integer *n, doublereal *xr, 
	doublereal *xi, ftnlen indi_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, n2;
    static doublereal ai, bi, ar, br, wi, wr, pi2;
    static logical lindi;
    static doublereal helpi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal helpr, whelp, wstpi, wstpr;


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

/*     For efficiency, no tests of the input scalar parameters are */
/*     performed. */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --xi;
    --xr;

    /* Function Body */
    lindi = lsame_(indi, "D", (ftnlen)1, (ftnlen)1);

/*     Initialisation. */

    pi2 = atan(1.) * 8.;
    if (lindi) {
	pi2 = -pi2;
    }

    whelp = pi2 / (doublereal) (*n << 1);
    wstpi = sin(whelp);
    whelp = sin(whelp * .5);
    wstpr = whelp * -2. * whelp;
    wi = 0.;

    if (lindi) {
	wr = 1.;
	xr[*n + 1] = xr[1];
	xi[*n + 1] = xi[1];
    } else {
	wr = -1.;
    }

/*     Recursion. */

    n2 = *n / 2 + 1;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n + 2 - i__;
	ar = xr[i__] + xr[j];
	ai = xi[i__] - xi[j];
	br = xi[i__] + xi[j];
	bi = xr[j] - xr[i__];
	if (lindi) {
	    ar *= .5;
	    ai *= .5;
	    br *= .5;
	    bi *= .5;
	}
	helpr = wr * br - wi * bi;
	helpi = wr * bi + wi * br;
	xr[i__] = ar + helpr;
	xi[i__] = ai + helpi;
	xr[j] = ar - helpr;
	xi[j] = helpi - ai;
	whelp = wr;
	wr = wr + wr * wstpr - wi * wstpi;
	wi = wi + wi * wstpr + whelp * wstpi;
/* L10: */
    }

    return 0;
/* *** Last line of DG01NY *** */
} /* dg01ny_ */

