/* MB03MY.f -- translated by f2c (version 20100827).
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

doublereal mb03my_(integer *nx, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal dx;


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

/*     To compute the absolute minimal value of NX elements in an array. */
/*     The function returns the value zero if NX < 1. */

/*     ARGUMENTS */

/*     NX      (input) INTEGER */
/*             The number of elements in X to be examined. */

/*     X       (input) DOUBLE PRECISION array, dimension (NX * INCX) */
/*             The one-dimensional array of which the absolute minimal */
/*             value of the elements is to be computed. */
/*             This array is not referenced if NX < 1. */

/*     INCX    (input) INTEGER */
/*             The increment to be taken in the array X, defining the */
/*             distance between two consecutive elements.  INCX >= 1. */
/*             INCX = 1, if all elements are contiguous in memory. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MB03AZ by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     June 16, 1997. */

/*     KEYWORDS */

/*     None. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*nx <= 0) {
	ret_val = 0.;
	return ret_val;
    }

    ret_val = abs(x[1]);

    i__1 = *nx * *incx;
    i__2 = *incx;
    for (i__ = *incx + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx = (d__1 = x[i__], abs(d__1));
	if (dx < ret_val) {
	    ret_val = dx;
	}
/* L20: */
    }

    return ret_val;
/* *** Last line of MB03MY *** */
} /* mb03my_ */

