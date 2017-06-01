/* MB04TU.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04tu_(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
    static doublereal dtemp;


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

/*     To perform the Givens transformation, defined by C (cos) and S */
/*     (sin), and interchange the vectors involved, i.e. */

/*        |X(i)|    | 0   1 |   | C   S |   |X(i)| */
/*        |    | := |       | x |       | x |    |, i = 1,...N. */
/*        |Y(i)|    | 1   0 |   |-S   C |   |Y(i)| */

/*     REMARK. This routine is a modification of DROT from BLAS. */
/*             This routine is called only by the SLICOT routines MB04TX */
/*             and MB04VX. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB04FU by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     January 26, 1998. */

/*     KEYWORDS */

/*     Othogonal transformation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx != 1 || *incy != 1) {

/*        Code for unequal increments or equal increments not equal to 1. */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    iy = (-(*n) + 1) * *incy + 1;
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dtemp = *c__ * y[iy] - *s * x[ix];
	    y[iy] = *c__ * x[ix] + *s * y[iy];
	    x[ix] = dtemp;
	    ix += *incx;
	    iy += *incy;
/* L20: */
	}

    } else {

/*        Code for both increments equal to 1. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dtemp = *c__ * y[i__] - *s * x[i__];
	    y[i__] = *c__ * x[i__] + *s * y[i__];
	    x[i__] = dtemp;
/* L40: */
	}

    }

    return 0;
/* *** Last line of MB04TU *** */
} /* mb04tu_ */

