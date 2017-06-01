/* SB01BX.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb01bx_(logical *reig, integer *n, doublereal *xr, 
	doublereal *xi, doublereal *wr, doublereal *wi, doublereal *s, 
	doublereal *p)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal x, y;


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

/*     To choose a real eigenvalue or a pair of complex conjugate */
/*     eigenvalues at "minimal" distance to a given real or complex */
/*     value. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     REIG    LOGICAL */
/*             Specifies the type of eigenvalues as follows: */
/*             = .TRUE.,  a real eigenvalue is to be selected; */
/*             = .FALSE., a pair of complex eigenvalues is to be */
/*                        selected. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of eigenvalues contained in the arrays WR */
/*             and WI.  N >= 1. */

/*     XR,XI   (input) DOUBLE PRECISION */
/*             If REIG = .TRUE., XR must contain the real value and XI */
/*             is assumed zero and therefore not referenced. */
/*             If REIG = .FALSE., XR must contain the real part and XI */
/*             the imaginary part, respectively, of the complex value. */

/*     WR,WI   (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if REIG = .TRUE., WR must contain the real */
/*             eigenvalues from which an eigenvalue at minimal distance */
/*             to XR is to be selected. In this case, WI is considered */
/*             zero and therefore not referenced. */
/*             On entry, if REIG = .FALSE., WR and WI must contain the */
/*             real and imaginary parts, respectively, of the eigenvalues */
/*             from which a pair of complex conjugate eigenvalues at */
/*             minimal "distance" to XR + jXI is to be selected. */
/*             The eigenvalues of each pair of complex conjugate */
/*             eigenvalues must appear consecutively. */
/*             On exit, the elements of these arrays are reordered such */
/*             that the selected eigenvalue(s) is (are) found in the */
/*             last element(s) of these arrays. */

/*     S,P     (output) DOUBLE PRECISION */
/*             If REIG = .TRUE., S (and also P) contains the value of */
/*             the selected real eigenvalue. */
/*             If REIG = .FALSE., S and P contain the sum and product, */
/*             respectively, of the selected complex conjugate pair of */
/*             eigenvalues. */

/*     FURTHER COMMENTS */

/*     For efficiency reasons, |x| + |y| is used for a complex number */
/*     x + jy, instead of its modulus. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     February 1999. Based on the RASP routine PMDIST. */

/*     REVISIONS */

/*     March 30, 1999, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     Feb. 15, 2004, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --wi;
    --wr;

    /* Function Body */
    j = 1;
    if (*reig) {
	y = (d__1 = wr[1] - *xr, abs(d__1));
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    x = (d__1 = wr[i__] - *xr, abs(d__1));
	    if (x < y) {
		y = x;
		j = i__;
	    }
/* L10: */
	}
	*s = wr[j];
	k = *n - j;
	if (k > 0) {
	    i__1 = j + k - 1;
	    for (i__ = j; i__ <= i__1; ++i__) {
		wr[i__] = wr[i__ + 1];
/* L20: */
	    }
	    wr[*n] = *s;
	}
	*p = *s;
    } else {
	y = (d__1 = wr[1] - *xr, abs(d__1)) + (d__2 = wi[1] - *xi, abs(d__2));
	i__1 = *n;
	for (i__ = 3; i__ <= i__1; i__ += 2) {
	    x = (d__1 = wr[i__] - *xr, abs(d__1)) + (d__2 = wi[i__] - *xi, 
		    abs(d__2));
	    if (x < y) {
		y = x;
		j = i__;
	    }
/* L30: */
	}
	x = wr[j];
	y = wi[j];
	k = *n - j - 1;
	if (k > 0) {
	    i__1 = j + k - 1;
	    for (i__ = j; i__ <= i__1; ++i__) {
		wr[i__] = wr[i__ + 2];
		wi[i__] = wi[i__ + 2];
/* L40: */
	    }
	    wr[*n - 1] = x;
	    wi[*n - 1] = y;
	    wr[*n] = x;
	    wi[*n] = -y;
	}
	*s = x + x;
	*p = x * x + y * y;
    }

    return 0;
/* *** End of SB01BX *** */
} /* sb01bx_ */

