/* SB02CX.f -- translated by f2c (version 20100827).
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

logical sb02cx_(doublereal *reig, doublereal *ieig)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static doublereal eps, tol;
    extern doublereal dlamch_(char *, ftnlen);


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

/*     To select the purely imaginary eigenvalues in computing the */
/*     H-infinity norm of a system. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     REIG    (input) DOUBLE PRECISION */
/*             The real part of the current eigenvalue considered. */

/*     IEIG    (input) DOUBLE PRECISION */
/*             The imaginary part of the current eigenvalue considered. */

/*     METHOD */

/*     The function value SB02CX is set to .TRUE. for a purely imaginary */
/*     eigenvalue and to .FALSE., otherwise. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     P. Hr. Petkov, Technical University of Sofia, May, 1999. */

/*     REVISIONS */

/*     P. Hr. Petkov, Technical University of Sofia, Oct. 2000. */

/*     KEYWORDS */

/*     H-infinity norm, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Get the machine precision. */

    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Set the tolerance in the determination of the purely */
/*     imaginary eigenvalues. */

    tol = eps * 100.;
    ret_val = abs(*reig) < tol;

    return ret_val;
/* *** Last line of SB02CX *** */
} /* sb02cx_ */

