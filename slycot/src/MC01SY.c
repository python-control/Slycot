/* MC01SY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mc01sy_(doublereal *m, integer *e, integer *b, 
	doublereal *a, logical *ovflow)
{
    static integer et;
    static doublereal mt, base;
    static integer emin, emax, expon;
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

/*     To find a real number A from its mantissa M and its exponent E, */
/*     i.e., */
/*        A = M * B**E. */
/*     M and E need not be the standard floating-point values. */
/*     If ABS(A) < B**(EMIN-1), i.e. the smallest positive model number, */
/*     then the routine returns A = 0. */
/*     If M = 0, then the routine returns A = 0 regardless of the value */
/*     of E. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) DOUBLE PRECISION */
/*             The mantissa of the floating-point representation of A. */

/*     E       (input) INTEGER */
/*             The exponent of the floating-point representation of A. */

/*     B       (input) INTEGER */
/*             The base of the floating-point arithmetic. */

/*     A       (output) DOUBLE PRECISION */
/*             The value of M * B**E. */

/*     OVFLOW  (output) LOGICAL */
/*             The value .TRUE., if ABS(M) * B**E >= B**EMAX (where EMAX */
/*             is the largest possible exponent) and .FALSE. otherwise. */
/*             A is not defined if OVFLOW = .TRUE.. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GY by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    *ovflow = FALSE_;

    if (*m == 0. || *e == 0) {
	*a = *m;
	return 0;
    }

/*     Determination of the mantissa MT and the exponent ET of the */
/*     standard floating-point representation. */

    emin = (integer) dlamch_("Minimum exponent", (ftnlen)16);
    emax = (integer) dlamch_("Largest exponent", (ftnlen)16);
    mt = *m;
    et = *e;
/*     WHILE ( ABS( MT ) >= B ) DO */
L20:
    if (abs(mt) >= (doublereal) (*b)) {
	mt /= *b;
	++et;
	goto L20;
    }
/*     END WHILE 20 */
/*     WHILE ( ABS( MT ) < 1 ) DO */
L40:
    if (abs(mt) < 1.) {
	mt *= *b;
	--et;
	goto L40;
    }
/*     END WHILE 40 */

    if (et < emin) {
	*a = 0.;
	return 0;
    }

    if (et >= emax) {
	*ovflow = TRUE_;
	return 0;
    }

/*     Computation of the value of A by the relation */
/*     M * B**E = A * (BASE)**EXPON */

    expon = abs(et);
    *a = mt;
    base = (doublereal) (*b);
    if (et < 0) {
	base = 1. / base;
    }
/*     WHILE ( not EXPON = 0 ) DO */
L60:
    if (expon != 0) {
	if (expon % 2 == 0) {
	    base *= base;
	    expon /= 2;
	} else {
	    *a *= base;
	    --expon;
	}
	goto L60;
    }
/*     END WHILE 60 */

    return 0;
/* *** Last line of MC01SY *** */
} /* mc01sy_ */

