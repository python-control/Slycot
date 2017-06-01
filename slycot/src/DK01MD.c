/* DK01MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int dk01md_(char *type__, integer *n, doublereal *a, integer 
	*info, ftnlen type_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__, n1;
    static doublereal fn, buf, temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical mtype, ntype;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical mntype;


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

/*     To apply an anti-aliasing window to a real signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             Indicates the type of window to be applied to the signal */
/*             as follows: */
/*             = 'M':  Hamming window; */
/*             = 'N':  Hann window; */
/*             = 'Q':  Quadratic window. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N >= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the signal to be */
/*             processed. */
/*             On exit, this array contains the windowing function. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If TYPE = 'M', then a Hamming window is applied to A(1),...,A(N), */
/*     which yields */
/*       _ */
/*       A(i) = (0.54 + 0.46*cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N. */

/*     If TYPE = 'N', then a Hann window is applied to A(1),...,A(N), */
/*     which yields */
/*       _ */
/*       A(i) = 0.5*(1 + cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N. */

/*     If TYPE = 'Q', then a quadratic window is applied to A(1),..., */
/*     A(N), which yields */
/*       _ */
/*       A(i) = (1 - 2*((i-1)/(N-1))**2)*(1 - (i-1)/(N-1))*A(i), */
/*                                             i = 1,2,...,(N-1)/2+1; */
/*       _ */
/*       A(i) = 2*(1 - ((i-1)/(N-1))**3)*A(i), i = (N-1)/2+2,...,N. */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DK01AD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Digital signal processing, Hamming window, Hann window, real */
/*     signals, windowing. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --a;

    /* Function Body */
    *info = 0;
    mtype = lsame_(type__, "M", (ftnlen)1, (ftnlen)1);
    ntype = lsame_(type__, "N", (ftnlen)1, (ftnlen)1);
    mntype = mtype || ntype;

/*     Test the input scalar arguments. */

    if (! mntype && ! lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n <= 0) {
	*info = -2;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("DK01MD", &i__1, (ftnlen)6);
	return 0;
    }

    fn = (doublereal) (*n - 1);
    if (mntype) {
	temp = atan(1.) * 4. / fn;
    }

    if (mtype) {

/*        Hamming window. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    a[i__] *= cos(temp * (doublereal) (i__ - 1)) * .46 + .54;
/* L10: */
	}

    } else if (ntype) {

/*        Hann window. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    a[i__] = a[i__] * .5 * (cos(temp * (doublereal) (i__ - 1)) + 1.);
/* L20: */
	}

    } else {

/*        Quadratic window. */

	n1 = (*n - 1) / 2 + 1;

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    buf = (doublereal) (i__ - 1) / fn;
/* Computing 2nd power */
	    d__1 = buf;
	    temp = d__1 * d__1;
	    if (i__ <= n1) {
		a[i__] = a[i__] * (1. - temp * 2.) * (1. - buf);
	    } else {
		a[i__] = a[i__] * 2. * (1. - buf * temp);
	    }
/* L30: */
	}

    }

    return 0;
/* *** Last line of DK01MD *** */
} /* dk01md_ */

