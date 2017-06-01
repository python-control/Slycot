/* TD05AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b13 = 90.;

/* Subroutine */ int td05ad_(char *unitf, char *output, integer *np1, integer 
	*mp1, doublereal *w, doublereal *a, doublereal *b, doublereal *valr, 
	doublereal *vali, integer *info, ftnlen unitf_len, ftnlen output_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double atan(doublereal), pow_di(doublereal *, integer *), d_imag(
	    doublecomplex *), d_sign(doublereal *, doublereal *), d_lg10(
	    doublereal *);

    /* Local variables */
    static doublereal g;
    static integer i__, m, n, m2, n2;
    static doublereal w2, wc, bimag, breal, timag;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal treal;
    static doublecomplex ztemp;
    static doublereal twopi;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static integer iphase;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static logical lunitf;
    static integer npzero, nzzero;
    static logical loutpu;


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

/*     Given a complex valued rational function of frequency (transfer */
/*     function) G(jW) this routine will calculate its complex value or */
/*     its magnitude and phase for a specified frequency value. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UNITF   CHARACTER*1 */
/*             Indicates the choice of frequency unit as follows: */
/*             = 'R':  Input frequency W in radians/second; */
/*             = 'H':  Input frequency W in hertz. */

/*     OUTPUT  CHARACTER*1 */
/*             Indicates the choice of co-ordinates for output as folows: */
/*             = 'C':  Cartesian co-ordinates (output real and imaginary */
/*                     parts of G(jW)); */
/*             = 'P':  Polar co-ordinates (output magnitude and phase */
/*                     of G(jW)). */

/*     Input/Output Parameters */

/*     NP1     (input) INTEGER */
/*             The order of the denominator + 1, i.e. N + 1.  NP1 >= 1. */

/*     MP1     (input) INTEGER */
/*             The order of the numerator + 1, i.e. M + 1.  MP1 >= 1. */

/*     W       (input) DOUBLE PRECISION */
/*             The frequency value W for which the transfer function is */
/*             to be evaluated. */

/*     A       (input) DOUBLE PRECISION array, dimension (NP1) */
/*             This array must contain the vector of denominator */
/*             coefficients in ascending order of powers. That is, A(i) */
/*             must contain the coefficient of (jW)**(i-1) for i = 1, */
/*             2,...,NP1. */

/*     B       (input) DOUBLE PRECISION array, dimension (MP1) */
/*             This array must contain the vector of numerator */
/*             coefficients in ascending order of powers. That is, B(i) */
/*             must contain the coefficient of (jW)**(i-1) for i = 1, */
/*             2,...,MP1. */

/*     VALR    (output) DOUBLE PRECISION */
/*             If OUTPUT = 'C', VALR contains the real part of G(jW). */
/*             If OUTPUT = 'P', VALR contains the magnitude of G(jW) */
/*                              in dBs. */

/*     VALI    (output) DOUBLE PRECISION */
/*             If OUTPUT = 'C', VALI contains the imaginary part of */
/*                              G(jW). */
/*             If OUTPUT = 'P', VALI contains the phase of G(jW) in */
/*                              degrees. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the frequency value W is a pole of G(jW), or all */
/*                   the coefficients of the A polynomial are zero. */

/*     METHOD */

/*     By substituting the values of A, B and W in the following */
/*     formula: */

/*            B(1)+B(2)*(jW)+B(3)*(jW)**2+...+B(MP1)*(jW)**(MP1-1) */
/*     G(jW) = ---------------------------------------------------. */
/*            A(1)+A(2)*(jW)+A(3)*(jW)**2+...+A(NP1)*(jW)**(NP1-1) */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N+M) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TD01AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, March 1981. */

/*     REVISIONS */

/*     February 1997. */
/*     February 22, 1998 (changed the name of TD01MD). */

/*     KEYWORDS */

/*     Elementary polynomial operations, frequency response, matrix */
/*     fraction, polynomial matrix, state-space representation, transfer */
/*     matrix. */

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
    --b;
    --a;

    /* Function Body */
    *info = 0;
    lunitf = lsame_(unitf, "H", (ftnlen)1, (ftnlen)1);
    loutpu = lsame_(output, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lunitf && ! lsame_(unitf, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! loutpu && ! lsame_(output, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*np1 < 1) {
	*info = -3;
    } else if (*mp1 < 1) {
	*info = -4;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TD05AD", &i__1, (ftnlen)6);
	return 0;
    }

    m = *mp1 - 1;
    n = *np1 - 1;
    wc = *w;
    twopi = atan(1.) * 8.;
    if (lunitf) {
	wc *= twopi;
    }
/* Computing 2nd power */
    d__1 = wc;
    w2 = d__1 * d__1;

/*     Determine the orders z (NZZERO) and p (NPZERO) of the factors */
/*     (jW)**k in the numerator and denominator polynomials, by counting */
/*     the zero trailing coefficients.  The value of G(jW) will then be */
/*     computed as (jW)**(z-p)*m(jW)/n(jW), for appropriate m and n. */

    i__ = 0;

L10:
    ++i__;
    if (i__ <= m) {
	if (b[i__] == 0.) {
	    goto L10;
	}
    }

    nzzero = i__ - 1;
    i__ = 0;

L20:
    ++i__;
    if (i__ <= n) {
	if (a[i__] == 0.) {
	    goto L20;
	}
    }

    npzero = i__ - 1;
    iphase = nzzero - npzero;

    m2 = (m - nzzero) % 2;

/*     Add real parts of the numerator m(jW). */

    treal = b[*mp1 - m2];

    i__1 = nzzero + 1;
    for (i__ = m - 1 - m2; i__ >= i__1; i__ += -2) {
	treal = b[i__] - w2 * treal;
/* L30: */
    }

/*     Add imaginary parts of the numerator m(jW). */

    if (m == 0) {
	timag = 0.;
    } else {
	timag = b[m + m2];

	i__1 = nzzero + 2;
	for (i__ = m + m2 - 2; i__ >= i__1; i__ += -2) {
	    timag = b[i__] - w2 * timag;
/* L40: */
	}

	timag *= wc;
    }

    n2 = (n - npzero) % 2;

/*     Add real parts of the denominator n(jW). */

    breal = a[*np1 - n2];

    i__1 = npzero + 1;
    for (i__ = n - 1 - n2; i__ >= i__1; i__ += -2) {
	breal = a[i__] - w2 * breal;
/* L50: */
    }

/*     Add imaginary parts of the denominator n(jW). */

    if (n == 0) {
	bimag = 0.;
    } else {
	bimag = a[n + n2];

	i__1 = npzero + 2;
	for (i__ = n + n2 - 2; i__ >= i__1; i__ += -2) {
	    bimag = a[i__] - w2 * bimag;
/* L60: */
	}

	bimag *= wc;
    }

/* Computing MAX */
    d__1 = abs(breal), d__2 = abs(bimag);
    if (max(d__1,d__2) == 0. || *w == 0. && iphase < 0) {

/*        Error return:  The specified frequency W is a pole of G(jW), */
/*              or all the coefficients of the A polynomial are zero. */

	*info = 1;
    } else {

/*        Evaluate the complex number W**(z-p)*m(jW)/n(jW). */

	z__2.r = treal, z__2.i = timag;
	z__3.r = breal, z__3.i = bimag;
	zladiv_(&z__1, &z__2, &z__3);
	ztemp.r = z__1.r, ztemp.i = z__1.i;
	*valr = ztemp.r * pow_di(&wc, &iphase);
	*vali = d_imag(&ztemp) * pow_di(&wc, &iphase);

	if (! loutpu) {

/*           Cartesian co-ordinates: Update the result for j**(z-p). */

	    i__ = abs(iphase) % 4;
	    if (iphase > 0 && i__ > 1 || iphase < 0 && (i__ == 1 || i__ == 2))
		     {
		*valr = -(*valr);
		*vali = -(*vali);
	    }

	    if (i__ % 2 != 0) {
		g = *valr;
		*valr = -(*vali);
		*vali = g;
	    }

	} else {

/*           Polar co-ordinates: Compute the magnitude and phase. */

	    g = dlapy2_(valr, vali);

	    if (*valr == 0.) {
		*vali = d_sign(&c_b13, vali);
	    } else {
		*vali = atan(*vali / *valr) / twopi * 360.;
		if (*vali == 0. && nzzero == m && npzero == n && b[nzzero + 1]
			 * a[npzero + 1] < 0.) {
		    *vali = 180.;
		}
	    }

	    *valr = d_lg10(&g) * 20.;

	    if (iphase != 0) {
		*vali += (doublereal) (nzzero - npzero) * 90.;
	    }
	}

    }

    return 0;
/* *** Last line of TD05AD *** */
} /* td05ad_ */

