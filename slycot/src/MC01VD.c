/* MC01VD.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 1.;

/* Subroutine */ int mc01vd_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *z1re, doublereal *z1im, doublereal *z2re, doublereal *
	z2im, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal w, m1, m2;
    static integer ea, eb, ec, ed;
    static doublereal ma, mb, mc, md;
    static integer eb2;
    static doublereal big, absa, absb, absc;
    static integer beta;
    extern /* Subroutine */ int mc01sw_(doublereal *, integer *, doublereal *,
	     integer *);
    static doublereal sfmin;
    extern /* Subroutine */ int mc01sy_(doublereal *, integer *, integer *, 
	    doublereal *, logical *);
    extern doublereal dlamch_(char *, ftnlen);
    static integer eaplec;
    static logical ovflow;


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

/*     To compute the roots of a quadratic equation with real */
/*     coefficients. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the quadratic term. */

/*     B       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the linear term. */

/*     C       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the constant term. */

/*     Z1RE    (output) DOUBLE PRECISION */
/*     Z1IM    (output) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of the largest */
/*             root in magnitude. */

/*     Z2RE    (output) DOUBLE PRECISION */
/*     Z2IM    (output) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of the */
/*             smallest root in magnitude. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if on entry, either A = B = 0.0 or A = 0.0 and the */
/*                   root -C/B overflows; in this case Z1RE, Z1IM, Z2RE */
/*                   and Z2IM are unassigned; */
/*             = 2:  if on entry, A = 0.0; in this case Z1RE contains */
/*                   BIG and Z1IM contains zero, where BIG is a */
/*                   representable number near the overflow threshold */
/*                   of the machine (see LAPACK Library Routine DLAMCH); */
/*             = 3:  if on entry, either C = 0.0 and the root -B/A */
/*                   overflows or A, B and C are non-zero and the largest */
/*                   real root in magnitude cannot be computed without */
/*                   overflow; in this case Z1RE contains BIG and Z1IM */
/*                   contains zero; */
/*             = 4:  if the roots cannot be computed without overflow; in */
/*                   this case Z1RE, Z1IM, Z2RE and Z2IM are unassigned. */

/*     METHOD */

/*     The routine computes the roots (r1 and r2) of the real quadratic */
/*     equation */
/*             2 */
/*        a * x  + b * x + c = 0 */

/*     as */
/*             - b - SIGN(b) * SQRT(b * b - 4 * a * c)             c */
/*        r1 = ---------------------------------------  and r2 = ------ */
/*                              2 * a                            a * r1 */

/*     unless a = 0, in which case */

/*             -c */
/*        r1 = --. */
/*              b */

/*     Precautions are taken to avoid overflow and underflow wherever */
/*     possible. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01JD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Quadratic equation, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Detect special cases. */

    *info = 0;
    beta = (integer) dlamch_("Base", (ftnlen)4);
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
    big = 1. / sfmin;
    if (*a == 0.) {
	if (*b == 0.) {
	    *info = 1;
	} else {
	    ovflow = FALSE_;
	    *z2re = 0.;
	    if (*c__ != 0.) {
		absb = abs(*b);
		if (absb >= 1.) {
		    if (abs(*c__) >= absb * sfmin) {
			*z2re = -(*c__) / *b;
		    }
		} else {
		    if (abs(*c__) <= absb * big) {
			*z2re = -(*c__) / *b;
		    } else {
			ovflow = TRUE_;
			*z2re = big;
			if (d_sign(&c_b4, b) * d_sign(&c_b4, c__) > 0.) {
			    *z2re = -big;
			}
		    }
		}
	    }
	    if (ovflow) {
		*info = 1;
	    } else {
		*z1re = big;
		*z1im = 0.;
		*z2im = 0.;
		*info = 2;
	    }
	}
	return 0;
    }

    if (*c__ == 0.) {
	ovflow = FALSE_;
	*z1re = 0.;
	if (*b != 0.) {
	    absa = abs(*a);
	    if (absa >= 1.) {
		if (abs(*b) >= absa * sfmin) {
		    *z1re = -(*b) / *a;
		}
	    } else {
		if (abs(*b) <= absa * big) {
		    *z1re = -(*b) / *a;
		} else {
		    ovflow = TRUE_;
		    *z1re = big;
		}
	    }
	}
	if (ovflow) {
	    *info = 3;
	}
	*z1im = 0.;
	*z2re = 0.;
	*z2im = 0.;
	return 0;
    }

/*     A and C are non-zero. */

    if (*b == 0.) {
	ovflow = FALSE_;
	absc = sqrt((abs(*c__)));
	absa = sqrt((abs(*a)));
	w = 0.;
	if (absa >= 1.) {
	    if (absc >= absa * sfmin) {
		w = absc / absa;
	    }
	} else {
	    if (absc <= absa * big) {
		w = absc / absa;
	    } else {
		ovflow = TRUE_;
		w = big;
	    }
	}
	if (ovflow) {
	    *info = 4;
	} else {
	    if (d_sign(&c_b4, a) * d_sign(&c_b4, c__) > 0.) {
		*z1re = 0.;
		*z2re = 0.;
		*z1im = w;
		*z2im = -w;
	    } else {
		*z1re = w;
		*z2re = -w;
		*z1im = 0.;
		*z2im = 0.;
	    }
	}
	return 0;
    }

/*     A, B and C are non-zero. */

    mc01sw_(a, &beta, &ma, &ea);
    mc01sw_(b, &beta, &mb, &eb);
    mc01sw_(c__, &beta, &mc, &ec);

/*     Compute a 'near' floating-point representation of the discriminant */
/*     D = MD * BETA**ED. */

    eaplec = ea + ec;
    eb2 = eb << 1;
    if (eaplec > eb2) {
	d__1 = mb * mb;
	i__1 = eb2 - eaplec;
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
	w -= ma * 4. * mc;
	mc01sw_(&w, &beta, &md, &ed);
	ed += eaplec;
    } else {
	d__1 = ma * 4. * mc;
	i__1 = eaplec - eb2;
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
	w = mb * mb - w;
	mc01sw_(&w, &beta, &md, &ed);
	ed += eb2;
    }

    if (ed % 2 != 0) {
	++ed;
	md /= beta;
    }

/*     Complex roots. */

    if (md < 0.) {
	d__1 = -mb / (ma * 2);
	i__1 = eb - ea;
	mc01sy_(&d__1, &i__1, &beta, z1re, &ovflow);
	if (ovflow) {
	    *info = 4;
	} else {
	    d__1 = sqrt(-md) / (ma * 2);
	    i__1 = ed / 2 - ea;
	    mc01sy_(&d__1, &i__1, &beta, z1im, &ovflow);
	    if (ovflow) {
		*info = 4;
	    } else {
		*z2re = *z1re;
		*z2im = -(*z1im);
	    }
	}
	return 0;
    }

/*     Real roots. */

    md = sqrt(md);
    ed /= 2;
    if (ed > eb) {
	d__1 = abs(mb);
	i__1 = eb - ed;
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
	w += md;
	m1 = -d_sign(&c_b4, &mb) * w / (ma * 2);
	i__1 = ed - ea;
	mc01sy_(&m1, &i__1, &beta, z1re, &ovflow);
	if (ovflow) {
	    *z1re = big;
	    *info = 3;
	}
	m2 = -d_sign(&c_b4, &mb) * 2 * mc / w;
	i__1 = ec - ed;
	mc01sy_(&m2, &i__1, &beta, z2re, &ovflow);
    } else {
	i__1 = ed - eb;
	mc01sy_(&md, &i__1, &beta, &w, &ovflow);
	w += abs(mb);
	m1 = -d_sign(&c_b4, &mb) * w / (ma * 2);
	i__1 = eb - ea;
	mc01sy_(&m1, &i__1, &beta, z1re, &ovflow);
	if (ovflow) {
	    *z1re = big;
	    *info = 3;
	}
	m2 = -d_sign(&c_b4, &mb) * 2 * mc / w;
	i__1 = ec - eb;
	mc01sy_(&m2, &i__1, &beta, z2re, &ovflow);
    }
    *z1im = 0.;
    *z2im = 0.;

    return 0;
/* *** Last line of MC01VD *** */
} /* mc01vd_ */

