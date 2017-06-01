/* MC01SD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mc01sd_(integer *dp, doublereal *p, integer *s, integer *
	t, doublereal *mant, integer *e, integer *iwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, m, v0, v1, lb, ub, dv, inc, beta;
    extern /* Subroutine */ int mc01sw_(doublereal *, integer *, doublereal *,
	     integer *);
    extern integer mc01sx_(integer *, integer *, integer *, doublereal *);
    extern /* Subroutine */ int mc01sy_(doublereal *, integer *, integer *, 
	    doublereal *, logical *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

/*     To scale the coefficients of the real polynomial P(x) such that */
/*     the coefficients of the scaled polynomial Q(x) = sP(tx) have */
/*     minimal variation, where s and t are real scalars. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial P(x).  DP >= 0. */

/*     P       (input/output) DOUBLE PRECISION array, dimension (DP+1) */
/*             On entry, this array must contain the coefficients of P(x) */
/*             in increasing powers of x. */
/*             On exit, this array contains the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     S       (output) INTEGER */
/*             The exponent of the floating-point representation of the */
/*             scaling factor s = BASE**S, where BASE is the base of the */
/*             machine representation of floating-point numbers (see */
/*             LAPACK Library Routine DLAMCH). */

/*     T       (output) INTEGER */
/*             The exponent of the floating-point representation of the */
/*             scaling factor t = BASE**T. */

/*     MANT    (output) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array contains the mantissas of the standard */
/*             floating-point representation of the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     E       (output) INTEGER array, dimension (DP+1) */
/*             This array contains the exponents of the standard */
/*             floating-point representation of the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (DP+1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, P(x) is the zero polynomial. */

/*     METHOD */

/*     Define the variation of the coefficients of the real polynomial */

/*                                         2                DP */
/*        P(x) = p(0) + p(1) * x + p(2) * x  + ... + p(DP) x */

/*     whose non-zero coefficients can be represented as */
/*                          e(i) */
/*        p(i) = m(i) * BASE     (where 1 <= ABS(m(i)) < BASE) */

/*     by */

/*        V = max(e(i)) - min(e(i)), */

/*     where max and min are taken over the indices i for which p(i) is */
/*     non-zero. */
/*                                        DP         i    i */
/*     For the scaled polynomial P(cx) = SUM p(i) * c  * x  with */
/*                                       i=0 */
/*                j */
/*     c  = (BASE) , the variation V(j) is given by */

/*       V(j) = max(e(i) + j * i) - min(e(i) + j * i). */

/*     Using the fact that V(j) is a convex function of j, the routine */
/*     determines scaling factors s = (BASE)**S and t = (BASE)**T such */
/*     that the coefficients of the scaled polynomial Q(x) = sP(tx) */
/*     satisfy the following conditions: */

/*       (a) 1 <= q(0) < BASE and */

/*       (b) the variation of the coefficients of Q(x) is minimal. */

/*     Further details can be found in [1]. */

/*     REFERENCES */

/*     [1] Dunaway, D.K. */
/*         Calculation of Zeros of a Real Polynomial through */
/*         Factorization using Euclid's Algorithm. */
/*         SIAM J. Numer. Anal., 11, pp. 1087-1104, 1974. */

/*     NUMERICAL ASPECTS */

/*     Since the scaling is performed on the exponents of the floating- */
/*     point representation of the coefficients of P(x), no rounding */
/*     errors occur during the computation of the coefficients of Q(x). */

/*     FURTHER COMMENTS */

/*     The scaling factors s and t are BASE dependent. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    --iwork;
    --e;
    --mant;
    --p;

    /* Function Body */
    if (*dp < 0) {
	*info = -1;

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC01SD", &i__1, (ftnlen)6);
	return 0;
    }

    *info = 0;
    lb = 1;
/*     WHILE ( LB <= DP+1 and P(LB) = 0 ) DO */
L20:
    if (lb <= *dp + 1) {
	if (p[lb] == 0.) {
	    ++lb;
	    goto L20;
	}
    }
/*     END WHILE 20 */

/*     LB = MIN( i: P(i) non-zero). */

    if (lb == *dp + 2) {
	*info = 1;
	return 0;
    }

    ub = *dp + 1;
/*     WHILE ( P(UB) = 0 ) DO */
L40:
    if (p[ub] == 0.) {
	--ub;
	goto L40;
    }
/*     END WHILE 40 */

/*     UB = MAX(i: P(i) non-zero). */

    beta = (integer) dlamch_("Base", (ftnlen)4);

    i__1 = *dp + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mc01sw_(&p[i__], &beta, &mant[i__], &e[i__]);
/* L60: */
    }

/*     First prescaling. */

    m = e[lb];
    if (m != 0) {

	i__1 = ub;
	for (i__ = lb; i__ <= i__1; ++i__) {
	    if (mant[i__] != 0.) {
		e[i__] -= m;
	    }
/* L80: */
	}

    }
    *s = -m;

/*     Second prescaling. */

    if (ub > 1) {
	d__1 = (doublereal) e[ub] / (doublereal) (ub - 1);
	m = i_dnnt(&d__1);
    }

    i__1 = ub;
    for (i__ = lb; i__ <= i__1; ++i__) {
	if (mant[i__] != 0.) {
	    e[i__] -= m * (i__ - 1);
	}
/* L100: */
    }

    *t = -m;

    v0 = mc01sx_(&lb, &ub, &e[1], &mant[1]);
    j = 1;

    i__1 = ub;
    for (i__ = lb; i__ <= i__1; ++i__) {
	if (mant[i__] != 0.) {
	    iwork[i__] = e[i__] + (i__ - 1);
	}
/* L120: */
    }

    v1 = mc01sx_(&lb, &ub, &iwork[1], &mant[1]);
    dv = v1 - v0;
    if (dv != 0) {
	if (dv > 0) {
	    j = 0;
	    inc = -1;
	    v1 = v0;
	    dv = -dv;

	    i__1 = ub;
	    for (i__ = lb; i__ <= i__1; ++i__) {
		iwork[i__] = e[i__];
/* L130: */
	    }

	} else {
	    inc = 1;
	}
/*        WHILE ( DV < 0 ) DO */
L140:
	if (dv < 0) {
	    v0 = v1;

	    i__1 = ub;
	    for (i__ = lb; i__ <= i__1; ++i__) {
		e[i__] = iwork[i__];
/* L150: */
	    }

	    j += inc;

	    i__1 = ub;
	    for (i__ = lb; i__ <= i__1; ++i__) {
		iwork[i__] = e[i__] + inc * (i__ - 1);
/* L160: */
	    }

	    v1 = mc01sx_(&lb, &ub, &iwork[1], &mant[1]);
	    dv = v1 - v0;
	    goto L140;
	}
/*        END WHILE 140 */
	*t = *t + j - inc;
    }

/*     Evaluation of the output parameters. */

    i__1 = ub;
    for (i__ = lb; i__ <= i__1; ++i__) {
	mc01sy_(&mant[i__], &e[i__], &beta, &p[i__], &ovflow);
/* L180: */
    }

    return 0;
/* *** Last line of MC01SD *** */
} /* mc01sd_ */

