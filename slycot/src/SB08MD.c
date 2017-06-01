/* SB08MD.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b30 = -1.;

/* Subroutine */ int sb08md_(char *acona, integer *da, doublereal *a, 
	doublereal *res, doublereal *e, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen acona_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal w, a0;
    static integer i0, nc;
    static doublereal si;
    static integer lq;
    static doublereal mu;
    static integer da1;
    static doublereal xda;
    static integer lay;
    static doublereal eps, muj;
    static integer binc, ldif, lphi;
    static logical conv;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal signi, signj;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb08my_(integer *, doublereal *, 
	    doublereal *, doublereal *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal signi0, simin1, sqrta0;
    static integer lambda;
    extern doublereal dlamch_(char *, ftnlen);
    static logical lacona;
    extern integer idamax_(integer *, doublereal *, integer *);
    static logical stable;
    static integer lphend, layend;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal tolphi, sqrtmj, sqrtmu;


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

/*     To compute a real polynomial E(s) such that */

/*        (a)  E(-s) * E(s) = A(-s) * A(s) and */
/*        (b)  E(s) is stable - that is, all the zeros of E(s) have */
/*             non-positive real parts, */

/*     which corresponds to computing the spectral factorization of the */
/*     real polynomial A(s) arising from continuous optimality problems. */

/*     The input polynomial may be supplied either in the form */

/*        A(s) = a(0) + a(1) * s + ... + a(DA) * s**DA */

/*     or as */

/*        B(s) = A(-s) * A(s) */
/*             = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA)        (1) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ACONA   CHARACTER*1 */
/*             Indicates whether the coefficients of A(s) or B(s) = */
/*             A(-s) * A(s) are to be supplied as follows: */
/*             = 'A':  The coefficients of A(s) are to be supplied; */
/*             = 'B':  The coefficients of B(s) are to be supplied. */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(s) and E(s).  DA >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (DA+1) */
/*             On entry, this array must contain either the coefficients */
/*             of the polynomial A(s) in increasing powers of s if */
/*             ACONA = 'A', or the coefficients of the polynomial B(s) in */
/*             increasing powers of s**2 (see equation (1)) if ACONA = */
/*             'B'. */
/*             On exit, this array contains the coefficients of the */
/*             polynomial B(s) in increasing powers of s**2. */

/*     RES     (output) DOUBLE PRECISION */
/*             An estimate of the accuracy with which the coefficients of */
/*             the polynomial E(s) have been computed (see also METHOD */
/*             and NUMERICAL ASPECTS). */

/*     E       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             The coefficients of the spectral factor E(s) in increasing */
/*             powers of s. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 5*DA+5. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, A(I) = 0.0, for I = 1,2,...,DA+1. */
/*             = 2:  if on entry, ACONA = 'B' but the supplied */
/*                   coefficients of the polynomial B(s) are not the */
/*                   coefficients of A(-s) * A(s) for some real A(s); */
/*                   in this case, RES and E are unassigned; */
/*             = 3:  if the iterative process (see METHOD) has failed to */
/*                   converge in 30 iterations; */
/*             = 4:  if the last computed iterate (see METHOD) is */
/*                   unstable. If ACONA = 'B', then the supplied */
/*                   coefficients of the polynomial B(s) may not be the */
/*                   coefficients of A(-s) * A(s) for some real A(s). */

/*     METHOD */
/*         _                                               _ */
/*     Let A(s) be the conjugate polynomial of A(s), i.e., A(s) = A(-s). */

/*     The method used by the routine is based on applying the */
/*     Newton-Raphson iteration to the function */
/*               _       _ */
/*        F(e) = A * A - e * e, */

/*     which leads to the iteration formulae (see [1]): */

/*        _(i)   (i)  _(i)   (i)     _      ) */
/*        q   * x   + x   * q    = 2 A * A  ) */
/*                                          )   for i = 0, 1, 2,... */
/*         (i+1)    (i)   (i)               ) */
/*        q     = (q   + x   )/2            ) */

/*                    (0)         DA */
/*     Starting from q   = (1 + s)   (which has no zeros in the closed */
/*                                                  (1)   (2)   (3) */
/*     right half-plane), the sequence of iterates q   , q   , q   ,... */
/*     converges to a solution of F(e) = 0 which has no zeros in the */
/*     open right half-plane. */

/*     The iterates satisfy the following conditions: */

/*              (i) */
/*        (a)  q   is a stable polynomial (no zeros in the closed right */
/*             half-plane) and */

/*              (i)        (i-1) */
/*        (b)  q   (1) <= q     (1). */

/*                                       (i-1)                       (i) */
/*     The iterative process stops with q     , (where i <= 30)  if q */
/*     violates either (a) or (b), or if the condition */
/*                       _(i) (i)  _ */
/*        (c)  RES  = ||(q   q   - A A)|| < tol, */

/*     is satisfied, where || . || denotes the largest coefficient of */
/*                     _(i) (i)  _ */
/*     the polynomial (q   q   - A A) and tol is an estimate of the */
/*                                                    _(i)  (i) */
/*     rounding error in the computed coefficients of q    q   . If there */
/*     is no convergence after 30 iterations then the routine returns */
/*     with the Error Indicator (INFO) set to 3, and the value of RES may */
/*     indicate whether or not the last computed iterate is close to the */
/*     solution. */

/*     If ACONA = 'B', then it is possible that the equation e(-s) * */
/*     e(s) = B(s) has no real solution, which will be the case if A(1) */
/*     < 0 or if ( -1)**DA * A(DA+1) < 0. */

/*     REFERENCES */

/*     [1] Vostry, Z. */
/*         New Algorithm for Polynomial Spectral Factorization with */
/*         Quadratic Convergence II. */
/*         Kybernetika, 12, pp. 248-259, 1976. */

/*     NUMERICAL ASPECTS */

/*     The conditioning of the problem depends upon the distance of the */
/*     zeros of A(s) from the imaginary axis and on their multiplicity. */
/*     For a well-conditioned problem the accuracy of the computed */
/*     coefficients of E(s) is of the order of RES. However, for problems */
/*     with zeros near the imaginary axis or with multiple zeros, the */
/*     value of RES may be an overestimate of the true accuracy. */

/*     FURTHER COMMENTS */

/*     In order for the problem e(-s) * e(s) = B(s) to have a real */
/*     solution e(s), it is necessary and sufficient that B(j*omega) */
/*     >= 0 for any purely imaginary argument j*omega (see [1]). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08AD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Factorization, Laplace transform, optimal control, optimal */
/*     filtering, polynomial operations, spectral factorization, zeros. */

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
    --dwork;
    --e;
    --a;

    /* Function Body */
    *info = 0;
    lacona = lsame_(acona, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lacona && ! lsame_(acona, "B", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*da < 0) {
	*info = -2;
    } else if (*ldwork < *da * 5 + 5) {
	*info = -7;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB08MD", &i__1, (ftnlen)6);
	return 0;
    }

    if (! lacona) {
	i__1 = *da + 1;
	dcopy_(&i__1, &a[1], &c__1, &e[1], &c__1);
    } else {
	w = 0.;
	sb08my_(da, &a[1], &e[1], &w);
    }

/*     Reduce E such that the first and the last element are non-zero. */

    da1 = *da + 1;

/*     WHILE ( DA1 >= 1 and E(DA1) = 0 ) DO */
L20:
    if (da1 >= 1) {
	if (e[da1] == 0.) {
	    --da1;
	    goto L20;
	}
    }
/*     END WHILE 20 */

    --da1;
    if (da1 < 0) {
	*info = 1;
	return 0;
    }

    i0 = 1;

/*     WHILE ( E(I0) = 0 ) DO */
L40:
    if (e[i0] == 0.) {
	++i0;
	goto L40;
    }
/*     END WHILE 40 */

    --i0;
    if (i0 != 0) {
	if (i0 % 2 == 0) {
	    signi0 = 1.;
	} else {
	    signi0 = -1.;
	}

	i__1 = da1 - i0 + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    e[i__] = signi0 * e[i__ + i0];
/* L60: */
	}

	da1 -= i0;
    }
    if (da1 % 2 == 0) {
	signi = 1.;
    } else {
	signi = -1.;
    }
    nc = da1 + 1;
    if (e[1] < 0. || e[nc] * signi < 0.) {
	*info = 2;
	return 0;
    }

/*     Initialization. */

    eps = dlamch_("Epsilon", (ftnlen)7);
    si = 1. / dlamch_("Safe minimum", (ftnlen)12);
    lq = 1;
    lay = lq + nc;
    lambda = lay + nc;
    lphi = lambda + nc;
    ldif = lphi + nc;

    a0 = e[1];
    binc = 1;

/*     Computation of the starting polynomial and scaling of the input */
/*     polynomial. */

    d__2 = a0 / (d__1 = e[nc], abs(d__1));
    d__3 = 1. / (doublereal) da1;
    mu = pow_dd(&d__2, &d__3);
    muj = 1.;

    i__1 = nc;
    for (j = 1; j <= i__1; ++j) {
	w = e[j] * muj / a0;
	a[j] = w;
	e[j] = (doublereal) binc;
	dwork[lq + j - 1] = (doublereal) binc;
	muj *= mu;
	binc = binc * (nc - j) / j;
/* L80: */
    }

    conv = FALSE_;
    stable = TRUE_;

/*     The contents of the arrays is, cf [1], */

/*     E : the last computed stable polynomial q   ; */
/*                                              i-1 */
/*     DWORK(LAY+1,...,LAY+DA1-1)  : a'(1), ..., a'(DA1-1), these values */
/*                                   are changed during the computation */
/*                                   into y; */
/*          (LAMBDA+1,...,LAMBDA+DA1-2) : lambda(1), ..., lambda(DA1-2), */
/*                                        the factors of the Routh */
/*                                        stability test, (lambda(i) is */
/*                                        P(i) in [1]); */
/*          (LPHI+1,...,LPHI+DA1-1) : phi(1), ..., phi(DA1-1), the values */
/*                                    phi(i,j), see [1], scheme (11); */
/*          (LDIF,...,LDIF+DA1) : the coeffs of q (-s) * q (s) - b(s). */
/*                                               i        i */
/*     DWORK(LQ,...,LQ+DA1) : the last computed polynomial q . */
/*                                                          i */
    i__ = 0;

/*     WHILE ( I < 30 and CONV = FALSE and STABLE = TRUE ) DO */
L100:
    if (i__ < 30 && ! conv && stable) {
	++i__;
	dcopy_(&nc, &a[1], &c__1, &dwork[lay], &c__1);
	dcopy_(&nc, &dwork[lq], &c__1, &dwork[lphi], &c__1);
	m = da1 / 2;
	layend = lay + da1;
	lphend = lphi + da1;
	xda = a[nc] / dwork[lq + da1];

	i__1 = m;
	for (k = 1; k <= i__1; ++k) {
	    dwork[lay + k] -= dwork[lphi + (k << 1)];
	    dwork[layend - k] -= dwork[lphend - (k << 1)] * xda;
/* L120: */
	}

/*        Computation of lambda(k) and y(k). */

	k = 1;

/*        WHILE ( K <= DA1 - 2 and STABLE = TRUE ) DO */
L140:
	if (k <= da1 - 2 && stable) {
	    if (dwork[lphi + k] <= 0.) {
		stable = FALSE_;
	    }
	    if (stable) {
		w = dwork[lphi + k - 1] / dwork[lphi + k];
		dwork[lambda + k] = w;
		i__1 = (da1 - k) / 2;
		d__1 = -w;
		daxpy_(&i__1, &d__1, &dwork[lphi + k + 2], &c__2, &dwork[lphi 
			+ k + 1], &c__2);
		w = dwork[lay + k] / dwork[lphi + k];
		dwork[lay + k] = w;
		i__1 = (da1 - k) / 2;
		d__1 = -w;
		daxpy_(&i__1, &d__1, &dwork[lphi + k + 2], &c__2, &dwork[lay 
			+ k + 1], &c__1);
		++k;
	    }
	    goto L140;
	}
/*        END WHILE 140 */

	if (dwork[lphi + da1 - 1] <= 0.) {
	    stable = FALSE_;
	} else {
	    dwork[lay + da1 - 1] /= dwork[lphi + da1 - 1];
	}

/*        STABLE = The polynomial q    is stable. */
/*                                 i-1 */
	if (stable) {

/*           Computation of x  and q . */
/*                           i      i */

	    for (k = da1 - 2; k >= 1; --k) {
		w = dwork[lambda + k];
		i__1 = (da1 - k) / 2;
		d__1 = -w;
		daxpy_(&i__1, &d__1, &dwork[lay + k + 1], &c__2, &dwork[lay + 
			k], &c__2);
/* L160: */
	    }

	    dwork[lay + da1] = xda;

	    dcopy_(&nc, &dwork[lq], &c__1, &e[1], &c__1);
	    simin1 = si;
	    si = dwork[lq];
	    signj = -1.;

	    i__1 = da1;
	    for (j = 1; j <= i__1; ++j) {
		w = (dwork[lq + j] + signj * dwork[lay + j]) * .5;
		dwork[lq + j] = w;
		si += w;
		signj = -signj;
/* L180: */
	    }

	    tolphi = eps;
	    sb08my_(&da1, &e[1], &dwork[ldif], &tolphi);
	    daxpy_(&nc, &c_b30, &a[1], &c__1, &dwork[ldif], &c__1);
	    *res = (d__1 = dwork[idamax_(&nc, &dwork[ldif], &c__1) + ldif - 1]
		    , abs(d__1));

/*           Convergency test. */

	    if (si > simin1 || *res < tolphi) {
		conv = TRUE_;
	    }
	    goto L100;
	}
    }
/*     END WHILE 100 */

/*     Backscaling. */

    mu = 1. / mu;
    sqrta0 = sqrt(a0);
    sqrtmu = sqrt(mu);
    muj = 1.;
    sqrtmj = 1.;

    i__1 = nc;
    for (j = 1; j <= i__1; ++j) {
	e[j] = e[j] * sqrta0 * sqrtmj;
	a[j] = a[j] * a0 * muj;
	muj *= mu;
	sqrtmj *= sqrtmu;
/* L200: */
    }

    if (i0 != 0) {

	for (j = nc; j >= 1; --j) {
	    e[i0 + j] = e[j];
	    a[i0 + j] = signi0 * a[j];
/* L220: */
	}

	i__1 = i0;
	for (j = 1; j <= i__1; ++j) {
	    e[j] = 0.;
	    a[j] = 0.;
/* L240: */
	}

    }

    if (! conv) {
	if (stable) {
	    *info = 3;
	} else {
	    *info = 4;
	}
    }

    return 0;
/* *** Last line of SB08MD *** */
} /* sb08md_ */

