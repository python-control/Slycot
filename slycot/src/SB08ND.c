/* SB08ND.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 2.;
static integer c_n1 = -1;
static doublereal c_b24 = -1.;

/* Subroutine */ int sb08nd_(char *acona, integer *da, doublereal *a, 
	doublereal *res, doublereal *e, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen acona_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, w, a0;
    static integer nc, lq;
    static doublereal sa0;
    static integer nck, lro;
    static doublereal res0;
    static integer leta;
    static logical conv;
    static doublereal tolq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), sb08ny_(integer *, doublereal *, 
	    doublereal *, doublereal *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer lambda;
    static logical lacona;
    static integer lalpha;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical hurwtz;


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

/*     To compute a real polynomial E(z) such that */

/*        (a)  E(1/z) * E(z) = A(1/z) * A(z) and */
/*        (b)  E(z) is stable - that is, E(z) has no zeros with modulus */
/*             greater than 1, */

/*     which corresponds to computing the spectral factorization of the */
/*     real polynomial A(z) arising from discrete optimality problems. */

/*     The input polynomial may be supplied either in the form */

/*     A(z) = a(0) + a(1) * z + ... + a(DA) * z**DA */

/*     or as */

/*     B(z) = A(1/z) * A(z) */
/*          = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA) */
/*                                                                    (1) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ACONA   CHARACTER*1 */
/*             Indicates whether the coefficients of A(z) or B(z) = */
/*             A(1/z) * A(z) are to be supplied as follows: */
/*             = 'A':  The coefficients of A(z) are to be supplied; */
/*             = 'B':  The coefficients of B(z) are to be supplied. */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(z) and E(z).  DA >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (DA+1) */
/*             On entry, if ACONA = 'A', this array must contain the */
/*             coefficients of the polynomial A(z) in increasing powers */
/*             of z, and if ACONA = 'B', this array must contain the */
/*             coefficients b ,b ,...,b   of the polynomial B(z) in */
/*                           0  1      DA */
/*             equation (1). That is, A(i) = b    for i = 1,2,...,DA+1. */
/*                                            i-1 */
/*             On exit, this array contains the coefficients of the */
/*             polynomial B(z) in eqation (1). Specifically, A(i) */
/*             contains b   ,  for i = 1,2,...DA+1. */
/*                       i-1 */

/*     RES     (output) DOUBLE PRECISION */
/*             An estimate of the accuracy with which the coefficients of */
/*             the polynomial E(z) have been computed (see also METHOD */
/*             and NUMERICAL ASPECTS). */

/*     E       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             The coefficients of the spectral factor E(z) in increasing */
/*             powers of z. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 5*DA+5. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  if on entry, ACONA = 'B' but the supplied */
/*                   coefficients of the polynomial B(z) are not the */
/*                   coefficients of A(1/z) * A(z) for some real A(z); */
/*                   in this case, RES and E are unassigned; */
/*             = 3:  if the iterative process (see METHOD) has failed to */
/*                   converge in 30 iterations; */
/*             = 4:  if the last computed iterate (see METHOD) is */
/*                   unstable. If ACONA = 'B', then the supplied */
/*                   coefficients of the polynomial B(z) may not be the */
/*                   coefficients of A(1/z) * A(z) for some real A(z). */

/*     METHOD */
/*         _                                               _ */
/*     Let A(z) be the conjugate polynomial of A(z), i.e., A(z) = A(1/z). */

/*     The method used by the routine is based on applying the */
/*     Newton-Raphson iteration to the function */
/*               _       _ */
/*        F(e) = A * A - e * e, */

/*     which leads to the iteration formulae (see [1] and [2]) */

/*        _(i)   (i)  _(i)   (i)     _      ) */
/*        q   * x   + x   * q    = 2 A * A  ) */
/*                                          )   for i = 0, 1, 2,... */
/*         (i+1)    (i)   (i)               ) */
/*        q     = (q   + x   )/2            ) */

/*     The iteration starts from */

/*         (0)                                        DA */
/*        q   (z) = (b(0) + b(1) * z + ... + b(DA) * z  ) / SQRT( b(0)) */

/*     which is a Hurwitz polynomial that has no zeros in the closed unit */
/*                                            (i) */
/*     circle (see [2], Theorem 3). Then lim q   = e, the convergence is */
/*     uniform and e is a Hurwitz polynomial. */

/*     The iterates satisfy the following conditions: */
/*              (i) */
/*        (a)  q    has no zeros in the closed unit circle, */
/*              (i)     (i-1) */
/*        (b)  q    <= q     and */
/*              0       0 */
/*              DA   (i) 2    DA     2 */
/*        (c)  SUM (q   )  - SUM (A )  >= 0. */
/*             k=0   k       k=0   k */
/*                                     (i) */
/*     The iterative process stops if q    violates (a), (b) or (c), */
/*     or if the condition */
/*                       _(i) (i)  _ */
/*        (d)  RES  = ||(q   q   - A A)|| < tol, */

/*     is satisfied, where || . || denotes the largest coefficient of */
/*                     _(i) (i)  _ */
/*     the polynomial (q   q   - A A) and tol is an estimate of the */
/*                                                    _(i)  (i) */
/*     rounding error in the computed coefficients of q    q   . If */
/*                                            (i-1) */
/*     condition (a) or (b) is violated then q      is taken otherwise */
/*      (i) */
/*     q    is used. Thus the computed reciprocal polynomial E(z) = z**DA */
/*     * q(1/z) is stable. If there is no convergence after 30 iterations */
/*     then the routine returns with the Error Indicator (INFO) set to 3, */
/*     and the value of RES may indicate whether or not the last computed */
/*     iterate is close to the solution. */
/*                                               (0) */
/*     If ACONA = 'B', then it is possible that q    is not a Hurwitz */
/*     polynomial, in which case the equation e(1/z) * e(z) = B(z) has no */
/*     real solution (see [2], Theorem 3). */

/*     REFERENCES */

/*     [1] Kucera, V. */
/*         Discrete Linear Control, The polynomial Approach. */
/*         John Wiley & Sons, Chichester, 1979. */

/*     [2] Vostry, Z. */
/*         New Algorithm for Polynomial Spectral Factorization with */
/*         Quadratic Convergence I. */
/*         Kybernetika, 11, pp. 415-422, 1975. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08BD by F. Delebecque and */
/*     A.J. Geurts. */

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
	xerbla_("SB08ND", &i__1, (ftnlen)6);
	return 0;
    }

    nc = *da + 1;
    if (! lacona) {
	if (a[1] <= 0.) {
	    *info = 2;
	    return 0;
	}
	dcopy_(&nc, &a[1], &c__1, &e[1], &c__1);
    } else {
	sb08ny_(da, &a[1], &e[1], &w);
    }

/*     Initialization. */

    lalpha = 1;
    lro = lalpha + nc;
    leta = lro + nc;
    lambda = leta + nc;
    lq = lambda + nc;

    a0 = e[1];
    sa0 = sqrt(a0);
    s = 0.;

    i__1 = nc;
    for (j = 1; j <= i__1; ++j) {
	w = e[j];
	a[j] = w;
	w /= sa0;
	e[j] = w;
	dwork[lq - 1 + j] = w;
/* Computing 2nd power */
	d__1 = w;
	s += d__1 * d__1;
/* L20: */
    }

    res0 = s - a0;

/*     The contents of the arrays is, cf [1], Section 7.6, */

/*     E : the last computed Hurwitz polynomial q   ; */
/*                                               i-1 */
/*     DWORK(LALPHA,..,LALPHA+DA-K)  : alpha(k,0),...alpha(k,n-k); */
/*          (LRO,...,LRO+DA-K)       : alpha(k,n-k),...,alpha(k); */
/*          (LETA,...,LETA+DA)       : eta(0),...,eta(n); */
/*          (LAMBDA,...,LAMBDA+DA-1) : lambda(0),...,lambda(n-1) */

/*     DWORK(LQ,...,LQ+DA) : the last computed polynomial q . */
/*                                                         i */
    i__ = 0;
    conv = FALSE_;
    hurwtz = TRUE_;

/*     WHILE ( I < 30 and CONV = FALSE and HURWTZ = TRUE ) DO */
L40:
    if (i__ < 30 && ! conv && hurwtz) {
	++i__;
	dcopy_(&nc, &a[1], &c__1, &dwork[leta], &c__1);
	dscal_(&nc, &c_b11, &dwork[leta], &c__1);
	dcopy_(&nc, &dwork[lq], &c__1, &dwork[lalpha], &c__1);

/*        Computation of lambda(k) and eta(k). */

	k = 1;

/*        WHILE ( K <= DA and HURWTZ = TRUE ) DO */
L60:
	if (k <= *da && hurwtz) {
	    nck = nc - k;
	    i__1 = nck + 1;
	    dcopy_(&i__1, &dwork[lalpha], &c_n1, &dwork[lro], &c__1);
	    w = dwork[lalpha + nck] / dwork[lro + nck];
	    if (abs(w) >= 1.) {
		hurwtz = FALSE_;
	    }
	    if (hurwtz) {
		dwork[lambda + k - 1] = w;
		d__1 = -w;
		daxpy_(&nck, &d__1, &dwork[lro], &c__1, &dwork[lalpha], &c__1)
			;
		w = dwork[leta + nck] / dwork[lalpha];
		dwork[leta + nck] = w;
		i__1 = nck - 1;
		d__1 = -w;
		daxpy_(&i__1, &d__1, &dwork[lalpha + 1], &c_n1, &dwork[leta + 
			1], &c__1);
		++k;
	    }
	    goto L60;
	}
/*        END WHILE 60 */

/*        HURWTZ = The polynomial q    is a Hurwitz polynomial. */
/*                                 i-1 */
	if (hurwtz) {
	    dcopy_(&nc, &dwork[lq], &c__1, &e[1], &c__1);

/*           Accuracy test. */

	    sb08ny_(da, &e[1], &dwork[lq], &tolq);
	    daxpy_(&nc, &c_b24, &a[1], &c__1, &dwork[lq], &c__1);
	    *res = (d__1 = dwork[idamax_(&nc, &dwork[lq], &c__1) + lq - 1], 
		    abs(d__1));
	    conv = *res < tolq || res0 < 0.;

	    if (! conv) {
		dwork[leta] = dwork[leta] * .5 / dwork[lalpha];

/*              Computation of x  and q . */
/*                              i      i */
/*              DWORK(LETA,...,LETA+DA)   : eta(k,0),...,eta(k,n) */
/*                   (LRO,...,LRO+DA-K+1) : eta(k,n-k+1),...,eta(k,0) */

		for (k = *da; k >= 1; --k) {
		    nck = nc - k + 1;
		    dcopy_(&nck, &dwork[leta], &c_n1, &dwork[lro], &c__1);
		    w = dwork[lambda + k - 1];
		    d__1 = -w;
		    daxpy_(&nck, &d__1, &dwork[lro], &c__1, &dwork[leta], &
			    c__1);
/* L80: */
		}

		s = 0.;

		i__1 = *da;
		for (j = 0; j <= i__1; ++j) {
		    w = (dwork[leta + j] + e[j + 1]) * .5;
		    dwork[lq + j] = w;
/* Computing 2nd power */
		    d__1 = w;
		    s += d__1 * d__1;
/* L100: */
		}

		res0 = s - a0;

/*              Test on the monotonicity of q . */
/*                                           0 */
		conv = dwork[lq] > e[1];
		goto L40;
	    }
	}
    }
/*     END WHILE 40 */

/*     Reverse the order of the coefficients in the array E. */

    dswap_(&nc, &e[1], &c__1, &dwork[1], &c_n1);
    dswap_(&nc, &dwork[1], &c__1, &e[1], &c__1);

    if (! conv) {
	if (hurwtz) {
	    *info = 3;
	} else if (i__ == 1) {
	    *info = 2;
	} else {
	    *info = 4;
	}
    }

    return 0;
/* *** Last line of SB08ND *** */
} /* sb08nd_ */

