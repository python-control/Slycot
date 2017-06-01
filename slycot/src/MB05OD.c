/* MB05OD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b25 = 2.;

/* Subroutine */ int mb05od_(char *balanc, integer *n, integer *ndiag, 
	doublereal *delta, doublereal *a, integer *lda, integer *mdig, 
	integer *idig, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen balanc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_lg10(doublereal *), sqrt(doublereal), exp(doublereal), log(
	    doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal p, s, u, bd;
    static integer ij, ik;
    static doublereal fn, gn, ss, tr, xn;
    static integer im1;
    static doublereal sd2, big, eps, var, tmp1;
    static integer ndec, base;
    static doublereal eabs, rerl, temp, rerr, size;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal sum2d;
    extern /* Subroutine */ int mb04md_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    static integer ifail;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical lbals;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal avgev, small;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb05oy_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ndagm1, ndagm2, ndecm1, jwora1, jwora2, jwora3;
    static doublereal ovrth2;
    static char actbal[1];
    static integer jworv1, jworv2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal eavgev, factor;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal maxred;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal underf;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal emnorm, vareps;
    static integer mpower;
    static doublereal ovrthr;


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

/*     To compute exp(A*delta) where A is a real N-by-N matrix and delta */
/*     is a scalar value. The routine also returns the minimal number of */
/*     accurate digits in the 1-norm of exp(A*delta) and the number of */
/*     accurate digits in the 1-norm of exp(A*delta) at 95% confidence */
/*     level. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Specifies whether or not a balancing transformation (done */
/*             by SLICOT Library routine MB04MD) is required, as follows: */
/*             = 'N', do not use balancing; */
/*             = 'S', use balancing (scaling). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     NDIAG   (input) INTEGER */
/*             The specified order of the diagonal Pade approximant. */
/*             In the absence of further information NDIAG should */
/*             be set to 9.  NDIAG should not exceed 15.  NDIAG >= 1. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             The scalar value delta of the problem. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On input, the leading N-by-N part of this array must */
/*             contain the matrix A of the problem. (This is not needed */
/*             if DELTA = 0.) */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains the solution matrix exp(A*delta). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     MDIG    (output) INTEGER */
/*             The minimal number of accurate digits in the 1-norm of */
/*             exp(A*delta). */

/*     IDIG    (output) INTEGER */
/*             The number of accurate digits in the 1-norm of */
/*             exp(A*delta) at 95% confidence level. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= N*(2*N+NDIAG+1)+NDIAG, if N >  1. */
/*             LDWORK >= 1,                     if N <= 1. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  if MDIG = 0 and IDIG > 0, warning for possible */
/*                   inaccuracy (the exponential has been computed); */
/*             = 2:  if MDIG = 0 and IDIG = 0, warning for severe */
/*                   inaccuracy (the exponential has been computed); */
/*             = 3:  if balancing has been requested, but it failed to */
/*                   reduce the matrix norm and was not actually used. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the norm of matrix A*delta (after a possible */
/*                   balancing) is too large to obtain an accurate */
/*                   result; */
/*             = 2:  if the coefficient matrix (the denominator of the */
/*                   Pade approximant) is exactly singular; try a */
/*                   different value of NDIAG; */
/*             = 3:  if the solution exponential would overflow, possibly */
/*                   due to a too large value DELTA; the calculations */
/*                   stopped prematurely. This error is not likely to */
/*                   appear. */

/*     METHOD */

/*     The exponential of the matrix A is evaluated from a diagonal Pade */
/*     approximant. This routine is a modification of the subroutine */
/*     PADE, described in reference [1]. The routine implements an */
/*     algorithm which exploits the identity */

/*         (exp[(2**-m)*A]) ** (2**m) = exp(A), */

/*     where m is an integer determined by the algorithm, to improve the */
/*     accuracy for matrices with large norms. */

/*     REFERENCES */

/*     [1] Ward, R.C. */
/*         Numerical computation of the matrix exponential with accuracy */
/*         estimate. */
/*         SIAM J. Numer. Anal., 14, pp. 600-610, 1977. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05CD by T.W.C. Williams, Kingston */
/*     Polytechnic, March 1982. */

/*     REVISIONS */

/*     June 14, 1997, April 25, 2003, December 12, 2004. */

/*     KEYWORDS */

/*     Continuous-time system, matrix algebra, matrix exponential, */
/*     matrix operations, Pade approximation. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *iwarn = 0;
    *info = 0;
    lbals = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lbals)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ndiag < 1) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldwork < 1 || *ldwork < *n * ((*n << 1) + *ndiag + 1) + *
	    ndiag && *n > 1) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB05OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    eps = dlamch_("Epsilon", (ftnlen)7);
    d__1 = 1. / eps;
    ndec = (integer) (d_lg10(&d__1) + 1.);

    if (*n == 0) {
	*mdig = ndec;
	*idig = ndec;
	return 0;
    }

/*     Set some machine parameters. */

    base = (integer) dlamch_("Base", (ftnlen)4);
    ndecm1 = ndec - 1;
    underf = dlamch_("Underflow", (ftnlen)9);
    ovrthr = dlamch_("Overflow", (ftnlen)8);
    ovrth2 = sqrt(ovrthr);

    if (*delta == 0.) {

/*        The DELTA = 0 case. */

	dlaset_("Full", n, n, &c_b10, &c_b11, &a[a_offset], lda, (ftnlen)4);
	*mdig = ndecm1;
	*idig = ndecm1;
	return 0;
    }

    if (*n == 1) {

/*        The 1-by-1 case. */

	a[a_dim1 + 1] = exp(a[a_dim1 + 1] * *delta);
	*mdig = ndecm1;
	*idig = ndecm1;
	return 0;
    }

/*     Set pointers for the workspace. */

    jwora1 = 1;
    jwora2 = jwora1 + *n * *n;
    jwora3 = jwora2 + *n * *ndiag;
    jworv1 = jwora3 + *n * *n;
    jworv2 = jworv1 + *n;

/*     Compute Pade coefficients in DWORK(JWORV2:JWORV2+NDIAG-1). */

    dwork[jworv2] = .5;

    i__1 = *ndiag;
    for (i__ = 2; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	dwork[jworv2 + im1] = dwork[jworv2 + i__ - 2] * (doublereal) (*ndiag 
		- im1) / (doublereal) (i__ * ((*ndiag << 1) - im1));
/* L20: */
    }

/* Computing 2nd power */
    d__1 = eps;
/* Computing 2nd power */
    d__2 = (doublereal) base;
    vareps = d__1 * d__1 * ((d__2 * d__2 - 1.) / (log((doublereal) base) * 
	    24.));
    xn = (doublereal) (*n);
    tr = 0.;

/*     Apply a translation with the mean of the eigenvalues of A*DELTA. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dscal_(n, delta, &a[i__ * a_dim1 + 1], &c__1);
	tr += a[i__ + i__ * a_dim1];
/* L40: */
    }

    avgev = tr / xn;
    if (avgev > log(ovrthr) || avgev < log(underf)) {
	avgev = 0.;
    }
    if (avgev != 0.) {
	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
		ftnlen)6);

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    a[i__ + i__ * a_dim1] -= avgev;
/* L60: */
	}

	temp = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
		ftnlen)6);
	if (temp > anorm * .5) {

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + i__ * a_dim1] += avgev;
/* L80: */
	    }

	    avgev = 0.;
	}
    }
    *(unsigned char *)actbal = *(unsigned char *)balanc;
    if (lbals) {

/*        Balancing (scaling) has been requested.  First, save A. */

	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[jwora1], n, (ftnlen)4)
		;
	maxred = 200.;
	mb04md_(n, &maxred, &a[a_offset], lda, &dwork[jworv1], info);
	if (maxred < 1.) {

/*           Recover the matrix and reset DWORK(JWORV1,...,JWORV1+N-1) */
/*           to 1, as no reduction of the norm occured (unlikely event). */

	    dlacpy_("Full", n, n, &dwork[jwora1], n, &a[a_offset], lda, (
		    ftnlen)4);
	    *(unsigned char *)actbal = 'N';
	    dwork[jworv1] = 1.;
	    i__1 = *n - 1;
	    dcopy_(&i__1, &dwork[jworv1], &c__0, &dwork[jworv1 + 1], &c__1);
	    *iwarn = 3;
	}
    }

/*     Scale the matrix by 2**(-M), where M is the minimum integer */
/*     so that the resulted matrix has the 1-norm less than 0.5. */

    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);
    m = 0;
    if (anorm >= .5) {
	mpower = (integer) (log(ovrthr) / log(2.));
	m = (integer) (log(anorm) / log(2.)) + 1;
	if (m > mpower) {

/*           Error return: The norm of A*DELTA is too large. */

	    *info = 1;
	    return 0;
	}
	factor = pow_di(&c_b25, &m);
	if (m + 1 < mpower) {
	    ++m;
	    factor *= 2.;
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = 1. / factor;
	    dscal_(n, &d__1, &a[i__ * a_dim1 + 1], &c__1);
/* L120: */
	}

    }
    ndagm1 = *ndiag - 1;
    ndagm2 = ndagm1 - 1;
    ij = 0;

/*     Compute the factors of the diagonal Pade approximant. */
/*     The loop 200 takes the accuracy requirements into account: */
/*     Pade coefficients decrease with K, so the calculations should */
/*     be performed in backward order, one column at a time. */
/*     (A BLAS 3 implementation in forward order, using DGEMM, could */
/*     possibly be less accurate.) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dgemv_("No transpose", n, n, &c_b11, &a[a_offset], lda, &a[j * a_dim1 
		+ 1], &c__1, &c_b10, &dwork[jwora2], &c__1, (ftnlen)12);
	ik = 0;

	i__2 = ndagm2;
	for (k = 1; k <= i__2; ++k) {
	    dgemv_("No transpose", n, n, &c_b11, &a[a_offset], lda, &dwork[
		    jwora2 + ik], &c__1, &c_b10, &dwork[jwora2 + ik + *n], &
		    c__1, (ftnlen)12);
	    ik += *n;
/* L140: */
	}

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = 0.;
	    u = 0.;
	    ik = ndagm2 * *n + i__ - 1;

	    for (k = ndagm1; k >= 1; --k) {
		p = dwork[jworv2 + k] * dwork[jwora2 + ik];
		ik -= *n;
		s += p;
		if ((k + 1) % 2 == 0) {
		    u += p;
		} else {
		    u -= p;
		}
/* L160: */
	    }

	    p = dwork[jworv2] * a[i__ + j * a_dim1];
	    s += p;
	    u -= p;
	    if (i__ == j) {
		s += 1.;
		u += 1.;
	    }
	    dwork[jwora3 + ij] = s;
	    dwork[jwora1 + ij] = u;
	    ++ij;
/* L180: */
	}

/* L200: */
    }

/*     Compute the exponential of the scaled matrix, using diagonal Pade */
/*     approximants.  As, in theory [1], the denominator of the Pade */
/*     approximant should be very well conditioned, no condition estimate */
/*     is computed. */

    dgetrf_(n, n, &dwork[jwora1], n, &iwork[1], &ifail);
    if (ifail > 0) {

/*        Error return: The matrix is exactly singular. */

	*info = 2;
	return 0;
    }

    dlacpy_("Full", n, n, &dwork[jwora3], n, &a[a_offset], lda, (ftnlen)4);
    dgetrs_("No transpose", n, n, &dwork[jwora1], n, &iwork[1], &a[a_offset], 
	    lda, &ifail, (ftnlen)12);

/*     Prepare for the calculation of the accuracy estimates. */
/*     Note that ANORM here is in the range [1, e]. */

    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);
    if (anorm >= 1.) {
	eabs = (xn * 19. + 47.) * (eps * anorm);
    } else {
	eabs = (xn * 19. + 47.) * eps * anorm;
    }
    if (m != 0) {
	var = xn * vareps;
	fn = xn * 4. / ((xn + 2.) * (xn + 1.));
/* Computing 2nd power */
	d__1 = xn + 2.;
/* Computing 2nd power */
	d__2 = xn + 1.;
	gn = ((xn * 2. + 10.) * xn - 4.) / (d__1 * d__1 * (d__2 * d__2));

/*        Square-up the computed exponential matrix M times, with caution */
/*        for avoiding overflows. */

	i__1 = m;
	for (k = 1; k <= i__1; ++k) {
	    if (anorm > ovrth2) {

/*              The solution could overflow. */

		d__1 = 1. / anorm;
		dgemm_("No transpose", "No transpose", n, n, n, &d__1, &a[
			a_offset], lda, &a[a_offset], lda, &c_b10, &dwork[
			jwora1], n, (ftnlen)12, (ftnlen)12);
		s = dlange_("1-norm", n, n, &dwork[jwora1], n, &dwork[jwora1],
			 (ftnlen)6);
		if (anorm <= ovrthr / s) {
		    dlascl_("General", n, n, &c_b11, &anorm, n, n, &dwork[
			    jwora1], n, info, (ftnlen)7);
		    temp = ovrthr;
		} else {

/*                 Error return: The solution would overflow. */
/*                 This will not happen on most machines, due to the */
/*                 selection of M. */

		    *info = 3;
		    return 0;
		}
	    } else {
		dgemm_("No transpose", "No transpose", n, n, n, &c_b11, &a[
			a_offset], lda, &a[a_offset], lda, &c_b10, &dwork[
			jwora1], n, (ftnlen)12, (ftnlen)12);
/* Computing 2nd power */
		d__1 = anorm;
		temp = d__1 * d__1;
	    }
	    if (eabs < 1.) {
		eabs = (anorm * 2. + eabs) * eabs + xn * (eps * temp);
	    } else if (eabs < sqrt(1. - xn * eps + ovrthr / temp) * anorm - 
		    anorm) {
/* Computing 2nd power */
		d__1 = eabs;
		eabs = xn * (eps * temp) + anorm * eabs * 2. + d__1 * d__1;
	    } else {
		eabs = ovrthr;
	    }

	    tmp1 = fn * var + gn * (temp * vareps);
	    if (tmp1 > ovrthr / temp) {
		var = ovrthr;
	    } else {
		var = tmp1 * temp;
	    }

	    dlacpy_("Full", n, n, &dwork[jwora1], n, &a[a_offset], lda, (
		    ftnlen)4);
	    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1],
		     (ftnlen)6);
/* L220: */
	}

    } else {
	var = xn * 12. * vareps;
    }

/*     Apply back transformations, if balancing was effectively used. */

    mb05oy_(actbal, n, &c__1, n, &a[a_offset], lda, &dwork[jworv1], info, (
	    ftnlen)1);
    eavgev = exp(avgev);
    emnorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);

/*     Compute auxiliary quantities needed for the accuracy estimates. */

    big = 1.;
    small = 1.;
    if (lbals) {

/*        Compute norms of the diagonal scaling matrix and its inverse. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = dwork[jworv1 + i__ - 1];
	    if (big < u) {
		big = u;
	    }
	    if (small > u) {
		small = u;
	    }
/* L240: */
	}

	sum2d = dnrm2_(n, &dwork[jworv1], &c__1);
    } else {
	sum2d = sqrt(xn);
    }

/*     Update the exponential for the initial translation, and update the */
/*     auxiliary quantities needed for the accuracy estimates. */

    sd2 = sqrt(xn * 8. * vareps) * anorm;
    bd = sqrt(var);
    ss = max(bd,sd2);
    bd = min(bd,sd2);
/* Computing 2nd power */
    d__1 = bd / ss;
    sd2 = ss * sqrt(d__1 * d__1 + 1.);
    if (sd2 <= 1.) {
	sd2 = 2. / xn * sum2d * sd2;
    } else if (sum2d / xn < ovrthr / 2. / sd2) {
	sd2 = 2. / xn * sum2d * sd2;
    } else {
	sd2 = ovrthr;
    }
    if (lbals) {
	size = 0.;
    } else {
	if (sd2 < ovrthr - emnorm) {
	    size = emnorm + sd2;
	} else {
	    size = ovrthr;
	}
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ss = dasum_(n, &a[j * a_dim1 + 1], &c__1);
	dscal_(n, &eavgev, &a[j * a_dim1 + 1], &c__1);
	if (lbals) {
	    bd = dwork[jworv1 + j - 1];
/* Computing MAX */
	    d__1 = size, d__2 = ss + sd2 / bd;
	    size = max(d__1,d__2);
	}
/* L260: */
    }

/*     Set the accuracy estimates and warning errors, if any. */

    rerr = d_lg10(&big) + d_lg10(&eabs) - d_lg10(&small) - d_lg10(&emnorm) - 
	    d_lg10(&eps);
    if (size > emnorm) {
	d__1 = (size / emnorm - 1.) / eps;
	rerl = d_lg10(&d__1);
    } else {
	rerl = 0.;
    }
/* Computing MIN */
    i__1 = ndec - (integer) (rerr + .5);
    *mdig = min(i__1,ndecm1);
/* Computing MIN */
    i__1 = ndec - (integer) (rerl + .5);
    *idig = min(i__1,ndecm1);

    if (*mdig <= 0) {
	*mdig = 0;
	*iwarn = 1;
    }
    if (*idig <= 0) {
	*idig = 0;
	*iwarn = 2;
    }

    return 0;
/* *** Last line of MB05OD *** */
} /* mb05od_ */

