/* MB05ND.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb05nd_(integer *n, doublereal *delta, doublereal *a, 
	integer *lda, doublereal *ex, integer *ldex, doublereal *exint, 
	integer *ldexin, doublereal *tol, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ex_dim1, ex_offset, exint_dim1, exint_offset, 
	    i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, l, ij, kk, iq, nn;
    static doublereal eps, err, f2iq1;
    static integer i2iq1;
    static doublereal qmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal delsc;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer jscal;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal small;
    extern /* Subroutine */ int dgesv_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal fnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal fnorm2, coeffd;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal coeffn;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute */

/*     (a)    F(delta) =  exp(A*delta) and */

/*     (b)    H(delta) =  Int[F(s) ds] from s = 0 to s = delta, */

/*     where A is a real N-by-N matrix and delta is a scalar value. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             The scalar value delta of the problem. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix A of the problem. (Array A need not be set if */
/*             DELTA = 0.) */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,N). */

/*     EX      (output) DOUBLE PRECISION array, dimension (LDEX,N) */
/*             The leading N-by-N part of this array contains an */
/*             approximation to F(delta). */

/*     LDEX    INTEGER */
/*             The leading dimension of array EX.  LDEX >= MAX(1,N). */

/*     EXINT   (output) DOUBLE PRECISION array, dimension (LDEXIN,N) */
/*             The leading N-by-N part of this array contains an */
/*             approximation to H(delta). */

/*     LDEXIN  INTEGER */
/*             The leading dimension of array EXINT.  LDEXIN >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the order of the */
/*             Pade approximation to H(t), where t is a scale factor */
/*             determined by the routine. A reasonable value for TOL may */
/*             be SQRT(EPS), where EPS is the machine precision (see */
/*             LAPACK Library routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. LDWORK >= MAX(1,N*(N+1)). */
/*             For optimum performance LDWORK should be larger (2*N*N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, the (i,i) element of the denominator of */
/*                   the Pade approximation is zero, so the denominator */
/*                   is exactly singular; */
/*             = N+1:  if DELTA = (delta * frobenius norm of matrix A) is */
/*                   probably too large to permit meaningful computation. */
/*                   That is, DELTA > SQRT(BIG), where BIG is a */
/*                   representable number near the overflow threshold of */
/*                   the machine (see LAPACK Library Routine DLAMCH). */

/*     METHOD */

/*     This routine uses a Pade approximation to H(t) for some small */
/*     value of t (where 0 < t <= delta) and then calculates F(t) from */
/*     H(t). Finally, the results are re-scaled to give F(delta) and */
/*     H(delta). For a detailed description of the implementation of this */
/*     algorithm see [1]. */

/*     REFERENCES */

/*     [1] Benson, C.J. */
/*         The numerical evaluation of the matrix exponential and its */
/*         integral. */
/*         Report 82/03, Control Systems Research Group, */
/*         School of Electronic Engineering and Computer */
/*         Science, Kingston Polytechnic, January 1982. */

/*     [2] Ward, R.C. */
/*         Numerical computation of the matrix exponential with accuracy */
/*         estimate. */
/*         SIAM J. Numer. Anal., 14, pp. 600-610, 1977. */

/*     [3] Moler, C.B. and Van Loan, C.F. */
/*         Nineteen Dubious Ways to Compute the Exponential of a Matrix. */
/*         SIAM Rev., 20, pp. 801-836, 1978. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine MB05BD by C.J. Benson, Kingston */
/*     Polytechnic, January 1982. */

/*     REVISIONS */

/*     - */

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
    ex_dim1 = *ldex;
    ex_offset = 1 + ex_dim1;
    ex -= ex_offset;
    exint_dim1 = *ldexin;
    exint_offset = 1 + exint_dim1;
    exint -= exint_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    nn = *n * *n;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldex < max(1,*n)) {
	*info = -6;
    } else if (*ldexin < max(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = nn + *n;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -12;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB05ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    dwork[1] = 1.;
    if (*n == 0) {
	return 0;
    }

    dlaset_("Full", n, n, &c_b4, &c_b4, &ex[ex_offset], ldex, (ftnlen)4);
    dlaset_("Full", n, n, &c_b4, &c_b4, &exint[exint_offset], ldexin, (ftnlen)
	    4);

    if (*delta == 0.) {
	dlaset_("Upper", n, n, &c_b4, &c_b11, &ex[ex_offset], ldex, (ftnlen)5)
		;
	return 0;
    }

    if (*n == 1) {
	ex[ex_dim1 + 1] = exp(*delta * a[a_dim1 + 1]);
	if (a[a_dim1 + 1] == 0.) {
	    exint[exint_dim1 + 1] = *delta;
	} else {
	    exint[exint_dim1 + 1] = 1. / a[a_dim1 + 1] * ex[ex_dim1 + 1] - 1. 
		    / a[a_dim1 + 1];
	}
	return 0;
    }

/*     Set some machine parameters. */

    eps = dlamch_("Epsilon", (ftnlen)7);
    small = dlamch_("Safe minimum", (ftnlen)12) / eps;

/*     First calculate the Frobenius norm of A, and the scaling factor. */

    fnorm = *delta * dlange_("Frobenius", n, n, &a[a_offset], lda, &dwork[1], 
	    (ftnlen)9);

    if (fnorm > sqrt(1. / small)) {
	*info = *n + 1;
	return 0;
    }

    jscal = 0;
    delsc = *delta;
/*     WHILE ( FNORM >= HALF ) DO */
L20:
    if (fnorm >= .5) {
	++jscal;
	delsc *= .5;
	fnorm *= .5;
	goto L20;
    }
/*     END WHILE 20 */

/*     Calculate the order of the Pade approximation needed to satisfy */
/*     the requested relative error  TOL. */

/* Computing 2nd power */
    d__1 = fnorm;
    fnorm2 = d__1 * d__1;
    iq = 1;
    qmax = fnorm / 3.;
/* Computing 2nd power */
    d__1 = fnorm2;
    err = *delta / delsc * (d__1 * d__1) / 4.8;
/*     WHILE ( ERR > TOL*( 2*IQ + 3 - FNORM )/1.64 and QMAX >= EPS ) DO */
L40:
    if (err > *tol * ((doublereal) ((iq << 1) + 3) - fnorm) / 1.64) {
	++iq;
	qmax = qmax * (doublereal) (iq + 1) * fnorm / (doublereal) ((iq << 1) 
		* ((iq << 1) + 1));
	if (qmax >= eps) {
/* Computing 2nd power */
	    i__1 = (iq << 1) + 3;
	    err = err * fnorm2 * (doublereal) ((iq << 1) + 5) / (doublereal) (
		    i__1 * i__1 * ((iq << 1) + 4));
	    goto L40;
	}
    }
/*     END WHILE 40 */

/*     Initialise DWORK (to contain succesive powers of A), */
/*                EXINT (to contain the numerator) and */
/*                EX    (to contain the denominator). */

    i2iq1 = (iq << 1) + 1;
    f2iq1 = (doublereal) i2iq1;
    coeffd = -((doublereal) iq) / f2iq1;
    coeffn = .5 / f2iq1;
    ij = 1;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[ij] = delsc * a[i__ + j * a_dim1];
	    exint[i__ + j * exint_dim1] = coeffn * dwork[ij];
	    ex[i__ + j * ex_dim1] = coeffd * dwork[ij];
	    ++ij;
/* L60: */
	}

	exint[j + j * exint_dim1] += 1.;
	ex[j + j * ex_dim1] += 1.;
/* L80: */
    }

    i__1 = iq;
    for (kk = 2; kk <= i__1; ++kk) {

/*        Calculate the next power of  A*DELSC,  and update the numerator */
/*        and denominator. */

	coeffd = -coeffd * (doublereal) (iq + 1 - kk) / (doublereal) (kk * (
		i2iq1 + 1 - kk));
	if (kk % 2 == 0) {
	    coeffn = coeffd / (doublereal) (kk + 1);
	} else {
	    coeffn = -coeffd / (doublereal) (i2iq1 - kk);
	}
	ij = 1;

	if (*ldwork >= nn << 1) {

/*           Enough space for a BLAS 3 calculation. */

	    dgemm_("No transpose", "No transpose", n, n, n, &delsc, &a[
		    a_offset], lda, &dwork[1], n, &c_b4, &dwork[nn + 1], n, (
		    ftnlen)12, (ftnlen)12);
	    dcopy_(&nn, &dwork[nn + 1], &c__1, &dwork[1], &c__1);

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		daxpy_(n, &coeffn, &dwork[ij], &c__1, &exint[j * exint_dim1 + 
			1], &c__1);
		daxpy_(n, &coeffd, &dwork[ij], &c__1, &ex[j * ex_dim1 + 1], &
			c__1);
		ij += *n;
/* L100: */
	    }

	} else {

/*           Not enough space for a BLAS 3 calculation. Use BLAS 2. */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		dgemv_("No transpose", n, n, &c_b11, &a[a_offset], lda, &
			dwork[ij], &c__1, &c_b4, &dwork[nn + 1], &c__1, (
			ftnlen)12);
		dcopy_(n, &dwork[nn + 1], &c__1, &dwork[ij], &c__1);
		dscal_(n, &delsc, &dwork[ij], &c__1);
		daxpy_(n, &coeffn, &dwork[ij], &c__1, &exint[j * exint_dim1 + 
			1], &c__1);
		daxpy_(n, &coeffd, &dwork[ij], &c__1, &ex[j * ex_dim1 + 1], &
			c__1);
		ij += *n;
/* L120: */
	    }

	}
/* L140: */
    }

/*     We now have numerator in EXINT, denominator in EX. */

/*     Solve the set of N systems of linear equations for the columns of */
/*     EXINT  using the LU factorization of EX. */

    dgesv_(n, n, &ex[ex_offset], ldex, &iwork[1], &exint[exint_offset], 
	    ldexin, info);
    if (*info != 0) {
	return 0;
    }

/*     Now we can form EX from EXINT using the formula: */
/*     EX = EXINT * A  +  I */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(n, &delsc, &exint[j * exint_dim1 + 1], &c__1);
/* L160: */
    }

    dgemm_("No transpose", "No transpose", n, n, n, &c_b11, &exint[
	    exint_offset], ldexin, &a[a_offset], lda, &c_b4, &ex[ex_offset], 
	    ldex, (ftnlen)12, (ftnlen)12);

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ex[j + j * ex_dim1] += 1.;
/* L180: */
    }

/*     EX  and  EXINT  have been evaluated at  DELSC,  so the results */
/*     must be re-scaled to give the function values at  DELTA. */

/*     EXINT(2t) = EXINT(t) * ^ EX(t) + I [ */
/*     EX(2t) = EX(t) * EX(t) */

/*     DWORK  is used to accumulate products. */

    i__1 = jscal;
    for (l = 1; l <= i__1; ++l) {
	dlacpy_("Full", n, n, &exint[exint_offset], ldexin, &dwork[1], n, (
		ftnlen)4);
	dgemm_("No transpose", "No transpose", n, n, n, &c_b11, &dwork[1], n, 
		&ex[ex_offset], ldex, &c_b11, &exint[exint_offset], ldexin, (
		ftnlen)12, (ftnlen)12);
	dlacpy_("Full", n, n, &ex[ex_offset], ldex, &dwork[1], n, (ftnlen)4);
	dgemm_("No transpose", "No transpose", n, n, n, &c_b11, &dwork[1], n, 
		&dwork[1], n, &c_b4, &ex[ex_offset], ldex, (ftnlen)12, (
		ftnlen)12);
/* L200: */
    }

    dwork[1] = (doublereal) (nn << 1);
    return 0;
/* *** Last line of MB05ND *** */
} /* mb05nd_ */

