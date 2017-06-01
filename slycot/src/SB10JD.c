/* SB10JD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 1.;
static doublereal c_b10 = 0.;
static integer c__1 = 1;
static doublereal c_b48 = -1.;

/* Subroutine */ int sb10jd_(integer *n, integer *m, integer *np, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *e, integer *
	lde, integer *nsys, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, is, iu, iv, ib2, ic2, ns1, ia12, ia21, isa, lwa;
    static doublereal eps, tol;
    static integer iwrk, info2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer lwamax, minwrk;


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

/*     To convert the descriptor state-space system */

/*     E*dx/dt = A*x + B*u */
/*           y = C*x + D*u */

/*     into regular state-space form */

/*      dx/dt = Ad*x + Bd*u */
/*          y = Cd*x + Dd*u . */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the descriptor system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state matrix A of the descriptor system. */
/*             On exit, the leading NSYS-by-NSYS part of this array */
/*             contains the state matrix Ad of the converted system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the descriptor system. */
/*             On exit, the leading NSYS-by-M part of this array */
/*             contains the input matrix Bd of the converted system. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading NP-by-N part of this array must */
/*             contain the output matrix C of the descriptor system. */
/*             On exit, the leading NP-by-NSYS part of this array */
/*             contains the output matrix Cd of the converted system. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading NP-by-M part of this array must */
/*             contain the matrix D of the descriptor system. */
/*             On exit, the leading NP-by-M part of this array contains */
/*             the matrix Dd of the converted system. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix E of the descriptor system. */
/*             On exit, this array contains no useful information. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= max(1,N). */

/*     NSYS    (output) INTEGER */
/*             The order of the converted state-space system. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max( 1, 2*N*N + 2*N + N*MAX( 5, N + M + NP ) ). */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the iteration for computing singular value */
/*                   decomposition did not converge. */

/*     METHOD */

/*     The routine performs the transformations described in [1]. */

/*     REFERENCES */

/*     [1] Chiang, R.Y. and Safonov, M.G. */
/*         Robust Control Toolbox User's Guide. */
/*         The MathWorks Inc., Natick, Mass., 1992. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000, */
/*     Feb. 2001. */

/*     KEYWORDS */

/*     Descriptor systems, state-space models. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*np < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldc < max(1,*np)) {
	*info = -9;
    } else if (*ldd < max(1,*np)) {
	*info = -11;
    } else if (*lde < max(1,*n)) {
	*info = -13;
    }

/*     Compute workspace. */

/* Computing MAX */
/* Computing MAX */
    i__3 = 5, i__4 = *n + *m + *np;
    i__1 = 1, i__2 = (*n << 1) * (*n + 1) + *n * max(i__3,i__4);
    minwrk = max(i__1,i__2);
    if (*ldwork < minwrk) {
	*info = -16;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10JD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*nsys = 0;
	dwork[1] = 1.;
	return 0;
    }

/*     Set tol. */

    eps = dlamch_("Epsilon", (ftnlen)7);
    tol = sqrt(eps);

/*     Workspace usage. */

    is = 0;
    iu = is + *n;
    iv = iu + *n * *n;

    iwrk = iv + *n * *n;

/*     Compute the SVD of E. */
/*     Additional workspace:  need   5*N; prefer larger. */

    i__1 = *ldwork - iwrk;
    dgesvd_("S", "S", n, n, &e[e_offset], lde, &dwork[is + 1], &dwork[iu + 1],
	     n, &dwork[iv + 1], n, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1,
	     (ftnlen)1);
    if (info2 != 0) {
	*info = 1;
	return 0;
    }
/* Computing MAX */
    i__1 = minwrk, i__2 = (integer) (dwork[iwrk + 1] + iwrk);
    lwamax = max(i__1,i__2);

/*     Determine the rank of E. */

    ns1 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dwork[is + i__] > tol) {
	    ++ns1;
	}
/* L10: */
    }
    if (ns1 > 0) {

/*        Transform A. */
/*        Additional workspace:  need   N*max(N,M,NP). */

	dgemm_("T", "N", n, n, n, &c_b9, &dwork[iu + 1], n, &a[a_offset], lda,
		 &c_b10, &dwork[iwrk + 1], n, (ftnlen)1, (ftnlen)1);
	dgemm_("N", "T", n, n, n, &c_b9, &dwork[iwrk + 1], n, &dwork[iv + 1], 
		n, &c_b10, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);

/*        Transform B. */

	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[iwrk + 1], n, (ftnlen)
		4);
	dgemm_("T", "N", n, m, n, &c_b9, &dwork[iu + 1], n, &dwork[iwrk + 1], 
		n, &c_b10, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1);

/*        Transform C. */

	dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwrk + 1], np, (
		ftnlen)4);
	dgemm_("N", "T", np, n, n, &c_b9, &dwork[iwrk + 1], np, &dwork[iv + 1]
		, n, &c_b10, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);

	k = *n - ns1;
	if (k > 0) {
	    isa = iu + k * k;
	    iv = isa + k;
	    iwrk = iv + k * max(k,ns1);

/*           Compute the SVD of A22. */
/*           Additional workspace:  need   5*K; prefer larger. */

	    i__1 = *ldwork - iwrk;
	    dgesvd_("S", "S", &k, &k, &a[ns1 + 1 + (ns1 + 1) * a_dim1], lda, &
		    dwork[isa + 1], &dwork[iu + 1], &k, &dwork[iv + 1], &k, &
		    dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (ftnlen)1);
	    if (info2 != 0) {
		*info = 1;
		return 0;
	    }
	    ia12 = iwrk;
	    ib2 = ia12 + ns1 * k;
	    ic2 = ib2 + k * *m;

	    lwa = (integer) dwork[iwrk + 1] + iwrk;
/* Computing MAX */
	    i__1 = max(lwa,lwamax), i__2 = ic2 + k * *np;
	    lwamax = max(i__1,i__2);

/*           Compute the transformed A12. */

	    dgemm_("N", "T", &ns1, &k, &k, &c_b9, &a[(ns1 + 1) * a_dim1 + 1], 
		    lda, &dwork[iv + 1], &k, &c_b10, &dwork[ia12 + 1], &ns1, (
		    ftnlen)1, (ftnlen)1);

/*           Compute CC2. */

	    dgemm_("N", "T", np, &k, &k, &c_b9, &c__[(ns1 + 1) * c_dim1 + 1], 
		    ldc, &dwork[iv + 1], &k, &c_b10, &dwork[ic2 + 1], np, (
		    ftnlen)1, (ftnlen)1);

/*           Compute the transformed A21. */

	    ia21 = iv;
	    dgemm_("T", "N", &k, &ns1, &k, &c_b9, &dwork[iu + 1], &k, &a[ns1 
		    + 1 + a_dim1], lda, &c_b10, &dwork[ia21 + 1], &k, (ftnlen)
		    1, (ftnlen)1);

/*           Compute BB2. */

	    dgemm_("T", "N", &k, m, &k, &c_b9, &dwork[iu + 1], &k, &b[ns1 + 1 
		    + b_dim1], ldb, &c_b10, &dwork[ib2 + 1], &k, (ftnlen)1, (
		    ftnlen)1);

/*           Compute A12*pinv(A22) and CC2*pinv(A22). */

	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		scale = 0.;
		if (dwork[isa + j] > tol) {
		    scale = 1. / dwork[isa + j];
		}
		dscal_(&ns1, &scale, &dwork[ia12 + (j - 1) * ns1 + 1], &c__1);
		dscal_(np, &scale, &dwork[ic2 + (j - 1) * *np + 1], &c__1);
/* L20: */
	    }

/*           Compute Ad. */

	    dgemm_("N", "N", &ns1, &ns1, &k, &c_b48, &dwork[ia12 + 1], &ns1, &
		    dwork[ia21 + 1], &k, &c_b9, &a[a_offset], lda, (ftnlen)1, 
		    (ftnlen)1);

/*           Compute Bd. */

	    dgemm_("N", "N", &ns1, m, &k, &c_b48, &dwork[ia12 + 1], &ns1, &
		    dwork[ib2 + 1], &k, &c_b9, &b[b_offset], ldb, (ftnlen)1, (
		    ftnlen)1);

/*           Compute Cd. */

	    dgemm_("N", "N", np, &ns1, &k, &c_b48, &dwork[ic2 + 1], np, &
		    dwork[ia21 + 1], &k, &c_b9, &c__[c_offset], ldc, (ftnlen)
		    1, (ftnlen)1);

/*           Compute Dd. */

	    dgemm_("N", "N", np, m, &k, &c_b48, &dwork[ic2 + 1], np, &dwork[
		    ib2 + 1], &k, &c_b9, &d__[d_offset], ldd, (ftnlen)1, (
		    ftnlen)1);
	}
	i__1 = ns1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scale = 1. / sqrt(dwork[is + i__]);
	    dscal_(&ns1, &scale, &a[i__ + a_dim1], lda);
	    dscal_(m, &scale, &b[i__ + b_dim1], ldb);
/* L30: */
	}
	i__1 = ns1;
	for (j = 1; j <= i__1; ++j) {
	    scale = 1. / sqrt(dwork[is + j]);
	    dscal_(&ns1, &scale, &a[j * a_dim1 + 1], &c__1);
	    dscal_(np, &scale, &c__[j * c_dim1 + 1], &c__1);
/* L40: */
	}
	*nsys = ns1;
    } else {
	d__1 = -1. / eps;
	dlaset_("F", n, n, &c_b10, &d__1, &a[a_offset], lda, (ftnlen)1);
	dlaset_("F", n, m, &c_b10, &c_b10, &b[b_offset], ldb, (ftnlen)1);
	dlaset_("F", np, n, &c_b10, &c_b10, &c__[c_offset], ldc, (ftnlen)1);
	*nsys = *n;
    }
    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10JD *** */
} /* sb10jd_ */

