/* FB01VD.f -- translated by f2c (version 20100827).
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

static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b15 = 2.;
static doublereal c_b26 = 0.;
static doublereal c_b40 = -1.;

/* Subroutine */ int fb01vd_(integer *n, integer *m, integer *l, doublereal *
	p, integer *ldp, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *q, integer *ldq, 
	doublereal *r__, integer *ldr, doublereal *k, integer *ldk, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, k_dim1, 
	    k_offset, p_dim1, p_offset, q_dim1, q_offset, r_dim1, r_offset, 
	    i__1, i__2;

    /* Local variables */
    static integer j, n1, ldw;
    extern /* Subroutine */ int mb01rd_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dscal_(integer *, doublereal *, 
	    doublereal *, integer *), dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrsm_(char *, char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    static doublereal rnorm;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpocon_(char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);


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

/*     To compute one recursion of the conventional Kalman filter */
/*     equations. This is one update of the Riccati difference equation */
/*     and the Kalman filter gain. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices P      and A .  N >= 0. */
/*                       i|i-1      i */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e., the order of the matrix */
/*             Q .  M >= 0. */
/*              i */

/*     L       (input) INTEGER */
/*             The actual output dimension, i.e., the order of the matrix */
/*             R .  L >= 0. */
/*              i */

/*     P       (input/output) DOUBLE PRECISION array, dimension (LDP,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain P     , the state covariance matrix at instant */
/*                      i|i-1 */
/*             (i-1). The upper triangular part only is needed. */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains P     , the state covariance matrix at */
/*                             i+1|i */
/*             instant i. The strictly lower triangular part is not set. */
/*             Otherwise, the leading N-by-N part of this array contains */
/*             P     , its input value. */
/*              i|i-1 */

/*     LDP     INTEGER */
/*             The leading dimension of array P.  LDP >= MAX(1,N). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain A , */
/*                                                                 i */
/*             the state transition matrix of the discrete system at */
/*             instant i. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain B , */
/*                                                                 i */
/*             the input weight matrix of the discrete system at */
/*             instant i. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array must contain C , */
/*                                                                 i */
/*             the output weight matrix of the discrete system at */
/*             instant i. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,L). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,M) */
/*             The leading M-by-M part of this array must contain Q , */
/*                                                                 i */
/*             the input (process) noise covariance matrix at instant i. */
/*             The diagonal elements of this array are modified by the */
/*             routine, but are restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,M). */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,L) */
/*             On entry, the leading L-by-L part of this array must */
/*             contain R , the output (measurement) noise covariance */
/*                      i */
/*             matrix at instant i. */
/*             On exit, if INFO = 0, or INFO = L+1, the leading L-by-L */
/*                                                                  1/2 */
/*             upper triangular part of this array contains (RINOV )   , */
/*                                                                i */
/*             the square root (left Cholesky factor) of the covariance */
/*             matrix of the innovations at instant i. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,L). */

/*     K       (output) DOUBLE PRECISION array, dimension (LDK,L) */
/*             If INFO = 0, the leading N-by-L part of this array */
/*             contains K , the Kalman filter gain matrix at instant i. */
/*                       i */
/*             If INFO > 0, the leading N-by-L part of this array */
/*             contains the matrix product P     C'. */
/*                                          i|i-1 i */

/*     LDK     INTEGER */
/*             The leading dimension of array K.  LDK >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the matrix RINOV . If the user sets TOL > 0, then the */
/*                             i */
/*             given value of TOL is used as a lower bound for the */
/*             reciprocal condition number of that matrix; a matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be nonsingular. If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = L*L*EPS, is used instead, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, or INFO = L+1, DWORK(1) returns an */
/*             estimate of the reciprocal of the condition number (in the */
/*             1-norm) of the matrix RINOV . */
/*                                        i */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,L*N+3*L,N*N,N*M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -k, the k-th argument had an illegal */
/*                   value; */
/*             = k:  if INFO = k, 1 <= k <= L, the leading minor of order */
/*                   k of the matrix RINOV  is not positive-definite, and */
/*                                        i */
/*                   its Cholesky factorization could not be completed; */
/*             = L+1: the matrix RINOV  is singular, i.e., the condition */
/*                                    i */
/*                   number estimate of RINOV  (in the 1-norm) exceeds */
/*                                           i */
/*                   1/TOL. */

/*     METHOD */

/*     The conventional Kalman filter gain used at the i-th recursion */
/*     step is of the form */

/*                            -1 */
/*        K  = P     C'  RINOV  , */
/*         i    i|i-1 i       i */

/*     where RINOV  = C P     C' + R , and the state covariance matrix */
/*                i    i i|i-1 i    i */

/*     P      is updated by the discrete-time difference Riccati equation */
/*      i|i-1 */

/*        P      = A  (P      - K C P     ) A'  + B Q B'. */
/*         i+1|i    i   i|i-1    i i i|i-1   i     i i i */

/*     Using these two updates, the combined time and measurement update */
/*     of the state X      is given by */
/*                   i|i-1 */

/*        X      = A X      + A K (Y  - C X     ), */
/*         i+1|i    i i|i-1    i i  i    i i|i-1 */

/*     where Y  is the new observation at step i. */
/*            i */

/*     REFERENCES */

/*     [1] Anderson, B.D.O. and Moore, J.B. */
/*         Optimal Filtering, */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1979. */

/*     [2] Verhaegen, M.H.G. and Van Dooren, P. */
/*         Numerical Aspects of Different Kalman Filter Implementations. */
/*         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, 1986. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*             3   2 */
/*      3/2 x N + N  x (3 x L + M/2) */

/*     operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine FB01JD by M.H.G. Verhaegen, */
/*     M. Vanbegin, and P. Van Dooren. */

/*     REVISIONS */

/*     February 20, 1998, November 20, 2003, April 20, 2004. */

/*     KEYWORDS */

/*     Kalman filtering, optimal filtering, recursive estimation. */

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
    p_dim1 = *ldp;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    k_dim1 = *ldk;
    k_offset = 1 + k_dim1;
    k -= k_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    n1 = max(1,*n);
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*ldp < n1) {
	*info = -5;
    } else if (*lda < n1) {
	*info = -7;
    } else if (*ldb < n1) {
	*info = -9;
    } else if (*ldc < max(1,*l)) {
	*info = -11;
    } else if (*ldq < max(1,*m)) {
	*info = -13;
    } else if (*ldr < max(1,*l)) {
	*info = -15;
    } else if (*ldk < n1) {
	*info = -17;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *l * *n + *l * 3, i__1 = max(i__1,i__2), i__2 = *n * 
		*n, i__1 = max(i__1,i__2), i__2 = *n * *m;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -21;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("FB01VD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (max(*n,*l) == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Efficiently compute RINOV = CPC' + R in R and put CP in DWORK and */
/*     PC' in K. (The content of DWORK on exit from MB01RD is used.) */
/*     Workspace: need L*N. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code.) */

    mb01rd_("Upper", "No transpose", l, n, &c_b5, &c_b5, &r__[r_offset], ldr, 
	    &c__[c_offset], ldc, &p[p_offset], ldp, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);
    ldw = max(1,*l);

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	dcopy_(n, &dwork[j], &ldw, &k[j * k_dim1 + 1], &c__1);
/* L10: */
    }

    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[1], &ldw, (ftnlen)4);
    dtrmm_("Right", "Upper", "Transpose", "Non-unit", l, n, &c_b5, &p[
	    p_offset], ldp, &dwork[1], &ldw, (ftnlen)5, (ftnlen)5, (ftnlen)9, 
	    (ftnlen)8);
    i__1 = *ldp + 1;
    dscal_(n, &c_b15, &p[p_offset], &i__1);

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	daxpy_(n, &c_b5, &k[j * k_dim1 + 1], &c__1, &dwork[j], &ldw);
	dcopy_(n, &dwork[j], &ldw, &k[j * k_dim1 + 1], &c__1);
/* L20: */
    }

/*     Calculate the Cholesky decomposition U'U of the innovation */
/*     covariance matrix RINOV, and its reciprocal condition number. */
/*     Workspace: need L*N + 3*L. */

    jwork = *l * *n + 1;
    rnorm = dlansy_("1-norm", "Upper", l, &r__[r_offset], ldr, &dwork[jwork], 
	    (ftnlen)6, (ftnlen)5);

    toldef = *tol;
    if (toldef <= 0.) {
	toldef = (doublereal) (*l * *l) * dlamch_("Epsilon", (ftnlen)7);
    }
    dpotrf_("Upper", l, &r__[r_offset], ldr, info, (ftnlen)5);
    if (*info != 0) {
	return 0;
    }

    dpocon_("Upper", l, &r__[r_offset], ldr, &rnorm, &rcond, &dwork[jwork], &
	    iwork[1], info, (ftnlen)5);

    if (rcond < toldef) {

/*        Error return: RINOV is numerically singular. */

	*info = *l + 1;
	dwork[1] = rcond;
	return 0;
    }

    if (*l > 1) {
	i__1 = *l - 1;
	i__2 = *l - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b26, &c_b26, &r__[r_dim1 + 2], ldr, 
		(ftnlen)5);
    }
/*                                                          -1 */
/*     Calculate the Kalman filter gain matrix  K = PC'RINOV . */
/*     Workspace: need L*N. */

    dtrsm_("Right", "Upper", "No transpose", "Non-unit", n, l, &c_b5, &r__[
	    r_offset], ldr, &k[k_offset], ldk, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
    dtrsm_("Right", "Upper", "Transpose", "Non-unit", n, l, &c_b5, &r__[
	    r_offset], ldr, &k[k_offset], ldk, (ftnlen)5, (ftnlen)5, (ftnlen)
	    9, (ftnlen)8);

/*     First part of the Riccati equation update: compute A(P-KCP)A'. */
/*     The upper triangular part of the symmetric matrix P-KCP is formed. */
/*     Workspace: need max(L*N,N*N). */

    jwork = 1;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dgemv_("No transpose", &j, l, &c_b40, &k[k_offset], ldk, &dwork[jwork]
		, &c__1, &c_b5, &p[j * p_dim1 + 1], &c__1, (ftnlen)12);
	jwork += *l;
/* L30: */
    }

    mb01rd_("Upper", "No transpose", n, n, &c_b26, &c_b5, &p[p_offset], ldp, &
	    a[a_offset], lda, &p[p_offset], ldp, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);

/*     Second part of the Riccati equation update: add BQB'. */
/*     Workspace: need N*M. */

    mb01rd_("Upper", "No transpose", n, m, &c_b5, &c_b5, &p[p_offset], ldp, &
	    b[b_offset], ldb, &q[q_offset], ldq, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);
    i__1 = *ldq + 1;
    dscal_(m, &c_b15, &q[q_offset], &i__1);

/*     Set the reciprocal of the condition number estimate. */

    dwork[1] = rcond;

    return 0;
/* *** Last line of FB01VD *** */
} /* fb01vd_ */

