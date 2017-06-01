/* TB01VD.f -- translated by f2c (version 20100827).
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
static doublereal c_b18 = 1.;
static doublereal c_b42 = 0.;
static doublereal c_b93 = -.5;
static doublereal c_b96 = -1.;

/* Subroutine */ int tb01vd_(char *apply, integer *n, integer *m, integer *l, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *x0, 
	doublereal *theta, integer *ltheta, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen apply_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), tan(doublereal);

    /* Local variables */
    static integer i__, j, k, ca, ia, in, iq;
    static doublereal ri;
    static integer ir;
    static doublereal ti;
    static integer it, iu, ldt, iwi, iwr, ldca;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal piby2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb03od_(char *, char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dtrmm_(char *, 
	    char *, char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrsm_(char *, char *, char *
	    , char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lapply;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer wrkopt;


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

/*     To convert the linear discrete-time system given as (A, B, C, D), */
/*     with initial state x0, into the output normal form [1], with */
/*     parameter vector THETA. The matrix A is assumed to be stable. */
/*     The matrices A, B, C, D and the vector x0 are converted, so that */
/*     on exit they correspond to the system defined by THETA. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     APPLY   CHARACTER*1 */
/*             Specifies whether or not the parameter vector should be */
/*             transformed using a bijective mapping, as follows: */
/*             = 'A' : apply the bijective mapping to the N vectors in */
/*                     THETA corresponding to the matrices A and C; */
/*             = 'N' : do not apply the bijective mapping. */
/*             The transformation performed when APPLY = 'A' allows */
/*             to get rid of the constraints norm(THETAi) < 1, i = 1:N. */
/*             A call of the SLICOT Library routine TB01VY associated to */
/*             a call of TB01VD must use the same value of APPLY. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A, assumed to be stable. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed system state matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed system input matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed system output matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,L). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading L-by-M part of this array must contain the */
/*             system input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,L). */

/*     X0      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state of the */
/*             system, x0. */
/*             On exit, this array contains the transformed initial state */
/*             of the system, corresponding to the output normal form */
/*             with parameter vector THETA. */

/*     THETA   (output) DOUBLE PRECISION array, dimension (LTHETA) */
/*             The leading N*(L+M+1)+L*M part of this array contains the */
/*             parameter vector that defines a system (A, B, C, D, x0) */
/*             which is equivalent up to a similarity transformation to */
/*             the system given on entry. The parameters are: */

/*             THETA(1:N*L)                      : parameters for A, C; */
/*             THETA(N*L+1:N*(L+M))              : parameters for B; */
/*             THETA(N*(L+M)+1:N*(L+M)+L*M)      : parameters for D; */
/*             THETA(N*(L+M)+L*M+1:N*(L+M+1)+L*M): parameters for x0. */

/*     LTHETA  INTEGER */
/*             The length of array THETA.  LTHETA >= N*(L+M+1)+L*M. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*N*L + N*L + N, */
/*                           N*N + MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), */
/*                                     N*M)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the Lyapunov equation A'*Q*A - Q = -scale^2*C'*C */
/*                   could only be solved with scale = 0; */
/*             = 2:  if matrix A is not discrete-time stable; */
/*             = 3:  if the QR algorithm failed to converge for */
/*                   matrix A. */

/*     METHOD */

/*     The matrices A and C are converted to output normal form. */
/*     First, the Lyapunov equation */

/*        A'*Q*A - Q = -scale^2*C'*C, */

/*     is solved in the Cholesky factor T, T'*T = Q, and then T is used */
/*     to get the transformation matrix. */

/*     The matrix B and the initial state x0 are transformed accordingly. */

/*     Then, the QR factorization of the transposed observability matrix */
/*     is computed, and the matrix Q is used to further transform the */
/*     system matrices. The parameters characterizing A and C are finally */
/*     obtained by applying a set of N orthogonal transformations. */

/*     REFERENCES */

/*     [1] Peeters, R.L.M., Hanzon, B., and Olivi, M. */
/*         Balanced realizations of discrete-time stable all-pass */
/*         systems and the tangential Schur algorithm. */
/*         Proceedings of the European Control Conference, */
/*         31 August - 3 September 1999, Karlsruhe, Germany. */
/*         Session CP-6, Discrete-time Systems, 1999. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Feb. 2002, Feb. 2004. */

/*     KEYWORDS */

/*     Asymptotically stable, Lyapunov equation, output normal form, */
/*     parameter estimation, similarity transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

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
    --x0;
    --theta;
    --dwork;

    /* Function Body */
    lapply = lsame_(apply, "A", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! (lapply || lsame_(apply, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldc < max(1,*l)) {
	*info = -10;
    } else if (*ldd < max(1,*l)) {
	*info = -12;
    } else if (*ltheta < *n * (*m + *l + 1) + *l * *m) {
	*info = -15;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n * (*n + max(*n,*l) + 6) + min(*n,*l), i__4 = *n * *m;
	i__1 = 1, i__2 = *n * *n * *l + *n * *l + *n, i__1 = max(i__1,i__2), 
		i__2 = *n * *n + max(i__3,i__4);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -17;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TB01VD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*n,*m);
    if (max(i__1,*l) == 0) {
	dwork[1] = 1.;
	return 0;
    } else if (*n == 0) {
	i__1 = max(1,*l);
	dlacpy_("Full", l, m, &d__[d_offset], ldd, &theta[1], &i__1, (ftnlen)
		4);
	dwork[1] = 1.;
	return 0;
    } else if (*l == 0) {
	dlacpy_("Full", n, m, &b[b_offset], ldb, &theta[1], n, (ftnlen)4);
	dcopy_(n, &x0[1], &c__1, &theta[*n * *m + 1], &c__1);
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    wrkopt = 1;
    piby2 = atan(1.) * 2.;

/*     Convert A and C to output normal form. */
/*     First, solve the Lyapunov equation */
/*        A'*Q*A - Q = -scale^2*C'*C, */
/*     in the Cholesky factor T, T'*T = Q, and use T to get the */
/*     transformation matrix. Copy A and C, to preserve them. */

/*     Workspace: need   N*(2*N + MAX(N,L) + 6) + MIN(N,L). */
/*                prefer larger. */

/*     Initialize the indices in the workspace. */

    ldt = max(*n,*l);
    ca = 1;
    ia = 1;
    it = ia + *n * *n;
    iu = it + ldt * *n;
    iwr = iu + *n * *n;
    iwi = iwr + *n;

    jwork = iwi + *n;

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[it], &ldt, (ftnlen)4);

    i__1 = *ldwork - jwork + 1;
    sb03od_("Discrete", "NotFactored", "NoTranspose", n, l, &dwork[ia], n, &
	    dwork[iu], n, &dwork[it], &ldt, &scale, &dwork[iwr], &dwork[iwi], 
	    &dwork[jwork], &i__1, info, (ftnlen)8, (ftnlen)11, (ftnlen)11);
    if (*info != 0) {
	if (*info == 6) {
	    *info = 3;
	} else {
	    *info = 2;
	}
	return 0;
    }
    wrkopt = (integer) dwork[jwork] + jwork - 1;

    if (scale == 0.) {
	*info = 1;
	return 0;
    }

/*     Compute A = T*A*T^(-1). */

    dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, n, &c_b18, &dwork[it]
	    , &ldt, &a[a_offset], lda, (ftnlen)4, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);

    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", n, n, &c_b18, &dwork[
	    it], &ldt, &a[a_offset], lda, (ftnlen)5, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);
    if (*m > 0) {

/*        Compute B = (1/scale)*T*B. */

	d__1 = 1. / scale;
	dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, m, &d__1, &dwork[
		it], &ldt, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
    }

/*     Compute x0 = (1/scale)*T*x0. */

    dtrmv_("Upper", "NoTranspose", "NonUnit", n, &dwork[it], &ldt, &x0[1], &
	    c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
    d__1 = 1. / scale;
    dscal_(n, &d__1, &x0[1], &c__1);

/*     Compute C = scale*C*T^(-1). */

    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", l, n, &scale, &dwork[
	    it], &ldt, &c__[c_offset], ldc, (ftnlen)5, (ftnlen)5, (ftnlen)11, 
	    (ftnlen)7);

/*     Now, the system has been transformed to the output normal form. */
/*     Build the transposed observability matrix in DWORK(CA) and compute */
/*     its QR factorization. */

    ma02ad_("Full", l, n, &c__[c_offset], ldc, &dwork[ca], n, (ftnlen)4);

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dgemm_("Transpose", "NoTranspose", n, l, n, &c_b18, &a[a_offset], lda,
		 &dwork[ca + (i__ - 1) * *n * *l], n, &c_b42, &dwork[ca + i__ 
		* *n * *l], n, (ftnlen)9, (ftnlen)11);
/* L10: */
    }

/*     Compute the QR factorization. */

/*     Workspace: need   N*N*L + N + L*N. */
/*                prefer N*N*L + N + NB*L*N. */

    itau = ca + *n * *n * *l;
    jwork = itau + *n;
    i__1 = *l * *n;
    i__2 = *ldwork - jwork + 1;
    dgeqrf_(n, &i__1, &dwork[ca], n, &dwork[itau], &dwork[jwork], &i__2, info)
	    ;
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Compute Q such that R has all diagonal elements nonnegative. */
/*     Only the first N*N part of R is needed. Move the details */
/*     of the QR factorization process, to gain memory and efficiency. */

/*     Workspace: need   2*N*N + 2*N. */
/*                prefer 2*N*N + N + NB*N. */

    ir = *n * *n + 1;
    if (*l != 2) {
	dcopy_(n, &dwork[itau], &c__1, &dwork[ir + *n * *n], &c__1);
    }
    dlacpy_("Lower", n, n, &dwork[ca], n, &dwork[ir], n, (ftnlen)5);
    itau = ir + *n * *n;
    jwork = itau + *n;

    iq = 1;
    dlaset_("Full", n, n, &c_b42, &c_b18, &dwork[iq], n, (ftnlen)4);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dwork[ir + (i__ - 1) * (*n + 1)] < 0.) {
	    dwork[iq + (i__ - 1) * (*n + 1)] = -1.;
	}
/* L20: */
    }

    i__1 = *ldwork - jwork + 1;
    dormqr_("Left", "NoTranspose", n, n, n, &dwork[ir], n, &dwork[itau], &
	    dwork[iq], n, &dwork[jwork], &i__1, info, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);
    jwork = ir;

/*     Now, the transformation matrix Q is in DWORK(IQ). */

/*     Compute A = Q'*A*Q. */

    dgemm_("Transpose", "NoTranspose", n, n, n, &c_b18, &dwork[iq], n, &a[
	    a_offset], lda, &c_b42, &dwork[jwork], n, (ftnlen)9, (ftnlen)11);
    dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b18, &dwork[jwork], n, &
	    dwork[iq], n, &c_b42, &a[a_offset], lda, (ftnlen)11, (ftnlen)11);

    if (*m > 0) {

/*        Compute B = Q'*B. */
/*        Workspace: need   N*N + N*M. */

	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[jwork], n, (ftnlen)4);
	dgemm_("Transpose", "NoTranspose", n, m, n, &c_b18, &dwork[iq], n, &
		dwork[jwork], n, &c_b42, &b[b_offset], ldb, (ftnlen)9, (
		ftnlen)11);
    }

/*     Compute C = C*Q. */
/*     Workspace: need   N*N + N*L. */

    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[jwork], l, (ftnlen)4);
    dgemm_("NoTranspose", "NoTranspose", l, n, n, &c_b18, &dwork[jwork], l, &
	    dwork[iq], n, &c_b42, &c__[c_offset], ldc, (ftnlen)11, (ftnlen)11)
	    ;

/*     Compute x0 = Q'*x0. */

    dcopy_(n, &x0[1], &c__1, &dwork[jwork], &c__1);
    dgemv_("Transpose", n, n, &c_b18, &dwork[iq], n, &dwork[jwork], &c__1, &
	    c_b42, &x0[1], &c__1, (ftnlen)9);

/*     Now, copy C and A into the workspace to make it easier to read out */
/*     the corresponding part of THETA, and to apply the transformations. */

    ldca = *n + *l;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(l, &c__[i__ * c_dim1 + 1], &c__1, &dwork[ca + (i__ - 1) * ldca]
		, &c__1);
	dcopy_(n, &a[i__ * a_dim1 + 1], &c__1, &dwork[ca + *l + (i__ - 1) * 
		ldca], &c__1);
/* L30: */
    }

    jwork = ca + ldca * *n;

/*     The parameters characterizing A and C are extracted in this loop. */
/*     Workspace: need   N*(N + L + 1). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(l, &dwork[ca + 1 + (*n - i__) * (ldca + 1)], &c__1, &theta[(
		i__ - 1) * *l + 1], &c__1);
	ri = dwork[ca + (*n - i__) * (ldca + 1)];
	ti = dnrm2_(l, &theta[(i__ - 1) * *l + 1], &c__1);

/*        Multiply the part of [C; A] which will be currently transformed */
/*        with Ui = [ -THETAi, Si; RI, THETAi' ] from the left, without */
/*        storing Ui. Ui has the size (L+1)-by-(L+1). */

	dgemv_("Transpose", l, n, &c_b18, &dwork[ca + *n - i__ + 1], &ldca, &
		theta[(i__ - 1) * *l + 1], &c__1, &c_b42, &dwork[jwork], &
		c__1, (ftnlen)9);

	if (ti > 0.) {
	    d__1 = (ri - 1.) / ti / ti;
	    dger_(l, n, &d__1, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[
		    jwork], &c__1, &dwork[ca + *n - i__ + 1], &ldca);
	} else {

/*           The call below is for the limiting case. */

	    dger_(l, n, &c_b93, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[
		    jwork], &c__1, &dwork[ca + *n - i__ + 1], &ldca);
	}

	dger_(l, n, &c_b96, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[ca + *n 
		- i__], &ldca, &dwork[ca + *n - i__ + 1], &ldca);
	daxpy_(n, &ri, &dwork[ca + *n - i__], &ldca, &dwork[jwork], &c__1);

/*        Move these results to their appropriate locations. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    in = ca + *n - i__ + (j - 1) * ldca;
	    i__3 = in + *l;
	    for (k = in + 1; k <= i__3; ++k) {
		dwork[k - 1] = dwork[k];
/* L40: */
	    }
	    dwork[in + *l] = dwork[jwork + j - 1];
/* L50: */
	}

/*        Now, apply the bijective mapping, which allows to get rid */
/*        of the constraint norm(THETAi) < 1. */

	if (lapply && ti != 0.) {
	    d__1 = tan(ti * piby2) / ti;
	    dscal_(l, &d__1, &theta[(i__ - 1) * *l + 1], &c__1);
	}

/* L60: */
    }

    if (*m > 0) {

/*        The next part of THETA is B. */

	dlacpy_("Full", n, m, &b[b_offset], ldb, &theta[*n * *l + 1], n, (
		ftnlen)4);

/*        Copy the matrix D. */

	dlacpy_("Full", l, m, &d__[d_offset], ldd, &theta[*n * (*l + *m) + 1],
		 l, (ftnlen)4);
    }

/*     Copy the initial state x0. */

    dcopy_(n, &x0[1], &c__1, &theta[*n * (*l + *m) + *l * *m + 1], &c__1);

    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of TB01VD *** */
} /* tb01vd_ */

