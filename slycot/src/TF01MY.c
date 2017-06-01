/* TF01MY.f -- translated by f2c (version 20100827).
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
static doublereal c_b8 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int tf01my_(integer *n, integer *m, integer *p, integer *ny, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *u, 
	integer *ldu, doublereal *x, doublereal *y, integer *ldy, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer nb, ik, is, ns;
    static doublereal upd;
    static integer iyl, irem, maxn;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the output sequence of a linear time-invariant */
/*     open-loop system given by its discrete-time state-space model */
/*     (A,B,C,D), where A is an N-by-N general matrix. */

/*     The initial state vector x(1) must be supplied by the user. */

/*     This routine differs from SLICOT Library routine TF01MD in the */
/*     way the input and output trajectories are stored. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NY      (input) INTEGER */
/*             The number of output vectors y(k) to be computed. */
/*             NY >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input matrix B of the system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             output matrix C of the system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             direct link matrix D of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             The leading NY-by-M part of this array must contain the */
/*             input vector sequence u(k), for k = 1,2,...,NY. */
/*             Specifically, the k-th row of U must contain u(k)'. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NY). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state vector */
/*             x(1) which consists of the N initial states of the system. */
/*             On exit, this array contains the final state vector */
/*             x(NY+1) of the N states of the system at instant NY+1. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,P) */
/*             The leading NY-by-P part of this array contains the output */
/*             vector sequence y(1),y(2),...,y(NY) such that the k-th */
/*             row of Y contains y(k)' (the outputs at instant k), */
/*             for k = 1,2,...,NY. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,NY). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= N. */
/*             For better performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an initial state vector x(1), the output vector sequence */
/*     y(1), y(2),..., y(NY) is obtained via the formulae */

/*        x(k+1) = A x(k) + B u(k) */
/*        y(k)   = C x(k) + D u(k), */

/*     where each element y(k) is a vector of length P containing the */
/*     outputs at instant k and k = 1,2,...,NY. */

/*     REFERENCES */

/*     [1] Luenberger, D.G. */
/*         Introduction to Dynamic Systems: Theory, Models and */
/*         Applications. */
/*         John Wiley & Sons, New York, 1979. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately (N + M) x (N + P) x NY */
/*     multiplications and additions. */

/*     FURTHER COMMENTS */

/*     The implementation exploits data locality and uses BLAS 3 */
/*     operations as much as possible, given the workspace length. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, multivariable system, state-space model, */
/*     state-space representation, time response. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --x;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    maxn = max(1,*n);
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*ny < 0) {
	*info = -4;
    } else if (*lda < maxn) {
	*info = -6;
    } else if (*ldb < maxn) {
	*info = -8;
    } else if (*ldc < max(1,*p)) {
	*info = -10;
    } else if (*ldd < max(1,*p)) {
	*info = -12;
    } else if (*ldu < max(1,*ny)) {
	*info = -14;
    } else if (*ldy < max(1,*ny)) {
	*info = -17;
    } else if (*ldwork < *n) {
	*info = -19;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TF01MY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*ny,*p) == 0) {
	return 0;
    } else if (*n == 0) {

/*        Non-dynamic system: compute the output vectors. */

	if (*m == 0) {
	    dlaset_("Full", ny, p, &c_b4, &c_b4, &y[y_offset], ldy, (ftnlen)4)
		    ;
	} else {
	    dgemm_("No transpose", "Transpose", ny, p, m, &c_b8, &u[u_offset],
		     ldu, &d__[d_offset], ldd, &c_b4, &y[y_offset], ldy, (
		    ftnlen)12, (ftnlen)9);
	}
	return 0;
    }

/*     Determine the block size (taken as for LAPACK routine DGETRF). */

    i__1 = max(*m,*p);
    nb = ilaenv_(&c__1, "DGETRF", " ", ny, &i__1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     Find the number of state vectors that can be accommodated in */
/*     the provided workspace and initialize. */

/* Computing MIN */
    i__1 = *ldwork / *n, i__2 = nb * nb / *n, i__1 = min(i__1,i__2);
    ns = min(i__1,*ny);

    if (ns <= 1 || *ny * max(*m,*p) <= nb * nb) {

/*        LDWORK < 2*N or small problem: */
/*                     only BLAS 2 calculations are used in the loop */
/*                     for computing the output corresponding to D = 0. */
/*        One row of the array Y is computed for each loop index value. */

	i__1 = *ny;
	for (ik = 1; ik <= i__1; ++ik) {
	    dgemv_("No transpose", p, n, &c_b8, &c__[c_offset], ldc, &x[1], &
		    c__1, &c_b4, &y[ik + y_dim1], ldy, (ftnlen)12);

	    dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &x[1], &
		    c__1, &c_b4, &dwork[1], &c__1, (ftnlen)12);
	    dgemv_("No transpose", n, m, &c_b8, &b[b_offset], ldb, &u[ik + 
		    u_dim1], ldu, &c_b8, &dwork[1], &c__1, (ftnlen)12);

	    dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);
/* L10: */
	}

    } else {

/*        LDWORK >= 2*N and large problem: */
/*        some BLAS 3 calculations can also be used. */

	iyl = *ny / ns * ns;
	if (*m == 0) {
	    upd = 0.;
	} else {
	    upd = 1.;
	}

	dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

	i__1 = iyl;
	i__2 = ns;
	for (ik = 1; i__2 < 0 ? ik >= i__1 : ik <= i__1; ik += i__2) {

/*           Compute the current NS-1 state vectors in the workspace. */

	    i__3 = ns - 1;
	    dgemm_("No transpose", "Transpose", n, &i__3, m, &c_b8, &b[
		    b_offset], ldb, &u[ik + u_dim1], ldu, &c_b4, &dwork[*n + 
		    1], &maxn, (ftnlen)12, (ftnlen)9);

	    i__3 = ns - 1;
	    for (is = 1; is <= i__3; ++is) {
		dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &dwork[
			(is - 1) * *n + 1], &c__1, &upd, &dwork[is * *n + 1], 
			&c__1, (ftnlen)12);
/* L20: */
	    }

/*           Initialize the current NS output vectors. */

	    dgemm_("Transpose", "Transpose", &ns, p, n, &c_b8, &dwork[1], &
		    maxn, &c__[c_offset], ldc, &c_b4, &y[ik + y_dim1], ldy, (
		    ftnlen)9, (ftnlen)9);

/*           Prepare the next iteration. */

	    dgemv_("No transpose", n, m, &c_b8, &b[b_offset], ldb, &u[ik + ns 
		    - 1 + u_dim1], ldu, &c_b4, &dwork[1], &c__1, (ftnlen)12);
	    dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &dwork[(ns 
		    - 1) * *n + 1], &c__1, &upd, &dwork[1], &c__1, (ftnlen)12)
		    ;
/* L30: */
	}

	irem = *ny - iyl;

	if (irem > 1) {

/*           Compute the last IREM output vectors. */
/*           First, compute the current IREM-1 state vectors. */

	    ik = iyl + 1;
	    i__2 = irem - 1;
	    dgemm_("No transpose", "Transpose", n, &i__2, m, &c_b8, &b[
		    b_offset], ldb, &u[ik + u_dim1], ldu, &c_b4, &dwork[*n + 
		    1], &maxn, (ftnlen)12, (ftnlen)9);

	    i__2 = irem - 1;
	    for (is = 1; is <= i__2; ++is) {
		dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &dwork[
			(is - 1) * *n + 1], &c__1, &upd, &dwork[is * *n + 1], 
			&c__1, (ftnlen)12);
/* L40: */
	    }

/*           Initialize the last IREM output vectors. */

	    dgemm_("Transpose", "Transpose", &irem, p, n, &c_b8, &dwork[1], &
		    maxn, &c__[c_offset], ldc, &c_b4, &y[ik + y_dim1], ldy, (
		    ftnlen)9, (ftnlen)9);

/*           Prepare the final state vector. */

	    dgemv_("No transpose", n, m, &c_b8, &b[b_offset], ldb, &u[ik + 
		    irem - 1 + u_dim1], ldu, &c_b4, &dwork[1], &c__1, (ftnlen)
		    12);
	    dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &dwork[(
		    irem - 1) * *n + 1], &c__1, &upd, &dwork[1], &c__1, (
		    ftnlen)12);

	} else if (irem == 1) {

/*           Compute the last 1 output vectors. */

	    dgemv_("No transpose", p, n, &c_b8, &c__[c_offset], ldc, &dwork[1]
		    , &c__1, &c_b4, &y[ik + y_dim1], ldy, (ftnlen)12);

/*           Prepare the final state vector. */

	    dcopy_(n, &dwork[1], &c__1, &dwork[*n + 1], &c__1);
	    dgemv_("No transpose", n, m, &c_b8, &b[b_offset], ldb, &u[ik + 
		    u_dim1], ldu, &c_b4, &dwork[1], &c__1, (ftnlen)12);
	    dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &dwork[*n 
		    + 1], &c__1, &upd, &dwork[1], &c__1, (ftnlen)12);
	}

/*        Set the final state vector. */

	dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);

    }

/*     Add the direct contribution of the input to the output vectors. */

    dgemm_("No transpose", "Transpose", ny, p, m, &c_b8, &u[u_offset], ldu, &
	    d__[d_offset], ldd, &c_b8, &y[y_offset], ldy, (ftnlen)12, (ftnlen)
	    9);

    return 0;
/* *** Last line of TF01MY *** */
} /* tf01my_ */

