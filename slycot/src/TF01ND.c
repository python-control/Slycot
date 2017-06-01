/* TF01ND.f -- translated by f2c (version 20100827).
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

static doublereal c_b6 = 0.;
static doublereal c_b10 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tf01nd_(char *uplo, integer *n, integer *m, integer *p, 
	integer *ny, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	 doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *u, integer *ldu, doublereal *x, doublereal *y, integer *
	ldy, doublereal *dwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ik;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical luplo;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);


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
/*     (A,B,C,D), where A is an N-by-N upper or lower Hessenberg matrix. */

/*     The initial state vector x(1) must be supplied by the user. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Indicates whether the user wishes to use an upper or lower */
/*             Hessenberg matrix as follows: */
/*             = 'U':  Upper Hessenberg matrix; */
/*             = 'L':  Lower Hessenberg matrix. */

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
/*             If UPLO = 'U', the leading N-by-N upper Hessenberg part */
/*             of this array must contain the state matrix A of the */
/*             system. */
/*             If UPLO = 'L', the leading N-by-N lower Hessenberg part */
/*             of this array must contain the state matrix A of the */
/*             system. */
/*             The remainder of the leading N-by-N part is not */
/*             referenced. */

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

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,NY) */
/*             The leading M-by-NY part of this array must contain the */
/*             input vector sequence u(k), for k = 1,2,...,NY. */
/*             Specifically, the k-th column of U must contain u(k). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,M). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state vector */
/*             x(1) which consists of the N initial states of the system. */
/*             On exit, this array contains the final state vector */
/*             x(NY+1) of the N states of the system at instant NY. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,NY) */
/*             The leading P-by-NY part of this array contains the output */
/*             vector sequence y(1),y(2),...,y(NY) such that the k-th */
/*             column of Y contains y(k) (the outputs at instant k), */
/*             for k = 1,2,...,NY. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,P). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

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

/*     The algorithm requires approximately ((N+M)xP + (N/2+M)xN) x NY */
/*     multiplications and additions. */

/*     FURTHER COMMENTS */

/*     The processing time required by this routine will be approximately */
/*     half that required by the SLICOT Library routine TF01MD, which */
/*     treats A as a general matrix. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TF01BD by S. Van Huffel, Katholieke */
/*     Univ. Leuven, Belgium. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003. */

/*     KEYWORDS */

/*     Discrete-time system, Hessenberg form, multivariable system, */
/*     state-space model, state-space representation, time response. */

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
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*ny < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*p)) {
	*info = -11;
    } else if (*ldd < max(1,*p)) {
	*info = -13;
    } else if (*ldu < max(1,*m)) {
	*info = -15;
    } else if (*ldy < max(1,*p)) {
	*info = -18;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TF01ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*p,*ny) == 0) {
	return 0;
    } else if (*n == 0) {

/*        Non-dynamic system: compute the output vectors. */

	if (*m == 0) {
	    dlaset_("Full", p, ny, &c_b6, &c_b6, &y[y_offset], ldy, (ftnlen)4)
		    ;
	} else {
	    dgemm_("No transpose", "No transpose", p, ny, m, &c_b10, &d__[
		    d_offset], ldd, &u[u_offset], ldu, &c_b6, &y[y_offset], 
		    ldy, (ftnlen)12, (ftnlen)12);
	}
	return 0;
    }

    dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

    i__1 = *ny;
    for (ik = 1; ik <= i__1; ++ik) {
	dgemv_("No transpose", p, n, &c_b10, &c__[c_offset], ldc, &dwork[1], &
		c__1, &c_b6, &y[ik * y_dim1 + 1], &c__1, (ftnlen)12);

	dtrmv_(uplo, "No transpose", "Non-unit", n, &a[a_offset], lda, &dwork[
		1], &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);

	if (luplo) {

	    i__2 = *n;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		dwork[i__] += a[i__ + (i__ - 1) * a_dim1] * x[i__ - 1];
/* L10: */
	    }

	} else {

	    i__2 = *n - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[i__] += a[i__ + (i__ + 1) * a_dim1] * x[i__ + 1];
/* L20: */
	    }

	}

	dgemv_("No transpose", n, m, &c_b10, &b[b_offset], ldb, &u[ik * 
		u_dim1 + 1], &c__1, &c_b10, &dwork[1], &c__1, (ftnlen)12);

	dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);
/* L30: */
    }

    dgemm_("No transpose", "No transpose", p, ny, m, &c_b10, &d__[d_offset], 
	    ldd, &u[u_offset], ldu, &c_b10, &y[y_offset], ldy, (ftnlen)12, (
	    ftnlen)12);

    return 0;
/* *** Last line of TF01ND *** */
} /* tf01nd_ */

