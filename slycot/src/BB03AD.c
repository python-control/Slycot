/* BB03AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b84 = -1.;
static doublereal c_b132 = 2.;

/* Subroutine */ int bb03ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *vec, integer *n, integer *m, doublereal *e, 
	integer *lde, doublereal *a, integer *lda, doublereal *y, integer *
	ldy, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *u, integer *ldu, char *note, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen def_len, ftnlen note_len)
{
    /* Initialized data */

    static logical vecdef[8] = { TRUE_,TRUE_,FALSE_,TRUE_,TRUE_,FALSE_,FALSE_,
	    FALSE_ };

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double pow_di(doublereal *, integer *), pow_dd(doublereal *, doublereal *)
	    ;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ttm1, ttp1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal twobyn;


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

/*     To generate benchmark examples of (generalized) continuous-time */
/*     Lyapunov equations */

/*        T           T */
/*       A  X  E  +  E  X A  =  Y . */

/*     In some examples, the right hand side has the form */

/*                T */
/*       Y  =  - B  B */

/*     and the solution can be represented as a product of Cholesky */
/*     factors */

/*              T */
/*       X  =  U  U . */

/*     E, A, Y, X, and U are real N-by-N matrices, and B is M-by-N. Note */
/*     that E can be the identity matrix. For some examples, B, X, or U */
/*     are not provided. */

/*     This routine is an implementation of the benchmark library */
/*     CTLEX (Version 1.0) described in [1]. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER*1 */
/*             Specifies the kind of values used as parameters when */
/*             generating parameter-dependent and scalable examples */
/*             (i.e., examples with NR(1) = 2, 3, or 4): */
/*             DEF = 'D' or 'd': Default values are used. */
/*             DEF = 'N' or 'n': Values set in DPAR and IPAR are used. */
/*             This parameter is not referenced if NR(1) = 1. */
/*             Note that the scaling parameter of examples with */
/*             NR(1) = 3 or 4 is considered as a regular parameter in */
/*             this context. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension 2 */
/*             Specifies the index of the desired example according */
/*             to [1]. */
/*             NR(1) defines the group: */
/*                   1 : parameter-free problems of fixed size */
/*                   2 : parameter-dependent problems of fixed size */
/*                   3 : parameter-free problems of scalable size */
/*                   4 : parameter-dependent problems of scalable size */
/*             NR(2) defines the number of the benchmark example */
/*             within a certain group according to [1]. */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension 2 */
/*             On entry, if DEF = 'N' or 'n' and the desired example */
/*             depends on real parameters, then the array DPAR must */
/*             contain the values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Example 4.1, DPAR(1) and DPAR(2) define 'r' and 's', */
/*             respectively. */
/*             For Example 4.2, DPAR(1) and DPAR(2) define 'lambda' and */
/*             's', respectively. */
/*             For Examples 4.3 and 4.4, DPAR(1) defines the parameter */
/*             't'. */
/*             On exit, if DEF = 'D' or 'd' and the desired example */
/*             depends on real parameters, then the array DPAR is */
/*             overwritten by the default values given in [1]. */

/*     IPAR    (input/output) INTEGER array of DIMENSION at least 1 */
/*             On entry, if DEF = 'N' or 'n' and the desired example */
/*             depends on integer parameters, then the array IPAR must */
/*             contain the values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 4.1, 4.2, and 4.3, IPAR(1) defines 'n'. */
/*             For Example 4.4, IPAR(1) defines 'q'. */
/*             On exit, if DEF = 'D' or 'd' and the desired example */
/*             depends on integer parameters, then the array IPAR is */
/*             overwritten by the default values given in [1]. */

/*     VEC     (output) LOGICAL array, dimension 8 */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(1) and VEC(2) refer to N and M, respectively, and are */
/*             always .TRUE. */
/*             VEC(3) is .TRUE. iff E is NOT the identity matrix. */
/*             VEC(4) and VEC(5) refer to A and Y, respectively, and are */
/*             always .TRUE. */
/*             VEC(6) is .TRUE. iff B is provided. */
/*             VEC(7) is .TRUE. iff the solution matrix X is provided. */
/*             VEC(8) is .TRUE. iff the Cholesky factor U is provided. */

/*     N       (output) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices E and A. */

/*     M       (output) INTEGER */
/*             The number of rows in the matrix B. If B is not provided */
/*             for the desired example, M = 0 is returned. */

/*     E       (output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix E. */
/*             NOTE that this array is overwritten (by the identity */
/*             matrix), if VEC(3) = .FALSE. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= N. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix Y. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading M-by-N part of this array contains the */
/*             matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= M. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix X. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= N. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix U. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= N. */

/*     NOTE    (output) CHARACTER*70 */
/*             String containing short information about the chosen */
/*             example. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             For Examples 4.1 and 4.2., LDWORK >= 2*IPAR(1) is */
/*             required. */
/*             For the other examples, no workspace is needed, i.e., */
/*             LDWORK >= 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; in particular, INFO = -3 or -4 indicates */
/*                   that at least one of the parameters in DPAR or */
/*                   IPAR, respectively, has an illegal value. */

/*     REFERENCES */

/*     [1]  D. Kressner, V. Mehrmann, and T. Penzl. */
/*          CTLEX - a Collection of Benchmark Examples for Continuous- */
/*          Time Lyapunov Equations. */
/*          SLICOT Working Note 1999-6, 1999. */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz) */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please contact Volker Mehrmann */
/*     (Email: volker.mehrmann@mathematik.tu-chemnitz.de). */

/*     REVISIONS */

/*     June 1999, V. Sima. */

/*     KEYWORDS */

/*     continuous-time Lyapunov equations */

/*     ******************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. Intrinsic Functions .. */
/*     .. Data Statements .. */
/*     . default values for availabilities . */
    /* Parameter adjustments */
    --nr;
    --dpar;
    --ipar;
    --vec;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --dwork;

    /* Function Body */

/*     .. Executable Statements .. */

    *info = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	vec[i__] = vecdef[i__ - 1];
/* L10: */
    }

    if (nr[1] == 4) {
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
	    *info = -1;
	    return 0;
	}

	if (nr[2] == 1) {
	    s_copy(note, "CTLEX: Example 4.1", (ftnlen)70, (ftnlen)18);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 10;
		dpar[1] = 1.5;
		dpar[2] = 1.5;
	    }
	    if (dpar[1] <= 1. || dpar[2] <= 1.) {
		*info = -3;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 1;
	    if (*lde < *n) {
		*info = -9;
	    }
	    if (*lda < *n) {
		*info = -11;
	    }
	    if (*ldy < *n) {
		*info = -13;
	    }
	    if (*ldb < *m) {
		*info = -15;
	    }
	    if (*ldx < *n) {
		*info = -17;
	    }
	    if (*ldwork < *n << 1) {
		*info = -22;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[6] = TRUE_;
	    vec[7] = TRUE_;
	    twobyn = 2. / (doublereal) (*n);
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
	    dlaset_("A", m, n, &c_b8, &c_b8, &b[b_offset], ldb, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &c_b8, &x[x_offset], ldx, (ftnlen)1);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		temp = pow_di(&dpar[1], &i__2);
		a[j + j * a_dim1] = -temp;
		dwork[j] = 1.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__ - 1;
		    x[i__ + j * x_dim1] = (doublereal) (i__ * j) / (temp + 
			    pow_di(&dpar[1], &i__3));
/* L20: */
		}
/* L30: */
	    }
/*         H1 * A */
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         H1 * X */
	    dgemv_("T", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &x[
		    x_offset], ldx);
/*         X * H1 */
	    dgemv_("N", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &x[
		    x_offset], ldx);
/*         S A INV(S), INV(S) X INV(S), B INV(S) */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		b[j * b_dim1 + 1] = (doublereal) (j - *n - 1) / pow_di(&dpar[
			2], &i__2);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__ + j - 2;
		    x[i__ + j * x_dim1] /= pow_di(&dpar[2], &i__3);
		    i__3 = i__ - j;
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
/* L40: */
		}
		dwork[j] = 1. - j % 2 * 2.;
/* L50: */
	    }
/*         H2 * A */
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         H2 * X */
	    dgemv_("T", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &x[
		    x_offset], ldx);
/*         X * H2 */
	    dgemv_("N", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &x[
		    x_offset], ldx);
/*         B * H2 */
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
	    dger_(n, n, &c_b84, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);

	} else if (nr[2] == 2) {
	    s_copy(note, "CTLEX: Example 4.2", (ftnlen)70, (ftnlen)18);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 10;
		dpar[1] = -.5;
		dpar[2] = 1.5;
	    }
	    if (dpar[1] >= 0. || dpar[2] <= 1.) {
		*info = -3;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 1;
	    if (*lde < *n) {
		*info = -9;
	    }
	    if (*lda < *n) {
		*info = -11;
	    }
	    if (*ldy < *n) {
		*info = -13;
	    }
	    if (*ldb < *m) {
		*info = -15;
	    }
	    if (*ldwork < *n << 1) {
		*info = -22;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[6] = TRUE_;
	    twobyn = 2. / (doublereal) (*n);
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &dpar[1], &a[a_offset], lda, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
	    d__1 = -twobyn;
	    d__2 = 1. - twobyn;
	    dlaset_("A", m, n, &d__1, &d__2, &b[b_offset], ldb, (ftnlen)1);
	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dwork[i__] = 1.;
		a[i__ + (i__ + 1) * a_dim1] = 1.;
/* L60: */
	    }
	    dwork[*n] = 1.;
/*         H1 * A */
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         S A INV(S), B INV(S) */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		b[j * b_dim1 + 1] /= pow_di(&dpar[2], &i__2);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__ - j;
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
/* L70: */
		}
		dwork[j] = 1. - j % 2 * 2.;
/* L80: */
	    }
/*         H2 * A */
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
	    d__1 = -twobyn;
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         B * H2 */
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
	    dger_(n, n, &c_b84, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);

	} else if (nr[2] == 3) {
	    s_copy(note, "CTLEX: Example 4.3", (ftnlen)70, (ftnlen)18);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 10;
		dpar[1] = 10.;
	    }
	    if (dpar[1] < 0.) {
		*info = -3;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 0;
	    if (*lde < *n) {
		*info = -9;
	    }
	    if (*lda < *n) {
		*info = -11;
	    }
	    if (*ldy < *n) {
		*info = -13;
	    }
	    if (*ldx < *n) {
		*info = -17;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[3] = TRUE_;
	    vec[7] = TRUE_;
	    d__1 = -dpar[1];
	    temp = pow_dd(&c_b132, &d__1);
	    dlaset_("U", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("L", n, n, &temp, &c_b9, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("L", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
	    dlaset_("U", n, n, &c_b9, &c_b8, &a[a_offset], lda, (ftnlen)1);
	    dlaset_("A", n, n, &c_b9, &c_b9, &x[x_offset], ldx, (ftnlen)1);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + i__ * a_dim1] = (doublereal) (i__ - 1) + temp;
/* L90: */
	    }
/* Computing 2nd power */
	    d__1 = temp;
	    y[y_dim1 + 1] = temp * 2. + (doublereal) (*n - 1) * 2. * (d__1 * 
		    d__1);
/* Computing 2nd power */
	    d__1 = temp;
	    ttp1 = (doublereal) (*n + 1) * 2. * temp + 2. - d__1 * d__1;
/* Computing 2nd power */
	    d__1 = temp;
	    ttm1 = (doublereal) (*n - 1) * 2. * temp + 2. - d__1 * d__1;
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		y[i__ + y_dim1] = y[y_dim1 + 1] + (doublereal) (i__ - 1) * 
			ttm1;
/* L100: */
	    }
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__ + j * y_dim1] = y[i__ + y_dim1] + (doublereal) (j - 
			    1) * (ttp1 - i__ * 4. * temp);
/* L110: */
		}
/* L120: */
	    }

	} else if (nr[2] == 4) {
	    s_copy(note, "CTLEX: Example 4.4", (ftnlen)70, (ftnlen)18);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 10;
		dpar[1] = 1.5;
	    }
	    if (dpar[1] < 1.) {
		*info = -3;
	    }
	    if (ipar[1] < 1) {
		*info = -4;
	    }
	    *n = ipar[1] * 3;
	    *m = 1;
	    if (*lde < *n) {
		*info = -9;
	    }
	    if (*lda < *n) {
		*info = -11;
	    }
	    if (*ldy < *n) {
		*info = -13;
	    }
	    if (*ldb < *m) {
		*info = -15;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[3] = TRUE_;
	    vec[6] = TRUE_;
	    dlaset_("A", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp = -pow_di(&dpar[1], &i__);
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    for (k = 0; k <= 2; ++k) {
			a[*n - i__ * 3 + 3 + (j * 3 - k) * a_dim1] = temp;
			a[*n - i__ * 3 + 2 + (j * 3 - k) * a_dim1] = temp * 
				2.;
/* L130: */
		    }
/* L140: */
		}
		a[*n - i__ * 3 + 3 + (i__ * 3 - 2) * a_dim1] = temp;
		a[*n - i__ * 3 + 2 + (i__ * 3 - 2) * a_dim1] = temp * 2.;
		a[*n - i__ * 3 + 2 + (i__ * 3 - 1) * a_dim1] = temp * 2.;
		a[*n - i__ * 3 + 2 + i__ * 3 * a_dim1] = temp;
		a[*n - i__ * 3 + 1 + i__ * 3 * a_dim1] = temp;
/* L150: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (j > 1) {
		    daxpy_(n, &c_b9, &a[j - 1 + a_dim1], lda, &a[j + a_dim1], 
			    lda);
		}
		b[j * b_dim1 + 1] = (doublereal) j;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    e[i__ + (*n - j + 1) * e_dim1] = (doublereal) min(i__,j);
		    y[i__ + j * y_dim1] = -((doublereal) (i__ * j));
/* L160: */
		}
/* L170: */
	    }

	} else {
	    *info = -2;
	}
    } else {
	*info = -2;
    }

    return 0;
/* *** Last Line of BB03AD *** */
} /* bb03ad_ */

