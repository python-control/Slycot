/* MB02ID.f -- translated by f2c (version 20100827).
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
static doublereal c_b18 = 1.;
static doublereal c_b42 = -1.;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02id_(char *job, integer *k, integer *l, integer *m, 
	integer *n, integer *rb, integer *rc, doublereal *tc, integer *ldtc, 
	doublereal *tr, integer *ldtr, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, 
	    i__10, i__11, i__12;

    /* Local variables */
    static integer i__, x, y, nb, kk, pt, pdi, len, pni, ppi, pdw, rnk, pnr, 
	    ppr, ierr, ipvt[1];
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02kd_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    mb02cu_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dgemm_(char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb02cv_(char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dgels_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    static logical compo, compu;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), dgeqrf_(integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin;
    extern /* Subroutine */ int dtrtri_(char *, char *, integer *, doublereal 
	    *, integer *, integer *, ftnlen, ftnlen);
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

/*     To solve the overdetermined or underdetermined real linear systems */
/*     involving an M*K-by-N*L block Toeplitz matrix T that is specified */
/*     by its first block column and row. It is assumed that T has full */
/*     rank. */
/*     The following options are provided: */

/*     1. If JOB = 'O' or JOB = 'A' :  find the least squares solution of */
/*        an overdetermined system, i.e., solve the least squares problem */

/*                  minimize || B - T*X ||.                           (1) */

/*     2. If JOB = 'U' or JOB = 'A' :  find the minimum norm solution of */
/*        the undetermined system */
/*                   T */
/*                  T * X = C.                                        (2) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the problem to be solved as follows */
/*             = 'O':  solve the overdetermined system (1); */
/*             = 'U':  solve the underdetermined system (2); */
/*             = 'A':  solve (1) and (2). */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T. */
/*             0 <= N <= M*K / L. */

/*     RB      (input) INTEGER */
/*             If JOB = 'O' or 'A', the number of columns in B.  RB >= 0. */

/*     RC      (input) INTEGER */
/*             If JOB = 'U' or 'A', the number of columns in C.  RC >= 0. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry, the leading M*K-by-L part of this array must */
/*             contain the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC.  LDTC >= MAX(1,M*K) */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L) */
/*             On entry, the leading K-by-(N-1)*L part of this array must */
/*             contain the 2nd to the N-th blocks of the first block row */
/*             of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     B       (input/output)  DOUBLE PRECISION array, dimension (LDB,RB) */
/*             On entry, if JOB = 'O' or JOB = 'A', the leading M*K-by-RB */
/*             part of this array must contain the right hand side */
/*             matrix B of the overdetermined system (1). */
/*             On exit, if JOB = 'O' or JOB = 'A', the leading N*L-by-RB */
/*             part of this array contains the solution of the */
/*             overdetermined system (1). */
/*             This array is not referenced if JOB = 'U'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,M*K),  if JOB = 'O'  or  JOB = 'A'; */
/*             LDB >= 1,           if JOB = 'U'. */

/*     C       (input)  DOUBLE PRECISION array, dimension (LDC,RC) */
/*             On entry, if JOB = 'U' or JOB = 'A', the leading N*L-by-RC */
/*             part of this array must contain the right hand side */
/*             matrix C of the underdetermined system (2). */
/*             On exit, if JOB = 'U' or JOB = 'A', the leading M*K-by-RC */
/*             part of this array contains the solution of the */
/*             underdetermined system (2). */
/*             This array is not referenced if JOB = 'O'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDB >= 1,           if JOB = 'O'; */
/*             LDB >= MAX(1,M*K),  if JOB = 'U'  or  JOB = 'A'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Let x = MAX( 2*N*L*(L+K) + (6+N)*L,(N*L+M*K+1)*L + M*K ) */
/*             and y = N*M*K*L + N*L, then */
/*             if MIN( M,N ) = 1 and JOB = 'O', */
/*                         LDWORK >= MAX( y + MAX( M*K,RB ),1 ); */
/*             if MIN( M,N ) = 1 and JOB = 'U', */
/*                         LDWORK >= MAX( y + MAX( M*K,RC ),1 ); */
/*             if MIN( M,N ) = 1 and JOB = 'A', */
/*                         LDWORK >= MAX( y +MAX( M*K,MAX( RB,RC ),1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'O', */
/*                         LDWORK >= MAX( x,N*L*RB + 1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'U', */
/*                         LDWORK >= MAX( x,N*L*RC + 1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'A', */
/*                         LDWORK >= MAX( x,N*L*MAX( RB,RC ) + 1 ). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The Toeplitz matrix */
/*                   associated with T is (numerically) not of full rank. */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O( L*L*K*(N+M)*log(N+M) + N*N*L*L*(L+K) ) */
/*     and additionally */

/*     if JOB = 'O' or JOB = 'A', */
/*                  O( (K*L+RB*L+K*RB)*(N+M)*log(N+M) + N*N*L*L*RB ); */
/*     if JOB = 'U' or JOB = 'A', */
/*                  O( (K*L+RC*L+K*RC)*(N+M)*log(N+M) + N*N*L*L*RC ); */

/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */
/*     D. Kressner, Technical Univ. Berlin, Germany, July 2002. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    tc_dim1 = *ldtc;
    tc_offset = 1 + tc_dim1;
    tc -= tc_offset;
    tr_dim1 = *ldtr;
    tr_offset = 1 + tr_dim1;
    tr -= tr_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    compo = lsame_(job, "O", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
    compu = lsame_(job, "U", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
/* Computing MAX */
    i__1 = (*n << 1) * *l * (*l + *k) + (*n + 6) * *l, i__2 = (*n * *l + *m * 
	    *k + 1) * *l + *m * *k;
    x = max(i__1,i__2);
    y = *n * *m * *k * *l + *n * *l;
    if (min(*m,*n) == 1) {
/* Computing MAX */
	i__1 = *m * *k;
	wrkmin = max(i__1,1);
	if (compo) {
	    wrkmin = max(wrkmin,*rb);
	}
	if (compu) {
	    wrkmin = max(wrkmin,*rc);
	}
/* Computing MAX */
	i__1 = y + wrkmin;
	wrkmin = max(i__1,1);
    } else {
	wrkmin = x;
	if (compo) {
/* Computing MAX */
	    i__1 = wrkmin, i__2 = *n * *l * *rb + 1;
	    wrkmin = max(i__1,i__2);
	}
	if (compu) {
/* Computing MAX */
	    i__1 = wrkmin, i__2 = *n * *l * *rc + 1;
	    wrkmin = max(i__1,i__2);
	}
    }
    wrkopt = 1;

/*     Check the scalar input parameters. */

    if (! (compo || compu)) {
	*info = -1;
    } else if (*k < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0 || *n * *l > *m * *k) {
	*info = -5;
    } else if (compo && *rb < 0) {
	*info = -6;
    } else if (compu && *rc < 0) {
	*info = -7;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m * *k;
	if (*ldtc < max(i__1,i__2)) {
	    *info = -9;
	} else if (*ldtr < max(1,*k)) {
	    *info = -11;
	} else if (*ldb < 1 || compo && *ldb < *m * *k) {
	    *info = -13;
	} else if (*ldc < 1 || compu && *ldc < *m * *k) {
	    *info = -15;
	} else if (*ldwork < wrkmin) {
	    dwork[1] = (doublereal) wrkmin;
	    *info = -17;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02ID", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = *n * *l;
    if (compo && min(i__1,*rb) == 0) {
	compo = FALSE_;
    }
/* Computing MIN */
    i__1 = *n * *l;
    if (compu && min(i__1,*rc) == 0) {
	i__1 = *m * *k;
	dlaset_("Full", &i__1, rc, &c_b8, &c_b8, &c__[c_offset], ldc, (ftnlen)
		4);
	compu = FALSE_;
    }
    if (! (compo || compu)) {
	dwork[1] = 1.;
	return 0;
    }

/*     Check cases M = 1 or N = 1. */

    if (min(*m,*n) == 1) {
	pdw = *k * *l * *m * *n;
	if (compo) {
	    i__1 = *m * *k;
	    i__2 = *m * *k;
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		    ftnlen)3);
	    i__1 = (*n - 1) * *l;
	    i__2 = *m * *k;
	    dlacpy_("All", k, &i__1, &tr[tr_offset], ldtr, &dwork[*k * *l + 1]
		    , &i__2, (ftnlen)3);
	    i__1 = *m * *k;
	    i__2 = *n * *l;
	    i__3 = *m * *k;
	    i__4 = *ldwork - pdw;
	    dgels_("NonTranspose", &i__1, &i__2, rb, &dwork[1], &i__3, &b[
		    b_offset], ldb, &dwork[pdw + 1], &i__4, &ierr, (ftnlen)12)
		    ;
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + 1] + pdw;
	    wrkopt = max(i__1,i__2);
	}
	if (compu) {
	    i__1 = *m * *k;
	    i__2 = *m * *k;
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		    ftnlen)3);
	    i__1 = (*n - 1) * *l;
	    i__2 = *m * *k;
	    dlacpy_("All", k, &i__1, &tr[tr_offset], ldtr, &dwork[*k * *l + 1]
		    , &i__2, (ftnlen)3);
	    i__1 = *m * *k;
	    i__2 = *n * *l;
	    i__3 = *m * *k;
	    i__4 = *ldwork - pdw;
	    dgels_("Transpose", &i__1, &i__2, rc, &dwork[1], &i__3, &c__[
		    c_offset], ldc, &dwork[pdw + 1], &i__4, &ierr, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + 1] + pdw;
	    wrkopt = max(i__1,i__2);
	}
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

/*     Step 1:  Compute the generator. */

    if (compo) {
	i__1 = *n * *l;
	i__2 = *ldwork - *n * *l * *rb;
	mb02kd_("Column", "Transpose", k, l, m, n, rb, &c_b18, &c_b8, &tc[
		tc_offset], ldtc, &tr[tr_offset], ldtr, &b[b_offset], ldb, &
		dwork[1], &i__1, &dwork[*n * *l * *rb + 1], &i__2, &ierr, (
		ftnlen)6, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[*n * *l * *rb + 1] + *n * *l * *
		rb;
	wrkopt = max(i__1,i__2);
	i__1 = *n * *l;
	i__2 = *n * *l;
	dlacpy_("All", &i__1, rb, &dwork[1], &i__2, &b[b_offset], ldb, (
		ftnlen)3);
    }

    pdw = *n * *l * *l + 1;
    i__1 = *m * *k;
    i__2 = *m * *k;
    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[pdw], &i__2, (
	    ftnlen)3);
    i__1 = *m * *k;
    i__2 = *m * *k;
    i__3 = *ldwork - pdw - (*m * *k + 1) * *l - 1;
    dgeqrf_(&i__1, l, &dwork[pdw], &i__2, &dwork[pdw + *m * *k * *l], &dwork[
	    pdw + (*m * *k + 1) * *l], &i__3, &ierr);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*m * *k + 1) * *l] + pdw + (*
	    m * *k + 1) * *l - 1;
    wrkopt = max(i__1,i__2);

    i__1 = pdw + *m * *k * *l - 1;
    i__2 = *m * *k + 1;
    for (i__ = pdw; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (dwork[i__] == 0.) {
	    *info = 1;
	    return 0;
	}
/* L10: */
    }

    i__2 = *m * *k;
    i__1 = *n * *l;
    ma02ad_("Upper", l, l, &dwork[pdw], &i__2, &dwork[1], &i__1, (ftnlen)5);
    i__2 = *m * *k;
    i__1 = *m * *k;
    i__3 = *ldwork - pdw - (*m * *k + 1) * *l - 1;
    dorgqr_(&i__2, l, l, &dwork[pdw], &i__1, &dwork[pdw + *m * *k * *l], &
	    dwork[pdw + (*m * *k + 1) * *l], &i__3, &ierr);
/* Computing MAX */
    i__2 = wrkopt, i__1 = (integer) dwork[pdw + (*m * *k + 1) * *l] + pdw + (*
	    m * *k + 1) * *l - 1;
    wrkopt = max(i__2,i__1);
    i__2 = *n - 1;
    i__1 = *m * *k;
    i__3 = *n * *l;
    i__4 = *ldwork - pdw - *m * *k * *l + 1;
    mb02kd_("Row", "Transpose", k, l, m, &i__2, l, &c_b18, &c_b8, &tc[
	    tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[pdw], &i__1, &
	    dwork[*l + 1], &i__3, &dwork[pdw + *m * *k * *l], &i__4, &ierr, (
	    ftnlen)3, (ftnlen)9);
/* Computing MAX */
    i__2 = wrkopt, i__1 = (integer) dwork[pdw + *m * *k * *l] + pdw + *m * *k 
	    * *l - 1;
    wrkopt = max(i__2,i__1);
    ppr = *n * *l * *l + 1;
    pnr = *n * *l * (*l + *k) + 1;
    i__2 = (*n - 1) * *l;
    i__1 = *n * *l;
    ma02ad_("All", k, &i__2, &tr[tr_offset], ldtr, &dwork[ppr + *l], &i__1, (
	    ftnlen)3);
    i__2 = (*n - 1) * *l;
    i__1 = *n * *l;
    i__3 = *n * *l;
    dlacpy_("All", &i__2, l, &dwork[*l + 1], &i__1, &dwork[pnr + *l], &i__3, (
	    ftnlen)3);
    pt = (*m - 1) * *k + 1;
    pdw = pnr + *n * *l * *l + *l;

/* Computing MIN */
    i__1 = *m, i__3 = *n - 1;
    i__2 = min(i__1,i__3);
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *n * *l;
	ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw], &i__1, (
		ftnlen)3);
	pt -= *k;
	pdw += *l;
/* L30: */
    }

    pt = 1;

    i__2 = *n - 1;
    for (i__ = *m + 1; i__ <= i__2; ++i__) {
	i__1 = *n * *l;
	ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw], &i__1, 
		(ftnlen)3);
	pt += *l;
	pdw += *l;
/* L40: */
    }

    if (compo) {

/*        Apply the first reduction step to T'*B. */

	i__2 = *n * *l;
	dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rb, &c_b18, &
		dwork[1], &i__2, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
	i__2 = (*n - 1) * *l;
	i__1 = *n * *l;
	dgemm_("NoTranspose", "NoTranspose", &i__2, rb, l, &c_b18, &dwork[*l 
		+ 1], &i__1, &b[b_offset], ldb, &c_b42, &b[*l + 1 + b_dim1], 
		ldb, (ftnlen)11, (ftnlen)11);
	i__2 = *n * *l;
	dtrsm_("Left", "Lower", "Transpose", "NonUnit", l, rb, &c_b18, &dwork[
		1], &i__2, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9,
		 (ftnlen)7);
    }

    if (compu) {

/*        Apply the first reduction step to C. */

	i__2 = *n * *l;
	dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rc, &c_b18, &
		dwork[1], &i__2, &c__[c_offset], ldc, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
	i__2 = (*n - 1) * *l;
	i__1 = *n * *l;
	dgemm_("NoTranspose", "NoTranspose", &i__2, rc, l, &c_b18, &dwork[*l 
		+ 1], &i__1, &c__[c_offset], ldc, &c_b42, &c__[*l + 1 + 
		c_dim1], ldc, (ftnlen)11, (ftnlen)11);
	i__2 = *n * *l;
	dtrsm_("Left", "Lower", "Transpose", "NonUnit", l, rc, &c_b18, &dwork[
		1], &i__2, &c__[c_offset], ldc, (ftnlen)4, (ftnlen)5, (ftnlen)
		9, (ftnlen)7);
    }

    pdi = (*n - 1) * *l + 1;
    i__2 = *n * *l;
    i__1 = *n * *l;
    dlacpy_("Lower", l, l, &dwork[1], &i__2, &dwork[pdi], &i__1, (ftnlen)5);
    i__2 = *n * *l;
    dtrtri_("Lower", "NonUnit", l, &dwork[pdi], &i__2, &ierr, (ftnlen)5, (
	    ftnlen)7);
    i__2 = *l - 1;
    i__1 = *n * *l;
    i__3 = *n * *l;
    ma02ad_("Lower", &i__2, l, &dwork[pdi + 1], &i__1, &dwork[((*n << 1) - 1) 
	    * *l + 1], &i__3, (ftnlen)5);
    i__2 = *l - 1;
    i__1 = *n * *l;
    dlaset_("Lower", &i__2, l, &c_b8, &c_b8, &dwork[pdi + 1], &i__1, (ftnlen)
	    5);
    i__2 = *n * *l;
    i__1 = *n * *l;
    dlacpy_("Upper", l, l, &dwork[pdi], &i__2, &dwork[pnr], &i__1, (ftnlen)5);
    i__2 = *l - 1;
    i__1 = *n * *l;
    dlaset_("Lower", &i__2, l, &c_b8, &c_b8, &dwork[pnr + 1], &i__1, (ftnlen)
	    5);
    i__2 = *n * *l;
    dlaset_("All", l, k, &c_b8, &c_b8, &dwork[ppr], &i__2, (ftnlen)3);
    i__2 = *n * *l;
    dlaset_("All", l, k, &c_b8, &c_b8, &dwork[pnr + *n * *l * *l], &i__2, (
	    ftnlen)3);

    ppi = ppr;
    ppr += *l;
    pni = pnr;
    pnr += *l;
    pdw = (*n << 1) * *l * (*l + *k) + 1;
    len = (*n - 1) * *l;

/*     Determine block size for the involved block Householder */
/*     transformations. */

/* Computing MIN */
    i__1 = *n * *l;
    i__2 = ilaenv_(&c__1, "DGELQF", " ", &i__1, l, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
    nb = min(i__2,*l);
    kk = pdw + *l * 6 - 1;
/* Computing MAX */
    i__2 = wrkopt, i__1 = kk + *n * *l * nb;
    wrkopt = max(i__2,i__1);
    kk = *ldwork - kk;
    if (kk < *n * *l * nb) {
	nb = kk / (*n * *l);
    }
/* Computing MAX */
    i__3 = *n * *l;
    i__2 = 2, i__1 = ilaenv_(&c__2, "DGELQF", " ", &i__3, l, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
    nbmin = max(i__2,i__1);
    if (nb < nbmin) {
	nb = 0;
    }

    i__2 = *n * *l;
    i__1 = *l;
    for (i__ = *l + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	i__3 = *l + *k;
	i__4 = *l + *k;
	i__5 = *n * *l;
	i__6 = *n * *l;
	i__7 = *n * *l;
	i__8 = *ldwork - pdw - *l * 6 + 1;
	mb02cu_("Column", l, &i__3, &i__4, &nb, &dwork[1], &i__5, &dwork[ppr],
		 &i__6, &dwork[pnr], &i__7, &rnk, ipvt, &dwork[pdw], &c_b8, &
		dwork[pdw + *l * 6], &i__8, &ierr, (ftnlen)6);
	if (ierr != 0) {

/*           Error return:  The rank condition is (numerically) not */
/*                          satisfied. */

	    *info = 1;
	    return 0;
	}
	i__3 = len - *l;
	i__4 = *l + *k;
	i__5 = *l + *k;
	i__6 = *n * *l;
	i__7 = *n * *l;
	i__8 = *n * *l;
	i__9 = *n * *l;
	i__10 = *n * *l;
	i__11 = *n * *l;
	i__12 = *ldwork - pdw - *l * 6 + 1;
	mb02cv_("Column", "NoStructure", l, &i__3, &i__4, &i__5, &nb, &c_n1, &
		dwork[1], &i__6, &dwork[ppr], &i__7, &dwork[pnr], &i__8, &
		dwork[*l + 1], &i__9, &dwork[ppr + *l], &i__10, &dwork[pnr + *
		l], &i__11, &dwork[pdw], &dwork[pdw + *l * 6], &i__12, &ierr, 
		(ftnlen)6, (ftnlen)11);
	pdi -= *l;
	if (compo) {

/*           Block Gaussian elimination to B. */

	    i__3 = *n * *l;
	    dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rb, &c_b42, 
		    &dwork[1], &i__3, &b[i__ + b_dim1], ldb, (ftnlen)4, (
		    ftnlen)5, (ftnlen)12, (ftnlen)7);
	    if (len > *l) {
		i__3 = len - *l;
		i__4 = *n * *l;
		dgemm_("NonTranspose", "NonTranspose", &i__3, rb, l, &c_b18, &
			dwork[*l + 1], &i__4, &b[i__ + b_dim1], ldb, &c_b18, &
			b[i__ + *l + b_dim1], ldb, (ftnlen)12, (ftnlen)12);
	    }
	}
	if (compu) {

/*           Block Gaussian elimination to C. */

	    i__3 = *n * *l;
	    dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rc, &c_b42, 
		    &dwork[1], &i__3, &c__[i__ + c_dim1], ldc, (ftnlen)4, (
		    ftnlen)5, (ftnlen)12, (ftnlen)7);
	    if (len > *l) {
		i__3 = len - *l;
		i__4 = *n * *l;
		dgemm_("NonTranspose", "NonTranspose", &i__3, rc, l, &c_b18, &
			dwork[*l + 1], &i__4, &c__[i__ + c_dim1], ldc, &c_b18,
			 &c__[i__ + *l + c_dim1], ldc, (ftnlen)12, (ftnlen)12)
			;
	    }
	}
	i__3 = *n * *l;
	dlaset_("All", l, l, &c_b8, &c_b8, &dwork[pdi], &i__3, (ftnlen)3);
	i__3 = i__ + *l - 1;
	i__4 = *l + *k;
	i__5 = *l + *k;
	i__6 = *n * *l;
	i__7 = *n * *l;
	i__8 = *n * *l;
	i__9 = *n * *l;
	i__10 = *n * *l;
	i__11 = *n * *l;
	i__12 = *ldwork - pdw - *l * 6 + 1;
	mb02cv_("Column", "Triangular", l, &i__3, &i__4, &i__5, &nb, &c_n1, &
		dwork[1], &i__6, &dwork[ppr], &i__7, &dwork[pnr], &i__8, &
		dwork[pdi], &i__9, &dwork[ppi], &i__10, &dwork[pni], &i__11, &
		dwork[pdw], &dwork[pdw + *l * 6], &i__12, &ierr, (ftnlen)6, (
		ftnlen)10);
	if (compo) {

/*           Apply block Gaussian elimination to B. */

	    i__3 = i__ - 1;
	    i__4 = *n * *l;
	    dgemm_("NoTranspose", "NoTranspose", &i__3, rb, l, &c_b18, &dwork[
		    pdi], &i__4, &b[i__ + b_dim1], ldb, &c_b18, &b[b_offset], 
		    ldb, (ftnlen)11, (ftnlen)11);
	    i__3 = *n * *l;
	    dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", l, rb, &c_b18, 
		    &dwork[(*n - 1) * *l + 1], &i__3, &b[i__ + b_dim1], ldb, (
		    ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);
	}
	if (compu) {

/*           Apply block Gaussian elimination to C. */

	    i__3 = i__ - 1;
	    i__4 = *n * *l;
	    dgemm_("NonTranspose", "NonTranspose", &i__3, rc, l, &c_b18, &
		    dwork[pdi], &i__4, &c__[i__ + c_dim1], ldc, &c_b18, &c__[
		    c_offset], ldc, (ftnlen)12, (ftnlen)12);
	    i__3 = *n * *l;
	    dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", l, rc, &c_b18, 
		    &dwork[(*n - 1) * *l + 1], &i__3, &c__[i__ + c_dim1], ldc,
		     (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);
	}
	len -= *l;
	pnr += *l;
	ppr += *l;
/* L50: */
    }

    if (compu) {
	i__1 = *m * *k;
	i__2 = *ldwork - *m * *k * *rc;
	mb02kd_("Column", "NonTranspose", k, l, m, n, rc, &c_b18, &c_b8, &tc[
		tc_offset], ldtc, &tr[tr_offset], ldtr, &c__[c_offset], ldc, &
		dwork[1], &i__1, &dwork[*m * *k * *rc + 1], &i__2, &ierr, (
		ftnlen)6, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[*m * *k * *rc + 1] + *m * *k * *
		rc;
	wrkopt = max(i__1,i__2);
	i__1 = *m * *k;
	i__2 = *m * *k;
	dlacpy_("All", &i__1, rc, &dwork[1], &i__2, &c__[c_offset], ldc, (
		ftnlen)3);
    }
    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of MB02ID *** */
} /* mb02id_ */

