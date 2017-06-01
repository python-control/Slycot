/* MB02JX.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int mb02jx_(char *job, integer *k, integer *l, integer *m, 
	integer *n, doublereal *tc, integer *ldtc, doublereal *tr, integer *
	ldtr, integer *rnk, doublereal *q, integer *ldq, doublereal *r__, 
	integer *ldr, integer *jpvt, doublereal *tol1, doublereal *tol2, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jj, kk, mk, pp, pt, gap, len, pdp, pdq, nzc, pdw;
    static doublereal nrm;
    static integer pnq, pnr, ppr, rdef, rrdf, ierr;
    static logical last;
    static doublereal temp;
    static integer rrnk;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal ltol1, ltol2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02kd_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dscal_(integer *, doublereal *, doublereal *, integer *), mb02cu_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen), mb02cv_(char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer cpcol;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical compq;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dgeqp3_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), dorgqr_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    static integer wrkmin, wrkopt;


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

/*     To compute a low rank QR factorization with column pivoting of a */
/*     K*M-by-L*N block Toeplitz matrix T with blocks of size (K,L); */
/*     specifically, */
/*                                     T */
/*                           T P =  Q R , */

/*     where R is lower trapezoidal, P is a block permutation matrix */
/*     and Q^T Q = I. The number of columns in R is equivalent to the */
/*     numerical rank of T with respect to the given tolerance TOL1. */
/*     Note that the pivoting scheme is local, i.e., only columns */
/*     belonging to the same block in T are permuted. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the output of the routine as follows: */
/*             = 'Q':  computes Q and R; */
/*             = 'R':  only computes R. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in one block of T.  K >= 0. */

/*     L       (input)  INTEGER */
/*             The number of columns in one block of T.  L >= 0. */

/*     M       (input)  INTEGER */
/*             The number of blocks in one block column of T.  M >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in one block row of T.  N >= 0. */

/*     TC      (input) DOUBLE PRECISION array, dimension (LDTC, L) */
/*             The leading M*K-by-L part of this array must contain */
/*             the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,M*K). */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L) */
/*             The leading K-by-(N-1)*L part of this array must contain */
/*             the first block row of T without the leading K-by-L */
/*             block. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     RNK     (output)  INTEGER */
/*             The number of columns in R, which is equivalent to the */
/*             numerical rank of T. */

/*     Q       (output)  DOUBLE PRECISION array, dimension (LDQ,RNK) */
/*             If JOB = 'Q', then the leading M*K-by-RNK part of this */
/*             array contains the factor Q. */
/*             If JOB = 'R', then this array is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= MAX(1,M*K),  if JOB = 'Q'; */
/*             LDQ >= 1,           if JOB = 'R'. */

/*     R       (output)  DOUBLE PRECISION array, dimension (LDR,RNK) */
/*             The leading N*L-by-RNK part of this array contains the */
/*             lower trapezoidal factor R. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1,N*L) */

/*     JPVT    (output)  INTEGER array, dimension (MIN(M*K,N*L)) */
/*             This array records the column pivoting performed. */
/*             If JPVT(j) = k, then the j-th column of T*P was */
/*             the k-th column of T. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If TOL1 >= 0.0, the user supplied diagonal tolerance; */
/*             if TOL1 < 0.0, a default diagonal tolerance is used. */

/*     TOL2    DOUBLE PRECISION */
/*             If TOL2 >= 0.0, the user supplied offdiagonal tolerance; */
/*             if TOL2 < 0.0, a default offdiagonal tolerance is used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK;  DWORK(2) and DWORK(3) return the used values */
/*             for TOL1 and TOL2, respectively. */
/*             On exit, if INFO = -19,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 3, ( M*K + ( N - 1 )*L )*( L + 2*K ) + 9*L */
/*                                 + MAX(M*K,(N-1)*L) ),    if JOB = 'Q'; */
/*             LDWORK >= MAX( 3, ( N - 1 )*L*( L + 2*K + 1 ) + 9*L, */
/*                                 M*K*( L + 1 ) + L ),     if JOB = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  due to perturbations induced by roundoff errors, or */
/*                   removal of nearly linearly dependent columns of the */
/*                   generator, the Schur algorithm encountered a */
/*                   situation where a diagonal element in the negative */
/*                   generator is larger in magnitude than the */
/*                   corresponding diagonal element in the positive */
/*                   generator (modulo TOL1); */
/*             = 2:  due to perturbations induced by roundoff errors, or */
/*                   removal of nearly linearly dependent columns of the */
/*                   generator, the Schur algorithm encountered a */
/*                   situation where diagonal elements in the positive */
/*                   and negative generator are equal in magnitude */
/*                   (modulo TOL1), but the offdiagonal elements suggest */
/*                   that these columns are not linearly dependent */
/*                   (modulo TOL2*ABS(diagonal element)). */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */
/*     If, during the process, the hyperbolic norm of a row in the */
/*     leading part of the generator is found to be less than or equal */
/*     to TOL1, then this row is not reduced. If the difference of the */
/*     corresponding columns has a norm less than or equal to TOL2 times */
/*     the magnitude of the leading element, then this column is removed */
/*     from the generator, as well as from R. Otherwise, the algorithm */
/*     breaks down. TOL1 is set to norm(TC)*sqrt(eps) and TOL2 is set */
/*     to N*L*sqrt(eps) by default. */
/*     If M*K > L, the columns of T are permuted so that the diagonal */
/*     elements in one block column of R have decreasing magnitudes. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(K*RNK*L*M*N) floating point operations. */

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
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --jpvt;
    --dwork;

    /* Function Body */
    *info = 0;
    wrkopt = 3;
    mk = *m * *k;
    compq = lsame_(job, "Q", (ftnlen)1, (ftnlen)1);
    if (compq) {
/* Computing MAX */
/* Computing MAX */
	i__3 = mk, i__4 = (*n - 1) * *l;
	i__1 = 3, i__2 = (mk + (*n - 1) * *l) * (*l + (*k << 1)) + *l * 9 + 
		max(i__3,i__4);
	wrkmin = max(i__1,i__2);
    } else {
/* Computing MAX */
/* Computing MAX */
	i__3 = (*n - 1) * *l * (*l + (*k << 1) + 1) + *l * 9, i__4 = mk * (*l 
		+ 1) + *l;
	i__1 = 3, i__2 = max(i__3,i__4);
	wrkmin = max(i__1,i__2);
    }

/*     Check the scalar input parameters. */

    if (! (compq || lsame_(job, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*k < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*ldtc < max(1,mk)) {
	*info = -7;
    } else if (*ldtr < max(1,*k)) {
	*info = -9;
    } else if (*ldq < 1 || compq && *ldq < mk) {
	*info = -12;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * *l;
	if (*ldr < max(i__1,i__2)) {
	    *info = -14;
	} else if (*ldwork < wrkmin) {
	    dwork[1] = (doublereal) wrkmin;
	    *info = -19;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02JX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*m,*n), i__1 = min(i__1,*k);
    if (min(i__1,*l) == 0) {
	*rnk = 0;
	dwork[1] = (doublereal) wrkopt;
	dwork[2] = 0.;
	dwork[3] = 0.;
	return 0;
    }

    wrkopt = wrkmin;

    if (mk <= *l) {

/*        Catch M*K <= L. */

	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &dwork[1], &mk, (ftnlen)
		3);
	pdw = mk * *l + 1;
	jwork = pdw + mk;
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(&mk, l, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &i__1, &
		ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	ma02ad_("Upper part", &mk, l, &dwork[1], &mk, &r__[r_offset], ldr, (
		ftnlen)10);
	i__1 = *ldwork - jwork + 1;
	dorgqr_(&mk, &mk, &mk, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	if (compq) {
	    dlacpy_("All", &mk, &mk, &dwork[1], &mk, &q[q_offset], ldq, (
		    ftnlen)3);
	}
	pdw = mk * mk + 1;
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *ldwork - pdw + 1;
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &mk, &c_b10, &c_b11, &
		    tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1], &mk,
		     &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &i__2, &ierr, (
		    ftnlen)3, (ftnlen)9);
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);

	i__1 = mk;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jpvt[i__] = i__;
/* L10: */
	}

	*rnk = mk;
	dwork[1] = (doublereal) wrkopt;
	dwork[2] = 0.;
	dwork[3] = 0.;
	return 0;
    }

/*     Compute the generator: */

/*     1st column of the generator. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jpvt[i__] = 0;
/* L20: */
    }

    ltol1 = *tol1;
    ltol2 = *tol2;

    if (compq) {
	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &q[q_offset], ldq, (
		ftnlen)3);
	i__1 = *ldwork - *l;
	dgeqp3_(&mk, l, &q[q_offset], ldq, &jpvt[1], &dwork[1], &dwork[*l + 1]
		, &i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
	wrkopt = max(i__1,i__2);

	if (ltol1 < 0.) {

/*           Compute default tolerance LTOL1. */

/*           Estimate the 2-norm of the first block column of the */
/*           matrix with 5 power iterations. */

	    temp = 1. / sqrt((doublereal) (*l));
	    dlaset_("All", l, &c__1, &temp, &temp, &dwork[*l + 1], &c__1, (
		    ftnlen)3);

	    for (i__ = 1; i__ <= 5; ++i__) {
		dtrmv_("Upper", "NonTranspose", "NonUnit", l, &q[q_offset], 
			ldq, &dwork[*l + 1], &c__1, (ftnlen)5, (ftnlen)12, (
			ftnlen)7);
		dtrmv_("Upper", "Transpose", "NonUnit", l, &q[q_offset], ldq, 
			&dwork[*l + 1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)
			7);
		nrm = dnrm2_(l, &dwork[*l + 1], &c__1);
		d__1 = 1. / nrm;
		dscal_(l, &d__1, &dwork[*l + 1], &c__1);
/* L30: */
	    }

	    ltol1 = sqrt(nrm * dlamch_("Epsilon", (ftnlen)7));
	}

	i__ = *l;

L40:
	if ((d__1 = q[i__ + i__ * q_dim1], abs(d__1)) <= ltol1) {
	    --i__;
	    if (i__ > 0) {
		goto L40;
	    }
	}

	rrnk = i__;
	rrdf = *l - rrnk;
	ma02ad_("Upper", &rrnk, l, &q[q_offset], ldq, &r__[r_offset], ldr, (
		ftnlen)5);
	if (rrnk > 1) {
	    i__1 = *l - 1;
	    i__2 = rrnk - 1;
	    dlaset_("Upper", &i__1, &i__2, &c_b11, &c_b11, &r__[(r_dim1 << 1) 
		    + 1], ldr, (ftnlen)5);
	}
	i__1 = *ldwork - *l;
	dorgqr_(&mk, l, &rrnk, &q[q_offset], ldq, &dwork[1], &dwork[*l + 1], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
	wrkopt = max(i__1,i__2);
	if (*n > 1) {
	    i__1 = *n - 1;
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &rrnk, &c_b10, &c_b11,
		     &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &q[q_offset],
		     ldq, &r__[*l + 1 + r_dim1], ldr, &dwork[1], ldwork, &
		    ierr, (ftnlen)3, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	}

    } else {

	pdw = mk * *l + 1;
	jwork = pdw + *l;
	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &dwork[1], &mk, (ftnlen)
		3);
	i__1 = *ldwork - jwork + 1;
	dgeqp3_(&mk, l, &dwork[1], &mk, &jpvt[1], &dwork[pdw], &dwork[jwork], 
		&i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

	if (ltol1 < 0.) {

/*           Compute default tolerance LTOL1. */

/*           Estimate the 2-norm of the first block column of the */
/*           matrix with 5 power iterations. */

	    temp = 1. / sqrt((doublereal) (*l));
	    dlaset_("All", l, &c__1, &temp, &temp, &dwork[jwork], &c__1, (
		    ftnlen)3);

	    for (i__ = 1; i__ <= 5; ++i__) {
		dtrmv_("Upper", "NonTranspose", "NonUnit", l, &dwork[1], &mk, 
			&dwork[jwork], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)
			7);
		dtrmv_("Upper", "Transpose", "NonUnit", l, &dwork[1], &mk, &
			dwork[jwork], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)7);
		nrm = dnrm2_(l, &dwork[jwork], &c__1);
		d__1 = 1. / nrm;
		dscal_(l, &d__1, &dwork[jwork], &c__1);
/* L50: */
	    }

	    ltol1 = sqrt(nrm * dlamch_("Epsilon", (ftnlen)7));
	}

	rrnk = *l;
	i__ = (*l - 1) * mk + *l;

L60:
	if ((d__1 = dwork[i__], abs(d__1)) <= ltol1) {
	    --rrnk;
	    i__ = i__ - mk - 1;
	    if (i__ > 0) {
		goto L60;
	    }
	}

	rrdf = *l - rrnk;
	ma02ad_("Upper part", &rrnk, l, &dwork[1], &mk, &r__[r_offset], ldr, (
		ftnlen)10);
	if (rrnk > 1) {
	    i__1 = *l - 1;
	    i__2 = rrnk - 1;
	    dlaset_("Upper", &i__1, &i__2, &c_b11, &c_b11, &r__[(r_dim1 << 1) 
		    + 1], ldr, (ftnlen)5);
	}
	i__1 = *ldwork - jwork + 1;
	dorgqr_(&mk, l, &rrnk, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *ldwork - pdw + 1;
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &rrnk, &c_b10, &c_b11,
		     &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1], &
		    mk, &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &i__2, &ierr,
		     (ftnlen)3, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	    wrkopt = max(i__1,i__2);
	}
    }

/*     Quick return if N = 1. */

    if (*n == 1) {
	*rnk = rrnk;
	dwork[1] = (doublereal) wrkopt;
	dwork[2] = ltol1;
	dwork[3] = 0.;
	return 0;
    }

/*     Compute default tolerance LTOL2. */

    if (ltol2 < 0.) {
	ltol2 = (doublereal) (*n * *l) * sqrt(dlamch_("Epsilon", (ftnlen)7));
    }

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	dcopy_(&rrnk, &r__[j + r_dim1], ldr, &r__[*l + jpvt[j] + (rrnk + 1) * 
		r_dim1], ldr);
/* L70: */
    }

    if (*n > 2) {
	i__1 = (*n - 2) * *l;
	dlacpy_("All", &i__1, &rrnk, &r__[*l + 1 + r_dim1], ldr, &r__[(*l << 
		1) + 1 + (rrnk + 1) * r_dim1], ldr, (ftnlen)3);
    }

/*     2nd column of the generator. */

    if (rrdf > 0) {
	i__1 = min(rrdf,*k);
	i__2 = (*n - 1) * *l;
	ma02ad_("All", &i__1, &i__2, &tr[tr_offset], ldtr, &r__[*l + 1 + ((
		rrnk << 1) + 1) * r_dim1], ldr, (ftnlen)3);
    }
    if (*k > rrdf) {
	i__1 = *k - rrdf;
	i__2 = (*n - 1) * *l;
	i__3 = (*n - 1) * *l;
	ma02ad_("All", &i__1, &i__2, &tr[rrdf + 1 + tr_dim1], ldtr, &dwork[1],
		 &i__3, (ftnlen)3);
    }

/*     3rd column of the generator. */

/* Computing MAX */
    i__1 = 0, i__2 = *k - rrdf;
    pnr = (*n - 1) * *l * max(i__1,i__2) + 1;
    i__1 = (*n - 1) * *l;
    i__2 = (*n - 1) * *l;
    dlacpy_("All", &i__1, &rrnk, &r__[*l + 1 + r_dim1], ldr, &dwork[pnr], &
	    i__2, (ftnlen)3);

/*     4th column of the generator. */

    pdw = pnr + (*n - 1) * *l * rrnk;
    pt = (*m - 1) * *k + 1;

/* Computing MIN */
    i__2 = *m, i__3 = *n - 1;
    i__1 = min(i__2,i__3);
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = (*n - 1) * *l;
	ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw], &i__2, (
		ftnlen)3);
	pt -= *k;
	pdw += *l;
/* L80: */
    }

    pt = 1;

    i__1 = *n - 1;
    for (i__ = *m + 1; i__ <= i__1; ++i__) {
	i__2 = (*n - 1) * *l;
	ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw], &i__2, 
		(ftnlen)3);
	pt += *l;
	pdw += *l;
/* L90: */
    }

    if (compq) {
	pdq = pnr + (*n - 1) * *l * (rrnk + *k);
/* Computing MAX */
	i__1 = 0, i__2 = *k - rrdf;
	pnq = pdq + mk * max(i__1,i__2);
	pdw = pnq + mk * (rrnk + *k);
	dlacpy_("All", &mk, &rrnk, &q[q_offset], ldq, &dwork[pnq], &mk, (
		ftnlen)3);
	if (*m > 1) {
	    i__1 = (*m - 1) * *k;
	    dlacpy_("All", &i__1, &rrnk, &q[q_offset], ldq, &q[*k + 1 + (rrnk 
		    + 1) * q_dim1], ldq, (ftnlen)3);
	}
	dlaset_("All", k, &rrnk, &c_b11, &c_b11, &q[(rrnk + 1) * q_dim1 + 1], 
		ldq, (ftnlen)3);
	if (rrdf > 0) {
	    dlaset_("All", &mk, &rrdf, &c_b11, &c_b10, &q[((rrnk << 1) + 1) * 
		    q_dim1 + 1], ldq, (ftnlen)3);
	}
/* Computing MAX */
	i__2 = 0, i__3 = *k - rrdf;
	i__1 = max(i__2,i__3);
	dlaset_("All", &rrdf, &i__1, &c_b11, &c_b11, &dwork[pdq], &mk, (
		ftnlen)3);
	i__1 = *m * *k - rrdf;
/* Computing MAX */
	i__3 = 0, i__4 = *k - rrdf;
	i__2 = max(i__3,i__4);
	dlaset_("All", &i__1, &i__2, &c_b11, &c_b10, &dwork[pdq + rrdf], &mk, 
		(ftnlen)3);
	dlaset_("All", &mk, k, &c_b11, &c_b11, &dwork[pnq + mk * rrnk], &mk, (
		ftnlen)3);
    } else {
	pdw = pnr + (*n - 1) * *l * (rrnk + *k);
    }
    ppr = 1;
    *rnk = rrnk;
    rdef = rrdf;
    len = *n * *l;
/* Computing MIN */
    i__1 = *n * *l;
    gap = *n * *l - min(i__1,mk);

/*     KK is the number of columns in the leading part of the */
/*     generator. After sufficiently many rank drops or if */
/*     M*K < N*L it may be less than L. */

/* Computing MIN */
    i__1 = *l + *k - rdef;
    kk = min(i__1,*l);
/* Computing MIN */
    i__1 = kk, i__2 = mk - *l;
    kk = min(i__1,i__2);

/*     Generator reduction process. */

/* Computing MIN */
    i__2 = mk, i__3 = *n * *l;
    i__1 = min(i__2,i__3);
    i__4 = *l;
    for (i__ = *l + 1; i__4 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__4) {
/* Computing MIN */
	i__2 = mk, i__3 = *n * *l;
	if (i__ + *l <= min(i__2,i__3)) {
	    last = FALSE_;
	} else {
	    last = TRUE_;
	}
/* Computing MAX */
	i__2 = *k - rdef;
	pp = kk + max(i__2,0);
	len -= *l;
	i__2 = *l + *k - rdef;
	i__3 = (*n - 1) * *l;
	i__5 = (*n - 1) * *l;
	i__6 = *ldwork - pdw - *l * 5 + 1;
	mb02cu_("Deficient", &kk, &pp, &i__2, &c_n1, &r__[i__ + (*rnk + 1) * 
		r_dim1], ldr, &dwork[ppr], &i__3, &dwork[pnr], &i__5, &rrnk, &
		jpvt[i__], &dwork[pdw], &ltol1, &dwork[pdw + *l * 5], &i__6, &
		ierr, (ftnlen)9);
	if (ierr != 0) {

/*           Error return:  The current generator is indefinite. */

	    *info = 1;
	    return 0;
	}

/*        Apply pivoting to other columns of R. */

	pdp = pdw + *l * 6 - i__;

	i__2 = i__ + kk - 1;
	for (j = i__; j <= i__2; ++j) {
	    jpvt[j] = jpvt[j] + i__ - 1;
	    dwork[pdp + jpvt[j]] = (doublereal) j;
/* L100: */
	}

	i__2 = i__ + kk - 1;
	for (j = i__; j <= i__2; ++j) {
	    temp = (doublereal) j;
	    jj = j - 1;

L110:
	    ++jj;
	    if (dwork[pdp + jj] != temp) {
		goto L110;
	    }

	    if (jj != j) {
		dwork[pdp + jj] = dwork[pdp + j];
		dswap_(rnk, &r__[j + r_dim1], ldr, &r__[jj + r_dim1], ldr);
	    }
/* L120: */
	}

	i__2 = i__ + *l - 1;
	for (j = i__ + kk; j <= i__2; ++j) {
	    jpvt[j] = j;
/* L130: */
	}

/*        Apply reduction to other rows of R. */

	if (len > kk) {
	    i__2 = len - kk;
	    i__3 = *l + *k - rdef;
	    i__5 = (*n - 1) * *l;
	    i__6 = (*n - 1) * *l;
	    i__7 = (*n - 1) * *l;
	    i__8 = (*n - 1) * *l;
	    i__9 = *ldwork - pdw - *l * 5 + 1;
	    mb02cv_("Deficient", "NoStructure", &kk, &i__2, &pp, &i__3, &c_n1,
		     &rrnk, &r__[i__ + (*rnk + 1) * r_dim1], ldr, &dwork[ppr],
		     &i__5, &dwork[pnr], &i__6, &r__[i__ + kk + (*rnk + 1) * 
		    r_dim1], ldr, &dwork[ppr + kk], &i__7, &dwork[pnr + kk], &
		    i__8, &dwork[pdw], &dwork[pdw + *l * 5], &i__9, &ierr, (
		    ftnlen)9, (ftnlen)11);
	}

/*        Apply reduction to Q. */

	if (compq) {
	    i__2 = *l + *k - rdef;
	    i__3 = (*n - 1) * *l;
	    i__5 = (*n - 1) * *l;
	    i__6 = *ldwork - pdw - *l * 5 + 1;
	    mb02cv_("Deficient", "NoStructure", &kk, &mk, &pp, &i__2, &c_n1, &
		    rrnk, &r__[i__ + (*rnk + 1) * r_dim1], ldr, &dwork[ppr], &
		    i__3, &dwork[pnr], &i__5, &q[(*rnk + 1) * q_dim1 + 1], 
		    ldq, &dwork[pdq], &mk, &dwork[pnq], &mk, &dwork[pdw], &
		    dwork[pdw + *l * 5], &i__6, &ierr, (ftnlen)9, (ftnlen)11);
	}

/*        Inspection of the rank deficient columns: */
/*        Look for small diagonal entries. */

	nzc = 0;

	i__2 = rrnk + 1;
	for (j = kk; j >= i__2; --j) {
	    if ((d__1 = r__[i__ + j - 1 + (*rnk + j) * r_dim1], abs(d__1)) <= 
		    ltol1) {
		++nzc;
	    }
/* L140: */
	}

/*        The last NZC columns of the generator cannot be removed. */
/*        Now, decide whether for the other rank deficient columns */
/*        it is safe to remove. */

	pt = pnr;

	i__2 = kk - nzc;
	for (j = rrnk + 1; j <= i__2; ++j) {
	    temp = r__[i__ + j - 1 + (*rnk + j) * r_dim1];
	    i__3 = len - j - gap;
	    dscal_(&i__3, &temp, &r__[i__ + j + (*rnk + j) * r_dim1], &c__1);
	    i__3 = len - j - gap;
	    d__1 = -dwork[pt + j - 1];
	    daxpy_(&i__3, &d__1, &dwork[pt + j], &c__1, &r__[i__ + j + (*rnk 
		    + j) * r_dim1], &c__1);
	    i__3 = len - j - gap;
	    if (dnrm2_(&i__3, &r__[i__ + j + (*rnk + j) * r_dim1], &c__1) > 
		    ltol2 * abs(temp)) {

/*              Unlucky case: */
/*              It is neither advisable to remove the whole column nor */
/*              possible to remove the diagonal entries by Hyperbolic */
/*              rotations. */

		*info = 2;
		return 0;
	    }
	    pt += (*n - 1) * *l;
/* L150: */
	}

/*        Annihilate unwanted elements in the factor R. */

	rrdf = kk - rrnk;
	i__2 = i__ - 1;
	dlaset_("All", &i__2, &rrnk, &c_b11, &c_b11, &r__[(*rnk + 1) * r_dim1 
		+ 1], ldr, (ftnlen)3);
	i__2 = *l - 1;
	i__3 = rrnk - 1;
	dlaset_("Upper", &i__2, &i__3, &c_b11, &c_b11, &r__[i__ + (*rnk + 2) *
		 r_dim1], ldr, (ftnlen)5);

/*        Construct the generator for the next step. */

	if (! last) {

/*           Compute KK for the next step. */

/* Computing MIN */
	    i__2 = *l + *k - rdef - rrdf + nzc;
	    kk = min(i__2,*l);
/* Computing MIN */
	    i__2 = kk, i__3 = mk - i__ - *l + 1;
	    kk = min(i__2,i__3);

	    if (kk <= 0) {
		*rnk += rrnk;
		goto L200;
	    }

	    dlaset_("All", l, &rrdf, &c_b11, &c_b11, &r__[i__ + (*rnk + rrnk 
		    + 1) * r_dim1], ldr, (ftnlen)3);

/*           The columns with small diagonal entries form parts of the */
/*           new positive generator. */

	    if (rrdf - nzc > 0 && nzc > 0) {
		cpcol = min(nzc,kk);

		i__2 = *rnk + rrnk + cpcol;
		for (j = *rnk + rrnk + 1; j <= i__2; ++j) {
		    i__3 = len - *l;
		    dcopy_(&i__3, &r__[i__ + *l + (j + rrdf - nzc) * r_dim1], 
			    &c__1, &r__[i__ + *l + j * r_dim1], &c__1);
/* L160: */
		}

	    }

/*           Construct the leading parts of the positive generator. */

/* Computing MIN */
	    i__2 = rrnk, i__3 = kk - nzc;
	    cpcol = min(i__2,i__3);
	    if (cpcol > 0) {

		i__2 = i__ + *l - 1;
		for (j = i__; j <= i__2; ++j) {
		    dcopy_(&cpcol, &r__[j + (*rnk + 1) * r_dim1], ldr, &r__[
			    jpvt[j] + *l + (*rnk + rrnk + nzc + 1) * r_dim1], 
			    ldr);
/* L170: */
		}

		if (len > *l << 1) {
		    i__2 = len - (*l << 1);
		    dlacpy_("All", &i__2, &cpcol, &r__[i__ + *l + (*rnk + 1) *
			     r_dim1], ldr, &r__[i__ + (*l << 1) + (*rnk + 
			    rrnk + nzc + 1) * r_dim1], ldr, (ftnlen)3);
		}
	    }
	    ppr += *l;

/*           Refill the leading parts of the positive generator. */

/* Computing MIN */
	    i__2 = *k - rdef, i__3 = kk - rrnk - nzc;
	    cpcol = min(i__2,i__3);
	    if (cpcol > 0) {
		i__2 = len - *l;
		i__3 = (*n - 1) * *l;
		dlacpy_("All", &i__2, &cpcol, &dwork[ppr], &i__3, &r__[i__ + *
			l + (*rnk + (rrnk << 1) + nzc + 1) * r_dim1], ldr, (
			ftnlen)3);
		ppr += cpcol * (*n - 1) * *l;
	    }
	    pnr = pnr + (rrdf - nzc) * (*n - 1) * *l + *l;

/*           Do the same things for Q. */

	    if (compq) {
		if (rrdf - nzc > 0 && nzc > 0) {
		    cpcol = min(nzc,kk);

		    i__2 = *rnk + rrnk + cpcol;
		    for (j = *rnk + rrnk + 1; j <= i__2; ++j) {
			dcopy_(&mk, &q[(j + rrdf - nzc) * q_dim1 + 1], &c__1, 
				&q[j * q_dim1 + 1], &c__1);
/* L180: */
		    }

		}
/* Computing MIN */
		i__2 = rrnk, i__3 = kk - nzc;
		cpcol = min(i__2,i__3);
		if (cpcol > 0) {
		    dlaset_("All", k, &cpcol, &c_b11, &c_b11, &q[(*rnk + rrnk 
			    + nzc + 1) * q_dim1 + 1], ldq, (ftnlen)3);
		    if (*m > 1) {
			i__2 = (*m - 1) * *k;
			dlacpy_("All", &i__2, &cpcol, &q[(*rnk + 1) * q_dim1 
				+ 1], ldq, &q[*k + 1 + (*rnk + rrnk + nzc + 1)
				 * q_dim1], ldq, (ftnlen)3);
		    }
		}
/* Computing MIN */
		i__2 = *k - rdef, i__3 = kk - rrnk - nzc;
		cpcol = min(i__2,i__3);
		if (cpcol > 0) {
		    dlacpy_("All", &mk, &cpcol, &dwork[pdq], &mk, &q[(*rnk + (
			    rrnk << 1) + nzc + 1) * q_dim1 + 1], ldq, (ftnlen)
			    3);
		    pdq += cpcol * mk;
		}
		pnq += (rrdf - nzc) * mk;
	    }
	}
	*rnk += rrnk;
	rdef = rdef + rrdf - nzc;
/* L190: */
    }

L200:
    dwork[1] = (doublereal) wrkopt;
    dwork[2] = ltol1;
    dwork[3] = ltol2;

/* *** Last line of MB02JX *** */
    return 0;
} /* mb02jx_ */

