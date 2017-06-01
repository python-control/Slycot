/* MB02KD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 0.;
static integer c__1 = 1;
static doublereal c_b23 = 1.;

/* Subroutine */ int mb02kd_(char *ldblk, char *trans, integer *k, integer *l,
	 integer *m, integer *n, integer *r__, doublereal *alpha, doublereal *
	beta, doublereal *tc, integer *ldtc, doublereal *tr, integer *ldtr, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen ldblk_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, p, p1, p2, q1, q2, r1, r2, s1, s2;
    static doublereal t1, t2, cf;
    static integer pb, pc, jj, kk, ll, mk, ln, ir, nl;
    static doublereal sf, th;
    static integer pp, pt, icp, icq, len, pdw, dimb, dimc;
    static doublereal coef, scal;
    static integer meth, ierr, shft;
    static char wght[1];
    static integer wpos;
    extern /* Subroutine */ int dg01od_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static doublereal param;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical fullc;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical ltran;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lmult;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
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

/*     To compute the matrix product */

/*               C = alpha*op( T )*B + beta*C, */

/*     where alpha and beta are scalars and T is a block Toeplitz matrix */
/*     specified by its first block column TC and first block row TR; */
/*     B and C are general matrices of appropriate dimensions. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LDBLK   CHARACTER*1 */
/*             Specifies where the (1,1)-block of T is stored, as */
/*             follows: */
/*             = 'C':  in the first block of TC; */
/*             = 'R':  in the first block of TR. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( T ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( T ) = T; */
/*             = 'T':  op( T ) = T'; */
/*             = 'C':  op( T ) = T'. */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T.  N >= 0. */

/*     R       (input) INTEGER */
/*             The number of columns in B and C.  R >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then TC, TR and B */
/*             are not referenced. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then C need not be set */
/*             before entry. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry with LDBLK = 'C', the leading M*K-by-L part of */
/*             this array must contain the first block column of T. */
/*             On entry with LDBLK = 'R', the leading (M-1)*K-by-L part */
/*             of this array must contain the 2nd to the M-th blocks of */
/*             the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,M*K),      if LDBLK = 'C'; */
/*             LDTC >= MAX(1,(M-1)*K),  if LDBLK = 'R'. */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,k) */
/*             where k is (N-1)*L when LDBLK = 'C' and is N*L when */
/*             LDBLK = 'R'. */
/*             On entry with LDBLK = 'C', the leading K-by-(N-1)*L part */
/*             of this array must contain the 2nd to the N-th blocks of */
/*             the first block row of T. */
/*             On entry with LDBLK = 'R', the leading K-by-N*L part of */
/*             this array must contain the first block row of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     B       (input)  DOUBLE PRECISION array, dimension (LDB,R) */
/*             On entry with TRANS = 'N', the leading N*L-by-R part of */
/*             this array must contain the matrix B. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             M*K-by-R part of this array must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N*L),  if TRANS = 'N'; */
/*             LDB >= MAX(1,M*K),  if TRANS = 'T' or TRANS = 'C'. */

/*     C       (input/output)  DOUBLE PRECISION array, dimension (LDC,R) */
/*             On entry with TRANS = 'N', the leading M*K-by-R part of */
/*             this array must contain the matrix C. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             N*L-by-R part of this array must contain the matrix C. */
/*             On exit with TRANS = 'N', the leading M*K-by-R part of */
/*             this array contains the updated matrix C. */
/*             On exit with TRANS = 'T' or TRANS = 'C', the leading */
/*             N*L-by-R part of this array contains the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= MAX(1,M*K),  if TRANS = 'N'; */
/*             LDC >= MAX(1,N*L),  if TRANS = 'T' or TRANS = 'C'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 1. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     For point Toeplitz matrices or sufficiently large block Toeplitz */
/*     matrices, this algorithm uses convolution algorithms based on */
/*     the fast Hartley transforms [1]. Otherwise, TC is copied in */
/*     reversed order into the workspace such that C can be computed from */
/*     barely M matrix-by-matrix multiplications. */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O( (K*L+R*L+K*R)*(N+M)*log(N+M) + K*L*R ) */
/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     March 2004. */

/*     KEYWORDS */

/*     Convolution, elementary matrix operations, */
/*     fast Hartley transform, Toeplitz matrix. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    fullc = lsame_(ldblk, "C", (ftnlen)1, (ftnlen)1);
    ltran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);
    lmult = *alpha != 0.;
    mk = *m * *k;
    nl = *n * *l;

/*     Check the scalar input parameters. */

    if (! (fullc || lsame_(ldblk, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (ltran || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*k < 0) {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*r__ < 0) {
	*info = -7;
    } else if (lmult && fullc && *ldtc < max(1,mk)) {
	*info = -11;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = (*m - 1) * *k;
	if (lmult && ! fullc && *ldtc < max(i__1,i__2)) {
	    *info = -11;
	} else if (lmult && *ldtr < max(1,*k)) {
	    *info = -13;
	} else if (lmult && ! ltran && *ldb < max(1,nl)) {
	    *info = -15;
	} else if (lmult && ltran && *ldb < max(1,mk)) {
	    *info = -15;
	} else if (! ltran && *ldc < max(1,mk)) {
	    *info = -17;
	} else if (ltran && *ldc < max(1,nl)) {
	    *info = -17;
	} else if (*ldwork < 1) {
	    dwork[1] = 1.;
	    *info = -19;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02KD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Scale C beforehand. */

    if (*beta == 0.) {
	if (ltran) {
	    dlaset_("All", &nl, r__, &c_b9, &c_b9, &c__[c_offset], ldc, (
		    ftnlen)3);
	} else {
	    dlaset_("All", &mk, r__, &c_b9, &c_b9, &c__[c_offset], ldc, (
		    ftnlen)3);
	}
    } else if (*beta != 1.) {
	if (ltran) {

	    i__1 = *r__;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dscal_(&nl, beta, &c__[i__ * c_dim1 + 1], &c__1);
/* L10: */
	    }

	} else {

	    i__1 = *r__;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dscal_(&mk, beta, &c__[i__ * c_dim1 + 1], &c__1);
/* L20: */
	    }

	}
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(mk,nl);
    if (! lmult || min(i__1,*r__) == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     The parameter PARAM is the watershed between conventional */
/*     multiplication and convolution. This is of course depending */
/*     on the used computer architecture. The lower this value is set */
/*     the more likely the routine will use convolution to compute */
/*     op( T )*B. Note that if there is enough workspace available, */
/*     convolution is always used for point Toeplitz matrices. */

    param = 950.;

/*     Decide which method to choose, based on the block sizes and */
/*     the available workspace. */

    len = 1;
    p = 0;

L30:
    if (len < *m + *n - 1) {
	len <<= 1;
	++p;
	goto L30;
    }

    coef = (doublereal) (*m * *n) * 3. * (doublereal) (*k * *l) * (doublereal)
	     (*r__) / (doublereal) (len * (*k * *l + *l * *r__ + *k * *r__));

    if (fullc) {
	p1 = mk * *l;
	shft = 0;
    } else {
	p1 = (*m - 1) * *k * *l;
	shft = 1;
    }
    if (*k * *l == 1 && min(*m,*n) > 1) {
	wrkopt = len * (*r__ + 2) - p;
	meth = 3;
    } else if (len < *m * *n && coef >= param) {
	wrkopt = len * (*k * *l + *k * *r__ + *l * *r__ + 1) - p;
	meth = 3;
    } else {
	meth = 2;
	wrkopt = p1;
    }

    if (*ldwork < wrkopt) {
	--meth;
    }
    if (*ldwork < p1) {
	meth = 1;
    }

/*     Start computations. */

    if (meth == 1 && ! ltran) {

/*        Method 1 is the most unlucky way to multiply Toeplitz matrices */
/*        with vectors. Due to the memory restrictions it is not */
/*        possible to flip TC. */

	pc = 1;

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pt = (i__ - 1 - shft) * *k + 1;
	    pb = 1;

	    i__2 = i__;
	    for (j = shft + 1; j <= i__2; ++j) {
		dgemm_("No Transpose", "No Transpose", k, r__, l, alpha, &tc[
			pt + tc_dim1], ldtc, &b[pb + b_dim1], ldb, &c_b23, &
			c__[pc + c_dim1], ldc, (ftnlen)12, (ftnlen)12);
		pt -= *k;
		pb += *l;
/* L40: */
	    }

	    if (*n > i__ - shft) {
		i__2 = (*n - i__ + shft) * *l;
		dgemm_("No Transpose", "No Transpose", k, r__, &i__2, alpha, &
			tr[tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &
			c__[pc + c_dim1], ldc, (ftnlen)12, (ftnlen)12);
	    }
	    pc += *k;
/* L50: */
	}

    } else if (meth == 1 && ltran) {

	pb = 1;

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pt = (i__ - 1 - shft) * *k + 1;
	    pc = 1;

	    i__2 = i__;
	    for (j = shft + 1; j <= i__2; ++j) {
		dgemm_("Transpose", "No Transpose", l, r__, k, alpha, &tc[pt 
			+ tc_dim1], ldtc, &b[pb + b_dim1], ldb, &c_b23, &c__[
			pc + c_dim1], ldc, (ftnlen)9, (ftnlen)12);
		pt -= *k;
		pc += *l;
/* L60: */
	    }

	    if (*n > i__ - shft) {
		i__2 = (*n - i__ + shft) * *l;
		dgemm_("Transpose", "No Transpose", &i__2, r__, k, alpha, &tr[
			tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &c__[
			pc + c_dim1], ldc, (ftnlen)9, (ftnlen)12);
	    }
	    pb += *k;
/* L70: */
	}

    } else if (meth == 2 && ! ltran) {

/*        In method 2 TC is flipped resulting in less calls to the BLAS */
/*        routine DGEMM. Actually this seems often to be the best way to */
/*        multiply with Toeplitz matrices except the point Toeplitz */
/*        case. */

	pt = (*m - 1 - shft) * *k + 1;

	i__1 = (*m - shft) * *k * *l;
	i__2 = *k * *l;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
	    pt -= *k;
/* L80: */
	}

	pt = (*m - 1) * *k * *l + 1;
	pc = 1;

	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
	    i__3 = i__ - shft;
	    i__1 = min(i__3,*n) * *l;
	    dgemm_("No Transpose", "No Transpose", k, r__, &i__1, alpha, &
		    dwork[pt], k, &b[b_offset], ldb, &c_b23, &c__[pc + c_dim1]
		    , ldc, (ftnlen)12, (ftnlen)12);
	    if (*n > i__ - shft) {
		i__1 = (*n - i__ + shft) * *l;
		dgemm_("No Transpose", "No Transpose", k, r__, &i__1, alpha, &
			tr[tr_offset], ldtr, &b[(i__ - shft) * *l + 1 + 
			b_dim1], ldb, &c_b23, &c__[pc + c_dim1], ldc, (ftnlen)
			12, (ftnlen)12);
	    }
	    pc += *k;
	    pt -= *k * *l;
/* L90: */
	}

    } else if (meth == 2 && ltran) {

	pt = (*m - 1 - shft) * *k + 1;

	i__2 = (*m - shft) * *k * *l;
	i__1 = *k * *l;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
	    pt -= *k;
/* L100: */
	}

	pt = (*m - 1) * *k * *l + 1;
	pb = 1;

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	    i__3 = i__ - shft;
	    i__2 = min(i__3,*n) * *l;
	    dgemm_("Tranpose", "No Transpose", &i__2, r__, k, alpha, &dwork[
		    pt], k, &b[pb + b_dim1], ldb, &c_b23, &c__[c_offset], ldc,
		     (ftnlen)8, (ftnlen)12);
	    if (*n > i__ - shft) {
		i__2 = (*n - i__ + shft) * *l;
		dgemm_("Transpose", "No Transpose", &i__2, r__, k, alpha, &tr[
			tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &c__[(
			i__ - shft) * *l + 1 + c_dim1], ldc, (ftnlen)9, (
			ftnlen)12);
	    }
	    pb += *k;
	    pt -= *k * *l;
/* L110: */
	}

    } else if (meth == 3) {

/*        In method 3 the matrix-vector product is computed by a suitable */
/*        block convolution via fast Hartley transforms similar to the */
/*        SLICOT routine DE01PD. */

/*        Step 1: Copy input data into the workspace arrays. */

	pdw = 1;
	if (ltran) {
	    dimb = *k;
	    dimc = *l;
	} else {
	    dimb = *l;
	    dimc = *k;
	}
	pb = len * *k * *l;
	pc = len * (*k * *l + dimb * *r__);
	if (ltran) {
	    if (fullc) {
		i__1 = len * *k;
		dlacpy_("All", k, l, &tc[tc_offset], ldtc, &dwork[1], &i__1, (
			ftnlen)3);
	    }

	    i__1 = *n - 1 + shft;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = len * *k;
		dlacpy_("All", k, l, &tr[((i__ - 1) * *l + 1) * tr_dim1 + 1], 
			ldtr, &dwork[(i__ - shft) * *k + 1], &i__2, (ftnlen)3)
			;
/* L120: */
	    }

	    pdw = *n * *k + 1;
	    r1 = (len - *m - *n + 1) * *k;
	    i__1 = len * *k;
	    dlaset_("All", &r1, l, &c_b9, &c_b9, &dwork[pdw], &i__1, (ftnlen)
		    3);
	    pdw += r1;

	    i__1 = *k - shft * *k + 1;
	    i__2 = -(*k);
	    for (i__ = (*m - 1 - shft) * *k + 1; i__2 < 0 ? i__ >= i__1 : i__ 
		    <= i__1; i__ += i__2) {
		i__3 = len * *k;
		dlacpy_("All", k, l, &tc[i__ + tc_dim1], ldtc, &dwork[pdw], &
			i__3, (ftnlen)3);
		pdw += *k;
/* L130: */
	    }

	    pdw = pb + 1;
	    i__2 = len * *k;
	    dlacpy_("All", &mk, r__, &b[b_offset], ldb, &dwork[pdw], &i__2, (
		    ftnlen)3);
	    pdw += mk;
	    i__2 = (len - *m) * *k;
	    i__1 = len * *k;
	    dlaset_("All", &i__2, r__, &c_b9, &c_b9, &dwork[pdw], &i__1, (
		    ftnlen)3);
	} else {
	    if (! fullc) {
		i__2 = len * *k;
		dlacpy_("All", k, l, &tr[tr_offset], ldtr, &dwork[1], &i__2, (
			ftnlen)3);
	    }
	    i__2 = (*m - shft) * *k;
	    i__1 = len * *k;
	    dlacpy_("All", &i__2, l, &tc[tc_offset], ldtc, &dwork[shft * *k + 
		    1], &i__1, (ftnlen)3);
	    pdw = mk + 1;
	    r1 = (len - *m - *n + 1) * *k;
	    i__2 = len * *k;
	    dlaset_("All", &r1, l, &c_b9, &c_b9, &dwork[pdw], &i__2, (ftnlen)
		    3);
	    pdw += r1;

	    i__2 = shft * *l + 1;
	    i__1 = -(*l);
	    for (i__ = (*n - 2 + shft) * *l + 1; i__1 < 0 ? i__ >= i__2 : i__ 
		    <= i__2; i__ += i__1) {
		i__3 = len * *k;
		dlacpy_("All", k, l, &tr[i__ * tr_dim1 + 1], ldtr, &dwork[pdw]
			, &i__3, (ftnlen)3);
		pdw += *k;
/* L140: */
	    }

	    pdw = pb + 1;
	    i__1 = len * *l;
	    dlacpy_("All", &nl, r__, &b[b_offset], ldb, &dwork[pdw], &i__1, (
		    ftnlen)3);
	    pdw += nl;
	    i__1 = (len - *n) * *l;
	    i__2 = len * *l;
	    dlaset_("All", &i__1, r__, &c_b9, &c_b9, &dwork[pdw], &i__2, (
		    ftnlen)3);
	}

/*        Take point Toeplitz matrices into extra consideration. */

	if (*k * *l == 1) {
	    *(unsigned char *)wght = 'N';
	    dg01od_("OutputScrambled", wght, &len, &dwork[1], &dwork[pc + 1], 
		    &ierr, (ftnlen)15, (ftnlen)1);

	    i__1 = pb + len * *r__ - 1;
	    i__2 = len;
	    for (i__ = pb; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) 
		    {
		dg01od_("OutputScrambled", wght, &len, &dwork[i__ + 1], &
			dwork[pc + 1], &ierr, (ftnlen)15, (ftnlen)1);
		scal = *alpha / (doublereal) len;
		dwork[i__ + 1] = scal * dwork[i__ + 1] * dwork[1];
		dwork[i__ + 2] = scal * dwork[i__ + 2] * dwork[2];
		scal /= 2.;

		ln = 1;

		i__3 = p - 1;
		for (ll = 1; ll <= i__3; ++ll) {
		    ln <<= 1;
		    r1 = ln << 1;

		    i__4 = ln + ln / 2;
		    for (p1 = ln + 1; p1 <= i__4; ++p1) {
			t1 = dwork[p1] + dwork[r1];
			t2 = dwork[p1] - dwork[r1];
			th = t2 * dwork[i__ + p1];
			dwork[i__ + p1] = scal * (t1 * dwork[i__ + p1] + t2 * 
				dwork[i__ + r1]);
			dwork[i__ + r1] = scal * (t1 * dwork[i__ + r1] - th);
			--r1;
/* L150: */
		    }

/* L160: */
		}

		dg01od_("InputScrambled", wght, &len, &dwork[i__ + 1], &dwork[
			pc + 1], &ierr, (ftnlen)14, (ftnlen)1);
/* L170: */
	    }

	    pc = pb;
	    goto L420;
	}

/*        Step 2: Compute the weights for the Hartley transforms. */

	pdw = pc;
	r1 = 1;
	ln = 1;
	th = atan(1.) * 4. / (doublereal) len;

	i__2 = p - 2;
	for (ll = 1; ll <= i__2; ++ll) {
	    ln <<= 1;
	    th *= 2.;
	    cf = cos(th);
	    sf = sin(th);
	    dwork[pdw + r1] = cf;
	    dwork[pdw + r1 + 1] = sf;
	    r1 += 2;

	    i__1 = ln - 2;
	    for (i__ = 1; i__ <= i__1; i__ += 2) {
		dwork[pdw + r1] = cf * dwork[pdw + i__] - sf * dwork[pdw + 
			i__ + 1];
		dwork[pdw + r1 + 1] = sf * dwork[pdw + i__] + cf * dwork[pdw 
			+ i__ + 1];
		r1 += 2;
/* L180: */
	    }

/* L190: */
	}

	p1 = 3;
	q1 = r1 - 2;

	for (ll = p - 2; ll >= 1; --ll) {

	    i__2 = q1;
	    for (i__ = p1; i__ <= i__2; i__ += 4) {
		dwork[pdw + r1] = dwork[pdw + i__];
		dwork[pdw + r1 + 1] = dwork[pdw + i__ + 1];
		r1 += 2;
/* L200: */
	    }

	    p1 = q1 + 4;
	    q1 = r1 - 2;
/* L210: */
	}

/*        Step 3: Compute the Hartley transforms with scrambled output. */

	j = 0;
	kk = *k;

/*        WHILE   J < (L*LEN*K + R*LEN*DIMB), */

L220:

	ln = len;
	wpos = pdw + 1;

	for (pp = p - 1; pp >= 1; --pp) {
	    ln /= 2;
	    p2 = 1;
	    q2 = ln * kk + 1;
	    r2 = ln / 2 * kk + 1;
	    s2 = r2 + q2 - 1;

	    i__2 = len / (ln << 1) - 1;
	    for (i__ = 0; i__ <= i__2; ++i__) {

		i__1 = kk - 1;
		for (ir = 0; ir <= i__1; ++ir) {
		    t1 = dwork[q2 + ir + j];
		    dwork[q2 + ir + j] = dwork[p2 + ir + j] - t1;
		    dwork[p2 + ir + j] += t1;
		    t1 = dwork[s2 + ir + j];
		    dwork[s2 + ir + j] = dwork[r2 + ir + j] - t1;
		    dwork[r2 + ir + j] += t1;
/* L230: */
		}

		p1 = p2 + kk;
		q1 = p1 + ln * kk;
		r1 = q1 - (kk << 1);
		s1 = r1 + ln * kk;

		i__1 = wpos + ln - 3;
		for (jj = wpos; jj <= i__1; jj += 2) {
		    cf = dwork[jj];
		    sf = dwork[jj + 1];

		    i__3 = kk - 1;
		    for (ir = 0; ir <= i__3; ++ir) {
			t1 = dwork[p1 + ir + j] - dwork[q1 + ir + j];
			t2 = dwork[r1 + ir + j] - dwork[s1 + ir + j];
			dwork[p1 + ir + j] += dwork[q1 + ir + j];
			dwork[r1 + ir + j] += dwork[s1 + ir + j];
			dwork[q1 + ir + j] = cf * t1 + sf * t2;
			dwork[s1 + ir + j] = -cf * t2 + sf * t1;
/* L240: */
		    }

		    p1 += kk;
		    q1 += kk;
		    r1 -= kk;
		    s1 -= kk;
/* L250: */
		}

		p2 += (kk << 1) * ln;
		q2 += (kk << 1) * ln;
		r2 += (kk << 1) * ln;
		s2 += (kk << 1) * ln;
/* L260: */
	    }

	    wpos = wpos + ln - 2;
/* L270: */
	}

	i__2 = len * kk;
	i__1 = kk << 1;
	for (icp = kk + 1; i__1 < 0 ? icp >= i__2 : icp <= i__2; icp += i__1) 
		{
	    icq = icp - kk;

	    i__3 = kk - 1;
	    for (ir = 0; ir <= i__3; ++ir) {
		t1 = dwork[icp + ir + j];
		dwork[icp + ir + j] = dwork[icq + ir + j] - t1;
		dwork[icq + ir + j] += t1;
/* L280: */
	    }

/* L290: */
	}

	j += len * kk;
	if (j == *l * len * *k) {
	    kk = dimb;
	}
	if (j < pc) {
	    goto L220;
	}
/*        END WHILE 220 */

/*        Step 4: Compute a Hadamard like product. */

	i__1 = len - p;
	dcopy_(&i__1, &dwork[pdw + 1], &c__1, &dwork[pdw + 1 + *r__ * len * 
		dimc], &c__1);
	pdw += *r__ * len * dimc;
	scal = *alpha / (doublereal) len;
	p1 = 1;
	r1 = len * *k * *l + 1;
	s1 = r1 + len * dimb * *r__;
	if (ltran) {
	    kk = *l;
	    ll = *k;
	} else {
	    kk = *k;
	    ll = *l;
	}
	i__1 = len * *k;
	i__2 = len * dimb;
	i__3 = len * dimc;
	dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1], &i__1,
		 &dwork[r1], &i__2, &c_b9, &dwork[s1], &i__3, (ftnlen)1, (
		ftnlen)12);
	p1 += *k;
	r1 += dimb;
	s1 += dimc;
	i__1 = len * *k;
	i__2 = len * dimb;
	i__3 = len * dimc;
	dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1], &i__1,
		 &dwork[r1], &i__2, &c_b9, &dwork[s1], &i__3, (ftnlen)1, (
		ftnlen)12);
	scal /= 2.;
	ln = 1;

	i__1 = p - 1;
	for (pp = 1; pp <= i__1; ++pp) {
	    ln <<= 1;
	    p2 = ((ln << 1) - 1) * *k + 1;
	    r1 = pb + ln * dimb + 1;
	    r2 = pb + ((ln << 1) - 1) * dimb + 1;
	    s1 = pc + ln * dimc + 1;
	    s2 = pc + ((ln << 1) - 1) * dimc + 1;

	    i__2 = (ln + ln / 2) * *k;
	    i__3 = *k;
	    for (p1 = ln * *k + 1; i__3 < 0 ? p1 >= i__2 : p1 <= i__2; p1 += 
		    i__3) {

		i__4 = len * *k * (*l - 1);
		i__5 = len * *k;
		for (j = 0; i__5 < 0 ? j >= i__4 : j <= i__4; j += i__5) {

		    i__6 = p1 + *k - 1;
		    for (i__ = p1; i__ <= i__6; ++i__) {
			t1 = dwork[p2];
			dwork[p2] = dwork[j + i__] - t1;
			dwork[j + i__] += t1;
			++p2;
/* L300: */
		    }

		    p2 += (len - 1) * *k;
/* L310: */
		}

		p2 -= len * *k * *l;
		i__5 = len * *k;
		i__4 = len * dimb;
		i__6 = len * dimc;
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1]
			, &i__5, &dwork[r1], &i__4, &c_b9, &dwork[s1], &i__6, 
			(ftnlen)1, (ftnlen)12);
		i__5 = len * *k;
		i__4 = len * dimb;
		i__6 = len * dimc;
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p2]
			, &i__5, &dwork[r2], &i__4, &c_b23, &dwork[s1], &i__6,
			 (ftnlen)1, (ftnlen)12);
		i__5 = len * *k;
		i__4 = len * dimb;
		i__6 = len * dimc;
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1]
			, &i__5, &dwork[r2], &i__4, &c_b9, &dwork[s2], &i__6, 
			(ftnlen)1, (ftnlen)12);
		d__1 = -scal;
		i__5 = len * *k;
		i__4 = len * dimb;
		i__6 = len * dimc;
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &d__1, &dwork[p2]
			, &i__5, &dwork[r1], &i__4, &c_b23, &dwork[s2], &i__6,
			 (ftnlen)1, (ftnlen)12);
		p2 -= *k;
		r1 += dimb;
		r2 -= dimb;
		s1 += dimc;
		s2 -= dimc;
/* L320: */
	    }

/* L330: */
	}

/*        Step 5: Hartley transform with scrambled input. */

	i__1 = pc + len * dimc * *r__;
	i__3 = len * dimc;
	for (j = pc; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {

	    i__2 = len * dimc;
	    i__5 = dimc << 1;
	    for (icp = dimc + 1; i__5 < 0 ? icp >= i__2 : icp <= i__2; icp += 
		    i__5) {
		icq = icp - dimc;

		i__4 = dimc - 1;
		for (ir = 0; ir <= i__4; ++ir) {
		    t1 = dwork[icp + ir + j];
		    dwork[icp + ir + j] = dwork[icq + ir + j] - t1;
		    dwork[icq + ir + j] += t1;
/* L340: */
		}

/* L350: */
	    }

	    ln = 1;
	    wpos = pdw + len - (p << 1) + 1;

	    i__5 = p - 1;
	    for (pp = 1; pp <= i__5; ++pp) {
		ln <<= 1;
		p2 = 1;
		q2 = ln * dimc + 1;
		r2 = ln / 2 * dimc + 1;
		s2 = r2 + q2 - 1;

		i__2 = len / (ln << 1) - 1;
		for (i__ = 0; i__ <= i__2; ++i__) {

		    i__4 = dimc - 1;
		    for (ir = 0; ir <= i__4; ++ir) {
			t1 = dwork[q2 + ir + j];
			dwork[q2 + ir + j] = dwork[p2 + ir + j] - t1;
			dwork[p2 + ir + j] += t1;
			t1 = dwork[s2 + ir + j];
			dwork[s2 + ir + j] = dwork[r2 + ir + j] - t1;
			dwork[r2 + ir + j] += t1;
/* L360: */
		    }

		    p1 = p2 + dimc;
		    q1 = p1 + ln * dimc;
		    r1 = q1 - (dimc << 1);
		    s1 = r1 + ln * dimc;

		    i__4 = wpos + ln - 3;
		    for (jj = wpos; jj <= i__4; jj += 2) {
			cf = dwork[jj];
			sf = dwork[jj + 1];

			i__6 = dimc - 1;
			for (ir = 0; ir <= i__6; ++ir) {
			    t1 = cf * dwork[q1 + ir + j] + sf * dwork[s1 + ir 
				    + j];
			    t2 = -cf * dwork[s1 + ir + j] + sf * dwork[q1 + 
				    ir + j];
			    dwork[q1 + ir + j] = dwork[p1 + ir + j] - t1;
			    dwork[p1 + ir + j] += t1;
			    dwork[s1 + ir + j] = dwork[r1 + ir + j] - t2;
			    dwork[r1 + ir + j] += t2;
/* L370: */
			}

			p1 += dimc;
			q1 += dimc;
			r1 -= dimc;
			s1 -= dimc;
/* L380: */
		    }

		    p2 += (dimc << 1) * ln;
		    q2 += (dimc << 1) * ln;
		    r2 += (dimc << 1) * ln;
		    s2 += (dimc << 1) * ln;
/* L390: */
		}

		wpos = wpos - (ln << 1) + 2;
/* L400: */
	    }

/* L410: */
	}

/*        Step 6: Copy data from workspace to output. */

L420:

	if (ltran) {
	    i__ = nl;
	} else {
	    i__ = mk;
	}

	i__3 = *r__ - 1;
	for (j = 0; j <= i__3; ++j) {
	    daxpy_(&i__, &c_b23, &dwork[pc + j * len * dimc + 1], &c__1, &c__[
		    (j + 1) * c_dim1 + 1], &c__1);
/* L430: */
	}

    }
    dwork[1] = (doublereal) max(1,wrkopt);
    return 0;

/* *** Last line of MB02KD *** */
} /* mb02kd_ */

