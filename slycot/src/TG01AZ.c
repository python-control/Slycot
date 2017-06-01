/* TG01AZ.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b18 = 10.;
static doublereal c_b56 = .5;

/* Subroutine */ int tg01az_(char *job, integer *l, integer *n, integer *m, 
	integer *p, doublereal *thresh, doublecomplex *a, integer *lda, 
	doublecomplex *e, integer *lde, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublereal *lscale, doublereal *
	rscale, doublereal *dwork, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_lg10(doublereal *), d_imag(doublecomplex *), z_abs(doublecomplex 
	    *), d_sign(doublereal *, doublereal *), pow_di(doublereal *, 
	    integer *);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer jc;
    static doublereal ta, tb, tc, te;
    static integer ir;
    static doublereal ew;
    static integer it, kw1, kw2, kw3, kw4, kw5;
    static doublereal cab, rab, ewc, cor, dum[1], sum;
    static integer nrp2, icab, lcab;
    static doublereal beta, coef;
    static integer irab, lrab;
    static doublereal basl, cmax;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5, gamma, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sfmin;
    static logical withb, withc;
    static doublereal sfmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kount;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal pgamma;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static integer lsfmin;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static integer lsfmax;


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

/*     To balance the matrices of the system pencil */

/*             S =  ( A  B ) - lambda ( E  0 ) :=  Q - lambda Z, */
/*                  ( C  0 )          ( 0  0 ) */

/*     corresponding to the descriptor triple (A-lambda E,B,C), */
/*     by balancing. This involves diagonal similarity transformations */
/*     (Dl*A*Dr - lambda Dl*E*Dr, Dl*B, C*Dr) applied to the system */
/*     (A-lambda E,B,C) to make the rows and columns of system pencil */
/*     matrices */

/*                  diag(Dl,I) * S * diag(Dr,I) */

/*     as close in norm as possible. Balancing may reduce the 1-norms */
/*     of the matrices of the system pencil S. */

/*     The balancing can be performed optionally on the following */
/*     particular system pencils */

/*              S = A-lambda E, */

/*              S = ( A-lambda E  B ),    or */

/*              S = ( A-lambda E ). */
/*                  (     C      ) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates which matrices are involved in balancing, as */
/*             follows: */
/*             = 'A':  All matrices are involved in balancing; */
/*             = 'B':  B, A and E matrices are involved in balancing; */
/*             = 'C':  C, A and E matrices are involved in balancing; */
/*             = 'N':  B and C matrices are not involved in balancing. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     THRESH  (input) DOUBLE PRECISION */
/*             Threshold value for magnitude of elements: */
/*             elements with magnitude less than or equal to */
/*             THRESH are ignored for balancing. THRESH >= 0. */
/*             The magnitude is computed as the sum of the absolute */
/*             values of the real and imaginary parts. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*A*Dr. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*E*Dr. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, if M > 0, the leading L-by-M part of this array */
/*             contains the balanced matrix Dl*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, if P > 0, the leading P-by-N part of this array */
/*             contains the balanced matrix C*Dr. */
/*             The array C is not referenced if P = 0. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     LSCALE  (output) DOUBLE PRECISION array, dimension (L) */
/*             The scaling factors applied to S from left.  If Dl(j) is */
/*             the scaling factor applied to row j, then */
/*             SCALE(j) = Dl(j), for j = 1,...,L. */

/*     RSCALE  (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to S from right.  If Dr(j) is */
/*             the scaling factor applied to column j, then */
/*             SCALE(j) = Dr(j), for j = 1,...,N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*(L+N)) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation */
/*                            -1 */
/*                  diag(Dl,I)  * S * diag(Dr,I) */

/*     to make the 1-norms of each row of the first L rows of S and its */
/*     corresponding N columns nearly equal. */

/*     Information about the diagonal matrices Dl and Dr are returned in */
/*     the vectors LSCALE and RSCALE, respectively. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     [2] R.C. Ward, R. C. */
/*         Balancing the generalized eigenvalue problem. */
/*         SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balancing, eigenvalue, matrix algebra, matrix operations, */
/*     similarity transformation. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */

/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --lscale;
    --rscale;
    --dwork;

    /* Function Body */
    *info = 0;
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*l < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*thresh < 0.) {
	*info = -6;
    } else if (*lda < max(1,*l)) {
	*info = -8;
    } else if (*lde < max(1,*l)) {
	*info = -10;
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
	*info = -12;
    } else if (*ldc < max(1,*p)) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TG01AZ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*l == 0 || *n == 0) {
	dum[0] = 1.;
	if (*l > 0) {
	    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
	} else if (*n > 0) {
	    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
	}
	return 0;
    }

/*     Initialize balancing and allocate work storage. */

    kw1 = *n;
    kw2 = kw1 + *l;
    kw3 = kw2 + *l;
    kw4 = kw3 + *n;
    kw5 = kw4 + *l;
    dum[0] = 0.;
    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
    i__1 = (*l + *n) * 3;
    dcopy_(&i__1, dum, &c__0, &dwork[1], &c__1);

/*     Compute right side vector in resulting linear equations. */

    basl = d_lg10(&c_b18);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = i__ + j * e_dim1;
	    te = (d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2));
	    i__3 = i__ + j * a_dim1;
	    ta = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2));
	    if (ta > *thresh) {
		ta = d_lg10(&ta) / basl;
	    } else {
		ta = 0.;
	    }
	    if (te > *thresh) {
		te = d_lg10(&te) / basl;
	    } else {
		te = 0.;
	    }
	    dwork[i__ + kw4] = dwork[i__ + kw4] - ta - te;
	    dwork[j + kw5] = dwork[j + kw5] - ta - te;
/* L10: */
	}
/* L20: */
    }

    if (*m == 0) {
	withb = FALSE_;
	tb = 0.;
    }
    if (*p == 0) {
	withc = FALSE_;
	tc = 0.;
    }

    if (withb) {
	i__1 = *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = izamax_(m, &b[i__ + b_dim1], ldb);
	    i__2 = i__ + j * b_dim1;
	    tb = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j * 
		    b_dim1]), abs(d__2));
	    if (tb > *thresh) {
		tb = d_lg10(&tb) / basl;
		dwork[i__ + kw4] -= tb;
	    }
/* L30: */
	}
    }

    if (withc) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__ = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
	    i__2 = i__ + j * c_dim1;
	    tc = (d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[i__ + 
		    j * c_dim1]), abs(d__2));
	    if (tc > *thresh) {
		tc = d_lg10(&tc) / basl;
		dwork[j + kw5] -= tc;
	    }
/* L40: */
	}
    }

    coef = 1. / (doublereal) (*l + *n);
    coef2 = coef * coef;
    coef5 = coef2 * .5;
    nrp2 = max(*l,*n) + 2;
    beta = 0.;
    it = 1;

/*     Start generalized conjugate gradient iteration. */

L50:

    gamma = ddot_(l, &dwork[kw4 + 1], &c__1, &dwork[kw4 + 1], &c__1) + ddot_(
	    n, &dwork[kw5 + 1], &c__1, &dwork[kw5 + 1], &c__1);

    ew = 0.;
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ew += dwork[i__ + kw4];
/* L60: */
    }

    ewc = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ewc += dwork[i__ + kw5];
/* L70: */
    }

/* Computing 2nd power */
    d__1 = ew;
/* Computing 2nd power */
    d__2 = ewc;
/* Computing 2nd power */
    d__3 = ew - ewc;
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
    if (gamma == 0.) {
	goto L160;
    }
    if (it != 1) {
	beta = gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.);
    tc = coef5 * (ew - ewc * 3.);

    i__1 = *n + *l;
    dscal_(&i__1, &beta, &dwork[1], &c__1);

    daxpy_(l, &coef, &dwork[kw4 + 1], &c__1, &dwork[kw1 + 1], &c__1);
    daxpy_(n, &coef, &dwork[kw5 + 1], &c__1, &dwork[1], &c__1);

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dwork[j] += tc;
/* L80: */
    }

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[i__ + kw1] += t;
/* L90: */
    }

/*     Apply matrix to vector. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kount = 0;
	sum = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = i__ + j * a_dim1;
	    if ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2)) > *thresh) {
		++kount;
		sum += dwork[j];
	    }
	    i__3 = i__ + j * e_dim1;
	    if ((d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2)) > *thresh) {
		++kount;
		sum += dwork[j];
	    }
/* L100: */
	}
	if (withb) {
	    j = izamax_(m, &b[i__ + b_dim1], ldb);
	    i__2 = i__ + j * b_dim1;
	    if ((d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j * 
		    b_dim1]), abs(d__2)) > *thresh) {
		++kount;
	    }
	}
	dwork[i__ + kw2] = (doublereal) kount * dwork[i__ + kw1] + sum;
/* L110: */
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	kount = 0;
	sum = 0.;
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * a_dim1;
	    if ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2)) > *thresh) {
		++kount;
		sum += dwork[i__ + kw1];
	    }
	    i__3 = i__ + j * e_dim1;
	    if ((d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2)) > *thresh) {
		++kount;
		sum += dwork[i__ + kw1];
	    }
/* L120: */
	}
	if (withc) {
	    i__ = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
	    i__2 = i__ + j * c_dim1;
	    if ((d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[i__ + j 
		    * c_dim1]), abs(d__2)) > *thresh) {
		++kount;
	    }
	}
	dwork[j + kw3] = (doublereal) kount * dwork[j] + sum;
/* L130: */
    }

    sum = ddot_(l, &dwork[kw1 + 1], &c__1, &dwork[kw2 + 1], &c__1) + ddot_(n, 
	    &dwork[1], &c__1, &dwork[kw3 + 1], &c__1);
    alpha = gamma / sum;

/*     Determine correction to current iteration. */

    cmax = 0.;
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cor = alpha * dwork[i__ + kw1];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	lscale[i__] += cor;
/* L140: */
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	cor = alpha * dwork[j];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	rscale[j] += cor;
/* L150: */
    }
    if (cmax < .5) {
	goto L160;
    }

    d__1 = -alpha;
    daxpy_(l, &d__1, &dwork[kw2 + 1], &c__1, &dwork[kw4 + 1], &c__1);
    d__1 = -alpha;
    daxpy_(n, &d__1, &dwork[kw3 + 1], &c__1, &dwork[kw5 + 1], &c__1);

    pgamma = gamma;
    ++it;
    if (it <= nrp2) {
	goto L50;
    }

/*     End generalized conjugate gradient iteration. */

L160:
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
    sfmax = 1. / sfmin;
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
    lsfmax = (integer) (d_lg10(&sfmax) / basl);

/*     Compute left diagonal scaling matrix. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	irab = izamax_(n, &a[i__ + a_dim1], lda);
	rab = z_abs(&a[i__ + irab * a_dim1]);
	irab = izamax_(n, &e[i__ + e_dim1], lde);
/* Computing MAX */
	d__1 = rab, d__2 = z_abs(&e[i__ + irab * e_dim1]);
	rab = max(d__1,d__2);
	if (withb) {
	    irab = izamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
	    d__1 = rab, d__2 = z_abs(&b[i__ + irab * b_dim1]);
	    rab = max(d__1,d__2);
	}
	d__1 = rab + sfmin;
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
	ir = (integer) (lscale[i__] + d_sign(&c_b56, &lscale[i__]));
/* Computing MIN */
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
	ir = min(i__2,i__3);
	lscale[i__] = pow_di(&c_b18, &ir);
/* L170: */
    }

/*     Compute right diagonal scaling matrix. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	icab = izamax_(l, &a[j * a_dim1 + 1], &c__1);
	cab = z_abs(&a[icab + j * a_dim1]);
	icab = izamax_(l, &e[j * e_dim1 + 1], &c__1);
/* Computing MAX */
	d__1 = cab, d__2 = z_abs(&e[icab + j * e_dim1]);
	cab = max(d__1,d__2);
	if (withc) {
	    icab = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
/* Computing MAX */
	    d__1 = cab, d__2 = z_abs(&c__[icab + j * c_dim1]);
	    cab = max(d__1,d__2);
	}
	d__1 = cab + sfmin;
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
	jc = (integer) (rscale[j] + d_sign(&c_b56, &rscale[j]));
/* Computing MIN */
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
	jc = min(i__2,i__3);
	rscale[j] = pow_di(&c_b18, &jc);
/* L180: */
    }

/*     Row scaling of matrices A, E and B. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zdscal_(n, &lscale[i__], &a[i__ + a_dim1], lda);
	zdscal_(n, &lscale[i__], &e[i__ + e_dim1], lde);
	if (withb) {
	    zdscal_(m, &lscale[i__], &b[i__ + b_dim1], ldb);
	}
/* L190: */
    }

/*     Column scaling of matrices A, E and C. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	zdscal_(l, &rscale[j], &a[j * a_dim1 + 1], &c__1);
	zdscal_(l, &rscale[j], &e[j * e_dim1 + 1], &c__1);
	if (withc) {
	    zdscal_(p, &rscale[j], &c__[j * c_dim1 + 1], &c__1);
	}
/* L200: */
    }

    return 0;
/* *** Last line of TG01AZ *** */
} /* tg01az_ */

