/* MB02CV.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02cv_(char *typeg, char *strucg, integer *k, integer *
	n, integer *p, integer *q, integer *nb, integer *rnk, doublereal *a1, 
	integer *lda1, doublereal *a2, integer *lda2, doublereal *b, integer *
	ldb, doublereal *f1, integer *ldf1, doublereal *f2, integer *ldf2, 
	doublereal *g, integer *ldg, doublereal *cs, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen typeg_len, ftnlen strucg_len)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, f1_dim1,
	     f1_offset, f2_dim1, f2_offset, g_dim1, g_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer ib, jj, nbl, len;
    static doublereal tau;
    static integer pos, col2, pst2;
    static doublereal beta;
    static logical lcol;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltri;
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical lrdef;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlarfb_(char *, char *, char 
	    *, char *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarft_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    static integer wrkmin;


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

/*     To apply the transformations created by the SLICOT Library routine */
/*     MB02CU on other columns / rows of the generator, contained in the */
/*     arrays F1, F2 and G. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPEG   CHARACTER*1 */
/*             Specifies the type of the generator, as follows: */
/*             = 'D':  generator is column oriented and rank */
/*                     deficient; */
/*             = 'C':  generator is column oriented and not rank */
/*                     deficient; */
/*             = 'R':  generator is row oriented and not rank */
/*                     deficient. */
/*             Note that this parameter must be equivalent with the */
/*             used TYPEG in the call of MB02CU. */

/*     STRUCG  CHARACTER*1 */
/*             Information about the structure of the generators, */
/*             as follows: */
/*             = 'T':  the trailing block of the positive generator */
/*                     is upper / lower triangular, and the trailing */
/*                     block of the negative generator is zero; */
/*             = 'N':  no special structure to mention. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in A1 to be processed.  K >= 0. */

/*     N       (input)  INTEGER */
/*             If TYPEG = 'D'  or  TYPEG = 'C', the number of rows in F1; */
/*             if TYPEG = 'R', the number of columns in F1.  N >= 0. */

/*     P       (input)  INTEGER */
/*             The number of columns of the positive generator.  P >= K. */

/*     Q       (input)  INTEGER */
/*             The number of columns in B. */
/*             If TYPEG = 'D',        Q >= K; */
/*             If TYPEG = 'C' or 'R', Q >= 0. */

/*     NB      (input)  INTEGER */
/*             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies */
/*             the block size to be used in the blocked parts of the */
/*             algorithm. NB must be equivalent with the used block size */
/*             in the routine MB02CU. */

/*     RNK     (input)  INTEGER */
/*             If TYPEG = 'D', the number of linearly independent columns */
/*             in the generator as returned by MB02CU.  0 <= RNK <= K. */
/*             If TYPEG = 'C' or 'R', the value of this parameter is */
/*             irrelevant. */

/*     A1      (input)  DOUBLE PRECISION array, dimension */
/*             (LDA1, K) */
/*             On entry, if TYPEG = 'D', the leading K-by-K part of this */
/*             array must contain the matrix A1 as returned by MB02CU. */
/*             If TYPEG = 'C' or 'R', this array is not referenced. */

/*     LDA1    INTEGER */
/*             The leading dimension of the array A1. */
/*             If TYPEG = 'D',                   LDA1 >= MAX(1,K); */
/*             if TYPEG = 'C'  or  TYPEG = 'R',  LDA1 >= 1. */

/*     A2      (input)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDA2, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array must contain the matrix */
/*             A2 as returned by MB02CU. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array must contain the matrix A2 as returned by */
/*             MB02CU. */

/*     LDA2    INTEGER */
/*             The leading dimension of the array A2. */
/*             If P = K,                  LDA2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDA2 >= MAX(1,K); */
/*             if P > K and TYPEG = 'R',  LDA2 >= P-K. */

/*     B       (input)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q); */
/*             if TYPEG = 'R',                   dimension (LDB, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array must contain the matrix B as */
/*             returned by MB02CU. */
/*             On entry, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array must contain the matrix B as returned by MB02CU. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             If Q = 0,                  LDB >= 1; */
/*             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDB >= MAX(1,K); */
/*             if Q > 0 and TYPEG = 'R',  LDB >= Q. */

/*     F1      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF1, K); */
/*             if TYPEG = 'R',                   dimension (LDF1, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-K part of this array must contain the first part */
/*             of the positive generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading K-by-N part of this */
/*             array must contain the first part of the positive */
/*             generator to be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-K part of this array contains the first part of the */
/*             transformed positive generator. */
/*             On exit, if TYPEG = 'R', the leading K-by-N part of this */
/*             array contains the first part of the transformed positive */
/*             generator. */

/*     LDF1    INTEGER */
/*             The leading dimension of the array F1. */
/*             If TYPEG = 'D'  or  TYPEG = 'C',   LDF1 >= MAX(1,N); */
/*             if TYPEG = 'R',                    LDF1 >= MAX(1,K). */

/*     F2      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDF2, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-(P-K) part of this array must contain the second part */
/*             of the positive generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-N part of */
/*             this array must contain the second part of the positive */
/*             generator to be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-(P-K) part of this array contains the second part of */
/*             the transformed positive generator. */
/*             On exit, if TYPEG = 'R', the leading (P-K)-by-N part of */
/*             this array contains the second part of the transformed */
/*             positive generator. */

/*     LDF2    INTEGER */
/*             The leading dimension of the array F2. */
/*             If P = K,                  LDF2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDF2 >= MAX(1,N); */
/*             if P > K and TYPEG = 'R',  LDF2 >= P-K. */

/*     G       (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDG, Q); */
/*             if TYPEG = 'R',                   dimension (LDG, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-Q part of this array must contain the negative part */
/*             of the generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading Q-by-N part of this */
/*             array must contain the negative part of the generator to */
/*             be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-Q part of this array contains the transformed */
/*             negative generator. */
/*             On exit, if TYPEG = 'R', the leading Q-by-N part of this */
/*             array contains the transformed negative generator. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G. */
/*             If Q = 0,                  LDG >= 1; */
/*             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDG >= MAX(1,N); */
/*             if Q > 0 and TYPEG = 'R',  LDG >= Q. */

/*     CS      (input)  DOUBLE PRECISION array, dimension (x) */
/*             If TYPEG = 'D' and P = K,                   x = 3*K; */
/*             If TYPEG = 'D' and P > K,                   x = 5*K; */
/*             If (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K; */
/*             If (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K. */
/*             On entry, the first x elements of this array must contain */
/*             Givens and modified hyperbolic rotation parameters, and */
/*             scalar factors of the Householder transformations as */
/*             returned by MB02CU. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = -23,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             TYPEG = 'D':               LDWORK >= MAX(1,N); */
/*             (TYPEG = 'C' or TYPEG = 'R')  and  NB <= 0: */
/*                                        LDWORK >= MAX(1,N); */
/*             (TYPEG = 'C' or TYPEG = 'R')  and  NB >= 1: */
/*                                        LDWORK >= MAX(1,( N + K )*NB). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N*K*( P + Q )) floating point operations. */

/*     METHOD */

/*     The Householder transformations and modified hyperbolic rotations */
/*     computed by SLICOT Library routine MB02CU are applied to the */
/*     corresponding parts of the matrices F1, F2 and G. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     March 2004, March 2007. */

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
    a1_dim1 = *lda1;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    a2_dim1 = *lda2;
    a2_offset = 1 + a2_dim1;
    a2 -= a2_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    f1_dim1 = *ldf1;
    f1_offset = 1 + f1_dim1;
    f1 -= f1_offset;
    f2_dim1 = *ldf2;
    f2_offset = 1 + f2_dim1;
    f2 -= f2_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --cs;
    --dwork;

    /* Function Body */
    *info = 0;
/* Computing MAX */
    i__1 = 0, i__2 = *p - *k;
    col2 = max(i__1,i__2);
    lrdef = lsame_(typeg, "D", (ftnlen)1, (ftnlen)1);
    lcol = lsame_(typeg, "C", (ftnlen)1, (ftnlen)1);
    ltri = lsame_(strucg, "T", (ftnlen)1, (ftnlen)1);
    if (lrdef) {
	wrkmin = max(1,*n);
    } else {
	if (*nb >= 1) {
/* Computing MAX */
	    i__1 = 1, i__2 = (*n + *k) * *nb;
	    wrkmin = max(i__1,i__2);
	} else {
	    wrkmin = max(1,*n);
	}
    }

/*     Check the scalar input parameters. */

    if (! (lcol || lrdef || lsame_(typeg, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (ltri || lsame_(strucg, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*k < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*p < *k) {
	*info = -5;
    } else if (*q < 0 || lrdef && *q < *k) {
	*info = -6;
    } else if (lrdef && (*rnk < 0 || *rnk > *k)) {
	*info = -8;
    } else if (*lda1 < 1 || lrdef && *lda1 < *k) {
	*info = -10;
    } else if (*p == *k && *lda2 < 1 || *p > *k && (lrdef || lcol) && *lda2 < 
	    max(1,*k) || *p > *k && ! (lrdef || lcol) && *lda2 < *p - *k) {
	*info = -12;
    } else if (*q == 0 && *ldb < 1 || *q > 0 && (lrdef || lcol) && *ldb < max(
	    1,*k) || *q > 0 && ! (lrdef || lcol) && *ldb < *q) {
	*info = -14;
    } else if ((lrdef || lcol) && *ldf1 < max(1,*n)) {
	*info = -16;
    } else if (! (lrdef || lcol) && *ldf1 < max(1,*k)) {
	*info = -16;
    } else if (*p == *k && *ldf2 < 1 || *p > *k && (lrdef || lcol) && *ldf2 < 
	    max(1,*n) || *p > *k && ! (lrdef || lcol) && *ldf2 < *p - *k) {
	*info = -18;
    } else if (*q == 0 && *ldg < 1 || *q > 0 && (lrdef || lcol) && *ldg < max(
	    1,*n) || *q > 0 && ! (lrdef || lcol) && *ldg < *q) {
	*info = -20;
    } else if (*ldwork < wrkmin) {
	dwork[1] = (doublereal) wrkmin;
	*info = -23;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02CV", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*k,*n) == 0 || ! lrdef && *q == 0 && *p == *k) {
	return 0;
    }

    if (lrdef) {

/*        Deficient generator. */

	if (col2 == 0) {
	    pst2 = *k << 1;
	} else {
	    pst2 = *k << 2;
	}

	i__1 = *rnk;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Apply elementary reflectors. */

	    if (col2 > 1) {
		tau = a2[i__ + a2_dim1];
		a2[i__ + a2_dim1] = 1.;
		dlarf_("Right", n, &col2, &a2[i__ + a2_dim1], lda2, &tau, &f2[
			f2_offset], ldf2, &dwork[1], (ftnlen)5);
		a2[i__ + a2_dim1] = tau;
	    }

	    if (*k > i__) {
		alpha = a1[i__ + i__ * a1_dim1];
		a1[i__ + i__ * a1_dim1] = 1.;
		i__2 = *k - i__ + 1;
		dlarf_("Right", n, &i__2, &a1[i__ + i__ * a1_dim1], lda1, &cs[
			pst2 + i__], &f1[i__ * f1_dim1 + 1], ldf1, &dwork[1], 
			(ftnlen)5);
		a1[i__ + i__ * a1_dim1] = alpha;
	    }

	    if (col2 > 0) {
		c__ = cs[(*k << 1) + (i__ << 1) - 1];
		s = cs[(*k << 1) + (i__ << 1)];
		drot_(n, &f1[i__ * f1_dim1 + 1], &c__1, &f2[f2_offset], &c__1,
			 &c__, &s);
	    }

	    if (*q > 1) {
		tau = b[i__ + b_dim1];
		b[i__ + b_dim1] = 1.;
		dlarf_("Right", n, q, &b[i__ + b_dim1], ldb, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)5);
		b[i__ + b_dim1] = tau;
	    }

/*           Apply hyperbolic rotation. */

	    c__ = cs[(i__ << 1) - 1];
	    s = cs[i__ * 2];
	    d__1 = 1. / c__;
	    dscal_(n, &d__1, &f1[i__ * f1_dim1 + 1], &c__1);
	    d__1 = -s / c__;
	    daxpy_(n, &d__1, &g[g_dim1 + 1], &c__1, &f1[i__ * f1_dim1 + 1], &
		    c__1);
	    dscal_(n, &c__, &g[g_dim1 + 1], &c__1);
	    d__1 = -s;
	    daxpy_(n, &d__1, &f1[i__ * f1_dim1 + 1], &c__1, &g[g_dim1 + 1], &
		    c__1);
/* L10: */
	}

	len = *q;
	pos = 1;

	i__1 = *k;
	for (j = *rnk + 1; j <= i__1; ++j) {

/*           Apply the reductions working on singular rows. */

	    if (col2 > 1) {
		tau = a2[j + a2_dim1];
		a2[j + a2_dim1] = 1.;
		dlarf_("Right", n, &col2, &a2[j + a2_dim1], lda2, &tau, &f2[
			f2_offset], ldf2, &dwork[1], (ftnlen)5);
		a2[j + a2_dim1] = tau;
	    }
	    if (*k > j) {
		alpha = a1[j + j * a1_dim1];
		a1[j + j * a1_dim1] = 1.;
		i__2 = *k - j + 1;
		dlarf_("Right", n, &i__2, &a1[j + j * a1_dim1], lda1, &cs[
			pst2 + j], &f1[j * f1_dim1 + 1], ldf1, &dwork[1], (
			ftnlen)5);
		a1[j + j * a1_dim1] = alpha;
	    }
	    if (col2 > 0) {
		c__ = cs[(*k << 1) + (j << 1) - 1];
		s = cs[(*k << 1) + (j << 1)];
		drot_(n, &f1[j * f1_dim1 + 1], &c__1, &f2[f2_offset], &c__1, &
			c__, &s);
	    }
	    if (len > 1) {
		beta = b[j + pos * b_dim1];
		b[j + pos * b_dim1] = 1.;
		dlarf_("Right", n, &len, &b[j + pos * b_dim1], ldb, &cs[(j << 
			1) - 1], &g[pos * g_dim1 + 1], ldg, &dwork[1], (
			ftnlen)5);
		b[j + pos * b_dim1] = beta;
	    }
	    --len;
	    ++pos;
/* L20: */
	}

    } else if (lcol) {

/*        Column oriented and not deficient generator. */

/*        Apply an LQ like hyperbolic/orthogonal blocked decomposition. */

	if (ltri) {
/* Computing MAX */
	    i__1 = *n - *k;
	    len = max(i__1,0);
	} else {
	    len = *n;
	}
	if (col2 > 0) {

	    nbl = min(col2,*nb);
	    if (nbl > 0) {

/*              Blocked version. */

		i__1 = *k - nbl + 1;
		i__2 = nbl;
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
		    i__3 = *k - i__ + 1;
		    ib = min(i__3,nbl);
		    i__3 = *n + *k;
		    dlarft_("Forward", "Rowwise", &col2, &ib, &a2[i__ + 
			    a2_dim1], lda2, &cs[(*k << 2) + i__], &dwork[1], &
			    i__3, (ftnlen)7, (ftnlen)7);
		    i__3 = *n + *k;
		    i__4 = *n + *k;
		    dlarfb_("Right", "No Transpose", "Forward", "Rowwise", &
			    len, &col2, &ib, &a2[i__ + a2_dim1], lda2, &dwork[
			    1], &i__3, &f2[f2_offset], ldf2, &dwork[ib + 1], &
			    i__4, (ftnlen)5, (ftnlen)12, (ftnlen)7, (ftnlen)7)
			    ;

		    i__3 = i__ + ib - 1;
		    for (j = i__; j <= i__3; ++j) {
			tau = a2[j + a2_dim1];
			a2[j + a2_dim1] = 1.;
/* Computing MIN */
			i__5 = col2, i__6 = j - i__ + 1;
			i__4 = min(i__5,i__6);
			dlarf_("Right", &len, &i__4, &a2[j + a2_dim1], lda2, &
				tau, &f2[f2_offset], ldf2, &dwork[1], (ftnlen)
				5);
			a2[j + a2_dim1] = tau;
			c__ = cs[(*k << 1) + (j << 1) - 1];
			s = cs[(*k << 1) + (j << 1)];
			drot_(&len, &f1[j * f1_dim1 + 1], &c__1, &f2[
				f2_offset], &c__1, &c__, &s);
			if (ltri) {
			    ++len;
			    temp = f1[len + j * f1_dim1];
			    f1[len + j * f1_dim1] = c__ * temp;
			    f2[len + f2_dim1] = -s * temp;

			    i__4 = col2;
			    for (jj = 2; jj <= i__4; ++jj) {
				f2[len + jj * f2_dim1] = 0.;
/* L30: */
			    }

			}
/* L40: */
		    }

/* L50: */
		}

	    } else {
		i__ = 1;
	    }

/*           Unblocked version for the last or only block. */

	    i__2 = *k;
	    for (j = i__; j <= i__2; ++j) {
		if (col2 > 1) {
		    tau = a2[j + a2_dim1];
		    a2[j + a2_dim1] = 1.;
		    dlarf_("Right", &len, &col2, &a2[j + a2_dim1], lda2, &tau,
			     &f2[f2_offset], ldf2, &dwork[1], (ftnlen)5);
		    a2[j + a2_dim1] = tau;
		}

		c__ = cs[(*k << 1) + (j << 1) - 1];
		s = cs[(*k << 1) + (j << 1)];
		drot_(&len, &f1[j * f1_dim1 + 1], &c__1, &f2[f2_offset], &
			c__1, &c__, &s);
		if (ltri) {
		    ++len;
		    temp = f1[len + j * f1_dim1];
		    f1[len + j * f1_dim1] = c__ * temp;
		    f2[len + f2_dim1] = -s * temp;

		    i__1 = col2;
		    for (jj = 2; jj <= i__1; ++jj) {
			f2[len + jj * f2_dim1] = 0.;
/* L60: */
		    }

		}
/* L70: */
	    }

	    pst2 = *k * 5;
	} else {
	    pst2 = *k << 1;
	}

	if (ltri) {
	    len = *n - *k;
	} else {
	    len = *n;
	}

	nbl = min(*q,*nb);
	if (nbl > 0) {

/*           Blocked version. */

	    i__2 = *k - nbl + 1;
	    i__1 = nbl;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *k - i__ + 1;
		ib = min(i__3,nbl);
		i__3 = *n + *k;
		dlarft_("Forward", "Rowwise", q, &ib, &b[i__ + b_dim1], ldb, &
			cs[pst2 + i__], &dwork[1], &i__3, (ftnlen)7, (ftnlen)
			7);
		i__3 = *n + *k;
		i__4 = *n + *k;
		dlarfb_("Right", "NonTranspose", "Forward", "Rowwise", &len, 
			q, &ib, &b[i__ + b_dim1], ldb, &dwork[1], &i__3, &g[
			g_offset], ldg, &dwork[ib + 1], &i__4, (ftnlen)5, (
			ftnlen)12, (ftnlen)7, (ftnlen)7);

		i__3 = i__ + ib - 1;
		for (j = i__; j <= i__3; ++j) {
		    tau = b[j + b_dim1];
		    b[j + b_dim1] = 1.;
		    i__4 = j - i__ + 1;
		    dlarf_("Right", &len, &i__4, &b[j + b_dim1], ldb, &tau, &
			    g[g_offset], ldg, &dwork[1], (ftnlen)5);
		    b[j + b_dim1] = tau;

/*                 Apply hyperbolic rotation. */

		    c__ = cs[(j << 1) - 1];
		    s = cs[j * 2];
		    d__1 = 1. / c__;
		    dscal_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1);
		    d__1 = -s / c__;
		    daxpy_(&len, &d__1, &g[g_offset], &c__1, &f1[j * f1_dim1 
			    + 1], &c__1);
		    dscal_(&len, &c__, &g[g_offset], &c__1);
		    d__1 = -s;
		    daxpy_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1, &g[
			    g_offset], &c__1);
		    if (ltri) {
			++len;
			g[len + g_dim1] = -s / c__ * f1[len + j * f1_dim1];
			f1[len + j * f1_dim1] /= c__;

			i__4 = *q;
			for (jj = 2; jj <= i__4; ++jj) {
			    g[len + jj * g_dim1] = 0.;
/* L80: */
			}

		    }
/* L90: */
		}

/* L100: */
	    }

	} else {
	    i__ = 1;
	}

/*        Unblocked version for the last or only block. */

	i__1 = *k;
	for (j = i__; j <= i__1; ++j) {
	    if (*q > 1) {
		tau = b[j + b_dim1];
		b[j + b_dim1] = 1.;
		dlarf_("Right", &len, q, &b[j + b_dim1], ldb, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)5);
		b[j + b_dim1] = tau;
	    }
	    if (*q > 0) {

/*              Apply hyperbolic rotation. */

		c__ = cs[(j << 1) - 1];
		s = cs[j * 2];
		d__1 = 1. / c__;
		dscal_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1);
		d__1 = -s / c__;
		daxpy_(&len, &d__1, &g[g_offset], &c__1, &f1[j * f1_dim1 + 1],
			 &c__1);
		dscal_(&len, &c__, &g[g_offset], &c__1);
		d__1 = -s;
		daxpy_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1, &g[g_offset],
			 &c__1);
		if (ltri) {
		    ++len;
		    g[len + g_dim1] = -s / c__ * f1[len + j * f1_dim1];
		    f1[len + j * f1_dim1] /= c__;

		    i__2 = *q;
		    for (jj = 2; jj <= i__2; ++jj) {
			g[len + jj * g_dim1] = 0.;
/* L110: */
		    }

		}
	    }
/* L120: */
	}

    } else {

/*        Row oriented and not deficient generator. */

	if (ltri) {
/* Computing MAX */
	    i__1 = *n - *k;
	    len = max(i__1,0);
	} else {
	    len = *n;
	}

	if (col2 > 0) {
	    nbl = min(*nb,col2);
	    if (nbl > 0) {

/*              Blocked version. */

		i__1 = *k - nbl + 1;
		i__2 = nbl;
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
		    i__3 = *k - i__ + 1;
		    ib = min(i__3,nbl);
		    i__3 = *n + *k;
		    dlarft_("Forward", "Columnwise", &col2, &ib, &a2[i__ * 
			    a2_dim1 + 1], lda2, &cs[(*k << 2) + i__], &dwork[
			    1], &i__3, (ftnlen)7, (ftnlen)10);
		    i__3 = *n + *k;
		    i__4 = *n + *k;
		    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &
			    col2, &len, &ib, &a2[i__ * a2_dim1 + 1], lda2, &
			    dwork[1], &i__3, &f2[f2_offset], ldf2, &dwork[ib 
			    + 1], &i__4, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			    ftnlen)10);

		    i__3 = i__ + ib - 1;
		    for (j = i__; j <= i__3; ++j) {
			tau = a2[j * a2_dim1 + 1];
			a2[j * a2_dim1 + 1] = 1.;
/* Computing MIN */
			i__5 = col2, i__6 = j - i__ + 1;
			i__4 = min(i__5,i__6);
			dlarf_("Left", &i__4, &len, &a2[j * a2_dim1 + 1], &
				c__1, &tau, &f2[f2_offset], ldf2, &dwork[1], (
				ftnlen)4);
			a2[j * a2_dim1 + 1] = tau;
			c__ = cs[(*k << 1) + (j << 1) - 1];
			s = cs[(*k << 1) + (j << 1)];
			drot_(&len, &f1[j + f1_dim1], ldf1, &f2[f2_offset], 
				ldf2, &c__, &s);
			if (ltri) {
			    ++len;
			    temp = f1[j + len * f1_dim1];
			    f1[j + len * f1_dim1] = c__ * temp;
			    f2[len * f2_dim1 + 1] = -s * temp;

			    i__4 = col2;
			    for (jj = 2; jj <= i__4; ++jj) {
				f2[jj + len * f2_dim1] = 0.;
/* L130: */
			    }

			}
/* L140: */
		    }

/* L150: */
		}

	    } else {
		i__ = 1;
	    }

/*           Unblocked version for the last or only block. */

	    i__2 = *k;
	    for (j = i__; j <= i__2; ++j) {
		if (col2 > 1) {
		    tau = a2[j * a2_dim1 + 1];
		    a2[j * a2_dim1 + 1] = 1.;
		    dlarf_("Left", &col2, &len, &a2[j * a2_dim1 + 1], &c__1, &
			    tau, &f2[f2_offset], ldf2, &dwork[1], (ftnlen)4);
		    a2[j * a2_dim1 + 1] = tau;
		}

		c__ = cs[(*k << 1) + (j << 1) - 1];
		s = cs[(*k << 1) + (j << 1)];
		drot_(&len, &f1[j + f1_dim1], ldf1, &f2[f2_offset], ldf2, &
			c__, &s);
		if (ltri) {
		    ++len;
		    temp = f1[j + len * f1_dim1];
		    f1[j + len * f1_dim1] = c__ * temp;
		    f2[len * f2_dim1 + 1] = -s * temp;

		    i__1 = col2;
		    for (jj = 2; jj <= i__1; ++jj) {
			f2[jj + len * f2_dim1] = 0.;
/* L160: */
		    }

		}
/* L170: */
	    }

	    pst2 = *k * 5;
	} else {
	    pst2 = *k << 1;
	}

	if (ltri) {
	    len = *n - *k;
	} else {
	    len = *n;
	}

	nbl = min(*q,*nb);
	if (nbl > 0) {

/*           Blocked version. */

	    i__2 = *k - nbl + 1;
	    i__1 = nbl;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *k - i__ + 1;
		ib = min(i__3,nbl);
		i__3 = *n + *k;
		dlarft_("Forward", "Columnwise", q, &ib, &b[i__ * b_dim1 + 1],
			 ldb, &cs[pst2 + i__], &dwork[1], &i__3, (ftnlen)7, (
			ftnlen)10);
		i__3 = *n + *k;
		i__4 = *n + *k;
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", q, &len,
			 &ib, &b[i__ * b_dim1 + 1], ldb, &dwork[1], &i__3, &g[
			g_offset], ldg, &dwork[ib + 1], &i__4, (ftnlen)4, (
			ftnlen)9, (ftnlen)7, (ftnlen)10);

		i__3 = i__ + ib - 1;
		for (j = i__; j <= i__3; ++j) {
		    tau = b[j * b_dim1 + 1];
		    b[j * b_dim1 + 1] = 1.;
		    i__4 = j - i__ + 1;
		    dlarf_("Left", &i__4, &len, &b[j * b_dim1 + 1], &c__1, &
			    tau, &g[g_offset], ldg, &dwork[1], (ftnlen)4);
		    b[j * b_dim1 + 1] = tau;

/*                 Apply hyperbolic rotation. */

		    c__ = cs[(j << 1) - 1];
		    s = cs[j * 2];
		    d__1 = 1. / c__;
		    dscal_(&len, &d__1, &f1[j + f1_dim1], ldf1);
		    d__1 = -s / c__;
		    daxpy_(&len, &d__1, &g[g_offset], ldg, &f1[j + f1_dim1], 
			    ldf1);
		    dscal_(&len, &c__, &g[g_offset], ldg);
		    d__1 = -s;
		    daxpy_(&len, &d__1, &f1[j + f1_dim1], ldf1, &g[g_offset], 
			    ldg);
		    if (ltri) {
			++len;
			g[len * g_dim1 + 1] = -s / c__ * f1[j + len * f1_dim1]
				;
			f1[j + len * f1_dim1] /= c__;

			i__4 = *q;
			for (jj = 2; jj <= i__4; ++jj) {
			    g[jj + len * g_dim1] = 0.;
/* L180: */
			}

		    }
/* L190: */
		}

/* L200: */
	    }

	} else {
	    i__ = 1;
	}

/*        Unblocked version for the last or only block. */

	i__1 = *k;
	for (j = i__; j <= i__1; ++j) {
	    if (*q > 1) {
		tau = b[j * b_dim1 + 1];
		b[j * b_dim1 + 1] = 1.;
		dlarf_("Left", q, &len, &b[j * b_dim1 + 1], &c__1, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)4);
		b[j * b_dim1 + 1] = tau;
	    }
	    if (*q > 0) {

/*              Apply hyperbolic rotation. */

		c__ = cs[(j << 1) - 1];
		s = cs[j * 2];
		d__1 = 1. / c__;
		dscal_(&len, &d__1, &f1[j + f1_dim1], ldf1);
		d__1 = -s / c__;
		daxpy_(&len, &d__1, &g[g_offset], ldg, &f1[j + f1_dim1], ldf1)
			;
		dscal_(&len, &c__, &g[g_offset], ldg);
		d__1 = -s;
		daxpy_(&len, &d__1, &f1[j + f1_dim1], ldf1, &g[g_offset], ldg)
			;
		if (ltri) {
		    ++len;
		    g[len * g_dim1 + 1] = -s / c__ * f1[j + len * f1_dim1];
		    f1[j + len * f1_dim1] /= c__;

		    i__2 = *q;
		    for (jj = 2; jj <= i__2; ++jj) {
			g[jj + len * g_dim1] = 0.;
/* L210: */
		    }

		}
	    }
/* L220: */
	}

    }

/* *** Last line of MB02CV *** */
    return 0;
} /* mb02cv_ */

