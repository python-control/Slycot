/* MB02CU.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02cu_(char *typeg, integer *k, integer *p, integer *q, 
	integer *nb, doublereal *a1, integer *lda1, doublereal *a2, integer *
	lda2, doublereal *b, integer *ldb, integer *rnk, integer *ipvt, 
	doublereal *cs, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen typeg_len)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer ib, jj, nbl, len, pdw, phv, pos, pvt, col2;
    static doublereal tau1, tau2;
    static integer pst2;
    static doublereal beta;
    static logical lcol;
    static doublereal dmax__;
    static integer imax, ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int ma02fd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical lrdef;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal alpha2;
    extern /* Subroutine */ int dgelq2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dgeqr2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), dlartg_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);
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

/*     To bring the first blocks of a generator to proper form. */
/*     The positive part of the generator is contained in the arrays A1 */
/*     and A2. The negative part of the generator is contained in B. */
/*     Transformation information will be stored and can be applied */
/*     via SLICOT Library routine MB02CV. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPEG   CHARACTER*1 */
/*             Specifies the type of the generator, as follows: */
/*             = 'D':  generator is column oriented and rank */
/*                     deficiencies are expected; */
/*             = 'C':  generator is column oriented and rank */
/*                     deficiencies are not expected; */
/*             = 'R':  generator is row oriented and rank */
/*                     deficiencies are not expected. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in A1 to be processed.  K >= 0. */

/*     P       (input)  INTEGER */
/*             The number of columns of the positive generator.  P >= K. */

/*     Q       (input)  INTEGER */
/*             The number of columns in B containing the negative */
/*             generators. */
/*             If TYPEG = 'D',        Q >= K; */
/*             If TYPEG = 'C' or 'R', Q >= 0. */

/*     NB      (input)  INTEGER */
/*             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies */
/*             the block size to be used in the blocked parts of the */
/*             algorithm. If NB <= 0, an unblocked algorithm is used. */

/*     A1      (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDA1, K) */
/*             On entry, the leading K-by-K part of this array must */
/*             contain the leading submatrix of the positive part of the */
/*             generator. If TYPEG = 'C', A1 is assumed to be lower */
/*             triangular and the strictly upper triangular part is not */
/*             referenced. If TYPEG = 'R', A1 is assumed to be upper */
/*             triangular and the strictly lower triangular part is not */
/*             referenced. */
/*             On exit, if TYPEG = 'D', the leading K-by-RNK part of this */
/*             array contains the lower trapezoidal part of the proper */
/*             generator and information for the Householder */
/*             transformations applied during the reduction process. */
/*             On exit, if TYPEG = 'C', the leading K-by-K part of this */
/*             array contains the leading lower triangular part of the */
/*             proper generator. */
/*             On exit, if TYPEG = 'R', the leading K-by-K part of this */
/*             array contains the leading upper triangular part of the */
/*             proper generator. */

/*     LDA1    INTEGER */
/*             The leading dimension of the array A1.  LDA1 >= MAX(1,K). */

/*     A2      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDA2, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array must contain the (K+1)-st */
/*             to P-th columns of the positive part of the generator. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array must contain the (K+1)-st to P-th rows of the */
/*             positive part of the generator. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array contains information for */
/*             Householder transformations. */
/*             On exit, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array contains information for Householder */
/*             transformations. */

/*     LDA2    INTEGER */
/*             The leading dimension of the array A2. */
/*             If P = K,                   LDA2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                         LDA2 >= MAX(1,K); */
/*             if P > K and TYPEG = 'R',   LDA2 >= P-K. */

/*     B       (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q); */
/*             if TYPEG = 'R',                   dimension (LDB, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array must contain the negative part */
/*             of the generator. */
/*             On entry, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array must contain the negative part of the generator. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array contains information for */
/*             Householder transformations. */
/*             On exit, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array contains information for Householder transformations. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             If Q = 0,                  LDB >= 1; */
/*             if Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDB >= MAX(1,K); */
/*             if Q > 0 and TYPEG = 'R',  LDB >= Q. */

/*     RNK     (output)  INTEGER */
/*             If TYPEG = 'D', the number of columns in the reduced */
/*             generator which are found to be linearly independent. */
/*             If TYPEG = 'C' or TYPEG = 'R', then RNK is not set. */

/*     IPVT    (output)  INTEGER array, dimension (K) */
/*             If TYPEG = 'D', then if IPVT(i) = k, the k-th row of the */
/*             proper generator is the reduced i-th row of the input */
/*             generator. */
/*             If TYPEG = 'C' or TYPEG = 'R', this array is not */
/*             referenced. */

/*     CS      (output)  DOUBLE PRECISION array, dimension (x) */
/*             If TYPEG = 'D' and P = K,                   x = 3*K; */
/*             if TYPEG = 'D' and P > K,                   x = 5*K; */
/*             if (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K; */
/*             if (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K. */
/*             On exit, the first x elements of this array contain */
/*             necessary information for the SLICOT library routine */
/*             MB02CV (Givens and modified hyperbolic rotation */
/*             parameters, scalar factors of the Householder */
/*             transformations). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If TYPEG = 'D', this number specifies the used tolerance */
/*             for handling deficiencies. If the hyperbolic norm */
/*             of two diagonal elements in the positive and negative */
/*             generators appears to be less than or equal to TOL, then */
/*             the corresponding columns are not reduced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = -17,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,4*K),         if TYPEG = 'D'; */
/*             LDWORK >= MAX(1,MAX(NB,1)*K), if TYPEG = 'C' or 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if TYPEG = 'D', the generator represents a */
/*                   (numerically) indefinite matrix; and if TYPEG = 'C' */
/*                   or TYPEG = 'R', the generator represents a */
/*                   (numerically) semidefinite matrix. */

/*     METHOD */

/*     If TYPEG = 'C' or TYPEG = 'R', blocked Householder transformations */
/*     and modified hyperbolic rotations are used to downdate the */
/*     matrix [ A1  A2  sqrt(-1)*B ], cf. [1], [2]. */
/*     If TYPEG = 'D', then an algorithm with row pivoting is used. In */
/*     the first stage it maximizes the hyperbolic norm of the active */
/*     row. As soon as the hyperbolic norm is below the threshold TOL, */
/*     the strategy is changed. Now, in the second stage, the algorithm */
/*     applies an LQ decomposition with row pivoting on B such that */
/*     the Euclidean norm of the active row is maximized. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(K *( P + Q )) floating point operations. */

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
    a1_dim1 = *lda1;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    a2_dim1 = *lda2;
    a2_offset = 1 + a2_dim1;
    a2 -= a2_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --ipvt;
    --cs;
    --dwork;

    /* Function Body */
    *info = 0;
    col2 = *p - *k;
    lrdef = lsame_(typeg, "D", (ftnlen)1, (ftnlen)1);
    lcol = lsame_(typeg, "C", (ftnlen)1, (ftnlen)1);
    if (lrdef) {
/* Computing MAX */
	i__1 = 1, i__2 = *k << 2;
	wrkmin = max(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = *nb * *k, i__1 = max(i__1,i__2);
	wrkmin = max(i__1,*k);
    }

/*     Check the scalar input parameters. */

    if (! (lcol || lrdef || lsame_(typeg, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*k < 0) {
	*info = -2;
    } else if (*p < *k) {
	*info = -3;
    } else if (*q < 0 || lrdef && *q < *k) {
	*info = -4;
    } else if (*lda1 < max(1,*k)) {
	*info = -7;
    } else if (*p == *k && *lda2 < 1 || *p > *k && (lrdef || lcol) && *lda2 < 
	    max(1,*k) || *p > *k && ! (lrdef || lcol) && *lda2 < *p - *k) {
	*info = -9;
    } else if (*q == 0 && *ldb < 1 || *q > 0 && (lrdef || lcol) && *ldb < max(
	    1,*k) || *q > 0 && ! (lrdef || lcol) && *ldb < *q) {
	*info = -11;
    } else if (*ldwork < wrkmin) {
	dwork[1] = (doublereal) wrkmin;
	*info = -17;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02CU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*k == 0 || ! lrdef && *q == 0 && *p == *k) {
	if (lrdef) {
	    *rnk = 0;
	}
	return 0;
    }

    if (lrdef) {

/*        Deficient generator. */

	if (col2 == 0) {
	    pst2 = *k << 1;
	} else {
	    pst2 = *k << 2;
	}

/*        Initialize partial hyperbolic row norms. */

	*rnk = 0;
	phv = *k * 3;

	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ipvt[i__] = i__;
	    dwork[i__] = dnrm2_(k, &a1[i__ + a1_dim1], lda1);
/* L10: */
	}

	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = dnrm2_(&col2, &a2[i__ + a2_dim1], lda2);
	    dwork[i__] = dlapy2_(&dwork[i__], &d__1);
	    dwork[i__ + *k] = dwork[i__];
/* L20: */
	}

	pdw = *k << 1;

	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++pdw;
	    dwork[pdw] = dnrm2_(q, &b[i__ + b_dim1], ldb);
/* L30: */
	}

/*        Compute factorization. */

	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Determine pivot row and swap if necessary. */

	    pdw = i__;
	    alpha = (d__1 = dwork[pdw], abs(d__1));
	    beta = (d__1 = dwork[pdw + (*k << 1)], abs(d__1));
	    d__2 = sqrt((d__1 = alpha - beta, abs(d__1))) * sqrt(alpha + beta)
		    ;
	    d__3 = alpha - beta;
	    dmax__ = d_sign(&d__2, &d__3);
	    imax = i__;

	    i__2 = *k - i__;
	    for (j = 1; j <= i__2; ++j) {
		++pdw;
		alpha = (d__1 = dwork[pdw], abs(d__1));
		beta = (d__1 = dwork[pdw + (*k << 1)], abs(d__1));
		d__2 = sqrt((d__1 = alpha - beta, abs(d__1))) * sqrt(alpha + 
			beta);
		d__3 = alpha - beta;
		temp = d_sign(&d__2, &d__3);
		if (temp > dmax__) {
		    imax = i__ + j;
		    dmax__ = temp;
		}
/* L40: */
	    }

/*           Proceed with the reduction if the hyperbolic norm is */
/*           beyond the threshold. */

	    if (dmax__ > *tol) {

		pvt = imax;
		if (pvt != i__) {
		    dswap_(k, &a1[pvt + a1_dim1], lda1, &a1[i__ + a1_dim1], 
			    lda1);
		    dswap_(&col2, &a2[pvt + a2_dim1], lda2, &a2[i__ + a2_dim1]
			    , lda2);
		    dswap_(q, &b[pvt + b_dim1], ldb, &b[i__ + b_dim1], ldb);
		    itemp = ipvt[pvt];
		    ipvt[pvt] = ipvt[i__];
		    ipvt[i__] = itemp;
		    dwork[pvt] = dwork[i__];
		    dwork[*k + pvt] = dwork[*k + i__];
		    dwork[(*k << 1) + pvt] = dwork[(*k << 1) + i__];
		}

/*              Generate and apply elementary reflectors. */

		if (col2 > 1) {
		    dlarfg_(&col2, &a2[i__ + a2_dim1], &a2[i__ + (a2_dim1 << 
			    1)], lda2, &tau2);
		    alpha2 = a2[i__ + a2_dim1];
		    if (*k > i__) {
			a2[i__ + a2_dim1] = 1.;
			i__2 = *k - i__;
			dlarf_("Right", &i__2, &col2, &a2[i__ + a2_dim1], 
				lda2, &tau2, &a2[i__ + 1 + a2_dim1], lda2, &
				dwork[phv + 1], (ftnlen)5);
		    }
		    a2[i__ + a2_dim1] = tau2;
		} else if (col2 > 0) {
		    alpha2 = a2[i__ + a2_dim1];
		    a2[i__ + a2_dim1] = 0.;
		}

		if (*k > i__) {
		    i__2 = *k - i__ + 1;
		    dlarfg_(&i__2, &a1[i__ + i__ * a1_dim1], &a1[i__ + (i__ + 
			    1) * a1_dim1], lda1, &tau1);
		    alpha = a1[i__ + i__ * a1_dim1];
		    a1[i__ + i__ * a1_dim1] = 1.;
		    i__2 = *k - i__;
		    i__3 = *k - i__ + 1;
		    dlarf_("Right", &i__2, &i__3, &a1[i__ + i__ * a1_dim1], 
			    lda1, &tau1, &a1[i__ + 1 + i__ * a1_dim1], lda1, &
			    dwork[phv + 1], (ftnlen)5);
		    cs[pst2 + i__] = tau1;
		} else {
		    alpha = a1[i__ + i__ * a1_dim1];
		}

		if (col2 > 0) {
		    temp = alpha;
		    dlartg_(&temp, &alpha2, &c__, &s, &alpha);
		    if (*k > i__) {
			i__2 = *k - i__;
			drot_(&i__2, &a1[i__ + 1 + i__ * a1_dim1], &c__1, &a2[
				i__ + 1 + a2_dim1], &c__1, &c__, &s);
		    }
		    cs[(*k << 1) + (i__ << 1) - 1] = c__;
		    cs[(*k << 1) + (i__ << 1)] = s;
		}
		a1[i__ + i__ * a1_dim1] = alpha;

		if (*q > 1) {
		    dlarfg_(q, &b[i__ + b_dim1], &b[i__ + (b_dim1 << 1)], ldb,
			     &tau2);
		    beta = b[i__ + b_dim1];
		    if (*k > i__) {
			b[i__ + b_dim1] = 1.;
			i__2 = *k - i__;
			dlarf_("Right", &i__2, q, &b[i__ + b_dim1], ldb, &
				tau2, &b[i__ + 1 + b_dim1], ldb, &dwork[phv + 
				1], (ftnlen)5);
		    }
		    b[i__ + b_dim1] = tau2;
		} else if (*q > 0) {
		    beta = b[i__ + b_dim1];
		    b[i__ + b_dim1] = 0.;
		} else {
		    beta = 0.;
		}

/*              Create hyperbolic Givens rotation. */

		ma02fd_(&a1[i__ + i__ * a1_dim1], &beta, &c__, &s, &ierr);
		if (ierr != 0) {

/*                 Error return:  This should not happen. */

		    *info = 1;
		    return 0;
		}

/*              Apply hyperbolic rotation. */

		if (*k > i__) {
		    i__2 = *k - i__;
		    d__1 = 1. / c__;
		    dscal_(&i__2, &d__1, &a1[i__ + 1 + i__ * a1_dim1], &c__1);
		    i__2 = *k - i__;
		    d__1 = -s / c__;
		    daxpy_(&i__2, &d__1, &b[i__ + 1 + b_dim1], &c__1, &a1[i__ 
			    + 1 + i__ * a1_dim1], &c__1);
		    i__2 = *k - i__;
		    dscal_(&i__2, &c__, &b[i__ + 1 + b_dim1], &c__1);
		    i__2 = *k - i__;
		    d__1 = -s;
		    daxpy_(&i__2, &d__1, &a1[i__ + 1 + i__ * a1_dim1], &c__1, 
			    &b[i__ + 1 + b_dim1], &c__1);
		}
		cs[(i__ << 1) - 1] = c__;
		cs[i__ * 2] = s;

/*              Downdate the norms in A1. */

		i__2 = *k;
		for (j = i__ + 1; j <= i__2; ++j) {
/* Computing 2nd power */
		    d__2 = (d__1 = a1[j + i__ * a1_dim1], abs(d__1)) / dwork[
			    j];
		    temp = 1. - d__2 * d__2;
/* Computing 2nd power */
		    d__1 = dwork[j] / dwork[*k + j];
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			i__3 = *k - i__;
			d__1 = dnrm2_(&i__3, &a1[j + (i__ + 1) * a1_dim1], 
				lda1);
			d__2 = dnrm2_(&col2, &a2[j + a2_dim1], lda2);
			dwork[j] = dlapy2_(&d__1, &d__2);
			dwork[*k + j] = dwork[j];
			dwork[(*k << 1) + j] = dnrm2_(q, &b[j + b_dim1], ldb);
		    } else {
			if (temp >= 0.) {
			    dwork[j] *= sqrt(temp);
			} else {
			    dwork[j] = -dwork[j] * sqrt(-temp);
			}
		    }
/* L50: */
		}

		++(*rnk);
	    } else if (abs(dmax__) < *tol) {

/*              Displacement is positive semidefinite. */
/*              Do an LQ decomposition with pivoting of the leftover */
/*              negative part to find diagonal elements with almost zero */
/*              norm. These columns cannot be removed from the */
/*              generator. */

/*              Initialize norms. */

		i__2 = *k;
		for (j = i__; j <= i__2; ++j) {
		    dwork[j] = dnrm2_(q, &b[j + b_dim1], ldb);
		    dwork[j + *k] = dwork[j];
/* L60: */
		}

		len = *q;
		pos = 1;

		i__2 = *k;
		for (j = i__; j <= i__2; ++j) {

/*                 Generate and apply elementary reflectors. */

		    i__3 = *k - j + 1;
		    pvt = j - 1 + idamax_(&i__3, &dwork[j], &c__1);

/*                 Swap rows if necessary. */

		    if (pvt != j) {
			dswap_(k, &a1[pvt + a1_dim1], lda1, &a1[j + a1_dim1], 
				lda1);
			dswap_(&col2, &a2[pvt + a2_dim1], lda2, &a2[j + 
				a2_dim1], lda2);
			dswap_(q, &b[pvt + b_dim1], ldb, &b[j + b_dim1], ldb);
			itemp = ipvt[pvt];
			ipvt[pvt] = ipvt[j];
			ipvt[j] = itemp;
			dwork[pvt] = dwork[j];
			dwork[*k + pvt] = dwork[*k + j];
		    }

/*                 Annihilate second part of the positive generators. */

		    if (col2 > 1) {
			dlarfg_(&col2, &a2[j + a2_dim1], &a2[j + (a2_dim1 << 
				1)], lda2, &tau2);
			alpha2 = a2[j + a2_dim1];
			if (*k > j) {
			    a2[j + a2_dim1] = 1.;
			    i__3 = *k - j;
			    dlarf_("Right", &i__3, &col2, &a2[j + a2_dim1], 
				    lda2, &tau2, &a2[j + 1 + a2_dim1], lda2, &
				    dwork[phv + 1], (ftnlen)5);
			}
			a2[j + a2_dim1] = tau2;
		    } else if (col2 > 0) {
			alpha2 = a2[j + a2_dim1];
			a2[j + a2_dim1] = 0.;
		    }

/*                 Transform first part of the positive generators to */
/*                 lower triangular form. */

		    if (*k > j) {
			i__3 = *k - j + 1;
			dlarfg_(&i__3, &a1[j + j * a1_dim1], &a1[j + (j + 1) *
				 a1_dim1], lda1, &tau1);
			alpha = a1[j + j * a1_dim1];
			a1[j + j * a1_dim1] = 1.;
			i__3 = *k - j;
			i__4 = *k - j + 1;
			dlarf_("Right", &i__3, &i__4, &a1[j + j * a1_dim1], 
				lda1, &tau1, &a1[j + 1 + j * a1_dim1], lda1, &
				dwork[phv + 1], (ftnlen)5);
			cs[pst2 + j] = tau1;
		    } else {
			alpha = a1[j + j * a1_dim1];
		    }

		    if (col2 > 0) {
			temp = alpha;
			dlartg_(&temp, &alpha2, &c__, &s, &alpha);
			if (*k > j) {
			    i__3 = *k - j;
			    drot_(&i__3, &a1[j + 1 + j * a1_dim1], &c__1, &a2[
				    j + 1 + a2_dim1], &c__1, &c__, &s);
			}
			cs[(*k << 1) + (j << 1) - 1] = c__;
			cs[(*k << 1) + (j << 1)] = s;
		    }
		    a1[j + j * a1_dim1] = alpha;

/*                 Transform negative part to lower triangular form. */

		    if (len > 1) {
			dlarfg_(&len, &b[j + pos * b_dim1], &b[j + (pos + 1) *
				 b_dim1], ldb, &tau2);
			beta = b[j + pos * b_dim1];
			if (*k > j) {
			    b[j + pos * b_dim1] = 1.;
			    i__3 = *k - j;
			    dlarf_("Right", &i__3, &len, &b[j + pos * b_dim1],
				     ldb, &tau2, &b[j + 1 + pos * b_dim1], 
				    ldb, &dwork[phv + 1], (ftnlen)5);
			}
			b[j + pos * b_dim1] = beta;
			cs[(j << 1) - 1] = tau2;
		    }

/*                 Downdate the norms of the rows in the negative part. */

		    i__3 = *k;
		    for (jj = j + 1; jj <= i__3; ++jj) {
			if (dwork[jj] != 0.) {
/* Computing 2nd power */
			    d__2 = (d__1 = b[jj + pos * b_dim1], abs(d__1)) / 
				    dwork[jj];
			    temp = 1. - d__2 * d__2;
			    temp = max(temp,0.);
/* Computing 2nd power */
			    d__1 = dwork[jj] / dwork[*k + jj];
			    temp2 = temp * .05 * (d__1 * d__1) + 1.;
			    if (temp2 == 1.) {
				i__4 = len - 1;
				dwork[jj] = dnrm2_(&i__4, &b[jj + (pos + 1) * 
					b_dim1], ldb);
				dwork[*k + jj] = dwork[jj];
			    } else {
				if (temp >= 0.) {
				    dwork[jj] *= sqrt(temp);
				} else {
				    dwork[jj] = -dwork[jj] * sqrt(-temp);
				}
			    }
			}
/* L70: */
		    }

		    --len;
		    ++pos;
/* L80: */
		}

		return 0;
	    } else {

/*              Error return: */

/*              Displacement is indefinite. */
/*              Due to roundoff error, positive semidefiniteness is */
/*              violated. This is a rather bad situation. There is no */
/*              meaningful way to continue the computations from this */
/*              point. */

		*info = 1;
		return 0;
	    }
/* L90: */
	}

    } else if (lcol) {

/*        Column oriented and not deficient generator. */

/*        Apply an LQ like hyperbolic/orthogonal blocked decomposition. */

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
		    dgelq2_(&ib, &col2, &a2[i__ + a2_dim1], lda2, &cs[(*k << 
			    2) + i__], &dwork[1], &ierr);
		    if (i__ + ib <= *k) {
			dlarft_("Forward", "Rowwise", &col2, &ib, &a2[i__ + 
				a2_dim1], lda2, &cs[(*k << 2) + i__], &dwork[
				1], k, (ftnlen)7, (ftnlen)7);
			i__3 = *k - i__ - ib + 1;
			dlarfb_("Right", "No Transpose", "Forward", "Rowwise",
				 &i__3, &col2, &ib, &a2[i__ + a2_dim1], lda2, 
				&dwork[1], k, &a2[i__ + ib + a2_dim1], lda2, &
				dwork[ib + 1], k, (ftnlen)5, (ftnlen)12, (
				ftnlen)7, (ftnlen)7);
		    }

/*                 Annihilate the remaining parts of A2. */

		    i__3 = i__ + ib - 1;
		    for (j = i__; j <= i__3; ++j) {
			if (col2 > 1) {
/* Computing MIN */
			    i__4 = col2, i__5 = j - i__ + 1;
			    len = min(i__4,i__5);
			    dlarfg_(&len, &a2[j + a2_dim1], &a2[j + (a2_dim1 
				    << 1)], lda2, &tau2);
			    alpha2 = a2[j + a2_dim1];
			    if (*k > j) {
				a2[j + a2_dim1] = 1.;
				i__4 = *k - j;
				dlarf_("Right", &i__4, &len, &a2[j + a2_dim1],
					 lda2, &tau2, &a2[j + 1 + a2_dim1], 
					lda2, &dwork[1], (ftnlen)5);
			    }
			    a2[j + a2_dim1] = tau2;
			} else {
			    alpha2 = a2[j + a2_dim1];
			    a2[j + a2_dim1] = 0.;
			}
			alpha = a1[j + j * a1_dim1];
			dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * 
				a1_dim1]);
			if (*k > j) {
			    i__4 = *k - j;
			    drot_(&i__4, &a1[j + 1 + j * a1_dim1], &c__1, &a2[
				    j + 1 + a2_dim1], &c__1, &c__, &s);
			}
			cs[(*k << 1) + (j << 1) - 1] = c__;
			cs[(*k << 1) + (j << 1)] = s;
/* L100: */
		    }

/* L110: */
		}

	    } else {
		i__ = 1;
	    }

/*           Unblocked version for the last or only block. */

	    i__2 = *k;
	    for (j = i__; j <= i__2; ++j) {
		if (col2 > 1) {
		    dlarfg_(&col2, &a2[j + a2_dim1], &a2[j + (a2_dim1 << 1)], 
			    lda2, &tau2);
		    alpha2 = a2[j + a2_dim1];
		    if (*k > j) {
			a2[j + a2_dim1] = 1.;
			i__1 = *k - j;
			dlarf_("Right", &i__1, &col2, &a2[j + a2_dim1], lda2, 
				&tau2, &a2[j + 1 + a2_dim1], lda2, &dwork[1], 
				(ftnlen)5);
		    }
		    a2[j + a2_dim1] = tau2;
		} else {
		    alpha2 = a2[j + a2_dim1];
		    a2[j + a2_dim1] = 0.;
		}
		alpha = a1[j + j * a1_dim1];
		dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * a1_dim1]);
		if (*k > j) {
		    i__1 = *k - j;
		    drot_(&i__1, &a1[j + 1 + j * a1_dim1], &c__1, &a2[j + 1 + 
			    a2_dim1], &c__1, &c__, &s);
		}
		cs[(*k << 1) + (j << 1) - 1] = c__;
		cs[(*k << 1) + (j << 1)] = s;
/* L120: */
	    }

	    pst2 = *k * 5;
	} else {
	    pst2 = *k << 1;
	}

/*        Annihilate B with hyperbolic transformations. */

	nbl = min(*nb,*q);
	if (nbl > 0) {

/*           Blocked version. */

	    i__2 = *k - nbl + 1;
	    i__1 = nbl;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *k - i__ + 1;
		ib = min(i__3,nbl);
		dgelq2_(&ib, q, &b[i__ + b_dim1], ldb, &cs[pst2 + i__], &
			dwork[1], &ierr);
		if (i__ + ib <= *k) {
		    dlarft_("Forward", "Rowwise", q, &ib, &b[i__ + b_dim1], 
			    ldb, &cs[pst2 + i__], &dwork[1], k, (ftnlen)7, (
			    ftnlen)7);
		    i__3 = *k - i__ - ib + 1;
		    dlarfb_("Right", "No Transpose", "Forward", "Rowwise", &
			    i__3, q, &ib, &b[i__ + b_dim1], ldb, &dwork[1], k,
			     &b[i__ + ib + b_dim1], ldb, &dwork[ib + 1], k, (
			    ftnlen)5, (ftnlen)12, (ftnlen)7, (ftnlen)7);
		}

/*              Annihilate the remaining parts of B. */

		i__3 = i__ + ib - 1;
		for (j = i__; j <= i__3; ++j) {
		    if (*q > 1) {
			i__4 = j - i__ + 1;
			dlarfg_(&i__4, &b[j + b_dim1], &b[j + (b_dim1 << 1)], 
				ldb, &tau2);
			alpha2 = b[j + b_dim1];
			if (*k > j) {
			    b[j + b_dim1] = 1.;
			    i__4 = *k - j;
			    i__5 = j - i__ + 1;
			    dlarf_("Right", &i__4, &i__5, &b[j + b_dim1], ldb,
				     &tau2, &b[j + 1 + b_dim1], ldb, &dwork[1]
				    , (ftnlen)5);
			}
			b[j + b_dim1] = tau2;
		    } else {
			alpha2 = b[j + b_dim1];
			b[j + b_dim1] = 0.;
		    }

/*                 Create hyperbolic rotation. */

		    ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
		    if (ierr != 0) {

/*                    Error return:  The matrix is not positive definite. */

			*info = 1;
			return 0;
		    }

/*                 Apply hyperbolic rotation. */

		    if (*k > j) {
			i__4 = *k - j;
			d__1 = 1. / c__;
			dscal_(&i__4, &d__1, &a1[j + 1 + j * a1_dim1], &c__1);
			i__4 = *k - j;
			d__1 = -s / c__;
			daxpy_(&i__4, &d__1, &b[j + 1 + b_dim1], &c__1, &a1[j 
				+ 1 + j * a1_dim1], &c__1);
			i__4 = *k - j;
			dscal_(&i__4, &c__, &b[j + 1 + b_dim1], &c__1);
			i__4 = *k - j;
			d__1 = -s;
			daxpy_(&i__4, &d__1, &a1[j + 1 + j * a1_dim1], &c__1, 
				&b[j + 1 + b_dim1], &c__1);
		    }
		    cs[(j << 1) - 1] = c__;
		    cs[j * 2] = s;
/* L130: */
		}

/* L140: */
	    }

	} else {
	    i__ = 1;
	}

/*        Unblocked version for the last or only block. */

	i__1 = *k;
	for (j = i__; j <= i__1; ++j) {
	    if (*q > 1) {
		dlarfg_(q, &b[j + b_dim1], &b[j + (b_dim1 << 1)], ldb, &tau2);
		alpha2 = b[j + b_dim1];
		if (*k > j) {
		    b[j + b_dim1] = 1.;
		    i__2 = *k - j;
		    dlarf_("Right", &i__2, q, &b[j + b_dim1], ldb, &tau2, &b[
			    j + 1 + b_dim1], ldb, &dwork[1], (ftnlen)5);
		}
		b[j + b_dim1] = tau2;
	    } else if (*q > 0) {
		alpha2 = b[j + b_dim1];
		b[j + b_dim1] = 0.;
	    }
	    if (*q > 0) {

/*              Create hyperbolic rotation. */

		ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/*              Apply hyperbolic rotation. */

		if (*k > j) {
		    i__2 = *k - j;
		    d__1 = 1. / c__;
		    dscal_(&i__2, &d__1, &a1[j + 1 + j * a1_dim1], &c__1);
		    i__2 = *k - j;
		    d__1 = -s / c__;
		    daxpy_(&i__2, &d__1, &b[j + 1 + b_dim1], &c__1, &a1[j + 1 
			    + j * a1_dim1], &c__1);
		    i__2 = *k - j;
		    dscal_(&i__2, &c__, &b[j + 1 + b_dim1], &c__1);
		    i__2 = *k - j;
		    d__1 = -s;
		    daxpy_(&i__2, &d__1, &a1[j + 1 + j * a1_dim1], &c__1, &b[
			    j + 1 + b_dim1], &c__1);
		}
		cs[(j << 1) - 1] = c__;
		cs[j * 2] = s;
	    }
/* L150: */
	}

    } else {

/*        Row oriented and not deficient generator. */

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
		    dgeqr2_(&col2, &ib, &a2[i__ * a2_dim1 + 1], lda2, &cs[(*k 
			    << 2) + i__], &dwork[1], &ierr);
		    if (i__ + ib <= *k) {
			dlarft_("Forward", "Columnwise", &col2, &ib, &a2[i__ *
				 a2_dim1 + 1], lda2, &cs[(*k << 2) + i__], &
				dwork[1], k, (ftnlen)7, (ftnlen)10);
			i__3 = *k - i__ - ib + 1;
			dlarfb_("Left", "Transpose", "Forward", "Columnwise", 
				&col2, &i__3, &ib, &a2[i__ * a2_dim1 + 1], 
				lda2, &dwork[1], k, &a2[(i__ + ib) * a2_dim1 
				+ 1], lda2, &dwork[ib + 1], k, (ftnlen)4, (
				ftnlen)9, (ftnlen)7, (ftnlen)10);
		    }

/*                 Annihilate the remaining parts of A2. */

		    i__3 = i__ + ib - 1;
		    for (j = i__; j <= i__3; ++j) {
			if (col2 > 1) {
/* Computing MIN */
			    i__4 = col2, i__5 = j - i__ + 1;
			    len = min(i__4,i__5);
			    dlarfg_(&len, &a2[j * a2_dim1 + 1], &a2[j * 
				    a2_dim1 + 2], &c__1, &tau2);
			    alpha2 = a2[j * a2_dim1 + 1];
			    if (*k > j) {
				a2[j * a2_dim1 + 1] = 1.;
				i__4 = *k - j;
				dlarf_("Left", &len, &i__4, &a2[j * a2_dim1 + 
					1], &c__1, &tau2, &a2[(j + 1) * 
					a2_dim1 + 1], lda2, &dwork[1], (
					ftnlen)4);
			    }
			    a2[j * a2_dim1 + 1] = tau2;
			} else {
			    alpha2 = a2[j * a2_dim1 + 1];
			    a2[j * a2_dim1 + 1] = 0.;
			}
			alpha = a1[j + j * a1_dim1];
			dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * 
				a1_dim1]);
			if (*k > j) {
			    i__4 = *k - j;
			    drot_(&i__4, &a1[j + (j + 1) * a1_dim1], lda1, &
				    a2[(j + 1) * a2_dim1 + 1], lda2, &c__, &s)
				    ;
			}
			cs[(*k << 1) + (j << 1) - 1] = c__;
			cs[(*k << 1) + (j << 1)] = s;
/* L160: */
		    }

/* L170: */
		}

	    } else {
		i__ = 1;
	    }

/*           Unblocked version for the last or only block. */

	    i__2 = *k;
	    for (j = i__; j <= i__2; ++j) {
		if (col2 > 1) {
		    dlarfg_(&col2, &a2[j * a2_dim1 + 1], &a2[j * a2_dim1 + 2],
			     &c__1, &tau2);
		    alpha2 = a2[j * a2_dim1 + 1];
		    if (*k > j) {
			a2[j * a2_dim1 + 1] = 1.;
			i__1 = *k - j;
			dlarf_("Left", &col2, &i__1, &a2[j * a2_dim1 + 1], &
				c__1, &tau2, &a2[(j + 1) * a2_dim1 + 1], lda2,
				 &dwork[1], (ftnlen)4);
		    }
		    a2[j * a2_dim1 + 1] = tau2;
		} else {
		    alpha2 = a2[j * a2_dim1 + 1];
		    a2[j * a2_dim1 + 1] = 0.;
		}
		alpha = a1[j + j * a1_dim1];
		dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * a1_dim1]);
		if (*k > j) {
		    i__1 = *k - j;
		    drot_(&i__1, &a1[j + (j + 1) * a1_dim1], lda1, &a2[(j + 1)
			     * a2_dim1 + 1], lda2, &c__, &s);
		}
		cs[(*k << 1) + (j << 1) - 1] = c__;
		cs[(*k << 1) + (j << 1)] = s;
/* L180: */
	    }

	    pst2 = *k * 5;
	} else {
	    pst2 = *k << 1;
	}

/*        Annihilate B with hyperbolic transformations. */

	nbl = min(*nb,*q);
	if (nbl > 0) {

/*           Blocked version. */

	    i__2 = *k - nbl + 1;
	    i__1 = nbl;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *k - i__ + 1;
		ib = min(i__3,nbl);
		dgeqr2_(q, &ib, &b[i__ * b_dim1 + 1], ldb, &cs[pst2 + i__], &
			dwork[1], &ierr);
		if (i__ + ib <= *k) {
		    dlarft_("Forward", "Columnwise", q, &ib, &b[i__ * b_dim1 
			    + 1], ldb, &cs[pst2 + i__], &dwork[1], k, (ftnlen)
			    7, (ftnlen)10);
		    i__3 = *k - i__ - ib + 1;
		    dlarfb_("Left", "Transpose", "Forward", "Columnwise", q, &
			    i__3, &ib, &b[i__ * b_dim1 + 1], ldb, &dwork[1], 
			    k, &b[(i__ + ib) * b_dim1 + 1], ldb, &dwork[ib + 
			    1], k, (ftnlen)4, (ftnlen)9, (ftnlen)7, (ftnlen)
			    10);
		}

/*              Annihilate the remaining parts of B. */

		i__3 = i__ + ib - 1;
		for (j = i__; j <= i__3; ++j) {
		    if (*q > 1) {
			i__4 = j - i__ + 1;
			dlarfg_(&i__4, &b[j * b_dim1 + 1], &b[j * b_dim1 + 2],
				 &c__1, &tau2);
			alpha2 = b[j * b_dim1 + 1];
			if (*k > j) {
			    b[j * b_dim1 + 1] = 1.;
			    i__4 = j - i__ + 1;
			    i__5 = *k - j;
			    dlarf_("Left", &i__4, &i__5, &b[j * b_dim1 + 1], &
				    c__1, &tau2, &b[(j + 1) * b_dim1 + 1], 
				    ldb, &dwork[1], (ftnlen)4);
			}
			b[j * b_dim1 + 1] = tau2;
		    } else {
			alpha2 = b[j * b_dim1 + 1];
			b[j * b_dim1 + 1] = 0.;
		    }

/*                 Create hyperbolic rotation. */

		    ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
		    if (ierr != 0) {

/*                    Error return:  The matrix is not positive definite. */

			*info = 1;
			return 0;
		    }

/*                 Apply hyperbolic rotation. */

		    if (*k > j) {
			i__4 = *k - j;
			d__1 = 1. / c__;
			dscal_(&i__4, &d__1, &a1[j + (j + 1) * a1_dim1], lda1)
				;
			i__4 = *k - j;
			d__1 = -s / c__;
			daxpy_(&i__4, &d__1, &b[(j + 1) * b_dim1 + 1], ldb, &
				a1[j + (j + 1) * a1_dim1], lda1);
			i__4 = *k - j;
			dscal_(&i__4, &c__, &b[(j + 1) * b_dim1 + 1], ldb);
			i__4 = *k - j;
			d__1 = -s;
			daxpy_(&i__4, &d__1, &a1[j + (j + 1) * a1_dim1], lda1,
				 &b[(j + 1) * b_dim1 + 1], ldb);
		    }
		    cs[(j << 1) - 1] = c__;
		    cs[j * 2] = s;
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
		dlarfg_(q, &b[j * b_dim1 + 1], &b[j * b_dim1 + 2], &c__1, &
			tau2);
		alpha2 = b[j * b_dim1 + 1];
		if (*k > j) {
		    b[j * b_dim1 + 1] = 1.;
		    i__2 = *k - j;
		    dlarf_("Left", q, &i__2, &b[j * b_dim1 + 1], &c__1, &tau2,
			     &b[(j + 1) * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)4);
		}
		b[j * b_dim1 + 1] = tau2;
	    } else if (*q > 0) {
		alpha2 = b[j * b_dim1 + 1];
		b[j * b_dim1 + 1] = 0.;
	    }
	    if (*q > 0) {

/*              Create hyperbolic rotation. */

		ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/*              Apply hyperbolic rotation. */

		if (*k > j) {
		    i__2 = *k - j;
		    d__1 = 1. / c__;
		    dscal_(&i__2, &d__1, &a1[j + (j + 1) * a1_dim1], lda1);
		    i__2 = *k - j;
		    d__1 = -s / c__;
		    daxpy_(&i__2, &d__1, &b[(j + 1) * b_dim1 + 1], ldb, &a1[j 
			    + (j + 1) * a1_dim1], lda1);
		    i__2 = *k - j;
		    dscal_(&i__2, &c__, &b[(j + 1) * b_dim1 + 1], ldb);
		    i__2 = *k - j;
		    d__1 = -s;
		    daxpy_(&i__2, &d__1, &a1[j + (j + 1) * a1_dim1], lda1, &b[
			    (j + 1) * b_dim1 + 1], ldb);
		}
		cs[(j << 1) - 1] = c__;
		cs[j * 2] = s;
	    }
/* L210: */
	}

    }

/* *** Last line of MB02CU *** */
    return 0;
} /* mb02cu_ */

