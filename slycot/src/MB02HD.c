/* MB02HD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static doublereal c_b15 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02hd_(char *triu, integer *k, integer *l, integer *m, 
	integer *ml, integer *n, integer *nu, integer *p, integer *s, 
	doublereal *tc, integer *ldtc, doublereal *tr, integer *ldtr, 
	doublereal *rb, integer *ldrb, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen triu_len)
{
    /* System generated locals */
    integer rb_dim1, rb_offset, tc_dim1, tc_offset, tr_dim1, tr_offset, i__1, 
	    i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, x, nb, kk, pt, pdc, len, pdr, pre, pfr, pdw, rnk, 
	    pnr, col2, len2, head, lenc, lenl, lenr, ierr;
    static logical ltri;
    static integer ipvt[1], posr, sizr, stps;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
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
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin;
    static char struct__[1];
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

/*     To compute, for a banded K*M-by-L*N block Toeplitz matrix T with */
/*     block size (K,L), specified by the nonzero blocks of its first */
/*     block column TC and row TR, a LOWER triangular matrix R (in band */
/*     storage scheme) such that */
/*                          T          T */
/*                         T  T  =  R R .                             (1) */

/*     It is assumed that the first MIN(M*K, N*L) columns of T are */
/*     linearly independent. */

/*     By subsequent calls of this routine, the matrix R can be computed */
/*     block column by block column. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRIU    CHARACTER*1 */
/*             Specifies the structure, if any, of the last blocks in TC */
/*             and TR, as follows: */
/*             = 'N':  TC and TR have no special structure; */
/*             = 'T':  TC and TR are upper and lower triangular, */
/*                     respectively. Depending on the block sizes, two */
/*                     different shapes of the last blocks in TC and TR */
/*                     are possible, as illustrated below: */

/*                     1)    TC       TR     2)   TC         TR */

/*                          x x x    x 0 0      x x x x    x 0 0 0 */
/*                          0 x x    x x 0      0 x x x    x x 0 0 */
/*                          0 0 x    x x x      0 0 x x    x x x 0 */
/*                          0 0 0    x x x */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 1. */

/*     ML      (input) INTEGER */
/*             The lower block bandwidth, i.e., ML + 1 is the number of */
/*             nonzero blocks in the first block column of T. */
/*             0 <= ML < M and (ML + 1)*K >= L and */
/*             if ( M*K <= N*L ),  ML >= M - INT( ( M*K - 1 )/L ) - 1; */
/*                                 ML >= M - INT( M*K/L ) or */
/*                                 MOD( M*K, L ) >= K; */
/*             if ( M*K >= N*L ),  ML*K >= N*( L - K ). */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T. */
/*             N >= 1. */

/*     NU      (input) INTEGER */
/*             The upper block bandwidth, i.e., NU + 1 is the number of */
/*             nonzero blocks in the first block row of T. */
/*             If TRIU = 'N',   0 <= NU < N and */
/*                              (M + NU)*L >= MIN( M*K, N*L ); */
/*             if TRIU = 'T',   MAX(1-ML,0) <= NU < N and */
/*                              (M + NU)*L >= MIN( M*K, N*L ). */

/*     P       (input)  INTEGER */
/*             The number of previously computed block columns of R. */
/*             P*L < MIN( M*K,N*L ) + L and P >= 0. */

/*     S       (input)  INTEGER */
/*             The number of block columns of R to compute. */
/*             (P+S)*L < MIN( M*K,N*L ) + L and S >= 0. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry, if P = 0, the leading (ML+1)*K-by-L part of this */
/*             array must contain the nonzero blocks in the first block */
/*             column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,(ML+1)*K),  if P = 0. */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,NU*L) */
/*             On entry, if P = 0, the leading K-by-NU*L part of this */
/*             array must contain the 2nd to the (NU+1)-st blocks of */
/*             the first block row of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR. */
/*             LDTR >= MAX(1,K),  if P = 0. */

/*     RB      (output)  DOUBLE PRECISION array, dimension */
/*             (LDRB,MIN( S*L,MIN( M*K,N*L )-P*L )) */
/*             On exit, if INFO = 0 and TRIU = 'N', the leading */
/*             MIN( ML+NU+1,N )*L-by-MIN( S*L,MIN( M*K,N*L )-P*L ) part */
/*             of this array contains the (P+1)-th to (P+S)-th block */
/*             column of the lower R factor (1) in band storage format. */
/*             On exit, if INFO = 0 and TRIU = 'T', the leading */
/*             MIN( (ML+NU)*L+1,N*L )-by-MIN( S*L,MIN( M*K,N*L )-P*L ) */
/*             part of this array contains the (P+1)-th to (P+S)-th block */
/*             column of the lower R factor (1) in band storage format. */
/*             For further details regarding the band storage scheme see */
/*             the documentation of the LAPACK routine DPBTF2. */

/*     LDRB    INTEGER */
/*             The leading dimension of the array RB. */
/*             LDRB >= MAX( MIN( ML+NU+1,N )*L,1 ),      if TRIU = 'N'; */
/*             LDRB >= MAX( MIN( (ML+NU)*L+1,N*L ),1 ),  if TRIU = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first 1 + 2*MIN( ML+NU+1,N )*L*(K+L) elements of DWORK */
/*             should be preserved during successive calls of the routine. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Let x = MIN( ML+NU+1,N ), then */
/*             LDWORK >= 1 + MAX( x*L*L + (2*NU+1)*L*K, */
/*                                2*x*L*(K+L) + (6+x)*L ),  if P = 0; */
/*             LDWORK >= 1 + 2*x*L*(K+L) + (6+x)*L,         if P > 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the full rank condition for the first MIN(M*K, N*L) */
/*                   columns of T is (numerically) violated. */

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

/*     The implemented method yields a factor R which has comparable */
/*     accuracy with the Cholesky factor of T^T * T. */
/*     The algorithm requires */
/*               2                                  2 */
/*           O( L *K*N*( ML + NU ) + N*( ML + NU )*L *( L + K ) ) */

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
    rb_dim1 = *ldrb;
    rb_offset = 1 + rb_dim1;
    rb -= rb_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    ltri = lsame_(triu, "T", (ftnlen)1, (ftnlen)1);
/* Computing MIN */
    i__1 = *ml + *nu + 1;
    x = min(i__1,*n);
    lenr = x * *l;
    if (ltri) {
/* Computing MIN */
	i__1 = (*ml + *nu) * *l + 1, i__2 = *n * *l;
	sizr = min(i__1,i__2);
    } else {
	sizr = lenr;
    }
    if (*p == 0) {
/* Computing MAX */
	i__1 = lenr * *l + ((*nu << 1) + 1) * *l * *k, i__2 = (lenr << 1) * (*
		k + *l) + (x + 6) * *l;
	wrkmin = max(i__1,i__2) + 1;
    } else {
	wrkmin = (lenr << 1) * (*k + *l) + 1 + (x + 6) * *l;
    }
    posr = 1;

/*     Check the scalar input parameters. */

    if (! (ltri || lsame_(triu, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*k < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*m < 1) {
	*info = -4;
    } else if (*ml >= *m || (*ml + 1) * *k < *l || *m * *k <= *n * *l && (*ml 
	    < *m - (*m * *k - 1) / *l - 1 || *ml < *m - *m * *k / *l && *m * *
	    k % *l < *k) || *m * *k >= *n * *l && *ml * *k < *n * (*l - *k)) {
	*info = -5;
    } else if (*n < 1) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MIN */
	i__1 = *m * *k, i__2 = *n * *l;
	if (*nu >= *n || *nu < 0 || ltri && *nu < 1 - *ml || (*m + *nu) * *l <
		 min(i__1,i__2)) {
	    *info = -7;
	} else /* if(complicated condition) */ {
/* Computing MIN */
	    i__1 = *m * *k, i__2 = *n * *l;
	    if (*p < 0 || *p * *l - *l >= min(i__1,i__2)) {
		*info = -8;
	    } else /* if(complicated condition) */ {
/* Computing MIN */
		i__1 = *m * *k, i__2 = *n * *l;
		if (*s < 0 || (*p + *s - 1) * *l >= min(i__1,i__2)) {
		    *info = -9;
		} else /* if(complicated condition) */ {
/* Computing MAX */
		    i__1 = 1, i__2 = (*ml + 1) * *k;
		    if (*p == 0 && *ldtc < max(i__1,i__2)) {
			*info = -11;
		    } else if (*p == 0 && *ldtr < max(1,*k)) {
			*info = -13;
		    } else if (*ldrb < max(sizr,1)) {
			*info = 15;
		    } else if (*ldwork < wrkmin) {
			dwork[1] = (doublereal) wrkmin;
			*info = -17;
		    }
		}
	    }
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02HD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*l * *k * *s == 0) {
	dwork[1] = 1.;
	return 0;
    }

    wrkopt = 1;

/*     Compute the generator if P = 0. */

    if (*p == 0) {

/*        1st column of the generator. */

	lenc = (*ml + 1) * *k;
/* Computing MAX */
/* Computing MIN */
	i__2 = *nu, i__3 = *n - *m;
	i__1 = *ml + 1 + min(i__2,i__3);
	lenl = max(i__1,0);
	pdc = lenr * *l + 1;
	pdw = pdc + lenc * *l;

/*        QR decomposition of the nonzero blocks in TC. */

	dlacpy_("All", &lenc, l, &tc[tc_offset], ldtc, &dwork[pdc + 1], &lenc,
		 (ftnlen)3);
	i__1 = *ldwork - pdw - *l;
	dgeqrf_(&lenc, l, &dwork[pdc + 1], &lenc, &dwork[pdw + 1], &dwork[pdw 
		+ *l + 1], &i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l + 1] + pdw + *l;
	wrkopt = max(i__1,i__2);

/*        The R factor is the transposed of the first block in the */
/*        generator. */

	ma02ad_("Upper part", l, l, &dwork[pdc + 1], &lenc, &dwork[2], &lenr, 
		(ftnlen)10);

/*        Get the first block column of the Q factor. */

	i__1 = *ldwork - pdw - *l;
	dorgqr_(&lenc, l, l, &dwork[pdc + 1], &lenc, &dwork[pdw + 1], &dwork[
		pdw + *l + 1], &i__1, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l + 1] + pdw + *l;
	wrkopt = max(i__1,i__2);

/*        Construct a flipped copy of TC for faster multiplication. */

	pt = lenc - (*k << 1) + 1;

	i__1 = pdw + *ml * *k * *l;
	i__2 = *k * *l;
	for (i__ = pdw + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
		 {
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
	    pt -= *k;
/* L10: */
	}

/*        Multiply T^T with the first block column of Q. */

	pdw = i__;
	pdr = *l + 2;
	len = *nu * *l;
	i__2 = lenr - *l;
	dlaset_("All", &i__2, l, &c_b10, &c_b10, &dwork[pdr], &lenr, (ftnlen)
		3);

	i__2 = *ml + 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
	    i__3 = i__ - 1, i__4 = *n - 1;
	    i__1 = min(i__3,i__4) * *l;
	    dgemm_("Transpose", "NonTranspose", &i__1, l, k, &c_b15, &dwork[
		    pdw], k, &dwork[pdc + 1], &lenc, &c_b15, &dwork[pdr], &
		    lenr, (ftnlen)9, (ftnlen)12);
	    if (len > 0) {
		dgemm_("Transpose", "NonTranspose", &len, l, k, &c_b15, &tr[
			tr_offset], ldtr, &dwork[pdc + 1], &lenc, &c_b15, &
			dwork[pdr + (i__ - 1) * *l], &lenr, (ftnlen)9, (
			ftnlen)12);
	    }
	    pdw -= *k * *l;
	    pdc += *k;
	    if (i__ >= *n - *nu) {
		len -= *l;
	    }
/* L20: */
	}

/*        Copy the first block column to R. */

	if (ltri) {

	    i__2 = *l;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
		i__3 = sizr, i__4 = *n * *l - i__ + 1;
		i__1 = min(i__3,i__4);
		dcopy_(&i__1, &dwork[(i__ - 1) * lenr + i__ + 1], &c__1, &rb[
			posr * rb_dim1 + 1], &c__1);
		++posr;
/* L30: */
	    }

	} else {

	    i__2 = *l;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = lenr - i__ + 1;
		dcopy_(&i__1, &dwork[(i__ - 1) * lenr + i__ + 1], &c__1, &rb[
			posr * rb_dim1 + 1], &c__1);
		if (lenr < *n * *l && i__ > 1) {
		    i__1 = i__ - 1;
		    dlaset_("All", &i__1, &c__1, &c_b10, &c_b10, &rb[lenr - 
			    i__ + 2 + posr * rb_dim1], ldrb, (ftnlen)3);
		}
		++posr;
/* L40: */
	    }

	}

/*        Quick return if N = 1. */

	if (*n == 1) {
	    dwork[1] = (doublereal) wrkopt;
	    return 0;
	}

/*        2nd column of the generator. */

	pdr = lenr * *l + 1;
	i__2 = *nu * *l;
	ma02ad_("All", k, &i__2, &tr[tr_offset], ldtr, &dwork[pdr + 1], &lenr,
		 (ftnlen)3);
	i__2 = lenr - *nu * *l;
	dlaset_("All", &i__2, k, &c_b10, &c_b10, &dwork[pdr + *nu * *l + 1], &
		lenr, (ftnlen)3);

/*        3rd column of the generator. */

	pnr = pdr + lenr * *k;
	i__2 = lenr - *l;
	dlacpy_("All", &i__2, l, &dwork[*l + 2], &lenr, &dwork[pnr + 1], &
		lenr, (ftnlen)3);
	dlaset_("All", l, l, &c_b10, &c_b10, &dwork[pnr + lenr - *l + 1], &
		lenr, (ftnlen)3);

/*        4th column of the generator. */

	pfr = pnr + lenr * *l;

	pdw = pfr + (*m - *ml - 1) * *l % lenr;
	pt = *ml * *k + 1;
/* Computing MIN */
	i__1 = *ml + 1;
	i__2 = min(i__1,lenl);
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw + 1], &
		    lenr, (ftnlen)3);
	    pt -= *k;
	    pdw = pfr + (pdw + *l - pfr) % lenr;
/* L50: */
	}
	pt = 1;
	i__2 = lenl;
	for (i__ = *ml + 2; i__ <= i__2; ++i__) {
	    ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw + 1],
		     &lenr, (ftnlen)3);
	    pt += *l;
	    pdw = pfr + (pdw + *l - pfr) % lenr;
/* L60: */
	}
	pre = 1;
	stps = *s - 1;
    } else {
	pdr = lenr * *l + 1;
	pnr = pdr + lenr * *k;
	pfr = pnr + lenr * *l;
	pre = *p;
	stps = *s;
    }

    pdw = pfr + lenr * *k;
    head = (pre - 1) * *l % lenr;

/*     Determine block size for the involved block Householder */
/*     transformations. */

/* Computing MIN */
    i__2 = ilaenv_(&c__1, "DGELQF", " ", &lenr, l, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
    nb = min(i__2,*l);
    kk = pdw + *l * 6;
/* Computing MAX */
    i__2 = wrkopt, i__1 = kk + lenr * nb;
    wrkopt = max(i__2,i__1);
    kk = *ldwork - kk;
    if (kk < lenr * nb) {
	nb = kk / lenr;
    }
/* Computing MAX */
    i__2 = 2, i__1 = ilaenv_(&c__2, "DGELQF", " ", &lenr, l, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
    nbmin = max(i__2,i__1);
    if (nb < nbmin) {
	nb = 0;
    }

/*     Generator reduction process. */

    i__2 = pre + stps - 1;
    for (i__ = pre; i__ <= i__2; ++i__) {

/*        The 4th generator column is not used in the first (M-ML) steps. */

	if (i__ < *m - *ml) {
	    col2 = *l;
	} else {
	    col2 = *k + *l;
	}

/* Computing MIN */
	i__1 = *l, i__3 = *m * *k - i__ * *l;
	kk = min(i__1,i__3);
	i__1 = kk + *k;
	i__3 = *ldwork - pdw - *l * 6;
	mb02cu_("Column", &kk, &i__1, &col2, &nb, &dwork[2], &lenr, &dwork[
		pdr + head + 1], &lenr, &dwork[pnr + head + 1], &lenr, &rnk, 
		ipvt, &dwork[pdw + 1], &c_b10, &dwork[pdw + *l * 6 + 1], &
		i__3, &ierr, (ftnlen)6);
	if (ierr != 0) {

/*           Error return:  The rank condition is (numerically) not */
/*                          satisfied. */

	    *info = 1;
	    return 0;
	}

/* Computing MAX */
/* Computing MIN */
	i__3 = (*n - i__) * *l - kk, i__4 = lenr - head - kk;
	i__1 = min(i__3,i__4);
	len = max(i__1,0);
/* Computing MAX */
/* Computing MIN */
	i__3 = (*n - i__) * *l - len - kk;
	i__1 = min(i__3,head);
	len2 = max(i__1,0);
	if (len == lenr - kk) {
	    *(unsigned char *)struct__ = *(unsigned char *)triu;
	} else {
	    *(unsigned char *)struct__ = 'N';
	}
	i__1 = kk + *k;
	i__3 = *ldwork - pdw - *l * 6;
	mb02cv_("Column", struct__, &kk, &len, &i__1, &col2, &nb, &c_n1, &
		dwork[2], &lenr, &dwork[pdr + head + 1], &lenr, &dwork[pnr + 
		head + 1], &lenr, &dwork[kk + 2], &lenr, &dwork[pdr + head + 
		kk + 1], &lenr, &dwork[pnr + head + kk + 1], &lenr, &dwork[
		pdw + 1], &dwork[pdw + *l * 6 + 1], &i__3, &ierr, (ftnlen)6, (
		ftnlen)1);

	if ((*n - i__) * *l >= lenr) {
	    *(unsigned char *)struct__ = *(unsigned char *)triu;
	} else {
	    *(unsigned char *)struct__ = 'N';
	}

	i__1 = kk + *k;
	i__3 = *ldwork - pdw - *l * 6;
	mb02cv_("Column", struct__, &kk, &len2, &i__1, &col2, &nb, &c_n1, &
		dwork[2], &lenr, &dwork[pdr + head + 1], &lenr, &dwork[pnr + 
		head + 1], &lenr, &dwork[kk + len + 2], &lenr, &dwork[pdr + 1]
		, &lenr, &dwork[pnr + 1], &lenr, &dwork[pdw + 1], &dwork[pdw 
		+ *l * 6 + 1], &i__3, &ierr, (ftnlen)6, (ftnlen)1);

	i__1 = *k + col2;
	dlaset_("All", l, &i__1, &c_b10, &c_b10, &dwork[pdr + head + 1], &
		lenr, (ftnlen)3);

/*        Copy current block column to R. */

	if (ltri) {

	    i__1 = kk;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__4 = sizr, i__5 = (*n - i__) * *l - j + 1;
		i__3 = min(i__4,i__5);
		dcopy_(&i__3, &dwork[(j - 1) * lenr + j + 1], &c__1, &rb[posr 
			* rb_dim1 + 1], &c__1);
		++posr;
/* L70: */
	    }

	} else {

	    i__1 = kk;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__4 = sizr - j + 1, i__5 = (*n - i__) * *l - j + 1;
		i__3 = min(i__4,i__5);
		dcopy_(&i__3, &dwork[(j - 1) * lenr + j + 1], &c__1, &rb[posr 
			* rb_dim1 + 1], &c__1);
		if (lenr < (*n - i__) * *l && j > 1) {
		    i__3 = j - 1;
/* Computing MIN */
		    i__4 = sizr - j + 1, i__5 = (*n - i__) * *l - j + 1;
		    dlaset_("All", &i__3, &c__1, &c_b10, &c_b10, &rb[min(i__4,
			    i__5) + 1 + posr * rb_dim1], ldrb, (ftnlen)3);
		}
		++posr;
/* L80: */
	    }

	}

	head = (head + *l) % lenr;
/* L90: */
    }

    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of MB02HD *** */
} /* mb02hd_ */

