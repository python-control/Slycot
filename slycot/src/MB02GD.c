/* MB02GD.f -- translated by f2c (version 20100827).
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

static doublereal c_b12 = 1.;
static integer c__1 = 1;
static doublereal c_b24 = 0.;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02gd_(char *typet, char *triu, integer *k, integer *n, 
	integer *nl, integer *p, integer *s, doublereal *t, integer *ldt, 
	doublereal *rb, integer *ldrb, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen typet_len, ftnlen triu_len)
{
    /* System generated locals */
    integer rb_dim1, rb_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, j, nb, jj, kk, len;
    static doublereal dum[1];
    static integer pre, pdw, rnk, len2, head, lenr, ierr;
    static logical ltri;
    static integer ipvt[1], posr, sizr, stps;
    extern /* Subroutine */ int mb02cu_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), 
	    mb02cv_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpotrf_(
	    char *, integer *, doublereal *, integer *, integer *, ftnlen);
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

/*     To compute the Cholesky factor of a banded symmetric positive */
/*     definite (s.p.d.) block Toeplitz matrix, defined by either its */
/*     first block row, or its first block column, depending on the */
/*     routine parameter TYPET. */

/*     By subsequent calls of this routine the Cholesky factor can be */
/*     computed block column by block column. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; the Cholesky factor is upper */
/*                     triangular; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; the Cholesky factor is */
/*                     lower triangular. This choice results in a column */
/*                     oriented algorithm which is usually faster. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     TRIU    CHARACTER*1 */
/*             Specifies the structure of the last block in T, as */
/*             follows: */
/*             = 'N':  the last block has no special structure; */
/*             = 'T':  the last block is lower / upper triangular. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 1. */
/*             If TRIU = 'N',   N >= 1; */
/*             if TRIU = 'T',   N >= 2. */

/*     NL      (input)  INTEGER */
/*             The lower block bandwidth, i.e., NL + 1 is the number of */
/*             nonzero blocks in the first block column of the block */
/*             Toeplitz matrix. */
/*             If TRIU = 'N',   0 <= NL < N; */
/*             if TRIU = 'T',   1 <= NL < N. */

/*     P       (input)  INTEGER */
/*             The number of previously computed block rows / columns of */
/*             the Cholesky factor.  0 <= P <= N. */

/*     S       (input)  INTEGER */
/*             The number of block rows / columns of the Cholesky factor */
/*             to compute.  0 <= S <= N - P. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,(NL+1)*K) / (LDT,K) */
/*             On entry, if P = 0, the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array must contain the first */
/*             block row / column of an s.p.d. block Toeplitz matrix. */
/*             On entry, if P > 0, the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array must contain the P-th */
/*             block row / column of the Cholesky factor. */
/*             On exit, if INFO = 0, then the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array contains the (P+S)-th */
/*             block row / column of the Cholesky factor. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K) / MAX(1,(NL+1)*K). */

/*     RB      (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDRB,MIN(P+NL+S,N)*K) / (LDRB,MIN(P+S,N)*K) */
/*             On entry, if TYPET = 'R'  and  TRIU = 'N'  and  P > 0, */
/*             the leading (NL+1)*K-by-MIN(NL,N-P)*K part of this array */
/*             must contain the (P*K+1)-st to ((P+NL)*K)-th columns */
/*             of the upper Cholesky factor in banded format from a */
/*             previous call of this routine. */
/*             On entry, if TYPET = 'R'  and  TRIU = 'T'  and  P > 0, */
/*             the leading (NL*K+1)-by-MIN(NL,N-P)*K part of this array */
/*             must contain the (P*K+1)-st to (MIN(P+NL,N)*K)-th columns */
/*             of the upper Cholesky factor in banded format from a */
/*             previous call of this routine. */
/*             On exit, if TYPET = 'R'  and  TRIU = 'N', the leading */
/*             (NL+1)*K-by-MIN(NL+S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the */
/*             upper Cholesky factor in banded format. */
/*             On exit, if TYPET = 'R'  and  TRIU = 'T', the leading */
/*             (NL*K+1)-by-MIN(NL+S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the */
/*             upper Cholesky factor in banded format. */
/*             On exit, if TYPET = 'C'  and  TRIU = 'N', the leading */
/*             (NL+1)*K-by-MIN(S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower */
/*             Cholesky factor in banded format. */
/*             On exit, if TYPET = 'C'  and  TRIU = 'T', the leading */
/*             (NL*K+1)-by-MIN(S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower */
/*             Cholesky factor in banded format. */
/*             For further details regarding the band storage scheme see */
/*             the documentation of the LAPACK routine DPBTF2. */

/*     LDRB    INTEGER */
/*             The leading dimension of the array RB. */
/*             If TRIU = 'N',   LDRB >= MAX( (NL+1)*K,1 ); */
/*             if TRIU = 'T',   LDRB >= NL*K+1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -13,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first 1 + ( NL + 1 )*K*K elements of DWORK should be */
/*             preserved during successive calls of the routine. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1 + ( NL + 1 )*K*K + NL*K. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The Toeplitz matrix */
/*                   associated with T is not (numerically) positive */
/*                   definite. */

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

/*     The implemented method is numerically stable. */
/*                                3 */
/*     The algorithm requires O( K *N*NL ) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     Mar. 2004. */

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
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    rb_dim1 = *ldrb;
    rb_offset = 1 + rb_dim1;
    rb -= rb_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    ltri = lsame_(triu, "T", (ftnlen)1, (ftnlen)1);
    lenr = (*nl + 1) * *k;
    if (ltri) {
	sizr = *nl * *k + 1;
    } else {
	sizr = lenr;
    }
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);
    wrkmin = (lenr + *nl) * *k + 1;

/*     Check the scalar input parameters. */

    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (ltri || lsame_(triu, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*k < 0) {
	*info = -3;
    } else if (ltri && *n < 2 || ! ltri && *n < 1) {
	*info = -4;
    } else if (*nl >= *n || ltri && *nl < 1 || ! ltri && *nl < 0) {
	*info = -5;
    } else if (*p < 0 || *p > *n) {
	*info = -6;
    } else if (*s < 0 || *s > *n - *p) {
	*info = -7;
    } else if (isrow && *ldt < max(1,*k) || ! isrow && *ldt < max(1,lenr)) {
	*info = -9;
    } else if (ltri && *ldrb < sizr || ! ltri && *ldrb < max(1,lenr)) {
	*info = -11;
    } else if (*ldwork < wrkmin) {
	dwork[1] = (doublereal) wrkmin;
	*info = -13;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02GD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*s * *k == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Compute the generator if P = 0. */

    if (*p == 0) {
	if (isrow) {
	    dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }
	    if (*nl > 0) {
		i__1 = *nl * *k;
		dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &
			c_b12, &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], 
			ldt, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
	    }

/*           Copy the first block row to RB. */

	    if (ltri) {

		i__1 = lenr - *k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = min(i__,*k);
/* Computing MAX */
		    i__3 = sizr - i__ + 1;
		    dcopy_(&i__2, &t[i__ * t_dim1 + 1], &c__1, &rb[max(i__3,1)
			     + i__ * rb_dim1], &c__1);
/* L10: */
		}

		for (i__ = *k; i__ >= 1; --i__) {
		    dcopy_(&i__, &t[*k - i__ + 1 + (lenr - i__ + 1) * t_dim1],
			     &c__1, &rb[(lenr - i__ + 1) * rb_dim1 + 1], &
			    c__1);
/* L20: */
		}

	    } else {

		i__1 = lenr;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = min(i__,*k);
/* Computing MAX */
		    i__3 = sizr - i__ + 1;
		    dcopy_(&i__2, &t[i__ * t_dim1 + 1], &c__1, &rb[max(i__3,1)
			     + i__ * rb_dim1], &c__1);
/* L30: */
		}

	    }

/*           Quick return if N = 1. */

	    if (*n == 1) {
		dwork[1] = 1.;
		return 0;
	    }

	    i__1 = *nl * *k;
	    dlacpy_("All", k, &i__1, &t[(*k + 1) * t_dim1 + 1], ldt, &dwork[2]
		    , k, (ftnlen)3);
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[*nl * *k * *k + 2], k,
		     (ftnlen)3);
	    posr = *k + 1;
	} else {
	    dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }
	    if (*nl > 0) {
		i__1 = *nl * *k;
		dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &
			c_b12, &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (
			ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)7);
	    }

/*           Copy the first block column to RB. */

	    posr = 1;
	    if (ltri) {

		i__1 = *k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dcopy_(&sizr, &t[i__ + i__ * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
		    ++posr;
/* L40: */
		}

	    } else {

		i__1 = *k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = lenr - i__ + 1;
		    dcopy_(&i__2, &t[i__ + i__ * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
		    if (lenr < *n * *k && i__ > 1) {
			i__2 = i__ - 1;
			dlaset_("All", &i__2, &c__1, &c_b24, &c_b24, &rb[lenr 
				- i__ + 2 + posr * rb_dim1], ldrb, (ftnlen)3);
		    }
		    ++posr;
/* L50: */
		}

	    }

/*           Quick return if N = 1. */

	    if (*n == 1) {
		dwork[1] = 1.;
		return 0;
	    }

	    i__1 = *nl * *k;
	    dlacpy_("All", &i__1, k, &t[*k + 1 + t_dim1], ldt, &dwork[2], &
		    lenr, (ftnlen)3);
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[*nl * *k + 2], &lenr, 
		    (ftnlen)3);
	}
	pre = 1;
	stps = *s - 1;
    } else {
	pre = *p;
	stps = *s;
	posr = 1;
    }

    pdw = lenr * *k + 1;
    head = (pre - 1) * *k % lenr;

/*     Determine block size for the involved block Householder */
/*     transformations. */

    if (isrow) {
/* Computing MIN */
	i__1 = ilaenv_(&c__1, "DGEQRF", " ", k, &lenr, &c_n1, &c_n1, (ftnlen)
		6, (ftnlen)1);
	nb = min(i__1,*k);
    } else {
/* Computing MIN */
	i__1 = ilaenv_(&c__1, "DGELQF", " ", &lenr, k, &c_n1, &c_n1, (ftnlen)
		6, (ftnlen)1);
	nb = min(i__1,*k);
    }
    kk = pdw + (*k << 2);
    wrkopt = kk + lenr * nb;
    kk = *ldwork - kk;
    if (kk < lenr * nb) {
	nb = kk / lenr;
    }
    if (isrow) {
/* Computing MAX */
	i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", k, &lenr, &c_n1, &c_n1,
		 (ftnlen)6, (ftnlen)1);
	nbmin = max(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = 2, i__2 = ilaenv_(&c__2, "DGELQF", " ", &lenr, k, &c_n1, &c_n1,
		 (ftnlen)6, (ftnlen)1);
	nbmin = max(i__1,i__2);
    }
    if (nb < nbmin) {
	nb = 0;
    }

/*     Generator reduction process. */

    if (isrow) {

	i__1 = pre + stps - 1;
	for (i__ = pre; i__ <= i__1; ++i__) {
	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cu_("Row", k, k, k, &nb, &t[t_offset], ldt, dum, &c__1, &
		    dwork[head * *k + 2], k, &rnk, ipvt, &dwork[pdw + 1], &
		    c_b24, &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)
		    3);

	    if (ierr != 0) {

/*              Error return:  The positive definiteness is (numerically) */
/*                             not satisfied. */

		*info = 1;
		return 0;
	    }

/* Computing MAX */
/* Computing MIN */
	    i__3 = (*n - i__) * *k - *k, i__4 = lenr - head - *k;
	    i__2 = min(i__3,i__4);
	    len = max(i__2,0);
/* Computing MAX */
/* Computing MIN */
	    i__3 = (*n - i__) * *k - len - *k;
	    i__2 = min(i__3,head);
	    len2 = max(i__2,0);
	    if (len == lenr - *k) {
		*(unsigned char *)struct__ = *(unsigned char *)triu;
	    } else {
		*(unsigned char *)struct__ = 'N';
	    }
	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cv_("Row", struct__, k, &len, k, k, &nb, &c_n1, dum, &c__1, 
		    dum, &c__1, &dwork[head * *k + 2], k, &t[(*k + 1) * 
		    t_dim1 + 1], ldt, dum, &c__1, &dwork[(head + *k) * *k + 2]
		    , k, &dwork[pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, 
		    &ierr, (ftnlen)3, (ftnlen)1);

	    if ((*n - i__) * *k >= lenr) {
		*(unsigned char *)struct__ = *(unsigned char *)triu;
	    } else {
		*(unsigned char *)struct__ = 'N';
	    }
	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cv_("Row", struct__, k, &len2, k, k, &nb, &c_n1, dum, &c__1, 
		    dum, &c__1, &dwork[head * *k + 2], k, &t[(*k + len + 1) * 
		    t_dim1 + 1], ldt, dum, &c__1, &dwork[2], k, &dwork[pdw + 
		    1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)3, 
		    (ftnlen)1);

	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[head * *k + 2], k, (
		    ftnlen)3);

/*           Copy current block row to RB. */

	    if (ltri) {

/* Computing MIN */
		i__3 = len + len2 + *k, i__4 = lenr - *k;
		i__2 = min(i__3,i__4);
		for (j = 1; j <= i__2; ++j) {
		    i__3 = min(j,*k);
/* Computing MAX */
		    i__4 = sizr - j + 1;
		    dcopy_(&i__3, &t[j * t_dim1 + 1], &c__1, &rb[max(i__4,1) 
			    + (posr + j - 1) * rb_dim1], &c__1);
/* L60: */
		}

		if (len + len2 + *k >= lenr) {

		    for (jj = *k; jj >= 1; --jj) {
			dcopy_(&jj, &t[*k - jj + 1 + (lenr - jj + 1) * t_dim1]
				, &c__1, &rb[(posr + lenr - jj) * rb_dim1 + 1]
				, &c__1);
/* L70: */
		    }

		}
		posr += *k;

	    } else {

		i__2 = len + len2 + *k;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = min(j,*k);
/* Computing MAX */
		    i__4 = sizr - j + 1;
		    dcopy_(&i__3, &t[j * t_dim1 + 1], &c__1, &rb[max(i__4,1) 
			    + (posr + j - 1) * rb_dim1], &c__1);
		    if (j > lenr - *k) {
			i__3 = sizr - j;
			dlaset_("All", &i__3, &c__1, &c_b24, &c_b24, &rb[(
				posr + j - 1) * rb_dim1 + 1], &c__1, (ftnlen)
				3);
		    }
/* L80: */
		}

		posr += *k;
	    }
	    head = (head + *k) % lenr;
/* L90: */
	}

    } else {

	i__1 = pre + stps - 1;
	for (i__ = pre; i__ <= i__1; ++i__) {

	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cu_("Column", k, k, k, &nb, &t[t_offset], ldt, dum, &c__1, &
		    dwork[head + 2], &lenr, &rnk, ipvt, &dwork[pdw + 1], &
		    c_b24, &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)
		    6);

	    if (ierr != 0) {

/*              Error return:  The positive definiteness is (numerically) */
/*                             not satisfied. */

		*info = 1;
		return 0;
	    }

/* Computing MAX */
/* Computing MIN */
	    i__3 = (*n - i__) * *k - *k, i__4 = lenr - head - *k;
	    i__2 = min(i__3,i__4);
	    len = max(i__2,0);
/* Computing MAX */
/* Computing MIN */
	    i__3 = (*n - i__) * *k - len - *k;
	    i__2 = min(i__3,head);
	    len2 = max(i__2,0);
	    if (len == lenr - *k) {
		*(unsigned char *)struct__ = *(unsigned char *)triu;
	    } else {
		*(unsigned char *)struct__ = 'N';
	    }
	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cv_("Column", struct__, k, &len, k, k, &nb, &c_n1, dum, &c__1,
		     dum, &c__1, &dwork[head + 2], &lenr, &t[*k + 1 + t_dim1],
		     ldt, dum, &c__1, &dwork[head + *k + 2], &lenr, &dwork[
		    pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (
		    ftnlen)6, (ftnlen)1);

	    if ((*n - i__) * *k >= lenr) {
		*(unsigned char *)struct__ = *(unsigned char *)triu;
	    } else {
		*(unsigned char *)struct__ = 'N';
	    }
	    i__2 = *ldwork - pdw - (*k << 2);
	    mb02cv_("Column", struct__, k, &len2, k, k, &nb, &c_n1, dum, &
		    c__1, dum, &c__1, &dwork[head + 2], &lenr, &t[*k + len + 
		    1 + t_dim1], ldt, dum, &c__1, &dwork[2], &lenr, &dwork[
		    pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (
		    ftnlen)6, (ftnlen)1);

	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[head + 2], &lenr, (
		    ftnlen)3);

/*           Copy current block column to RB. */

	    if (ltri) {

		i__2 = *k;
		for (j = 1; j <= i__2; ++j) {
/* Computing MIN */
		    i__4 = sizr, i__5 = (*n - i__) * *k - j + 1;
		    i__3 = min(i__4,i__5);
		    dcopy_(&i__3, &t[j + j * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
		    ++posr;
/* L100: */
		}

	    } else {

		i__2 = *k;
		for (j = 1; j <= i__2; ++j) {
/* Computing MIN */
		    i__4 = sizr - j + 1, i__5 = (*n - i__) * *k - j + 1;
		    i__3 = min(i__4,i__5);
		    dcopy_(&i__3, &t[j + j * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
		    if (lenr < (*n - i__) * *k) {
			i__3 = j - 1;
/* Computing MIN */
			i__4 = sizr - j + 1, i__5 = (*n - i__) * *k - j + 1;
			dlaset_("All", &i__3, &c__1, &c_b24, &c_b24, &rb[min(
				i__4,i__5) + 1 + posr * rb_dim1], ldrb, (
				ftnlen)3);
		    }
		    ++posr;
/* L110: */
		}

	    }
	    head = (head + *k) % lenr;
/* L120: */
	}

    }
    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of MB02GD *** */
} /* mb02gd_ */

