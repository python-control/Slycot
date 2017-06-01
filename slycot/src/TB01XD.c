/* TB01XD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb01xd_(char *jobd, integer *n, integer *m, integer *p, 
	integer *kl, integer *ku, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, integer *info, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, j1, nm1, lda1;
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer minmp, maxmp;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To apply a special transformation to a system given as a triple */
/*     (A,B,C), */

/*        A <-- P * A' * P,  B <-- P * C',  C <-- B' * P, */

/*     where P is a matrix with 1 on the secondary diagonal, and with 0 */
/*     in the other entries. Matrix A can be specified as a band matrix. */
/*     Optionally, matrix D of the system can be transposed. This */
/*     transformation is actually a special similarity transformation of */
/*     the dual system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of matrix B */
/*             and the number of columns of matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER. */
/*             The number of columns of matrix B. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER. */
/*             The number of rows of matrix C. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     KL      (input) INTEGER */
/*             The number of subdiagonals of A to be transformed. */
/*             MAX( 0, N-1 ) >= KL >= 0. */

/*     KU      (input) INTEGER */
/*             The number of superdiagonals of A to be transformed. */
/*             MAX( 0, N-1 ) >= KU >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed (pertransposed) matrix P*A'*P. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, the leading N-by-P part of this array contains */
/*             the dual input/state matrix P*C'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0 or  P > 0. */
/*             LDB >= 1        if M = 0 and P = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the dual state/output matrix B'*P. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1          if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDD,MAX(M,P)) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the original direct transmission */
/*             matrix D. */
/*             On exit, if JOBD = 'D', the leading M-by-P part of this */
/*             array contains the transposed direct transmission matrix */
/*             D'. The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,M,P) if JOBD = 'D'. */
/*             LDD >= 1          if JOBD = 'Z'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The rows and/or columns of the matrices of the triplet (A,B,C) */
/*     and, optionally, of the matrix D are swapped in a special way. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */
/*     Partly based on routine DMPTR (A. Varga, German Aerospace */
/*     Research Establishment, DLR, Aug. 1992). */


/*     REVISIONS */

/*     07-31-1998, 04-25-1999, A. Varga. */
/*     03-16-2004, V. Sima. */

/*     KEYWORDS */

/*     Matrix algebra, matrix operations, similarity transformation. */

/*  ********************************************************************* */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

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

    /* Function Body */
    *info = 0;
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
    maxmp = max(*m,*p);
    minmp = min(*m,*p);
    nm1 = *n - 1;

    if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*kl < 0 || *kl > max(0,nm1)) {
	*info = -5;
    } else if (*ku < 0 || *ku > max(0,nm1)) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (maxmp > 0 && *ldb < max(1,*n) || minmp == 0 && *ldb < 1) {
	*info = -10;
    } else if (*ldc < 1 || *n > 0 && *ldc < maxmp) {
	*info = -12;
    } else if (*ldd < 1 || ljobd && *ldd < maxmp) {
	*info = -14;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB01XD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (ljobd) {

/*        Replace D by D', if non-scalar. */

	i__1 = maxmp;
	for (j = 1; j <= i__1; ++j) {
	    if (j < minmp) {
		i__2 = minmp - j;
		dswap_(&i__2, &d__[j + 1 + j * d_dim1], &c__1, &d__[j + (j + 
			1) * d_dim1], ldd);
	    } else if (j > *p) {
		dcopy_(p, &d__[j * d_dim1 + 1], &c__1, &d__[j + d_dim1], ldd);
	    } else if (j > *m) {
		dcopy_(m, &d__[j + d_dim1], ldd, &d__[j * d_dim1 + 1], &c__1);
	    }
/* L5: */
	}

    }

    if (*n == 0) {
	return 0;
    }

/*     Replace matrix A by P*A'*P. */

    if (*kl == nm1 && *ku == nm1) {

/*        Full matrix A. */

	i__1 = nm1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n - j;
	    i__3 = -(*lda);
	    dswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[*n - j + 1 + (j + 1) *
		     a_dim1], &i__3);
/* L10: */
	}

    } else {

/*        Band matrix A. */

	lda1 = *lda + 1;

/*        Pertranspose the KL subdiagonals. */

/* Computing MIN */
	i__2 = *kl, i__3 = *n - 2;
	i__1 = min(i__2,i__3);
	for (j = 1; j <= i__1; ++j) {
	    j1 = (*n - j) / 2;
	    i__2 = -lda1;
	    dswap_(&j1, &a[j + 1 + a_dim1], &lda1, &a[*n - j1 + 1 + (*n - j1 
		    + 1 - j) * a_dim1], &i__2);
/* L20: */
	}

/*        Pertranspose the KU superdiagonals. */

/* Computing MIN */
	i__2 = *ku, i__3 = *n - 2;
	i__1 = min(i__2,i__3);
	for (j = 1; j <= i__1; ++j) {
	    j1 = (*n - j) / 2;
	    i__2 = -lda1;
	    dswap_(&j1, &a[(j + 1) * a_dim1 + 1], &lda1, &a[*n - j1 + 1 - j + 
		    (*n - j1 + 1) * a_dim1], &i__2);
/* L30: */
	}

/*        Pertranspose the diagonal. */

	j1 = *n / 2;
	i__1 = -lda1;
	dswap_(&j1, &a[a_dim1 + 1], &lda1, &a[*n - j1 + 1 + (*n - j1 + 1) * 
		a_dim1], &i__1);

    }

/*     Replace matrix B by P*C' and matrix C by B'*P. */

    i__1 = maxmp;
    for (j = 1; j <= i__1; ++j) {
	if (j <= minmp) {
	    i__2 = -(*ldc);
	    dswap_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], &i__2);
	} else if (j > *p) {
	    i__2 = -(*ldc);
	    dcopy_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], &i__2);
	} else {
	    i__2 = -(*ldc);
	    dcopy_(n, &c__[j + c_dim1], &i__2, &b[j * b_dim1 + 1], &c__1);
	}
/* L40: */
    }

    return 0;
/* *** Last line of TB01XD *** */
} /* tb01xd_ */

