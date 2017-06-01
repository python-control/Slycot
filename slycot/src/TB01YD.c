/* TB01YD.f -- translated by f2c (version 20100827).
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

static integer c_n1 = -1;
static integer c__1 = 1;

/* Subroutine */ int tb01yd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer j, nby2;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To apply a special similarity transformation to a system given as */
/*     a triple (A,B,C), */

/*        A <-- P * A * P,  B <-- P * B,  C <-- C * P, */

/*     where P is a matrix with 1 on the secondary diagonal, and with 0 */
/*     in the other entries. */

/*     ARGUMENTS */

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

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed matrix P*A*P. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed matrix P*B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0. */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*P. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The rows and/or columns of the matrices of the triplet (A,B,C) */
/*     are swapped in a special way. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */


/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

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

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < 1 || *m > 0 && *ldb < *n) {
	*info = -7;
    } else if (*ldc < max(1,*p)) {
	*info = -9;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB01YD", &i__1, (ftnlen)6);
	return 0;
    }

    if (*n <= 1) {
	return 0;
    }

/*     Transform the matrix A. */

    nby2 = *n / 2;

    i__1 = nby2;
    for (j = 1; j <= i__1; ++j) {
	dswap_(n, &a[j * a_dim1 + 1], &c_n1, &a[(*n - j + 1) * a_dim1 + 1], &
		c__1);
/* L10: */
    }

    if (*n % 2 != 0 && *n > 2) {
	dswap_(&nby2, &a[nby2 + 2 + (nby2 + 1) * a_dim1], &c_n1, &a[(nby2 + 1)
		 * a_dim1 + 1], &c__1);
    }

    if (*m > 0) {

/*        Transform the matrix B. */

	i__1 = nby2;
	for (j = 1; j <= i__1; ++j) {
	    dswap_(m, &b[j + b_dim1], ldb, &b[*n - j + 1 + b_dim1], ldb);
/* L20: */
	}

    }

    if (*p > 0) {

/*        Transform the matrix C. */

	i__1 = nby2;
	for (j = 1; j <= i__1; ++j) {
	    dswap_(p, &c__[j * c_dim1 + 1], &c__1, &c__[(*n - j + 1) * c_dim1 
		    + 1], &c__1);
/* L30: */
	}

    }

    return 0;
/* *** Last line of TB01YD *** */
} /* tb01yd_ */

