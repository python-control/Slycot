/* MB04NY.f -- translated by f2c (version 20100827).
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
static doublereal c_b15 = 1.;

/* Subroutine */ int mb04ny_(integer *m, integer *n, doublereal *v, integer *
	incv, doublereal *tau, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9;
    static integer iv;
    static doublereal sum;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *);


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

/*     To apply a real elementary reflector H to a real m-by-(n+1) */
/*     matrix C = [ A  B ], from the right, where A has one column. H is */
/*     represented in the form */
/*                                        ( 1 ) */
/*           H = I - tau * u *u',    u  = (   ), */
/*                                        ( v ) */
/*     where tau is a real scalar and v is a real n-vector. */

/*     If tau = 0, then H is taken to be the unit matrix. */

/*     In-line code is used if H has order < 11. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices A and B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix B.  N >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*ABS( INCV )) */
/*             The vector v in the representation of H. */

/*     INCV    (input) INTEGER */
/*             The increment between the elements of v.  INCV <> 0. */

/*     TAU     (input) DOUBLE PRECISION */
/*             The scalar factor of the elementary reflector H. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,1) */
/*             On entry, the leading M-by-1 part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading M-by-1 part of this array contains */
/*             the updated matrix A (the first column of C * H). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the updated matrix B (the last n columns of C * H). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (M) */
/*             DWORK is not referenced if H has order less than 11. */

/*     METHOD */

/*     The routine applies the elementary reflector H, taking the special */
/*     structure of C into account. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1998. */
/*     Based on LAPACK routines DLARFX and DLATZM. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, elementary reflector, orthogonal */
/*     transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --v;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dwork;

    /* Function Body */
    if (*tau == 0.) {
	return 0;
    }

/*     Form  C * H, where H has order n+1. */

    switch (*n + 1) {
	case 1:  goto L10;
	case 2:  goto L30;
	case 3:  goto L50;
	case 4:  goto L70;
	case 5:  goto L90;
	case 6:  goto L110;
	case 7:  goto L130;
	case 8:  goto L150;
	case 9:  goto L170;
	case 10:  goto L190;
    }

/*     Code for general N. Compute */

/*     w := C*u,  C := C - tau * w * u'. */

    dcopy_(m, &a[a_offset], &c__1, &dwork[1], &c__1);
    dgemv_("No transpose", m, n, &c_b15, &b[b_offset], ldb, &v[1], incv, &
	    c_b15, &dwork[1], &c__1, (ftnlen)12);
    d__1 = -(*tau);
    daxpy_(m, &d__1, &dwork[1], &c__1, &a[a_offset], &c__1);
    d__1 = -(*tau);
    dger_(m, n, &d__1, &dwork[1], &c__1, &v[1], incv, &b[b_offset], ldb);
    goto L210;
L10:

/*     Special code for 1 x 1 Householder */

    t1 = 1. - *tau;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	a[j + a_dim1] = t1 * a[j + a_dim1];
/* L20: */
    }
    goto L210;
L30:

/*     Special code for 2 x 2 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
/* L40: */
    }
    goto L210;
L50:

/*     Special code for 3 x 3 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
/* L60: */
    }
    goto L210;
L70:

/*     Special code for 4 x 4 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
/* L80: */
    }
    goto L210;
L90:

/*     Special code for 5 x 5 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
/* L100: */
    }
    goto L210;
L110:

/*     Special code for 6 x 6 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    iv += *incv;
    v5 = v[iv];
    t5 = *tau * v5;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
	b[j + b_dim1 * 5] -= sum * t5;
/* L120: */
    }
    goto L210;
L130:

/*     Special code for 7 x 7 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    iv += *incv;
    v5 = v[iv];
    t5 = *tau * v5;
    iv += *incv;
    v6 = v[iv];
    t6 = *tau * v6;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
	b[j + b_dim1 * 5] -= sum * t5;
	b[j + b_dim1 * 6] -= sum * t6;
/* L140: */
    }
    goto L210;
L150:

/*     Special code for 8 x 8 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    iv += *incv;
    v5 = v[iv];
    t5 = *tau * v5;
    iv += *incv;
    v6 = v[iv];
    t6 = *tau * v6;
    iv += *incv;
    v7 = v[iv];
    t7 = *tau * v7;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
	b[j + b_dim1 * 5] -= sum * t5;
	b[j + b_dim1 * 6] -= sum * t6;
	b[j + b_dim1 * 7] -= sum * t7;
/* L160: */
    }
    goto L210;
L170:

/*     Special code for 9 x 9 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    iv += *incv;
    v5 = v[iv];
    t5 = *tau * v5;
    iv += *incv;
    v6 = v[iv];
    t6 = *tau * v6;
    iv += *incv;
    v7 = v[iv];
    t7 = *tau * v7;
    iv += *incv;
    v8 = v[iv];
    t8 = *tau * v8;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7] + v8 * b[j + (b_dim1 << 3)];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
	b[j + b_dim1 * 5] -= sum * t5;
	b[j + b_dim1 * 6] -= sum * t6;
	b[j + b_dim1 * 7] -= sum * t7;
	b[j + (b_dim1 << 3)] -= sum * t8;
/* L180: */
    }
    goto L210;
L190:

/*     Special code for 10 x 10 Householder */

    iv = 1;
    if (*incv < 0) {
	iv = (-(*n) + 1) * *incv + 1;
    }
    v1 = v[iv];
    t1 = *tau * v1;
    iv += *incv;
    v2 = v[iv];
    t2 = *tau * v2;
    iv += *incv;
    v3 = v[iv];
    t3 = *tau * v3;
    iv += *incv;
    v4 = v[iv];
    t4 = *tau * v4;
    iv += *incv;
    v5 = v[iv];
    t5 = *tau * v5;
    iv += *incv;
    v6 = v[iv];
    t6 = *tau * v6;
    iv += *incv;
    v7 = v[iv];
    t7 = *tau * v7;
    iv += *incv;
    v8 = v[iv];
    t8 = *tau * v8;
    iv += *incv;
    v9 = v[iv];
    t9 = *tau * v9;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7] + v8 * b[j + (b_dim1 << 3)] + v9 * b[j + b_dim1 * 9];
	a[j + a_dim1] -= sum * *tau;
	b[j + b_dim1] -= sum * t1;
	b[j + (b_dim1 << 1)] -= sum * t2;
	b[j + b_dim1 * 3] -= sum * t3;
	b[j + (b_dim1 << 2)] -= sum * t4;
	b[j + b_dim1 * 5] -= sum * t5;
	b[j + b_dim1 * 6] -= sum * t6;
	b[j + b_dim1 * 7] -= sum * t7;
	b[j + (b_dim1 << 3)] -= sum * t8;
	b[j + b_dim1 * 9] -= sum * t9;
/* L200: */
    }
L210:
    return 0;
/* *** Last line of MB04NY *** */
} /* mb04ny_ */

