/* MB04PY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04py_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *dwork, 
	ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, sum;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;


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

/*     To apply a real elementary reflector H to a real m-by-n matrix */
/*     C, from either the left or the right. H is represented in the form */
/*                                        ( 1 ) */
/*           H = I - tau * u *u',    u  = (   ), */
/*                                        ( v ) */
/*     where tau is a real scalar and v is a real vector. */

/*     If tau = 0, then H is taken to be the unit matrix. */

/*     In-line code is used if H has order < 11. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Indicates whether the elementary reflector should be */
/*             applied from the left or from the right, as follows: */
/*             = 'L':  Compute H * C; */
/*             = 'R':  Compute C * H. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix C.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix C.  N >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*             (M-1), if SIDE = 'L', or */
/*             (N-1), if SIDE = 'R'. */
/*             The vector v in the representation of H. */

/*     TAU     (input) DOUBLE PRECISION */
/*             The scalar factor of the elementary reflector H. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix C. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix H * C, if SIDE = 'L', or C * H, if SIDE = 'R'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N), if SIDE = 'L', or */
/*                                               (M), if SIDE = 'R'. */
/*             DWORK is not referenced if H has order less than 11. */

/*     METHOD */

/*     The routine applies the elementary reflector H, taking its special */
/*     structure into account. The multiplications by the first component */
/*     of u (which is 1) are avoided, to increase the efficiency. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999. */
/*     This is a modification of LAPACK Library routine DLARFX. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, elementary reflector, orthogonal */
/*     transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    if (*tau == 0.) {
	return 0;
    }
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

	switch (*m) {
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

/*        Code for general M. */

/*        w := C'*u. */

	dcopy_(n, &c__[c_offset], ldc, &dwork[1], &c__1);
	i__1 = *m - 1;
	dgemv_("Transpose", &i__1, n, &c_b15, &c__[c_dim1 + 2], ldc, &v[1], &
		c__1, &c_b15, &dwork[1], &c__1, (ftnlen)9);

/*        C := C - tau * u * w'. */

	d__1 = -(*tau);
	daxpy_(n, &d__1, &dwork[1], &c__1, &c__[c_offset], ldc);
	i__1 = *m - 1;
	d__1 = -(*tau);
	dger_(&i__1, n, &d__1, &v[1], &c__1, &dwork[1], &c__1, &c__[c_dim1 + 
		2], ldc);
	goto L410;
L10:

/*        Special code for 1 x 1 Householder. */

	t1 = 1. - *tau;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
/* L20: */
	}
	goto L410;
L30:

/*        Special code for 2 x 2 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
/* L40: */
	}
	goto L410;
L50:

/*        Special code for 3 x 3 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
/* L60: */
	}
	goto L410;
L70:

/*        Special code for 4 x 4 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
/* L80: */
	}
	goto L410;
L90:

/*        Special code for 5 x 5 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
/* L100: */
	}
	goto L410;
L110:

/*        Special code for 6 x 6 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
	    c__[j * c_dim1 + 6] -= sum * t5;
/* L120: */
	}
	goto L410;
L130:

/*        Special code for 7 x 7 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
	    c__[j * c_dim1 + 6] -= sum * t5;
	    c__[j * c_dim1 + 7] -= sum * t6;
/* L140: */
	}
	goto L410;
L150:

/*        Special code for 8 x 8 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
	    c__[j * c_dim1 + 6] -= sum * t5;
	    c__[j * c_dim1 + 7] -= sum * t6;
	    c__[j * c_dim1 + 8] -= sum * t7;
/* L160: */
	}
	goto L410;
L170:

/*        Special code for 9 x 9 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8] + v8 * c__[j * 
		    c_dim1 + 9];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
	    c__[j * c_dim1 + 6] -= sum * t5;
	    c__[j * c_dim1 + 7] -= sum * t6;
	    c__[j * c_dim1 + 8] -= sum * t7;
	    c__[j * c_dim1 + 9] -= sum * t8;
/* L180: */
	}
	goto L410;
L190:

/*        Special code for 10 x 10 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8] + v8 * c__[j * 
		    c_dim1 + 9] + v9 * c__[j * c_dim1 + 10];
	    c__[j * c_dim1 + 1] -= sum * *tau;
	    c__[j * c_dim1 + 2] -= sum * t1;
	    c__[j * c_dim1 + 3] -= sum * t2;
	    c__[j * c_dim1 + 4] -= sum * t3;
	    c__[j * c_dim1 + 5] -= sum * t4;
	    c__[j * c_dim1 + 6] -= sum * t5;
	    c__[j * c_dim1 + 7] -= sum * t6;
	    c__[j * c_dim1 + 8] -= sum * t7;
	    c__[j * c_dim1 + 9] -= sum * t8;
	    c__[j * c_dim1 + 10] -= sum * t9;
/* L200: */
	}
	goto L410;
    } else {

/*        Form  C * H, where H has order n. */

	switch (*n) {
	    case 1:  goto L210;
	    case 2:  goto L230;
	    case 3:  goto L250;
	    case 4:  goto L270;
	    case 5:  goto L290;
	    case 6:  goto L310;
	    case 7:  goto L330;
	    case 8:  goto L350;
	    case 9:  goto L370;
	    case 10:  goto L390;
	}

/*        Code for general N. */

/*        w := C * u. */

	dcopy_(m, &c__[c_offset], &c__1, &dwork[1], &c__1);
	i__1 = *n - 1;
	dgemv_("No transpose", m, &i__1, &c_b15, &c__[(c_dim1 << 1) + 1], ldc,
		 &v[1], &c__1, &c_b15, &dwork[1], &c__1, (ftnlen)12);

/*        C := C - tau * w * u'. */

	d__1 = -(*tau);
	daxpy_(m, &d__1, &dwork[1], &c__1, &c__[c_offset], &c__1);
	i__1 = *n - 1;
	d__1 = -(*tau);
	dger_(m, &i__1, &d__1, &dwork[1], &c__1, &v[1], &c__1, &c__[(c_dim1 <<
		 1) + 1], ldc);
	goto L410;
L210:

/*        Special code for 1 x 1 Householder. */

	t1 = 1. - *tau;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    c__[j + c_dim1] = t1 * c__[j + c_dim1];
/* L220: */
	}
	goto L410;
L230:

/*        Special code for 2 x 2 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
/* L240: */
	}
	goto L410;
L250:

/*        Special code for 3 x 3 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
/* L260: */
	}
	goto L410;
L270:

/*        Special code for 4 x 4 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
/* L280: */
	}
	goto L410;
L290:

/*        Special code for 5 x 5 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
/* L300: */
	}
	goto L410;
L310:

/*        Special code for 6 x 6 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
	    c__[j + c_dim1 * 6] -= sum * t5;
/* L320: */
	}
	goto L410;
L330:

/*        Special code for 7 x 7 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
	    c__[j + c_dim1 * 6] -= sum * t5;
	    c__[j + c_dim1 * 7] -= sum * t6;
/* L340: */
	}
	goto L410;
L350:

/*        Special code for 8 x 8 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
	    c__[j + c_dim1 * 6] -= sum * t5;
	    c__[j + c_dim1 * 7] -= sum * t6;
	    c__[j + (c_dim1 << 3)] -= sum * t7;
/* L360: */
	}
	goto L410;
L370:

/*        Special code for 9 x 9 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)] + v8 * c__[j + 
		    c_dim1 * 9];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
	    c__[j + c_dim1 * 6] -= sum * t5;
	    c__[j + c_dim1 * 7] -= sum * t6;
	    c__[j + (c_dim1 << 3)] -= sum * t7;
	    c__[j + c_dim1 * 9] -= sum * t8;
/* L380: */
	}
	goto L410;
L390:

/*        Special code for 10 x 10 Householder. */

	v1 = v[1];
	t1 = *tau * v1;
	v2 = v[2];
	t2 = *tau * v2;
	v3 = v[3];
	t3 = *tau * v3;
	v4 = v[4];
	t4 = *tau * v4;
	v5 = v[5];
	t5 = *tau * v5;
	v6 = v[6];
	t6 = *tau * v6;
	v7 = v[7];
	t7 = *tau * v7;
	v8 = v[8];
	t8 = *tau * v8;
	v9 = v[9];
	t9 = *tau * v9;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)] + v8 * c__[j + 
		    c_dim1 * 9] + v9 * c__[j + c_dim1 * 10];
	    c__[j + c_dim1] -= sum * *tau;
	    c__[j + (c_dim1 << 1)] -= sum * t1;
	    c__[j + c_dim1 * 3] -= sum * t2;
	    c__[j + (c_dim1 << 2)] -= sum * t3;
	    c__[j + c_dim1 * 5] -= sum * t4;
	    c__[j + c_dim1 * 6] -= sum * t5;
	    c__[j + c_dim1 * 7] -= sum * t6;
	    c__[j + (c_dim1 << 3)] -= sum * t7;
	    c__[j + c_dim1 * 9] -= sum * t8;
	    c__[j + c_dim1 * 10] -= sum * t9;
/* L400: */
	}
	goto L410;
    }
L410:
    return 0;

/* *** Last line of MB04PY *** */
} /* mb04py_ */

