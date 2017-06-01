/* DE01PD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int de01pd_(char *conv, char *wght, integer *n, doublereal *
	a, doublereal *b, doublereal *w, integer *info, ftnlen conv_len, 
	ftnlen wght_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, l, m, p1, r1;
    static doublereal t1, t2, t3;
    static integer len;
    extern /* Subroutine */ int dg01od_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lconv, lwght;
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen);


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

/*     To compute the convolution or deconvolution of two real signals */
/*     A and B using the Hartley transform. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CONV    CHARACTER*1 */
/*             Indicates whether convolution or deconvolution is to be */
/*             performed as follows: */
/*             = 'C':  Convolution; */
/*             = 'D':  Deconvolution. */

/*     WGHT    CHARACTER*1 */
/*             Indicates whether the precomputed weights are available */
/*             or not, as follows: */
/*             = 'A':  available; */
/*             = 'N':  not available. */
/*             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is */
/*             set to 'A' on exit. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the first signal. */
/*             On exit, this array contains the convolution (if */
/*             CONV = 'C') or deconvolution (if CONV = 'D') of the two */
/*             signals. */

/*     B       (input) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the second signal. */
/*             NOTE that this array is overwritten. */

/*     W       (input/output) DOUBLE PRECISION array, */
/*                            dimension (N - LOG2(N)) */
/*             On entry with WGHT = 'A', this array must contain the long */
/*             weight vector computed by a previous call of this routine */
/*             or of the SLICOT Library routine DG01OD.f, with the same */
/*             value of N. If WGHT = 'N', the contents of this array on */
/*             entry is ignored. */
/*             On exit, this array contains the long weight vector. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine computes the convolution or deconvolution of two */
/*     real signals A and B using three scrambled Hartley transforms */
/*     (SLICOT Library routine DG01OD). */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N log(N)) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     KEYWORDS */

/*     Convolution, deconvolution, digital signal processing, */
/*     fast Hartley transform, real signals. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --w;
    --b;
    --a;

    /* Function Body */
    *info = 0;
    lconv = lsame_(conv, "C", (ftnlen)1, (ftnlen)1);
    lwght = lsame_(wght, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lconv && ! lsame_(conv, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lwght && ! lsame_(wght, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else {
	m = 0;
	j = 0;
	if (*n >= 1) {
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
L10:
	    if (j % 2 == 0) {
		j /= 2;
		++m;
		goto L10;
	    }
/*           END WHILE 10 */
	    if (j != 1) {
		*info = -3;
	    }
	} else if (*n < 0) {
	    *info = -3;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("DE01PD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n <= 0) {
	return 0;
    } else if (*n == 1) {
	if (lconv) {
	    a[1] *= b[1];
	} else {
	    a[1] /= b[1];
	}
	return 0;
    }

/*     Scrambled Hartley transforms of A and B. */

    dg01od_("OutputScrambled", wght, n, &a[1], &w[1], info, (ftnlen)15, (
	    ftnlen)1);
    dg01od_("OutputScrambled", wght, n, &b[1], &w[1], info, (ftnlen)15, (
	    ftnlen)1);

/*     Something similar to a Hadamard product/quotient. */

    len = 1;
    if (lconv) {
	a[1] = a[1] * 2. * b[1];
	a[2] = a[2] * 2. * b[2];

	i__1 = m - 1;
	for (l = 1; l <= i__1; ++l) {
	    len <<= 1;
	    r1 = len << 1;

	    i__2 = len + len / 2;
	    for (p1 = len + 1; p1 <= i__2; ++p1) {
		t1 = b[p1] + b[r1];
		t2 = b[p1] - b[r1];
		t3 = t2 * a[p1];
		a[p1] = t1 * a[p1] + t2 * a[r1];
		a[r1] = t1 * a[r1] - t3;
		--r1;
/* L20: */
	    }

/* L30: */
	}

    } else {

	a[1] = a[1] * .5 / b[1];
	a[2] = a[2] * .5 / b[2];

	i__1 = m - 1;
	for (l = 1; l <= i__1; ++l) {
	    len <<= 1;
	    r1 = len << 1;

	    i__2 = len + len / 2;
	    for (p1 = len + 1; p1 <= i__2; ++p1) {
		d__1 = b[p1] + b[r1];
		d__2 = b[r1] - b[p1];
		dladiv_(&a[p1], &a[r1], &d__1, &d__2, &t1, &t2);
		a[p1] = t1;
		a[r1] = t2;
		--r1;
/* L40: */
	    }

/* L50: */
	}

    }

/*     Transposed Hartley transform of A. */

    dg01od_("InputScrambled", wght, n, &a[1], &w[1], info, (ftnlen)14, (
	    ftnlen)1);
    if (lconv) {
	d__1 = .5 / (doublereal) (*n);
	dscal_(n, &d__1, &a[1], &c__1);
    } else {
	d__1 = 2. / (doublereal) (*n);
	dscal_(n, &d__1, &a[1], &c__1);
    }

    return 0;
/* *** Last line of DE01PD *** */
} /* de01pd_ */

