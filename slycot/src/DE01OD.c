/* DE01OD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int de01od_(char *conv, integer *n, doublereal *a, 
	doublereal *b, integer *info, ftnlen conv_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal ac, bc, ci, as;
    static integer kj;
    static doublereal bs, cr, ast;
    static integer nd2p1;
    extern /* Subroutine */ int dg01md_(char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dscal_(integer *, doublereal *, 
	    doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lconv;
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
/*     A and B. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CONV    CHARACTER*1 */
/*             Indicates whether convolution or deconvolution is to be */
/*             performed as follows: */
/*             = 'C':  Convolution; */
/*             = 'D':  Deconvolution. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2.  N >= 2. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the first signal. */
/*             On exit, this array contains the convolution (if */
/*             CONV = 'C') or deconvolution (if CONV = 'D') of the two */
/*             signals. */

/*     B       (input) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the second signal. */
/*             NOTE that this array is overwritten. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine computes the convolution or deconvolution of two real */
/*     signals A and B using an FFT algorithm (SLICOT Library routine */
/*     DG01MD). */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DE01CD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Convolution, deconvolution, digital signal processing, fast */
/*     Fourier transform, real signals. */

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
    --b;
    --a;

    /* Function Body */
    *info = 0;
    lconv = lsame_(conv, "C", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lconv && ! lsame_(conv, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else {
	j = 0;
	if (*n >= 2) {
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
L10:
	    if (j % 2 == 0) {
		j /= 2;
		goto L10;
	    }
/*           END WHILE 10 */
	}
	if (j != 1) {
	    *info = -2;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("DE01OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Fourier transform. */

    dg01md_("Direct", n, &a[1], &b[1], info, (ftnlen)6);

    if (lconv) {
	ast = a[1] * b[1];
    } else {
	if (b[1] == 0.) {
	    ast = 0.;
	} else {
	    ast = a[1] / b[1];
	}
    }

    nd2p1 = *n / 2 + 1;
    j = nd2p1;

    i__1 = *n;
    for (kj = nd2p1; kj <= i__1; ++kj) {

/*        Components of the transform of function A. */

	ac = (a[j] + a[kj]) * .5;
	as = (b[j] - b[kj]) * .5;

/*        Components of the transform of function B. */

	bc = (b[kj] + b[j]) * .5;
	bs = (a[kj] - a[j]) * .5;

/*        Deconvolution by complex division if CONV = 'D'; */
/*        Convolution by complex multiplication if CONV = 'C'. */

	if (lconv) {
	    cr = ac * bc - as * bs;
	    ci = as * bc + ac * bs;
	} else {
/* Computing MAX */
	    d__1 = abs(bc), d__2 = abs(bs);
	    if (max(d__1,d__2) == 0.) {
		cr = 0.;
		ci = 0.;
	    } else {
		dladiv_(&ac, &as, &bc, &bs, &cr, &ci);
	    }
	}

	a[j] = cr;
	b[j] = ci;
	a[kj] = cr;
	b[kj] = -ci;
	--j;
/* L20: */
    }
    a[1] = ast;
    b[1] = 0.;

/*     Inverse Fourier transform. */

    dg01md_("Inverse", n, &a[1], &b[1], info, (ftnlen)7);

    d__1 = 1. / (doublereal) (*n);
    dscal_(n, &d__1, &a[1], &c__1);

    return 0;
/* *** Last line of DE01OD *** */
} /* de01od_ */

