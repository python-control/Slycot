/* DG01MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int dg01md_(char *indi, integer *n, doublereal *xr, 
	doublereal *xi, integer *info, ftnlen indi_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal ti, wi, tr, wr, pi2;
    static logical lindi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal whelp, wstpi, wstpr;
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

/*     To compute the discrete Fourier transform, or inverse transform, */
/*     of a complex signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     INDI    CHARACTER*1 */
/*             Indicates whether a Fourier transform or inverse Fourier */
/*             transform is to be performed as follows: */
/*             = 'D':  (Direct) Fourier transform; */
/*             = 'I':  Inverse Fourier transform. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of complex samples.  N must be a power of 2. */
/*             N >= 2. */

/*     XR      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the real part of either */
/*             the complex signal z if INDI = 'D', or f(z) if INDI = 'I'. */
/*             On exit, this array contains either the real part of the */
/*             computed Fourier transform f(z) if INDI = 'D', or the */
/*             inverse Fourier transform z of f(z) if INDI = 'I'. */

/*     XI      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the imaginary part of */
/*             either z if INDI = 'D', or f(z) if INDI = 'I'. */
/*             On exit, this array contains either the imaginary part of */
/*             f(z) if INDI = 'D', or z if INDI = 'I'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If INDI = 'D', then the routine performs a discrete Fourier */
/*     transform on the complex signal Z(i), i = 1,2,...,N. If the result */
/*     is denoted by FZ(k), k = 1,2,...,N, then the relationship between */
/*     Z and FZ is given by the formula: */

/*                     N            ((k-1)*(i-1)) */
/*            FZ(k) = SUM ( Z(i) * V              ), */
/*                    i=1 */
/*                                     2 */
/*     where V = exp( -2*pi*j/N ) and j  = -1. */

/*     If INDI = 'I', then the routine performs an inverse discrete */
/*     Fourier transform on the complex signal FZ(k), k = 1,2,...,N. If */
/*     the result is denoted by Z(i), i = 1,2,...,N, then the */
/*     relationship between Z and FZ is given by the formula: */

/*                    N             ((k-1)*(i-1)) */
/*            Z(i) = SUM ( FZ(k) * W              ), */
/*                   k=1 */

/*     where W = exp( 2*pi*j/N ). */

/*     Note that a discrete Fourier transform, followed by an inverse */
/*     discrete Fourier transform, will result in a signal which is a */
/*     factor N larger than the original input signal. */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DG01AD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Complex signals, digital signal processing, fast Fourier */
/*     transform. */

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
    --xi;
    --xr;

    /* Function Body */
    *info = 0;
    lindi = lsame_(indi, "D", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lindi && ! lsame_(indi, "I", (ftnlen)1, (ftnlen)1)) {
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
	xerbla_("DG01MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Inplace shuffling of data. */

    j = 1;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (j > i__) {
	    tr = xr[i__];
	    ti = xi[i__];
	    xr[i__] = xr[j];
	    xi[i__] = xi[j];
	    xr[j] = tr;
	    xi[j] = ti;
	}
	k = *n / 2;
/*        REPEAT */
L20:
	if (j > k) {
	    j -= k;
	    k /= 2;
	    if (k >= 2) {
		goto L20;
	    }
	}
/*        UNTIL ( K.LT.2 ) */
	j += k;
/* L30: */
    }

/*     Transform by decimation in time. */

    pi2 = atan(1.) * 8.;
    if (lindi) {
	pi2 = -pi2;
    }

    i__ = 1;

/*     WHILE ( I.LT.N ) DO */

L40:
    if (i__ < *n) {
	l = i__ << 1;
	whelp = pi2 / (doublereal) l;
	wstpi = sin(whelp);
	whelp = sin(whelp * .5);
	wstpr = whelp * -2. * whelp;
	wr = 1.;
	wi = 0.;

	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *n;
	    i__3 = l;
	    for (k = j; i__3 < 0 ? k >= i__2 : k <= i__2; k += i__3) {
		m = k + i__;
		tr = wr * xr[m] - wi * xi[m];
		ti = wr * xi[m] + wi * xr[m];
		xr[m] = xr[k] - tr;
		xi[m] = xi[k] - ti;
		xr[k] += tr;
		xi[k] += ti;
/* L50: */
	    }

	    whelp = wr;
	    wr = wr + wr * wstpr - wi * wstpi;
	    wi = wi + whelp * wstpi + wi * wstpr;
/* L60: */
	}

	i__ = l;
	goto L40;
/*        END WHILE 40 */
    }

    return 0;
/* *** Last line of DG01MD *** */
} /* dg01md_ */

