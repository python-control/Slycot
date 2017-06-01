/* DG01ND.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int dg01nd_(char *indi, integer *n, doublereal *xr, 
	doublereal *xi, integer *info, ftnlen indi_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int dg01md_(char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static logical lindi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dg01ny_(char *, integer *, doublereal *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen);


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

/*     To compute the discrete Fourier transform, or inverse Fourier */
/*     transform, of a real signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     INDI    CHARACTER*1 */
/*             Indicates whether a Fourier transform or inverse Fourier */
/*             transform is to be performed as follows: */
/*             = 'D':  (Direct) Fourier transform; */
/*             = 'I':  Inverse Fourier transform. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             Half the number of real samples.  N must be a power of 2. */
/*             N >= 2. */

/*     XR      (input/output) DOUBLE PRECISION array, dimension (N+1) */
/*             On entry with INDI = 'D', the first N elements of this */
/*             array must contain the odd part of the input signal; for */
/*             example, XR(I) = A(2*I-1) for I = 1,2,...,N. */
/*             On entry with INDI = 'I', the first N+1 elements of this */
/*             array must contain the the real part of the input discrete */
/*             Fourier transform (computed, for instance, by a previous */
/*             call of the routine). */
/*             On exit with INDI = 'D', the first N+1 elements of this */
/*             array contain the real part of the output signal, that is */
/*             of the computed discrete Fourier transform. */
/*             On exit with INDI = 'I', the first N elements of this */
/*             array contain the odd part of the output signal, that is */
/*             of the computed inverse discrete Fourier transform. */

/*     XI      (input/output) DOUBLE PRECISION array, dimension (N+1) */
/*             On entry with INDI = 'D', the first N elements of this */
/*             array must contain the even part of the input signal; for */
/*             example, XI(I) = A(2*I) for I = 1,2,...,N. */
/*             On entry with INDI = 'I', the first N+1 elements of this */
/*             array must contain the the imaginary part of the input */
/*             discrete Fourier transform (computed, for instance, by a */
/*             previous call of the routine). */
/*             On exit with INDI = 'D', the first N+1 elements of this */
/*             array contain the imaginary part of the output signal, */
/*             that is of the computed discrete Fourier transform. */
/*             On exit with INDI = 'I', the first N elements of this */
/*             array contain the even part of the output signal, that is */
/*             of the computed inverse discrete Fourier transform. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let A(1),....,A(2*N) be a real signal of 2*N samples. Then the */
/*     first N+1 samples of the discrete Fourier transform of this signal */
/*     are given by the formula: */

/*                  2*N           ((m-1)*(i-1)) */
/*          FA(m) = SUM ( A(i) * W              ), */
/*                  i=1 */
/*                                                  2 */
/*     where m = 1,2,...,N+1, W = exp(-pi*j/N) and j = -1. */

/*     This transform can be computed as follows. First, transform A(i), */
/*     i = 1,2,...,2*N, into the complex signal Z(i) = (X(i),Y(i)), */
/*     i = 1,2,...,N. That is, X(i) = A(2*i-1) and Y(i) = A(2*i). Next, */
/*     perform a discrete Fourier transform on Z(i) by calling SLICOT */
/*     Library routine DG01MD. This gives a new complex signal FZ(k), */
/*     such that */

/*                   N            ((k-1)*(i-1)) */
/*          FZ(k) = SUM ( Z(i) * V              ), */
/*                  i=1 */

/*     where k = 1,2,...,N, V = exp(-2*pi*j/N).  Using the values of */
/*     FZ(k), the components of the discrete Fourier transform FA can be */
/*     computed by simple linear relations, implemented in the DG01NY */
/*     subroutine. */

/*     Finally, let */

/*          XR(k) = Re(FZ(k)), XI(k) = Im(FZ(k)),   k = 1,2,...,N, */

/*     be the contents of the arrays XR and XI on entry to DG01NY with */
/*     INDI = 'D', then on exit XR and XI contain the real and imaginary */
/*     parts of the Fourier transform of the original real signal A. */
/*     That is, */

/*          XR(m) = Re(FA(m)),  XI(m) = Im(FA(m)), */

/*     where m = 1,2,...,N+1. */

/*     If INDI = 'I', then the routine evaluates the inverse Fourier */
/*     transform of a complex signal which may itself be the discrete */
/*     Fourier transform of a real signal. */

/*     Let FA(m), m = 1,2,...,2*N, denote the full discrete Fourier */
/*     transform of a real signal A(i), i=1,2,...,2*N. The relationship */
/*     between FA and A is given by the formula: */

/*                 2*N            ((m-1)*(i-1)) */
/*          A(i) = SUM ( FA(m) * W              ), */
/*                 m=1 */

/*     where W = exp(pi*j/N). */

/*     Let */

/*          XR(m) = Re(FA(m)) and XI(m) = Im(FA(m)) for m = 1,2,...,N+1, */

/*     be the contents of the arrays XR and XI on entry to the routine */
/*     DG01NY with INDI = 'I', then on exit the first N samples of the */
/*     complex signal FZ are returned in XR and XI such that */

/*          XR(k) = Re(FZ(k)), XI(k) = Im(FZ(k)) and k = 1,2,...,N. */

/*     Next, an inverse Fourier transform is performed on FZ (e.g. by */
/*     calling SLICOT Library routine DG01MD), to give the complex signal */
/*     Z, whose i-th component is given by the formula: */

/*                  N             ((k-1)*(i-1)) */
/*          Z(i) = SUM ( FZ(k) * V              ), */
/*                 k=1 */

/*     where i = 1,2,...,N and V = exp(2*pi*j/N). */

/*     Finally, the 2*N samples of the real signal A can then be obtained */
/*     directly from Z. That is, */

/*          A(2*i-1) = Re(Z(i)) and A(2*i) = Im(Z(i)), for i = 1,2,...N. */

/*     Note that a discrete Fourier transform, followed by an inverse */
/*     transform will result in a signal which is a factor 2*N larger */
/*     than the original input signal. */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DG01BD by R. Dekeyser, and */
/*     F. Dumortier, State University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Complex signals, digital signal processing, fast Fourier */
/*     transform, real signals. */

/*     ****************************************************************** */

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
	xerbla_("DG01ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compute the Fourier transform of Z = (XR,XI). */

    if (! lindi) {
	dg01ny_(indi, n, &xr[1], &xi[1], (ftnlen)1);
    }

    dg01md_(indi, n, &xr[1], &xi[1], info, (ftnlen)1);

    if (lindi) {
	dg01ny_(indi, n, &xr[1], &xi[1], (ftnlen)1);
    }

    return 0;
/* *** Last line of DG01ND *** */
} /* dg01nd_ */

