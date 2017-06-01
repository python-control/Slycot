/* DF01MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int df01md_(char *sico, integer *n, doublereal *dt, 
	doublereal *a, doublereal *dwork, integer *info, ftnlen sico_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, m;
    static doublereal a0;
    static integer i2;
    static doublereal w1, w2, w3;
    static integer md2, ind1, ind2;
    static logical lsig;
    extern /* Subroutine */ int dg01nd_(char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lsico;
    static doublereal pibym;
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

/*     To compute the sine transform or cosine transform of a real */
/*     signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SICO    CHARACTER*1 */
/*             Indicates whether the sine transform or cosine transform */
/*             is to be computed as follows: */
/*             = 'S':  The sine transform is computed; */
/*             = 'C':  The cosine transform is computed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2 plus 1. */
/*             N >= 5. */

/*     DT      (input) DOUBLE PRECISION */
/*             The sampling time of the signal. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the signal to be */
/*             processed. */
/*             On exit, this array contains either the sine transform, if */
/*             SICO = 'S', or the cosine transform, if SICO = 'C', of the */
/*             given signal. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N+1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let A(1), A(2),..., A(N) be a real signal of N samples. */

/*     If SICO = 'S', the routine computes the sine transform of A as */
/*     follows. First, transform A(i), i = 1,2,...,N, into the complex */
/*     signal B(i), i = 1,2,...,(N+1)/2, where */

/*        B(1) = -2*A(2), */
/*        B(i) = {A(2i-2) - A(2i)} - j*A(2i-1) for i = 2,3,...,(N-1)/2, */
/*        B((N+1)/2) = 2*A(N-1) and j**2 = -1. */

/*     Next, perform a discrete inverse Fourier transform on B(i) by */
/*     calling SLICOT Library Routine DG01ND, to give the complex signal */
/*     Z(i), i = 1,2,...,(N-1)/2, from which the real signal C(i) may be */
/*     obtained as follows: */

/*        C(2i-1) = Re(Z(i)),  C(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2. */

/*     Finally, compute the sine transform coefficients S ,S ,...,S */
/*                                                       1  2      N */
/*     given by */

/*        S  = 0, */
/*         1 */
/*                {                     [C(k) + C(N+1-k)]     } */
/*        S  = DT*{[C(k) - C(N+1-k)] - -----------------------}, */
/*         k      {                    [2*sin(pi*(k-1)/(N-1))]} */

/*           for k = 2,3,...,N-1, and */

/*        S = 0. */
/*         N */

/*     If SICO = 'C', the routine computes the cosine transform of A as */
/*     follows. First, transform A(i), i = 1,2,...,N, into the complex */
/*     signal B(i), i = 1,2,...,(N+1)/2, where */

/*        B(1) = 2*A(1), */
/*        B(i) = 2*A(2i-1) + 2*j*{[A(2i-2) - A(2i)]} */
/*        for i = 2,3,...,(N-1)/2 and B((N+1)/2) = 2*A(N). */

/*     Next, perform a discrete inverse Fourier transform on B(i) by */
/*     calling SLICOT Library Routine DG01ND, to give the complex signal */
/*     Z(i), i = 1,2,...,(N-1)/2, from which the real signal D(i) may be */
/*     obtained as follows: */

/*        D(2i-1) = Re(Z(i)),  D(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2. */

/*     Finally, compute the cosine transform coefficients S ,S ,...,S */
/*                                                         1  2      N */
/*     given by */

/*        S  = 2*DT*[D(1) + A0], */
/*         1 */
/*                {                     [D(k) - D(N+1-k)]     } */
/*        S  = DT*{[D(k) + D(N+1-k)] - -----------------------}, */
/*         k      {                    [2*sin(pi*(k-1)/(N-1))]} */


/*           for k = 2,3,...,N-1, and */

/*        S  = 2*DT*[D(1) - A0], */
/*         N */
/*                 (N-1)/2 */
/*     where A0 = 2*SUM   A(2i). */
/*                  i=1 */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     [2] Oppenheim, A.V. and Schafer, R.W. */
/*         Discrete-Time Signal Processing. */
/*         Prentice-Hall Signal Processing Series, 1989. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DF01AD by F. Dumortier, and */
/*     R.M.C. Dekeyser, State University of Gent, Belgium. */

/*     REVISIONS */

/*     V. Sima, Jan. 2003. */

/*     KEYWORDS */

/*     Digital signal processing, fast Fourier transform, complex */
/*     signals. */

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
    --dwork;
    --a;

    /* Function Body */
    *info = 0;
    lsico = lsame_(sico, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lsico && ! lsame_(sico, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else {
	m = 0;
	if (*n > 4) {
	    m = *n - 1;
/*           WHILE ( MOD( M, 2 ).EQ.0 ) DO */
L10:
	    if (m % 2 == 0) {
		m /= 2;
		goto L10;
	    }
/*           END WHILE 10 */
	}
	if (m != 1) {
	    *info = -2;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("DF01MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialisation. */

    m = *n - 1;
    md2 = (*n + 1) / 2;
    pibym = atan(1.) * 4. / (doublereal) m;
    i2 = 1;
    dwork[md2 + 1] = 0.;
    dwork[md2 * 2] = 0.;

    if (lsico) {

/*        Sine transform. */

	lsig = TRUE_;
	dwork[1] = a[2] * -2.;
	dwork[md2] = a[m] * 2.;

	i__1 = m;
	for (i__ = 4; i__ <= i__1; i__ += 2) {
	    ++i2;
	    dwork[i2] = a[i__ - 2] - a[i__];
	    dwork[md2 + i2] = -a[i__ - 1];
/* L20: */
	}

    } else {

/*        Cosine transform. */

	lsig = FALSE_;
	dwork[1] = a[1] * 2.;
	dwork[md2] = a[*n] * 2.;
	a0 = a[2];

	i__1 = m;
	for (i__ = 4; i__ <= i__1; i__ += 2) {
	    ++i2;
	    dwork[i2] = a[i__ - 1] * 2.;
	    dwork[md2 + i2] = (a[i__ - 2] - a[i__]) * 2.;
	    a0 += a[i__];
/* L30: */
	}

	a0 *= 2.;
    }

/*     Inverse Fourier transform. */

    i__1 = md2 - 1;
    dg01nd_("Inverse", &i__1, &dwork[1], &dwork[md2 + 1], info, (ftnlen)7);

/*     Sine or cosine coefficients. */

    if (lsico) {
	a[1] = 0.;
	a[*n] = 0.;
    } else {
	a[1] = *dt * 2. * (dwork[1] + a0);
	a[*n] = *dt * 2. * (dwork[1] - a0);
    }

    ind1 = md2 + 1;
    ind2 = *n;

    i__1 = m - 1;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	w1 = dwork[ind1];
	w2 = dwork[ind2];
	if (lsig) {
	    w2 = -w2;
	}
	w3 = sin(pibym * (doublereal) i__) * 2.;
	a[i__ + 1] = *dt * (w1 + w2 - (w1 - w2) / w3);
	++ind1;
	--ind2;
/* L40: */
    }

    ind1 = 2;
    ind2 = md2 - 1;

    i__1 = m - 2;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
	w1 = dwork[ind1];
	w2 = dwork[ind2];
	if (lsig) {
	    w2 = -w2;
	}
	w3 = sin(pibym * (doublereal) i__) * 2.;
	a[i__ + 1] = *dt * (w1 + w2 - (w1 - w2) / w3);
	++ind1;
	--ind2;
/* L50: */
    }

    return 0;
/* *** Last line of DF01MD *** */
} /* df01md_ */

