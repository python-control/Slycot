/* DG01OD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int dg01od_(char *scr, char *wght, integer *n, doublereal *a,
	 doublereal *w, integer *info, ftnlen scr_len, ftnlen wght_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, l, m, p1, p2, q1, q2, r1, r2, s1, s2;
    static doublereal t1, t2, cf, sf, th;
    static integer len;
    static logical lfwd, lscr;
    static integer wpos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lwght;
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

/*     To compute the (scrambled) discrete Hartley transform of */
/*     a real signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SCR     CHARACTER*1 */
/*             Indicates whether the signal is scrambled on input or */
/*             on output as follows: */
/*             = 'N':  the signal is not scrambled at all; */
/*             = 'I':  the input signal is bit-reversed; */
/*             = 'O':  the output transform is bit-reversed. */

/*     WGHT    CHARACTER*1 */
/*             Indicates whether the precomputed weights are available */
/*             or not, as follows: */
/*             = 'A':  available; */
/*             = 'N':  not available. */
/*             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is */
/*             set to 'A' on exit. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             Number of real samples. N must be a power of 2. */
/*             N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry with SCR = 'N' or SCR = 'O', this array must */
/*             contain the input signal. */
/*             On entry with SCR = 'I', this array must contain the */
/*             bit-reversed input signal. */
/*             On exit with SCR = 'N' or SCR = 'I', this array contains */
/*             the Hartley transform of the input signal. */
/*             On exit with SCR = 'O', this array contains the */
/*             bit-reversed Hartley transform. */

/*     W       (input/output) DOUBLE PRECISION array, */
/*                            dimension (N - LOG2(N)) */
/*             On entry with WGHT = 'A', this array must contain the long */
/*             weight vector computed by a previous call of this routine */
/*             with the same value of N. If WGHT = 'N', the contents of */
/*             this array on entry is ignored. */
/*             On exit, this array contains the long weight vector. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine uses a Hartley butterfly algorithm as described */
/*     in [1]. */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable and requires O(N log(N)) */
/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     KEYWORDS */

/*     Digital signal processing, fast Hartley transform, real signals. */

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
    --a;

    /* Function Body */
    *info = 0;
    lfwd = lsame_(scr, "N", (ftnlen)1, (ftnlen)1) || lsame_(scr, "I", (ftnlen)
	    1, (ftnlen)1);
    lscr = lsame_(scr, "I", (ftnlen)1, (ftnlen)1) || lsame_(scr, "O", (ftnlen)
	    1, (ftnlen)1);
    lwght = lsame_(wght, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! (lfwd || lscr)) {
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
	xerbla_("DG01OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n <= 1) {
	return 0;
    }

    if (! lwght) {

/*        Compute the long weight vector via subvector scaling. */

	r1 = 1;
	len = 1;
	th = atan(1.) * 4. / (doublereal) (*n);

	i__1 = m - 2;
	for (l = 1; l <= i__1; ++l) {
	    len <<= 1;
	    th *= 2.;
	    cf = cos(th);
	    sf = sin(th);
	    w[r1] = cf;
	    w[r1 + 1] = sf;
	    r1 += 2;

	    i__2 = len - 2;
	    for (i__ = 1; i__ <= i__2; i__ += 2) {
		w[r1] = cf * w[i__] - sf * w[i__ + 1];
		w[r1 + 1] = sf * w[i__] + cf * w[i__ + 1];
		r1 += 2;
/* L20: */
	    }

/* L30: */
	}

	p1 = 3;
	q1 = r1 - 2;

	for (l = m - 2; l >= 1; --l) {

	    i__1 = q1;
	    for (i__ = p1; i__ <= i__1; i__ += 4) {
		w[r1] = w[i__];
		w[r1 + 1] = w[i__ + 1];
		r1 += 2;
/* L40: */
	    }

	    p1 = q1 + 4;
	    q1 = r1 - 2;
/* L50: */
	}

	*(unsigned char *)wght = 'A';

    }

    if (lfwd && ! lscr) {

/*        Inplace shuffling of data. */

	j = 1;

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (j > i__) {
		t1 = a[i__];
		a[i__] = a[j];
		a[j] = t1;
	    }
	    l = *n / 2;
/*           REPEAT */
L60:
	    if (j > l) {
		j -= l;
		l /= 2;
		if (l >= 2) {
		    goto L60;
		}
	    }
/*           UNTIL ( L.LT.2 ) */
	    j += l;
/* L70: */
	}

    }

    if (lfwd) {

/*        Compute Hartley transform with butterfly operators. */

	i__1 = *n;
	for (j = 2; j <= i__1; j += 2) {
	    t1 = a[j];
	    a[j] = a[j - 1] - t1;
	    a[j - 1] += t1;
/* L110: */
	}

	len = 1;
	wpos = *n - (m << 1) + 1;

	i__1 = m - 1;
	for (l = 1; l <= i__1; ++l) {
	    len <<= 1;
	    p2 = 1;
	    q2 = len + 1;
	    r2 = len / 2 + 1;
	    s2 = r2 + q2 - 1;

	    i__2 = *n / (len << 1) - 1;
	    for (i__ = 0; i__ <= i__2; ++i__) {
		t1 = a[q2];
		a[q2] = a[p2] - t1;
		a[p2] += t1;
		t1 = a[s2];
		a[s2] = a[r2] - t1;
		a[r2] += t1;

		p1 = p2 + 1;
		q1 = p1 + len;
		r1 = q1 - 2;
		s1 = r1 + len;

		i__3 = wpos + len - 3;
		for (j = wpos; j <= i__3; j += 2) {
		    cf = w[j];
		    sf = w[j + 1];
		    t1 = cf * a[q1] + sf * a[s1];
		    t2 = -cf * a[s1] + sf * a[q1];
		    a[q1] = a[p1] - t1;
		    a[p1] += t1;
		    a[s1] = a[r1] - t2;
		    a[r1] += t2;
		    ++p1;
		    ++q1;
		    --r1;
		    --s1;
/* L120: */
		}

		p2 += len << 1;
		q2 += len << 1;
		r2 += len << 1;
		s2 += len << 1;
/* L130: */
	    }

	    wpos = wpos - (len << 1) + 2;
/* L140: */
	}

    } else {

/*        Compute Hartley transform with transposed butterfly operators. */

	wpos = 1;
	len = *n;

	for (l = m - 1; l >= 1; --l) {
	    len /= 2;
	    p2 = 1;
	    q2 = len + 1;
	    r2 = len / 2 + 1;
	    s2 = r2 + q2 - 1;

	    i__1 = *n / (len << 1) - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		t1 = a[q2];
		a[q2] = a[p2] - t1;
		a[p2] += t1;
		t1 = a[s2];
		a[s2] = a[r2] - t1;
		a[r2] += t1;

		p1 = p2 + 1;
		q1 = p1 + len;
		r1 = q1 - 2;
		s1 = r1 + len;

		i__2 = wpos + len - 3;
		for (j = wpos; j <= i__2; j += 2) {
		    cf = w[j];
		    sf = w[j + 1];
		    t1 = a[p1] - a[q1];
		    t2 = a[r1] - a[s1];
		    a[p1] += a[q1];
		    a[r1] += a[s1];
		    a[q1] = cf * t1 + sf * t2;
		    a[s1] = -cf * t2 + sf * t1;
		    ++p1;
		    ++q1;
		    --r1;
		    --s1;
/* L210: */
		}

		p2 += len << 1;
		q2 += len << 1;
		r2 += len << 1;
		s2 += len << 1;
/* L220: */
	    }

	    wpos = wpos + len - 2;
/* L230: */
	}

	i__1 = *n;
	for (j = 2; j <= i__1; j += 2) {
	    t1 = a[j];
	    a[j] = a[j - 1] - t1;
	    a[j - 1] += t1;
/* L240: */
	}

    }
    return 0;
/* *** Last line of DG01OD *** */
} /* dg01od_ */

