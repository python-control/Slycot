/* TF01QD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tf01qd_(integer *nc, integer *nb, integer *n, integer *
	iord, doublereal *ar, doublereal *ma, doublereal *h__, integer *ldh, 
	integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, jj, jk, ki, nl;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nord, ldhnb;
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

/*     To compute N Markov parameters M(1), M(2),..., M(N) from a */
/*     multivariable system whose transfer function matrix G(z) is given. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NC      (input) INTEGER */
/*             The number of system outputs, i.e. the number of rows in */
/*             the transfer function matrix G(z).  NC >= 0. */

/*     NB      (input) INTEGER */
/*             The number of system inputs, i.e. the number of columns in */
/*             the transfer function matrix G(z).  NB >= 0. */

/*     N       (input) INTEGER */
/*             The number of Markov parameters M(k) to be computed. */
/*             N >= 0. */

/*     IORD    (input) INTEGER array, dimension (NC*NB) */
/*             This array must contain the order r of the elements of the */
/*             transfer function matrix G(z), stored row by row. */
/*             For example, the order of the (i,j)-th element of G(z) is */
/*             given by IORD((i-1)xNB+j). */

/*     AR      (input) DOUBLE PRECISION array, dimension (NA), where */
/*             NA = IORD(1) + IORD(2) + ... + IORD(NC*NB). */
/*             The leading NA elements of this array must contain the */
/*             denominator coefficients AR(1),...,AR(r) in equation (1) */
/*             of the (i,j)-th element of the transfer function matrix */
/*             G(z), stored row by row, i.e. in the order */
/*             (1,1),(1,2),...,(1,NB), (2,1),(2,2),...,(2,NB), ..., */
/*             (NC,1),(NC,2),...,(NC,NB). The coefficients must be given */
/*             in decreasing order of powers of z; the coefficient of the */
/*             highest order term is assumed to be equal to 1. */

/*     MA      (input) DOUBLE PRECISION array, dimension (NA) */
/*             The leading NA elements of this array must contain the */
/*             numerator coefficients MA(1),...,MA(r) in equation (1) */
/*             of the (i,j)-th element of the transfer function matrix */
/*             G(z), stored row by row, i.e. in the order */
/*             (1,1),(1,2),...,(1,NB), (2,1),(2,2),...,(2,NB), ..., */
/*             (NC,1),(NC,2),...,(NC,NB). The coefficients must be given */
/*             in decreasing order of powers of z. */

/*     H       (output) DOUBLE PRECISION array, dimension (LDH,N*NB) */
/*             The leading NC-by-N*NB part of this array contains the */
/*             multivariable Markov parameter sequence M(k), where each */
/*             parameter M(k) is an NC-by-NB matrix and k = 1,2,...,N. */
/*             The Markov parameters are stored such that H(i,(k-1)xNB+j) */
/*             contains the (i,j)-th element of M(k) for i = 1,2,...,NC */
/*             and j = 1,2,...,NB. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= MAX(1,NC). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The (i,j)-th element of G(z), defining the particular I/O transfer */
/*     between output i and input j, has the following form: */

/*                          -1         -2               -r */
/*                    MA(1)z   + MA(2)z   + ... + MA(r)z */
/*         G  (z) = ----------------------------------------.         (1) */
/*          ij                -1         -2               -r */
/*                  1 + AR(1)z   + AR(2)z   + ... + AR(r)z */

/*     The (i,j)-th element of G(z) is defined by its order r, its r */
/*     moving average coefficients (= numerator) MA(1),...,MA(r) and its */
/*     r autoregressive coefficients (= denominator) AR(1),...,AR(r). The */
/*     coefficient of the constant term in the denominator is assumed to */
/*     be equal to 1. */

/*     The relationship between the (i,j)-th element of the Markov */
/*     parameters M(1),M(2),...,M(N) and the corresponding element of the */
/*     transfer function matrix G(z) is given by: */

/*                               -1          -2                -k */
/*      G  (z) = M  (0) + M  (1)z   + M  (2)z   + ... + M  (k)z  + ...(2) */
/*       ij       ij       ij          ij                ij */

/*     Equating (1) and (2), we find that the relationship between the */
/*     (i,j)-th element of the Markov parameters M(k) and the ARMA */
/*     parameters AR(1),...,AR(r) and MA(1),...,MA(r) of the (i,j)-th */
/*     element of the transfer function matrix G(z) is as follows: */

/*        M  (1)   = MA(1), */
/*         ij */
/*                           k-1 */
/*        M  (k)   = MA(k) - SUM AR(p) x M  (k-p) for 1 < k <= r and */
/*         ij                p=1          ij */
/*                      r */
/*        M  (k+r) = - SUM AR(p) x M  (k+r-p) for k > 0. */
/*         ij          p=1          ij */

/*     From these expressions the Markov parameters M(k) are computed */
/*     element by element. */

/*     REFERENCES */

/*     [1] Luenberger, D.G. */
/*         Introduction to Dynamic Systems: Theory, Models and */
/*         Applications. */
/*         John Wiley & Sons, New York, 1979. */

/*     NUMERICAL ASPECTS */

/*     The computation of the (i,j)-th element of M(k) requires: */
/*        (k-1) multiplications and k additions if k <= r; */
/*          r   multiplications and r additions if k > r. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TF01ED by S. Van Huffel, Katholieke */
/*     Univ. Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Markov parameters, multivariable system, transfer function, */
/*     transfer matrix. */

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
    --iord;
    --ar;
    --ma;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*nc < 0) {
	*info = -1;
    } else if (*nb < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldh < max(1,*nc)) {
	*info = -8;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TF01QD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*nc,*nb);
    if (max(i__1,*n) == 0) {
	return 0;
    }

    ldhnb = *ldh * *nb;
    nl = 1;
    k = 1;

    i__1 = *nc;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *nb;
	for (j = 1; j <= i__2; ++j) {
	    nord = iord[k];
	    h__[i__ + j * h_dim1] = ma[nl];
	    jk = j;

	    i__3 = nord - 1;
	    for (ki = 1; ki <= i__3; ++ki) {
		jk += *nb;
		i__4 = -ldhnb;
		h__[i__ + jk * h_dim1] = ma[nl + ki] - ddot_(&ki, &ar[nl], &
			c__1, &h__[i__ + j * h_dim1], &i__4);
/* L20: */
	    }

	    i__3 = j + (*n - nord - 1) * *nb;
	    i__4 = *nb;
	    for (jj = j; i__4 < 0 ? jj >= i__3 : jj <= i__3; jj += i__4) {
		jk += *nb;
		i__5 = -ldhnb;
		h__[i__ + jk * h_dim1] = -ddot_(&nord, &ar[nl], &c__1, &h__[
			i__ + jj * h_dim1], &i__5);
/* L40: */
	    }

	    nl += nord;
	    ++k;
/* L50: */
	}

/* L60: */
    }

    return 0;
/* *** Last line of TF01QD *** */
} /* tf01qd_ */

