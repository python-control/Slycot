/* MB03ND.f -- translated by f2c (version 20100827).
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

integer mb03nd_(integer *n, doublereal *theta, doublereal *q2, doublereal *e2,
	 doublereal *pivmin, integer *info)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer j;
    static doublereal r__, t;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer numeig;


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

/*     To find the number of singular values of the bidiagonal matrix */

/*              |q(1) e(1)  .    ...    0   | */
/*              | 0   q(2) e(2)         .   | */
/*          J = | .                     .   | */
/*              | .                   e(N-1)| */
/*              | 0   ...     ...   0  q(N) | */

/*     which are less than or equal to a given bound THETA. */

/*     This routine is intended to be called only by other SLICOT */
/*     routines. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the bidiagonal matrix J.  N >= 0. */

/*     THETA   (input) DOUBLE PRECISION */
/*             Given bound. */
/*             Note: If THETA < 0.0 on entry, then MB03ND is set to 0 */
/*                   as the singular values of J are non-negative. */

/*     Q2      (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the squares of the diagonal */
/*             elements q(1),q(2),...,q(N) of the bidiagonal matrix J. */
/*             That is, Q2(i) = J(i,i)**2 for i = 1,2,...,N. */

/*     E2      (input) DOUBLE PRECISION array, dimension (N-1) */
/*             This array must contain the squares of the superdiagonal */
/*             elements e(1),e(2),...,e(N-1) of the bidiagonal matrix J. */
/*             That is, E2(k) = J(k,k+1)**2 for k = 1,2,...,N-1. */

/*     PIVMIN  (input) DOUBLE PRECISION */
/*             The minimum absolute value of a "pivot" in the Sturm */
/*             sequence loop. */
/*             PIVMIN >= max( max( |q(i)|, |e(k)| )**2*sf_min, sf_min ), */
/*             where i = 1,2,...,N, k = 1,2,...,N-1, and sf_min is at */
/*             least the smallest number that can divide one without */
/*             overflow (see LAPACK Library routine DLAMCH). */
/*             Note that this condition is not checked by the routine. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The computation of the number of singular values s(i) of J which */
/*     are less than or equal to THETA is based on applying Sylvester's */
/*     Law of Inertia, or equivalently, Sturm sequences [1,p.52] to the */
/*     unreduced symmetric tridiagonal matrices associated with J as */
/*     follows. Let T be the following 2N-by-2N symmetric matrix */
/*     associated with J: */

/*               | 0   J'| */
/*          T =  |       |. */
/*               | J   0 | */

/*     (The eigenvalues of T are given by s(1),s(2),...,s(N),-s(1),-s(2), */
/*     ...,-s(N)). Then, by permuting the rows and columns of T into the */
/*     order 1, N+1, 2, N+2, ..., N, 2N it follows that T is orthogonally */
/*     similar to the tridiagonal matrix T" with zeros on its diagonal */
/*     and q(1), e(1), q(2), e(2), ..., e(N-1), q(N) on its offdiagonals */
/*     [3,4]. If q(1),q(2),...,q(N) and e(1),e(2),...,e(N-1) are nonzero, */
/*     Sylvester's Law of Inertia may be applied directly to T". */
/*     Otherwise, T" is block diagonal and each diagonal block (which is */
/*     then unreduced) must be analysed separately by applying */
/*     Sylvester's Law of Inertia. */

/*     REFERENCES */

/*     [1] Parlett, B.N. */
/*         The Symmetric Eigenvalue Problem. */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1980. */

/*     [2] Demmel, J. and Kahan, W. */
/*         Computing Small Singular Values of Bidiagonal Matrices with */
/*         Guaranteed High Relative Accuracy. */
/*         Technical Report, Courant Inst., New York, March 1988. */

/*     [3] Van Huffel, S. and Vandewalle, J. */
/*         The Partial Total Least-Squares Algorithm. */
/*         J. Comput. and Appl. Math., 21, pp. 333-341, 1988. */

/*     [4] Golub, G.H. and Kahan, W. */
/*         Calculating the Singular Values and Pseudo-inverse of a */
/*         Matrix. */
/*         SIAM J. Numer. Anal., Ser. B, 2, pp. 205-224, 1965. */

/*     [5] Demmel, J.W., Dhillon, I. and Ren, H. */
/*         On the Correctness of Parallel Bisection in Floating Point. */
/*         Computer Science Division Technical Report UCB//CSD-94-805, */
/*         University of California, Berkeley, CA 94720, March 1994. */

/*     NUMERICAL ASPECTS */

/*     The singular values s(i) could also be obtained with the use of */
/*     the symmetric tridiagonal matrix T = J'J, whose eigenvalues are */
/*     the squared singular values of J [4,p.213]. However, the method */
/*     actually used by the routine is more accurate and equally */
/*     efficient (see [2]). */

/*     To avoid overflow, matrix J should be scaled so that its largest */
/*     element is no greater than  overflow**(1/2) * underflow**(1/4) */
/*     in absolute value (and not much smaller than that, for maximal */
/*     accuracy). */

/*     With respect to accuracy the following condition holds (see [2]): */

/*     If the established value is denoted by p, then at least p */
/*     singular values of J are less than or equal to */
/*     THETA/(1 - (3 x N - 1.5) x EPS) and no more than p singular values */
/*     are less than or equal to */
/*     THETA x (1 - (6 x N-2) x EPS)/(1 - (3 x N - 1.5) x EPS). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB03BD by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     July 10, 1997. */

/*     KEYWORDS */

/*     Bidiagonal matrix, singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments.  PIVMIN is not checked. */

    /* Parameter adjustments */
    --e2;
    --q2;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB03ND", &i__1, (ftnlen)6);
	return ret_val;
    }

/*     Quick return if possible. */

    if (*n == 0 || *theta < 0.) {
	ret_val = 0;
	return ret_val;
    }

    numeig = *n;
    t = -(*theta);
    r__ = t;
    if (abs(r__) < *pivmin) {
	r__ = -(*pivmin);
    }

    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	r__ = t - q2[j] / r__;
	if (abs(r__) < *pivmin) {
	    r__ = -(*pivmin);
	}
	if (r__ > 0.) {
	    --numeig;
	}
	r__ = t - e2[j] / r__;
	if (abs(r__) < *pivmin) {
	    r__ = -(*pivmin);
	}
	if (r__ > 0.) {
	    --numeig;
	}
/* L20: */
    }

    r__ = t - q2[*n] / r__;
    if (abs(r__) < *pivmin) {
	r__ = -(*pivmin);
    }
    if (r__ > 0.) {
	--numeig;
    }
    ret_val = numeig;

    return ret_val;
/* *** Last line of MB03ND *** */
} /* mb03nd_ */

