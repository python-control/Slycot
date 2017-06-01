/* MC01QD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mc01qd_(integer *da, integer *db, doublereal *a, 
	doublereal *b, doublereal *rq, integer *iwarn, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer n;
    static doublereal q;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), xerbla_(char *,
	     integer *, ftnlen);


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

/*     To compute, for two given real polynomials A(x) and B(x), the */
/*     quotient polynomial Q(x) and the remainder polynomial R(x) of */
/*     A(x) divided by B(x). */

/*     The polynomials Q(x) and R(x) satisfy the relationship */

/*        A(x) = B(x) * Q(x) + R(x), */

/*     where the degree of R(x) is less than the degree of B(x). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the numerator polynomial A(x).  DA >= -1. */

/*     DB      (input/output) INTEGER */
/*             On entry, the degree of the denominator polynomial B(x). */
/*             DB >= 0. */
/*             On exit, if B(DB+1) = 0.0 on entry, then DB contains the */
/*             index of the highest power of x for which B(DB+1) <> 0.0. */

/*     A       (input) DOUBLE PRECISION array, dimension (DA+1) */
/*             This array must contain the coefficients of the */
/*             numerator polynomial A(x) in increasing powers of x */
/*             unless DA = -1 on entry, in which case A(x) is taken */
/*             to be the zero polynomial. */

/*     B       (input) DOUBLE PRECISION array, dimension (DB+1) */
/*             This array must contain the coefficients of the */
/*             denominator polynomial B(x) in increasing powers of x. */

/*     RQ      (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             If DA < DB on exit, then this array contains the */
/*             coefficients of the remainder polynomial R(x) in */
/*             increasing powers of x; Q(x) is the zero polynomial. */
/*             Otherwise, the leading DB elements of this array contain */
/*             the coefficients of R(x) in increasing powers of x, and */
/*             the next (DA-DB+1) elements contain the coefficients of */
/*             Q(x) in increasing powers of x. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = k:  if the degree of the denominator polynomial B(x) has */
/*                   been reduced to (DB - k) because B(DB+1-j) = 0.0 on */
/*                   entry for j = 0, 1, ..., k-1 and B(DB+1-k) <> 0.0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, DB >= 0 and B(i) = 0.0, where */
/*                   i = 1, 2, ..., DB+1. */

/*     METHOD */

/*     Given real polynomials */
/*                                                  DA */
/*        A(x) = a(1) + a(2) * x + ... + a(DA+1) * x */

/*     and */
/*                                                  DB */
/*        B(x) = b(1) + b(2) * x + ... + b(DB+1) * x */

/*     where b(DB+1) is non-zero, the routine computes the coeffcients of */
/*     the quotient polynomial */
/*                                                     DA-DB */
/*        Q(x) = q(1) + q(2) * x + ... + q(DA-DB+1) * x */

/*     and the remainder polynomial */
/*                                                DB-1 */
/*        R(x) = r(1) + r(2) * x + ... + r(DB) * x */

/*     such that A(x) = B(x) * Q(x) + R(x). */

/*     The algorithm used is synthetic division of polynomials (see [1]), */
/*     which involves the following steps: */

/*        (a) compute q(k+1) = a(DB+k+1) / b(DB+1) */

/*     and */

/*        (b) set a(j) = a(j) - q(k+1) * b(j-k) for j = k+1, ..., DB+k. */

/*     Steps (a) and (b) are performed for k = DA-DB, DA-DB-1, ..., 0 and */
/*     the algorithm terminates with r(i) = a(i) for i = 1, 2, ..., DB. */

/*     REFERENCES */

/*     [1] Knuth, D.E. */
/*         The Art of Computer Programming, (Vol. 2, Seminumerical */
/*         Algorithms). */
/*         Addison-Wesley, Reading, Massachusetts (2nd Edition), 1981. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01ED by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    --rq;
    --b;
    --a;

    /* Function Body */
    *iwarn = 0;
    *info = 0;
    if (*da < -1) {
	*info = -1;
    } else if (*db < 0) {
	*info = -2;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC01QD", &i__1, (ftnlen)6);
	return 0;
    }

/*     WHILE ( DB >= 0 and B(DB+1) = 0 ) DO */
L20:
    if (*db >= 0) {
	if (b[*db + 1] == 0.) {
	    --(*db);
	    ++(*iwarn);
	    goto L20;
	}
    }
/*     END WHILE 20 */
    if (*db == -1) {
	*info = 1;
	return 0;
    }

/*     B(x) is non-zero. */

    if (*da >= 0) {
	n = *da;
	i__1 = n + 1;
	dcopy_(&i__1, &a[1], &c__1, &rq[1], &c__1);
/*        WHILE ( N >= DB ) DO */
L40:
	if (n >= *db) {
	    if (rq[n + 1] != 0.) {
		q = rq[n + 1] / b[*db + 1];
		d__1 = -q;
		daxpy_(db, &d__1, &b[1], &c__1, &rq[n - *db + 1], &c__1);
		rq[n + 1] = q;
	    }
	    --n;
	    goto L40;
	}
/*        END WHILE 40 */
    }

    return 0;
/* *** Last line of MC01QD *** */
} /* mc01qd_ */

