/* MB02SD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02sd_(integer *n, doublereal *h__, integer *ldh, 
	integer *ipiv, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, jp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
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

/*     To compute an LU factorization of an n-by-n upper Hessenberg */
/*     matrix H using partial pivoting with row interchanges. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
/*             On entry, the n-by-n upper Hessenberg matrix to be */
/*             factored. */
/*             On exit, the factors L and U from the factorization */
/*             H = P*L*U; the unit diagonal elements of L are not stored, */
/*             and L is lower bidiagonal. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the matrix */
/*             was interchanged with row IPIV(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, U(i,i) is exactly zero. The */
/*                   factorization has been completed, but the factor U */
/*                   is exactly singular, and division by zero will occur */
/*                   if it is used to solve a system of equations. */

/*     METHOD */

/*     The factorization has the form */
/*        H = P * L * U */
/*     where P is a permutation matrix, L is lower triangular with unit */
/*     diagonal elements (and one nonzero subdiagonal), and U is upper */
/*     triangular. */

/*     This is the right-looking Level 1 BLAS version of the algorithm */
/*     (adapted after DGETF2). */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000, */
/*     Jan. 2005. */

/*     KEYWORDS */

/*     Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ldh < max(1,*n)) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02SD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	jp = j;
	if (j < *n) {
	    if ((d__1 = h__[j + 1 + j * h_dim1], abs(d__1)) > (d__2 = h__[j + 
		    j * h_dim1], abs(d__2))) {
		jp = j + 1;
	    }
	}
	ipiv[j] = jp;
	if (h__[jp + j * h_dim1] != 0.) {

/*           Apply the interchange to columns J:N. */

	    if (jp != j) {
		i__2 = *n - j + 1;
		dswap_(&i__2, &h__[j + j * h_dim1], ldh, &h__[jp + j * h_dim1]
			, ldh);
	    }

/*           Compute element J+1 of J-th column. */

	    if (j < *n) {
		h__[j + 1 + j * h_dim1] /= h__[j + j * h_dim1];
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < *n) {

/*           Update trailing submatrix. */

	    i__2 = *n - j;
	    d__1 = -h__[j + 1 + j * h_dim1];
	    daxpy_(&i__2, &d__1, &h__[j + (j + 1) * h_dim1], ldh, &h__[j + 1 
		    + (j + 1) * h_dim1], ldh);
	}
/* L10: */
    }
    return 0;
/* *** Last line of MB02SD *** */
} /* mb02sd_ */

