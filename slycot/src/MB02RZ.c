/* MB02RZ.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b1 = {1.,0.};

/* Subroutine */ int mb02rz_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *h__, integer *ldh, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1, i__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, jp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztrsm_(
	    char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static logical notran;


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

/*     To solve a system of linear equations */
/*        H * X = B,  H' * X = B  or  H**H * X = B */
/*     with a complex upper Hessenberg N-by-N matrix H using the LU */
/*     factorization computed by MB02SZ. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of the system of equations: */
/*             = 'N':  H * X = B  (No transpose) */
/*             = 'T':  H'* X = B  (Transpose) */
/*             = 'C':  H**H * X = B  (Conjugate transpose) */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of right hand sides, i.e., the number of */
/*             columns of the matrix B.  NRHS >= 0. */

/*     H       (input) COMPLEX*16 array, dimension (LDH,N) */
/*             The factors L and U from the factorization H = P*L*U */
/*             as computed by MB02SZ. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices from MB02SZ; for 1<=i<=N, row i of the */
/*             matrix was interchanged with row IPIV(i). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS) */
/*             On entry, the right hand side matrix B. */
/*             On exit, the solution matrix X. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     INFO    (output) INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses the factorization */
/*        H = P * L * U */
/*     where P is a permutation matrix, L is lower triangular with unit */
/*     diagonal elements (and one nonzero subdiagonal), and U is upper */
/*     triangular. */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N x NRHS ) complex operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01FW by A.J. Laub, University of */
/*     Southern California, United States of America, May 1980. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Frequency response, Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldh < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02RZ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve H * X = B. */

/*        Solve L * X = B, overwriting B with X. */

/*        L is represented as a product of permutations and unit lower */
/*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
/*        where each transformation L(i) is a rank-one modification of */
/*        the identity matrix. */

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    jp = ipiv[j];
	    if (jp != j) {
		zswap_(nrhs, &b[jp + b_dim1], ldb, &b[j + b_dim1], ldb);
	    }
	    i__2 = j + 1 + j * h_dim1;
	    z__1.r = -h__[i__2].r, z__1.i = -h__[i__2].i;
	    zaxpy_(nrhs, &z__1, &b[j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
/* L10: */
	}

/*        Solve U * X = B, overwriting B with X. */

	ztrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b1, &
		h__[h_offset], ldh, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8);

    } else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*        Solve H' * X = B. */

/*        Solve U' * X = B, overwriting B with X. */

	ztrsm_("Left", "Upper", trans, "Non-unit", n, nrhs, &c_b1, &h__[
		h_offset], ldh, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        Solve L' * X = B, overwriting B with X. */

	for (j = *n - 1; j >= 1; --j) {
	    i__1 = j + 1 + j * h_dim1;
	    z__1.r = -h__[i__1].r, z__1.i = -h__[i__1].i;
	    zaxpy_(nrhs, &z__1, &b[j + 1 + b_dim1], ldb, &b[j + b_dim1], ldb);
	    jp = ipiv[j];
	    if (jp != j) {
		zswap_(nrhs, &b[jp + b_dim1], ldb, &b[j + b_dim1], ldb);
	    }
/* L20: */
	}

    } else {

/*        Solve H**H * X = B. */

/*        Solve U**H * X = B, overwriting B with X. */

	ztrsm_("Left", "Upper", trans, "Non-unit", n, nrhs, &c_b1, &h__[
		h_offset], ldh, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        Solve L**H * X = B, overwriting B with X. */

	for (j = *n - 1; j >= 1; --j) {
	    d_cnjg(&z__2, &h__[j + 1 + j * h_dim1]);
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zaxpy_(nrhs, &z__1, &b[j + 1 + b_dim1], ldb, &b[j + b_dim1], ldb);
	    jp = ipiv[j];
	    if (jp != j) {
		zswap_(nrhs, &b[jp + b_dim1], ldb, &b[j + b_dim1], ldb);
	    }
/* L30: */
	}

    }

    return 0;
/* *** Last line of MB02RZ *** */
} /* mb02rz_ */

