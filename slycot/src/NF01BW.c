/* NF01BW.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int nf01bw_(integer *n, integer *ipar, integer *lipar, 
	doublereal *dpar, integer *ldpar, doublereal *j, integer *ldj, 
	doublereal *x, integer *incx, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer m, bn, jl, ix, xl, st, bsm, bsn, ibsm, ibsn, nths;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute the matrix-vector product x <-- (J'*J + c*I)*x, for the */
/*     Jacobian J as received from SLICOT Library routine NF01BD: */

/*          /  dy(1)/dwb(1)  |  dy(1)/dtheta  \ */
/*     Jc = |       :        |       :        | . */
/*          \  dy(L)/dwb(L)  |  dy(L)/dtheta  / */

/*     This is a compressed representation of the actual structure */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_L  |  L_L  / */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the vector x. */
/*             N = BN*BSN + ST >= 0.  (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain ST, the number of parameters */
/*                     corresponding to the linear part.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, BN = L, */
/*                     for the parameters corresponding to the nonlinear */
/*                     part.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     J_k = dy(k)/dwb(k), k = 1:BN, if BN > 0, or the */
/*                     number of rows of the matrix J, if BN <= 1. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             The leading NR-by-NC part of this array must contain */
/*             the (compressed) representation (Jc) of the Jacobian */
/*             matrix J, where NR = BSM if BN <= 1, and NR = BN*BSM, */
/*             if BN > 1. */

/*     LDJ     (input) INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NR). */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*INCX) */
/*             On entry, this incremented array must contain the */
/*             vector x. */
/*             On exit, this incremented array contains the value of the */
/*             matrix-vector product (J'*J + c*I)*x. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X.  INCX >= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= NR. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The associativity of matrix multiplications is used; the result */
/*     is obtained as:  x_out = J'*( J*x ) + c*x. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Mar. 2001, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001, */
/*     Mar. 2002. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations, */
/*     Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipar;
    --dpar;
    j_dim1 = *ldj;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --x;
    --dwork;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*lipar < 4) {
	*info = -3;
    } else if (*ldpar < 1) {
	*info = -5;
    } else if (*incx < 1) {
	*info = -9;
    } else {
	st = ipar[1];
	bn = ipar[2];
	bsm = ipar[3];
	bsn = ipar[4];
	nths = bn * bsn;
	if (bn > 1) {
	    m = bn * bsm;
	} else {
	    m = bsm;
	}
/* Computing MIN */
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
	if (min(i__1,bsn) < 0) {
	    *info = -2;
	} else if (*n != nths + st) {
	    *info = -1;
	} else if (*ldj < max(1,m)) {
	    *info = -7;
	} else if (*ldwork < m) {
	    *info = -11;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BW", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    c__ = dpar[1];

    if (m == 0) {

/*        Special case, void Jacobian: x <-- c*x. */

	dscal_(n, &c__, &x[1], incx);
	return 0;
    }

    if (bn <= 1 || bsn == 0) {

/*        Special case, l <= 1 or BSN = 0: the Jacobian is represented */
/*        as a full matrix. Adapted code from NF01BX is included in-line. */

	dgemv_("NoTranspose", &m, n, &c_b4, &j[j_offset], ldj, &x[1], incx, &
		c_b5, &dwork[1], &c__1, (ftnlen)11);
	dgemv_("Transpose", &m, n, &c_b4, &j[j_offset], ldj, &dwork[1], &c__1,
		 &c__, &x[1], incx, (ftnlen)9);
	return 0;
    }

/*     General case: l > 1, BSN > 0, BSM > 0. */

    jl = bsn + 1;
    ix = bsn * *incx;
    xl = bn * ix + 1;

    if (st > 0) {
	dgemv_("NoTranspose", &m, &st, &c_b4, &j[jl * j_dim1 + 1], ldj, &x[xl]
		, incx, &c_b5, &dwork[1], &c__1, (ftnlen)11);
    } else {
	dwork[1] = 0.;
	dcopy_(&m, &dwork[1], &c__0, &dwork[1], &c__1);
    }
    ibsn = 1;

    i__1 = m;
    i__2 = bsm;
    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += i__2) {
	dgemv_("NoTranspose", &bsm, &bsn, &c_b4, &j[ibsm + j_dim1], ldj, &x[
		ibsn], incx, &c_b4, &dwork[ibsm], &c__1, (ftnlen)11);
	dgemv_("Transpose", &bsm, &bsn, &c_b4, &j[ibsm + j_dim1], ldj, &dwork[
		ibsm], &c__1, &c__, &x[ibsn], incx, (ftnlen)9);
	ibsn += ix;
/* L10: */
    }

    if (st > 0) {
	dgemv_("Transpose", &m, &st, &c_b4, &j[jl * j_dim1 + 1], ldj, &dwork[
		1], &c__1, &c__, &x[xl], incx, (ftnlen)9);
    }

    return 0;

/* *** Last line of NF01BW *** */
} /* nf01bw_ */

