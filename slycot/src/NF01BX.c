/* NF01BX.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int nf01bx_(integer *n, integer *ipar, integer *lipar, 
	doublereal *dpar, integer *ldpar, doublereal *j, integer *ldj, 
	doublereal *x, integer *incx, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1;

    /* Local variables */
    static doublereal c__;
    static integer m;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
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

/*     To compute (J'*J + c*I)*x, where J is an m-by-n real matrix, c is */
/*     a real scalar, I is the n-by-n identity matrix, and x is a real */
/*     n-vector. */

/*     NOTE: this routine must have the same arguments as SLICOT Library */
/*     routine NF01BW. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J.  N >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain the number of rows M of the Jacobian */
/*                     matrix J.  M >= 0. */
/*             IPAR is provided for compatibility with SLICOT Library */
/*             routine MD03AD. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 1. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ,N) */
/*             The leading M-by-N part of this array must contain the */
/*             Jacobian matrix J. */

/*     LDJ     INTEGER */
/*             The leading dimension of the array J.  LDJ >= MAX(1,M). */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*abs(INCX)) */
/*             On entry, this incremented array must contain the */
/*             vector x. */
/*             On exit, this incremented array contains the value of the */
/*             matrix-vector product (J'*J + c*I)*x. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X.  INCX <> 0. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= M. */

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
/*     Oct. 2000, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Mar. 2002, Oct. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations. */

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
    } else if (*lipar < 1) {
	*info = -3;
    } else if (*ldpar < 1) {
	*info = -5;
    } else if (*incx == 0) {
	*info = -9;
    } else {
	m = ipar[1];
	if (m < 0) {
	    *info = -2;
	} else if (*ldj < max(1,m)) {
	    *info = -7;
	} else if (*ldwork < m) {
	    *info = -11;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("NF01BX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    c__ = dpar[1];
    if (m == 0) {

/*        Special case, void J: x <-- c*x. */

	dscal_(n, &c__, &x[1], incx);
	return 0;
    }

    dgemv_("NoTranspose", &m, n, &c_b4, &j[j_offset], ldj, &x[1], incx, &c_b5,
	     &dwork[1], &c__1, (ftnlen)11);
    dgemv_("Transpose", &m, n, &c_b4, &j[j_offset], ldj, &dwork[1], &c__1, &
	    c__, &x[1], incx, (ftnlen)9);
    return 0;

/* *** Last line of NF01BX *** */
} /* nf01bx_ */

