/* TD03AY.f -- translated by f2c (version 20100827).
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

static doublereal c_b3 = 0.;
static doublereal c_b7 = 1.;

/* Subroutine */ int td03ay_(integer *mwork, integer *pwork, integer *index, 
	doublereal *dcoeff, integer *lddcoe, doublereal *ucoeff, integer *
	lduco1, integer *lduco2, integer *n, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, dcoeff_dim1, dcoeff_offset, ucoeff_dim1, ucoeff_dim2, 
	    ucoeff_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, ia, ja;
    static doublereal diag, temp;
    static integer jmax1;
    static doublereal umax1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ibias;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal absdia;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal absdmx, bignum;
    static integer indcur;
    static doublereal smlnum;


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

/*     Calculates a state-space representation for a (PWORK x MWORK) */
/*     transfer matrix given in the form of polynomial row vectors over */
/*     common denominators (not necessarily lcd's).  Such a description */
/*     is simply the polynomial matrix representation */

/*          T(s) = inv(D(s)) * U(s), */

/*     where D(s) is diagonal with (I,I)-th element D:I(s) of degree */
/*     INDEX(I); applying Wolovich's Observable Structure Theorem to */
/*     this left matrix fraction then yields an equivalent state-space */
/*     representation in observable companion form, of order */
/*     N = sum(INDEX(I)).  As D(s) is diagonal, the PWORK ordered */
/*     'non-trivial' columns of C and A are very simply calculated, these */
/*     submatrices being diagonal and (INDEX(I) x 1) - block diagonal, */
/*     respectively: finding B and D is also somewhat simpler than for */
/*     general P(s) as dealt with in TC04AD. Finally, the state-space */
/*     representation obtained here is not necessarily controllable */
/*     (as D(s) and U(s) are not necessarily relatively left prime), but */
/*     it is theoretically completely observable: however, its */
/*     observability matrix may be poorly conditioned, so it is safer */
/*     not to assume observability either. */

/*     REVISIONS */

/*     May 13, 1998. */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --index;
    dcoeff_dim1 = *lddcoe;
    dcoeff_offset = 1 + dcoeff_dim1;
    dcoeff -= dcoeff_offset;
    ucoeff_dim1 = *lduco1;
    ucoeff_dim2 = *lduco2;
    ucoeff_offset = 1 + ucoeff_dim1 * (1 + ucoeff_dim2);
    ucoeff -= ucoeff_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;

    /* Function Body */
    *info = 0;

/*     Initialize A and C to be zero, apart from 1's on the subdiagonal */
/*     of A. */

    dlaset_("Upper", n, n, &c_b3, &c_b3, &a[a_offset], lda, (ftnlen)5);
    if (*n > 1) {
	i__1 = *n - 1;
	i__2 = *n - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b3, &c_b7, &a[a_dim1 + 2], lda, (
		ftnlen)5);
    }

    dlaset_("Full", pwork, n, &c_b3, &c_b3, &c__[c_offset], ldc, (ftnlen)4);

/*     Calculate B and D, as well as 'non-trivial' elements of A and C. */
/*     Check if any leading coefficient of D(s) nearly zero: if so, exit. */
/*     Caution is taken to avoid overflow. */

    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
    bignum = 1. / smlnum;

    ibias = 2;
    ja = 0;

    i__1 = *pwork;
    for (i__ = 1; i__ <= i__1; ++i__) {
	absdia = (d__1 = dcoeff[i__ + dcoeff_dim1], abs(d__1));
	jmax1 = idamax_(mwork, &ucoeff[i__ + (ucoeff_dim2 + 1) * ucoeff_dim1],
		 lduco1);
	umax1 = (d__1 = ucoeff[i__ + (jmax1 + ucoeff_dim2) * ucoeff_dim1], 
		abs(d__1));
	if (absdia < smlnum || absdia < 1. && umax1 > absdia * bignum) {

/*           Error return. */

	    *info = i__;
	    return 0;
	}
	diag = 1. / dcoeff[i__ + dcoeff_dim1];
	indcur = index[i__];
	if (indcur != 0) {
	    ibias += indcur;
	    ja += indcur;
	    if (indcur >= 1) {
		jmax1 = idamax_(&indcur, &dcoeff[i__ + (dcoeff_dim1 << 1)], 
			lddcoe);
		absdmx = (d__1 = dcoeff[i__ + jmax1 * dcoeff_dim1], abs(d__1))
			;
		if (absdia >= 1.) {
		    if (umax1 > 1.) {
			if (absdmx / absdia > bignum / umax1) {

/*                       Error return. */

			    *info = i__;
			    return 0;
			}
		    }
		} else {
		    if (umax1 > 1.) {
			if (absdmx > bignum * absdia / umax1) {

/*                       Error return. */

			    *info = i__;
			    return 0;
			}
		    }
		}
	    }

/*           I-th 'non-trivial' sub-vector of A given from coefficients */
/*           of D:I(s), while I-th row block of B given from this and */
/*           row I of U(s). */

	    i__2 = indcur + 1;
	    for (k = 2; k <= i__2; ++k) {
		ia = ibias - k;
		temp = -diag * dcoeff[i__ + k * dcoeff_dim1];
		a[ia + ja * a_dim1] = temp;

		dcopy_(mwork, &ucoeff[i__ + (k * ucoeff_dim2 + 1) * 
			ucoeff_dim1], lduco1, &b[ia + b_dim1], ldb);
		daxpy_(mwork, &temp, &ucoeff[i__ + (ucoeff_dim2 + 1) * 
			ucoeff_dim1], lduco1, &b[ia + b_dim1], ldb);
/* L10: */
	    }

	    if (ja < *n) {
		a[ja + 1 + ja * a_dim1] = 0.;
	    }

/*           Finally, I-th 'non-trivial' entry of C and row of D obtained */
/*           also. */

	    c__[i__ + ja * c_dim1] = diag;
	}

	dcopy_(mwork, &ucoeff[i__ + (ucoeff_dim2 + 1) * ucoeff_dim1], lduco1, 
		&d__[i__ + d_dim1], ldd);
	dscal_(mwork, &diag, &d__[i__ + d_dim1], ldd);
/* L20: */
    }

    return 0;
/* *** Last line of TD03AY *** */
} /* td03ay_ */

