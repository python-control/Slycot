/* TB03AY.f -- translated by f2c (version 20100827).
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

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;
static doublereal c_b10 = -1.;
static integer c__1 = 1;

/* Subroutine */ int tb03ay_(integer *nr, doublereal *a, integer *lda, 
	integer *indblk, integer *nblk, doublereal *vcoeff, integer *ldvco1, 
	integer *ldvco2, doublereal *pcoeff, integer *ldpco1, integer *ldpco2,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, pcoeff_dim1, pcoeff_dim2, pcoeff_offset, 
	    vcoeff_dim1, vcoeff_dim2, vcoeff_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, ioff, joff, ncol, nrow;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrsm_(char *, char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer kplus, lwork, lstop;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer lstart, inplus;


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

/*     To calculate the (PWORK-by-NR) polynomial matrix V(s) one */
/*     (PWORK-by-NBLK(L-1)) block V:L-1(s) at a time, in reverse order */
/*     (L = INDBLK,...,1).  At each stage, the (NBLK(L)-by-NBLK(L)) poly- */
/*     nomial matrix W(s) = V2(s) * A2 is formed, where V2(s) is that */
/*     part of V(s) already computed and A2 is the subdiagonal (incl.) */
/*     part of the L-th column block of A; W(s) is temporarily stored in */
/*     the top left part of P(s), as is subsequently the further matrix */
/*     Wbar(s) = s * V:L(s) - W(s).  Then, except for the final stage */
/*     L = 1 (when the next step is to calculate P(s) itself, not here), */
/*     the top left part of V:L-1(s) is given by Wbar(s) * inv(R), where */
/*     R is the upper triangular part of the L-th superdiagonal block of */
/*     A.  Finally, note that the coefficient matrices W(.,.,K) can only */
/*     be non-zero for K = L + 1,...,INPLUS, with each of these matrices */
/*     having only its first NBLK(L-1) rows non-trivial.  Similarly, */
/*     Wbar(.,.,K) (and so clearly V:L-1(.,.,K) ) can only be non-zero */
/*     for K = L,...,INPLUS, with each of these having only its first */
/*     NBLK(K-1) rows non-trivial except for K = L, which has NBLK(L) */
/*     such rows. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     NOTE: In the interests of speed, this routine does not check the */
/*           inputs for errors. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --nblk;
    vcoeff_dim1 = *ldvco1;
    vcoeff_dim2 = *ldvco2;
    vcoeff_offset = 1 + vcoeff_dim1 * (1 + vcoeff_dim2);
    vcoeff -= vcoeff_offset;
    pcoeff_dim1 = *ldpco1;
    pcoeff_dim2 = *ldpco2;
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
    pcoeff -= pcoeff_offset;

    /* Function Body */
    *info = 0;
    inplus = *indblk + 1;
    joff = *nr;

/*     Calculate each column block V:LWORK-1(s) of V(s) in turn. */

    i__1 = *indblk;
    for (l = 1; l <= i__1; ++l) {
	lwork = inplus - l;

/*        Determine number of columns of V:LWORK(s) & its position in V. */

	ncol = nblk[lwork];
	joff -= ncol;

/*        Find limits for V2(s) * A2 calculation: skips zero rows */
/*        in V(s). */

	lstart = joff + 1;
	lstop = joff;

/*        Calculate W(s) and store (temporarily) in top left part */
/*        of P(s). */

	i__2 = inplus;
	for (k = lwork + 1; k <= i__2; ++k) {
	    nrow = nblk[k - 1];
	    lstop += nrow;
	    i__3 = lstop - lstart + 1;
	    dgemm_("No transpose", "No transpose", &nrow, &ncol, &i__3, &c_b6,
		     &vcoeff[(lstart + k * vcoeff_dim2) * vcoeff_dim1 + 1], 
		    ldvco1, &a[lstart + (joff + 1) * a_dim1], lda, &c_b7, &
		    pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (
		    ftnlen)12, (ftnlen)12);
/* L10: */
	}

/*        Replace W(s) by Wbar(s) = s * V:L(s) - W(s). */

	nrow = ncol;

	i__2 = *indblk;
	for (k = lwork; k <= i__2; ++k) {
	    kplus = k + 1;

	    i__3 = ncol;
	    for (j = 1; j <= i__3; ++j) {
		dscal_(&nrow, &c_b10, &pcoeff[(j + k * pcoeff_dim2) * 
			pcoeff_dim1 + 1], &c__1);
		daxpy_(&nrow, &c_b6, &vcoeff[(joff + j + kplus * vcoeff_dim2) 
			* vcoeff_dim1 + 1], &c__1, &pcoeff[(j + k * 
			pcoeff_dim2) * pcoeff_dim1 + 1], &c__1);
/* L20: */
	    }

	    nrow = nblk[k];
/* L30: */
	}

	i__2 = ncol;
	for (j = 1; j <= i__2; ++j) {
	    dscal_(&nrow, &c_b10, &pcoeff[(j + inplus * pcoeff_dim2) * 
		    pcoeff_dim1 + 1], &c__1);
/* L40: */
	}

	if (lwork != 1) {

/*           If not final stage, use the upper triangular R (from A) */
/*           to calculate V:L-1(s), finally storing this new block. */

	    ioff = joff - nblk[lwork - 1];

	    i__2 = ncol;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (a[ioff + i__ + (joff + i__) * a_dim1] == 0.) {

/*                 Error return. */

		    *info = i__;
		    return 0;
		}
/* L50: */
	    }

	    nrow = nblk[lwork];

	    i__2 = inplus;
	    for (k = lwork; k <= i__2; ++k) {
		dlacpy_("Full", &nrow, &ncol, &pcoeff[(k * pcoeff_dim2 + 1) * 
			pcoeff_dim1 + 1], ldpco1, &vcoeff[(ioff + 1 + k * 
			vcoeff_dim2) * vcoeff_dim1 + 1], ldvco1, (ftnlen)4);
		dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &nrow, &
			ncol, &c_b6, &a[ioff + 1 + (joff + 1) * a_dim1], lda, 
			&vcoeff[(ioff + 1 + k * vcoeff_dim2) * vcoeff_dim1 + 
			1], ldvco1, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)
			8);
		nrow = nblk[k];
/* L60: */
	    }

	}
/* L70: */
    }

    return 0;
/* *** Last line of TB03AY *** */
} /* tb03ay_ */

