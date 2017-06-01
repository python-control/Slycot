/* SB04RV.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;
static doublereal c_b9 = 1.;
static doublereal c_b11 = 0.;
static doublereal c_b19 = -1.;

/* Subroutine */ int sb04rv_(char *abschr, char *ul, integer *n, integer *m, 
	doublereal *c__, integer *ldc, integer *indx, doublereal *ab, integer 
	*ldab, doublereal *ba, integer *ldba, doublereal *d__, doublereal *
	dwork, ftnlen abschr_len, ftnlen ul_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, ba_dim1, ba_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);


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

/*     To construct the right-hand sides D for a system of equations in */
/*     quasi-Hessenberg form solved via SB04RX (case with 2 right-hand */
/*     sides). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ABSCHR  CHARACTER*1 */
/*             Indicates whether AB contains A or B, as follows: */
/*             = 'A':  AB contains A; */
/*             = 'B':  AB contains B. */

/*     UL      CHARACTER*1 */
/*             Indicates whether AB is upper or lower Hessenberg matrix, */
/*             as follows: */
/*             = 'U':  AB is upper Hessenberg; */
/*             = 'L':  AB is lower Hessenberg. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix B.  M >= 0. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,M) */
/*             The leading N-by-M part of this array must contain both */
/*             the not yet modified part of the coefficient matrix C of */
/*             the Sylvester equation X + AXB = C, and both the currently */
/*             computed part of the solution of the Sylvester equation. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     INDX    (input) INTEGER */
/*             The position of the first column/row of C to be used in */
/*             the construction of the right-hand side D. */

/*     AB      (input) DOUBLE PRECISION array, dimension (LDAB,*) */
/*             The leading N-by-N or M-by-M part of this array must */
/*             contain either A or B of the Sylvester equation */
/*             X + AXB = C. */

/*     LDAB    INTEGER */
/*             The leading dimension of array AB. */
/*             LDAB >= MAX(1,N) or LDAB >= MAX(1,M) (depending on */
/*             ABSCHR = 'A' or ABSCHR = 'B', respectively). */

/*     BA      (input) DOUBLE PRECISION array, dimension (LDBA,*) */
/*             The leading N-by-N or M-by-M part of this array must */
/*             contain either A or B of the Sylvester equation */
/*             X + AXB = C, the matrix not contained in AB. */

/*     LDBA    INTEGER */
/*             The leading dimension of array BA. */
/*             LDBA >= MAX(1,N) or LDBA >= MAX(1,M) (depending on */
/*             ABSCHR = 'B' or ABSCHR = 'A', respectively). */

/*     D       (output) DOUBLE PRECISION array, dimension (*) */
/*             The leading 2*N or 2*M part of this array (depending on */
/*             ABSCHR = 'B' or ABSCHR = 'A', respectively) contains the */
/*             right-hand side stored as a matrix with two rows. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             where LDWORK is equal to 2*N or 2*M (depending on */
/*             ABSCHR = 'B' or ABSCHR = 'A', respectively). */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     For speed, no tests on the input scalar arguments are made. */
/*     Quick return if possible. */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    ba_dim1 = *ldba;
    ba_offset = 1 + ba_dim1;
    ba -= ba_offset;
    --d__;
    --dwork;

    /* Function Body */
    if (*n == 0 || *m == 0) {
	return 0;
    }

    if (lsame_(abschr, "B", (ftnlen)1, (ftnlen)1)) {

/*        Construct the 2 columns of the right-hand side. */

	dcopy_(n, &c__[*indx * c_dim1 + 1], &c__1, &d__[1], &c__2);
	dcopy_(n, &c__[(*indx + 1) * c_dim1 + 1], &c__1, &d__[2], &c__2);
	if (lsame_(ul, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*indx > 1) {
		i__1 = *indx - 1;
		dgemv_("N", n, &i__1, &c_b9, &c__[c_offset], ldc, &ab[*indx * 
			ab_dim1 + 1], &c__1, &c_b11, &dwork[1], &c__1, (
			ftnlen)1);
		i__1 = *indx - 1;
		dgemv_("N", n, &i__1, &c_b9, &c__[c_offset], ldc, &ab[(*indx 
			+ 1) * ab_dim1 + 1], &c__1, &c_b11, &dwork[*n + 1], &
			c__1, (ftnlen)1);
		dgemv_("N", n, n, &c_b19, &ba[ba_offset], ldba, &dwork[1], &
			c__1, &c_b9, &d__[1], &c__2, (ftnlen)1);
		dgemv_("N", n, n, &c_b19, &ba[ba_offset], ldba, &dwork[*n + 1]
			, &c__1, &c_b9, &d__[2], &c__2, (ftnlen)1);
	    }
	} else {
	    if (*indx < *m - 1) {
		i__1 = *m - *indx - 1;
		dgemv_("N", n, &i__1, &c_b9, &c__[(*indx + 2) * c_dim1 + 1], 
			ldc, &ab[*indx + 2 + *indx * ab_dim1], &c__1, &c_b11, 
			&dwork[1], &c__1, (ftnlen)1);
		i__1 = *m - *indx - 1;
		dgemv_("N", n, &i__1, &c_b9, &c__[(*indx + 2) * c_dim1 + 1], 
			ldc, &ab[*indx + 2 + (*indx + 1) * ab_dim1], &c__1, &
			c_b11, &dwork[*n + 1], &c__1, (ftnlen)1);
		dgemv_("N", n, n, &c_b19, &ba[ba_offset], ldba, &dwork[1], &
			c__1, &c_b9, &d__[1], &c__2, (ftnlen)1);
		dgemv_("N", n, n, &c_b19, &ba[ba_offset], ldba, &dwork[*n + 1]
			, &c__1, &c_b9, &d__[2], &c__2, (ftnlen)1);
	    }
	}
    } else {

/*        Construct the 2 rows of the right-hand side. */

	dcopy_(m, &c__[*indx + c_dim1], ldc, &d__[1], &c__2);
	dcopy_(m, &c__[*indx + 1 + c_dim1], ldc, &d__[2], &c__2);
	if (lsame_(ul, "U", (ftnlen)1, (ftnlen)1)) {
	    if (*indx < *n - 1) {
		i__1 = *n - *indx - 1;
		dgemv_("T", &i__1, m, &c_b9, &c__[*indx + 2 + c_dim1], ldc, &
			ab[*indx + (*indx + 2) * ab_dim1], ldab, &c_b11, &
			dwork[1], &c__1, (ftnlen)1);
		i__1 = *n - *indx - 1;
		dgemv_("T", &i__1, m, &c_b9, &c__[*indx + 2 + c_dim1], ldc, &
			ab[*indx + 1 + (*indx + 2) * ab_dim1], ldab, &c_b11, &
			dwork[*m + 1], &c__1, (ftnlen)1);
		dgemv_("T", m, m, &c_b19, &ba[ba_offset], ldba, &dwork[1], &
			c__1, &c_b9, &d__[1], &c__2, (ftnlen)1);
		dgemv_("T", m, m, &c_b19, &ba[ba_offset], ldba, &dwork[*m + 1]
			, &c__1, &c_b9, &d__[2], &c__2, (ftnlen)1);
	    }
	} else {
	    if (*indx > 1) {
		i__1 = *indx - 1;
		dgemv_("T", &i__1, m, &c_b9, &c__[c_offset], ldc, &ab[*indx + 
			ab_dim1], ldab, &c_b11, &dwork[1], &c__1, (ftnlen)1);
		i__1 = *indx - 1;
		dgemv_("T", &i__1, m, &c_b9, &c__[c_offset], ldc, &ab[*indx + 
			1 + ab_dim1], ldab, &c_b11, &dwork[*m + 1], &c__1, (
			ftnlen)1);
		dgemv_("T", m, m, &c_b19, &ba[ba_offset], ldba, &dwork[1], &
			c__1, &c_b9, &d__[1], &c__2, (ftnlen)1);
		dgemv_("T", m, m, &c_b19, &ba[ba_offset], ldba, &dwork[*m + 1]
			, &c__1, &c_b9, &d__[2], &c__2, (ftnlen)1);
	    }
	}
    }

    return 0;
/* *** Last line of SB04RV *** */
} /* sb04rv_ */

