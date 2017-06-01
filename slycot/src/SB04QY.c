/* SB04QY.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int sb04qy_(integer *n, integer *m, integer *ind, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ipr, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer i__, j, k, i2, k1, k2, m1;
    static doublereal dum[1];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), sb04mw_(integer *, doublereal *, integer *, integer *)
	    , dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *), dtrmv_(char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);


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

/*     To construct and solve a linear algebraic system of order M whose */
/*     coefficient matrix is in upper Hessenberg form. Such systems */
/*     appear when solving discrete-time Sylvester equations using the */
/*     Hessenberg-Schur method. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix A.  M >= 0. */

/*     IND     (input) INTEGER */
/*             The index of the column in C to be computed.  IND >= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain an */
/*             upper Hessenberg matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading N-by-N part of this array must contain a */
/*             matrix in real Schur form. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the coefficient matrix C of the equation. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix C with column IND updated. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Workspace */

/*     D       DOUBLE PRECISION array, dimension (M*(M+1)/2+2*M) */

/*     IPR     INTEGER array, dimension (2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             > 0:  if INFO = IND, a singular matrix was encountered. */

/*     METHOD */

/*     A special linear algebraic system of order M, with coefficient */
/*     matrix in upper Hessenberg form is constructed and solved. The */
/*     coefficient matrix is stored compactly, row-wise. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     [2] Sima, V. */
/*         Algorithms for Linear-quadratic Optimization. */
/*         Marcel Dekker, Inc., New York, 1996. */

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
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --d__;
    --ipr;

    /* Function Body */
    if (*ind < *n) {
	dum[0] = 0.;
	dcopy_(m, dum, &c__0, &d__[1], &c__1);
	i__1 = *n;
	for (i__ = *ind + 1; i__ <= i__1; ++i__) {
	    daxpy_(m, &b[*ind + i__ * b_dim1], &c__[i__ * c_dim1 + 1], &c__1, 
		    &d__[1], &c__1);
/* L10: */
	}
	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    c__[i__ + *ind * c_dim1] -= a[i__ + (i__ - 1) * a_dim1] * d__[i__ 
		    - 1];
/* L20: */
	}
	dtrmv_("Upper", "No Transpose", "Non Unit", m, &a[a_offset], lda, &
		d__[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__ + *ind * c_dim1] -= d__[i__];
/* L30: */
	}
    }

    m1 = *m + 1;
    i2 = *m * m1 / 2 + m1;
    k2 = 1;
    k = *m;

/*     Construct the linear algebraic system of order M. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = m1 - k;
	dcopy_(&k, &a[i__ + j * a_dim1], lda, &d__[k2], &c__1);
	dscal_(&k, &b[*ind + *ind * b_dim1], &d__[k2], &c__1);
	k1 = k2;
	k2 += k;
	if (i__ > 1) {
	    ++k1;
	    --k;
	}
	d__[k1] += 1.;

/*        Store the right hand side. */

	d__[i2] = c__[i__ + *ind * c_dim1];
	++i2;
/* L40: */
    }

/*     Solve the linear algebraic system and store the solution in C. */

    sb04mw_(m, &d__[1], &ipr[1], info);

    if (*info != 0) {
	*info = *ind;
    } else {

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__ + *ind * c_dim1] = d__[ipr[i__]];
/* L50: */
	}

    }

    return 0;
/* *** Last line of SB04QY *** */
} /* sb04qy_ */

