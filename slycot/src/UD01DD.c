/* UD01DD.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 0.;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;

/* Subroutine */ int ud01dd_(integer *m, integer *n, integer *nin, doublereal 
	*a, integer *lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Builtin functions */
    integer s_rsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsle(void);

    /* Local variables */
    static integer i__, j;
    static doublereal aij;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 1, 0, 0 };



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

/*     To read the elements of a sparse matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     NIN     (input) INTEGER */
/*             The input channel from which the elements of A are read. */
/*             NIN >= 0. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading M-by-N part of this array contains the sparse */
/*             matrix A. The not assigned elements are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : if a row index i is read with i < 1 or i > M or */
/*                   a column index j is read with j < 1 or j > N. */
/*                   This is a warning. */

/*     METHOD */

/*     First, the elements A(i,j) with 1 <= i <= M and 1 <= j <= N are */
/*     set to zero. Next the nonzero elements are read from the input */
/*     file NIN. Each line of NIN must contain consecutively the values */
/*     i, j, A(i,j). The routine terminates after the last line has been */
/*     read. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on routine RDSPAR by A.J. Geurts, Eindhoven University of */
/*     Technology, Holland. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;

/*     Check the input scalar arguments. */

    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nin < 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("UD01DD", &i__1, (ftnlen)6);
	return 0;
    }

    dlaset_("Full", m, n, &c_b4, &c_b4, &a[a_offset], lda, (ftnlen)4);

/*     Read (i, j, A(i,j)) of the nonzero elements one by one. */

L10:
    io___1.ciunit = *nin;
    i__1 = s_rsle(&io___1);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&aij, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L20;
    }
    if (i__ < 1 || i__ > *m || j < 1 || j > *n) {
	*info = 1;
    } else {
	a[i__ + j * a_dim1] = aij;
    }
    goto L10;
L20:

    return 0;
/* *** Last line of UD01DD *** */
} /* ud01dd_ */

