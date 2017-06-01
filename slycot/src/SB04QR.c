/* SB04QR.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb04qr_(integer *m, doublereal *d__, integer *ipr, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal d1, d2, d3;
    static integer i1, i2, m1, mpi, mpi1, mpi2;
    static doublereal dmax__;
    static integer iprm, iprm1;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


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

/*     To solve a linear algebraic system of order M whose coefficient */
/*     matrix has zeros below the third subdiagonal and zero elements on */
/*     the third subdiagonal with even column indices. The matrix is */
/*     stored compactly, row-wise. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the system.  M >= 0, M even. */
/*             Note that parameter M should have twice the value in the */
/*             original problem (see SLICOT Library routine SB04QU). */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (M*M/2+4*M) */
/*             On entry, the first M*M/2 + 3*M elements of this array */
/*             must contain the coefficient matrix, stored compactly, */
/*             row-wise, and the next M elements must contain the right */
/*             hand side of the linear system, as set by SLICOT Library */
/*             routine SB04QU. */
/*             On exit, the content of this array is updated, the last M */
/*             elements containing the solution with components */
/*             interchanged (see IPR). */

/*     IPR     (output) INTEGER array, dimension (2*M) */
/*             The leading M elements contain information about the */
/*             row interchanges performed for solving the system. */
/*             Specifically, the i-th component of the solution is */
/*             specified by IPR(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if a singular matrix was encountered. */

/*     METHOD */

/*     Gaussian elimination with partial pivoting is used. The rows of */
/*     the matrix are not actually permuted, only their indices are */
/*     interchanged in array IPR. */

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

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipr;
    --d__;

    /* Function Body */
    *info = 0;
    i2 = *m * *m / 2 + *m * 3;
    mpi = *m;
    iprm = i2;
    m1 = *m;
    i1 = 1;

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++mpi;
	++iprm;
	ipr[mpi] = i1;
	ipr[i__] = iprm;
	i1 += m1;
	if (i__ >= 4 && i__ % 2 == 0) {
	    m1 += -2;
	}
/* L20: */
    }

    m1 = *m - 1;
    mpi1 = *m + 1;

/*     Reduce to upper triangular form. */

    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mpi = mpi1;
	++mpi1;
	iprm = ipr[mpi];
	d1 = d__[iprm];
	i1 = 3;
	if (i__ % 2 == 0) {
	    i1 = 2;
	}
	if (i__ == m1) {
	    i1 = 1;
	}
	mpi2 = mpi + i1;
	l = 0;
	dmax__ = abs(d1);

	i__2 = mpi2;
	for (j = mpi1; j <= i__2; ++j) {
	    d2 = d__[ipr[j]];
	    d3 = abs(d2);
	    if (d3 > dmax__) {
		dmax__ = d3;
		d1 = d2;
		l = j - mpi;
	    }
/* L40: */
	}

/*        Check singularity. */

	if (dmax__ == 0.) {
	    *info = 1;
	    return 0;
	}

	if (l > 0) {

/*           Permute the row indices. */

	    k = iprm;
	    j = mpi + l;
	    iprm = ipr[j];
	    ipr[j] = k;
	    ipr[mpi] = iprm;
	    k = ipr[i__];
	    i2 = i__ + l;
	    ipr[i__] = ipr[i2];
	    ipr[i2] = k;
	}
	++iprm;

/*        Annihilate the subdiagonal elements of the matrix. */

	i2 = i__;
	d3 = d__[ipr[i__]];

	i__2 = mpi2;
	for (j = mpi1; j <= i__2; ++j) {
	    ++i2;
	    iprm1 = ipr[j];
	    dmax__ = -d__[iprm1] / d1;
	    d__[ipr[i2]] += dmax__ * d3;
	    i__3 = *m - i__;
	    daxpy_(&i__3, &dmax__, &d__[iprm], &c__1, &d__[iprm1 + 1], &c__1);
	    ++ipr[j];
/* L60: */
	}

/* L80: */
    }

    mpi = *m + *m;
    iprm = ipr[mpi];

/*     Check singularity. */

    if (d__[iprm] == 0.) {
	*info = 1;
	return 0;
    }

/*     Back substitution. */

    d__[ipr[*m]] /= d__[iprm];

    for (i__ = m1; i__ >= 1; --i__) {
	--mpi;
	iprm = ipr[mpi];
	iprm1 = iprm;
	dmax__ = 0.;

	i__1 = *m;
	for (k = i__ + 1; k <= i__1; ++k) {
	    ++iprm1;
	    dmax__ += d__[ipr[k]] * d__[iprm1];
/* L100: */
	}

	d__[ipr[i__]] = (d__[ipr[i__]] - dmax__) / d__[iprm];
/* L120: */
    }

    return 0;
/* *** Last line of SB04QR *** */
} /* sb04qr_ */

