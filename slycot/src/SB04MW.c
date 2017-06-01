/* SB04MW.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb04mw_(integer *m, doublereal *d__, integer *ipr, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal d1, d2;
    static integer i1, m1, m2, mpi, iprm;
    static doublereal mult;
    static integer iprm1;
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
/*     matrix is in upper Hessenberg form, stored compactly, row-wise. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the system.  M >= 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (M*(M+1)/2+2*M) */
/*             On entry, the first M*(M+1)/2 + M elements of this array */
/*             must contain an upper Hessenberg matrix, stored compactly, */
/*             row-wise, and the next M elements must contain the right */
/*             hand side of the linear system, as set by SLICOT Library */
/*             routine SB04MY. */
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

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB04AW by G. Golub, S. Nash, and */
/*     C. Van Loan, Stanford University, California, United States of */
/*     America, January 1982. */

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
    m1 = *m * (*m + 3) / 2;
    m2 = *m + *m;
    mpi = *m;
    iprm = m1;
    m1 = *m;
    i1 = 1;

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++mpi;
	++iprm;
	ipr[mpi] = i1;
	ipr[i__] = iprm;
	i1 += m1;
	if (i__ > 1) {
	    --m1;
	}
/* L20: */
    }

    m1 = *m - 1;
    mpi = *m;

/*     Reduce to upper triangular form. */

    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = i__ + 1;
	++mpi;
	iprm = ipr[mpi];
	iprm1 = ipr[mpi + 1];
	d1 = d__[iprm];
	d2 = d__[iprm1];
	if (abs(d1) <= abs(d2)) {

/*           Permute the row indices. */

	    k = iprm;
	    ipr[mpi] = iprm1;
	    iprm = iprm1;
	    iprm1 = k;
	    k = ipr[i__];
	    ipr[i__] = ipr[i1];
	    ipr[i1] = k;
	    d1 = d2;
	}

/*        Check singularity. */

	if (d1 == 0.) {
	    *info = 1;
	    return 0;
	}

	mult = -d__[iprm1] / d1;
	++iprm1;
	ipr[mpi + 1] = iprm1;

/*        Annihilate the subdiagonal elements of the matrix. */

	d__[ipr[i1]] += mult * d__[ipr[i__]];
	i__2 = *m - i__;
	daxpy_(&i__2, &mult, &d__[iprm + 1], &c__1, &d__[iprm1], &c__1);
/* L40: */
    }

/*     Check singularity. */

    if (d__[ipr[m2]] == 0.) {
	*info = 1;
	return 0;
    }

/*     Back substitution. */

    d__[ipr[*m]] /= d__[ipr[m2]];
    mpi = m2;

    for (i__ = m1; i__ >= 1; --i__) {
	--mpi;
	iprm = ipr[mpi];
	iprm1 = iprm;
	mult = 0.;

	i__1 = *m;
	for (i1 = i__ + 1; i1 <= i__1; ++i1) {
	    ++iprm1;
	    mult += d__[ipr[i1]] * d__[iprm1];
/* L60: */
	}

	d__[ipr[i__]] = (d__[ipr[i__]] - mult) / d__[iprm];
/* L80: */
    }

    return 0;
/* *** Last line of SB04MW *** */
} /* sb04mw_ */

