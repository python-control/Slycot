/* UD01BD.f -- translated by f2c (version 20100827).
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

static integer c__5 = 5;
static integer c__1 = 1;

/* Subroutine */ int ud01bd_(integer *mp, integer *np, integer *dp, integer *
	nin, doublereal *p, integer *ldp1, integer *ldp2, integer *info)
{
    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rsfe(cilist *), e_rsfe(void), s_rsle(cilist *), do_lio(integer *
	    , integer *, char *, ftnlen), e_rsle(void);

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, "()", 0 };
    static cilist io___4 = { 0, 0, 0, 0, 0 };



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

/*     To read the coefficients of a matrix polynomial */
/*                                                    dp-1           dp */
/*        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  . */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     MP      (input) INTEGER */
/*             The number of rows of the matrix polynomial P(s). */
/*             MP >= 1. */

/*     NP      (input) INTEGER */
/*             The number of columns of the matrix polynomial P(s). */
/*             NP >= 1. */

/*     DP      (input) INTEGER */
/*             The degree of the matrix polynomial P(s).  DP >= 0. */

/*     NIN     (input) INTEGER */
/*             The input channel from which the elements of P(s) are */
/*             read.  NIN >= 0. */

/*     P       (output) DOUBLE PRECISION array, dimension */
/*             (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array contains */
/*             the coefficients of the matrix polynomial P(s). */
/*             Specifically, P(i,j,k) contains the coefficient of */
/*             s**(k-1) of the polynomial which is the (i,j)-th element */
/*             of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and */
/*             k = 1,2,...,DP+1. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MP. */

/*     LDP2    INTEGER */
/*             The second dimension of array P.  LDP2 >= NP. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The coefficients P(i), i = 0, ..., DP, which are MP-by-NP */
/*     matrices, are read from the input file NIN row by row. Each P(i) */
/*     must be preceded by a text line. This text line can be used to */
/*     indicate the coefficient matrices. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on routine RDMAPO by A.J. Geurts, Eindhoven University of */
/*     Technology, Holland. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    p_dim1 = *ldp1;
    p_dim2 = *ldp2;
    p_offset = 1 + p_dim1 * (1 + p_dim2);
    p -= p_offset;

    /* Function Body */
    *info = 0;

/*     Check the input scalar arguments. */

    if (*mp < 1) {
	*info = -1;
    } else if (*np < 1) {
	*info = -2;
    } else if (*dp < 0) {
	*info = -3;
    } else if (*nin < 0) {
	*info = -4;
    } else if (*ldp1 < *mp) {
	*info = -6;
    } else if (*ldp2 < *np) {
	*info = -7;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("UD01BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Skip the text line preceding P(i) and read P(i), i = 0, ..., DP, */
/*     row after row. */

    i__1 = *dp + 1;
    for (k = 1; k <= i__1; ++k) {
	io___2.ciunit = *nin;
	s_rsfe(&io___2);
	e_rsfe();

	i__2 = *mp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___4.ciunit = *nin;
	    s_rsle(&io___4);
	    i__3 = *np;
	    for (j = 1; j <= i__3; ++j) {
		do_lio(&c__5, &c__1, (char *)&p[i__ + (j + k * p_dim2) * 
			p_dim1], (ftnlen)sizeof(doublereal));
	    }
	    e_rsle();
/* L10: */
	}

/* L20: */
    }

    return 0;
/* *** Last line of UD01BD *** */
} /* ud01bd_ */

