/* UD01ND.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ud01nd_(integer *mp, integer *np, integer *dp, integer *
	l, integer *nout, doublereal *p, integer *ldp1, integer *ldp2, char *
	text, integer *info, ftnlen text_len)
{
    /* Format strings */
    static char fmt_99999[] = "(\002 \002)";
    static char fmt_99998[] = "(/,1x,a,\002(\002,i2,\002)\002,\002 (\002,i2"
	    ",\002X\002,i2,\002)\002)";
    static char fmt_99997[] = "(5x,5(6x,i2,7x))";
    static char fmt_99996[] = "(1x,i2,2x,5d15.7)";

    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_len(char *, ftnlen), s_wsfe(cilist *), e_wsfe(void), do_fio(
	    integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, j1, j2, n1, jj, ltext;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lentxt;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_99999, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_99998, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_99996, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_99996, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_99999, 0 };



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

/*     To print the MP-by-NP coefficient matrices of a matrix polynomial */
/*                                                    dp-1           dp */
/*        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  . */

/*     The elements of the matrices are output to 7 significant figures. */

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

/*     L       (input) INTEGER */
/*             The number of elements of the coefficient matrices to be */
/*             printed per line.  1 <= L <= 5. */

/*     NOUT    (input) INTEGER */
/*             The output channel to which the results are sent. */
/*             NOUT >= 0. */

/*     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array must */
/*             contain the coefficients of the matrix polynomial P(s). */
/*             Specifically, P(i,j,k) must contain the coefficient of */
/*             s**(k-1) of the polynomial which is the (i,j)-th element */
/*             of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and */
/*             k = 1,2,...,DP+1. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MP. */

/*     LDP2    INTEGER */
/*             The second dimension of array P.  LDP2 >= NP. */

/*     TEXT    (input) CHARACTER*72 */
/*             Title caption of the coefficient matrices to be printed. */
/*             TEXT is followed by the degree of the coefficient matrix, */
/*             within brackets. If TEXT = ' ', then the coefficient */
/*             matrices are separated by an empty line. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     For i = 1, 2, ..., DP + 1 the routine first prints the contents of */
/*     TEXT followed by (i-1) as a title, followed by the elements of the */
/*     MP-by-NP coefficient matrix P(i) such that */
/*     (i)  if NP < L, then the leading MP-by-NP part is printed; */
/*     (ii) if NP = k*L + p (where k, p > 0), then k MP-by-L blocks of */
/*          consecutive columns of P(i) are printed one after another */
/*          followed by one MP-by-p block containing the last p columns */
/*          of P(i). */
/*     Row numbers are printed on the left of each row and a column */
/*     number on top of each column. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on routine PRMAPO by A.J. Geurts, Eindhoven University of */
/*     Technology, Holland. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

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
    } else if (*l < 1 || *l > 5) {
	*info = -4;
    } else if (*nout < 0) {
	*info = -5;
    } else if (*ldp1 < *mp) {
	*info = -7;
    } else if (*ldp2 < *np) {
	*info = -8;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("UD01ND", &i__1, (ftnlen)6);
	return 0;
    }

    lentxt = i_len(text, text_len);
    ltext = min(72,lentxt);
/*     WHILE ( TEXT(LTEXT:LTEXT) =  ' ' ) DO */
L10:
    if (*(unsigned char *)&text[ltext - 1] == ' ') {
	--ltext;
	goto L10;
    }
/*     END WHILE 10 */

    i__1 = *dp + 1;
    for (k = 1; k <= i__1; ++k) {
	if (ltext == 0) {
	    io___4.ciunit = *nout;
	    s_wsfe(&io___4);
	    e_wsfe();
	} else {
	    io___5.ciunit = *nout;
	    s_wsfe(&io___5);
	    do_fio(&c__1, text, ltext);
	    i__2 = k - 1;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*mp), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	n1 = (*np - 1) / *l;
	j1 = 1;
	j2 = *l;

	i__2 = n1;
	for (j = 1; j <= i__2; ++j) {
	    io___10.ciunit = *nout;
	    s_wsfe(&io___10);
	    i__3 = j2;
	    for (jj = j1; jj <= i__3; ++jj) {
		do_fio(&c__1, (char *)&jj, (ftnlen)sizeof(integer));
	    }
	    e_wsfe();

	    i__3 = *mp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		io___13.ciunit = *nout;
		s_wsfe(&io___13);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		i__4 = j2;
		for (jj = j1; jj <= i__4; ++jj) {
		    do_fio(&c__1, (char *)&p[i__ + (jj + k * p_dim2) * p_dim1]
			    , (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
/* L20: */
	    }

	    j1 += *l;
	    j2 += *l;
/* L30: */
	}

	io___14.ciunit = *nout;
	s_wsfe(&io___14);
	i__2 = *np;
	for (j = j1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	}
	e_wsfe();

	i__2 = *mp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___15.ciunit = *nout;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    i__3 = *np;
	    for (jj = j1; jj <= i__3; ++jj) {
		do_fio(&c__1, (char *)&p[i__ + (jj + k * p_dim2) * p_dim1], (
			ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
/* L40: */
	}

/* L50: */
    }

    io___16.ciunit = *nout;
    s_wsfe(&io___16);
    e_wsfe();

    return 0;


/* *** Last line of UD01ND *** */
} /* ud01nd_ */

