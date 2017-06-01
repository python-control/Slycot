/* TF01PD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tf01pd_(integer *nh1, integer *nh2, integer *nr, integer 
	*nc, doublereal *h__, integer *ldh, doublereal *t, integer *ldt, 
	integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static integer ih, it, jt, ncol, nrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To construct the block Toeplitz expansion T of a multivariable */
/*     parameter sequence M(1),...,M(NR+NC-1), where each parameter M(k) */
/*     is an NH1-by-NH2 block matrix and k = 1,2,...,(NR+NC-1). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NH1     (input) INTEGER */
/*             The number of rows in each parameter M(k).  NH1 >= 0. */

/*     NH2     (input) INTEGER */
/*             The number of columns in each parameter M(k).  NH2 >= 0. */

/*     NR      (input) INTEGER */
/*             The number of parameters required in each column of the */
/*             block Toeplitz expansion matrix T.  NR >= 0. */

/*     NC      (input) INTEGER */
/*             The number of parameters required in each row of the */
/*             block Toeplitz expansion matrix T.  NC >= 0. */

/*     H       (input) DOUBLE PRECISION array, dimension */
/*             (LDH,(NR+NC-1)*NH2) */
/*             The leading NH1-by-(NR+NC-1)*NH2 part of this array must */
/*             contain the multivariable sequence M(k), where k = 1,2, */
/*             ...,(NR+NC-1). Specifically, each parameter M(k) is an */
/*             NH1-by-NH2 matrix whose (i,j)-th element must be stored in */
/*             H(i,(k-1)*NH2+j) for i = 1,2,...,NH1 and j = 1,2,...,NH2. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= MAX(1,NH1). */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,NH2*NC) */
/*             The leading NH1*NR-by-NH2*NC part of this array contains */
/*             the block Toeplitz expansion of the multivariable sequence */
/*             M(k). */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,NH1*NR). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The NH1-by-NH2 dimensional parameters M(k) of a multivariable */
/*     sequence are arranged into a matrix T in Toeplitz form such that */

/*                | M(NC)       M(NC-1)     M(NC-2)    . . .  M(1)  | */
/*                |                                                 | */
/*                | M(NC+1)     M(NC)       M(NC-1)    . . .  M(2)  | */
/*           T =  |  .           .           .                 .    |. */
/*                |  .           .           .                 .    | */
/*                |  .           .           .                 .    | */
/*                |                                                 | */
/*                | M(NR+NC-1)  M(NR+NC-2)  M(NR+NC-3) . . .  M(NR) | */

/*     REFERENCES */

/*     [1] Johvidov, J.S. */
/*         Hankel and Toeplitz Matrices and Forms: Algebraic Theory, */
/*         (translated by G.P.A. Thijsse, I. Gohberg, ed.). */
/*         Birkhaeuser, Boston, 1982. */

/*     NUMERICAL ASPECTS */

/*     The time taken is approximately proportional to */
/*     NH1 x NH2 x NR x NC. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TF01DD by S. Van Huffel, Katholieke */
/*     Univ. Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Multivariable system, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*nh1 < 0) {
	*info = -1;
    } else if (*nh2 < 0) {
	*info = -2;
    } else if (*nr < 0) {
	*info = -3;
    } else if (*nc < 0) {
	*info = -4;
    } else if (*ldh < max(1,*nh1)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *nh1 * *nr;
	if (*ldt < max(i__1,i__2)) {
	    *info = -8;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TF01PD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*nh1,*nh2), i__1 = max(i__1,*nr);
    if (max(i__1,*nc) == 0) {
	return 0;
    }

/*     Construct the last block column of T. */

    ih = 1;
    nrow = (*nr - 1) * *nh1;
    ncol = (*nc - 1) * *nh2 + 1;

    i__1 = nrow + *nh1;
    i__2 = *nh1;
    for (it = 1; i__2 < 0 ? it >= i__1 : it <= i__1; it += i__2) {
	dlacpy_("Full", nh1, nh2, &h__[ih * h_dim1 + 1], ldh, &t[it + ncol * 
		t_dim1], ldt, (ftnlen)4);
	ih += *nh2;
/* L10: */
    }

/*     Construct the remaining block columns of T in backward order. */

    i__2 = -(*nh2);
    for (jt = ncol - *nh2; i__2 < 0 ? jt >= 1 : jt <= 1; jt += i__2) {
	dlacpy_("Full", &nrow, nh2, &t[*nh1 + 1 + (jt + *nh2) * t_dim1], ldt, 
		&t[jt * t_dim1 + 1], ldt, (ftnlen)4);
	dlacpy_("Full", nh1, nh2, &h__[ih * h_dim1 + 1], ldh, &t[nrow + 1 + 
		jt * t_dim1], ldt, (ftnlen)4);
	ih += *nh2;
/* L20: */
    }

    return 0;
/* *** Last line of TF01PD *** */
} /* tf01pd_ */

