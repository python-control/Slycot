/* TB04BW.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb04bw_(char *order, integer *p, integer *m, integer *md,
	 integer *ign, integer *ldign, integer *igd, integer *ldigd, 
	doublereal *gn, doublereal *gd, doublereal *d__, integer *ldd, 
	integer *info, ftnlen order_len)
{
    /* System generated locals */
    integer d_dim1, d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, ii, nd, kk, km, nn;
    static doublereal dij;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical ascend;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To compute the sum of an P-by-M rational matrix G and a real */
/*     P-by-M matrix D. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             of the rational matrix are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     Input/Output Parameters */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1, i.e., */
/*             MD = MAX(IGN(I,J),IGD(I,J)) + 1. */
/*                  I,J */

/*     IGN     (input/output) INTEGER array, dimension (LDIGN,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the degrees of the numerator polynomials in G: */
/*             the (i,j) element of IGN must contain the degree of the */
/*             numerator polynomial of the polynomial ratio G(i,j). */
/*             On exit, the leading P-by-M part of this array contains */
/*             the degrees of the numerator polynomials in G + D. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (input) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             degrees of the denominator polynomials in G (and G + D): */
/*             the (i,j) element of IGD contains the degree of the */
/*             denominator polynomial of the polynomial ratio G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             On entry, this array must contain the coefficients of the */
/*             numerator polynomials, Num(i,j), of the rational matrix G. */
/*             The polynomials are stored in a column-wise order, i.e., */
/*             Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), Num(2,2), */
/*             ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., Num(P,M); */
/*             MD memory locations are reserved for each polynomial, */
/*             hence, the (i,j) polynomial is stored starting from the */
/*             location ((j-1)*P+i-1)*MD+1. The coefficients appear in */
/*             increasing or decreasing order of the powers of the */
/*             indeterminate, according to ORDER. */
/*             On exit, this array contains the coefficients of the */
/*             numerator polynomials of the rational matrix G + D, */
/*             stored similarly. */

/*     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array must contain the coefficients of the */
/*             denominator polynomials, Den(i,j), of the rational */
/*             matrix G. The polynomials are stored as for the */
/*             numerator polynomials. */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= max(1,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The (i,j) entry of the real matrix D is added to the (i,j) entry */
/*     of the matrix G, g(i,j), which is a ratio of two polynomials, */
/*     for i = 1 : P, and for j = 1 : M. If g(i,j) = 0, it is assumed */
/*     that its denominator is 1. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     FURTHER COMMENTS */

/*     Often, the rational matrix G is found from a state-space */
/*     representation (A,B,C), and D corresponds to the direct */
/*     feedthrough matrix of the system. The sum G + D gives the */
/*     transfer function matrix of the system (A,B,C,D). */
/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Based on the BIMASC Library routine TMCADD by A. Varga. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     State-space representation, transfer function. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

    /* Parameter adjustments */
    ign_dim1 = *ldign;
    ign_offset = 1 + ign_dim1;
    ign -= ign_offset;
    igd_dim1 = *ldigd;
    igd_offset = 1 + igd_dim1;
    igd -= igd_offset;
    --gn;
    --gd;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;

    /* Function Body */
    *info = 0;
    ascend = lsame_(order, "I", (ftnlen)1, (ftnlen)1);
    if (! ascend && ! lsame_(order, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*p < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*md < 1) {
	*info = -4;
    } else if (*ldign < max(1,*p)) {
	*info = -6;
    } else if (*ldigd < max(1,*p)) {
	*info = -8;
    } else if (*ldd < max(1,*p)) {
	*info = -12;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB04BW", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*p,*m) == 0) {
	return 0;
    }

    k = 1;

    if (ascend) {

/*        Polynomial coefficients are stored in increasing order. */

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dij = d__[i__ + j * d_dim1];
		if (dij != 0.) {
		    nn = ign[i__ + j * ign_dim1];
		    nd = igd[i__ + j * igd_dim1];
		    if (nn == 0 && nd == 0) {
			if (gn[k] == 0.) {
			    gn[k] = dij;
			} else {
			    gn[k] += dij * gd[k];
			}
		    } else {
			km = min(nn,nd) + 1;
			daxpy_(&km, &dij, &gd[k], &c__1, &gn[k], &c__1);
			if (nn < nd) {

			    i__3 = k + nd;
			    for (ii = k + km; ii <= i__3; ++ii) {
				gn[ii] = dij * gd[ii];
/* L10: */
			    }

			    ign[i__ + j * ign_dim1] = nd;
			}
		    }
		}
		k += *md;
/* L20: */
	    }

/* L30: */
	}

    } else {

/*        Polynomial coefficients are stored in decreasing order. */

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dij = d__[i__ + j * d_dim1];
		if (dij != 0.) {
		    nn = ign[i__ + j * ign_dim1];
		    nd = igd[i__ + j * igd_dim1];
		    if (nn == 0 && nd == 0) {
			if (gn[k] == 0.) {
			    gn[k] = dij;
			} else {
			    gn[k] += dij * gd[k];
			}
		    } else {
			km = min(nn,nd) + 1;
			if (nn < nd) {
			    kk = k + nd - nn;

			    i__3 = k;
			    for (ii = k + nn; ii >= i__3; --ii) {
				gn[ii + nd - nn] = gn[ii];
/* L35: */
			    }

			    i__3 = kk - 1;
			    for (ii = k; ii <= i__3; ++ii) {
				gn[ii] = dij * gd[ii];
/* L40: */
			    }

			    ign[i__ + j * ign_dim1] = nd;
			    daxpy_(&km, &dij, &gd[kk], &c__1, &gn[kk], &c__1);
			} else {
			    kk = k + nn - nd;
			    daxpy_(&km, &dij, &gd[k], &c__1, &gn[kk], &c__1);
			}
		    }
		}
		k += *md;
/* L50: */
	    }

/* L60: */
	}

    }

    return 0;
/* *** Last line of TB04BW *** */
} /* tb04bw_ */

