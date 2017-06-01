/* TB04BV.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb04bv_(char *order, integer *p, integer *m, integer *md,
	 integer *ign, integer *ldign, integer *igd, integer *ldigd, 
	doublereal *gn, doublereal *gd, doublereal *d__, integer *ldd, 
	doublereal *tol, integer *info, ftnlen order_len)
{
    /* System generated locals */
    integer d_dim1, d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, 
	    i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, ii, nd, kk, km, nn;
    static doublereal dij, eps;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical ascend;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal toldef;
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

/*     To separate the strictly proper part G0 from the constant part D */
/*     of an P-by-M proper transfer function matrix G. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             of the transfer function matrix are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     Input/Output Parameters */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1, i.e., */
/*             MD = MAX(IGD(I,J)) + 1. */
/*                  I,J */

/*     IGN     (input/output) INTEGER array, dimension (LDIGN,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the degrees of the numerator polynomials in G: */
/*             the (i,j) element of IGN must contain the degree of the */
/*             numerator polynomial of the polynomial ratio G(i,j). */
/*             On exit, the leading P-by-M part of this array contains */
/*             the degrees of the numerator polynomials in G0. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (input) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             degrees of the denominator polynomials in G (and G0): */
/*             the (i,j) element of IGD contains the degree of the */
/*             denominator polynomial of the polynomial ratio G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             On entry, this array must contain the coefficients of the */
/*             numerator polynomials, Num(i,j), of the transfer function */
/*             matrix G. The polynomials are stored in a column-wise */
/*             order, i.e., Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), */
/*             Num(2,2), ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., */
/*             Num(P,M); MD memory locations are reserved for each */
/*             polynomial, hence, the (i,j) polynomial is stored starting */
/*             from the location ((j-1)*P+i-1)*MD+1. The coefficients */
/*             appear in increasing or decreasing order of the powers */
/*             of the indeterminate, according to ORDER. */
/*             On exit, this array contains the coefficients of the */
/*             numerator polynomials of the strictly proper part G0 of */
/*             the transfer function matrix G, stored similarly. */

/*     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array must contain the coefficients of the */
/*             denominator polynomials, Den(i,j), of the transfer */
/*             function matrix G. The polynomials are stored as for the */
/*             numerator polynomials. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array contains the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= max(1,P). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the degrees of */
/*             the numerators Num0(i,j) of the strictly proper part of */
/*             the transfer function matrix G. If the user sets TOL > 0, */
/*             then the given value of TOL is used as an absolute */
/*             tolerance; the leading coefficients with absolute value */
/*             less than TOL are considered neglijible. If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by TOLDEF = IGN(i,j)*EPS*NORM( Num(i,j) ) is used */
/*             instead, where EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH), and NORM denotes the infinity */
/*             norm (the maximum coefficient in absolute value). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the transfer function matrix is not proper; */
/*             = 2:  if a denominator polynomial is null. */

/*     METHOD */

/*     The (i,j) entry of the real matrix D is zero, if the degree of */
/*     Num(i,j), IGN(i,j), is less than the degree of Den(i,j), IGD(i,j), */
/*     and it is given by the ratio of the leading coefficients of */
/*     Num(i,j) and Den(i,j), if IGN(i,j) is equal to IGD(i,j), */
/*     for i = 1 : P, and for j = 1 : M. */

/*     FURTHER COMMENTS */

/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Based on the BIMASC Library routine TMPRP by A. Varga. */

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
	xerbla_("TB04BV", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*p,*m) == 0) {
	return 0;
    }

/*     Prepare the computation of the default tolerance. */

    toldef = *tol;
    if (toldef <= 0.) {
	eps = dlamch_("Epsilon", (ftnlen)7);
    }

    k = 1;

    if (ascend) {

/*        Polynomial coefficients are stored in increasing order. */

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		nn = ign[i__ + j * ign_dim1];
		nd = igd[i__ + j * igd_dim1];
		if (nn > nd) {

/*                 Error return: the transfer function matrix is */
/*                               not proper. */

		    *info = 1;
		    return 0;
		} else if (nn < nd || nd == 0 && gn[k] == 0.) {
		    d__[i__ + j * d_dim1] = 0.;
		} else {

/*                 Here NN = ND. */

		    kk = k + nn;

		    if (gd[kk] == 0.) {

/*                    Error return: the denominator is null. */

			*info = 2;
			return 0;
		    }

		    dij = gn[kk] / gd[kk];
		    d__[i__ + j * d_dim1] = dij;
		    gn[kk] = 0.;
		    if (nn > 0) {
			d__1 = -dij;
			daxpy_(&nn, &d__1, &gd[k], &c__1, &gn[k], &c__1);
			if (*tol <= 0.) {
			    toldef = (doublereal) nn * eps * (d__1 = gn[
				    idamax_(&nn, &gn[k], &c__1)], abs(d__1));
			}
			km = nn;
			i__3 = km;
			for (ii = 1; ii <= i__3; ++ii) {
			    --kk;
			    --nn;
			    if ((d__1 = gn[kk], abs(d__1)) > toldef) {
				goto L20;
			    }
/* L10: */
			}

L20:

			ign[i__ + j * ign_dim1] = nn;
		    }
		}
		k += *md;
/* L30: */
	    }

/* L40: */
	}

    } else {

/*        Polynomial coefficients are stored in decreasing order. */

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		nn = ign[i__ + j * ign_dim1];
		nd = igd[i__ + j * igd_dim1];
		if (nn > nd) {

/*                 Error return: the transfer function matrix is */
/*                               not proper. */

		    *info = 1;
		    return 0;
		} else if (nn < nd || nd == 0 && gn[k] == 0.) {
		    d__[i__ + j * d_dim1] = 0.;
		} else {

/*                 Here NN = ND. */

		    kk = k;

		    if (gd[kk] == 0.) {

/*                    Error return: the denominator is null. */

			*info = 2;
			return 0;
		    }

		    dij = gn[kk] / gd[kk];
		    d__[i__ + j * d_dim1] = dij;
		    gn[kk] = 0.;
		    if (nn > 0) {
			d__1 = -dij;
			daxpy_(&nn, &d__1, &gd[k + 1], &c__1, &gn[k + 1], &
				c__1);
			if (*tol <= 0.) {
			    toldef = (doublereal) nn * eps * (d__1 = gn[
				    idamax_(&nn, &gn[k + 1], &c__1)], abs(
				    d__1));
			}
			km = nn;
			i__3 = km;
			for (ii = 1; ii <= i__3; ++ii) {
			    ++kk;
			    --nn;
			    if ((d__1 = gn[kk], abs(d__1)) > toldef) {
				goto L60;
			    }
/* L50: */
			}

L60:

			ign[i__ + j * ign_dim1] = nn;
			i__3 = nn;
			for (ii = 0; ii <= i__3; ++ii) {
			    gn[k + ii] = gn[kk + ii];
/* L70: */
			}

		    }
		}
		k += *md;
/* L80: */
	    }

/* L90: */
	}

    }

    return 0;
/* *** Last line of TB04BV *** */
} /* tb04bv_ */

