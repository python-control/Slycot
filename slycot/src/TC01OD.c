/* TC01OD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tc01od_(char *leri, integer *m, integer *p, integer *
	indlim, doublereal *pcoeff, integer *ldpco1, integer *ldpco2, 
	doublereal *qcoeff, integer *ldqco1, integer *ldqco2, integer *info, 
	ftnlen leri_len)
{
    /* System generated locals */
    integer pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, qcoeff_dim2,
	     qcoeff_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, porm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lleri;
    static integer mplim;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer minmp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To find the dual right (left) polynomial matrix representation of */
/*     a given left (right) polynomial matrix representation, where the */
/*     right and left polynomial matrix representations are of the form */
/*     Q(s)*inv(P(s)) and inv(P(s))*Q(s) respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether a left or right matrix fraction is input */
/*             as follows: */
/*             = 'L':  A left matrix fraction is input; */
/*             = 'R':  A right matrix fraction is input. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     INDLIM  (input) INTEGER */
/*             The highest value of K for which PCOEFF(.,.,K) and */
/*             QCOEFF(.,.,K) are to be transposed. */
/*             K = kpcoef + 1, where kpcoef is the maximum degree of the */
/*             polynomials in P(s).  INDLIM >= 1. */

/*     PCOEFF  (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,INDLIM) */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             On entry, the leading porm-by-porm-by-INDLIM part of this */
/*             array must contain the coefficients of the denominator */
/*             matrix P(s). */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of */
/*             polynomial (I,J) of P(s), where K = 1,2,...,INDLIM. */
/*             On exit, the leading porm-by-porm-by-INDLIM part of this */
/*             array contains the coefficients of the denominator matrix */
/*             P'(s) of the dual system. */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO1 >= MAX(1,M) if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO2 >= MAX(1,M) if LERI = 'R'. */

/*     QCOEFF  (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,INDLIM) */
/*             On entry, the leading P-by-M-by-INDLIM part of this array */
/*             must contain the coefficients of the numerator matrix */
/*             Q(s). */
/*             QCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of */
/*             polynomial (I,J) of Q(s), where K = 1,2,...,INDLIM. */
/*             On exit, the leading M-by-P-by-INDLIM part of the array */
/*             contains the coefficients of the numerator matrix Q'(s) */
/*             of the dual system. */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,M,P). */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If the given M-input/P-output left (right) polynomial matrix */
/*     representation has numerator matrix Q(s) and denominator matrix */
/*     P(s), its dual P-input/M-output right (left) polynomial matrix */
/*     representation simply has numerator matrix Q'(s) and denominator */
/*     matrix P'(s). */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TC01CD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    pcoeff_dim1 = *ldpco1;
    pcoeff_dim2 = *ldpco2;
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
    pcoeff -= pcoeff_offset;
    qcoeff_dim1 = *ldqco1;
    qcoeff_dim2 = *ldqco2;
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
    qcoeff -= qcoeff_offset;

    /* Function Body */
    *info = 0;
    lleri = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
    mplim = max(*m,*p);
    minmp = min(*m,*p);

/*     Test the input scalar arguments. */

    if (! lleri && ! lsame_(leri, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*indlim < 1) {
	*info = -4;
    } else if (lleri && *ldpco1 < max(1,*p) || ! lleri && *ldpco1 < max(1,*m))
	     {
	*info = -6;
    } else if (lleri && *ldpco2 < max(1,*p) || ! lleri && *ldpco2 < max(1,*m))
	     {
	*info = -7;
    } else if (*ldqco1 < max(1,mplim)) {
	*info = -9;
    } else if (*ldqco2 < max(1,mplim)) {
	*info = -10;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TC01OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *p == 0) {
	return 0;
    }

    if (mplim != 1) {

/*        Non-scalar system: transpose numerator matrix Q(s). */

	i__1 = *indlim;
	for (k = 1; k <= i__1; ++k) {

	    i__2 = mplim;
	    for (j = 1; j <= i__2; ++j) {
		if (j < minmp) {
		    i__3 = minmp - j;
		    dswap_(&i__3, &qcoeff[j + 1 + (j + k * qcoeff_dim2) * 
			    qcoeff_dim1], &c__1, &qcoeff[j + (j + 1 + k * 
			    qcoeff_dim2) * qcoeff_dim1], ldqco1);
		} else if (j > *p) {
		    dcopy_(p, &qcoeff[(j + k * qcoeff_dim2) * qcoeff_dim1 + 1]
			    , &c__1, &qcoeff[j + (k * qcoeff_dim2 + 1) * 
			    qcoeff_dim1], ldqco1);
		} else if (j > *m) {
		    dcopy_(m, &qcoeff[j + (k * qcoeff_dim2 + 1) * qcoeff_dim1]
			    , ldqco1, &qcoeff[(j + k * qcoeff_dim2) * 
			    qcoeff_dim1 + 1], &c__1);
		}
/* L10: */
	    }

/* L20: */
	}

/*        Find dimension of denominator matrix P(s): M (P) for */
/*        right (left) polynomial matrix representation. */

	porm = *m;
	if (lleri) {
	    porm = *p;
	}
	if (porm != 1) {

/*           Non-scalar P(s): transpose it. */

	    i__1 = *indlim;
	    for (k = 1; k <= i__1; ++k) {

		i__2 = porm - 1;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = porm - j;
		    dswap_(&i__3, &pcoeff[j + 1 + (j + k * pcoeff_dim2) * 
			    pcoeff_dim1], &c__1, &pcoeff[j + (j + 1 + k * 
			    pcoeff_dim2) * pcoeff_dim1], ldpco1);
/* L30: */
		}

/* L40: */
	    }

	}
    }

    return 0;
/* *** Last line of TC01OD *** */
} /* tc01od_ */

