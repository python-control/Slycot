/* MC03MD.f -- translated by f2c (version 20100827).
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
static doublereal c_b8 = 0.;

/* Subroutine */ int mc03md_(integer *rp1, integer *cp1, integer *cp2, 
	integer *dp1, integer *dp2, integer *dp3, doublereal *alpha, 
	doublereal *p1, integer *ldp11, integer *ldp12, doublereal *p2, 
	integer *ldp21, integer *ldp22, doublereal *p3, integer *ldp31, 
	integer *ldp32, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer p1_dim1, p1_dim2, p1_offset, p2_dim1, p2_dim2, p2_offset, p3_dim1,
	     p3_dim2, p3_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer e, h__, i__, j, k;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer dpol3;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static logical cfzero;


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

/*     To compute the coefficients of the real polynomial matrix */

/*        P(x) = P1(x) * P2(x) + alpha * P3(x), */

/*     where P1(x), P2(x) and P3(x) are given real polynomial matrices */
/*     and alpha is a real scalar. */

/*     Each of the polynomial matrices P1(x), P2(x) and P3(x) may be the */
/*     zero matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     RP1     (input) INTEGER */
/*             The number of rows of the matrices P1(x) and P3(x). */
/*             RP1 >= 0. */

/*     CP1     (input) INTEGER */
/*             The number of columns of matrix P1(x) and the number of */
/*             rows of matrix P2(x).  CP1 >= 0. */

/*     CP2     (input) INTEGER */
/*             The number of columns of the matrices P2(x) and P3(x). */
/*             CP2 >= 0. */

/*     DP1     (input) INTEGER */
/*             The degree of the polynomial matrix P1(x).  DP1 >= -1. */

/*     DP2     (input) INTEGER */
/*             The degree of the polynomial matrix P2(x).  DP2 >= -1. */

/*     DP3     (input/output) INTEGER */
/*             On entry, the degree of the polynomial matrix P3(x). */
/*             DP3 >= -1. */
/*             On exit, the degree of the polynomial matrix P(x). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar value alpha of the problem. */

/*     P1      (input) DOUBLE PRECISION array, dimension (LDP11,LDP12,*) */
/*             If DP1 >= 0, then the leading RP1-by-CP1-by-(DP1+1) part */
/*             of this array must contain the coefficients of the */
/*             polynomial matrix P1(x). Specifically, P1(i,j,k) must */
/*             contain the coefficient of x**(k-1) of the polynomial */
/*             which is the (i,j)-th element of P1(x), where i = 1,2,..., */
/*             RP1, j = 1,2,...,CP1 and k = 1,2,...,DP1+1. */
/*             If DP1 = -1, then P1(x) is taken to be the zero polynomial */
/*             matrix, P1 is not referenced and can be supplied as a */
/*             dummy array (i.e. set the parameters LDP11 = LDP12 = 1 and */
/*             declare this array to be P1(1,1,1) in the calling */
/*             program). */

/*     LDP11   INTEGER */
/*             The leading dimension of array P1. */
/*             LDP11 >= MAX(1,RP1) if DP1 >= 0, */
/*             LDP11 >= 1          if DP1 = -1. */

/*     LDP12   INTEGER */
/*             The second dimension of array P1. */
/*             LDP12 >= MAX(1,CP1) if DP1 >= 0, */
/*             LDP12 >= 1          if DP1 = -1. */

/*     P2      (input) DOUBLE PRECISION array, dimension (LDP21,LDP22,*) */
/*             If DP2 >= 0, then the leading CP1-by-CP2-by-(DP2+1) part */
/*             of this array must contain the coefficients of the */
/*             polynomial matrix P2(x). Specifically, P2(i,j,k) must */
/*             contain the coefficient of x**(k-1) of the polynomial */
/*             which is the (i,j)-th element of P2(x), where i = 1,2,..., */
/*             CP1, j = 1,2,...,CP2 and k = 1,2,...,DP2+1. */
/*             If DP2 = -1, then P2(x) is taken to be the zero polynomial */
/*             matrix, P2 is not referenced and can be supplied as a */
/*             dummy array (i.e. set the parameters LDP21 = LDP22 = 1 and */
/*             declare this array to be P2(1,1,1) in the calling */
/*             program). */

/*     LDP21   INTEGER */
/*             The leading dimension of array P2. */
/*             LDP21 >= MAX(1,CP1) if DP2 >= 0, */
/*             LDP21 >= 1          if DP2 = -1. */

/*     LDP22   INTEGER */
/*             The second dimension of array P2. */
/*             LDP22 >= MAX(1,CP2) if DP2 >= 0, */
/*             LDP22 >= 1          if DP2 = -1. */

/*     P3      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDP31,LDP32,n), where n = MAX(DP1+DP2,DP3,0)+1. */
/*             On entry, if DP3 >= 0, then the leading */
/*             RP1-by-CP2-by-(DP3+1) part of this array must contain the */
/*             coefficients of the polynomial matrix P3(x). Specifically, */
/*             P3(i,j,k) must contain the coefficient of x**(k-1) of the */
/*             polynomial which is the (i,j)-th element of P3(x), where */
/*             i = 1,2,...,RP1, j = 1,2,...,CP2 and k = 1,2,...,DP3+1. */
/*             If DP3 = -1, then P3(x) is taken to be the zero polynomial */
/*             matrix. */
/*             On exit, if DP3 >= 0 on exit (ALPHA <> 0.0 and DP3 <> -1, */
/*             on entry, or DP1 <> -1 and DP2 <> -1), then the leading */
/*             RP1-by-CP2-by-(DP3+1) part of this array contains the */
/*             coefficients of P(x). Specifically, P3(i,j,k) contains the */
/*             coefficient of x**(k-1) of the polynomial which is the */
/*             (i,j)-th element of P(x), where i = 1,2,...,RP1, j = 1,2, */
/*             ...,CP2 and k = 1,2,...,DP3+1. */
/*             If DP3 = -1 on exit, then the coefficients of P(x) (the */
/*             zero polynomial matrix) are not stored in the array. */

/*     LDP31   INTEGER */
/*             The leading dimension of array P3.  LDP31 >= MAX(1,RP1). */

/*     LDP32   INTEGER */
/*             The second dimension of array P3.   LDP32 >= MAX(1,CP2). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (CP1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given real polynomial matrices */

/*                DP1            i */
/*        P1(x) = SUM (A(i+1) * x ), */
/*                i=0 */

/*                DP2            i */
/*        P2(x) = SUM (B(i+1) * x ), */
/*                i=0 */

/*                DP3            i */
/*        P3(x) = SUM (C(i+1) * x ) */
/*                i=0 */

/*     and a real scalar alpha, the routine computes the coefficients */
/*     d ,d ,..., of the polynomial matrix */
/*      1  2 */

/*        P(x) = P1(x) * P2(x) + alpha * P3(x) */

/*     from the formula */

/*                 s */
/*        d    =  SUM (A(k+1) * B(i-k+1)) + alpha * C(i+1), */
/*         i+1    k=r */

/*     where i = 0,1,...,DP1+DP2 and r and s depend on the value of i */
/*     (e.g. if i <= DP1 and i <= DP2, then r = 0 and s = i). */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     FURTHER COMMENTS */

/*     Other elementary operations involving polynomial matrices can */
/*     easily be obtained by calling the appropriate BLAS routine(s). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC03AD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, input output description, */
/*     polynomial matrix, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    p1_dim1 = *ldp11;
    p1_dim2 = *ldp12;
    p1_offset = 1 + p1_dim1 * (1 + p1_dim2);
    p1 -= p1_offset;
    p2_dim1 = *ldp21;
    p2_dim2 = *ldp22;
    p2_offset = 1 + p2_dim1 * (1 + p2_dim2);
    p2 -= p2_offset;
    p3_dim1 = *ldp31;
    p3_dim2 = *ldp32;
    p3_offset = 1 + p3_dim1 * (1 + p3_dim2);
    p3 -= p3_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*rp1 < 0) {
	*info = -1;
    } else if (*cp1 < 0) {
	*info = -2;
    } else if (*cp2 < 0) {
	*info = -3;
    } else if (*dp1 < -1) {
	*info = -4;
    } else if (*dp2 < -1) {
	*info = -5;
    } else if (*dp3 < -1) {
	*info = -6;
    } else if (*dp1 == -1 && *ldp11 < 1 || *dp1 >= 0 && *ldp11 < max(1,*rp1)) 
	    {
	*info = -9;
    } else if (*dp1 == -1 && *ldp12 < 1 || *dp1 >= 0 && *ldp12 < max(1,*cp1)) 
	    {
	*info = -10;
    } else if (*dp2 == -1 && *ldp21 < 1 || *dp2 >= 0 && *ldp21 < max(1,*cp1)) 
	    {
	*info = -12;
    } else if (*dp2 == -1 && *ldp22 < 1 || *dp2 >= 0 && *ldp22 < max(1,*cp2)) 
	    {
	*info = -13;
    } else if (*ldp31 < max(1,*rp1)) {
	*info = -15;
    } else if (*ldp32 < max(1,*cp2)) {
	*info = -16;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC03MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*rp1 == 0 || *cp2 == 0) {
	return 0;
    }

    if (*alpha == 0.) {
	*dp3 = -1;
    }

    if (*dp3 >= 0) {

/*        P3(x) := ALPHA * P3(x). */

	i__1 = *dp3 + 1;
	for (k = 1; k <= i__1; ++k) {

	    i__2 = *cp2;
	    for (j = 1; j <= i__2; ++j) {
		dscal_(rp1, alpha, &p3[(j + k * p3_dim2) * p3_dim1 + 1], &
			c__1);
/* L20: */
	    }

/* L40: */
	}
    }

    if (*dp1 == -1 || *dp2 == -1 || *cp1 == 0) {
	return 0;
    }

/*     Neither of P1(x) and P2(x) is the zero polynomial. */

    dpol3 = *dp1 + *dp2;
    if (dpol3 > *dp3) {

/*        Initialize the additional part of P3(x) to zero. */

	i__1 = dpol3 + 1;
	for (k = *dp3 + 2; k <= i__1; ++k) {
	    dlaset_("Full", rp1, cp2, &c_b8, &c_b8, &p3[(k * p3_dim2 + 1) * 
		    p3_dim1 + 1], ldp31, (ftnlen)4);
/* L80: */
	}

	*dp3 = dpol3;
    }
/*                                                              k-1 */
/*     The inner product of the j-th row of the coefficient of x    of P1 */
/*                                                i-1 */
/*     and the h-th column of the coefficient of x    of P2(x) contribute */
/*                                                 k+i-2 */
/*     the (j,h)-th element of the coefficient of x      of P3(x). */

    i__1 = *dp1 + 1;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *rp1;
	for (j = 1; j <= i__2; ++j) {
	    dcopy_(cp1, &p1[j + (k * p1_dim2 + 1) * p1_dim1], ldp11, &dwork[1]
		    , &c__1);

	    i__3 = *dp2 + 1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		e = k + i__ - 1;

		i__4 = *cp2;
		for (h__ = 1; h__ <= i__4; ++h__) {
		    p3[j + (h__ + e * p3_dim2) * p3_dim1] = ddot_(cp1, &dwork[
			    1], &c__1, &p2[(h__ + i__ * p2_dim2) * p2_dim1 + 
			    1], &c__1) + p3[j + (h__ + e * p3_dim2) * p3_dim1]
			    ;
/* L100: */
		}

/* L120: */
	    }

/* L140: */
	}

/* L160: */
    }

/*     Computation of the exact degree of P3(x). */

    cfzero = TRUE_;
/*     WHILE ( DP3 >= 0 and CFZERO ) DO */
L180:
    if (*dp3 >= 0 && cfzero) {
	dpol3 = *dp3 + 1;

	i__1 = *cp2;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *rp1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (p3[i__ + (j + dpol3 * p3_dim2) * p3_dim1] != 0.) {
		    cfzero = FALSE_;
		}
/* L200: */
	    }

/* L220: */
	}

	if (cfzero) {
	    --(*dp3);
	}
	goto L180;
    }
/*     END WHILE 180 */

    return 0;
/* *** Last line of MC03MD *** */
} /* mc03md_ */

