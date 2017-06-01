/* MB04DI.f -- translated by f2c (version 20100827).
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

static doublereal c_b14 = -1.;

/* Subroutine */ int mb04di_(char *job, char *sgn, integer *n, integer *ilo, 
	doublereal *scale, integer *m, doublereal *v1, integer *ldv1, 
	doublereal *v2, integer *ldv2, integer *info, ftnlen job_len, ftnlen 
	sgn_len)
{
    /* System generated locals */
    integer v1_dim1, v1_offset, v2_dim1, v2_offset, i__1;

    /* Local variables */
    static integer i__, k;
    static logical lsgn, sysw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static logical lperm;
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

/*     To apply the inverse of a balancing transformation, computed by */
/*     the SLICOT Library routines MB04DD or MB04DS, to a 2*N-by-M matrix */

/*               [   V1   ] */
/*               [        ], */
/*               [ sgn*V2 ] */

/*     where sgn is either +1 or -1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the type of inverse transformation required: */
/*             = 'N':  do nothing, return immediately; */
/*             = 'P':  do inverse transformation for permutation only; */
/*             = 'S':  do inverse transformation for scaling only; */
/*             = 'B':  do inverse transformations for both permutation */
/*                     and scaling. */
/*             JOB must be the same as the argument JOB supplied to */
/*             MB04DD or MB04DS. */

/*     SGN     CHARACTER*1 */
/*             Specifies the sign to use for V2: */
/*             = 'P':  sgn = +1; */
/*             = 'N':  sgn = -1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrices V1 and V2. N >= 0. */

/*     ILO     (input) INTEGER */
/*             The integer ILO determined by MB04DD or MB04DS. */
/*             1 <= ILO <= N+1. */

/*     SCALE   (input) DOUBLE PRECISION array, dimension (N) */
/*             Details of the permutation and scaling factors, as */
/*             returned by MB04DD or MB04DS. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices V1 and V2.  M >= 0. */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix V1. */
/*             On exit, the leading N-by-M part of this array is */
/*             overwritten by the updated matrix V1 of the transformed */
/*             matrix. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1. LDV1 >= max(1,N). */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix V2. */
/*             On exit, the leading N-by-M part of this array is */
/*             overwritten by the updated matrix V2 of the transformed */
/*             matrix. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2. LDV2 >= max(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHABAK). */

/*     KEYWORDS */

/*     Balancing, Hamiltonian matrix, skew-Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --scale;
    v1_dim1 = *ldv1;
    v1_offset = 1 + v1_dim1;
    v1 -= v1_offset;
    v2_dim1 = *ldv2;
    v2_offset = 1 + v2_dim1;
    v2 -= v2_offset;

    /* Function Body */
    *info = 0;
    lperm = lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
    lscal = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
    lsgn = lsame_(sgn, "N", (ftnlen)1, (ftnlen)1);
    if (! lperm && ! lscal && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lsgn && ! lsame_(sgn, "P", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1 || *ilo > *n + 1) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*ldv1 < max(1,*n)) {
	*info = -8;
    } else if (*ldv2 < max(1,*n)) {
	*info = -10;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04DI", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
	return 0;
    }

/*     Inverse scaling. */

    if (lscal) {
	i__1 = *n;
	for (i__ = *ilo; i__ <= i__1; ++i__) {
	    drscl_(m, &scale[i__], &v1[i__ + v1_dim1], ldv1);
/* L20: */
	}
	i__1 = *n;
	for (i__ = *ilo; i__ <= i__1; ++i__) {
	    drscl_(m, &scale[i__], &v2[i__ + v2_dim1], ldv2);
/* L30: */
	}
    }

/*     Inverse permutation. */

    if (lperm) {
	for (i__ = *ilo - 1; i__ >= 1; --i__) {
	    k = (integer) scale[i__];
	    sysw = k > *n;
	    if (sysw) {
		k -= *n;
	    }

	    if (k != i__) {

/*              Exchange rows k <-> i. */

		dswap_(m, &v1[i__ + v1_dim1], ldv1, &v1[k + v1_dim1], ldv1);
		dswap_(m, &v2[i__ + v2_dim1], ldv2, &v2[k + v2_dim1], ldv2);
	    }

	    if (sysw) {

/*              Exchange V1(k,:) <-> V2(k,:). */

		dswap_(m, &v1[k + v1_dim1], ldv1, &v2[k + v2_dim1], ldv2);
		if (lsgn) {
		    dscal_(m, &c_b14, &v2[k + v2_dim1], ldv2);
		} else {
		    dscal_(m, &c_b14, &v1[k + v1_dim1], ldv1);
		}
	    }
/* L40: */
	}
    }

    return 0;
/* *** Last line of MB04DI *** */
} /* mb04di_ */

