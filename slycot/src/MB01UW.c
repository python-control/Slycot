/* MB01UW.f -- translated by f2c (version 20100827).
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
static doublereal c_b10 = 0.;
static doublereal c_b33 = 1.;
static integer c__0 = 0;

/* Subroutine */ int mb01uw_(char *side, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *h__, integer *ldh, doublereal *a, 
	integer *lda, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, jw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrmv_(char *, char *, char *, integer *, doublereal *
	    , integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical ltrans;


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

/*     To compute one of the matrix products */

/*        A : = alpha*op( H ) * A, or A : = alpha*A * op( H ), */

/*     where alpha is a scalar, A is an m-by-n matrix, H is an upper */
/*     Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*             = 'L':  A := alpha*op( H ) * A; */
/*             = 'R':  A := alpha*A * op( H ). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then H is not */
/*             referenced and A need not be set before entry. */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,k) */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */
/*             On entry with SIDE = 'L', the leading M-by-M upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             On entry with SIDE = 'R', the leading N-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             The elements below the subdiagonal are not referenced, */
/*             except possibly for those in the first column, which */
/*             could be overwritten, but are restored on exit. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,k), */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the computed product. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, alpha <> 0, and LDWORK >= M*N > 0, */
/*             DWORK contains a copy of the matrix A, having the leading */
/*             dimension M. */
/*             This array is not referenced when alpha = 0. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= 0,   if  alpha =  0 or MIN(M,N) = 0; */
/*             LDWORK >= M-1, if  SIDE  = 'L'; */
/*             LDWORK >= N-1, if  SIDE  = 'R'. */
/*             For maximal efficiency LDWORK should be at least M*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The required matrix product is computed in two steps. In the first */
/*     step, the upper triangle of H is used; in the second step, the */
/*     contribution of the subdiagonal is added. If the workspace can */
/*     accomodate a copy of A, a fast BLAS 3 DTRMM operation is used in */
/*     the first step. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, January 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldh < 1 || lside && *ldh < *m || ! lside && *ldh < *n) {
	*info = -7;
    } else if (*lda < max(1,*m)) {
	*info = -9;
    } else if (*ldwork < 0 || *alpha != 0. && min(*m,*n) > 0 && (lside && *
	    ldwork < *m - 1 || ! lside && *ldwork < *n - 1)) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01UW", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return, if possible. */

    if (min(*m,*n) == 0) {
	return 0;
    } else if (lside) {
	if (*m == 1) {
	    d__1 = *alpha * h__[h_dim1 + 1];
	    dscal_(n, &d__1, &a[a_offset], lda);
	    return 0;
	}
    } else {
	if (*n == 1) {
	    d__1 = *alpha * h__[h_dim1 + 1];
	    dscal_(m, &d__1, &a[a_offset], &c__1);
	    return 0;
	}
    }

    if (*alpha == 0.) {

/*        Set A to zero and return. */

	dlaset_("Full", m, n, &c_b10, &c_b10, &a[a_offset], lda, (ftnlen)4);
	return 0;
    }

    if (*ldwork >= *m * *n) {

/*        Enough workspace for a fast BLAS 3 calculation. */
/*        Save A in the workspace and compute one of the matrix products */
/*          A : = alpha*op( triu( H ) ) * A, or */
/*          A : = alpha*A * op( triu( H ) ), */
/*        involving the upper triangle of H. */

	dlacpy_("Full", m, n, &a[a_offset], lda, &dwork[1], m, (ftnlen)4);
	dtrmm_(side, "Upper", trans, "Non-unit", m, n, alpha, &h__[h_offset], 
		ldh, &a[a_offset], lda, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);

/*        Add the contribution of the subdiagonal of H. */
/*        If SIDE = 'L', the subdiagonal of H is swapped with the */
/*        corresponding elements in the first column of H, and the */
/*        calculations are organized for column operations. */

	if (lside) {
	    if (*m > 2) {
		i__1 = *m - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }
	    if (ltrans) {
		jw = 1;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    ++jw;
		    i__2 = *m - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += *alpha * h__[i__ + 1 + h_dim1] 
				* dwork[jw];
			++jw;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		jw = 0;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    ++jw;
		    i__2 = *m;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += *alpha * h__[i__ + h_dim1] * 
				dwork[jw];
			++jw;
/* L30: */
		    }
/* L40: */
		}
	    }
	    if (*m > 2) {
		i__1 = *m - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }

	} else {

	    if (ltrans) {
		jw = 1;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    if (h__[j + 1 + j * h_dim1] != 0.) {
			d__1 = *alpha * h__[j + 1 + j * h_dim1];
			daxpy_(m, &d__1, &dwork[jw], &c__1, &a[(j + 1) * 
				a_dim1 + 1], &c__1);
		    }
		    jw += *m;
/* L50: */
		}
	    } else {
		jw = *m + 1;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    if (h__[j + 1 + j * h_dim1] != 0.) {
			d__1 = *alpha * h__[j + 1 + j * h_dim1];
			daxpy_(m, &d__1, &dwork[jw], &c__1, &a[j * a_dim1 + 1]
				, &c__1);
		    }
		    jw += *m;
/* L60: */
		}
	    }
	}

    } else {

/*        Use a BLAS 2 calculation. */

	if (lside) {
	    if (*m > 2) {
		i__1 = *m - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }
	    if (ltrans) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

		    i__2 = *m - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			dwork[i__] = h__[i__ + 1 + h_dim1] * a[i__ + 1 + j * 
				a_dim1];
/* L70: */
		    }

/*                 Multiply the upper triangle of H by the j-th column */
/*                 of A, and add to the above result. */

		    dtrmv_("Upper", trans, "Non-unit", m, &h__[h_offset], ldh,
			     &a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)1, 
			    (ftnlen)8);
		    i__2 = *m - 1;
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[j * a_dim1 + 1]
			    , &c__1);
/* L80: */
		}

	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

		    i__2 = *m - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			dwork[i__] = h__[i__ + 1 + h_dim1] * a[i__ + j * 
				a_dim1];
/* L90: */
		    }

/*                 Multiply the upper triangle of H by the j-th column */
/*                 of A, and add to the above result. */

		    dtrmv_("Upper", trans, "Non-unit", m, &h__[h_offset], ldh,
			     &a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)1, 
			    (ftnlen)8);
		    i__2 = *m - 1;
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[j * a_dim1 + 2]
			    , &c__1);
/* L100: */
		}
	    }
	    if (*m > 2) {
		i__1 = *m - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }

	} else {

/*           Below, row-wise calculations are used for A. */

	    if (*n > 2) {
		i__1 = *n - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }
	    if (ltrans) {
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the i-th row of the product. */

		    i__2 = *n - 1;
		    for (j = 1; j <= i__2; ++j) {
			dwork[j] = a[i__ + j * a_dim1] * h__[j + 1 + h_dim1];
/* L110: */
		    }

/*                 Multiply the i-th row of A by the upper triangle of H, */
/*                 and add to the above result. */

		    dtrmv_("Upper", "NoTranspose", "Non-unit", n, &h__[
			    h_offset], ldh, &a[i__ + a_dim1], lda, (ftnlen)5, 
			    (ftnlen)11, (ftnlen)8);
		    i__2 = *n - 1;
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[i__ + (a_dim1 
			    << 1)], lda);
/* L120: */
		}

	    } else {
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the i-th row of the product. */

		    i__2 = *n - 1;
		    for (j = 1; j <= i__2; ++j) {
			dwork[j] = a[i__ + (j + 1) * a_dim1] * h__[j + 1 + 
				h_dim1];
/* L130: */
		    }

/*                 Multiply the i-th row of A by the upper triangle of H, */
/*                 and add to the above result. */

		    dtrmv_("Upper", "Transpose", "Non-unit", n, &h__[h_offset]
			    , ldh, &a[i__ + a_dim1], lda, (ftnlen)5, (ftnlen)
			    9, (ftnlen)8);
		    i__2 = *n - 1;
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[i__ + a_dim1], 
			    lda);
/* L140: */
		}
	    }
	    if (*n > 2) {
		i__1 = *n - 2;
		i__2 = *ldh + 1;
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
	    }

	}

/*        Scale the result by alpha. */

	if (*alpha != 1.) {
	    dlascl_("General", &c__0, &c__0, &c_b33, alpha, m, n, &a[a_offset]
		    , lda, info, (ftnlen)7);
	}
    }
    return 0;
/* *** Last line of MB01UW *** */
} /* mb01uw_ */

