/* MB01ZD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb01zd_(char *side, char *uplo, char *transt, char *diag,
	 integer *m, integer *n, integer *l, doublereal *alpha, doublereal *t,
	 integer *ldt, doublereal *h__, integer *ldh, integer *info, ftnlen 
	side_len, ftnlen uplo_len, ftnlen transt_len, ftnlen diag_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, m2;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical trans, upper;
    static integer nrowt;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

/*     To compute the matrix product */

/*        H := alpha*op( T )*H,   or   H := alpha*H*op( T ), */

/*     where alpha is a scalar, H is an m-by-n upper or lower */
/*     Hessenberg-like matrix (with l nonzero subdiagonals or */
/*     superdiagonals, respectively), T is a unit, or non-unit, */
/*     upper or lower triangular matrix, and op( T ) is one of */

/*        op( T ) = T   or   op( T ) = T'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the triangular matrix T appears on the */
/*             left or right in the matrix product, as follows: */
/*             = 'L':  the product alpha*op( T )*H is computed; */
/*             = 'R':  the product alpha*H*op( T ) is computed. */

/*     UPLO    CHARACTER*1 */
/*             Specifies the form of the matrices T and H, as follows: */
/*             = 'U':  the matrix T is upper triangular and the matrix H */
/*                     is upper Hessenberg-like; */
/*             = 'L':  the matrix T is lower triangular and the matrix H */
/*                     is lower Hessenberg-like. */

/*     TRANST  CHARACTER*1 */
/*             Specifies the form of op( T ) to be used, as follows: */
/*             = 'N':  op( T ) = T; */
/*             = 'T':  op( T ) = T'; */
/*             = 'C':  op( T ) = T'. */

/*     DIAG    CHARACTER*1. */
/*             Specifies whether or not T is unit triangular, as follows: */
/*             = 'U':  the matrix T is assumed to be unit triangular; */
/*             = 'N':  the matrix T is not assumed to be unit triangular. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of H.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of H.  N >= 0. */

/*     L       (input) INTEGER */
/*             If UPLO = 'U', matrix H has L nonzero subdiagonals. */
/*             If UPLO = 'L', matrix H has L nonzero superdiagonals. */
/*             MAX(0,M-1) >= L >= 0, if UPLO = 'U'; */
/*             MAX(0,N-1) >= L >= 0, if UPLO = 'L'. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then T is not */
/*             referenced and H need not be set before entry. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,k), where */
/*             k is m when SIDE = 'L' and is n when SIDE = 'R'. */
/*             If UPLO = 'U', the leading k-by-k upper triangular part */
/*             of this array must contain the upper triangular matrix T */
/*             and the strictly lower triangular part is not referenced. */
/*             If UPLO = 'L', the leading k-by-k lower triangular part */
/*             of this array must contain the lower triangular matrix T */
/*             and the strictly upper triangular part is not referenced. */
/*             Note that when DIAG = 'U', the diagonal elements of T are */
/*             not referenced either, but are assumed to be unity. */

/*     LDT     INTEGER */
/*             The leading dimension of array T. */
/*             LDT >= MAX(1,M), if SIDE = 'L'; */
/*             LDT >= MAX(1,N), if SIDE = 'R'. */

/*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
/*             On entry, if UPLO = 'U', the leading M-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg-like matrix H. */
/*             On entry, if UPLO = 'L', the leading M-by-N lower */
/*             Hessenberg part of this array must contain the lower */
/*             Hessenberg-like matrix H. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix product alpha*op( T )*H, if SIDE = 'L', */
/*             or alpha*H*op( T ), if SIDE = 'R'. If TRANST = 'N', this */
/*             product has the same pattern as the given matrix H; */
/*             the elements below the L-th subdiagonal (if UPLO = 'U'), */
/*             or above the L-th superdiagonal (if UPLO = 'L'), are not */
/*             referenced in this case. If TRANST = 'T', the elements */
/*             below the (N+L)-th row (if UPLO = 'U', SIDE = 'R', and */
/*             M > N+L), or at the right of the (M+L)-th column */
/*             (if UPLO = 'L', SIDE = 'L', and N > M+L), are not set to */
/*             zero nor referenced. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The calculations are efficiently performed taking the problem */
/*     structure into account. */

/*     FURTHER COMMENTS */

/*     The matrix H may have the following patterns, when m = 7, n = 6, */
/*     and l = 2 are used for illustration: */

/*               UPLO = 'U'                    UPLO = 'L' */

/*            [ x x x x x x ]               [ x x x 0 0 0 ] */
/*            [ x x x x x x ]               [ x x x x 0 0 ] */
/*            [ x x x x x x ]               [ x x x x x 0 ] */
/*        H = [ 0 x x x x x ],          H = [ x x x x x x ]. */
/*            [ 0 0 x x x x ]               [ x x x x x x ] */
/*            [ 0 0 0 x x x ]               [ x x x x x x ] */
/*            [ 0 0 0 0 x x ]               [ x x x x x x ] */

/*     The products T*H or H*T have the same pattern as H, but the */
/*     products T'*H or H*T' may be full matrices. */

/*     If m = n, the matrix H is upper or lower triangular, for l = 0, */
/*     and upper or lower Hessenberg, for l = 1. */

/*     This routine is a specialization of the BLAS 3 routine DTRMM. */
/*     BLAS 1 calls are used when appropriate, instead of in-line code, */
/*     in order to increase the efficiency. If the matrix H is full, or */
/*     its zero triangle has small order, an optimized DTRMM code could */
/*     be faster than MB01ZD. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    /* Function Body */
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    trans = lsame_(transt, "T", (ftnlen)1, (ftnlen)1) || lsame_(transt, "C", (
	    ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
    if (lside) {
	nrowt = *m;
    } else {
	nrowt = *n;
    }

    if (upper) {
	m2 = *m;
    } else {
	m2 = *n;
    }

    *info = 0;
    if (! (lside || lsame_(side, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (trans || lsame_(transt, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (nounit || lsame_(diag, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 0, i__2 = m2 - 1;
	if (*l < 0 || *l > max(i__1,i__2)) {
	    *info = -7;
	} else if (*ldt < max(1,nrowt)) {
	    *info = -10;
	} else if (*ldh < max(1,*m)) {
	    *info = -12;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01ZD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return, if possible. */

    if (min(*m,*n) == 0) {
	return 0;
    }

/*     Also, when alpha = 0. */

    if (*alpha == 0.) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (upper) {
		i1 = 1;
/* Computing MIN */
		i__2 = j + *l;
		i2 = min(i__2,*m);
	    } else {
/* Computing MAX */
		i__2 = 1, i__3 = j - *l;
		i1 = max(i__2,i__3);
		i2 = *m;
	    }

	    i__2 = i2;
	    for (i__ = i1; i__ <= i__2; ++i__) {
		h__[i__ + j * h_dim1] = 0.;
/* L10: */
	    }

/* L20: */
	}

	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (! trans) {

/*           Form  H := alpha*T*H. */

	    if (upper) {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
		    i__3 = j + *l;
		    i__2 = min(i__3,*m);
		    for (k = 1; k <= i__2; ++k) {
			if (h__[k + j * h_dim1] != 0.) {
			    temp = *alpha * h__[k + j * h_dim1];
			    i__3 = k - 1;
			    daxpy_(&i__3, &temp, &t[k * t_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);
			    if (nounit) {
				temp *= t[k + k * t_dim1];
			    }
			    h__[k + j * h_dim1] = temp;
			}
/* L30: */
		    }

/* L40: */
		}

	    } else {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

/* Computing MAX */
		    i__3 = 1, i__4 = j - *l;
		    i__2 = max(i__3,i__4);
		    for (k = *m; k >= i__2; --k) {
			if (h__[k + j * h_dim1] != 0.) {
			    temp = *alpha * h__[k + j * h_dim1];
			    h__[k + j * h_dim1] = temp;
			    if (nounit) {
				h__[k + j * h_dim1] *= t[k + k * t_dim1];
			    }
			    i__3 = *m - k;
			    daxpy_(&i__3, &temp, &t[k + 1 + k * t_dim1], &
				    c__1, &h__[k + 1 + j * h_dim1], &c__1);
			}
/* L50: */
		    }

/* L60: */
		}

	    }

	} else {

/*           Form  H := alpha*T'*H. */

	    if (upper) {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i1 = j + *l;

		    for (i__ = *m; i__ >= 1; --i__) {
			if (i__ > i1) {
			    temp = ddot_(&i1, &t[i__ * t_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);
			} else {
			    temp = h__[i__ + j * h_dim1];
			    if (nounit) {
				temp *= t[i__ + i__ * t_dim1];
			    }
			    i__2 = i__ - 1;
			    temp += ddot_(&i__2, &t[i__ * t_dim1 + 1], &c__1, 
				    &h__[j * h_dim1 + 1], &c__1);
			}
			h__[i__ + j * h_dim1] = *alpha * temp;
/* L70: */
		    }

/* L80: */
		}

	    } else {

/* Computing MIN */
		i__2 = *m + *l;
		i__1 = min(i__2,*n);
		for (j = 1; j <= i__1; ++j) {
		    i1 = j - *l;

		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			if (i__ < i1) {
			    i__3 = *m - i1 + 1;
			    temp = ddot_(&i__3, &t[i1 + i__ * t_dim1], &c__1, 
				    &h__[i1 + j * h_dim1], &c__1);
			} else {
			    temp = h__[i__ + j * h_dim1];
			    if (nounit) {
				temp *= t[i__ + i__ * t_dim1];
			    }
			    i__3 = *m - i__;
			    temp += ddot_(&i__3, &t[i__ + 1 + i__ * t_dim1], &
				    c__1, &h__[i__ + 1 + j * h_dim1], &c__1);
			}
			h__[i__ + j * h_dim1] = *alpha * temp;
/* L90: */
		    }

/* L100: */
		}

	    }

	}

    } else {

	if (! trans) {

/*           Form  H := alpha*H*T. */

	    if (upper) {

		for (j = *n; j >= 1; --j) {
/* Computing MIN */
		    i__1 = j + *l;
		    i2 = min(i__1,*m);
		    temp = *alpha;
		    if (nounit) {
			temp *= t[j + j * t_dim1];
		    }
		    dscal_(&i2, &temp, &h__[j * h_dim1 + 1], &c__1);

		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			d__1 = *alpha * t[k + j * t_dim1];
			daxpy_(&i2, &d__1, &h__[k * h_dim1 + 1], &c__1, &h__[
				j * h_dim1 + 1], &c__1);
/* L110: */
		    }

/* L120: */
		}

	    } else {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		    i__2 = 1, i__3 = j - *l;
		    i1 = max(i__2,i__3);
		    temp = *alpha;
		    if (nounit) {
			temp *= t[j + j * t_dim1];
		    }
		    i__2 = *m - i1 + 1;
		    dscal_(&i__2, &temp, &h__[i1 + j * h_dim1], &c__1);

		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			i__3 = *m - i1 + 1;
			d__1 = *alpha * t[k + j * t_dim1];
			daxpy_(&i__3, &d__1, &h__[i1 + k * h_dim1], &c__1, &
				h__[i1 + j * h_dim1], &c__1);
/* L130: */
		    }

/* L140: */
		}

	    }

	} else {

/*           Form  H := alpha*H*T'. */

	    if (upper) {
/* Computing MIN */
		i__1 = *n + *l;
		m2 = min(i__1,*m);

		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
		    i__2 = k + *l;
		    i1 = min(i__2,*m);
/* Computing MIN */
		    i__2 = k + *l;
		    i2 = min(i__2,m2);

		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (t[j + k * t_dim1] != 0.) {
			    temp = *alpha * t[j + k * t_dim1];
			    daxpy_(&i1, &temp, &h__[k * h_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);

			    i__3 = i2;
			    for (i__ = i1 + 1; i__ <= i__3; ++i__) {
				h__[i__ + j * h_dim1] = temp * h__[i__ + k * 
					h_dim1];
/* L150: */
			    }

			}
/* L160: */
		    }

		    temp = *alpha;
		    if (nounit) {
			temp *= t[k + k * t_dim1];
		    }
		    if (temp != 1.) {
			dscal_(&i2, &temp, &h__[k * h_dim1 + 1], &c__1);
		    }
/* L170: */
		}

	    } else {

		for (k = *n; k >= 1; --k) {
/* Computing MAX */
		    i__1 = 1, i__2 = k - *l;
		    i1 = max(i__1,i__2);
/* Computing MAX */
		    i__1 = 1, i__2 = k - *l + 1;
		    i2 = max(i__1,i__2);
/* Computing MIN */
		    i__1 = *m, i__2 = i2 - 1;
		    m2 = min(i__1,i__2);

		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (t[j + k * t_dim1] != 0.) {
			    temp = *alpha * t[j + k * t_dim1];
			    i__2 = *m - i2 + 1;
			    daxpy_(&i__2, &temp, &h__[i2 + k * h_dim1], &c__1,
				     &h__[i2 + j * h_dim1], &c__1);

			    i__2 = m2;
			    for (i__ = i1; i__ <= i__2; ++i__) {
				h__[i__ + j * h_dim1] = temp * h__[i__ + k * 
					h_dim1];
/* L180: */
			    }

			}
/* L190: */
		    }

		    temp = *alpha;
		    if (nounit) {
			temp *= t[k + k * t_dim1];
		    }
		    if (temp != 1.) {
			i__1 = *m - i1 + 1;
			dscal_(&i__1, &temp, &h__[i1 + k * h_dim1], &c__1);
		    }
/* L200: */
		}

	    }

	}

    }

    return 0;

/* *** Last line of MB01ZD *** */
} /* mb01zd_ */

