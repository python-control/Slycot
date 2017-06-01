/* MB01RY.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static integer c__0 = 0;
static doublereal c_b14 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01ry_(char *side, char *uplo, char *trans, integer *m, 
	doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr, 
	doublereal *h__, integer *ldh, doublereal *b, integer *ldb, 
	doublereal *dwork, integer *info, ftnlen side_len, ftnlen uplo_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, r_dim1, r_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    static logical luplo;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

/*     To compute either the upper or lower triangular part of one of the */
/*     matrix formulas */
/*        _ */
/*        R = alpha*R + beta*op( H )*B,                               (1) */
/*        _ */
/*        R = alpha*R + beta*B*op( H ),                               (2) */
/*                                                    _ */
/*     where alpha and beta are scalars, H, B, R, and R are m-by-m */
/*     matrices, H is an upper Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     The result is overwritten on R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*                     _ */
/*             = 'L':  R = alpha*R + beta*op( H )*B; */
/*                     _ */
/*             = 'R':  R = alpha*R + beta*B*op( H ). */

/*     UPLO    CHARACTER*1                               _ */
/*             Specifies which triangles of the matrices R and R are */
/*             computed and given, respectively, as follows: */
/*             = 'U':  the upper triangular part; */
/*             = 'L':  the lower triangular part. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER           _ */
/*             The order of the matrices R, R, H and B.  M >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then R need not be */
/*             set before entry. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then H and B are not */
/*             referenced. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry with UPLO = 'U', the leading M-by-M upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the matrix R; the strictly lower */
/*             triangular part of the array is not referenced. */
/*             On entry with UPLO = 'L', the leading M-by-M lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the matrix R; the strictly upper */
/*             triangular part of the array is not referenced. */
/*             On exit, the leading M-by-M upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L') of */
/*             this array contains the corresponding triangular part of */
/*                                 _ */
/*             the computed matrix R. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,M) */
/*             On entry, the leading M-by-M upper Hessenberg part of */
/*             this array must contain the upper Hessenberg part of the */
/*             matrix H. */
/*             The elements below the subdiagonal are not referenced, */
/*             except possibly for those in the first column, which */
/*             could be overwritten, but are restored on exit. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             LDWORK >= M, if  beta <> 0 and SIDE = 'L'; */
/*             LDWORK >= 0, if  beta =  0 or  SIDE = 'R'. */
/*             This array is not referenced when beta = 0 or SIDE = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix expression is efficiently evaluated taking the */
/*     Hessenberg/triangular structure into account. BLAS 2 operations */
/*     are used. A block algorithm can be constructed; it can use BLAS 3 */
/*     GEMM operations for most computations, and calls of this BLAS 2 */
/*     algorithm for computing the triangles. */

/*     FURTHER COMMENTS */

/*     The main application of this routine is when the result should */
/*     be a symmetric matrix, e.g., when B = X*op( H )', for (1), or */
/*     B = op( H )'*X, for (2), where B is already available and X = X'. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations. */

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
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*ldr < max(1,*m)) {
	*info = -8;
    } else if (*ldh < max(1,*m)) {
	*info = -10;
    } else if (*ldb < max(1,*m)) {
	*info = -12;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01RY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0) {
	return 0;
    }

    if (*beta == 0.) {
	if (*alpha == 0.) {

/*           Special case when both alpha = 0 and beta = 0. */

	    dlaset_(uplo, m, m, &c_b10, &c_b10, &r__[r_offset], ldr, (ftnlen)
		    1);
	} else {

/*           Special case beta = 0. */

	    if (*alpha != 1.) {
		dlascl_(uplo, &c__0, &c__0, &c_b14, alpha, m, m, &r__[
			r_offset], ldr, info, (ftnlen)1);
	    }
	}
	return 0;
    }

/*     General case: beta <> 0. */
/*     Compute the required triangle of (1) or (2) using BLAS 2 */
/*     operations. */

    if (lside) {

/*        To avoid repeated references to the subdiagonal elements of H, */
/*        these are swapped with the corresponding elements of H in the */
/*        first column, and are finally restored. */

	if (*m > 2) {
	    i__1 = *m - 2;
	    i__2 = *ldh + 1;
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
	}

	if (luplo) {
	    if (ltrans) {

		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {

/*                 Multiply the transposed upper triangle of the leading */
/*                 j-by-j submatrix of H by the leading part of the j-th */
/*                 column of B. */

		    dcopy_(&j, &b[j * b_dim1 + 1], &c__1, &dwork[1], &c__1);
		    dtrmv_("Upper", trans, "Non-unit", &j, &h__[h_offset], 
			    ldh, &dwork[1], &c__1, (ftnlen)5, (ftnlen)1, (
			    ftnlen)8);

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

/* Computing MIN */
		    i__3 = j, i__4 = *m - 1;
		    i__2 = min(i__3,i__4);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			r__[i__ + j * r_dim1] = *alpha * r__[i__ + j * r_dim1]
				 + *beta * (dwork[i__] + h__[i__ + 1 + h_dim1]
				 * b[i__ + 1 + j * b_dim1]);
/* L10: */
		    }

/* L20: */
		}

		r__[*m + *m * r_dim1] = *alpha * r__[*m + *m * r_dim1] + *
			beta * dwork[*m];

	    } else {

		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {

/*                 Multiply the upper triangle of the leading j-by-j */
/*                 submatrix of H by the leading part of the j-th column */
/*                 of B. */

		    dcopy_(&j, &b[j * b_dim1 + 1], &c__1, &dwork[1], &c__1);
		    dtrmv_("Upper", trans, "Non-unit", &j, &h__[h_offset], 
			    ldh, &dwork[1], &c__1, (ftnlen)5, (ftnlen)1, (
			    ftnlen)8);
		    if (j < *m) {

/*                    Multiply the remaining right part of the leading */
/*                    j-by-M submatrix of H by the trailing part of the */
/*                    j-th column of B. */

			i__2 = *m - j;
			dgemv_(trans, &j, &i__2, beta, &h__[(j + 1) * h_dim1 
				+ 1], ldh, &b[j + 1 + j * b_dim1], &c__1, 
				alpha, &r__[j * r_dim1 + 1], &c__1, (ftnlen)1)
				;
		    } else {
			dscal_(m, alpha, &r__[*m * r_dim1 + 1], &c__1);
		    }

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

		    r__[j * r_dim1 + 1] += *beta * dwork[1];

		    i__2 = j;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			r__[i__ + j * r_dim1] += *beta * (dwork[i__] + h__[
				i__ + h_dim1] * b[i__ - 1 + j * b_dim1]);
/* L30: */
		    }

/* L40: */
		}

	    }

	} else {

	    if (ltrans) {

		for (j = *m; j >= 1; --j) {

/*                 Multiply the transposed upper triangle of the trailing */
/*                 (M-j+1)-by-(M-j+1) submatrix of H by the trailing part */
/*                 of the j-th column of B. */

		    i__1 = *m - j + 1;
		    dcopy_(&i__1, &b[j + j * b_dim1], &c__1, &dwork[j], &c__1)
			    ;
		    i__1 = *m - j + 1;
		    dtrmv_("Upper", trans, "Non-unit", &i__1, &h__[j + j * 
			    h_dim1], ldh, &dwork[j], &c__1, (ftnlen)5, (
			    ftnlen)1, (ftnlen)8);
		    if (j > 1) {

/*                    Multiply the remaining left part of the trailing */
/*                    (M-j+1)-by-(j-1) submatrix of H' by the leading */
/*                    part of the j-th column of B. */

			i__1 = j - 1;
			i__2 = *m - j + 1;
			dgemv_(trans, &i__1, &i__2, beta, &h__[j * h_dim1 + 1]
				, ldh, &b[j * b_dim1 + 1], &c__1, alpha, &r__[
				j + j * r_dim1], &c__1, (ftnlen)1);
		    } else {
			dscal_(m, alpha, &r__[r_dim1 + 1], &c__1);
		    }

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

		    i__1 = *m - 1;
		    for (i__ = j; i__ <= i__1; ++i__) {
			r__[i__ + j * r_dim1] += *beta * (dwork[i__] + h__[
				i__ + 1 + h_dim1] * b[i__ + 1 + j * b_dim1]);
/* L50: */
		    }

		    r__[*m + j * r_dim1] += *beta * dwork[*m];
/* L60: */
		}

	    } else {

		for (j = *m; j >= 1; --j) {

/*                 Multiply the upper triangle of the trailing */
/*                 (M-j+1)-by-(M-j+1) submatrix of H by the trailing */
/*                 part of the j-th column of B. */

		    i__1 = *m - j + 1;
		    dcopy_(&i__1, &b[j + j * b_dim1], &c__1, &dwork[j], &c__1)
			    ;
		    i__1 = *m - j + 1;
		    dtrmv_("Upper", trans, "Non-unit", &i__1, &h__[j + j * 
			    h_dim1], ldh, &dwork[j], &c__1, (ftnlen)5, (
			    ftnlen)1, (ftnlen)8);

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

		    i__1 = *m;
		    for (i__ = max(j,2); i__ <= i__1; ++i__) {
			r__[i__ + j * r_dim1] = *alpha * r__[i__ + j * r_dim1]
				 + *beta * (dwork[i__] + h__[i__ + h_dim1] * 
				b[i__ - 1 + j * b_dim1]);
/* L70: */
		    }

/* L80: */
		}

		r__[r_dim1 + 1] = *alpha * r__[r_dim1 + 1] + *beta * dwork[1];

	    }
	}

	if (*m > 2) {
	    i__1 = *m - 2;
	    i__2 = *ldh + 1;
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
	}

    } else {

/*        Row-wise calculations are used for H, if SIDE = 'R' and */
/*        TRANS = 'T'. */

	if (luplo) {
	    if (ltrans) {
		r__[r_dim1 + 1] = *alpha * r__[r_dim1 + 1] + *beta * ddot_(m, 
			&b[b_offset], ldb, &h__[h_offset], ldh);

		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = *m - j + 2;
		    dgemv_("NoTranspose", &j, &i__2, beta, &b[(j - 1) * 
			    b_dim1 + 1], ldb, &h__[j + (j - 1) * h_dim1], ldh,
			     alpha, &r__[j * r_dim1 + 1], &c__1, (ftnlen)11);
/* L90: */
		}

	    } else {

		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j + 1;
		    dgemv_("NoTranspose", &j, &i__2, beta, &b[b_offset], ldb, 
			    &h__[j * h_dim1 + 1], &c__1, alpha, &r__[j * 
			    r_dim1 + 1], &c__1, (ftnlen)11);
/* L100: */
		}

		dgemv_("NoTranspose", m, m, beta, &b[b_offset], ldb, &h__[*m *
			 h_dim1 + 1], &c__1, alpha, &r__[*m * r_dim1 + 1], &
			c__1, (ftnlen)11);

	    }

	} else {

	    if (ltrans) {

		dgemv_("NoTranspose", m, m, beta, &b[b_offset], ldb, &h__[
			h_offset], ldh, alpha, &r__[r_dim1 + 1], &c__1, (
			ftnlen)11);

		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = *m - j + 1;
		    i__3 = *m - j + 2;
		    dgemv_("NoTranspose", &i__2, &i__3, beta, &b[j + (j - 1) *
			     b_dim1], ldb, &h__[j + (j - 1) * h_dim1], ldh, 
			    alpha, &r__[j + j * r_dim1], &c__1, (ftnlen)11);
/* L110: */
		}

	    } else {

		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m - j + 1;
		    i__3 = j + 1;
		    dgemv_("NoTranspose", &i__2, &i__3, beta, &b[j + b_dim1], 
			    ldb, &h__[j * h_dim1 + 1], &c__1, alpha, &r__[j + 
			    j * r_dim1], &c__1, (ftnlen)11);
/* L120: */
		}

		r__[*m + *m * r_dim1] = *alpha * r__[*m + *m * r_dim1] + *
			beta * ddot_(m, &b[*m + b_dim1], ldb, &h__[*m * 
			h_dim1 + 1], &c__1);

	    }
	}
    }

    return 0;
/* *** Last line of MB01RY *** */
} /* mb01ry_ */

