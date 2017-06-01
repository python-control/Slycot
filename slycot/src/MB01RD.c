/* MB01RD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = .5;
static doublereal c_b11 = 0.;
static integer c__0 = 0;
static doublereal c_b15 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01rd_(char *uplo, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr, 
	doublereal *a, integer *lda, doublereal *x, integer *ldx, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen uplo_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, x_dim1, x_offset, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer j, ldw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char ntran[12];
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nrowa;
    static logical luplo;
    static integer jwork;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
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

/*     To compute the matrix formula */
/*        _ */
/*        R = alpha*R + beta*op( A )*X*op( A )', */
/*                                                 _ */
/*     where alpha and beta are scalars, R, X, and R are symmetric */
/*     matrices, A is a general matrix, and op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A'. */

/*     The result is overwritten on R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1                                         _ */
/*             Specifies which triangles of the symmetric matrices R, R, */
/*             and X are given as follows: */
/*             = 'U':  the upper triangular part is given; */
/*             = 'L':  the lower triangular part is given. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( A ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER           _ */
/*             The order of the matrices R and R and the number of rows */
/*             of the matrix op( A ).  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix X and the number of columns of the */
/*             the matrix op( A ).  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then R need not be */
/*             set before entry, except when R is identified with X in */
/*             the call (which is possible only in this case). */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then A and X are not */
/*             referenced. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry with UPLO = 'U', the leading M-by-M upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix R; the strictly */
/*             lower triangular part of the array is used as workspace. */
/*             On entry with UPLO = 'L', the leading M-by-M lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix R; the strictly */
/*             upper triangular part of the array is used as workspace. */
/*             On exit, the leading M-by-M upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of */
/*             this array contains the corresponding triangular part of */
/*                                 _ */
/*             the computed matrix R. If beta <> 0, the remaining */
/*             strictly triangular part of this array contains the */
/*             corresponding part of the matrix expression */
/*             beta*op( A )*T*op( A )', where T is the triangular matrix */
/*             defined in the Method section. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,k) */
/*             where k is N when TRANS = 'N' and is M when TRANS = 'T' or */
/*             TRANS = 'C'. */
/*             On entry with TRANS = 'N', the leading M-by-N part of this */
/*             array must contain the matrix A. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             N-by-M part of this array must contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,l), */
/*             where l is M when TRANS = 'N' and is N when TRANS = 'T' or */
/*             TRANS = 'C'. */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix X and the strictly */
/*             lower triangular part of the array is not referenced. */
/*             On entry, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix X and the strictly */
/*             upper triangular part of the array is not referenced. */
/*             On exit, each diagonal element of this array has half its */
/*             input value, but the other elements are not modified. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, the leading M-by-N part of this */
/*             array (with the leading dimension MAX(1,M)) returns the */
/*             matrix product beta*op( A )*T, where T is the triangular */
/*             matrix defined in the Method section. */
/*             This array is not referenced when beta = 0. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,M*N), if  beta <> 0; */
/*             LDWORK >= 1,          if  beta =  0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -k, the k-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix expression is efficiently evaluated taking the symmetry */
/*     into account. Specifically, let X = T + T', with T an upper or */
/*     lower triangular matrix, defined by */

/*        T = triu( X ) - (1/2)*diag( X ),  if UPLO = 'U', */
/*        T = tril( X ) - (1/2)*diag( X ),  if UPLO = 'L', */

/*     where triu, tril, and diag denote the upper triangular part, lower */
/*     triangular part, and diagonal part of X, respectively. Then, */

/*        op( A )*X*op( A )' = B + B', */

/*     where B := op( A )*T*op( A )'. Matrix B is not symmetric, but it */
/*     can be written as tri( B ) + stri( B ), where tri denotes the */
/*     triangular part specified by UPLO, and stri denotes the remaining */
/*     strictly triangular part. Let R = V + V', with V defined as T */
/*     above. Then, the required triangular part of the result can be */
/*     written as */

/*        alpha*V + beta*tri( B )  + beta*(stri( B ))' + */
/*                 alpha*diag( V ) + beta*diag( tri( B ) ). */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*                   2         2 */
/*        3/2 x M x N + 1/2 x M */

/*     operations. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004, */
/*     Apr. 2004. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

    if (ltrans) {
	nrowa = *n;
	s_copy(ntran, "No transpose", (ftnlen)12, (ftnlen)12);
    } else {
	nrowa = *m;
	s_copy(ntran, "Transpose", (ftnlen)12, (ftnlen)9);
    }

    ldw = max(1,*m);

    if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldr < ldw) {
	*info = -8;
    } else if (*lda < max(1,nrowa)) {
	*info = -10;
    } else if (*ldx < max(1,*n)) {
	*info = -12;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m * *n;
	if (*beta != 0. && *ldwork < max(i__1,i__2) || *beta == 0. && *ldwork 
		< 1) {
	    *info = -14;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01RD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    i__1 = *ldx + 1;
    dscal_(n, &c_b10, &x[x_offset], &i__1);
    if (*m == 0) {
	return 0;
    }

    if (*beta == 0. || *n == 0) {
	if (*alpha == 0.) {

/*           Special case alpha = 0. */

	    dlaset_(uplo, m, m, &c_b11, &c_b11, &r__[r_offset], ldr, (ftnlen)
		    1);
	} else {

/*           Special case beta = 0 or N = 0. */

	    if (*alpha != 1.) {
		dlascl_(uplo, &c__0, &c__0, &c_b15, alpha, m, m, &r__[
			r_offset], ldr, info, (ftnlen)1);
	    }
	}
	return 0;
    }

/*     General case: beta <> 0. Efficiently compute */
/*        _ */
/*        R = alpha*R + beta*op( A )*X*op( A )', */

/*     as described in the Method section. */

/*     Compute W = beta*op( A )*T in DWORK. */
/*     Workspace: need M*N. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code.) */

    if (ltrans) {
	jwork = 1;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(m, &a[j + a_dim1], lda, &dwork[jwork], &c__1);
	    jwork += ldw;
/* L10: */
	}

    } else {
	dlacpy_("Full", m, n, &a[a_offset], lda, &dwork[1], &ldw, (ftnlen)4);
    }

    dtrmm_("Right", uplo, "No transpose", "Non-unit", m, n, beta, &x[x_offset]
	    , ldx, &dwork[1], &ldw, (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)
	    8);

/*     Compute Y = alpha*V + W*op( A )' in R. First, set to zero the */
/*     strictly triangular part of R not specified by UPLO. That part */
/*     will then contain beta*stri( B ). */

    if (*alpha != 0.) {
	if (*m > 1) {
	    if (luplo) {
		i__1 = *m - 1;
		i__2 = *m - 1;
		dlaset_("Lower", &i__1, &i__2, &c_b11, &c_b11, &r__[r_dim1 + 
			2], ldr, (ftnlen)5);
	    } else {
		i__1 = *m - 1;
		i__2 = *m - 1;
		dlaset_("Upper", &i__1, &i__2, &c_b11, &c_b11, &r__[(r_dim1 <<
			 1) + 1], ldr, (ftnlen)5);
	    }
	}
	i__1 = *ldr + 1;
	dscal_(m, &c_b10, &r__[r_offset], &i__1);
    }

    dgemm_("No transpose", ntran, m, m, n, &c_b15, &dwork[1], &ldw, &a[
	    a_offset], lda, alpha, &r__[r_offset], ldr, (ftnlen)12, (ftnlen)
	    12);

/*     Add the term corresponding to B', with B = op( A )*T*op( A )'. */

    if (luplo) {

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    daxpy_(&j, &c_b15, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1], &
		    c__1);
/* L20: */
	}

    } else {

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    daxpy_(&j, &c_b15, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], 
		    ldr);
/* L30: */
	}

    }

    return 0;
/* *** Last line of MB01RD *** */
} /* mb01rd_ */

