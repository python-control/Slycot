/* MB01RW.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 1.;
static doublereal c_b13 = 0.;

/* Subroutine */ int mb01rw_(char *uplo, char *trans, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *z__, integer *ldz, 
	doublereal *dwork, integer *info, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nottra;


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

/*     To compute the transformation of the symmetric matrix A by the */
/*     matrix Z in the form */

/*        A := op(Z)*A*op(Z)', */

/*     where op(Z) is either Z or its transpose, Z'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the upper or lower triangle of A */
/*             is stored: */
/*             = 'U':  Upper triangle of A is stored; */
/*             = 'L':  Lower triangle of A is stored. */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether op(Z) is Z or its transpose Z': */
/*             = 'N':  op(Z) = Z; */
/*             = 'T':  op(Z) = Z'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the resulting symmetric matrix op(Z)*A*op(Z)' */
/*             and the number of rows of the matrix Z, if TRANS = 'N', */
/*             or the number of columns of the matrix Z, if TRANS = 'T'. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the symmetric matrix A and the number of */
/*             columns of the matrix Z, if TRANS = 'N', or the number of */
/*             rows of the matrix Z, if TRANS = 'T'.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,MAX(M,N)) */
/*             On entry, the leading N-by-N upper or lower triangular */
/*             part of this array must contain the upper (UPLO = 'U') */
/*             or lower (UPLO = 'L') triangular part of the symmetric */
/*             matrix A. */
/*             On exit, the leading M-by-M upper or lower triangular */
/*             part of this array contains the upper (UPLO = 'U') or */
/*             lower (UPLO = 'L') triangular part of the symmetric */
/*             matrix op(Z)*A*op(Z)'. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M,N). */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDQ,K) */
/*             where K = N if TRANS = 'N' and K = M if TRANS = 'T'. */
/*             The leading M-by-N part, if TRANS = 'N', or N-by-M part, */
/*             if TRANS = 'T', of this array contains the matrix Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z. */
/*             LDZ >= MAX(1,M) if TRANS = 'N' and */
/*             LDZ >= MAX(1,N) if TRANS = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     FURTHER COMMENTS */

/*     This is a simpler, BLAS 2 version for MB01RD. */

/*     CONTRIBUTOR */

/*     A. Varga, DLR, Feb. 1995. */

/*     REVISIONS */

/*     April 1998 (T. Penzl). */
/*     Sep. 1998 (V. Sima). */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --dwork;

    /* Function Body */
    nottra = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (nottra || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*lda < max(i__1,*n)) {
	    *info = -6;
	} else if (nottra && *ldz < max(1,*m) || ! nottra && *ldz < max(1,*n))
		 {
	    *info = -8;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB01RW", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0) {
	return 0;
    }

    if (nottra) {

/*        Compute Z*A*Z'. */

	if (upper) {

/*           Compute Z*A in A (M-by-N). */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy_(&i__2, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
		i__2 = *n - j + 1;
		dcopy_(&i__2, &a[j + j * a_dim1], lda, &dwork[j], &c__1);
		dgemv_(trans, m, n, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[j * a_dim1 + 1], &c__1, (ftnlen)1);
/* L10: */
	    }

/*           Compute A*Z' in the upper triangular part of A. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
		i__2 = *m - i__ + 1;
		dgemv_(trans, &i__2, n, &c_b11, &z__[i__ + z_dim1], ldz, &
			dwork[1], &c__1, &c_b13, &a[i__ + i__ * a_dim1], lda, 
			(ftnlen)1);
/* L20: */
	    }

	} else {

/*           Compute A*Z' in A (N-by-M). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__ - 1;
		dcopy_(&i__2, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
		i__2 = *n - i__ + 1;
		dcopy_(&i__2, &a[i__ + i__ * a_dim1], &c__1, &dwork[i__], &
			c__1);
		dgemv_(trans, m, n, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[i__ + a_dim1], lda, (ftnlen)1);
/* L30: */
	    }

/*           Compute Z*A in the lower triangular part of A. */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
		i__2 = *m - j + 1;
		dgemv_(trans, &i__2, n, &c_b11, &z__[j + z_dim1], ldz, &dwork[
			1], &c__1, &c_b13, &a[j + j * a_dim1], &c__1, (ftnlen)
			1);
/* L40: */
	    }

	}
    } else {

/*        Compute Z'*A*Z. */

	if (upper) {

/*           Compute Z'*A in A (M-by-N). */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy_(&i__2, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
		i__2 = *n - j + 1;
		dcopy_(&i__2, &a[j + j * a_dim1], lda, &dwork[j], &c__1);
		dgemv_(trans, n, m, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[j * a_dim1 + 1], &c__1, (ftnlen)1);
/* L50: */
	    }

/*           Compute A*Z in the upper triangular part of A. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
		i__2 = *m - i__ + 1;
		dgemv_(trans, n, &i__2, &c_b11, &z__[i__ * z_dim1 + 1], ldz, &
			dwork[1], &c__1, &c_b13, &a[i__ + i__ * a_dim1], lda, 
			(ftnlen)1);
/* L60: */
	    }

	} else {

/*           Compute A*Z in A (N-by-M). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__ - 1;
		dcopy_(&i__2, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
		i__2 = *n - i__ + 1;
		dcopy_(&i__2, &a[i__ + i__ * a_dim1], &c__1, &dwork[i__], &
			c__1);
		dgemv_(trans, n, m, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[i__ + a_dim1], lda, (ftnlen)1);
/* L70: */
	    }

/*           Compute Z'*A in the lower triangular part of A. */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
		i__2 = *m - j + 1;
		dgemv_(trans, n, &i__2, &c_b11, &z__[j * z_dim1 + 1], ldz, &
			dwork[1], &c__1, &c_b13, &a[j + j * a_dim1], &c__1, (
			ftnlen)1);
/* L80: */
	    }

	}
    }

    return 0;
/* *** Last line of MB01RW *** */
} /* mb01rw_ */

