/* MB04ND.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04nd_(char *uplo, integer *n, integer *m, integer *p, 
	doublereal *r__, integer *ldr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *tau, doublereal *dwork, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, r_dim1, 
	    r_offset, i__1;

    /* Local variables */
    static integer i__, im, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04ny_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *);
    static logical luplo;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);


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

/*     To calculate an RQ factorization of the first block row and */
/*     apply the orthogonal transformations (from the right) also to the */
/*     second block row of a structured matrix, as follows */
/*                              _ */
/*       [ A   R ]        [ 0   R ] */
/*       [       ] * Q' = [ _   _ ] */
/*       [ C   B ]        [ C   B ] */
/*                 _ */
/*     where R and R are upper triangular. The matrix A can be full or */
/*     upper trapezoidal/triangular. The problem structure is exploited. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Indicates if the matrix A is or not triangular as follows: */
/*             = 'U':  Matrix A is upper trapezoidal/triangular; */
/*             = 'F':  Matrix A is full. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER                 _ */
/*             The order of the matrices R and R.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices B and C.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of columns of the matrices A and C.  P >= 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*                                                        _ */
/*             array contains the upper triangular matrix R. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,P) */
/*             On entry, if UPLO = 'F', the leading N-by-P part of this */
/*             array must contain the matrix A. For UPLO = 'U', if */
/*             N <= P, the upper triangle of the subarray A(1:N,P-N+1:P) */
/*             must contain the N-by-N upper triangular matrix A, and if */
/*             N >= P, the elements on and above the (N-P)-th subdiagonal */
/*             must contain the N-by-P upper trapezoidal matrix A. */
/*             On exit, if UPLO = 'F', the leading N-by-P part of this */
/*             array contains the trailing components (the vectors v, see */
/*             METHOD) of the elementary reflectors used in the */
/*             factorization. If UPLO = 'U', the upper triangle of the */
/*             subarray A(1:N,P-N+1:P) (if N <= P), or the elements on */
/*             and above the (N-P)-th subdiagonal (if N >= P), contain */
/*             the trailing components (the vectors v, see METHOD) of the */
/*             elementary reflectors used in the factorization. */
/*             The remaining elements are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading M-by-N part of this array contains */
/*                                 _ */
/*             the computed matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,P) */
/*             On entry, the leading M-by-P part of this array must */
/*             contain the matrix C. */
/*             On exit, the leading M-by-P part of this array contains */
/*                                 _ */
/*             the computed matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(N-1,M)) */

/*     METHOD */

/*     The routine uses N Householder transformations exploiting the zero */
/*     pattern of the block matrix.  A Householder matrix has the form */

/*                                     ( 1 ) */
/*        H  = I - tau *u *u',    u  = ( v ), */
/*         i          i  i  i      i   (  i) */

/*     where v  is a P-vector, if UPLO = 'F', or a min(N-i+1,P)-vector, */
/*            i */
/*     if UPLO = 'U'.  The components of v  are stored in the i-th row */
/*                                        i */
/*     of A, and tau  is stored in TAU(i), i = N,N-1,...,1. */
/*                  i */
/*     In-line code for applying Householder transformations is used */
/*     whenever possible (see MB04NY routine). */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1998. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary reflector, RQ factorization, orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --tau;
    --dwork;

    /* Function Body */
    if (min(*n,*p) == 0) {
	return 0;
    }

    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    if (luplo) {

	for (i__ = *n; i__ >= 1; --i__) {

/*           Annihilate the I-th row of A and apply the transformations */
/*           to the entire block matrix, exploiting its structure. */

/* Computing MIN */
	    i__1 = *n - i__ + 1;
	    im = min(i__1,*p);
/* Computing MAX */
	    i__1 = *p - *n + i__;
	    ip = max(i__1,1);
	    i__1 = im + 1;
	    dlarfg_(&i__1, &r__[i__ + i__ * r_dim1], &a[i__ + ip * a_dim1], 
		    lda, &tau[i__]);

/*           Compute */
/*                                                [ 1 ] */
/*           w := [ R(1:I-1,I)  A(1:I-1,IP:P) ] * [   ], */
/*                                                [ v ] */

/*           [ R(1:I-1,I)  A(1:I-1,IP:P) ] = */
/*           [ R(1:I-1,I)  A(1:I-1,IP:P) ] - tau * w * [ 1 v' ]. */

	    if (i__ > 0) {
		i__1 = i__ - 1;
		mb04ny_(&i__1, &im, &a[i__ + ip * a_dim1], lda, &tau[i__], &
			r__[i__ * r_dim1 + 1], ldr, &a[ip * a_dim1 + 1], lda, 
			&dwork[1]);
	    }


/*           Compute */
/*                                        [ 1 ] */
/*           w := [ B(:,I)  C(:,IP:P) ] * [   ], */
/*                                        [ v ] */

/*           [ B(:,I)  C(:,IP:P) ] = [ B(:,I)  C(:,IP:P) ] - */
/*                                   tau * w * [ 1 v' ]. */

	    if (*m > 0) {
		mb04ny_(m, &im, &a[i__ + ip * a_dim1], lda, &tau[i__], &b[i__ 
			* b_dim1 + 1], ldb, &c__[ip * c_dim1 + 1], ldc, &
			dwork[1]);
	    }
/* L10: */
	}

    } else {

	for (i__ = *n; i__ >= 2; --i__) {

/*           Annihilate the I-th row of A and apply the transformations */
/*           to the first block row, exploiting its structure. */

	    i__1 = *p + 1;
	    dlarfg_(&i__1, &r__[i__ + i__ * r_dim1], &a[i__ + a_dim1], lda, &
		    tau[i__]);

/*           Compute */
/*                                             [ 1 ] */
/*           w := [ R(1:I-1,I)  A(1:I-1,:) ] * [   ], */
/*                                             [ v ] */

/*           [ R(1:I-1,I)  A(1:I-1,:) ] = [ R(1:I-1,I)  A(1:I-1,:) ] - */
/*                                        tau * w * [ 1 v' ]. */

	    i__1 = i__ - 1;
	    mb04ny_(&i__1, p, &a[i__ + a_dim1], lda, &tau[i__], &r__[i__ * 
		    r_dim1 + 1], ldr, &a[a_offset], lda, &dwork[1]);
/* L20: */
	}

	i__1 = *p + 1;
	dlarfg_(&i__1, &r__[r_dim1 + 1], &a[a_dim1 + 1], lda, &tau[1]);
	if (*m > 0) {

/*           Apply the transformations to the second block row. */

	    for (i__ = *n; i__ >= 1; --i__) {

/*              Compute */
/*                                   [ 1 ] */
/*              w := [ B(:,I)  C ] * [   ], */
/*                                   [ v ] */

/*              [ B(:,I)  C ] = [ B(:,I)  C ] - tau * w * [ 1 v' ]. */

		mb04ny_(m, p, &a[i__ + a_dim1], lda, &tau[i__], &b[i__ * 
			b_dim1 + 1], ldb, &c__[c_offset], ldc, &dwork[1]);
/* L30: */
	    }

	}
    }
    return 0;
/* *** Last line of MB04ND *** */
} /* mb04nd_ */

