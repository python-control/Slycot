/* SB04PY.f -- translated by f2c (version 20100827).
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
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b28 = 1.;
static doublereal c_b31 = 0.;
static logical c_true = TRUE_;

/* Subroutine */ int sb04py_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, doublereal *
	dwork, integer *info, ftnlen trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer k1, k2, l1, l2;
    static doublereal a11, db, p11, p12, p21, p22, da11, vec[4]	/* was [2][2] 
	    */, dum[1], eps, sgn;
    static integer mnk1, mnk2, mnl1, mnl2;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, sumr;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sb04px_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer knext, lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notrna, notrnb;
    static doublereal smlnum;


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

/*     To solve for X the discrete-time Sylvester equation */

/*        op(A)*X*op(B) + ISGN*X = scale*C, */

/*     where op(A) = A or A**T, A and B are both upper quasi-triangular, */
/*     and ISGN = 1 or -1. A is M-by-M and B is N-by-N; the right hand */
/*     side C and the solution X are M-by-N; and scale is an output scale */
/*     factor, set less than or equal to 1 to avoid overflow in X. The */
/*     solution matrix X is overwritten onto C. */

/*     A and B must be in Schur canonical form (as returned by LAPACK */
/*     Library routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op(B) to be used, as follows: */
/*             = 'N':  op(B) = B    (No transpose); */
/*             = 'T':  op(B) = B**T (Transpose); */
/*             = 'C':  op(B) = B**T (Conjugate transpose = Transpose). */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A, and the number of rows in the */
/*             matrices X and C.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix B, and the number of columns in */
/*             the matrices X and C.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper quasi-triangular matrix A, in Schur canonical form. */
/*             The part of A below the first sub-diagonal is not */
/*             referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading N-by-N part of this array must contain the */
/*             upper quasi-triangular matrix B, in Schur canonical form. */
/*             The part of B below the first sub-diagonal is not */
/*             referenced. */

/*     LDB     (input) INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix C. */
/*             On exit, if INFO >= 0, the leading M-by-N part of this */
/*             array contains the solution matrix X. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  A and -ISGN*B have almost reciprocal eigenvalues; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrices A and B are unchanged). */

/*     METHOD */

/*     The solution matrix X is computed column-wise via a back */
/*     substitution scheme, an extension and refinement of the algorithm */
/*     in [1], similar to that used in [2] for continuous-time Sylvester */
/*     equations. A set of equivalent linear algebraic systems of */
/*     equations of order at most four are formed and solved using */
/*     Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is stable and reliable, since Gaussian elimination */
/*     with complete pivoting is used. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000. */
/*     D. Sima, University of Bucharest, April 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */
/*     Partly based on the routine SYLSV, A. Varga, 1992. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, matrix algebra, Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! notrnb && ! lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*isgn != 1 && *isgn != -1) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*m)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*m)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB04PY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *scale = 1.;
    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Set constants to control overflow. */

    eps = dlamch_("Precision", (ftnlen)9);
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
    bignum = 1. / smlnum;

/* Computing MAX */
    d__1 = smlnum, d__2 = eps * dlange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * dlange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
    smin = max(d__1,d__2);

    sgn = (doublereal) (*isgn);

    if (notrna && notrnb) {

/*        Solve    A*X*B + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*           A(K,K)*X(K,L)*B(L,L) + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                       M */
/*           R(K,L) = { SUM [A(K,J)*X(J,L)] } * B(L,L) + */
/*                     J=K+1 */
/*                       M             L-1 */
/*                      SUM { A(K,J) * SUM [X(J,I)*B(I,L)] }. */
/*                      J=K            I=1 */

/*        Start column loop (index = L) */
/*        L1 (L2) : column index of the first (last) row of X(K,L). */

	lnext = 1;

	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    if (l < lnext) {
		goto L60;
	    }
	    l1 = l;
	    if (l == *n) {
		l2 = l;
	    } else {
		if (b[l + 1 + l * b_dim1] != 0.) {
		    l2 = l + 1;
		} else {
		    l2 = l;
		}
		lnext = l2 + 1;
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

	    knext = *m;

	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L50;
		}
		k2 = k;
		if (k == 1) {
		    k1 = k;
		} else {
		    if (a[k + (k - 1) * a_dim1] != 0.) {
			k1 = k - 1;
		    } else {
			k1 = k;
		    }
		    knext = k1 - 1;
		}

/* Computing MIN */
		i__2 = k1 + 1;
		mnk1 = min(i__2,*m);
/* Computing MIN */
		i__2 = k2 + 1;
		mnk2 = min(i__2,*m);
		i__2 = *m - k2;
		p11 = ddot_(&i__2, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 + 
			l1 * c_dim1], &c__1);
		i__2 = l1 - 1;
		dwork[k1] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			b_dim1 + 1], &c__1);

		if (l1 == l2 && k1 == k2) {

		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
		    scaloc = 1.;

		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L10: */
			}

			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = *m - k2;
		    p21 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
		    i__2 = l1 - 1;
		    dwork[k2] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &b[l1 + l1 * 
			    b_dim1], &a[k1 + k1 * a_dim1], lda, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L20: */
			}

			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k2 + l1 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = *m - k1;
		    p12 = ddot_(&i__2, &a[k1 + mnk1 * a_dim1], lda, &c__[mnk1 
			    + l2 * c_dim1], &c__1);
		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

		    i__2 = l1 - 1;
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[k1 + k1 * a_dim1]
			    , &b[l1 + l1 * b_dim1], ldb, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L30: */
			}

			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1 + *m], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = *m - k2;
		    p21 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
		    i__2 = *m - k2;
		    p12 = ddot_(&i__2, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);
		    i__2 = *m - k2;
		    p22 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);

		    i__2 = l1 - 1;
		    dwork[k2] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

		    i__2 = l1 - 1;
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    i__2 = l1 - 1;
		    dwork[k2 + *m] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l2 + l1 * b_dim1]);

		    i__2 = *m - k1 + 1;
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l1 + l2 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

		    sb04px_(&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 + 
			    k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec,
			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L40: */
			}

			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
			i__2 = *m - k1 + 1;
			dscal_(&i__2, &scaloc, &dwork[k1 + *m], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[2];
		    c__[k2 + l1 * c_dim1] = x[1];
		    c__[k2 + l2 * c_dim1] = x[3];
		}

L50:
		;
	    }

L60:
	    ;
	}

    } else if (! notrna && notrnb) {

/*        Solve     A'*X*B + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*         A(K,K)'*X(K,L)*B(L,L) + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                      K-1 */
/*           R(K,L) = { SUM [A(J,K)'*X(J,L)] } * B(L,L) + */
/*                      J=1 */
/*                       K              L-1 */
/*                      SUM A(J,K)' * { SUM [X(J,I)*B(I,L)] }. */
/*                      J=1             I=1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

	lnext = 1;

	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    if (l < lnext) {
		goto L120;
	    }
	    l1 = l;
	    if (l == *n) {
		l2 = l;
	    } else {
		if (b[l + 1 + l * b_dim1] != 0.) {
		    l2 = l + 1;
		} else {
		    l2 = l;
		}
		lnext = l2 + 1;
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

	    knext = 1;

	    i__2 = *m;
	    for (k = 1; k <= i__2; ++k) {
		if (k < knext) {
		    goto L110;
		}
		k1 = k;
		if (k == *m) {
		    k2 = k;
		} else {
		    if (a[k + 1 + k * a_dim1] != 0.) {
			k2 = k + 1;
		    } else {
			k2 = k;
		    }
		    knext = k2 + 1;
		}

		i__3 = k1 - 1;
		p11 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			c_dim1 + 1], &c__1);
		i__3 = l1 - 1;
		dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			b_dim1 + 1], &c__1);

		if (l1 == l2 && k1 == k2) {

		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
		    scaloc = 1.;

		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {

			i__3 = *n;
			for (j = 1; j <= i__3; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L70: */
			}

			dscal_(&k1, &scaloc, &dwork[1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];

		} else if (l1 == l2 && k1 != k2) {

		    i__3 = k1 - 1;
		    p21 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
		    i__3 = l1 - 1;
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &b[l1 + l1 * b_dim1]
			    , &a[k1 + k1 * a_dim1], lda, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__3 = *n;
			for (j = 1; j <= i__3; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L80: */
			}

			dscal_(&k2, &scaloc, &dwork[1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k2 + l1 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 == k2) {

		    i__3 = k1 - 1;
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

		    i__3 = l1 - 1;
		    dwork[k1 + *m] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[k1 + k1 * a_dim1]
			    , &b[l1 + l1 * b_dim1], ldb, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__3 = *n;
			for (j = 1; j <= i__3; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L90: */
			}

			dscal_(&k1, &scaloc, &dwork[1], &c__1);
			dscal_(&k1, &scaloc, &dwork[*m + 1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 != k2) {

		    i__3 = k1 - 1;
		    p21 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
		    i__3 = k1 - 1;
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
		    i__3 = k1 - 1;
		    p22 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

		    i__3 = l1 - 1;
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l2 + l1 * b_dim1]);

		    i__3 = l1 - 1;
		    dwork[k1 + *m] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    i__3 = l1 - 1;
		    dwork[k2 + *m] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l1 + l2 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

		    sb04px_(&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__3 = *n;
			for (j = 1; j <= i__3; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L100: */
			}

			dscal_(&k2, &scaloc, &dwork[1], &c__1);
			dscal_(&k2, &scaloc, &dwork[*m + 1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[2];
		    c__[k2 + l1 * c_dim1] = x[1];
		    c__[k2 + l2 * c_dim1] = x[3];
		}

L110:
		;
	    }

L120:
	    ;
	}

    } else if (! notrna && ! notrnb) {

/*        Solve    A'*X*B' + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        top-right corner column by column by */

/*           A(K,K)'*X(K,L)*B(L,L)' + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                      K-1 */
/*           R(K,L) = { SUM [A(J,K)'*X(J,L)] } * B(L,L)' + */
/*                      J=1 */
/*                       K               N */
/*                      SUM A(J,K)' * { SUM [X(J,I)*B(L,I)'] }. */
/*                      J=1            I=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

	lnext = *n;

	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L180;
	    }
	    l2 = l;
	    if (l == 1) {
		l1 = l;
	    } else {
		if (b[l + (l - 1) * b_dim1] != 0.) {
		    l1 = l - 1;
		} else {
		    l1 = l;
		}
		lnext = l1 - 1;
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

	    knext = 1;

	    i__1 = *m;
	    for (k = 1; k <= i__1; ++k) {
		if (k < knext) {
		    goto L170;
		}
		k1 = k;
		if (k == *m) {
		    k2 = k;
		} else {
		    if (a[k + 1 + k * a_dim1] != 0.) {
			k2 = k + 1;
		    } else {
			k2 = k;
		    }
		    knext = k2 + 1;
		}

/* Computing MIN */
		i__2 = l1 + 1;
		mnl1 = min(i__2,*n);
/* Computing MIN */
		i__2 = l2 + 1;
		mnl2 = min(i__2,*n);
		i__2 = k1 - 1;
		p11 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			c_dim1 + 1], &c__1);
		i__2 = *n - l2;
		dwork[k1] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], ldc, &b[l1 
			+ mnl2 * b_dim1], ldb);

		if (l1 == l2 && k1 == k2) {
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
		    scaloc = 1.;

		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L130: */
			}

			dscal_(&k1, &scaloc, &dwork[1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = k1 - 1;
		    p21 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
		    i__2 = *n - l1;
		    dwork[k2] = ddot_(&i__2, &c__[k2 + mnl1 * c_dim1], ldc, &
			    b[l1 + mnl1 * b_dim1], ldb);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &b[l1 + l1 * b_dim1]
			    , &a[k1 + k1 * a_dim1], lda, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L140: */
			}

			dscal_(&k2, &scaloc, &dwork[1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k2 + l1 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = k1 - 1;
		    p12 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

		    i__2 = *n - l2;
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[k1 + k1 * 
			    a_dim1], &b[l1 + l1 * b_dim1], ldb, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L150: */
			}

			dscal_(&k1, &scaloc, &dwork[1], &c__1);
			dscal_(&k1, &scaloc, &dwork[*m + 1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = k1 - 1;
		    p21 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
		    i__2 = k1 - 1;
		    p12 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
		    i__2 = k1 - 1;
		    p22 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

		    i__2 = *n - l2;
		    dwork[k2] = ddot_(&i__2, &c__[k2 + mnl2 * c_dim1], ldc, &
			    b[l1 + mnl2 * b_dim1], ldb);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l1 + l2 * b_dim1]);

		    i__2 = *n - l2;
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    i__2 = *n - l2;
		    dwork[k2 + *m] = ddot_(&i__2, &c__[k2 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l2 + l1 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

		    sb04px_(&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 *
			     a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L160: */
			}

			dscal_(&k2, &scaloc, &dwork[1], &c__1);
			dscal_(&k2, &scaloc, &dwork[*m + 1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[2];
		    c__[k2 + l1 * c_dim1] = x[1];
		    c__[k2 + l2 * c_dim1] = x[3];
		}

L170:
		;
	    }

L180:
	    ;
	}

    } else {

/*        Solve    A*X*B' + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-right corner column by column by */

/*            A(K,K)*X(K,L)*B(L,L)' + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                       M */
/*           R(K,L) = { SUM [A(K,J)*X(J,L)] } * B(L,L)' + */
/*                     J=K+1 */
/*                       M              N */
/*                      SUM { A(K,J) * SUM [X(J,I)*B(L,I)'] }. */
/*                      J=K           I=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

	lnext = *n;

	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L240;
	    }
	    l2 = l;
	    if (l == 1) {
		l1 = l;
	    } else {
		if (b[l + (l - 1) * b_dim1] != 0.) {
		    l1 = l - 1;
		} else {
		    l1 = l;
		}
		lnext = l1 - 1;
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

	    knext = *m;

	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L230;
		}
		k2 = k;
		if (k == 1) {
		    k1 = k;
		} else {
		    if (a[k + (k - 1) * a_dim1] != 0.) {
			k1 = k - 1;
		    } else {
			k1 = k;
		    }
		    knext = k1 - 1;
		}

/* Computing MIN */
		i__1 = k1 + 1;
		mnk1 = min(i__1,*m);
/* Computing MIN */
		i__1 = k2 + 1;
		mnk2 = min(i__1,*m);
/* Computing MIN */
		i__1 = l1 + 1;
		mnl1 = min(i__1,*n);
/* Computing MIN */
		i__1 = l2 + 1;
		mnl2 = min(i__1,*n);
		i__1 = *m - k2;
		p11 = ddot_(&i__1, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 + 
			l1 * c_dim1], &c__1);
		i__1 = *n - l2;
		dwork[k1] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], ldc, &b[l1 
			+ mnl2 * b_dim1], ldb);

		if (l1 == l2 && k1 == k2) {

		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
		    scaloc = 1.;

		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {

			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L190: */
			}

			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];

		} else if (l1 == l2 && k1 != k2) {

		    i__1 = *m - k2;
		    p21 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
		    i__1 = *n - l1;
		    dwork[k2] = ddot_(&i__1, &c__[k2 + mnl1 * c_dim1], ldc, &
			    b[l1 + mnl1 * b_dim1], ldb);
		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &b[l1 + l1 * 
			    b_dim1], &a[k1 + k1 * a_dim1], lda, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L200: */
			}

			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k2 + l1 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 == k2) {

		    i__1 = *m - k1;
		    p12 = ddot_(&i__1, &a[k1 + mnk1 * a_dim1], lda, &c__[mnk1 
			    + l2 * c_dim1], &c__1);
		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

		    i__1 = *n - l2;
		    dwork[k1 + *m] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    d__1 = -sgn;
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[k1 + k1 * 
			    a_dim1], &b[l1 + l1 * b_dim1], ldb, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L210: */
			}

			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1 + *m], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[1];

		} else if (l1 != l2 && k1 != k2) {

		    i__1 = *m - k2;
		    p21 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
		    i__1 = *m - k2;
		    p12 = ddot_(&i__1, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);
		    i__1 = *m - k2;
		    p22 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);

		    i__1 = *n - l2;
		    dwork[k2] = ddot_(&i__1, &c__[k2 + mnl2 * c_dim1], ldc, &
			    b[l1 + mnl2 * b_dim1], ldb);
		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l1 + l2 * b_dim1]);

		    i__1 = *n - l2;
		    dwork[k1 + *m] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    i__1 = *n - l2;
		    dwork[k2 + *m] = ddot_(&i__1, &c__[k2 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

		    i__1 = *m - k1 + 1;
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l2 + l1 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

		    sb04px_(&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {

			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L220: */
			}

			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
			i__1 = *m - k1 + 1;
			dscal_(&i__1, &scaloc, &dwork[k1 + *m], &c__1);
			*scale *= scaloc;
		    }
		    c__[k1 + l1 * c_dim1] = x[0];
		    c__[k1 + l2 * c_dim1] = x[2];
		    c__[k2 + l1 * c_dim1] = x[1];
		    c__[k2 + l2 * c_dim1] = x[3];
		}

L230:
		;
	    }

L240:
	    ;
	}

    }

    return 0;
/* *** Last line of SB04PY *** */
} /* sb04py_ */

