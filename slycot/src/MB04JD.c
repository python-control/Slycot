/* MB04JD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04jd_(integer *n, integer *m, integer *p, integer *l, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	tau, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static doublereal first;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dgelqf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen), dormlq_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static doublereal wrkopt;


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

/*     To compute an LQ factorization of an n-by-m matrix A (A = L * Q), */
/*     having a min(n,p)-by-p zero triangle in the upper right-hand side */
/*     corner, as shown below, for n = 8, m = 7, and p = 2: */

/*            [ x x x x x 0 0 ] */
/*            [ x x x x x x 0 ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*        A = [ x x x x x x x ], */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */

/*     and optionally apply the transformations to an l-by-m matrix B */
/*     (from the right). The problem structure is exploited. This */
/*     computation is useful, for instance, in combined measurement and */
/*     time update of one iteration of the time-invariant Kalman filter */
/*     (square root covariance filter). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix A.  M >= 0. */

/*     P       (input) INTEGER */
/*             The order of the zero triagle.  P >= 0. */

/*     L       (input) INTEGER */
/*             The number of rows of the matrix B.  L >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix A. The elements corresponding to the */
/*             zero MIN(N,P)-by-P upper trapezoidal/triangular part */
/*             (if P > 0) are not referenced. */
/*             On exit, the elements on and below the diagonal of this */
/*             array contain the N-by-MIN(N,M) lower trapezoidal matrix */
/*             L (L is lower triangular, if N <= M) of the LQ */
/*             factorization, and the relevant elements above the */
/*             diagonal contain the trailing components (the vectors v, */
/*             see Method) of the elementary reflectors used in the */
/*             factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the updated matrix B. */
/*             If L = 0, this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,L). */

/*     TAU     (output) DOUBLE PRECISION array, dimension MIN(N,M) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,N-1,N-P,L). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses min(N,M) Householder transformations exploiting */
/*     the zero pattern of the matrix.  A Householder matrix has the form */

/*                                     ( 1 ), */
/*        H  = I - tau *u *u',    u  = ( v ) */
/*         i          i  i  i      i   (  i) */

/*     where v  is an (M-P+I-2)-vector.  The components of v  are stored */
/*            i                                             i */
/*     in the i-th row of A, beginning from the location i+1, and tau */
/*                                                                   i */
/*     is stored in TAU(i). */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary reflector, LQ factorization, orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*l)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1, i__1 = max(i__1,i__2), i__2 = *n - *p, i__1 =
		 max(i__1,i__2);
	if (*ldwork < max(i__1,*l)) {
	    *info = -11;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB04JD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*m,*n) == 0) {
	dwork[1] = 1.;
	return 0;
    } else if (*m <= *p + 1) {
	i__1 = min(*n,*m);
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tau[i__] = 0.;
/* L5: */
	}
	dwork[1] = 1.;
	return 0;
    }

/*     Annihilate the superdiagonal elements of A and apply the */
/*     transformations to B, if L > 0. */
/*     Workspace: need MAX(N-1,L). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    i__1 = min(*n,*p);
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Exploit the structure of the I-th row of A. */

	i__2 = *m - *p;
	dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + (i__ + 1) * a_dim1], 
		lda, &tau[i__]);
	if (tau[i__] != 0.) {

	    first = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.;

	    if (i__ < *n) {
		i__2 = *n - i__;
		i__3 = *m - *p;
		dlarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &dwork[1], 
			(ftnlen)5);
	    }
	    if (*l > 0) {
		i__2 = *m - *p;
		dlarf_("Right", l, &i__2, &a[i__ + i__ * a_dim1], lda, &tau[
			i__], &b[i__ * b_dim1 + 1], ldb, &dwork[1], (ftnlen)5)
			;
	    }

	    a[i__ + i__ * a_dim1] = first;
	}
/* L10: */
    }

/* Computing MAX */
    d__1 = 1., d__2 = (doublereal) (*n - 1), d__1 = max(d__1,d__2), d__2 = (
	    doublereal) (*l);
    wrkopt = max(d__1,d__2);

/*     Fast LQ factorization of the remaining trailing submatrix, if any. */
/*     Workspace: need N-P;  prefer (N-P)*NB. */

    if (*n > *p) {
	i__1 = *n - *p;
	i__2 = *m - *p;
	dgelqf_(&i__1, &i__2, &a[*p + 1 + (*p + 1) * a_dim1], lda, &tau[*p + 
		1], &dwork[1], ldwork, info);
	wrkopt = max(wrkopt,dwork[1]);

	if (*l > 0) {

/*           Apply the transformations to B. */
/*           Workspace: need L;  prefer L*NB. */

	    i__1 = *m - *p;
	    i__2 = min(*n,*m) - *p;
	    dormlq_("Right", "Transpose", l, &i__1, &i__2, &a[*p + 1 + (*p + 
		    1) * a_dim1], lda, &tau[*p + 1], &b[(*p + 1) * b_dim1 + 1]
		    , ldb, &dwork[1], ldwork, info, (ftnlen)5, (ftnlen)9);
	    wrkopt = max(wrkopt,dwork[1]);
	}
    }

    dwork[1] = wrkopt;
    return 0;
/* *** Last line of MB04JD *** */
} /* mb04jd_ */

