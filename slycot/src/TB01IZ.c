/* TB01IZ.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb01iz_(char *job, integer *n, integer *m, integer *p, 
	doublereal *maxred, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);

    /* Local variables */
    static doublereal f, g;
    static integer i__, j;
    static doublereal s, ca, co, ra, ro;
    static integer ica, ira;
    static doublereal sred;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical withb, withc;
    static doublereal snorm, sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical noconv;
    static doublereal maxnrm;
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


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

/*     To reduce the 1-norm of a system matrix */

/*             S =  ( A  B ) */
/*                  ( C  0 ) */

/*     corresponding to the triple (A,B,C), by balancing. This involves */
/*     a diagonal similarity transformation inv(D)*A*D applied */
/*     iteratively to A to make the rows and columns of */
/*                           -1 */
/*                  diag(D,I)  * S * diag(D,I) */

/*     as close in norm as possible. */

/*     The balancing can be performed optionally on the following */
/*     particular system matrices */

/*              S = A,    S = ( A  B )    or    S = ( A ) */
/*                                                  ( C ) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates which matrices are involved in balancing, as */
/*             follows: */
/*             = 'A':  All matrices are involved in balancing; */
/*             = 'B':  B and A matrices are involved in balancing; */
/*             = 'C':  C and A matrices are involved in balancing; */
/*             = 'N':  B and C matrices are not involved in balancing. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of matrix B */
/*             and the number of columns of matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER. */
/*             The number of columns of matrix B. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER. */
/*             The number of rows of matrix C. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     MAXRED  (input/output) DOUBLE PRECISION */
/*             On entry, the maximum allowed reduction in the 1-norm of */
/*             S (in an iteration) if zero rows or columns are */
/*             encountered. */
/*             If MAXRED > 0.0, MAXRED must be larger than one (to enable */
/*             the norm reduction). */
/*             If MAXRED <= 0.0, then the value 10.0 for MAXRED is */
/*             used. */
/*             On exit, if the 1-norm of the given matrix S is non-zero, */
/*             the ratio between the 1-norm of the given matrix and the */
/*             1-norm of the balanced matrix. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced matrix inv(D)*A*D. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, if M > 0, the leading N-by-M part of this array */
/*             must contain the system input matrix B. */
/*             On exit, if M > 0, the leading N-by-M part of this array */
/*             contains the balanced matrix inv(D)*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0. */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, if P > 0, the leading P-by-N part of this array */
/*             must contain the system output matrix C. */
/*             On exit, if P > 0, the leading P-by-N part of this array */
/*             contains the balanced matrix C*D. */
/*             The array C is not referenced if P = 0. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to S.  If D(j) is the scaling */
/*             factor applied to row and column j, then SCALE(j) = D(j), */
/*             for j = 1,...,N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation */
/*                           -1 */
/*                  diag(D,I)  * S * diag(D,I) */

/*     to make the 1-norms of each row of the first N rows of S and its */
/*     corresponding column nearly equal. */

/*     Information about the diagonal matrix D is returned in the vector */
/*     SCALE. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balancing, eigenvalue, matrix algebra, matrix operations, */
/*     similarity transformation. */

/*  ********************************************************************* */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

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
    --scale;

    /* Function Body */
    *info = 0;
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*maxred > 0. && *maxred < 1.) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*m > 0 && *ldb < max(1,*n) || *m == 0 && *ldb < 1) {
	*info = -9;
    } else if (*ldc < max(1,*p)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TB01IZ", &i__1, (ftnlen)6);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

/*     Compute the 1-norm of the required part of matrix S and exit if */
/*     it is zero. */

    snorm = 0.;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	scale[j] = 1.;
	co = dzasum_(n, &a[j * a_dim1 + 1], &c__1);
	if (withc && *p > 0) {
	    co += dzasum_(p, &c__[j * c_dim1 + 1], &c__1);
	}
	snorm = max(snorm,co);
/* L10: */
    }

    if (withb) {

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = snorm, d__2 = dzasum_(n, &b[j * b_dim1 + 1], &c__1);
	    snorm = max(d__1,d__2);
/* L20: */
	}

    }

    if (snorm == 0.) {
	return 0;
    }

/*     Set some machine parameters and the maximum reduction in the */
/*     1-norm of S if zero rows or columns are encountered. */

    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
    sfmax1 = 1. / sfmin1;
    sfmin2 = sfmin1 * 10.;
    sfmax2 = 1. / sfmin2;

    sred = *maxred;
    if (sred <= 0.) {
	sred = 10.;
    }

/* Computing MAX */
    d__1 = snorm / sred;
    maxnrm = max(d__1,sfmin1);

/*     Balance the matrix. */

/*     Iterative loop for norm reduction. */

L30:
    noconv = FALSE_;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	co = 0.;
	ro = 0.;

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (j == i__) {
		goto L40;
	    }
	    i__3 = j + i__ * a_dim1;
	    co += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + i__ * 
		    a_dim1]), abs(d__2));
	    i__3 = i__ + j * a_dim1;
	    ro += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2));
L40:
	    ;
	}

	ica = izamax_(n, &a[i__ * a_dim1 + 1], &c__1);
	ca = z_abs(&a[ica + i__ * a_dim1]);
	ira = izamax_(n, &a[i__ + a_dim1], lda);
	ra = z_abs(&a[i__ + ira * a_dim1]);

	if (withc && *p > 0) {
	    co += dzasum_(p, &c__[i__ * c_dim1 + 1], &c__1);
	    ica = izamax_(p, &c__[i__ * c_dim1 + 1], &c__1);
/* Computing MAX */
	    d__1 = ca, d__2 = z_abs(&c__[ica + i__ * c_dim1]);
	    ca = max(d__1,d__2);
	}

	if (withb && *m > 0) {
	    ro += dzasum_(m, &b[i__ + b_dim1], ldb);
	    ira = izamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
	    d__1 = ra, d__2 = z_abs(&b[i__ + ira * b_dim1]);
	    ra = max(d__1,d__2);
	}

/*        Special case of zero CO and/or RO. */

	if (co == 0. && ro == 0.) {
	    goto L90;
	}
	if (co == 0.) {
	    if (ro <= maxnrm) {
		goto L90;
	    }
	    co = maxnrm;
	}
	if (ro == 0.) {
	    if (co <= maxnrm) {
		goto L90;
	    }
	    ro = maxnrm;
	}

/*        Guard against zero CO or RO due to underflow. */

	g = ro / 10.;
	f = 1.;
	s = co + ro;
L50:
/* Computing MAX */
	d__1 = max(f,co);
/* Computing MIN */
	d__2 = min(ro,g);
	if (co >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
	    goto L60;
	}
	f *= 10.;
	co *= 10.;
	ca *= 10.;
	g /= 10.;
	ro /= 10.;
	ra /= 10.;
	goto L50;

L60:
	g = co / 10.;
L70:
/* Computing MIN */
	d__1 = min(f,co), d__1 = min(d__1,g);
	if (g < ro || max(ro,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
	    goto L80;
	}
	f /= 10.;
	co /= 10.;
	ca /= 10.;
	g /= 10.;
	ro *= 10.;
	ra *= 10.;
	goto L70;

/*        Now balance. */

L80:
	if (co + ro >= s * .95) {
	    goto L90;
	}
	if (f < 1. && scale[i__] < 1.) {
	    if (f * scale[i__] <= sfmin1) {
		goto L90;
	    }
	}
	if (f > 1. && scale[i__] > 1.) {
	    if (scale[i__] >= sfmax1 / f) {
		goto L90;
	    }
	}
	g = 1. / f;
	scale[i__] *= f;
	noconv = TRUE_;

	zdscal_(n, &g, &a[i__ + a_dim1], lda);
	zdscal_(n, &f, &a[i__ * a_dim1 + 1], &c__1);
	if (*m > 0) {
	    zdscal_(m, &g, &b[i__ + b_dim1], ldb);
	}
	if (*p > 0) {
	    zdscal_(p, &f, &c__[i__ * c_dim1 + 1], &c__1);
	}

L90:
	;
    }

    if (noconv) {
	goto L30;
    }

/*     Set the norm reduction parameter. */

    *maxred = snorm;
    snorm = 0.;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	co = dzasum_(n, &a[j * a_dim1 + 1], &c__1);
	if (withc && *p > 0) {
	    co += dzasum_(p, &c__[j * c_dim1 + 1], &c__1);
	}
	snorm = max(snorm,co);
/* L100: */
    }

    if (withb) {

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = snorm, d__2 = dzasum_(n, &b[j * b_dim1 + 1], &c__1);
	    snorm = max(d__1,d__2);
/* L110: */
	}

    }
    *maxred /= snorm;
    return 0;
/* *** Last line of TB01IZ *** */
} /* tb01iz_ */

