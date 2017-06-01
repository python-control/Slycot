/* MB04MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04md_(integer *n, doublereal *maxred, doublereal *a, 
	integer *lda, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j;
    static doublereal r__, s, ca, ra;
    static integer ica, ira;
    static doublereal sred;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm, sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconv;
    static doublereal maxnrm;


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

/*     To reduce the 1-norm of a general real matrix A by balancing. */
/*     This involves diagonal similarity transformations applied */
/*     iteratively to A to make the rows and columns as close in norm as */
/*     possible. */

/*     This routine can be used instead LAPACK Library routine DGEBAL, */
/*     when no reduction of the 1-norm of the matrix is possible with */
/*     DGEBAL, as for upper triangular matrices. LAPACK Library routine */
/*     DGEBAK, with parameters ILO = 1, IHI = N, and JOB = 'S', should */
/*     be used to apply the backward transformation. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     MAXRED  (input/output) DOUBLE PRECISION */
/*             On entry, the maximum allowed reduction in the 1-norm of */
/*             A (in an iteration) if zero rows or columns are */
/*             encountered. */
/*             If MAXRED > 0.0, MAXRED must be larger than one (to enable */
/*             the norm reduction). */
/*             If MAXRED <= 0.0, then the value 10.0 for MAXRED is */
/*             used. */
/*             On exit, if the 1-norm of the given matrix A is non-zero, */
/*             the ratio between the 1-norm of the given matrix and the */
/*             1-norm of the balanced matrix. Usually, this ratio will be */
/*             larger than one, but it can sometimes be one, or even less */
/*             than one (for instance, for some companion matrices). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the input matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to A.  If D(j) is the scaling */
/*             factor applied to row and column j, then SCALE(j) = D(j), */
/*             for j = 1,...,N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation inv(D) * A * D to make the 1-norms of each row */
/*     of A and its corresponding column nearly equal. */

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

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04AD by T.W.C. Williams, */
/*     Kingston Polytechnic, United Kingdom, October 1984. */
/*     This subroutine is based on LAPACK routine DGEBAL, and routine */
/*     BALABC (A. Varga, German Aerospace Research Establishment, DLR). */


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
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --scale;

    /* Function Body */
    *info = 0;

    if (*n < 0) {
	*info = -1;
    } else if (*maxred > 0. && *maxred < 1.) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04MD", &i__1, (ftnlen)6);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	scale[i__] = 1.;
/* L10: */
    }

/*     Compute the 1-norm of matrix A and exit if it is zero. */

    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &scale[1], (ftnlen)6);
    if (anorm == 0.) {
	return 0;
    }

/*     Set some machine parameters and the maximum reduction in the */
/*     1-norm of A if zero rows or columns are encountered. */

    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
    sfmax1 = 1. / sfmin1;
    sfmin2 = sfmin1 * 10.;
    sfmax2 = 1. / sfmin2;

    sred = *maxred;
    if (sred <= 0.) {
	sred = 10.;
    }

/* Computing MAX */
    d__1 = anorm / sred;
    maxnrm = max(d__1,sfmin1);

/*     Balance the matrix. */

/*     Iterative loop for norm reduction. */

L20:
    noconv = FALSE_;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__ = 0.;
	r__ = 0.;

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (j == i__) {
		goto L30;
	    }
	    c__ += (d__1 = a[j + i__ * a_dim1], abs(d__1));
	    r__ += (d__1 = a[i__ + j * a_dim1], abs(d__1));
L30:
	    ;
	}
	ica = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
	ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
	ira = idamax_(n, &a[i__ + a_dim1], lda);
	ra = (d__1 = a[i__ + ira * a_dim1], abs(d__1));

/*        Special case of zero C and/or R. */

	if (c__ == 0. && r__ == 0.) {
	    goto L80;
	}
	if (c__ == 0.) {
	    if (r__ <= maxnrm) {
		goto L80;
	    }
	    c__ = maxnrm;
	}
	if (r__ == 0.) {
	    if (c__ <= maxnrm) {
		goto L80;
	    }
	    r__ = maxnrm;
	}

/*        Guard against zero C or R due to underflow. */

	g = r__ / 10.;
	f = 1.;
	s = c__ + r__;
L40:
/* Computing MAX */
	d__1 = max(f,c__);
/* Computing MIN */
	d__2 = min(r__,g);
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
	    goto L50;
	}
	f *= 10.;
	c__ *= 10.;
	ca *= 10.;
	r__ /= 10.;
	g /= 10.;
	ra /= 10.;
	goto L40;

L50:
	g = c__ / 10.;
L60:
/* Computing MIN */
	d__1 = min(f,c__), d__1 = min(d__1,g);
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
	    goto L70;
	}
	f /= 10.;
	c__ /= 10.;
	g /= 10.;
	ca /= 10.;
	r__ *= 10.;
	ra *= 10.;
	goto L60;

/*        Now balance. */

L70:
	if (c__ + r__ >= s * .95) {
	    goto L80;
	}
	if (f < 1. && scale[i__] < 1.) {
	    if (f * scale[i__] <= sfmin1) {
		goto L80;
	    }
	}
	if (f > 1. && scale[i__] > 1.) {
	    if (scale[i__] >= sfmax1 / f) {
		goto L80;
	    }
	}
	g = 1. / f;
	scale[i__] *= f;
	noconv = TRUE_;

	dscal_(n, &g, &a[i__ + a_dim1], lda);
	dscal_(n, &f, &a[i__ * a_dim1 + 1], &c__1);

L80:
	;
    }

    if (noconv) {
	goto L20;
    }

/*     Set the norm reduction parameter. */

    *maxred = anorm / dlange_("1-norm", n, n, &a[a_offset], lda, &scale[1], (
	    ftnlen)6);

    return 0;
/* *** End of MB04MD *** */
} /* mb04md_ */

