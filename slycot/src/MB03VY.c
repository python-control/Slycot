/* MB03VY.f -- translated by f2c (version 20100827).
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

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;

/* Subroutine */ int mb03vy_(integer *n, integer *p, integer *ilo, integer *
	ihi, doublereal *a, integer *lda1, integer *lda2, doublereal *tau, 
	integer *ldtau, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, tau_dim1, tau_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, nh;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dorghr_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
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

/*     To generate the real orthogonal matrices Q_1, Q_2, ..., Q_p, */
/*     which are defined as the product of ihi-ilo elementary reflectors */
/*     of order n, as returned by SLICOT Library routine MB03VD: */

/*        Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices Q_1, Q_2, ..., Q_p.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number p of transformation matrices.  P >= 1. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             The values of the indices ilo and ihi, respectively, used */
/*             in the previous call of the SLICOT Library routine MB03VD. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA1,LDA2,N) */
/*             On entry, the leading N-by-N strictly lower triangular */
/*             part of A(*,*,j) must contain the vectors which define the */
/*             elementary reflectors used for reducing A_j, as returned */
/*             by SLICOT Library routine MB03VD, j = 1, ..., p. */
/*             On exit, the leading N-by-N part of A(*,*,j) contains the */
/*             N-by-N orthogonal matrix Q_j, j = 1, ..., p. */

/*     LDA1    INTEGER */
/*             The first leading dimension of the array A. */
/*             LDA1 >= max(1,N). */

/*     LDA2    INTEGER */
/*             The second leading dimension of the array A. */
/*             LDA2 >= max(1,N). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (LDTAU,P) */
/*             The leading N-1 elements in the j-th column must contain */
/*             the scalar factors of the elementary reflectors used to */
/*             form the matrix Q_j, as returned by SLICOT Library routine */
/*             MB03VD. */

/*     LDTAU   INTEGER */
/*             The leading dimension of the array TAU. */
/*             LDTAU >= max(1,N-1). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Each matrix Q_j is generated as the product of the elementary */
/*     reflectors used for reducing A_j. Standard LAPACK routines for */
/*     Hessenberg and QR decompositions are used. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P. */
/*         The periodic Schur decomposition: algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         Periodic Schur form and some matrix equations. */
/*         Proc. of the Symposium on the Mathematical Theory of Networks */
/*         and Systems (MTNS'93), Regensburg, Germany (U. Helmke, */
/*         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga, */
/*     German Aerospace Center, DLR Oberpfaffenhofen, February 1999. */
/*     Partly based on the routine PSHTR by A. Varga */
/*     (DLR Oberpfaffenhofen), November 26, 1995. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, periodic systems, */
/*     similarity transformation, triangular form. */

/*     ****************************************************************** */

/*     .. Parameters .. */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda1;
    a_dim2 = *lda2;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    tau_dim1 = *ldtau;
    tau_offset = 1 + tau_dim1;
    tau -= tau_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*p < 1) {
	*info = -2;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -3;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -4;
    } else if (*lda1 < max(1,*n)) {
	*info = -6;
    } else if (*lda2 < max(1,*n)) {
	*info = -7;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*ldtau < max(i__1,i__2)) {
	    *info = -9;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03VY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Generate the orthogonal matrix Q_1. */

    dorghr_(n, ilo, ihi, &a[a_offset], lda1, &tau[tau_offset], &dwork[1], 
	    ldwork, info);
    wrkopt = dwork[1];

    nh = *ihi - *ilo + 1;

    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {

/*        Generate the orthogonal matrix Q_j. */
/*        Set the first ILO-1 and the last N-IHI rows and columns of Q_j */
/*        to those of the unit matrix. */

	i__2 = *ilo - 1;
	dlaset_("Full", n, &i__2, &c_b5, &c_b6, &a[(j * a_dim2 + 1) * a_dim1 
		+ 1], lda1, (ftnlen)4);
	i__2 = *ilo - 1;
	dlaset_("Full", &i__2, &nh, &c_b5, &c_b5, &a[(*ilo + j * a_dim2) * 
		a_dim1 + 1], lda1, (ftnlen)4);
	if (nh > 1) {
	    i__2 = nh - 1;
	    dorgqr_(&nh, &nh, &i__2, &a[*ilo + (*ilo + j * a_dim2) * a_dim1], 
		    lda1, &tau[*ilo + j * tau_dim1], &dwork[1], ldwork, info);
	}
	if (*ihi < *n) {
	    i__2 = *n - *ihi;
	    dlaset_("Full", &i__2, &nh, &c_b5, &c_b5, &a[*ihi + 1 + (*ilo + j 
		    * a_dim2) * a_dim1], lda1, (ftnlen)4);
	    i__2 = *n - *ihi;
	    dlaset_("Full", ihi, &i__2, &c_b5, &c_b5, &a[(*ihi + 1 + j * 
		    a_dim2) * a_dim1 + 1], lda1, (ftnlen)4);
	    i__2 = *n - *ihi;
	    i__3 = *n - *ihi;
	    dlaset_("Full", &i__2, &i__3, &c_b5, &c_b6, &a[*ihi + 1 + (*ihi + 
		    1 + j * a_dim2) * a_dim1], lda1, (ftnlen)4);
	}
/* L20: */
    }

    dwork[1] = max(wrkopt,dwork[1]);
    return 0;

/* *** Last line of MB03VY *** */
} /* mb03vy_ */

