/* AB13CD.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b1 = {1.,0.};
static doublereal c_b16 = -1.;
static doublereal c_b17 = 1.;
static doublereal c_b38 = 0.;

doublereal ab13cd_(integer *n, integer *m, integer *np, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, doublecomplex *cwork, integer *
	lcwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, k, l, iw2, iw3, iw4, iw5, iw6, iw7, iw8, iw9;
    static doublereal den;
    static integer iw10, iw11, iw12;
    static doublereal rat;
    static integer icw2, icw3, icw4, sdim, iter;
    static doublereal temp;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal gamma, fpeak, omega;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    extern logical sb02cx_();
    extern /* Subroutine */ int dgesv_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    extern logical sb02mv_();
    extern /* Subroutine */ int mb01rx_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer icwrk;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static doublereal wimax, wrmin;
    extern /* Subroutine */ int dposv_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dsyrk_(char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), zgesv_(integer *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal gammal, gammau;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static integer lwamax, lcwamx;
    static doublereal ratmax;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer mincwr;
    static logical complx;
    extern /* Subroutine */ int zgesvd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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

/*     To compute the H-infinity norm of the continuous-time stable */
/*     system */

/*                          | A | B | */
/*                   G(s) = |---|---| . */
/*                          | C | D | */

/*     FUNCTION VALUE */

/*     AB13CD  DOUBLE PRECISION */
/*             If INFO = 0, the H-infinity norm of the system, HNORM, */
/*             i.e., the peak gain of the frequency response (as measured */
/*             by the largest singular value in the MIMO case). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used to set the accuracy in determining the */
/*             norm. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK, and DWORK(2) contains the frequency where the */
/*             gain of the frequency response achieves its peak value */
/*             HNORM. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(2,4*N*N+2*M*M+3*M*N+M*NP+2*(N+NP)*NP+10*N+ */
/*                             6*max(M,NP)). */
/*             For good performance, LDWORK must generally be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) contains the optimal value */
/*             of LCWORK. */

/*     LCWORK  INTEGER */
/*             The dimension of the array CWORK. */
/*             LCWORK >= max(1,(N+M)*(N+NP)+3*max(M,NP)). */
/*             For good performance, LCWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the system is unstable; */
/*             = 2:  the tolerance is too small (the algorithm for */
/*                   computing the H-infinity norm did not converge); */
/*             = 3:  errors in computing the eigenvalues of A or of the */
/*                   Hamiltonian matrix (the QR algorithm did not */
/*                   converge); */
/*             = 4:  errors in computing singular values. */

/*     METHOD */

/*     The routine implements the method presented in [1]. */

/*     REFERENCES */

/*     [1] Bruinsma, N.A. and Steinbuch, M. */
/*         A fast algorithm to compute the Hinfinity-norm of a transfer */
/*         function matrix. */
/*         Systems & Control Letters, vol. 14, pp. 287-293, 1990. */

/*     NUMERICAL ASPECTS */

/*     If the algorithm does not converge (INFO = 2), the tolerance must */
/*     be increased. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999, */
/*     Oct. 2000. */
/*     P.Hr. Petkov, October 2000. */
/*     A. Varga, October 2000. */
/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, system norm. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --iwork;
    --dwork;
    --cwork;
    --bwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*np < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldc < max(1,*np)) {
	*info = -9;
    } else if (*ldd < max(1,*np)) {
	*info = -11;
    }

/*     Compute workspace. */

/* Computing MAX */
    i__1 = 2, i__2 = (*n << 2) * *n + (*m << 1) * *m + *m * 3 * *n + *m * *np 
	    + (*n + *np << 1) * *np + *n * 10 + max(*m,*np) * 6;
    minwrk = max(i__1,i__2);
    if (*ldwork < minwrk) {
	*info = -15;
    }
/* Computing MAX */
    i__1 = 1, i__2 = (*n + *m) * (*n + *np) + max(*m,*np) * 3;
    mincwr = max(i__1,i__2);
    if (*lcwork < mincwr) {
	*info = -17;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("AB13CD", &i__1, (ftnlen)6);
	return ret_val;
    }

/*     Quick return if possible. */

    if (*m == 0 || *np == 0) {
	return ret_val;
    }

/*     Workspace usage. */

    iw2 = *n;
    iw3 = iw2 + *n;
    iw4 = iw3 + *n * *n;
    iw5 = iw4 + *n * *m;
    iw6 = iw5 + *np * *m;
    iwrk = iw6 + min(*np,*m);

/*     Determine the maximum singular value of G(infinity) = D . */

    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iw5 + 1], np, (ftnlen)
	    4);
    i__1 = *ldwork - iwrk;
    dgesvd_("N", "N", np, m, &dwork[iw5 + 1], np, &dwork[iw6 + 1], &dwork[1], 
	    np, &dwork[1], m, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
    if (info2 > 0) {
	*info = 4;
	return ret_val;
    }
    gammal = dwork[iw6 + 1];
    fpeak = 1e30;
    lwamax = (integer) dwork[iwrk + 1] + iwrk;

/*     Quick return if N = 0 . */

    if (*n == 0) {
	ret_val = gammal;
	dwork[1] = 2.;
	dwork[2] = 0.;
	cwork[1].r = 1., cwork[1].i = 0.;
	return ret_val;
    }

/*     Stability check. */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw3 + 1], n, (ftnlen)4);
    i__1 = *ldwork - iwrk;
    dgees_("N", "S", (L_fp)sb02mv_, n, &dwork[iw3 + 1], n, &sdim, &dwork[1], &
	    dwork[iw2 + 1], &dwork[1], n, &dwork[iwrk + 1], &i__1, &bwork[1], 
	    &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 3;
	return ret_val;
    }
    if (sdim < *n) {
	*info = 1;
	return ret_val;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(i__1,lwamax);

/*     Determine the maximum singular value of G(0) = -C*inv(A)*B + D . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw3 + 1], n, (ftnlen)4);
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[iw4 + 1], n, (ftnlen)4);
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iw5 + 1], np, (ftnlen)
	    4);
    dgesv_(n, m, &dwork[iw3 + 1], n, &iwork[1], &dwork[iw4 + 1], n, &info2);
    if (info2 > 0) {
	*info = 1;
	return ret_val;
    }
    dgemm_("N", "N", np, m, n, &c_b16, &c__[c_offset], ldc, &dwork[iw4 + 1], 
	    n, &c_b17, &dwork[iw5 + 1], np, (ftnlen)1, (ftnlen)1);
    i__1 = *ldwork - iwrk;
    dgesvd_("N", "N", np, m, &dwork[iw5 + 1], np, &dwork[iw6 + 1], &dwork[1], 
	    np, &dwork[1], m, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
    if (info2 > 0) {
	*info = 4;
	return ret_val;
    }
    if (gammal < dwork[iw6 + 1]) {
	gammal = dwork[iw6 + 1];
	fpeak = 0.;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(i__1,lwamax);

/*     Find a frequency which is close to the peak frequency. */

    complx = FALSE_;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dwork[iw2 + i__] != 0.) {
	    complx = TRUE_;
	}
/* L10: */
    }
    if (! complx) {
	wrmin = abs(dwork[1]);
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if (wrmin > (d__1 = dwork[i__], abs(d__1))) {
		wrmin = (d__2 = dwork[i__], abs(d__2));
	    }
/* L20: */
	}
	omega = wrmin;
    } else {
	ratmax = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    den = dlapy2_(&dwork[i__], &dwork[iw2 + i__]);
	    rat = (d__1 = dwork[iw2 + i__] / dwork[i__] / den, abs(d__1));
	    if (ratmax < rat) {
		ratmax = rat;
		wimax = den;
	    }
/* L30: */
	}
	omega = wimax;
    }

/*     Workspace usage. */

    icw2 = *n * *n;
    icw3 = icw2 + *n * *m;
    icw4 = icw3 + *np * *n;
    icwrk = icw4 + *np * *m;

/*     Determine the maximum singular value of */
/*     G(omega) = C*inv(j*omega*In - A)*B + D . */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + (j - 1) * *n;
	    d__1 = -a[i__ + j * a_dim1];
	    cwork[i__3].r = d__1, cwork[i__3].i = 0.;
/* L40: */
	}
	i__2 = j + (j - 1) * *n;
	z__2.r = omega * 0., z__2.i = omega * 1.;
	i__3 = j + j * a_dim1;
	z__1.r = z__2.r - a[i__3], z__1.i = z__2.i;
	cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
/* L50: */
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = icw2 + i__ + (j - 1) * *n;
	    i__4 = i__ + j * b_dim1;
	    cwork[i__3].r = b[i__4], cwork[i__3].i = 0.;
/* L60: */
	}
/* L70: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = icw3 + i__ + (j - 1) * *np;
	    i__4 = i__ + j * c_dim1;
	    cwork[i__3].r = c__[i__4], cwork[i__3].i = 0.;
/* L80: */
	}
/* L90: */
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = icw4 + i__ + (j - 1) * *np;
	    i__4 = i__ + j * d_dim1;
	    cwork[i__3].r = d__[i__4], cwork[i__3].i = 0.;
/* L100: */
	}
/* L110: */
    }
    zgesv_(n, m, &cwork[1], n, &iwork[1], &cwork[icw2 + 1], n, &info2);
    if (info2 > 0) {
	*info = 1;
	return ret_val;
    }
    zgemm_("N", "N", np, m, n, &c_b1, &cwork[icw3 + 1], np, &cwork[icw2 + 1], 
	    n, &c_b1, &cwork[icw4 + 1], np, (ftnlen)1, (ftnlen)1);
    i__1 = *lcwork - icwrk;
    zgesvd_("N", "N", np, m, &cwork[icw4 + 1], np, &dwork[iw6 + 1], &cwork[1],
	     np, &cwork[1], m, &cwork[icwrk + 1], &i__1, &dwork[iwrk + 1], &
	    info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 4;
	return ret_val;
    }
    if (gammal < dwork[iw6 + 1]) {
	gammal = dwork[iw6 + 1];
	fpeak = omega;
    }
    i__1 = icwrk + 1;
    lcwamx = (integer) cwork[i__1].r + icwrk;

/*     Workspace usage. */

    iw2 = *m * *n;
    iw3 = iw2 + *m * *m;
    iw4 = iw3 + *np * *np;
    iw5 = iw4 + *m * *m;
    iw6 = iw5 + *m * *n;
    iw7 = iw6 + *m * *n;
    iw8 = iw7 + *np * *np;
    iw9 = iw8 + *np * *n;
    iw10 = iw9 + (*n << 2) * *n;
    iw11 = iw10 + (*n << 1);
    iw12 = iw11 + (*n << 1);
    iwrk = iw12 + min(*np,*m);

/*     Compute D'*C . */

    dgemm_("T", "N", m, n, np, &c_b17, &d__[d_offset], ldd, &c__[c_offset], 
	    ldc, &c_b38, &dwork[1], m, (ftnlen)1, (ftnlen)1);

/*     Compute D'*D . */

    dsyrk_("U", "T", m, np, &c_b17, &d__[d_offset], ldd, &c_b38, &dwork[iw2 + 
	    1], m, (ftnlen)1, (ftnlen)1);

/*     Compute D*D' . */

    dsyrk_("U", "N", np, m, &c_b17, &d__[d_offset], ldd, &c_b38, &dwork[iw3 + 
	    1], np, (ftnlen)1, (ftnlen)1);

/*     Main iteration loop for gamma. */

    iter = 0;
L120:
    ++iter;
    if (iter > 10) {
	*info = 2;
	return ret_val;
    }
    gamma = (*tol * 2. + 1.) * gammal;

/*     Compute R = GAMMA^2*Im - D'*D . */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw4 + i__ + (j - 1) * *m] = -dwork[iw2 + i__ + (j - 1) * *m]
		    ;
/* L130: */
	}
/* Computing 2nd power */
	d__1 = gamma;
	dwork[iw4 + j + (j - 1) * *m] = d__1 * d__1 - dwork[iw2 + j + (j - 1) 
		* *m];
/* L140: */
    }

/*     Compute inv(R)*D'*C . */

    dlacpy_("Full", m, n, &dwork[1], m, &dwork[iw5 + 1], m, (ftnlen)4);
    dpotrf_("U", m, &dwork[iw4 + 1], m, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 2;
	return ret_val;
    }
    dpotrs_("U", m, n, &dwork[iw4 + 1], m, &dwork[iw5 + 1], m, &info2, (
	    ftnlen)1);

/*     Compute inv(R)*B' . */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw6 + i__ + (j - 1) * *m] = b[j + i__ * b_dim1];
/* L150: */
	}
/* L160: */
    }
    dpotrs_("U", m, n, &dwork[iw4 + 1], m, &dwork[iw6 + 1], m, &info2, (
	    ftnlen)1);

/*     Compute S = GAMMA^2*Ip - D*D' . */

    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw7 + i__ + (j - 1) * *np] = -dwork[iw3 + i__ + (j - 1) * *
		    np];
/* L170: */
	}
/* Computing 2nd power */
	d__1 = gamma;
	dwork[iw7 + j + (j - 1) * *np] = d__1 * d__1 - dwork[iw3 + j + (j - 1)
		 * *np];
/* L180: */
    }

/*     Compute inv(S)*C . */

    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iw8 + 1], np, (ftnlen)
	    4);
    dposv_("U", np, n, &dwork[iw7 + 1], np, &dwork[iw8 + 1], np, &info2, (
	    ftnlen)1);
    if (info2 > 0) {
	*info = 2;
	return ret_val;
    }

/*     Construct the Hamiltonian matrix . */

    i__1 = *n << 1;
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw9 + 1], &i__1, (ftnlen)
	    4);
    i__1 = *n << 1;
    dgemm_("N", "N", n, n, m, &c_b17, &b[b_offset], ldb, &dwork[iw5 + 1], m, &
	    c_b17, &dwork[iw9 + 1], &i__1, (ftnlen)1, (ftnlen)1);
    d__1 = -gamma;
    i__1 = *n << 1;
    mb01rx_("Left", "Upper", "Transpose", n, np, &c_b38, &d__1, &dwork[iw9 + *
	    n + 1], &i__1, &c__[c_offset], ldc, &dwork[iw8 + 1], np, &info2, (
	    ftnlen)4, (ftnlen)5, (ftnlen)9);
    i__1 = *n << 1;
    ma02ed_("Upper", n, &dwork[iw9 + *n + 1], &i__1, (ftnlen)5);
    i__1 = *n << 1;
    mb01rx_("Left", "Upper", "NoTranspose", n, m, &c_b38, &gamma, &dwork[iw9 
	    + (*n << 1) * *n + 1], &i__1, &b[b_offset], ldb, &dwork[iw6 + 1], 
	    m, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)11);
    i__1 = *n << 1;
    ma02ed_("Upper", n, &dwork[iw9 + (*n << 1) * *n + 1], &i__1, (ftnlen)5);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw9 + (*n << 1) * *n + *n + i__ + (j - 1 << 1) * *n] = 
		    -dwork[iw9 + j + (i__ - 1 << 1) * *n];
/* L190: */
	}
/* L200: */
    }

/*     Compute the eigenvalues of the Hamiltonian matrix. */

    i__1 = *n << 1;
    i__2 = *n << 1;
    i__3 = *n << 1;
    i__4 = *ldwork - iwrk;
    dgees_("N", "S", (L_fp)sb02cx_, &i__1, &dwork[iw9 + 1], &i__2, &sdim, &
	    dwork[iw10 + 1], &dwork[iw11 + 1], &dwork[1], &i__3, &dwork[iwrk 
	    + 1], &i__4, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 3;
	return ret_val;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(i__1,lwamax);

    if (sdim == 0) {
	gammau = gamma;
	goto L330;
    }

/*     Store the positive imaginary parts. */

    j = 0;
    i__1 = sdim - 1;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	++j;
	dwork[iw10 + j] = dwork[iw11 + i__];
/* L210: */
    }
    k = j;

    if (k >= 2) {

/*        Reorder the imaginary parts. */

	i__1 = k - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = k;
	    for (l = j + 1; l <= i__2; ++l) {
		if (dwork[iw10 + j] <= dwork[iw10 + l]) {
		    goto L220;
		}
		temp = dwork[iw10 + j];
		dwork[iw10 + j] = dwork[iw10 + l];
		dwork[iw10 + l] = temp;
L220:
		;
	    }
/* L230: */
	}

/*        Determine the next frequency. */

	i__1 = k - 1;
	for (l = 1; l <= i__1; ++l) {
	    omega = (dwork[iw10 + l] + dwork[iw10 + l + 1]) / 2.;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = i__ + (j - 1) * *n;
		    d__1 = -a[i__ + j * a_dim1];
		    cwork[i__4].r = d__1, cwork[i__4].i = 0.;
/* L240: */
		}
		i__3 = j + (j - 1) * *n;
		z__2.r = omega * 0., z__2.i = omega * 1.;
		i__4 = j + j * a_dim1;
		z__1.r = z__2.r - a[i__4], z__1.i = z__2.i;
		cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
/* L250: */
	    }
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = icw2 + i__ + (j - 1) * *n;
		    i__5 = i__ + j * b_dim1;
		    cwork[i__4].r = b[i__5], cwork[i__4].i = 0.;
/* L260: */
		}
/* L270: */
	    }
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *np;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = icw3 + i__ + (j - 1) * *np;
		    i__5 = i__ + j * c_dim1;
		    cwork[i__4].r = c__[i__5], cwork[i__4].i = 0.;
/* L280: */
		}
/* L290: */
	    }
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *np;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = icw4 + i__ + (j - 1) * *np;
		    i__5 = i__ + j * d_dim1;
		    cwork[i__4].r = d__[i__5], cwork[i__4].i = 0.;
/* L300: */
		}
/* L310: */
	    }
	    zgesv_(n, m, &cwork[1], n, &iwork[1], &cwork[icw2 + 1], n, &info2)
		    ;
	    if (info2 > 0) {
		*info = 1;
		return ret_val;
	    }
	    zgemm_("N", "N", np, m, n, &c_b1, &cwork[icw3 + 1], np, &cwork[
		    icw2 + 1], n, &c_b1, &cwork[icw4 + 1], np, (ftnlen)1, (
		    ftnlen)1);
	    i__2 = *lcwork - icwrk;
	    zgesvd_("N", "N", np, m, &cwork[icw4 + 1], np, &dwork[iw6 + 1], &
		    cwork[1], np, &cwork[1], m, &cwork[icwrk + 1], &i__2, &
		    dwork[iwrk + 1], &info2, (ftnlen)1, (ftnlen)1);
	    if (info2 > 0) {
		*info = 4;
		return ret_val;
	    }
	    if (gammal < dwork[iw6 + 1]) {
		gammal = dwork[iw6 + 1];
		fpeak = omega;
	    }
/* Computing MAX */
	    i__3 = icwrk + 1;
	    i__2 = (integer) cwork[i__3].r + icwrk;
	    lcwamx = max(i__2,lcwamx);
/* L320: */
	}
    }
    goto L120;
L330:
    ret_val = (gammal + gammau) / 2.;

    dwork[1] = (doublereal) lwamax;
    dwork[2] = fpeak;
    cwork[1].r = (doublereal) lcwamx, cwork[1].i = 0.;
    return ret_val;
/* *** End of AB13CD *** */
} /* ab13cd_ */

