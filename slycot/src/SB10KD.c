/* SB10KD.f -- translated by f2c (version 20100827).
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

static doublereal c_b5 = 1.;
static doublereal c_b6 = 0.;
static doublereal c_b17 = -1.;

/* Subroutine */ int sb10kd_(integer *n, integer *m, integer *np, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *factor, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *rcond, integer *iwork, 
	doublereal *dwork, integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, dk_dim1, 
	    dk_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, i2, i3, i4, i5, i6, i7, i8, i9, n2, i10, i11, 
	    i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, 
	    i25, i26, ns, lwa, sdim, iwrk, info2;
    static doublereal gamma;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), sb02od_(char *, char *, char *, char *, char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), 
	    dsyev_(char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen),
	     dsyrk_(char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal rnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern logical select_();
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dsytrf_(char *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), dsytrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);


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

/*     To compute the matrices of the positive feedback controller */

/*              | Ak | Bk | */
/*          K = |----|----| */
/*              | Ck | Dk | */

/*     for the shaped plant */

/*              | A | B | */
/*          G = |---|---| */
/*              | C | 0 | */

/*     in the Discrete-Time Loop Shaping Design Procedure. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the plant.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A of the shaped plant. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B of the shaped plant. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C of the shaped plant. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     FACTOR  (input) DOUBLE PRECISION */
/*             = 1  implies that an optimal controller is required; */
/*             > 1  implies that a suboptimal controller is required */
/*                  achieving a performance FACTOR less than optimal. */
/*             FACTOR >= 1. */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix Ak. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NP) */
/*             The leading N-by-NP part of this array contains the */
/*             controller input matrix Bk. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading M-by-N part of this array contains the */
/*             controller output matrix Ck. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK.  LDCK >= max(1,M). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NP) */
/*             The leading M-by-NP part of this array contains the */
/*             controller matrix Dk. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK.  LDDK >= max(1,M). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*             RCOND(1) contains an estimate of the reciprocal condition */
/*                      number of the linear system of equations from */
/*                      which the solution of the P-Riccati equation is */
/*                      obtained; */
/*             RCOND(2) contains an estimate of the reciprocal condition */
/*                      number of the linear system of equations from */
/*                      which the solution of the Q-Riccati equation is */
/*                      obtained; */
/*             RCOND(3) contains an estimate of the reciprocal condition */
/*                      number of the linear system of equations from */
/*                      which the solution of the X-Riccati equation is */
/*                      obtained; */
/*             RCOND(4) contains an estimate of the reciprocal condition */
/*                      number of the matrix Rx + Bx'*X*Bx (see the */
/*                      comments in the code). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*max(N,NP+M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 15*N*N + 6*N + */
/*                       max( 14*N+23, 16*N, 2*N+NP+M, 3*(NP+M) ) + */
/*                       max( N*N, 11*N*NP + 2*M*M + 8*NP*NP + 8*M*N + */
/*                                 4*M*NP + NP ). */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the P-Riccati equation is not solved successfully; */
/*             = 2:  the Q-Riccati equation is not solved successfully; */
/*             = 3:  the X-Riccati equation is not solved successfully; */
/*             = 4:  the iteration to compute eigenvalues failed to */
/*                   converge; */
/*             = 5:  the matrix Rx + Bx'*X*Bx is singular; */
/*             = 6:  the closed-loop system is unstable. */

/*     METHOD */

/*     The routine implements the method presented in [1]. */

/*     REFERENCES */

/*     [1] McFarlane, D. and Glover, K. */
/*         A loop shaping design procedure using H_infinity synthesis. */
/*         IEEE Trans. Automat. Control, vol. AC-37, no. 6, pp. 759-769, */
/*         1992. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the results depends on the conditioning of the */
/*     two Riccati equations solved in the controller design. For */
/*     better conditioning it is advised to take FACTOR > 1. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke University Leuven, January 2001, */
/*     February 2001. */

/*     KEYWORDS */

/*     H_infinity control, Loop-shaping design, Robust control. */

/*     ****************************************************************** */

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

/*     Decode and Test input parameters. */

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
    ak_dim1 = *ldak;
    ak_offset = 1 + ak_dim1;
    ak -= ak_offset;
    bk_dim1 = *ldbk;
    bk_offset = 1 + bk_dim1;
    bk -= bk_offset;
    ck_dim1 = *ldck;
    ck_offset = 1 + ck_dim1;
    ck -= ck_offset;
    dk_dim1 = *lddk;
    dk_offset = 1 + dk_dim1;
    dk -= dk_offset;
    --rcond;
    --iwork;
    --dwork;
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
    } else if (*factor < 1.) {
	*info = -10;
    } else if (*ldak < max(1,*n)) {
	*info = -12;
    } else if (*ldbk < max(1,*n)) {
	*info = -14;
    } else if (*ldck < max(1,*m)) {
	*info = -16;
    } else if (*lddk < max(1,*m)) {
	*info = -18;
    }

/*     Compute workspace. */

/* Computing MAX */
    i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2), i__2 = (*n << 
	    1) + *np + *m, i__1 = max(i__1,i__2), i__2 = (*np + *m) * 3;
/* Computing MAX */
    i__3 = *n * *n, i__4 = *n * 11 * *np + (*m << 1) * *m + (*np << 3) * *np 
	    + (*m << 3) * *n + (*m << 2) * *np + *np;
    minwrk = *n * 15 * *n + *n * 6 + max(i__1,i__2) + max(i__3,i__4);
    if (*ldwork < minwrk) {
	*info = -22;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10KD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0) {
	rcond[1] = 1.;
	rcond[2] = 1.;
	rcond[3] = 1.;
	rcond[4] = 1.;
	dwork[1] = 1.;
	return 0;
    }

/*     Workspace usage. */

    n2 = *n << 1;
    i1 = *n * *n;
    i2 = i1 + *n * *n;
    i3 = i2 + *n * *n;
    i4 = i3 + *n * *n;
    i5 = i4 + n2;
    i6 = i5 + n2;
    i7 = i6 + n2;
    i8 = i7 + n2 * n2;
    i9 = i8 + n2 * n2;

    iwrk = i9 + n2 * n2;
    lwamax = 0;

/*     Compute Cr = C'*C . */

    dsyrk_("U", "T", n, np, &c_b5, &c__[c_offset], ldc, &c_b6, &dwork[i2 + 1],
	     n, (ftnlen)1, (ftnlen)1);

/*     Compute Dr = B*B' . */

    dsyrk_("U", "N", n, m, &c_b5, &b[b_offset], ldb, &c_b6, &dwork[i3 + 1], n,
	     (ftnlen)1, (ftnlen)1);
/*                                                     -1 */
/*     Solution of the Riccati equation A'*P*(In + Dr*P) *A - P + Cr = 0. */

    i__1 = *ldwork - iwrk;
    sb02od_("D", "G", "N", "U", "Z", "S", n, m, np, &a[a_offset], lda, &dwork[
	    i3 + 1], n, &dwork[i2 + 1], n, &dwork[1], m, &dwork[1], n, &rcond[
	    1], &dwork[1], n, &dwork[i4 + 1], &dwork[i5 + 1], &dwork[i6 + 1], 
	    &dwork[i7 + 1], &n2, &dwork[i8 + 1], &n2, &dwork[i9 + 1], &n2, &
	    c_b17, &iwork[1], &dwork[iwrk + 1], &i__1, &bwork[1], &info2, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 1;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);

/*     Transpose A in AK (used as workspace). */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ak[i__ + j * ak_dim1] = a[j + i__ * a_dim1];
/* L30: */
	}
/* L40: */
    }
/*                                                    -1 */
/*     Solution of the Riccati equation A*Q*(In + Cr*Q) *A' - Q + Dr = 0. */

    i__1 = *ldwork - iwrk;
    sb02od_("D", "G", "N", "U", "Z", "S", n, m, np, &ak[ak_offset], ldak, &
	    dwork[i2 + 1], n, &dwork[i3 + 1], n, &dwork[1], m, &dwork[1], n, &
	    rcond[2], &dwork[i1 + 1], n, &dwork[i4 + 1], &dwork[i5 + 1], &
	    dwork[i6 + 1], &dwork[i7 + 1], &n2, &dwork[i8 + 1], &n2, &dwork[
	    i9 + 1], &n2, &c_b17, &iwork[1], &dwork[iwrk + 1], &i__1, &bwork[
	    1], &info2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1,
	     (ftnlen)1);
    if (info2 != 0) {
	*info = 2;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);

/*     Compute gamma. */

    dgemm_("N", "N", n, n, n, &c_b5, &dwork[i1 + 1], n, &dwork[1], n, &c_b6, &
	    ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
    i__1 = *ldwork - iwrk;
    dgees_("N", "N", (L_fp)select_, n, &ak[ak_offset], ldak, &sdim, &dwork[i6 
	    + 1], &dwork[i7 + 1], &dwork[iwrk + 1], n, &dwork[iwrk + 1], &
	    i__1, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 4;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);
    gamma = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = gamma, d__2 = dwork[i6 + i__];
	gamma = max(d__1,d__2);
/* L50: */
    }
    gamma = *factor * sqrt(gamma + 1.);

/*     Workspace usage. */

    i3 = i2 + *n * *np;
    i4 = i3 + *np * *np;
    i5 = i4 + *np * *np;
    i6 = i5 + *np * *np;
    i7 = i6 + *np;
    i8 = i7 + *np * *np;
    i9 = i8 + *np * *np;
    i10 = i9 + *np * *np;
    i11 = i10 + *n * *np;
    i12 = i11 + *n * *np;
    i13 = i12 + (*np + *m) * (*np + *m);
    i14 = i13 + *n * (*np + *m);
    i15 = i14 + *n * (*np + *m);
    i16 = i15 + *n * *n;
    i17 = i16 + n2;
    i18 = i17 + n2;
    i19 = i18 + n2;
    i20 = i19 + (n2 + *np + *m) * (n2 + *np + *m);
    i21 = i20 + (n2 + *np + *m) * n2;

    iwrk = i21 + n2 * n2;

/*     Compute Q*C' . */

    dgemm_("N", "T", n, np, n, &c_b5, &dwork[i1 + 1], n, &c__[c_offset], ldc, 
	    &c_b6, &dwork[i2 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Compute Ip + C*Q*C' . */

    dlaset_("Full", np, np, &c_b6, &c_b5, &dwork[i3 + 1], np, (ftnlen)4);
    dgemm_("N", "N", np, np, n, &c_b5, &c__[c_offset], ldc, &dwork[i2 + 1], n,
	     &c_b5, &dwork[i3 + 1], np, (ftnlen)1, (ftnlen)1);

/*     Compute the eigenvalues and eigenvectors of Ip + C'*Q*C */

    dlacpy_("U", np, np, &dwork[i3 + 1], np, &dwork[i5 + 1], np, (ftnlen)1);
    i__1 = *ldwork - iwrk;
    dsyev_("V", "U", np, &dwork[i5 + 1], np, &dwork[i6 + 1], &dwork[iwrk + 1],
	     &i__1, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 4;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);
/*                            -1 */
/*     Compute ( Ip + C'*Q*C )  . */

    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i9 + i__ + (j - 1) * *np] = dwork[i5 + j + (i__ - 1) * *np] 
		    / dwork[i6 + i__];
/* L60: */
	}
/* L70: */
    }
    dgemm_("N", "N", np, np, np, &c_b5, &dwork[i5 + 1], np, &dwork[i9 + 1], 
	    np, &c_b6, &dwork[i4 + 1], np, (ftnlen)1, (ftnlen)1);

/*     Compute Z2 . */

    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i9 + i__ + (j - 1) * *np] = dwork[i5 + j + (i__ - 1) * *np] 
		    / sqrt(dwork[i6 + i__]);
/* L80: */
	}
/* L90: */
    }
    dgemm_("N", "N", np, np, np, &c_b5, &dwork[i5 + 1], np, &dwork[i9 + 1], 
	    np, &c_b6, &dwork[i7 + 1], np, (ftnlen)1, (ftnlen)1);
/*               -1 */
/*     Compute Z2  . */

    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i9 + i__ + (j - 1) * *np] = dwork[i5 + j + (i__ - 1) * *np] 
		    * sqrt(dwork[i6 + i__]);
/* L100: */
	}
/* L110: */
    }
    dgemm_("N", "N", np, np, np, &c_b5, &dwork[i5 + 1], np, &dwork[i9 + 1], 
	    np, &c_b6, &dwork[i8 + 1], np, (ftnlen)1, (ftnlen)1);

/*     Compute A*Q*C' . */

    dgemm_("N", "N", n, np, n, &c_b5, &a[a_offset], lda, &dwork[i2 + 1], n, &
	    c_b6, &dwork[i10 + 1], n, (ftnlen)1, (ftnlen)1);
/*                                        -1 */
/*     Compute H = -A*Q*C'*( Ip + C*Q*C' )  . */

    dgemm_("N", "N", n, np, np, &c_b17, &dwork[i10 + 1], n, &dwork[i4 + 1], 
	    np, &c_b6, &dwork[i11 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Compute Rx . */

    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *np + *m;
    dlaset_("F", &i__1, &i__2, &c_b6, &c_b5, &dwork[i12 + 1], &i__3, (ftnlen)
	    1);
    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i12 + i__ + (j - 1) * (*np + *m)] = dwork[i3 + i__ + (j - 1)
		     * *np];
/* L120: */
	}
	dwork[i12 + j + (j - 1) * (*np + *m)] = dwork[i3 + j + (j - 1) * *np] 
		- gamma * gamma;
/* L130: */
    }

/*     Compute Bx . */

    dgemm_("N", "N", n, np, np, &c_b17, &dwork[i11 + 1], n, &dwork[i8 + 1], 
	    np, &c_b6, &dwork[i13 + 1], n, (ftnlen)1, (ftnlen)1);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i13 + *n * *np + i__ + (j - 1) * *n] = b[i__ + j * b_dim1];
/* L140: */
	}
/* L150: */
    }

/*     Compute Sx . */

    dgemm_("T", "N", n, np, np, &c_b5, &c__[c_offset], ldc, &dwork[i8 + 1], 
	    np, &c_b6, &dwork[i14 + 1], n, (ftnlen)1, (ftnlen)1);
    dlaset_("F", n, m, &c_b6, &c_b6, &dwork[i14 + *n * *np + 1], n, (ftnlen)1)
	    ;

/*     Solve the Riccati equation */
/*                                                      -1 */
/*       X = A'*X*A + Cx - (Sx + A'*X*Bx)*(Rx + Bx'*X*B ) *(Sx'+Bx'*X*A). */

    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = n2 + *np + *m;
    i__4 = n2 + *np + *m;
    i__5 = *ldwork - iwrk;
    sb02od_("D", "B", "C", "U", "N", "S", n, &i__1, np, &a[a_offset], lda, &
	    dwork[i13 + 1], n, &c__[c_offset], ldc, &dwork[i12 + 1], &i__2, &
	    dwork[i14 + 1], n, &rcond[3], &dwork[i15 + 1], n, &dwork[i16 + 1],
	     &dwork[i17 + 1], &dwork[i18 + 1], &dwork[i19 + 1], &i__3, &dwork[
	    i20 + 1], &i__4, &dwork[i21 + 1], &n2, &c_b17, &iwork[1], &dwork[
	    iwrk + 1], &i__5, &bwork[1], &info2, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 3;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);

    i22 = i16;
    i23 = i22 + (*np + *m) * *n;
    i24 = i23 + (*np + *m) * (*np + *m);
    i25 = i24 + (*np + *m) * *n;
    i26 = i25 + *m * *n;

    iwrk = i25;

/*     Compute Bx'*X . */

    i__1 = *np + *m;
    i__2 = *np + *m;
    dgemm_("T", "N", &i__1, n, n, &c_b5, &dwork[i13 + 1], n, &dwork[i15 + 1], 
	    n, &c_b6, &dwork[i22 + 1], &i__2, (ftnlen)1, (ftnlen)1);

/*     Compute Rx + Bx'*X*Bx . */

    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *np + *m;
    i__4 = *np + *m;
    dlacpy_("F", &i__1, &i__2, &dwork[i12 + 1], &i__3, &dwork[i23 + 1], &i__4,
	     (ftnlen)1);
    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *np + *m;
    i__4 = *np + *m;
    dgemm_("N", "N", &i__1, &i__2, n, &c_b5, &dwork[i22 + 1], &i__3, &dwork[
	    i13 + 1], n, &c_b5, &dwork[i23 + 1], &i__4, (ftnlen)1, (ftnlen)1);

/*     Compute -( Sx' + Bx'*X*A ) . */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *np + *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i24 + i__ + (j - 1) * (*np + *m)] = dwork[i14 + j + (i__ - 
		    1) * *n];
/* L160: */
	}
/* L170: */
    }
    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *np + *m;
    dgemm_("N", "N", &i__1, n, n, &c_b17, &dwork[i22 + 1], &i__2, &a[a_offset]
	    , lda, &c_b17, &dwork[i24 + 1], &i__3, (ftnlen)1, (ftnlen)1);

/*     Factorize Rx + Bx'*X*Bx . */

    i__1 = *np + *m;
    i__2 = *np + *m;
    rnorm = dlansy_("1", "U", &i__1, &dwork[i23 + 1], &i__2, &dwork[iwrk + 1],
	     (ftnlen)1, (ftnlen)1);
    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *ldwork - iwrk;
    dsytrf_("U", &i__1, &dwork[i23 + 1], &i__2, &iwork[1], &dwork[iwrk + 1], &
	    i__3, &info2, (ftnlen)1);
    if (info2 != 0) {
	*info = 5;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);
    i__1 = *np + *m;
    i__2 = *np + *m;
    dsycon_("U", &i__1, &dwork[i23 + 1], &i__2, &iwork[1], &rnorm, &rcond[4], 
	    &dwork[iwrk + 1], &iwork[*np + *m + 1], &info2, (ftnlen)1);
/*                                   -1 */
/*     Compute F = -( Rx + Bx'*X*Bx )  ( Sx' + Bx'*X*A ) . */

    i__1 = *np + *m;
    i__2 = *np + *m;
    i__3 = *np + *m;
    dsytrs_("U", &i__1, n, &dwork[i23 + 1], &i__2, &iwork[1], &dwork[i24 + 1],
	     &i__3, &info2, (ftnlen)1);

/*     Compute B'*X . */

    dgemm_("T", "N", m, n, n, &c_b5, &b[b_offset], ldb, &dwork[i15 + 1], n, &
	    c_b6, &dwork[i25 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute Im + B'*X*B . */

    dlaset_("F", m, m, &c_b6, &c_b5, &dwork[i23 + 1], m, (ftnlen)1);
    dgemm_("N", "N", m, m, n, &c_b5, &dwork[i25 + 1], m, &b[b_offset], ldb, &
	    c_b5, &dwork[i23 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Factorize Im + B'*X*B . */

    dpotrf_("U", m, &dwork[i23 + 1], m, &info2, (ftnlen)1);
/*                            -1 */
/*     Compute ( Im + B'*X*B )  B'*X . */

    dpotrs_("U", m, n, &dwork[i23 + 1], m, &dwork[i25 + 1], m, &info2, (
	    ftnlen)1);
/*                                 -1 */
/*     Compute Dk = ( Im + B'*X*B )  B'*X*H . */

    dgemm_("N", "N", m, np, n, &c_b5, &dwork[i25 + 1], m, &dwork[i11 + 1], n, 
	    &c_b6, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute Bk = -H + B*Dk . */

    dlacpy_("F", n, np, &dwork[i11 + 1], n, &bk[bk_offset], ldbk, (ftnlen)1);
    dgemm_("N", "N", n, np, m, &c_b5, &b[b_offset], ldb, &dk[dk_offset], lddk,
	     &c_b17, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);
/*                  -1 */
/*     Compute Dk*Z2  . */

    dgemm_("N", "N", m, np, np, &c_b5, &dk[dk_offset], lddk, &dwork[i8 + 1], 
	    np, &c_b6, &dwork[i26 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute F1 + Z2*C . */

    i__1 = *np + *m;
    dlacpy_("F", np, n, &dwork[i24 + 1], &i__1, &dwork[i12 + 1], np, (ftnlen)
	    1);
    dgemm_("N", "N", np, n, np, &c_b5, &dwork[i7 + 1], np, &c__[c_offset], 
	    ldc, &c_b5, &dwork[i12 + 1], np, (ftnlen)1, (ftnlen)1);
/*                            -1 */
/*     Compute Ck = F2 - Dk*Z2  *( F1 + Z2*C ) . */

    i__1 = *np + *m;
    dlacpy_("F", m, n, &dwork[i24 + *np + 1], &i__1, &ck[ck_offset], ldck, (
	    ftnlen)1);
    dgemm_("N", "N", m, n, np, &c_b17, &dwork[i26 + 1], m, &dwork[i12 + 1], 
	    np, &c_b5, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);

/*     Compute Ak = A + H*C + B*Ck . */

    dlacpy_("F", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)1);
    dgemm_("N", "N", n, n, np, &c_b5, &dwork[i11 + 1], n, &c__[c_offset], ldc,
	     &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", n, n, m, &c_b5, &b[b_offset], ldb, &ck[ck_offset], ldck, 
	    &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Workspace usage. */

    i1 = *m * *n;
    i2 = i1 + n2 * n2;
    i3 = i2 + n2;

    iwrk = i3 + n2;

/*     Compute Dk*C . */

    dgemm_("N", "N", m, n, np, &c_b5, &dk[dk_offset], lddk, &c__[c_offset], 
	    ldc, &c_b6, &dwork[1], m, (ftnlen)1, (ftnlen)1);

/*     Compute the closed-loop state matrix. */

    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i1 + 1], &n2, (ftnlen)1);
    dgemm_("N", "N", n, n, m, &c_b17, &b[b_offset], ldb, &dwork[1], m, &c_b5, 
	    &dwork[i1 + 1], &n2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", n, n, np, &c_b17, &bk[bk_offset], ldbk, &c__[c_offset], 
	    ldc, &c_b6, &dwork[i1 + *n + 1], &n2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", n, n, m, &c_b5, &b[b_offset], ldb, &ck[ck_offset], ldck, 
	    &c_b6, &dwork[i1 + n2 * *n + 1], &n2, (ftnlen)1, (ftnlen)1);
    dlacpy_("F", n, n, &ak[ak_offset], ldak, &dwork[i1 + n2 * *n + *n + 1], &
	    n2, (ftnlen)1);

/*     Compute the closed-loop poles. */

    i__1 = *ldwork - iwrk;
    dgees_("N", "N", (L_fp)select_, &n2, &dwork[i1 + 1], &n2, &sdim, &dwork[
	    i2 + 1], &dwork[i3 + 1], &dwork[iwrk + 1], n, &dwork[iwrk + 1], &
	    i__1, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 4;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1] + iwrk;
    lwamax = max(lwa,lwamax);

/*     Check the stability of the closed-loop system. */

    ns = 0;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dlapy2_(&dwork[i2 + i__], &dwork[i3 + i__]) > 1.) {
	    ++ns;
	}
/* L180: */
    }
    if (ns > 0) {
	*info = 6;
	return 0;
    }

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10KD *** */
} /* sb10kd_ */

