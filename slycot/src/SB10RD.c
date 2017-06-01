/* SB10RD.f -- translated by f2c (version 20100827).
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

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static doublereal c_b19 = -1.;

/* Subroutine */ int sb10rd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *gamma, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *f, integer *ldf, 
	doublereal *h__, integer *ldh, doublereal *tu, integer *ldtu, 
	doublereal *ty, integer *ldty, doublereal *x, integer *ldx, 
	doublereal *y, integer *ldy, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, f_dim1, f_offset, h_dim1, h_offset, tu_dim1, 
	    tu_offset, ty_dim1, ty_offset, x_dim1, x_offset, y_dim1, y_offset,
	     i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m1, m2, ij, nd1, nd2, np1, np2, iw1, iw2, iw3, iw4,
	     id11, id12, id21, iwb, iwc;
    static doublereal eps;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb01rx_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dgetri_(integer *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen), dgetrs_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dsytrs_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);


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

/*     To compute the matrices of an H-infinity (sub)optimal controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     from the state feedback matrix F and output injection matrix H as */
/*     determined by the SLICOT Library routine SB10QD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0. */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0. */
/*             M-NCON >= NMEAS. */

/*     GAMMA   (input) DOUBLE PRECISION */
/*             The value of gamma. It is assumed that gamma is */
/*             sufficiently large so that the controller is admissible. */
/*             GAMMA >= 0. */

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

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array must contain the */
/*             state feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= max(1,M). */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,NP) */
/*             The leading N-by-NP part of this array must contain the */
/*             output injection matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     TU      (input) DOUBLE PRECISION array, dimension (LDTU,M2) */
/*             The leading M2-by-M2 part of this array must contain the */
/*             control transformation matrix TU, as obtained by the */
/*             SLICOT Library routine SB10PD. */

/*     LDTU    INTEGER */
/*             The leading dimension of the array TU.  LDTU >= max(1,M2). */

/*     TY      (input) DOUBLE PRECISION array, dimension (LDTY,NP2) */
/*             The leading NP2-by-NP2 part of this array must contain the */
/*             measurement transformation matrix TY, as obtained by the */
/*             SLICOT Library routine SB10PD. */

/*     LDTY    INTEGER */
/*             The leading dimension of the array TY. */
/*             LDTY >= max(1,NP2). */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix X, solution of the X-Riccati equation, as obtained */
/*             by the SLICOT Library routine SB10QD. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix Y, solution of the Y-Riccati equation, as obtained */
/*             by the SLICOT Library routine SB10QD. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= max(1,N). */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array contains the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array contains the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array contains the */
/*             controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK = max(2*(max(NP,M)-M2-NP2,M2,N),NP2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(1, M2*NP2 + NP2*NP2 + M2*M2 + */
/*                           max(D1*D1 + max(2*D1, (D1+D2)*NP2), */
/*                               D2*D2 + max(2*D2, D2*M2), 3*N, */
/*                               N*(2*NP2 + M2) + */
/*                               max(2*N*M2, M2*NP2 + */
/*                                           max(M2*M2+3*M2, NP2*(2*NP2+ */
/*                                                  M2+max(NP2,N)))))) */
/*             where D1 = NP1 - M2, D2 = M1 - NP2, */
/*                  NP1 = NP - NP2, M1 = M - M2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             max( 1, Q*(3*Q + 3*N + max(2*N, 4*Q + max(Q, N)))). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 2:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is zero. */

/*     METHOD */

/*     The routine implements the Glover's and Doyle's formulas [1],[2]. */

/*     REFERENCES */

/*     [1] Glover, K. and Doyle, J.C. */
/*         State-space formulae for all stabilizing controllers that */
/*         satisfy an Hinf norm bound and relations to risk sensitivity. */
/*         Systems and Control Letters, vol. 11, pp. 167-172, 1988. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Oct. 2001. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H-infinity optimal control, robust */
/*     control. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    tu_dim1 = *ldtu;
    tu_offset = 1 + tu_dim1;
    tu -= tu_offset;
    ty_dim1 = *ldty;
    ty_offset = 1 + ty_dim1;
    ty -= ty_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
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
    --iwork;
    --dwork;

    /* Function Body */
    m1 = *m - *ncon;
    m2 = *ncon;
    np1 = *np - *nmeas;
    np2 = *nmeas;

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*np < 0) {
	*info = -3;
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
	*info = -4;
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
	*info = -5;
    } else if (*gamma < 0.) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldc < max(1,*np)) {
	*info = -12;
    } else if (*ldd < max(1,*np)) {
	*info = -14;
    } else if (*ldf < max(1,*m)) {
	*info = -16;
    } else if (*ldh < max(1,*n)) {
	*info = -18;
    } else if (*ldtu < max(1,m2)) {
	*info = -20;
    } else if (*ldty < max(1,np2)) {
	*info = -22;
    } else if (*ldx < max(1,*n)) {
	*info = -24;
    } else if (*ldy < max(1,*n)) {
	*info = -26;
    } else if (*ldak < max(1,*n)) {
	*info = -28;
    } else if (*ldbk < max(1,*n)) {
	*info = -30;
    } else if (*ldck < max(1,m2)) {
	*info = -32;
    } else if (*lddk < max(1,m2)) {
	*info = -34;
    } else {

/*        Compute workspace. */

	nd1 = np1 - m2;
	nd2 = m1 - np2;
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	i__5 = nd1 << 1, i__6 = (nd1 + nd2) * np2;
/* Computing MAX */
	i__7 = nd2 << 1, i__8 = nd2 * m2;
/* Computing MAX */
/* Computing MAX */
	i__11 = m2 * m2 + m2 * 3, i__12 = np2 * ((np2 << 1) + m2 + max(np2,*n)
		);
	i__9 = (*n << 1) * m2, i__10 = m2 * np2 + max(i__11,i__12);
	i__3 = nd1 * nd1 + max(i__5,i__6), i__4 = nd2 * nd2 + max(i__7,i__8), 
		i__3 = max(i__3,i__4), i__4 = *n * 3, i__3 = max(i__3,i__4), 
		i__4 = *n * ((np2 << 1) + m2) + max(i__9,i__10);
	i__1 = 1, i__2 = m2 * np2 + np2 * np2 + m2 * m2 + max(i__3,i__4);
	minwrk = max(i__1,i__2);
	if (*ldwork < minwrk) {
	    *info = -37;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10RD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Get the machine precision. */

    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Workspace usage. */

    id11 = 1;
    id21 = id11 + m2 * np2;
    id12 = id21 + np2 * np2;
    iw1 = id12 + m2 * m2;
    iw2 = iw1 + nd1 * nd1;
    iw3 = iw2 + nd1 * np2;
    iwrk = iw2;

/*     Set D11HAT := -D1122 . */

    ij = id11;
    i__1 = np2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = m2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[ij] = -d__[nd1 + i__ + (nd2 + j) * d_dim1];
	    ++ij;
/* L10: */
	}
/* L20: */
    }

/*     Set D21HAT := Inp2 . */

    dlaset_("Upper", &np2, &np2, &c_b7, &c_b8, &dwork[id21], &np2, (ftnlen)5);

/*     Set D12HAT := Im2 . */

    dlaset_("Lower", &m2, &m2, &c_b7, &c_b8, &dwork[id12], &m2, (ftnlen)5);

/*     Compute D11HAT, D21HAT, D12HAT . */

    lwamax = 0;
    if (nd1 > 0) {
	if (nd2 == 0) {

/*           Compute D21HAT'*D21HAT = Inp2 - D1112'*D1112/gamma^2 . */

/* Computing 2nd power */
	    d__2 = *gamma;
	    d__1 = -1. / (d__2 * d__2);
	    dsyrk_("U", "T", &np2, &nd1, &d__1, &d__[d_offset], ldd, &c_b8, &
		    dwork[id21], &np2, (ftnlen)1, (ftnlen)1);
	} else {

/*           Compute gdum = gamma^2*Ind1 - D1111*D1111' . */

/* Computing 2nd power */
	    d__2 = *gamma;
	    d__1 = d__2 * d__2;
	    dlaset_("U", &nd1, &nd1, &c_b7, &d__1, &dwork[iw1], &nd1, (ftnlen)
		    1);
	    dsyrk_("U", "N", &nd1, &nd2, &c_b19, &d__[d_offset], ldd, &c_b8, &
		    dwork[iw1], &nd1, (ftnlen)1, (ftnlen)1);
	    anorm = dlansy_("I", "U", &nd1, &dwork[iw1], &nd1, &dwork[iwrk], (
		    ftnlen)1, (ftnlen)1);
	    i__1 = *ldwork - iwrk + 1;
	    dsytrf_("U", &nd1, &dwork[iw1], &nd1, &iwork[1], &dwork[iwrk], &
		    i__1, &info2, (ftnlen)1);
	    if (info2 > 0) {
		*info = 1;
		return 0;
	    }
	    lwamax = (integer) dwork[iwrk] + iwrk - 1;
	    dsycon_("U", &nd1, &dwork[iw1], &nd1, &iwork[1], &anorm, &rcond, &
		    dwork[iwrk], &iwork[nd1 + 1], &info2, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

	    if (rcond < eps) {
		*info = 1;
		return 0;
	    }

/*           Compute inv(gdum)*D1112 . */

	    dlacpy_("Full", &nd1, &np2, &d__[(nd2 + 1) * d_dim1 + 1], ldd, &
		    dwork[iw2], &nd1, (ftnlen)4);
	    dsytrs_("U", &nd1, &np2, &dwork[iw1], &nd1, &iwork[1], &dwork[iw2]
		    , &nd1, &info2, (ftnlen)1);

/*           Compute D11HAT = -D1121*D1111'*inv(gdum)*D1112 - D1122 . */

	    dgemm_("T", "N", &nd2, &np2, &nd1, &c_b8, &d__[d_offset], ldd, &
		    dwork[iw2], &nd1, &c_b7, &dwork[iw3], &nd2, (ftnlen)1, (
		    ftnlen)1);
	    dgemm_("N", "N", &m2, &np2, &nd2, &c_b19, &d__[nd1 + 1 + d_dim1], 
		    ldd, &dwork[iw3], &nd2, &c_b8, &dwork[id11], &m2, (ftnlen)
		    1, (ftnlen)1);

/*           Compute D21HAT'*D21HAT = Inp2 - D1112'*inv(gdum)*D1112 . */

	    mb01rx_("Left", "Upper", "Transpose", &np2, &nd1, &c_b8, &c_b19, &
		    dwork[id21], &np2, &d__[(nd2 + 1) * d_dim1 + 1], ldd, &
		    dwork[iw2], &nd1, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)9)
		    ;

	    iw2 = iw1 + nd2 * nd2;
	    iwrk = iw2;

/*           Compute gdum = gamma^2*Ind2 - D1111'*D1111 . */

/* Computing 2nd power */
	    d__2 = *gamma;
	    d__1 = d__2 * d__2;
	    dlaset_("L", &nd2, &nd2, &c_b7, &d__1, &dwork[iw1], &nd2, (ftnlen)
		    1);
	    dsyrk_("L", "T", &nd2, &nd1, &c_b19, &d__[d_offset], ldd, &c_b8, &
		    dwork[iw1], &nd2, (ftnlen)1, (ftnlen)1);
	    anorm = dlansy_("I", "L", &nd2, &dwork[iw1], &nd2, &dwork[iwrk], (
		    ftnlen)1, (ftnlen)1);
	    i__1 = *ldwork - iwrk + 1;
	    dsytrf_("L", &nd2, &dwork[iw1], &nd2, &iwork[1], &dwork[iwrk], &
		    i__1, &info2, (ftnlen)1);
	    if (info2 > 0) {
		*info = 1;
		return 0;
	    }
/* Computing MAX */
	    i__1 = (integer) dwork[iwrk] + iwrk - 1;
	    lwamax = max(i__1,lwamax);
	    dsycon_("L", &nd2, &dwork[iw1], &nd2, &iwork[1], &anorm, &rcond, &
		    dwork[iwrk], &iwork[nd2 + 1], &info2, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

	    if (rcond < eps) {
		*info = 1;
		return 0;
	    }

/*           Compute inv(gdum)*D1121' . */

	    ma02ad_("Full", &m2, &nd2, &d__[nd1 + 1 + d_dim1], ldd, &dwork[
		    iw2], &nd2, (ftnlen)4);
	    dsytrs_("L", &nd2, &m2, &dwork[iw1], &nd2, &iwork[1], &dwork[iw2],
		     &nd2, &info2, (ftnlen)1);

/*           Compute D12HAT*D12HAT' = Im2 - D1121*inv(gdum)*D1121' . */

	    mb01rx_("Left", "Lower", "NoTranspose", &m2, &nd2, &c_b8, &c_b19, 
		    &dwork[id12], &m2, &d__[nd1 + 1 + d_dim1], ldd, &dwork[
		    iw2], &nd2, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)11);
	}
    } else {
	if (nd2 > 0) {

/*           Compute D12HAT*D12HAT' = Im2 - D1121*D1121'/gamma^2 . */

/* Computing 2nd power */
	    d__2 = *gamma;
	    d__1 = -1. / (d__2 * d__2);
	    dsyrk_("L", "N", &m2, &nd2, &d__1, &d__[d_offset], ldd, &c_b8, &
		    dwork[id12], &m2, (ftnlen)1, (ftnlen)1);
	}
    }

/*     Compute D21HAT using Cholesky decomposition. */

    dpotrf_("U", &np2, &dwork[id21], &np2, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }

/*     Compute D12HAT using Cholesky decomposition. */

    dpotrf_("L", &m2, &dwork[id12], &m2, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }
/*             _ */
/*     Compute Z = In - Y*X/gamma^2 and its LU factorization in AK . */

    iwrk = iw1;
    dlaset_("Full", n, n, &c_b7, &c_b8, &ak[ak_offset], ldak, (ftnlen)4);
/* Computing 2nd power */
    d__2 = *gamma;
    d__1 = -1. / (d__2 * d__2);
    dgemm_("N", "N", n, n, n, &d__1, &y[y_offset], ldy, &x[x_offset], ldx, &
	    c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
    anorm = dlange_("1", n, n, &ak[ak_offset], ldak, &dwork[iwrk], (ftnlen)1);
    dgetrf_(n, n, &ak[ak_offset], ldak, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }
    dgecon_("1", n, &ak[ak_offset], ldak, &anorm, &rcond, &dwork[iwrk], &
	    iwork[*n + 1], info, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

    if (rcond < eps) {
	*info = 1;
	return 0;
    }

    iwb = iw1;
    iwc = iwb + *n * np2;
    iw1 = iwc + (m2 + np2) * *n;
    iw2 = iw1 + *n * m2;

/*     Compute C2' + F12' in BK . */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = np2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    bk[j + i__ * bk_dim1] = c__[np1 + i__ + j * c_dim1] + f[nd2 + i__ 
		    + j * f_dim1];
/* L30: */
	}
/* L40: */
    }
/*                                                          _ */
/*     Compute the transpose of (C2 + F12)*Z , with Z = inv(Z) . */

    dgetrs_("Transpose", n, &np2, &ak[ak_offset], ldak, &iwork[1], &bk[
	    bk_offset], ldbk, &info2, (ftnlen)9);

/*     Compute the transpose of F2*Z . */

    ma02ad_("Full", &m2, n, &f[m1 + 1 + f_dim1], ldf, &dwork[iw1], n, (ftnlen)
	    4);
    dgetrs_("Transpose", n, &m2, &ak[ak_offset], ldak, &iwork[1], &dwork[iw1],
	     n, &info2, (ftnlen)9);

/*     Compute the transpose of C1HAT = F2*Z - D11HAT*(C2 + F12)*Z . */

    dgemm_("N", "T", n, &m2, &np2, &c_b19, &bk[bk_offset], ldbk, &dwork[id11],
	     &m2, &c_b8, &dwork[iw1], n, (ftnlen)1, (ftnlen)1);

/*     Compute CHAT . */

    i__1 = m2 + np2;
    dgemm_("N", "T", &m2, n, &m2, &c_b8, &tu[tu_offset], ldtu, &dwork[iw1], n,
	     &c_b7, &dwork[iwc], &i__1, (ftnlen)1, (ftnlen)1);
    i__1 = m2 + np2;
    ma02ad_("Full", n, &np2, &bk[bk_offset], ldbk, &dwork[iwc + m2], &i__1, (
	    ftnlen)4);
    i__1 = m2 + np2;
    dtrmm_("L", "U", "N", "N", &np2, n, &c_b19, &dwork[id21], &np2, &dwork[
	    iwc + m2], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute B2 + H12 . */

    ij = iw2;
    i__1 = m2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[ij] = b[i__ + (m1 + j) * b_dim1] + h__[i__ + (nd1 + j) * 
		    h_dim1];
	    ++ij;
/* L50: */
	}
/* L60: */
    }

/*     Compute A + HC in AK . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
    dgemm_("N", "N", n, n, np, &c_b8, &h__[h_offset], ldh, &c__[c_offset], 
	    ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute AHAT = A + HC + (B2 + H12)*C1HAT in AK . */

    dgemm_("N", "T", n, n, &m2, &c_b8, &dwork[iw2], n, &dwork[iw1], n, &c_b8, 
	    &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute B1HAT = -H2 + (B2 + H12)*D11HAT in BK . */

    dlacpy_("Full", n, &np2, &h__[(np1 + 1) * h_dim1 + 1], ldh, &bk[bk_offset]
	    , ldbk, (ftnlen)4);
    dgemm_("N", "N", n, &np2, &m2, &c_b8, &dwork[iw2], n, &dwork[id11], &m2, &
	    c_b19, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute the first block of BHAT, BHAT1 . */

    dgemm_("N", "N", n, &np2, &np2, &c_b8, &bk[bk_offset], ldbk, &ty[
	    ty_offset], ldty, &c_b7, &dwork[iwb], n, (ftnlen)1, (ftnlen)1);

/*     Compute Tu*D11HAT . */

    dgemm_("N", "N", &m2, &np2, &m2, &c_b8, &tu[tu_offset], ldtu, &dwork[id11]
	    , &m2, &c_b7, &dwork[iw1], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute Tu*D11HAT*Ty in DK . */

    dgemm_("N", "N", &m2, &np2, &np2, &c_b8, &dwork[iw1], &m2, &ty[ty_offset],
	     ldty, &c_b7, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute P = Im2 + Tu*D11HAT*Ty*D22 and its condition. */

    iw2 = iw1 + m2 * np2;
    iwrk = iw2 + m2 * m2;
    dlaset_("Full", &m2, &m2, &c_b7, &c_b8, &dwork[iw2], &m2, (ftnlen)4);
    dgemm_("N", "N", &m2, &m2, &np2, &c_b8, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b8, &dwork[iw2], &m2, (ftnlen)1, (
	    ftnlen)1);
    anorm = dlange_("1", &m2, &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1);
    dgetrf_(&m2, &m2, &dwork[iw2], &m2, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 2;
	return 0;
    }
    dgecon_("1", &m2, &dwork[iw2], &m2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    m2 + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

    if (rcond < eps) {
	*info = 2;
	return 0;
    }

/*     Find the controller matrix CK, CK = inv(P)*CHAT(1:M2,:) . */

    i__1 = m2 + np2;
    dlacpy_("Full", &m2, n, &dwork[iwc], &i__1, &ck[ck_offset], ldck, (ftnlen)
	    4);
    dgetrs_("NoTranspose", &m2, n, &dwork[iw2], &m2, &iwork[1], &ck[ck_offset]
	    , ldck, &info2, (ftnlen)11);

/*     Find the controller matrices AK, BK, and DK, exploiting the */
/*     special structure of the relations. */

/*     Compute Q = Inp2 + D22*Tu*D11HAT*Ty and its LU factorization. */

    iw3 = iw2 + np2 * np2;
    iw4 = iw3 + np2 * m2;
    iwrk = iw4 + np2 * np2;
    dlaset_("Full", &np2, &np2, &c_b7, &c_b8, &dwork[iw2], &np2, (ftnlen)4);
    dgemm_("N", "N", &np2, &np2, &m2, &c_b8, &d__[np1 + 1 + (m1 + 1) * d_dim1]
	    , ldd, &dk[dk_offset], lddk, &c_b8, &dwork[iw2], &np2, (ftnlen)1, 
	    (ftnlen)1);
    dgetrf_(&np2, &np2, &dwork[iw2], &np2, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 2;
	return 0;
    }

/*     Compute A1 = inv(Q)*D22 and inv(Q) . */

    dlacpy_("Full", &np2, &m2, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, &dwork[
	    iw3], &np2, (ftnlen)4);
    dgetrs_("NoTranspose", &np2, &m2, &dwork[iw2], &np2, &iwork[1], &dwork[
	    iw3], &np2, &info2, (ftnlen)11);
    i__1 = *ldwork - iwrk + 1;
    dgetri_(&np2, &dwork[iw2], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Compute A2 = ( inv(Ty) - inv(Q)*inv(Ty) - */
/*                    A1*Tu*D11HAT )*inv(D21HAT) . */

    dlacpy_("Full", &np2, &np2, &ty[ty_offset], ldty, &dwork[iw4], &np2, (
	    ftnlen)4);
    dgetrf_(&np2, &np2, &dwork[iw4], &np2, &iwork[1], &info2);
    i__1 = *ldwork - iwrk + 1;
    dgetri_(&np2, &dwork[iw4], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);

    dlacpy_("Full", &np2, &np2, &dwork[iw4], &np2, &dwork[iwrk], &np2, (
	    ftnlen)4);
    dgemm_("N", "N", &np2, &np2, &np2, &c_b19, &dwork[iw2], &np2, &dwork[iwrk]
	    , &np2, &c_b8, &dwork[iw4], &np2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &np2, &np2, &m2, &c_b19, &dwork[iw3], &np2, &dwork[iw1], 
	    &m2, &c_b8, &dwork[iw4], &np2, (ftnlen)1, (ftnlen)1);
    dtrmm_("R", "U", "N", "N", &np2, &np2, &c_b8, &dwork[id21], &np2, &dwork[
	    iw4], &np2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute [ A1  A2 ]*CHAT . */

    i__1 = m2 + np2;
    i__2 = m2 + np2;
    dgemm_("N", "N", &np2, n, &i__1, &c_b8, &dwork[iw3], &np2, &dwork[iwc], &
	    i__2, &c_b7, &dwork[iwrk], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute AK := AHAT - BHAT1*[ A1  A2 ]*CHAT . */

    dgemm_("N", "N", n, n, &np2, &c_b19, &dwork[iwb], n, &dwork[iwrk], &np2, &
	    c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute BK := BHAT1*inv(Q) . */

    dgemm_("N", "N", n, &np2, &np2, &c_b8, &dwork[iwb], n, &dwork[iw2], &np2, 
	    &c_b7, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute DK := Tu*D11HAT*Ty*inv(Q) . */

    dgemm_("N", "N", &m2, &np2, &np2, &c_b8, &dk[dk_offset], lddk, &dwork[iw2]
	    , &np2, &c_b7, &dwork[iw3], &m2, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", &m2, &np2, &dwork[iw3], &m2, &dk[dk_offset], lddk, (
	    ftnlen)4);

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10RD *** */
} /* sb10rd_ */

