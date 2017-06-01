/* SB10SD.f -- translated by f2c (version 20100827).
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

static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb10sd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *x, integer *ldx, doublereal *y, integer *ldy, doublereal *
	rcond, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, 
	    i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m1, m2, nd1, nd2, np1, np2, iw2, iwb, iwc, iwg, iwi, 
	    iwq, iwr, iws, iwt, iwu, iwv;
    static doublereal sepd, ferr, toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb02od_(char *, char *, char *, char *, char *, char *, integer *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), sb02sd_(char *, char *, 
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), 
	    mb01rx_(char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpocon_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
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

/*     To compute the matrices of the H2 optimal controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     for the normalized discrete-time system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 | D11 D12 |   | C | D | */
/*                   | C2 | D21  0  | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank with D12 = | 0 | and D21 is */
/*                                             | I | */
/*          full row rank with D21 = | 0 I | as obtained by the */
/*          SLICOT Library routine SB10PD, */

/*               j*Theta */
/*     (A3) | A-e       *I  B2  | has full column rank for all */
/*          |    C1         D12 | */

/*          0 <= Theta < 2*Pi , */


/*               j*Theta */
/*     (A4) | A-e       *I  B1  | has full row rank for all */
/*          |    C2         D21 | */

/*          0 <= Theta < 2*Pi . */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0, */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0, */
/*             M-NCON >= NMEAS. */

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
/*             system input/output matrix D. Only the leading */
/*             (NP-NP2)-by-(M-M2) submatrix D11 is used. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

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

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the matrix */
/*             X, solution of the X-Riccati equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains the matrix */
/*             Y, solution of the Y-Riccati equation. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*             RCOND contains estimates of the reciprocal condition */
/*             numbers of the matrices which are to be inverted and the */
/*             reciprocal condition numbers of the Riccati equations */
/*             which have to be solved during the computation of the */
/*             controller. (See the description of the algorithm in [2].) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      matrix Im2 + B2'*X2*B2; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      matrix Ip2 + C2*Y2*C2'; */
/*             RCOND(3) contains the reciprocal condition number of the */
/*                      X-Riccati equation; */
/*             RCOND(4) contains the reciprocal condition number of the */
/*                      Y-Riccati equation. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used in determining the nonsingularity of the */
/*             matrices which must be inverted. If TOL <= 0, then a */
/*             default value equal to sqrt(EPS) is used, where EPS is the */
/*             relative machine precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(M2,2*N,N*N,NP2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(1, 14*N*N+6*N+max(14*N+23,16*N), */
/*                              M2*(N+M2+max(3,M1)), NP2*(N+NP2+3)), */
/*             where M1 = M - M2. */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 2:  if the matrix Im2 + B2'*X2*B2 is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL); */
/*             = 3:  if the Y-Riccati equation was not solved */
/*                   successfully; */
/*             = 4:  if the matrix Ip2 + C2*Y2*C2' is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL). */

/*     METHOD */

/*     The routine implements the formulas given in [1]. The X- and */
/*     Y-Riccati equations are solved with condition estimates. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M. */
/*         Fortran 77 routines for Hinf and H2 design of linear */
/*         discrete-time control systems. */
/*         Report 99-8, Department of Engineering, Leicester University, */
/*         April 1999. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     matrices which are to be inverted and on the condition numbers of */
/*     the matrix Riccati equations which are to be solved in the */
/*     computation of the controller. (The corresponding reciprocal */
/*     condition numbers are given in the output array RCOND.) */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     January 2003. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, LQG, LQR, optimal */
/*     regulator, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External functions .. */
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
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --rcond;
    --iwork;
    --dwork;
    --bwork;

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
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*np)) {
	*info = -11;
    } else if (*ldd < max(1,*np)) {
	*info = -13;
    } else if (*ldak < max(1,*n)) {
	*info = -15;
    } else if (*ldbk < max(1,*n)) {
	*info = -17;
    } else if (*ldck < max(1,m2)) {
	*info = -19;
    } else if (*lddk < max(1,m2)) {
	*info = -21;
    } else if (*ldx < max(1,*n)) {
	*info = -23;
    } else if (*ldy < max(1,*n)) {
	*info = -25;
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
	i__3 = *n * 14 + 23, i__4 = *n << 4;
	i__1 = 1, i__2 = *n * 14 * *n + *n * 6 + max(i__3,i__4), i__1 = max(
		i__1,i__2), i__2 = m2 * (*n + m2 + max(3,m1)), i__1 = max(
		i__1,i__2), i__2 = np2 * (*n + np2 + 3);
	minwrk = max(i__1,i__2);
	if (*ldwork < minwrk) {
	    *info = -30;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10SD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
	rcond[1] = 1.;
	rcond[2] = 1.;
	rcond[3] = 1.;
	rcond[4] = 1.;
	dwork[1] = 1.;
	return 0;
    }

    nd1 = np1 - m2;
    nd2 = m1 - np2;
    toll = *tol;
    if (toll <= 0.) {

/*        Set the default value of the tolerance for nonsingularity test. */

	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
    }

/*     Workspace usage. */

    iwq = 1;
    iwg = iwq + *n * *n;
    iwr = iwg + *n * *n;
    iwi = iwr + (*n << 1);
    iwb = iwi + (*n << 1);
    iws = iwb + (*n << 1);
    iwt = iws + (*n << 2) * *n;
    iwu = iwt + (*n << 2) * *n;
    iwrk = iwu + (*n << 2) * *n;
    iwc = iwr;
    iwv = iwc + *n * *n;

/*     Compute Ax = A - B2*D12'*C1 in AK . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
    dgemm_("N", "N", n, n, &m2, &c_b7, &b[(m1 + 1) * b_dim1 + 1], ldb, &c__[
	    nd1 + 1 + c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (
	    ftnlen)1);

/*     Compute Cx = C1'*C1 - C1'*D12*D12'*C1 . */

    if (nd1 > 0) {
	dsyrk_("L", "T", n, &nd1, &c_b8, &c__[c_offset], ldc, &c_b12, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
    } else {
	dlaset_("L", n, n, &c_b12, &c_b12, &dwork[iwq], n, (ftnlen)1);
    }

/*     Compute Dx = B2*B2' . */

    dsyrk_("L", "N", n, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &c_b12, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the discrete-time Riccati equation */
/*        Ax'*inv(In + X2*Dx)*X2*Ax - X2 + Cx  = 0 . */
/*     Workspace:  need   14*N*N + 6*N + max(14*N+23,16*N); */
/*                 prefer larger. */

    i__1 = *n << 1;
    i__2 = *n << 1;
    i__3 = *n << 1;
    i__4 = *ldwork - iwrk + 1;
    sb02od_("D", "G", "N", "L", "Z", "S", n, &m2, &np1, &ak[ak_offset], ldak, 
	    &dwork[iwg], n, &dwork[iwq], n, &dwork[iwrk], m, &dwork[iwrk], n, 
	    &rcond2, &x[x_offset], ldx, &dwork[iwr], &dwork[iwi], &dwork[iwb],
	     &dwork[iws], &i__1, &dwork[iwt], &i__2, &dwork[iwu], &i__3, &
	    toll, &iwork[1], &dwork[iwrk], &i__4, &bwork[1], &info2, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Condition estimation. */
/*     Workspace:  need   4*N*N + max(N*N+5*N,max(3,2*N*N)+N*N); */
/*                 prefer larger. */

    iwrk = iwv + *n * *n;
    i__1 = *ldwork - iwrk + 1;
    sb02sd_("C", "N", "N", "L", "O", n, &ak[ak_offset], ldak, &dwork[iwc], n, 
	    &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &x[x_offset], ldx,
	     &sepd, &rcond[3], &ferr, &iwork[1], &dwork[iwrk], &i__1, &info2, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	rcond[3] = 0.;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

    iw2 = m2 * *n + 1;
    iwrk = iw2 + m2 * m2;

/*     Compute B2'*X2 . */

    dgemm_("T", "N", &m2, n, n, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &x[
	    x_offset], ldx, &c_b12, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute Im2 + B2'*X2*B2 . */

    dlaset_("L", &m2, &m2, &c_b12, &c_b8, &dwork[iw2], &m2, (ftnlen)1);
    mb01rx_("Left", "Lower", "N", &m2, n, &c_b8, &c_b8, &dwork[iw2], &m2, &
	    dwork[1], &m2, &b[(m1 + 1) * b_dim1 + 1], ldb, &info2, (ftnlen)4, 
	    (ftnlen)5, (ftnlen)1);

/*     Compute the Cholesky factorization of Im2 + B2'*X2*B2 . */
/*     Workspace:  need   M2*N + M2*M2 + max(3*M2,M2*M1); */
/*                 prefer larger. */

    anorm = dlansy_("I", "L", &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1, 
	    (ftnlen)1);
    dpotrf_("L", &m2, &dwork[iw2], &m2, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 2;
	return 0;
    }
    dpocon_("L", &m2, &dwork[iw2], &m2, &anorm, &rcond[1], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

    if (rcond[1] < toll) {
	*info = 2;
	return 0;
    }

/*     Compute -( B2'*X2*A + D12'*C1 ) in CK . */

    dlacpy_("Full", &m2, n, &c__[nd1 + 1 + c_dim1], ldc, &ck[ck_offset], ldck,
	     (ftnlen)4);
    dgemm_("N", "N", &m2, n, n, &c_b7, &dwork[1], &m2, &a[a_offset], lda, &
	    c_b7, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);

/*     Compute F2 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*A + D12'*C1 ) . */

    dpotrs_("L", &m2, n, &dwork[iw2], &m2, &ck[ck_offset], ldck, &info2, (
	    ftnlen)1);

/*     Compute -( B2'*X2*B1 + D12'*D11 ) . */

    dlacpy_("Full", &m2, &m1, &d__[nd1 + 1 + d_dim1], ldd, &dwork[iwrk], &m2, 
	    (ftnlen)4);
    dgemm_("N", "N", &m2, &m1, n, &c_b7, &dwork[1], &m2, &b[b_offset], ldb, &
	    c_b7, &dwork[iwrk], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute F0 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*B1 + D12'*D11 ) . */

    dpotrs_("L", &m2, &m1, &dwork[iw2], &m2, &dwork[iwrk], &m2, &info2, (
	    ftnlen)1);

/*     Save F0*D21' in DK . */

    dlacpy_("Full", &m2, &np2, &dwork[iwrk + nd2 * m2], &m2, &dk[dk_offset], 
	    lddk, (ftnlen)4);

/*     Workspace usage. */

    iwrk = iwu + (*n << 2) * *n;

/*     Compute Ay = A - B1*D21'*C2 in AK . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
    dgemm_("N", "N", n, n, &np2, &c_b7, &b[(nd2 + 1) * b_dim1 + 1], ldb, &c__[
	    np1 + 1 + c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (
	    ftnlen)1);

/*     Transpose Ay in-situ. */

    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	dswap_(&j, &ak[j + 1 + ak_dim1], ldak, &ak[(j + 1) * ak_dim1 + 1], &
		c__1);
/* L20: */
    }

/*     Compute Cy = B1*B1' - B1*D21'*D21*B1' . */

    if (nd2 > 0) {
	dsyrk_("U", "N", n, &nd2, &c_b8, &b[b_offset], ldb, &c_b12, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
    } else {
	dlaset_("U", n, n, &c_b12, &c_b12, &dwork[iwq], n, (ftnlen)1);
    }

/*     Compute Dy = C2'*C2 . */

    dsyrk_("U", "T", n, &np2, &c_b8, &c__[np1 + 1 + c_dim1], ldc, &c_b12, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the discrete-time Riccati equation */
/*        Ay*inv( In + Y2*Dy )*Y2*Ay' - Y2 + Cy = 0 . */

    i__1 = *n << 1;
    i__2 = *n << 1;
    i__3 = *n << 1;
    i__4 = *ldwork - iwrk + 1;
    sb02od_("D", "G", "N", "U", "Z", "S", n, &np2, &m1, &ak[ak_offset], ldak, 
	    &dwork[iwg], n, &dwork[iwq], n, &dwork[iwrk], m, &dwork[iwrk], n, 
	    &rcond2, &y[y_offset], ldy, &dwork[iwr], &dwork[iwi], &dwork[iwb],
	     &dwork[iws], &i__1, &dwork[iwt], &i__2, &dwork[iwu], &i__3, &
	    toll, &iwork[1], &dwork[iwrk], &i__4, &bwork[1], &info2, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 3;
	return 0;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Condition estimation. */

    iwrk = iwv + *n * *n;
    i__1 = *ldwork - iwrk + 1;
    sb02sd_("C", "N", "N", "U", "O", n, &ak[ak_offset], ldak, &dwork[iwc], n, 
	    &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &y[y_offset], ldy,
	     &sepd, &rcond[4], &ferr, &iwork[1], &dwork[iwrk], &i__1, &info2, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	rcond[4] = 0.;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

    iw2 = *n * np2 + 1;
    iwrk = iw2 + np2 * np2;

/*     Compute Y2*C2' . */

    dgemm_("N", "T", n, &np2, n, &c_b8, &y[y_offset], ldy, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b12, &dwork[1], n, (ftnlen)1, (ftnlen)1);

/*     Compute Ip2 + C2*Y2*C2' . */

    dlaset_("U", &np2, &np2, &c_b12, &c_b8, &dwork[iw2], &np2, (ftnlen)1);
    mb01rx_("Left", "Upper", "N", &np2, n, &c_b8, &c_b8, &dwork[iw2], &np2, &
	    c__[np1 + 1 + c_dim1], ldc, &dwork[1], n, &info2, (ftnlen)4, (
	    ftnlen)5, (ftnlen)1);

/*     Compute the Cholesky factorization of Ip2 + C2*Y2*C2' . */

    anorm = dlansy_("I", "U", &np2, &dwork[iw2], &np2, &dwork[iwrk], (ftnlen)
	    1, (ftnlen)1);
    dpotrf_("U", &np2, &dwork[iw2], &np2, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 4;
	return 0;
    }
    dpocon_("U", &np2, &dwork[iw2], &np2, &anorm, &rcond[2], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

    if (rcond[2] < toll) {
	*info = 4;
	return 0;
    }

/*     Compute A*Y2*C2' + B1*D21' in BK . */

    dlacpy_("Full", n, &np2, &b[(nd2 + 1) * b_dim1 + 1], ldb, &bk[bk_offset], 
	    ldbk, (ftnlen)4);
    dgemm_("N", "N", n, &np2, n, &c_b8, &a[a_offset], lda, &dwork[1], n, &
	    c_b8, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute L2 = -( A*Y2*C2' + B1*D21' )*inv( Ip2 + C2*Y2*C2' ) . */

    dtrsm_("R", "U", "N", "N", n, &np2, &c_b7, &dwork[iw2], &np2, &bk[
	    bk_offset], ldbk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    dtrsm_("R", "U", "T", "N", n, &np2, &c_b8, &dwork[iw2], &np2, &bk[
	    bk_offset], ldbk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute F2*Y2*C2' + F0*D21' . */

    dgemm_("N", "N", &m2, &np2, n, &c_b8, &ck[ck_offset], ldck, &dwork[1], n, 
	    &c_b8, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute DK = L0 = ( F2*Y2*C2' + F0*D21' )*inv( Ip2 + C2*Y2*C2' ) . */

    dtrsm_("R", "U", "N", "N", &m2, &np2, &c_b8, &dwork[iw2], &np2, &dk[
	    dk_offset], lddk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    dtrsm_("R", "U", "T", "N", &m2, &np2, &c_b8, &dwork[iw2], &np2, &dk[
	    dk_offset], lddk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute CK = F2 - L0*C2 . */

    dgemm_("N", "N", &m2, n, &np2, &c_b7, &dk[dk_offset], lddk, &c__[np1 + 1 
	    + c_dim1], ldc, &c_b8, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1)
	    ;

/*     Find AK = A + B2*( F2 - L0*C2 ) + L2*C2 . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
    dgemm_("N", "N", n, n, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &ck[
	    ck_offset], ldck, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)
	    1);
    dgemm_("N", "N", n, n, &np2, &c_b8, &bk[bk_offset], ldbk, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Find BK = -L2 + B2*L0 . */

    dgemm_("N", "N", n, &np2, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &dk[
	    dk_offset], lddk, &c_b7, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)
	    1);

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10SD *** */
} /* sb10sd_ */

