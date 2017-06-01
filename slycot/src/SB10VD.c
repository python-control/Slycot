/* SB10VD.f -- translated by f2c (version 20100827).
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

static doublereal c_b6 = -1.;
static doublereal c_b7 = 1.;
static doublereal c_b11 = 0.;

/* Subroutine */ int sb10vd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *f, integer *
	ldf, doublereal *h__, integer *ldh, doublereal *x, integer *ldx, 
	doublereal *y, integer *ldy, doublereal *xycond, integer *iwork, 
	doublereal *dwork, integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, f_dim1, 
	    f_offset, h_dim1, h_offset, x_dim1, x_offset, y_dim1, y_offset, 
	    i__1;

    /* Local variables */
    static integer m1, m2, n2, nd1, nd2, np1, np2, iwg, iwi;
    static doublereal sep;
    static integer iwq, iwr, iws, iwt, iwv;
    static doublereal ferr;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb02rd_(char *, char *, char *, char *, char *, char *, char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    logical *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(char *, char *, integer *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), dlacpy_(char *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static integer lwamax, minwrk;


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

/*     To compute the state feedback and the output injection */
/*     matrices for an H2 optimal n-state controller for the system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank with D12 = | 0 | and D21 is */
/*                                             | I | */
/*          full row rank with D21 = | 0 I | as obtained by the */
/*          SLICOT Library routine SB10UD. Matrix D is not used */
/*          explicitly. */

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

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading NCON-by-N part of this array contains the */
/*             state feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= max(1,NCON). */

/*     H       (output) DOUBLE PRECISION array, dimension (LDH,NMEAS) */
/*             The leading N-by-NMEAS part of this array contains the */
/*             output injection matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

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

/*     XYCOND  (output) DOUBLE PRECISION array, dimension (2) */
/*             XYCOND(1) contains an estimate of the reciprocal condition */
/*                       number of the X-Riccati equation; */
/*             XYCOND(2) contains an estimate of the reciprocal condition */
/*                       number of the Y-Riccati equation. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*N,N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 13*N*N + 12*N + 5. */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 2:  if the Y-Riccati equation was not solved */
/*                   successfully. */

/*     METHOD */

/*     The routine implements the formulas given in [1], [2]. The X- */
/*     and Y-Riccati equations are solved with condition and accuracy */
/*     estimates [3]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     [3] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V. */
/*         DGRSVX and DMSRIC: Fortan 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Tech. Univ. */
/*         Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */

/*     The precision of the solution of the matrix Riccati equations */
/*     can be controlled by the values of the condition numbers */
/*     XYCOND(1) and XYCOND(2) of these equations. */

/*     FURTHER COMMENTS */

/*     The Riccati equations are solved by the Schur approach */
/*     implementing condition and accuracy estimates. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999. */

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
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --xycond;
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
    } else if (*ldf < max(1,*ncon)) {
	*info = -13;
    } else if (*ldh < max(1,*n)) {
	*info = -15;
    } else if (*ldx < max(1,*n)) {
	*info = -17;
    } else if (*ldy < max(1,*n)) {
	*info = -19;
    } else {

/*        Compute workspace. */

	minwrk = *n * 13 * *n + *n * 12 + 5;
	if (*ldwork < minwrk) {
	    *info = -23;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10VD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
	dwork[1] = 1.;
	xycond[1] = 1.;
	xycond[2] = 1.;
	return 0;
    }

    nd1 = np1 - m2;
    nd2 = m1 - np2;
    n2 = *n << 1;

/*     Workspace usage. */

    iwq = *n * *n + 1;
    iwg = iwq + *n * *n;
    iwt = iwg + *n * *n;
    iwv = iwt + *n * *n;
    iwr = iwv + *n * *n;
    iwi = iwr + n2;
    iws = iwi + n2;
    iwrk = iws + (*n << 2) * *n;

/*     Compute Ax = A - B2*D12'*C1 . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
    dgemm_("N", "N", n, n, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &c__[
	    nd1 + 1 + c_dim1], ldc, &c_b7, &dwork[1], n, (ftnlen)1, (ftnlen)1)
	    ;

/*     Compute Cx = C1'*C1 - C1'*D12*D12'*C1 . */

    if (nd1 > 0) {
	dsyrk_("L", "T", n, &nd1, &c_b7, &c__[c_offset], ldc, &c_b11, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
    } else {
	dlaset_("L", n, n, &c_b11, &c_b11, &dwork[iwq], n, (ftnlen)1);
    }

/*     Compute Dx = B2*B2' . */

    dsyrk_("L", "N", n, &m2, &c_b7, &b[(m1 + 1) * b_dim1 + 1], ldb, &c_b11, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the Riccati equation Ax'*X + X*Ax + Cx - X*Dx*X = 0 . */
/*     Workspace:  need   13*N*N + 12*N + 5; */
/*                 prefer larger. */

    i__1 = *ldwork - iwrk + 1;
    sb02rd_("All", "Continuous", "NotUsed", "NoTranspose", "Lower", "General"
	    "Scaling", "Stable", "NotFactored", "Original", n, &dwork[1], n, &
	    dwork[iwt], n, &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &x[
	    x_offset], ldx, &sep, &xycond[1], &ferr, &dwork[iwr], &dwork[iwi],
	     &dwork[iws], &n2, &iwork[1], &dwork[iwrk], &i__1, &bwork[1], &
	    info2, (ftnlen)3, (ftnlen)10, (ftnlen)7, (ftnlen)11, (ftnlen)5, (
	    ftnlen)14, (ftnlen)6, (ftnlen)11, (ftnlen)8);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }

    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Compute F = -D12'*C1 - B2'*X . */

    dlacpy_("Full", &m2, n, &c__[nd1 + 1 + c_dim1], ldc, &f[f_offset], ldf, (
	    ftnlen)4);
    dgemm_("T", "N", &m2, n, n, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &x[
	    x_offset], ldx, &c_b6, &f[f_offset], ldf, (ftnlen)1, (ftnlen)1);

/*     Compute Ay = A - B1*D21'*C2 . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
    dgemm_("N", "N", n, n, &np2, &c_b6, &b[(nd2 + 1) * b_dim1 + 1], ldb, &c__[
	    np1 + 1 + c_dim1], ldc, &c_b7, &dwork[1], n, (ftnlen)1, (ftnlen)1)
	    ;

/*     Compute Cy = B1*B1' - B1*D21'*D21*B1' . */

    if (nd2 > 0) {
	dsyrk_("U", "N", n, &nd2, &c_b7, &b[b_offset], ldb, &c_b11, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
    } else {
	dlaset_("U", n, n, &c_b11, &c_b11, &dwork[iwq], n, (ftnlen)1);
    }

/*     Compute Dy = C2'*C2 . */

    dsyrk_("U", "T", n, &np2, &c_b7, &c__[np1 + 1 + c_dim1], ldc, &c_b11, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the Riccati equation Ay*Y + Y*Ay' + Cy - Y*Dy*Y = 0 . */
/*     Workspace:  need   13*N*N + 12*N + 5; */
/*                 prefer larger. */

    i__1 = *ldwork - iwrk + 1;
    sb02rd_("All", "Continuous", "NotUsed", "Transpose", "Upper", "GeneralSc"
	    "aling", "Stable", "NotFactored", "Original", n, &dwork[1], n, &
	    dwork[iwt], n, &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &y[
	    y_offset], ldy, &sep, &xycond[2], &ferr, &dwork[iwr], &dwork[iwi],
	     &dwork[iws], &n2, &iwork[1], &dwork[iwrk], &i__1, &bwork[1], &
	    info2, (ftnlen)3, (ftnlen)10, (ftnlen)7, (ftnlen)9, (ftnlen)5, (
	    ftnlen)14, (ftnlen)6, (ftnlen)11, (ftnlen)8);
    if (info2 > 0) {
	*info = 2;
	return 0;
    }

/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Compute H = -B1*D21' - Y*C2' . */

    dlacpy_("Full", n, &np2, &b[(nd2 + 1) * b_dim1 + 1], ldb, &h__[h_offset], 
	    ldh, (ftnlen)4);
    dgemm_("N", "T", n, &np2, n, &c_b6, &y[y_offset], ldy, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b6, &h__[h_offset], ldh, (ftnlen)1, (ftnlen)1);

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10VD *** */
} /* sb10vd_ */

