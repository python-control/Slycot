/* SB10HD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb10hd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *rcond, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, 
	    i__8;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m1, m2, np1, np2, iwc, iwd, iwf, iwh, iwy;
    static doublereal toll;
    static integer iwrk, iwtu, iwty, info2;
    extern /* Subroutine */ int sb10ud_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *), sb10vd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *), sb10wd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

/*     To compute the matrices of the H2 optimal n-state controller */

/*              | AK | BK | */
/*          K = |----|----| */
/*              | CK | DK | */

/*     for the system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| , */
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) The block D11 of D is zero, */

/*     (A3) D12 is full column rank and D21 is full row rank. */

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
/*             system input/output matrix D. */

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

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix; */
/*             RCOND(3) contains an estimate of the reciprocal condition */
/*                      number of the X-Riccati equation; */
/*             RCOND(4) contains an estimate of the reciprocal condition */
/*                      number of the Y-Riccati equation. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of the applied */
/*             transformations for computing the normalized form in */
/*             SLICOT Library routine SB10UD. Transformation matrices */
/*             whose reciprocal condition numbers are less than TOL are */
/*             not allowed. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*N,N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 + */
/*                       max(max(M2 + NP1*NP1 + */
/*                               max(NP1*N,3*M2+NP1,5*M2), */
/*                               NP2 + M1*M1 + */
/*                               max(M1*N,3*NP2+M1,5*NP2), */
/*                               N*M2,NP2*N,NP2*M2,1), */
/*                               N*(14*N+12+M2+NP2)+5), */
/*             where M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             2*Q*(3*Q+2*N)+max(1,Q*(Q+max(N,5)+1),N*(14*N+12+2*Q)+5). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance TOL; */
/*             = 2:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance TOL; */
/*             = 3:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices D12 or D21). */
/*             = 4:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 5:  if the Y-Riccati equation was not solved */
/*                   successfully. */

/*     METHOD */

/*     The routine implements the formulas given in [1], [2]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations and on the condition numbers of */
/*     the two Riccati equations, as given by the values of RCOND(1), */
/*     RCOND(2), RCOND(3) and RCOND(4), respectively. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, Oct. 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Jan. 2000, Feb. 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, optimal regulator, */
/*     robust control. */

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
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	i__5 = np1 * *n, i__6 = m2 * 3 + np1, i__5 = max(i__5,i__6), i__6 = 
		m2 * 5;
/* Computing MAX */
	i__7 = m1 * *n, i__8 = np2 * 3 + m1, i__7 = max(i__7,i__8), i__8 = 
		np2 * 5;
	i__3 = m2 + np1 * np1 + max(i__5,i__6), i__4 = np2 + m1 * m1 + max(
		i__7,i__8), i__3 = max(i__3,i__4), i__4 = *n * m2, i__3 = max(
		i__3,i__4), i__4 = np2 * *n, i__3 = max(i__3,i__4), i__4 = 
		np2 * m2, i__3 = max(i__3,i__4);
	i__1 = max(i__3,1), i__2 = *n * (*n * 14 + 12 + m2 + np2) + 5;
	minwrk = *n * *m + *np * (*n + *m) + m2 * m2 + np2 * np2 + max(i__1,
		i__2);
	if (*ldwork < minwrk) {
	    *info = -26;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10HD", &i__1, (ftnlen)6);
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

    toll = *tol;
    if (toll <= 0.) {

/*        Set the default value of the tolerance for rank tests. */

	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
    }

/*     Workspace usage. */

    iwc = *n * *m + 1;
    iwd = iwc + *np * *n;
    iwtu = iwd + *np * *m;
    iwty = iwtu + m2 * m2;
    iwrk = iwty + np2 * np2;

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwc], np, (ftnlen)4);
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iwd], np, (ftnlen)4);

/*     Transform the system so that D12 and D21 satisfy the formulas */
/*     in the computation of the H2 optimal controller. */

    i__1 = *ldwork - iwrk + 1;
    sb10ud_(n, m, np, ncon, nmeas, &dwork[1], n, &dwork[iwc], np, &dwork[iwd],
	     np, &dwork[iwtu], &m2, &dwork[iwty], &np2, &rcond[1], &toll, &
	    dwork[iwrk], &i__1, &info2);
    if (info2 > 0) {
	*info = info2;
	return 0;
    }
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

    iwy = iwrk;
    iwf = iwy + *n * *n;
    iwh = iwf + m2 * *n;
    iwrk = iwh + *n * np2;

/*     Compute the optimal state feedback and output injection matrices. */
/*     AK is used to store X. */

    i__1 = *ldwork - iwrk + 1;
    sb10vd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwf], &m2, &dwork[iwh], n, &ak[ak_offset], ldak, 
	    &dwork[iwy], n, &rcond[3], &iwork[1], &dwork[iwrk], &i__1, &bwork[
	    1], &info2);
    if (info2 > 0) {
	*info = info2 + 3;
	return 0;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Compute the H2 optimal controller. */

    sb10wd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwd], np, &dwork[iwf], &m2, &dwork[iwh], n, &
	    dwork[iwtu], &m2, &dwork[iwty], &np2, &ak[ak_offset], ldak, &bk[
	    bk_offset], ldbk, &ck[ck_offset], ldck, &dk[dk_offset], lddk, &
	    info2);

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10HD *** */
} /* sb10hd_ */

