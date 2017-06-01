/* SB10ED.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb10ed_(integer *n, integer *m, integer *np, integer *
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
	     dk_dim1, dk_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, m1, m2, nl, m2l, np1, np2, lw1, lw2, lw3, lw4, lw5, 
	    lw6, iwc, iwd, nlp, npl, iwx, iwy;
    static doublereal toll;
    static integer iwrk, iwtu, iwty, info2;
    extern /* Subroutine */ int sb10pd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), sb10sd_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *), sb10td_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
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

/*                           | AK | BK | */
/*                       K = |----|----| */
/*                           | CK | DK | */

/*     for the discrete-time system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| , */
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank and D21 is full row rank, */

/*               j*Theta */
/*     (A3) | A-e       *I  B2  | has full column rank for all */
/*          |    C1         D12 | */

/*          0 <= Theta < 2*Pi , */


/*               j*Theta */
/*     (A4) | A-e       *I  B1  |  has full row rank for all */
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

/*     A       (input/worksp.) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */
/*             This array is modified internally, but it is restored on */
/*             exit. */

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

/*     RCOND   (output) DOUBLE PRECISION array, dimension (7) */
/*             RCOND contains estimates the reciprocal condition */
/*             numbers of the matrices which are to be inverted and the */
/*             reciprocal condition numbers of the Riccati equations */
/*             which have to be solved during the computation of the */
/*             controller. (See the description of the algorithm in [2].) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix TU; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix TY; */
/*             RCOND(3) contains the reciprocal condition number of the */
/*                      matrix Im2 + B2'*X2*B2; */
/*             RCOND(4) contains the reciprocal condition number of the */
/*                      matrix Ip2 + C2*Y2*C2'; */
/*             RCOND(5) contains the reciprocal condition number of the */
/*                      X-Riccati equation; */
/*             RCOND(6) contains the reciprocal condition number of the */
/*                      Y-Riccati equation; */
/*             RCOND(7) contains the reciprocal condition number of the */
/*                      matrix Im2 + DKHAT*D22 . */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of the */
/*             transformations applied for diagonalizing D12 and D21, */
/*             and for checking the nonsingularity of the matrices to be */
/*             inverted. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*M2,2*N,N*N,NP2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 + */
/*                       max(1,LW1,LW2,LW3,LW4,LW5,LW6), where */
/*             LW1 = (N+NP1+1)*(N+M2) + max(3*(N+M2)+N+NP1,5*(N+M2)), */
/*             LW2 = (N+NP2)*(N+M1+1) + max(3*(N+NP2)+N+M1,5*(N+NP2)), */
/*             LW3 = M2 + NP1*NP1 + max(NP1*max(N,M1),3*M2+NP1,5*M2), */
/*             LW4 = NP2 + M1*M1 + max(max(N,NP1)*M1,3*NP2+M1,5*NP2), */
/*             LW5 = 2*N*N+max(1,14*N*N+6*N+max(14*N+23,16*N),M2*(N+M2+ */
/*                             max(3,M1)),NP2*(N+NP2+3)), */
/*             LW6 = max(N*M2,N*NP2,M2*NP2,M2*M2+4*M2), */
/*             with M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             2*Q*(3*Q+2*N)+max(1,(N+Q)*(N+Q+6),Q*(Q+max(N,Q,5)+1), */
/*                     2*N*N+max(1,14*N*N+6*N+max(14*N+23,16*N), */
/*                               Q*(N+Q+max(Q,3)))). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*                                      j*Theta */
/*             = 1:  if the matrix | A-e       *I  B2  | had not full */
/*                                 |      C1       D12 | */
/*                   column rank in respect to the tolerance EPS; */
/*                                      j*Theta */
/*             = 2:  if the matrix | A-e       *I  B1  |  had not full */
/*                                 |      C2       D21 | */
/*                   row rank in respect to the tolerance EPS; */
/*             = 3:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance TOL; */
/*             = 4:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance TOL; */
/*             = 5:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices |A-I  B2 |, |A-I  B1 |, D12 or D21). */
/*                                |C1   D12|  |C2   D21| */
/*             = 6:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 7:  if the matrix Im2 + B2'*X2*B2 is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL); */
/*             = 8:  if the Y-Riccati equation was not solved */
/*                   successfully; */
/*             = 9:  if the matrix Ip2 + C2*Y2*C2' is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL); */
/*             =10:  if the matrix Im2 + DKHAT*D22 is singular, or its */
/*                   estimated condition number is larger than or equal */
/*                   to 1/TOL. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

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

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Feb. 2000, Nov. 2005. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, optimal regulator, */
/*     robust control. */

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
    nl = max(1,*n);
    npl = max(1,*np);
    m2l = max(1,m2);
    nlp = max(1,np2);

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
    } else if (*lda < nl) {
	*info = -7;
    } else if (*ldb < nl) {
	*info = -9;
    } else if (*ldc < npl) {
	*info = -11;
    } else if (*ldd < npl) {
	*info = -13;
    } else if (*ldak < nl) {
	*info = -15;
    } else if (*ldbk < nl) {
	*info = -17;
    } else if (*ldck < m2l) {
	*info = -19;
    } else if (*lddk < m2l) {
	*info = -21;
    } else {

/*        Compute workspace. */

/* Computing MAX */
	i__1 = (*n + m2) * 3 + *n + np1, i__2 = (*n + m2) * 5;
	lw1 = (*n + np1 + 1) * (*n + m2) + max(i__1,i__2);
/* Computing MAX */
	i__1 = (*n + np2) * 3 + *n + m1, i__2 = (*n + np2) * 5;
	lw2 = (*n + np2) * (*n + m1 + 1) + max(i__1,i__2);
/* Computing MAX */
	i__1 = np1 * max(*n,m1), i__2 = m2 * 3 + np1, i__1 = max(i__1,i__2), 
		i__2 = m2 * 5;
	lw3 = m2 + np1 * np1 + max(i__1,i__2);
/* Computing MAX */
	i__1 = max(*n,np1) * m1, i__2 = np2 * 3 + m1, i__1 = max(i__1,i__2), 
		i__2 = np2 * 5;
	lw4 = np2 + m1 * m1 + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
	i__3 = *n * 14 + 23, i__4 = *n << 4;
	i__1 = 1, i__2 = *n * 14 * *n + *n * 6 + max(i__3,i__4), i__1 = max(
		i__1,i__2), i__2 = m2 * (*n + m2 + max(3,m1)), i__1 = max(
		i__1,i__2), i__2 = np2 * (*n + np2 + 3);
	lw5 = (*n << 1) * *n + max(i__1,i__2);
/* Computing MAX */
	i__1 = *n * m2, i__2 = *n * np2, i__1 = max(i__1,i__2), i__2 = m2 * 
		np2, i__1 = max(i__1,i__2), i__2 = m2 * m2 + (m2 << 2);
	lw6 = max(i__1,i__2);
/* Computing MAX */
	i__1 = max(1,lw1), i__1 = max(i__1,lw2), i__1 = max(i__1,lw3), i__1 = 
		max(i__1,lw4), i__1 = max(i__1,lw5);
	minwrk = *n * *m + *np * (*n + *m) + m2 * m2 + np2 * np2 + max(i__1,
		lw6);
	if (*ldwork < minwrk) {
	    *info = -26;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10ED", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 && max(m2,np2) == 0) {
	rcond[1] = 1.;
	rcond[2] = 1.;
	rcond[3] = 1.;
	rcond[4] = 1.;
	rcond[5] = 1.;
	rcond[6] = 1.;
	rcond[7] = 1.;
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

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], &nl, (ftnlen)4);
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwc], &npl, (ftnlen)4);
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iwd], &npl, (ftnlen)4);

/*     Transform the system so that D12 and D21 satisfy the formulas */
/*     in the computation of the H2 optimal controller. */
/*     Since SLICOT Library routine SB10PD performs the tests */
/*     corresponding to the continuous-time counterparts of the */
/*     assumptions (A3) and (A4), for the frequency w = 0, the */
/*     next SB10PD routine call uses A - I. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + i__ * a_dim1] += -1.;
/* L10: */
    }

    i__1 = *ldwork - iwrk + 1;
    sb10pd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], &nl, &dwork[
	    iwc], &npl, &dwork[iwd], &npl, &dwork[iwtu], &m2l, &dwork[iwty], &
	    nlp, &rcond[1], &toll, &dwork[iwrk], &i__1, &info2);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + i__ * a_dim1] += 1.;
/* L20: */
    }

    if (info2 > 0) {
	*info = info2;
	return 0;
    }
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

    iwx = iwrk;
    iwy = iwx + *n * *n;
    iwrk = iwy + *n * *n;

/*     Compute the optimal H2 controller for the normalized system. */

    i__1 = *ldwork - iwrk + 1;
    sb10sd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], &nl, &dwork[
	    iwc], &npl, &dwork[iwd], &npl, &ak[ak_offset], ldak, &bk[
	    bk_offset], ldbk, &ck[ck_offset], ldck, &dk[dk_offset], lddk, &
	    dwork[iwx], &nl, &dwork[iwy], &nl, &rcond[3], &toll, &iwork[1], &
	    dwork[iwrk], &i__1, &bwork[1], &info2);
    if (info2 > 0) {
	*info = info2 + 5;
	return 0;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

    iwrk = iwx;

/*     Compute the H2 optimal controller for the original system. */

    i__1 = *ldwork - iwrk + 1;
    sb10td_(n, m, np, ncon, nmeas, &dwork[iwd], &npl, &dwork[iwtu], &m2l, &
	    dwork[iwty], &nlp, &ak[ak_offset], ldak, &bk[bk_offset], ldbk, &
	    ck[ck_offset], ldck, &dk[dk_offset], lddk, &rcond[7], &toll, &
	    iwork[1], &dwork[iwrk], &i__1, &info2);
    if (info2 > 0) {
	*info = 10;
	return 0;
    }

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10ED *** */
} /* sb10ed_ */

