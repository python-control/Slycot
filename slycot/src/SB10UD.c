/* SB10UD.f -- translated by f2c (version 20100827).
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
static doublereal c_b15 = 1.;
static doublereal c_b16 = 0.;

/* Subroutine */ int sb10ud_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *tu, integer *
	ldtu, doublereal *ty, integer *ldty, doublereal *rcond, doublereal *
	tol, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, d_dim1, d_offset, tu_dim1, 
	    tu_offset, ty_dim1, ty_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m1, m2, iq, nd1, nd2, np1, np2;
    static doublereal toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dswap_(
	    integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
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

/*     To reduce the matrices D12 and D21 of the linear time-invariant */
/*     system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     to unit diagonal form, and to transform the matrices B and C to */
/*     satisfy the formulas in the computation of the H2 optimal */
/*     controller. */

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

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading NP-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading NP-by-N part of this array contains */
/*             the transformed system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading NP-by-M part of this array must */
/*             contain the system input/output matrix D. */
/*             The (NP-NMEAS)-by-(M-NCON) leading submatrix D11 is not */
/*             referenced. */
/*             On exit, the trailing NMEAS-by-NCON part (in the leading */
/*             NP-by-M part) of this array contains the transformed */
/*             submatrix D22. */
/*             The transformed submatrices D12 = [ 0  Im2 ]' and */
/*             D21 = [ 0  Inp2 ] are not stored. The corresponding part */
/*             of this array contains no useful information. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     TU      (output) DOUBLE PRECISION array, dimension (LDTU,M2) */
/*             The leading M2-by-M2 part of this array contains the */
/*             control transformation matrix TU. */

/*     LDTU    INTEGER */
/*             The leading dimension of the array TU.  LDTU >= max(1,M2). */

/*     TY      (output) DOUBLE PRECISION array, dimension (LDTY,NP2) */
/*             The leading NP2-by-NP2 part of this array contains the */
/*             measurement transformation matrix TY. */

/*     LDTY    INTEGER */
/*             The leading dimension of the array TY. */
/*             LDTY >= max(1,NP2). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (2) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix TU; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix TY. */
/*             RCOND is set even if INFO = 1 or INFO = 2; if INFO = 1, */
/*             then RCOND(2) was not computed, but it is set to 0. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of the applied */
/*             transformations. Transformation matrices TU and TY whose */
/*             reciprocal condition numbers are less than TOL are not */
/*             allowed. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= MAX( M2 + NP1*NP1 + MAX(NP1*N,3*M2+NP1,5*M2), */
/*                            NP2 + M1*M1  + MAX(M1*N,3*NP2+M1,5*NP2), */
/*                            N*M2, NP2*N, NP2*M2, 1 ) */
/*             where M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is */
/*             MAX(1,Q*(Q+MAX(N,5)+1)). */

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
/*                   did not converge (when computing the SVD of D12 or */
/*                   D21). */

/*     METHOD */

/*     The routine performs the transformations described in [1], [2]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The precision of the transformations can be controlled by the */
/*     condition numbers of the matrices TU and TY as given by the */
/*     values of RCOND(1) and RCOND(2), respectively. An error return */
/*     with INFO = 1 or INFO = 2 will be obtained if the condition */
/*     number of TU or TY, respectively, would exceed 1/TOL. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Feb. 2000. */

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
/*     .. External Functions */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

    /* Parameter adjustments */
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    tu_dim1 = *ldtu;
    tu_offset = 1 + tu_dim1;
    tu -= tu_offset;
    ty_dim1 = *ldty;
    ty_offset = 1 + ty_dim1;
    ty -= ty_offset;
    --rcond;
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
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldc < max(1,*np)) {
	*info = -9;
    } else if (*ldd < max(1,*np)) {
	*info = -11;
    } else if (*ldtu < max(1,m2)) {
	*info = -13;
    } else if (*ldty < max(1,np2)) {
	*info = -15;
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
	i__3 = np1 * *n, i__4 = m2 * 3 + np1, i__3 = max(i__3,i__4), i__4 = 
		m2 * 5;
/* Computing MAX */
	i__5 = m1 * *n, i__6 = np2 * 3 + m1, i__5 = max(i__5,i__6), i__6 = 
		np2 * 5;
	i__1 = 1, i__2 = m2 + np1 * np1 + max(i__3,i__4), i__1 = max(i__1,
		i__2), i__2 = np2 + m1 * m1 + max(i__5,i__6), i__1 = max(i__1,
		i__2), i__2 = *n * m2, i__1 = max(i__1,i__2), i__2 = np2 * *n,
		 i__1 = max(i__1,i__2), i__2 = np2 * m2;
	minwrk = max(i__1,i__2);
	if (*ldwork < minwrk) {
	    *info = -19;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10UD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
	rcond[1] = 1.;
	rcond[2] = 1.;
	dwork[1] = 1.;
	return 0;
    }

    nd1 = np1 - m2;
    nd2 = m1 - np2;
    toll = *tol;
    if (toll <= 0.) {

/*        Set the default value of the tolerance for condition tests. */

	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
    }

/*     Determine SVD of D12, D12 = U12 S12 V12', and check if D12 has */
/*     full column rank. V12' is stored in TU. */
/*     Workspace:  need   M2 + NP1*NP1 + max(3*M2+NP1,5*M2); */
/*                 prefer larger. */

    iq = m2 + 1;
    iwrk = iq + np1 * np1;

    i__1 = *ldwork - iwrk + 1;
    dgesvd_("A", "A", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[1], 
	    &dwork[iq], &np1, &tu[tu_offset], ldtu, &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 3;
	return 0;
    }

    rcond[1] = dwork[m2] / dwork[1];
    if (rcond[1] <= toll) {
	rcond[2] = 0.;
	*info = 1;
	return 0;
    }
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Determine Q12. */

    if (nd1 > 0) {
	dlacpy_("Full", &np1, &m2, &dwork[iq], &np1, &d__[(m1 + 1) * d_dim1 + 
		1], ldd, (ftnlen)4);
	dlacpy_("Full", &np1, &nd1, &dwork[iq + np1 * m2], &np1, &dwork[iq], &
		np1, (ftnlen)4);
	dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[
		iq + np1 * nd1], &np1, (ftnlen)4);
    }

/*     Determine Tu by transposing in-situ and scaling. */

    i__1 = m2 - 1;
    for (j = 1; j <= i__1; ++j) {
	dswap_(&j, &tu[j + 1 + tu_dim1], ldtu, &tu[(j + 1) * tu_dim1 + 1], &
		c__1);
/* L10: */
    }

    i__1 = m2;
    for (j = 1; j <= i__1; ++j) {
	d__1 = 1. / dwork[j];
	dscal_(&m2, &d__1, &tu[j * tu_dim1 + 1], &c__1);
/* L20: */
    }

/*     Determine C1 =: Q12'*C1. */
/*     Workspace:  M2 + NP1*NP1 + NP1*N. */

    dgemm_("T", "N", &np1, n, &np1, &c_b15, &dwork[iq], &np1, &c__[c_offset], 
	    ldc, &c_b16, &dwork[iwrk], &np1, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", &np1, n, &dwork[iwrk], &np1, &c__[c_offset], ldc, (ftnlen)
	    4);
/* Computing MAX */
    i__1 = iwrk + np1 * *n - 1;
    lwamax = max(i__1,lwamax);

/*     Determine SVD of D21, D21 = U21 S21 V21', and check if D21 has */
/*     full row rank. U21 is stored in TY. */
/*     Workspace:  need   NP2 + M1*M1 + max(3*NP2+M1,5*NP2); */
/*                 prefer larger. */

    iq = np2 + 1;
    iwrk = iq + m1 * m1;

    i__1 = *ldwork - iwrk + 1;
    dgesvd_("A", "A", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[1], &ty[
	    ty_offset], ldty, &dwork[iq], &m1, &dwork[iwrk], &i__1, &info2, (
	    ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 3;
	return 0;
    }

    rcond[2] = dwork[np2] / dwork[1];
    if (rcond[2] <= toll) {
	*info = 2;
	return 0;
    }
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Determine Q21. */

    if (nd2 > 0) {
	dlacpy_("Full", &np2, &m1, &dwork[iq], &m1, &d__[np1 + 1 + d_dim1], 
		ldd, (ftnlen)4);
	dlacpy_("Full", &nd2, &m1, &dwork[iq + np2], &m1, &dwork[iq], &m1, (
		ftnlen)4);
	dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[iq + 
		nd2], &m1, (ftnlen)4);
    }

/*     Determine Ty by scaling and transposing in-situ. */

    i__1 = np2;
    for (j = 1; j <= i__1; ++j) {
	d__1 = 1. / dwork[j];
	dscal_(&np2, &d__1, &ty[j * ty_dim1 + 1], &c__1);
/* L30: */
    }

    i__1 = np2 - 1;
    for (j = 1; j <= i__1; ++j) {
	dswap_(&j, &ty[j + 1 + ty_dim1], ldty, &ty[(j + 1) * ty_dim1 + 1], &
		c__1);
/* L40: */
    }

/*     Determine B1 =: B1*Q21'. */
/*     Workspace:  NP2 + M1*M1 + N*M1. */

    dgemm_("N", "T", n, &m1, &m1, &c_b15, &b[b_offset], ldb, &dwork[iq], &m1, 
	    &c_b16, &dwork[iwrk], n, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", n, &m1, &dwork[iwrk], n, &b[b_offset], ldb, (ftnlen)4);
/* Computing MAX */
    i__1 = iwrk + *n * m1 - 1;
    lwamax = max(i__1,lwamax);

/*     Determine B2 =: B2*Tu. */
/*     Workspace:  N*M2. */

    dgemm_("N", "N", n, &m2, &m2, &c_b15, &b[(m1 + 1) * b_dim1 + 1], ldb, &tu[
	    tu_offset], ldtu, &c_b16, &dwork[1], n, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", n, &m2, &dwork[1], n, &b[(m1 + 1) * b_dim1 + 1], ldb, (
	    ftnlen)4);

/*     Determine C2 =: Ty*C2. */
/*     Workspace:  NP2*N. */

    dgemm_("N", "N", &np2, n, &np2, &c_b15, &ty[ty_offset], ldty, &c__[np1 + 
	    1 + c_dim1], ldc, &c_b16, &dwork[1], &np2, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", &np2, n, &dwork[1], &np2, &c__[np1 + 1 + c_dim1], ldc, (
	    ftnlen)4);

/*     Determine D22 =: Ty*D22*Tu. */
/*     Workspace:  NP2*M2. */

    dgemm_("N", "N", &np2, &m2, &np2, &c_b15, &ty[ty_offset], ldty, &d__[np1 
	    + 1 + (m1 + 1) * d_dim1], ldd, &c_b16, &dwork[1], &np2, (ftnlen)1,
	     (ftnlen)1);
    dgemm_("N", "N", &np2, &m2, &m2, &c_b15, &dwork[1], &np2, &tu[tu_offset], 
	    ldtu, &c_b16, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, (ftnlen)1, (
	    ftnlen)1);

/* Computing MAX */
    i__1 = *n * max(m2,np2), i__2 = np2 * m2, i__1 = max(i__1,i__2);
    lwamax = max(i__1,lwamax);
    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB10UD *** */
} /* sb10ud_ */

