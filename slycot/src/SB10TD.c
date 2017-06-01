/* SB10TD.f -- translated by f2c (version 20100827).
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

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;
static doublereal c_b39 = -1.;

/* Subroutine */ int sb10td_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *d__, integer *ldd, doublereal *tu, 
	integer *ldtu, doublereal *ty, integer *ldty, doublereal *ak, integer 
	*ldak, doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *rcond, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer ak_dim1, ak_offset, bk_dim1, bk_offset, ck_dim1, ck_offset, 
	    d_dim1, d_offset, dk_dim1, dk_offset, tu_dim1, tu_offset, ty_dim1,
	     ty_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m1, m2, np1, np2;
    static doublereal toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal anorm;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer minwrk;


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

/*     To compute the matrices of the H2 optimal discrete-time controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     from the matrices of the controller for the normalized system, */
/*     as determined by the SLICOT Library routine SB10SD. */

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

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. Only the trailing */
/*             NMEAS-by-NCON submatrix D22 is used. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

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

/*     AK      (input/output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain controller state matrix for the normalized system */
/*             as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading N-by-N part of this array contains */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDBK,NMEAS) */
/*             On entry, the leading N-by-NMEAS part of this array must */
/*             contain controller input matrix for the normalized system */
/*             as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading N-by-NMEAS part of this array */
/*             contains controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (input/output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             On entry, the leading NCON-by-N part of this array must */
/*             contain controller output matrix for the normalized */
/*             system as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading NCON-by-N part of this array contains */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDDK,NMEAS) */
/*             On entry, the leading NCON-by-NMEAS part of this array */
/*             must contain controller matrix DK for the normalized */
/*             system as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading NCON-by-NMEAS part of this array */
/*             contains controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             RCOND contains an estimate of the reciprocal condition */
/*             number of the matrix Im2 + DKHAT*D22 which must be */
/*             inverted in the computation of the controller. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used in determining the nonsingularity of the */
/*             matrix which must be inverted. If TOL <= 0, then a default */
/*             value equal to sqrt(EPS) is used, where EPS is the */
/*             relative machine precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(N*M2,N*NP2,M2*NP2,M2*M2+4*M2). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix Im2 + DKHAT*D22 is singular, or the */
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
/*     input and output transformations and of the matrix Im2 + */
/*     DKHAT*D22. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Jan. 2000. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    tu_dim1 = *ldtu;
    tu_offset = 1 + tu_dim1;
    tu -= tu_offset;
    ty_dim1 = *ldty;
    ty_offset = 1 + ty_dim1;
    ty -= ty_offset;
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
    } else if (*ldd < max(1,*np)) {
	*info = -7;
    } else if (*ldtu < max(1,m2)) {
	*info = -9;
    } else if (*ldty < max(1,np2)) {
	*info = -11;
    } else if (*ldak < max(1,*n)) {
	*info = -13;
    } else if (*ldbk < max(1,*n)) {
	*info = -15;
    } else if (*ldck < max(1,m2)) {
	*info = -17;
    } else if (*lddk < max(1,m2)) {
	*info = -19;
    } else {

/*        Compute workspace. */

/* Computing MAX */
	i__1 = *n * m2, i__2 = *n * np2, i__1 = max(i__1,i__2), i__2 = m2 * 
		np2, i__1 = max(i__1,i__2), i__2 = m2 * (m2 + 4);
	minwrk = max(i__1,i__2);
	if (*ldwork < minwrk) {
	    *info = -24;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10TD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
	*rcond = 1.;
	return 0;
    }

    toll = *tol;
    if (toll <= 0.) {

/*        Set the default value of the tolerance for nonsingularity test. */

	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
    }

/*     Find BKHAT . */

    dgemm_("N", "N", n, &np2, &np2, &c_b6, &bk[bk_offset], ldbk, &ty[
	    ty_offset], ldty, &c_b7, &dwork[1], n, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", n, &np2, &dwork[1], n, &bk[bk_offset], ldbk, (ftnlen)4);

/*     Find CKHAT . */

    dgemm_("N", "N", &m2, n, &m2, &c_b6, &tu[tu_offset], ldtu, &ck[ck_offset],
	     ldck, &c_b7, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", &m2, n, &dwork[1], &m2, &ck[ck_offset], ldck, (ftnlen)4);

/*     Compute DKHAT . */

    dgemm_("N", "N", &m2, &np2, &m2, &c_b6, &tu[tu_offset], ldtu, &dk[
	    dk_offset], lddk, &c_b7, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &m2, &np2, &np2, &c_b6, &dwork[1], &m2, &ty[ty_offset], 
	    ldty, &c_b7, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute Im2 + DKHAT*D22 . */

    iwrk = m2 * m2 + 1;
    dlaset_("Full", &m2, &m2, &c_b7, &c_b6, &dwork[1], &m2, (ftnlen)4);
    dgemm_("N", "N", &m2, &m2, &np2, &c_b6, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b6, &dwork[1], &m2, (ftnlen)1, (
	    ftnlen)1);
    anorm = dlange_("1", &m2, &m2, &dwork[1], &m2, &dwork[iwrk], (ftnlen)1);
    dgetrf_(&m2, &m2, &dwork[1], &m2, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }
    dgecon_("1", &m2, &dwork[1], &m2, &anorm, rcond, &dwork[iwrk], &iwork[m2 
	    + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

    if (*rcond < toll) {
	*info = 1;
	return 0;
    }

/*     Compute CK . */

    dgetrs_("N", &m2, n, &dwork[1], &m2, &iwork[1], &ck[ck_offset], ldck, &
	    info2, (ftnlen)1);

/*     Compute DK . */

    dgetrs_("N", &m2, &np2, &dwork[1], &m2, &iwork[1], &dk[dk_offset], lddk, &
	    info2, (ftnlen)1);

/*     Compute AK . */

    dgemm_("N", "N", n, &m2, &np2, &c_b6, &bk[bk_offset], ldbk, &d__[np1 + 1 
	    + (m1 + 1) * d_dim1], ldd, &c_b7, &dwork[1], n, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("N", "N", n, n, &m2, &c_b39, &dwork[1], n, &ck[ck_offset], ldck, &
	    c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute BK . */

    dgemm_("N", "N", n, &np2, &m2, &c_b39, &dwork[1], n, &dk[dk_offset], lddk,
	     &c_b6, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);
    return 0;
/* *** Last line of SB10TD *** */
} /* sb10td_ */

