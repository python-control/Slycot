/* SB10LD.f -- translated by f2c (version 20100827).
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
static doublereal c_b9 = -1.;

/* Subroutine */ int sb10ld_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *ac, integer *ldac, doublereal *bc, integer *ldbc, 
	doublereal *cc, integer *ldcc, doublereal *dc, integer *lddc, integer 
	*iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ac_dim1, ac_offset, ak_dim1, ak_offset, b_dim1, 
	    b_offset, bc_dim1, bc_offset, bk_dim1, bk_offset, c_dim1, 
	    c_offset, cc_dim1, cc_offset, ck_dim1, ck_offset, d_dim1, 
	    d_offset, dc_dim1, dc_offset, dk_dim1, dk_offset, i__1;

    /* Local variables */
    static integer m1, m2, n2, np1, np2, iw2, iw3, iw4, iw5, iw6, iw7, iw8;
    static doublereal eps;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond, anorm;
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
	    integer *), xerbla_(char *, integer *, ftnlen);
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

/*     To compute the matrices of the closed-loop system */

/*              | AC | BC | */
/*          G = |----|----|, */
/*              | CC | DC | */

/*     from the matrices of the open-loop system */

/*               | A | B | */
/*           P = |---|---| */
/*               | C | D | */

/*     and the matrices of the controller */

/*              | AK | BK | */
/*          K = |----|----|. */
/*              | CK | DK | */

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

/*     AK      (input) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array must contain the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (input) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array must contain the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (input) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array must contain the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (input) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array must contain */
/*             the controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     AC      (output) DOUBLE PRECISION array, dimension (LDAC,2*N) */
/*             The leading 2*N-by-2*N part of this array contains the */
/*             closed-loop system state matrix AC. */

/*     LDAC    INTEGER */
/*             The leading dimension of the array AC. */
/*             LDAC >= max(1,2*N). */

/*     BC      (output) DOUBLE PRECISION array, dimension (LDBC,M-NCON) */
/*             The leading 2*N-by-(M-NCON) part of this array contains */
/*             the closed-loop system input matrix BC. */

/*     LDBC    INTEGER */
/*             The leading dimension of the array BC. */
/*             LDBC >= max(1,2*N). */

/*     CC      (output) DOUBLE PRECISION array, dimension (LDCC,2*N) */
/*             The leading (NP-NMEAS)-by-2*N part of this array contains */
/*             the closed-loop system output matrix CC. */

/*     LDCC    INTEGER */
/*             The leading dimension of the array CC. */
/*             LDCC >= max(1,NP-NMEAS). */

/*     DC      (output) DOUBLE PRECISION array, dimension (LDDC,M-NCON) */
/*             The leading (NP-NMEAS)-by-(M-NCON) part of this array */
/*             contains the closed-loop system input/output matrix DC. */

/*     LDDC    INTEGER */
/*             The leading dimension of the array DC. */
/*             LDDC >= max(1,NP-NMEAS). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*max(NCON,NMEAS) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 2*M*M+NP*NP+2*M*N+M*NP+2*N*NP. */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicactor */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix Inp2 - D22*DK is singular to working */
/*                   precision; */
/*             = 2:  if the matrix Im2 - DK*D22 is singular to working */
/*                   precision. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

/*     REFERENCES */

/*     [1] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     matrices  Inp2 - D22*DK  and  Im2 - DK*D22. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999. */
/*     A. Markovski, Technical University, Sofia, April, 2003. */

/*     KEYWORDS */

/*     Closed loop systems, feedback control, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */

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
    ac_dim1 = *ldac;
    ac_offset = 1 + ac_dim1;
    ac -= ac_offset;
    bc_dim1 = *ldbc;
    bc_offset = 1 + bc_dim1;
    bc -= bc_offset;
    cc_dim1 = *ldcc;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    dc_dim1 = *lddc;
    dc_offset = 1 + dc_dim1;
    dc -= dc_offset;
    --iwork;
    --dwork;

    /* Function Body */
    n2 = *n << 1;
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
    } else if (*ldac < max(1,n2)) {
	*info = -23;
    } else if (*ldbc < max(1,n2)) {
	*info = -25;
    } else if (*ldcc < max(1,np1)) {
	*info = -27;
    } else if (*lddc < max(1,np1)) {
	*info = -29;
    } else {

/*        Compute workspace. */

	minwrk = (*m << 1) * *m + *np * *np + (*m << 1) * *n + *m * *np + (*n 
		<< 1) * *np;
	if (*ldwork < minwrk) {
	    *info = -32;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB10LD", &i__1, (ftnlen)6);
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

    iw2 = np2 * np2 + 1;
    iw3 = iw2 + m2 * m2;
    iw4 = iw3 + np2 * *n;
    iw5 = iw4 + m2 * *n;
    iw6 = iw5 + np2 * m1;
    iw7 = iw6 + m2 * m1;
    iw8 = iw7 + m2 * *n;
    iwrk = iw8 + np2 * *n;

/*     Compute inv(Inp2 - D22*DK) . */

    dlaset_("Full", &np2, &np2, &c_b5, &c_b6, &dwork[1], &np2, (ftnlen)4);
    dgemm_("N", "N", &np2, &np2, &m2, &c_b9, &d__[np1 + 1 + (m1 + 1) * d_dim1]
	    , ldd, &dk[dk_offset], lddk, &c_b6, &dwork[1], &np2, (ftnlen)1, (
	    ftnlen)1);
    anorm = dlange_("1", &np2, &np2, &dwork[1], &np2, &dwork[iwrk], (ftnlen)1)
	    ;
    dgetrf_(&np2, &np2, &dwork[1], &np2, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 1;
	return 0;
    }
    dgecon_("1", &np2, &dwork[1], &np2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    np2 + 1], info, (ftnlen)1);
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Return if the matrix is singular to working precision. */

    if (rcond < eps) {
	*info = 1;
	return 0;
    }
    i__1 = *ldwork - iwrk + 1;
    dgetri_(&np2, &dwork[1], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Compute inv(Im2 - DK*D22) . */

    dlaset_("Full", &m2, &m2, &c_b5, &c_b6, &dwork[iw2], &m2, (ftnlen)4);
    dgemm_("N", "N", &m2, &m2, &np2, &c_b9, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b6, &dwork[iw2], &m2, (ftnlen)1, (
	    ftnlen)1);
    anorm = dlange_("1", &m2, &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1);
    dgetrf_(&m2, &m2, &dwork[iw2], &m2, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 2;
	return 0;
    }
    dgecon_("1", &m2, &dwork[iw2], &m2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    m2 + 1], info, (ftnlen)1);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Return if the matrix is singular to working precision. */

    if (rcond < eps) {
	*info = 2;
	return 0;
    }
    i__1 = *ldwork - iwrk + 1;
    dgetri_(&m2, &dwork[iw2], &m2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    lwamax = max(i__1,lwamax);

/*     Compute inv(Inp2 - D22*DK)*C2 . */

    dgemm_("N", "N", &np2, n, &np2, &c_b6, &dwork[1], &np2, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b5, &dwork[iw3], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute DK*inv(Inp2 - D22*DK)*C2 . */

    dgemm_("N", "N", &m2, n, &np2, &c_b6, &dk[dk_offset], lddk, &dwork[iw3], &
	    np2, &c_b5, &dwork[iw4], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute inv(Inp2 - D22*DK)*D21 . */

    dgemm_("N", "N", &np2, &m1, &np2, &c_b6, &dwork[1], &np2, &d__[np1 + 1 + 
	    d_dim1], ldd, &c_b5, &dwork[iw5], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute DK*inv(Inp2 - D22*DK)*D21 . */

    dgemm_("N", "N", &m2, &m1, &np2, &c_b6, &dk[dk_offset], lddk, &dwork[iw5],
	     &np2, &c_b5, &dwork[iw6], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute inv(Im2 - DK*D22)*CK . */

    dgemm_("N", "N", &m2, n, &m2, &c_b6, &dwork[iw2], &m2, &ck[ck_offset], 
	    ldck, &c_b5, &dwork[iw7], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute D22*inv(Im2 - DK*D22)*CK . */

    dgemm_("N", "N", &np2, n, &m2, &c_b6, &d__[np1 + 1 + (m1 + 1) * d_dim1], 
	    ldd, &dwork[iw7], &m2, &c_b5, &dwork[iw8], &np2, (ftnlen)1, (
	    ftnlen)1);

/*     Compute AC . */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ac[ac_offset], ldac, (ftnlen)4);
    dgemm_("N", "N", n, n, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[
	    iw4], &m2, &c_b6, &ac[ac_offset], ldac, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", n, n, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[
	    iw7], &m2, &c_b5, &ac[(*n + 1) * ac_dim1 + 1], ldac, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("N", "N", n, n, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw3], &
	    np2, &c_b5, &ac[*n + 1 + ac_dim1], ldac, (ftnlen)1, (ftnlen)1);
    dlacpy_("Full", n, n, &ak[ak_offset], ldak, &ac[*n + 1 + (*n + 1) * 
	    ac_dim1], ldac, (ftnlen)4);
    dgemm_("N", "N", n, n, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw8], &
	    np2, &c_b6, &ac[*n + 1 + (*n + 1) * ac_dim1], ldac, (ftnlen)1, (
	    ftnlen)1);

/*     Compute BC . */

    dlacpy_("Full", n, &m1, &b[b_offset], ldb, &bc[bc_offset], ldbc, (ftnlen)
	    4);
    dgemm_("N", "N", n, &m1, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &
	    dwork[iw6], &m2, &c_b6, &bc[bc_offset], ldbc, (ftnlen)1, (ftnlen)
	    1);
    dgemm_("N", "N", n, &m1, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw5], &
	    np2, &c_b5, &bc[*n + 1 + bc_dim1], ldbc, (ftnlen)1, (ftnlen)1);

/*     Compute CC . */

    dlacpy_("Full", &np1, n, &c__[c_offset], ldc, &cc[cc_offset], ldcc, (
	    ftnlen)4);
    dgemm_("N", "N", &np1, n, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, &
	    dwork[iw4], &m2, &c_b6, &cc[cc_offset], ldcc, (ftnlen)1, (ftnlen)
	    1);
    dgemm_("N", "N", &np1, n, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, &
	    dwork[iw7], &m2, &c_b5, &cc[(*n + 1) * cc_dim1 + 1], ldcc, (
	    ftnlen)1, (ftnlen)1);

/*     Compute DC . */

    dlacpy_("Full", &np1, &m1, &d__[d_offset], ldd, &dc[dc_offset], lddc, (
	    ftnlen)4);
    dgemm_("N", "N", &np1, &m1, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, 
	    &dwork[iw6], &m2, &c_b6, &dc[dc_offset], lddc, (ftnlen)1, (ftnlen)
	    1);

    return 0;
/* *** Last line of SB10LD *** */
} /* sb10ld_ */

