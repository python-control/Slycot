/* MB03TS.f -- translated by f2c (version 20100827).
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
static doublereal c_b8 = -1.;
static integer c__2 = 2;
static integer c__4 = 4;
static doublereal c_b21 = 0.;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__3 = 3;
static doublereal c_b203 = 1.;

/* Subroutine */ int mb03ts_(logical *isham, logical *wantu, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	u1, integer *ldu1, doublereal *u2, integer *ldu2, integer *j1, 
	integer *n1, integer *n2, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, u1_dim1, u1_offset, u2_dim1, 
	    u2_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal d__[16]	/* was [4][4] */;
    static integer k;
    static doublereal v[3], x[4]	/* was [2][2] */;
    static integer j2, j3, j4;
    static doublereal v1[3], v2[3], a11, a22, a33;
    static integer nd;
    static doublereal cs, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, tau2;
    static logical lblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr2_(char 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen), mb01md_(char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen), 
	    mb01nd_(char *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal dnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dsymv_(char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal xnorm;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasy2_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), dlarfx_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    static doublereal thresh, smlnum;


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

/*     To swap diagonal blocks A11 and A22 of order 1 or 2 in the upper */
/*     quasi-triangular matrix A contained in a skew-Hamiltonian matrix */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G = -G, */
/*                   [  0   A  ] */

/*     or in a Hamiltonian matrix */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G =  G. */
/*                   [  0  -A  ] */

/*     This routine is a modified version of the LAPACK subroutine */
/*     DLAEX2. */

/*     The matrix A must be in Schur canonical form (as returned by the */
/*     LAPACK routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ISHAM   LOGIGAL */
/*             Specifies the type of X: */
/*             = .TRUE.:   X is a Hamiltonian matrix; */
/*             = .FALSE.:  X is a skew-Hamiltonian matrix. */

/*     WANTU   LOGIGAL */
/*             = .TRUE.:   update the matrices U1 and U2 containing the */
/*                         Schur vectors; */
/*             = .FALSE.:  do not update U1 and U2. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A, in Schur */
/*             canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the reordered matrix A, again in Schur canonical form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular part of the symmetric */
/*             matrix G, if ISHAM = .TRUE., or the strictly upper */
/*             triangular part of the skew-symmetric matrix G, otherwise. */
/*             The rest of this array is not referenced. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper or strictly upper triangular part of the */
/*             symmetric or skew-symmetric matrix G, respectively, */
/*             updated by the orthogonal transformation which reorders A. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array must contain the matrix U1. */
/*             On exit, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array contains U1, postmultiplied by the orthogonal */
/*             transformation which reorders A. See the description in */
/*             the SLICOT subroutine MB03TD for further details. */
/*             If WANTU = .FALSE., this array is not referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1. */
/*             LDU1 >= MAX(1,N),  if WANTU = .TRUE.; */
/*             LDU1 >= 1,         otherwise. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array must contain the matrix U2. */
/*             On exit, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array contains U2, postmultiplied by the orthogonal */
/*             transformation which reorders A. */
/*             If WANTU = .FALSE., this array is not referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2. */
/*             LDU2 >= MAX(1,N),  if WANTU = .TRUE.; */
/*             LDU2 >= 1,         otherwise. */

/*     J1      (input) INTEGER */
/*             The index of the first row of the first block A11. */
/*             If J1+N1 < N, then A11 is swapped with the block starting */
/*             at (J1+N1+1)-th diagonal element. */
/*             If J1+N1 > N, then A11 is the last block in A and swapped */
/*             with -A11', if ISHAM = .TRUE., */
/*             or    A11', if ISHAM = .FALSE.. */

/*     N1      (input) INTEGER */
/*             The order of the first block A11. N1 = 0, 1 or 2. */

/*     N2      (input) INTEGER */
/*             The order of the second block A22. N2 = 0, 1 or 2. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  the transformed matrix A would be too far from Schur */
/*                   form; the blocks are not swapped and A, G, U1 and */
/*                   U2 are unchanged. */

/*     REFERENCES */

/*     [1] Bai, Z., and Demmel, J.W. */
/*        On swapping diagonal blocks in real Schur form. */
/*        Linear Algebra Appl., 186, pp. 73-95, 1993. */

/*     [2] Benner, P., Kressner, D., and Mehrmann, V. */
/*         Skew-Hamiltonian and Hamiltonian Eigenvalue Problems: Theory, */
/*         Algorithms and Applications. Techn. Report, TU Berlin, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAEX2). */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Quick return if possible. */

    if (*n == 0 || *n1 == 0 || *n2 == 0) {
	return 0;
    }
    lblk = *j1 + *n1 > *n;

    j2 = *j1 + 1;
    j3 = *j1 + 2;
    j4 = *j1 + 3;

    if (lblk && *n1 == 1) {

	if (*isham) {
	    a11 = a[*n + *n * a_dim1];
	    d__1 = a11 * -2.;
	    dlartg_(&g[*n + *n * g_dim1], &d__1, &cs, &sn, &temp);
	    i__1 = *n - 1;
	    drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], &
		    c__1, &cs, &sn);
	    a[*n + *n * a_dim1] = -a11;
	    if (*wantu) {
		drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1], 
			&c__1, &cs, &sn);
	    }
	} else {
	    i__1 = *n - 1;
	    dswap_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], &
		    c__1);
	    i__1 = *n - 1;
	    dscal_(&i__1, &c_b8, &a[*n * a_dim1 + 1], &c__1);
	    if (*wantu) {
		dswap_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1],
			 &c__1);
		dscal_(n, &c_b8, &u1[*n * u1_dim1 + 1], &c__1);
	    }
	}

    } else if (lblk && *n1 == 2) {

	if (*isham) {

/*           Reorder Hamiltonian matrix: */

/*                      [ A11  G11  ] */
/*                      [         T ]. */
/*                      [  0  -A11  ] */

	    nd = 4;
	    dlacpy_("Full", &c__2, &c__2, &a[*n - 1 + (*n - 1) * a_dim1], lda,
		     d__, &c__4, (ftnlen)4);
	    dlaset_("All", &c__2, &c__2, &c_b21, &c_b21, &d__[2], &c__4, (
		    ftnlen)3);
	    dlacpy_("Upper", &c__2, &c__2, &g[*n - 1 + (*n - 1) * g_dim1], 
		    ldg, &d__[8], &c__4, (ftnlen)5);
	    d__[9] = d__[12];
	    d__[10] = -d__[0];
	    d__[11] = -d__[4];
	    d__[14] = -d__[1];
	    d__[15] = -d__[5];
	    dnorm = dlange_("Max", &nd, &nd, d__, &c__4, &dwork[1], (ftnlen)3)
		    ;

/*           Compute machine-dependent threshold for test for accepting */
/*           swap. */

	    eps = dlamch_("P", (ftnlen)1);
	    smlnum = dlamch_("S", (ftnlen)1) / eps;
/* Computing MAX */
	    d__1 = eps * 40. * dnorm;
	    thresh = max(d__1,smlnum);

/*           Solve A11*X + X*A11' = scale*G11 for X. */

	    dlasy2_(&c_false, &c_false, &c_n1, &c__2, &c__2, d__, &c__4, &d__[
		    10], &c__4, &d__[8], &c__4, &scale, x, &c__2, &xnorm, &
		    ierr);

/*           Compute symplectic QR decomposition of */

/*                  (  -X11  -X12 ) */
/*                  (  -X21  -X22 ). */
/*                  ( scale    0  ) */
/*                  (    0  scale ) */

	    temp = -x[0];
	    dlartg_(&temp, &scale, v1, v2, x);
	    d__1 = -x[1];
	    dlartg_(x, &d__1, &v1[1], &v2[1], &temp);
	    x[2] = -x[2];
	    x[3] = -x[3];
	    x[0] = 0.;
	    x[1] = scale;
	    drot_(&c__1, &x[2], &c__1, x, &c__1, v1, v2);
	    drot_(&c__1, &x[2], &c__1, &x[3], &c__1, &v1[1], &v2[1]);
	    drot_(&c__1, x, &c__1, &x[1], &c__1, &v1[1], &v2[1]);
	    dlartg_(&x[3], &x[1], &v1[2], &v2[2], &temp);

/*           Perform swap provisionally on D. */

	    drot_(&c__4, d__, &c__4, &d__[2], &c__4, v1, v2);
	    drot_(&c__4, d__, &c__4, &d__[1], &c__4, &v1[1], &v2[1]);
	    drot_(&c__4, &d__[2], &c__4, &d__[3], &c__4, &v1[1], &v2[1]);
	    drot_(&c__4, &d__[1], &c__4, &d__[3], &c__4, &v1[2], &v2[2]);
	    drot_(&c__4, d__, &c__1, &d__[8], &c__1, v1, v2);
	    drot_(&c__4, d__, &c__1, &d__[4], &c__1, &v1[1], &v2[1]);
	    drot_(&c__4, &d__[8], &c__1, &d__[12], &c__1, &v1[1], &v2[1]);
	    drot_(&c__4, &d__[4], &c__1, &d__[12], &c__1, &v1[2], &v2[2]);

/*           Test whether to reject swap. */

/* Computing MAX */
	    d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), 
		    d__2 = abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[
		    7]);
	    if (max(d__1,d__2) > thresh) {
		goto L50;
	    }

	    dlacpy_("All", &c__2, &c__2, d__, &c__4, &a[*n - 1 + (*n - 1) * 
		    a_dim1], lda, (ftnlen)3);
	    dlacpy_("Upper", &c__2, &c__2, &d__[8], &c__4, &g[*n - 1 + (*n - 
		    1) * g_dim1], ldg, (ftnlen)5);

	    if (*n > 2) {
		i__1 = *n - 2;
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1, v1, v2);
		i__1 = *n - 2;
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &a[*n * a_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
		i__1 = *n - 2;
		drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
		i__1 = *n - 2;
		drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], 
			&c__1, &v1[2], &v2[2]);
	    }

	    if (*wantu) {
		drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1) * 
			u2_dim1 + 1], &c__1, v1, v2);
		drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u1[*n * u1_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
		drot_(n, &u2[(*n - 1) * u2_dim1 + 1], &c__1, &u2[*n * u2_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
		drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1], 
			&c__1, &v1[2], &v2[2]);
	    }

	} else {

	    if ((d__1 = a[*n - 1 + *n * a_dim1], abs(d__1)) > (d__2 = a[*n + (
		    *n - 1) * a_dim1], abs(d__2))) {
		temp = g[*n - 1 + *n * g_dim1];
		dlartg_(&temp, &a[*n - 1 + *n * a_dim1], &cs, &sn, &g[*n - 1 
			+ *n * g_dim1]);
		sn = -sn;
		i__1 = *n - 2;
		drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], 
			&c__1, &cs, &sn);

		a[*n - 1 + *n * a_dim1] = -sn * a[*n + (*n - 1) * a_dim1];
		temp = -cs * a[*n + (*n - 1) * a_dim1];
		a[*n + (*n - 1) * a_dim1] = g[*n - 1 + *n * g_dim1];
		g[*n - 1 + *n * g_dim1] = temp;
		if (*wantu) {
		    drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 
			    1], &c__1, &cs, &sn);
		}
		i__1 = *n - 2;
		dswap_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1);
		i__1 = *n - 2;
		dscal_(&i__1, &c_b8, &a[(*n - 1) * a_dim1 + 1], &c__1);
		if (*wantu) {
		    dswap_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1)
			     * u2_dim1 + 1], &c__1);
		    dscal_(n, &c_b8, &u1[(*n - 1) * u1_dim1 + 1], &c__1);
		}
	    } else {
		temp = g[*n - 1 + *n * g_dim1];
		dlartg_(&temp, &a[*n + (*n - 1) * a_dim1], &cs, &sn, &g[*n - 
			1 + *n * g_dim1]);
		i__1 = *n - 2;
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1, &cs, &sn);
		a[*n + (*n - 1) * a_dim1] = -sn * a[*n - 1 + *n * a_dim1];
		a[*n - 1 + *n * a_dim1] = cs * a[*n - 1 + *n * a_dim1];
		if (*wantu) {
		    drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1) 
			    * u2_dim1 + 1], &c__1, &cs, &sn);
		}
		i__1 = *n - 1;
		dswap_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1],
			 &c__1);
		i__1 = *n - 1;
		dscal_(&i__1, &c_b8, &a[*n * a_dim1 + 1], &c__1);
		if (*wantu) {
		    dswap_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 
			    + 1], &c__1);
		    dscal_(n, &c_b8, &u1[*n * u1_dim1 + 1], &c__1);
		}
	    }
	}

/*        Standardize new 2-by-2 block. */

	dlanv2_(&a[*n - 1 + (*n - 1) * a_dim1], &a[*n - 1 + *n * a_dim1], &a[*
		n + (*n - 1) * a_dim1], &a[*n + *n * a_dim1], &wr1, &wi1, &
		wr2, &wi2, &cs, &sn);
	i__1 = *n - 2;
	drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &a[*n * a_dim1 + 1], &
		c__1, &cs, &sn);
	if (*isham) {
	    temp = g[*n - 1 + *n * g_dim1];
	    i__1 = *n - 1;
	    drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 + 1]
		    , &c__1, &cs, &sn);
	    tau = cs * temp + sn * g[*n + *n * g_dim1];
	    g[*n + *n * g_dim1] = cs * g[*n + *n * g_dim1] - sn * temp;
	    g[*n - 1 + (*n - 1) * g_dim1] = cs * g[*n - 1 + (*n - 1) * g_dim1]
		     + sn * tau;
	    drot_(&c__1, &g[*n - 1 + *n * g_dim1], ldg, &g[*n + *n * g_dim1], 
		    ldg, &cs, &sn);
	} else {
	    i__1 = *n - 2;
	    drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 + 1]
		    , &c__1, &cs, &sn);
	}
	if (*wantu) {
	    drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u1[*n * u1_dim1 + 1]
		    , &c__1, &cs, &sn);
	    drot_(n, &u2[(*n - 1) * u2_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1]
		    , &c__1, &cs, &sn);
	}

    } else if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

	a11 = a[*j1 + *j1 * a_dim1];
	a22 = a[j2 + j2 * a_dim1];

/*        Determine the transformation to perform the interchange. */

	d__1 = a22 - a11;
	dlartg_(&a[*j1 + j2 * a_dim1], &d__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix A. */

	if (j3 <= *n) {
	    i__1 = *n - *j1 - 1;
	    drot_(&i__1, &a[*j1 + j3 * a_dim1], lda, &a[j2 + j3 * a_dim1], 
		    lda, &cs, &sn);
	}
	i__1 = *j1 - 1;
	drot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[j2 * a_dim1 + 1], &c__1, 
		&cs, &sn);

	a[*j1 + *j1 * a_dim1] = a22;
	a[j2 + j2 * a_dim1] = a11;

/*        Apply transformation to the matrix G. */

	if (*isham) {
	    temp = g[*j1 + j2 * g_dim1];
	    drot_(j1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &c__1,
		     &cs, &sn);
	    tau = cs * temp + sn * g[j2 + j2 * g_dim1];
	    g[j2 + j2 * g_dim1] = cs * g[j2 + j2 * g_dim1] - sn * temp;
	    g[*j1 + *j1 * g_dim1] = cs * g[*j1 + *j1 * g_dim1] + sn * tau;
	    i__1 = *n - *j1;
	    drot_(&i__1, &g[*j1 + j2 * g_dim1], ldg, &g[j2 + j2 * g_dim1], 
		    ldg, &cs, &sn);
	} else {
	    if (*n > *j1 + 1) {
		i__1 = *n - *j1 - 1;
		drot_(&i__1, &g[*j1 + (*j1 + 2) * g_dim1], ldg, &g[j2 + (*j1 
			+ 2) * g_dim1], ldg, &cs, &sn);
	    }
	    i__1 = *j1 - 1;
	    drot_(&i__1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

	    drot_(n, &u1[*j1 * u1_dim1 + 1], &c__1, &u1[j2 * u1_dim1 + 1], &
		    c__1, &cs, &sn);
	    drot_(n, &u2[*j1 * u2_dim1 + 1], &c__1, &u2[j2 * u2_dim1 + 1], &
		    c__1, &cs, &sn);
	}

    } else {

/*        Swapping involves at least one 2-by-2 block. */

/*        Copy the diagonal block of order N1+N2 to the local array D */
/*        and compute its norm. */

	nd = *n1 + *n2;
	dlacpy_("Full", &nd, &nd, &a[*j1 + *j1 * a_dim1], lda, d__, &c__4, (
		ftnlen)4);
	dnorm = dlange_("Max", &nd, &nd, d__, &c__4, &dwork[1], (ftnlen)3);

/*        Compute machine-dependent threshold for test for accepting */
/*        swap. */

	eps = dlamch_("P", (ftnlen)1);
	smlnum = dlamch_("S", (ftnlen)1) / eps;
/* Computing MAX */
	d__1 = eps * 30. * dnorm;
	thresh = max(d__1,smlnum);

/*        Solve A11*X - X*A22 = scale*A12 for X. */

	dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 + 
		(*n1 + 1 << 2) - 5], &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

	k = *n1 + *n1 + *n2 - 3;
	switch (k) {
	    case 1:  goto L10;
	    case 2:  goto L20;
	    case 3:  goto L30;
	}

L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that: */

/*        ( scale, X11, X12 ) H = ( 0, 0, * ). */

	v[0] = scale;
	v[1] = x[0];
	v[2] = x[2];
	dlarfg_(&c__3, &v[2], v, &c__1, &tau);
	v[2] = 1.;
	a11 = a[*j1 + *j1 * a_dim1];

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("Left", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)
		4);
	dlarfx_("Right", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (
		ftnlen)5);

/*        Test whether to reject swap. */

/* Computing MAX */
	d__2 = abs(d__[2]), d__3 = abs(d__[6]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[10] - a11, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix A. */

	i__1 = *n - *j1 + 1;
	dlarfx_("Left", &c__3, &i__1, v, &tau, &a[*j1 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)4);
	dlarfx_("Right", &j2, &c__3, v, &tau, &a[*j1 * a_dim1 + 1], lda, &
		dwork[1], (ftnlen)5);

	a[j3 + *j1 * a_dim1] = 0.;
	a[j3 + j2 * a_dim1] = 0.;
	a[j3 + j3 * a_dim1] = a11;

/*        Apply transformation to G. */

	if (*isham) {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
	    dsymv_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    temp = tau * -.5 * ddot_(&c__3, &dwork[1], &c__1, v, &c__1);
	    daxpy_(&c__3, &temp, v, &c__1, &dwork[1], &c__1);
	    dsyr2_("Upper", &c__3, &c_b8, v, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	} else {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
	    mb01md_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    mb01nd_("Upper", &c__3, &c_b203, v, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	}

	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

	    dlarfx_("R", n, &c__3, v, &tau, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
	    dlarfx_("R", n, &c__3, v, &tau, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
	}
	goto L40;

L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that: */

/*        H (  -X11 ) = ( * ) */
/*          (  -X21 ) = ( 0 ). */
/*          ( scale ) = ( 0 ) */

	v[0] = -x[0];
	v[1] = -x[1];
	v[2] = scale;
	dlarfg_(&c__3, v, &v[1], &c__1, &tau);
	v[0] = 1.;
	a33 = a[j3 + j3 * a_dim1];

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("L", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)1);
	dlarfx_("R", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
	d__2 = abs(d__[1]), d__3 = abs(d__[2]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[0] - a33, abs(d__1));
	if (max(d__2,d__3) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix A. */

	dlarfx_("Right", &j3, &c__3, v, &tau, &a[*j1 * a_dim1 + 1], lda, &
		dwork[1], (ftnlen)5);
	i__1 = *n - *j1;
	dlarfx_("Left", &c__3, &i__1, v, &tau, &a[*j1 + j2 * a_dim1], lda, &
		dwork[1], (ftnlen)4);

	a[*j1 + *j1 * a_dim1] = a33;
	a[j2 + *j1 * a_dim1] = 0.;
	a[j3 + *j1 * a_dim1] = 0.;

/*        Apply transformation to G. */

	if (*isham) {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
	    dsymv_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    temp = tau * -.5 * ddot_(&c__3, &dwork[1], &c__1, v, &c__1);
	    daxpy_(&c__3, &temp, v, &c__1, &dwork[1], &c__1);
	    dsyr2_("Upper", &c__3, &c_b8, v, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	} else {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
	    mb01md_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    mb01nd_("Upper", &c__3, &c_b203, v, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	}

	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

	    dlarfx_("R", n, &c__3, v, &tau, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
	    dlarfx_("R", n, &c__3, v, &tau, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
	}
	goto L40;

L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so */
/*        that: */

/*        H(2) H(1) (  -X11  -X12 ) = (  *  * ) */
/*                  (  -X21  -X22 )   (  0  * ). */
/*                  ( scale    0  )   (  0  0 ) */
/*                  (    0  scale )   (  0  0 ) */

	v1[0] = -x[0];
	v1[1] = -x[1];
	v1[2] = scale;
	dlarfg_(&c__3, v1, &v1[1], &c__1, &tau1);
	v1[0] = 1.;

	temp = -tau1 * (x[2] + v1[1] * x[3]);
	v2[0] = -temp * v1[1] - x[3];
	v2[1] = -temp * v1[2];
	v2[2] = scale;
	dlarfg_(&c__3, v2, &v2[1], &c__1, &tau2);
	v2[0] = 1.;

/*        Perform swap provisionally on diagonal block in D. */

	dlarfx_("L", &c__3, &c__4, v1, &tau1, d__, &c__4, &dwork[1], (ftnlen)
		1);
	dlarfx_("R", &c__4, &c__3, v1, &tau1, d__, &c__4, &dwork[1], (ftnlen)
		1);
	dlarfx_("L", &c__3, &c__4, v2, &tau2, &d__[1], &c__4, &dwork[1], (
		ftnlen)1);
	dlarfx_("R", &c__4, &c__3, v2, &tau2, &d__[4], &c__4, &dwork[1], (
		ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
	d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), d__2 = 
		abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[7]);
	if (max(d__1,d__2) > thresh) {
	    goto L50;
	}

/*        Accept swap: apply transformation to the entire matrix A. */

	i__1 = *n - *j1 + 1;
	dlarfx_("L", &c__3, &i__1, v1, &tau1, &a[*j1 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)1);
	dlarfx_("R", &j4, &c__3, v1, &tau1, &a[*j1 * a_dim1 + 1], lda, &dwork[
		1], (ftnlen)1);
	i__1 = *n - *j1 + 1;
	dlarfx_("L", &c__3, &i__1, v2, &tau2, &a[j2 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)1);
	dlarfx_("R", &j4, &c__3, v2, &tau2, &a[j2 * a_dim1 + 1], lda, &dwork[
		1], (ftnlen)1);

	a[j3 + *j1 * a_dim1] = 0.;
	a[j3 + j2 * a_dim1] = 0.;
	a[j4 + *j1 * a_dim1] = 0.;
	a[j4 + j2 * a_dim1] = 0.;

/*        Apply transformation to G. */

	if (*isham) {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v1, &tau1, &g[*j1 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
	    dsymv_("Upper", &c__3, &tau1, &g[*j1 + *j1 * g_dim1], ldg, v1, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    temp = tau1 * -.5 * ddot_(&c__3, &dwork[1], &c__1, v1, &c__1);
	    daxpy_(&c__3, &temp, v1, &c__1, &dwork[1], &c__1);
	    dsyr2_("Upper", &c__3, &c_b8, v1, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v1, &tau1, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }

	    i__1 = j2 - 1;
	    dlarfx_("Right", &i__1, &c__3, v2, &tau2, &g[j2 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
	    dsymv_("Upper", &c__3, &tau2, &g[j2 + j2 * g_dim1], ldg, v2, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    temp = tau2 * -.5 * ddot_(&c__3, &dwork[1], &c__1, v2, &c__1);
	    daxpy_(&c__3, &temp, v2, &c__1, &dwork[1], &c__1);
	    dsyr2_("Upper", &c__3, &c_b8, v2, &c__1, &dwork[1], &c__1, &g[j2 
		    + j2 * g_dim1], ldg, (ftnlen)5);
	    if (*n > j2 + 2) {
		i__1 = *n - j2 - 2;
		dlarfx_("Left", &c__3, &i__1, v2, &tau2, &g[j2 + (j2 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	} else {
	    i__1 = *j1 - 1;
	    dlarfx_("Right", &i__1, &c__3, v1, &tau1, &g[*j1 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
	    mb01md_("Upper", &c__3, &tau1, &g[*j1 + *j1 * g_dim1], ldg, v1, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    mb01nd_("Upper", &c__3, &c_b203, v1, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
	    if (*n > *j1 + 2) {
		i__1 = *n - *j1 - 2;
		dlarfx_("Left", &c__3, &i__1, v1, &tau1, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	    i__1 = j2 - 1;
	    dlarfx_("Right", &i__1, &c__3, v2, &tau2, &g[j2 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
	    mb01md_("Upper", &c__3, &tau2, &g[j2 + j2 * g_dim1], ldg, v2, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
	    mb01nd_("Upper", &c__3, &c_b203, v2, &c__1, &dwork[1], &c__1, &g[
		    j2 + j2 * g_dim1], ldg, (ftnlen)5);
	    if (*n > j2 + 2) {
		i__1 = *n - j2 - 2;
		dlarfx_("Left", &c__3, &i__1, v2, &tau2, &g[j2 + (j2 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
	    }
	}

	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

	    dlarfx_("R", n, &c__3, v1, &tau1, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
	    dlarfx_("R", n, &c__3, v2, &tau2, &u1[j2 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
	    dlarfx_("R", n, &c__3, v1, &tau1, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
	    dlarfx_("R", n, &c__3, v2, &tau2, &u2[j2 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
	}

L40:

	if (*n2 == 2) {

/*           Standardize new 2-by-2 block A11. */

	    dlanv2_(&a[*j1 + *j1 * a_dim1], &a[*j1 + j2 * a_dim1], &a[j2 + *
		    j1 * a_dim1], &a[j2 + j2 * a_dim1], &wr1, &wi1, &wr2, &
		    wi2, &cs, &sn);
	    i__1 = *n - *j1 - 1;
	    drot_(&i__1, &a[*j1 + (*j1 + 2) * a_dim1], lda, &a[j2 + (*j1 + 2) 
		    * a_dim1], lda, &cs, &sn);
	    i__1 = *j1 - 1;
	    drot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[j2 * a_dim1 + 1], &
		    c__1, &cs, &sn);
	    if (*isham) {
		temp = g[*j1 + j2 * g_dim1];
		drot_(j1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &
			c__1, &cs, &sn);
		tau = cs * temp + sn * g[j2 + j2 * g_dim1];
		g[j2 + j2 * g_dim1] = cs * g[j2 + j2 * g_dim1] - sn * temp;
		g[*j1 + *j1 * g_dim1] = cs * g[*j1 + *j1 * g_dim1] + sn * tau;
		i__1 = *n - *j1;
		drot_(&i__1, &g[*j1 + j2 * g_dim1], ldg, &g[j2 + j2 * g_dim1],
			 ldg, &cs, &sn);
	    } else {
		if (*n > *j1 + 1) {
		    i__1 = *n - *j1 - 1;
		    drot_(&i__1, &g[*j1 + (*j1 + 2) * g_dim1], ldg, &g[j2 + (*
			    j1 + 2) * g_dim1], ldg, &cs, &sn);
		}
		i__1 = *j1 - 1;
		drot_(&i__1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1],
			 &c__1, &cs, &sn);
	    }
	    if (*wantu) {
		drot_(n, &u1[*j1 * u1_dim1 + 1], &c__1, &u1[j2 * u1_dim1 + 1],
			 &c__1, &cs, &sn);
		drot_(n, &u2[*j1 * u2_dim1 + 1], &c__1, &u2[j2 * u2_dim1 + 1],
			 &c__1, &cs, &sn);
	    }
	}

	if (*n1 == 2) {

/*           Standardize new 2-by-2 block A22. */

	    j3 = *j1 + *n2;
	    j4 = j3 + 1;
	    dlanv2_(&a[j3 + j3 * a_dim1], &a[j3 + j4 * a_dim1], &a[j4 + j3 * 
		    a_dim1], &a[j4 + j4 * a_dim1], &wr1, &wi1, &wr2, &wi2, &
		    cs, &sn);
	    if (j3 + 2 <= *n) {
		i__1 = *n - j3 - 1;
		drot_(&i__1, &a[j3 + (j3 + 2) * a_dim1], lda, &a[j4 + (j3 + 2)
			 * a_dim1], lda, &cs, &sn);
	    }
	    i__1 = j3 - 1;
	    drot_(&i__1, &a[j3 * a_dim1 + 1], &c__1, &a[j4 * a_dim1 + 1], &
		    c__1, &cs, &sn);
	    if (*isham) {
		temp = g[j3 + j4 * g_dim1];
		drot_(&j3, &g[j3 * g_dim1 + 1], &c__1, &g[j4 * g_dim1 + 1], &
			c__1, &cs, &sn);
		tau = cs * temp + sn * g[j4 + j4 * g_dim1];
		g[j4 + j4 * g_dim1] = cs * g[j4 + j4 * g_dim1] - sn * temp;
		g[j3 + j3 * g_dim1] = cs * g[j3 + j3 * g_dim1] + sn * tau;
		i__1 = *n - j3;
		drot_(&i__1, &g[j3 + j4 * g_dim1], ldg, &g[j4 + j4 * g_dim1], 
			ldg, &cs, &sn);
	    } else {
		if (*n > j3 + 1) {
		    i__1 = *n - j3 - 1;
		    drot_(&i__1, &g[j3 + (j3 + 2) * g_dim1], ldg, &g[j4 + (j3 
			    + 2) * g_dim1], ldg, &cs, &sn);
		}
		i__1 = j3 - 1;
		drot_(&i__1, &g[j3 * g_dim1 + 1], &c__1, &g[j4 * g_dim1 + 1], 
			&c__1, &cs, &sn);
	    }
	    if (*wantu) {
		drot_(n, &u1[j3 * u1_dim1 + 1], &c__1, &u1[j4 * u1_dim1 + 1], 
			&c__1, &cs, &sn);
		drot_(n, &u2[j3 * u2_dim1 + 1], &c__1, &u2[j4 * u2_dim1 + 1], 
			&c__1, &cs, &sn);
	    }
	}

    }
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

L50:
    *info = 1;
    return 0;
/* *** Last line of MB03TS *** */
} /* mb03ts_ */

