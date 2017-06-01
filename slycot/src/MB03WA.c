/* MB03WA.f -- translated by f2c (version 20100827).
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

static integer c__4 = 4;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b42 = 1.;
static doublereal c_b48 = -1.;

/* Subroutine */ int mb03wa_(logical *wantq, logical *wantz, integer *n1, 
	integer *n2, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	 doublereal *q, integer *ldq, doublereal *z__, integer *ldz, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer i__, m;
    static doublereal s[16]	/* was [4][4] */, t[16]	/* was [4][4] */, be[
	    2], ai[2], ar[2], sa, sb, li[16]	/* was [4][4] */, ir[16]	
	    /* was [4][4] */, ss, ws, eps;
    static logical weak;
    static doublereal ddum, taul[4], dsum;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal taur[4], scpy[16]	/* was [4][4] */, tcpy[16]	/* 
	    was [4][4] */;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale, bqra21, brqa21;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal licop[16]	/* was [4][4] */;
    static integer linfo;
    static doublereal ircop[16]	/* was [4][4] */;
    extern /* Subroutine */ int mb03yt_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal dnorm;
    extern /* Subroutine */ int sb04ow_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal dwork[32];
    static integer iwork[4];
    extern /* Subroutine */ int dgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dgerq2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), dorg2r_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dorgr2_(integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dorm2r_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), dormr2_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal dscale;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), dlassq_(integer *
	    , doublereal *, integer *, doublereal *, doublereal *);
    static logical dtrong;
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

/*     To swap adjacent diagonal blocks A11*B11 and A22*B22 of size */
/*     1-by-1 or 2-by-2 in an upper (quasi) triangular matrix product */
/*     A*B by an orthogonal equivalence transformation. */

/*     (A, B) must be in periodic real Schur canonical form (as returned */
/*     by SLICOT Library routine MB03XP), i.e., A is block upper */
/*     triangular with 1-by-1 and 2-by-2 diagonal blocks, and B is upper */
/*     triangular. */

/*     Optionally, the matrices Q and Z of generalized Schur vectors are */
/*     updated. */

/*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)', */
/*         Z(in) * B(in) * Q(in)' = Z(out) * B(out) * Q(out)'. */

/*     This routine is largely based on the LAPACK routine DTGEX2 */
/*     developed by Bo Kagstrom and Peter Poromaa. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WANTQ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = .TRUE. :  The matrix Q is updated; */
/*             = .FALSE.:  the matrix Q is not required. */

/*     WANTZ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = .TRUE. :  The matrix Z is updated; */
/*             = .FALSE.:  the matrix Z is not required. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The order of the first block A11*B11. N1 = 0, 1 or 2. */

/*     N2      (input) INTEGER */
/*             The order of the second block A22*B22. N2 = 0, 1 or 2. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,N1+N2) */
/*             On entry, the leading (N1+N2)-by-(N1+N2) part of this */
/*             array must contain the matrix A. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the matrix A of the reordered pair. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. LDA >= MAX(1,N1+N2). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,N1+N2) */
/*             On entry, the leading (N1+N2)-by-(N1+N2) part of this */
/*             array must contain the matrix B. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the matrix B of the reordered pair. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. LDB >= MAX(1,N1+N2). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQ,N1+N2) */
/*             On entry, if WANTQ = .TRUE., the leading */
/*             (N1+N2)-by-(N1+N2) part of this array must contain the */
/*             orthogonal matrix Q. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the updated matrix Q. Q will be a rotation */
/*             matrix for N1=N2=1. */
/*             This array is not referenced if WANTQ = .FALSE.. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. LDQ >= 1. */
/*             If WANTQ = .TRUE., LDQ >= N1+N2. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDZ,N1+N2) */
/*             On entry, if WANTZ = .TRUE., the leading */
/*             (N1+N2)-by-(N1+N2) part of this array must contain the */
/*             orthogonal matrix Z. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the updated matrix Z. Z will be a rotation */
/*             matrix for N1=N2=1. */
/*             This array is not referenced if WANTZ = .FALSE.. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z. LDZ >= 1. */
/*             If WANTZ = .TRUE., LDZ >= N1+N2. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  the transformed matrix (A, B) would be */
/*                   too far from periodic Schur form; the blocks are */
/*                   not swapped and (A,B) and (Q,Z) are unchanged. */

/*     METHOD */

/*     In the current code both weak and strong stability tests are */
/*     performed. The user can omit the strong stability test by changing */
/*     the internal logical parameter WANDS to .FALSE.. See ref. [2] for */
/*     details. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A direct method for reordering eigenvalues in the generalized */
/*         real Schur form of a regular matrix pair (A,B), in M.S. Moonen */
/*         et al (eds.), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., 1993, pp. 195-218. */

/*     [2] Kagstrom, B., and Poromaa, P. */
/*         Computing eigenspaces with specified eigenvalues of a regular */
/*         matrix pair (A, B) and condition estimation: Theory, */
/*         algorithms and software, Numer. Algorithms, 1996, vol. 12, */
/*         pp. 369-407. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DTGPX2). */

/*     KEYWORDS */

/*     Eigenvalue, periodic Schur form, reordering */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;

/*     Quick return if possible. */
/*     For efficiency, the arguments are not checked. */

    if (*n1 <= 0 || *n2 <= 0) {
	return 0;
    }
    m = *n1 + *n2;

    weak = FALSE_;
    dtrong = FALSE_;

/*     Make a local copy of selected block. */

    dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, li, &c__4, (ftnlen)3);
    dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, ir, &c__4, (ftnlen)3);
    dlacpy_("Full", &m, &m, &a[a_offset], lda, s, &c__4, (ftnlen)4);
    dlacpy_("Full", &m, &m, &b[b_offset], ldb, t, &c__4, (ftnlen)4);

/*     Compute threshold for testing acceptance of swapping. */

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1) / eps;
    dscale = 0.;
    dsum = 1.;
    dlacpy_("Full", &m, &m, s, &c__4, dwork, &m, (ftnlen)4);
    i__1 = m * m;
    dlassq_(&i__1, dwork, &c__1, &dscale, &dsum);
    dlacpy_("Full", &m, &m, t, &c__4, dwork, &m, (ftnlen)4);
    i__1 = m * m;
    dlassq_(&i__1, dwork, &c__1, &dscale, &dsum);
    dnorm = dscale * sqrt(dsum);
/* Computing MAX */
    d__1 = eps * 10. * dnorm;
    thresh = max(d__1,smlnum);

    if (m == 2) {

/*        CASE 1: Swap 1-by-1 and 1-by-1 blocks. */

/*        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks */
/*        using Givens rotations and perform the swap tentatively. */

	f = s[5] * t[5] - t[0] * s[0];
	g = -s[5] * t[4] - t[0] * s[4];
	sb = abs(t[0]);
	sa = abs(s[5]);
	dlartg_(&f, &g, &ir[4], ir, &ddum);
	ir[1] = -ir[4];
	ir[5] = ir[0];
	drot_(&c__2, s, &c__1, &s[4], &c__1, ir, &ir[1]);
	drot_(&c__2, t, &c__4, &t[1], &c__4, ir, &ir[1]);
	if (sa >= sb) {
	    dlartg_(s, &s[1], li, &li[1], &ddum);
	} else {
	    dlartg_(&t[5], &t[1], li, &li[1], &ddum);
	    li[1] = -li[1];
	}
	drot_(&c__2, s, &c__4, &s[1], &c__4, li, &li[1]);
	drot_(&c__2, t, &c__1, &t[4], &c__1, li, &li[1]);
	li[5] = li[0];
	li[4] = -li[1];

/*        Weak stability test: */
/*           |S21| + |T21| <= O(EPS * F-norm((S, T))). */

	ws = abs(s[1]) + abs(t[1]);
	weak = ws <= thresh;
	if (! weak) {
	    goto L50;
	}

	if (TRUE_) {

/*           Strong stability test: */
/*             F-norm((A-QL'*S*QR, B-QR'*T*QL)) <= O(EPS*F-norm((A,B))). */

	    dlacpy_("Full", &m, &m, &a[a_offset], lda, &dwork[m * m], &m, (
		    ftnlen)4);
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &
		    c__4, s, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     ir, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
	    dscale = 0.;
	    dsum = 1.;
	    i__1 = m * m;
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);

	    dlacpy_("Full", &m, &m, &b[b_offset], ldb, &dwork[m * m], &m, (
		    ftnlen)4);
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &
		    c__4, t, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     li, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
	    i__1 = m * m;
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);
	    ss = dscale * sqrt(dsum);
	    dtrong = ss <= thresh;
	    if (! dtrong) {
		goto L50;
	    }
	}

/*        Update A and B. */

	dlacpy_("All", &m, &m, s, &c__4, &a[a_offset], lda, (ftnlen)3);
	dlacpy_("All", &m, &m, t, &c__4, &b[b_offset], ldb, (ftnlen)3);

/*        Set  N1-by-N2 (2,1) - blocks to ZERO. */

	a[a_dim1 + 2] = 0.;
	b[b_dim1 + 2] = 0.;

/*        Accumulate transformations into Q and Z if requested. */

	if (*wantq) {
	    drot_(&c__2, &q[q_dim1 + 1], &c__1, &q[(q_dim1 << 1) + 1], &c__1, 
		    li, &li[1]);
	}
	if (*wantz) {
	    drot_(&c__2, &z__[z_dim1 + 1], &c__1, &z__[(z_dim1 << 1) + 1], &
		    c__1, ir, &ir[1]);
	}

/*        Exit with INFO = 0 if swap was successfully performed. */

	return 0;

    } else {

/*        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2 */
/*                and 2-by-2 blocks. */

/*        Solve the periodic Sylvester equation */
/*                 S11 * R - L * S22 = SCALE * S12 */
/*                 T11 * L - R * T22 = SCALE * T12 */
/*        for R and L. Solutions in IR and LI. */

	dlacpy_("Full", n1, n2, &t[(*n1 + 1 << 2) - 4], &c__4, li, &c__4, (
		ftnlen)4);
	dlacpy_("Full", n1, n2, &s[(*n1 + 1 << 2) - 4], &c__4, &ir[*n2 + 1 + (
		*n1 + 1 << 2) - 5], &c__4, (ftnlen)4);
	sb04ow_(n1, n2, s, &c__4, &s[*n1 + 1 + (*n1 + 1 << 2) - 5], &c__4, &
		ir[*n2 + 1 + (*n1 + 1 << 2) - 5], &c__4, t, &c__4, &t[*n1 + 1 
		+ (*n1 + 1 << 2) - 5], &c__4, li, &c__4, &scale, iwork, &
		linfo);
	if (linfo != 0) {
	    goto L50;
	}

/*        Compute orthogonal matrix QL: */

/*                    QL' * LI = [ TL ] */
/*                               [ 0  ] */
/*        where */
/*                    LI =  [      -L              ]. */
/*                          [ SCALE * identity(N2) ] */

	i__1 = *n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dscal_(n1, &c_b48, &li[(i__ << 2) - 4], &c__1);
	    li[*n1 + i__ + (i__ << 2) - 5] = scale;
/* L10: */
	}
	dgeqr2_(&m, n2, li, &c__4, taul, dwork, &linfo);
	dorg2r_(&m, &m, n2, li, &c__4, taul, dwork, &linfo);

/*        Compute orthogonal matrix RQ: */

/*                    IR * RQ' =   [ 0  TR], */

/*         where IR = [ SCALE * identity(N1), R ]. */

	i__1 = *n1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ir[*n2 + i__ + (i__ << 2) - 5] = scale;
/* L20: */
	}
	dgerq2_(n1, &m, &ir[*n2], &c__4, taur, dwork, &linfo);
	dorgr2_(&m, &m, n1, ir, &c__4, taur, dwork, &linfo);

/*        Perform the swapping tentatively: */

	dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &c__4, s, 
		&c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
	dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b42, dwork, &m, ir,
		 &c__4, &c_b5, s, &c__4, (ftnlen)12, (ftnlen)9);
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &c__4, 
		t, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, dwork, &m, 
		li, &c__4, &c_b5, t, &c__4, (ftnlen)12, (ftnlen)12);
	dlacpy_("All", &m, &m, s, &c__4, scpy, &c__4, (ftnlen)3);
	dlacpy_("All", &m, &m, t, &c__4, tcpy, &c__4, (ftnlen)3);
	dlacpy_("All", &m, &m, ir, &c__4, ircop, &c__4, (ftnlen)3);
	dlacpy_("All", &m, &m, li, &c__4, licop, &c__4, (ftnlen)3);

/*        Triangularize the B-part by a QR factorization. */
/*        Apply transformation (from left) to A-part, giving S. */

	dgeqr2_(&m, &m, t, &c__4, taur, dwork, &linfo);
	dorm2r_("Right", "No Transpose", &m, &m, &m, t, &c__4, taur, s, &c__4,
		 dwork, &linfo, (ftnlen)5, (ftnlen)12);
	dorm2r_("Left", "Transpose", &m, &m, &m, t, &c__4, taur, ir, &c__4, 
		dwork, &linfo, (ftnlen)4, (ftnlen)9);

/*        Compute F-norm(S21) in BRQA21. (T21 is 0.) */

	dscale = 0.;
	dsum = 1.;
	i__1 = *n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dlassq_(n1, &s[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &dsum);
/* L30: */
	}
	brqa21 = dscale * sqrt(dsum);

/*        Triangularize the B-part by an RQ factorization. */
/*        Apply transformation (from right) to A-part, giving S. */

	dgerq2_(&m, &m, tcpy, &c__4, taul, dwork, &linfo);
	dormr2_("Left", "No Transpose", &m, &m, &m, tcpy, &c__4, taul, scpy, &
		c__4, dwork, &linfo, (ftnlen)4, (ftnlen)12);
	dormr2_("Right", "Transpose", &m, &m, &m, tcpy, &c__4, taul, licop, &
		c__4, dwork, &linfo, (ftnlen)5, (ftnlen)9);

/*        Compute F-norm(S21) in BQRA21. (T21 is 0.) */

	dscale = 0.;
	dsum = 1.;
	i__1 = *n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dlassq_(n1, &scpy[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &
		    dsum);
/* L40: */
	}
	bqra21 = dscale * sqrt(dsum);

/*        Decide which method to use. */
/*          Weak stability test: */
/*             F-norm(S21) <= O(EPS * F-norm((S, T))) */

	if (bqra21 <= brqa21 && bqra21 <= thresh) {
	    dlacpy_("All", &m, &m, scpy, &c__4, s, &c__4, (ftnlen)3);
	    dlacpy_("All", &m, &m, tcpy, &c__4, t, &c__4, (ftnlen)3);
	    dlacpy_("All", &m, &m, ircop, &c__4, ir, &c__4, (ftnlen)3);
	    dlacpy_("All", &m, &m, licop, &c__4, li, &c__4, (ftnlen)3);
	} else if (brqa21 >= thresh) {
	    goto L50;
	}

/*        Set lower triangle of B-part to zero */

	i__1 = m - 1;
	i__2 = m - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b5, &c_b5, &t[1], &c__4, (ftnlen)5);

	if (TRUE_) {

/*           Strong stability test: */
/*              F-norm((A-QL*S*QR', B-QR*T*QL')) <= O(EPS*F-norm((A,B))) */

	    dlacpy_("All", &m, &m, &a[a_offset], lda, &dwork[m * m], &m, (
		    ftnlen)3);
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &
		    c__4, s, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b48, dwork, 
		    &m, ir, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)12);
	    dscale = 0.;
	    dsum = 1.;
	    i__1 = m * m;
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);

	    dlacpy_("All", &m, &m, &b[b_offset], ldb, &dwork[m * m], &m, (
		    ftnlen)3);
	    dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &c__4,
		     t, &c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     li, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
	    i__1 = m * m;
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);
	    ss = dscale * sqrt(dsum);
	    dtrong = ss <= thresh;
	    if (! dtrong) {
		goto L50;
	    }

	}

/*        If the swap is accepted ("weakly" and "strongly"), apply the */
/*        transformations and set N1-by-N2 (2,1)-block to zero. */

	dlaset_("All", n1, n2, &c_b5, &c_b5, &s[*n2], &c__4, (ftnlen)3);

/*        Copy (S,T) to (A,B). */

	dlacpy_("All", &m, &m, s, &c__4, &a[a_offset], lda, (ftnlen)3);
	dlacpy_("All", &m, &m, t, &c__4, &b[b_offset], ldb, (ftnlen)3);
	dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, t, &c__4, (ftnlen)3);

/*        Standardize existing 2-by-2 blocks. */

	dlaset_("All", &m, &m, &c_b5, &c_b5, dwork, &m, (ftnlen)3);
	dwork[0] = 1.;
	t[0] = 1.;
	if (*n2 > 1) {
	    mb03yt_(&a[a_offset], lda, &b[b_offset], ldb, ar, ai, be, dwork, &
		    dwork[1], t, &t[1]);
	    dwork[m] = -dwork[1];
	    dwork[m + 1] = dwork[0];
	    t[*n2 + (*n2 << 2) - 5] = t[0];
	    t[4] = -t[1];
	}
	dwork[m * m - 1] = 1.;
	t[m + (m << 2) - 5] = 1.;

	if (*n1 > 1) {
	    mb03yt_(&a[*n2 + 1 + (*n2 + 1) * a_dim1], lda, &b[*n2 + 1 + (*n2 
		    + 1) * b_dim1], ldb, taur, taul, &dwork[m * m], &dwork[*
		    n2 * m + *n2], &dwork[*n2 * m + *n2 + 1], &t[*n2 + 1 + (*
		    n2 + 1 << 2) - 5], &t[m + (m - 1 << 2) - 5]);
	    dwork[m * m - 1] = dwork[*n2 * m + *n2];
	    dwork[m * m - 2] = -dwork[*n2 * m + *n2 + 1];
	    t[m + (m << 2) - 5] = t[*n2 + 1 + (*n2 + 1 << 2) - 5];
	    t[m - 1 + (m << 2) - 5] = -t[m + (m - 1 << 2) - 5];
	}

	dgemm_("Transpose", "No Transpose", n2, n1, n2, &c_b42, dwork, &m, &a[
		(*n2 + 1) * a_dim1 + 1], lda, &c_b5, &dwork[m * m], n2, (
		ftnlen)9, (ftnlen)12);
	dlacpy_("All", n2, n1, &dwork[m * m], n2, &a[(*n2 + 1) * a_dim1 + 1], 
		lda, (ftnlen)3);
	dgemm_("Transpose", "No Transpose", n2, n1, n2, &c_b42, t, &c__4, &b[(
		*n2 + 1) * b_dim1 + 1], ldb, &c_b5, &dwork[m * m], n2, (
		ftnlen)9, (ftnlen)12);
	dlacpy_("All", n2, n1, &dwork[m * m], n2, &b[(*n2 + 1) * b_dim1 + 1], 
		ldb, (ftnlen)3);
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &c__4, 
		dwork, &m, &c_b5, &dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
	dlacpy_("All", &m, &m, &dwork[m * m], &m, li, &c__4, (ftnlen)3);
	dgemm_("No Transpose", "No Transpose", n2, n1, n1, &c_b42, &a[(*n2 + 
		1) * a_dim1 + 1], lda, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &
		c__4, &c_b5, &dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
	dlacpy_("All", n2, n1, &dwork[m * m], &m, &a[(*n2 + 1) * a_dim1 + 1], 
		lda, (ftnlen)3);
	dgemm_("No Transpose", "No Transpose", n2, n1, n1, &c_b42, &b[(*n2 + 
		1) * b_dim1 + 1], ldb, &dwork[*n2 * m + *n2], &m, &c_b5, &
		dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
	dlacpy_("All", n2, n1, &dwork[m * m], &m, &b[(*n2 + 1) * b_dim1 + 1], 
		ldb, (ftnlen)3);
	dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, t, &c__4, ir, 
		&c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
	dlacpy_("All", &m, &m, dwork, &m, ir, &c__4, (ftnlen)3);

/*        Accumulate transformations into Q and Z if requested. */

	if (*wantq) {
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, &q[
		    q_offset], ldq, li, &c__4, &c_b5, dwork, &m, (ftnlen)12, (
		    ftnlen)12);
	    dlacpy_("All", &m, &m, dwork, &m, &q[q_offset], ldq, (ftnlen)3);
	}

	if (*wantz) {
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b42, &z__[
		    z_offset], ldz, ir, &c__4, &c_b5, dwork, &m, (ftnlen)12, (
		    ftnlen)9);
	    dlacpy_("Full", &m, &m, dwork, &m, &z__[z_offset], ldz, (ftnlen)4)
		    ;

	}

/*        Exit with INFO = 0 if swap was successfully performed. */

	return 0;

    }

/*     Exit with INFO = 1 if swap was rejected. */

L50:

    *info = 1;
    return 0;
/* *** Last line of MB03WA *** */
} /* mb03wa_ */

