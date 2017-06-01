/* SB04OW.f -- translated by f2c (version 20100827).
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

static integer c__8 = 8;
static integer c__1 = 1;
static doublereal c_b23 = -1.;
static doublereal c_b37 = 1.;
static doublereal c_b51 = 0.;

/* Subroutine */ int sb04ow_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *e, integer *lde, 
	doublereal *f, integer *ldf, doublereal *scale, integer *iwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, p, q;
    static doublereal z__[64]	/* was [8][8] */;
    static integer ie, je, mb, nb, ii, jj, is, js;
    static doublereal rhs[8];
    static integer isp1, jsp1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr, zdim, ipiv[8], jpiv[8];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dgemv_(
	    char *, integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dgesc2_(integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *), dgetc2_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To solve a periodic Sylvester equation */

/*              A * R - L * B = scale * C                           (1) */
/*              D * L - R * E = scale * F, */

/*     using Level 1 and 2 BLAS, where R and L are unknown M-by-N */
/*     matrices, (A, D), (B, E) and (C, F) are given matrix pairs of */
/*     size M-by-M, N-by-N and M-by-N, respectively, with real entries. */
/*     (A, D) and (B, E) must be in periodic Schur form, i.e. A, B are */
/*     upper quasi triangular and D, E are upper triangular. The solution */
/*     (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output scaling */
/*     factor chosen to avoid overflow. */

/*     This routine is largely based on the LAPACK routine DTGSY2 */
/*     developed by Bo Kagstrom and Peter Poromaa. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of A and D, and the row dimension of C, F, R */
/*             and L.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of B and E, and the column dimension of C, F, R */
/*             and L.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the upper quasi triangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi triangular matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand-side of the first matrix equation */
/*             in (1). */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution R. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,M). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the upper triangular matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,M). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand-side of the second matrix equation */
/*             in (1). */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution L. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the arrays */
/*             C and F will hold the solutions R and L, respectively, to */
/*             a slightly perturbed system but the input matrices A, B, D */
/*             and E have not been changed. If SCALE = 0, C and F will */
/*             hold solutions to the homogeneous system with C = F = 0. */
/*             Normally, SCALE = 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+N+2) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  the matrix products A*D and B*E have common or very */
/*                   close eigenvalues. */

/*     METHOD */

/*     In matrix notation solving equation (1) corresponds to solving */
/*     Z*x = scale*b, where Z is defined as */

/*         Z = [  kron(In, A)  -kron(B', Im) ]            (2) */
/*             [ -kron(E', Im)  kron(In, D)  ], */

/*     Ik is the identity matrix of size k and X' is the transpose of X. */
/*     kron(X, Y) is the Kronecker product between the matrices X and Y. */
/*     In the process of solving (1), we solve a number of such systems */
/*     where Dim(Im), Dim(In) = 1 or 2. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A Direct Method for Reordering Eigenvalues in the Generalized */
/*         Real Schur Form of a Regular Matrix Pair (A,B). M.S. Moonen */
/*         et al (eds.), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., pp. 195-218, 1993. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         A Schur approach for solving some periodic matrix equations. */
/*         U. Helmke et al (eds.), Systems and Networks: Mathematical */
/*         Theory and Applications, Akademie Verlag, Berlin, vol. 77, */
/*         pp. 339-362, 1994. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DTGPY2). */

/*     KEYWORDS */

/*     Matrix equation, periodic Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --iwork;

    /* Function Body */
    *info = 0;
    ierr = 0;
    if (*m <= 0) {
	*info = -1;
    } else if (*n <= 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -6;
    } else if (*ldc < max(1,*m)) {
	*info = -8;
    } else if (*ldd < max(1,*m)) {
	*info = -10;
    } else if (*lde < max(1,*n)) {
	*info = -12;
    } else if (*ldf < max(1,*m)) {
	*info = -14;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB04OW", &i__1, (ftnlen)6);
	return 0;
    }

/*     Determine block structure of A. */

    p = 0;
    i__ = 1;
L10:
    if (i__ > *m) {
	goto L20;
    }
    ++p;
    iwork[p] = i__;
    if (i__ == *m) {
	goto L20;
    }
    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
	i__ += 2;
    } else {
	++i__;
    }
    goto L10;
L20:
    iwork[p + 1] = *m + 1;

/*     Determine block structure of B. */

    q = p + 1;
    j = 1;
L30:
    if (j > *n) {
	goto L40;
    }
    ++q;
    iwork[q] = j;
    if (j == *n) {
	goto L40;
    }
    if (b[j + 1 + j * b_dim1] != 0.) {
	j += 2;
    } else {
	++j;
    }
    goto L30;
L40:
    iwork[q + 1] = *n + 1;

/*     Solve (I, J) - subsystem */
/*       A(I,I) * R(I,J) - L(I,J) * B(J,J) = C(I,J) */
/*       D(I,I) * L(I,J) - R(I,J) * E(J,J) = F(I,J) */
/*     for I = P, P - 1, ..., 1; J = 1, 2, ..., Q. */

    *scale = 1.;
    scaloc = 1.;
    i__1 = q;
    for (j = p + 2; j <= i__1; ++j) {
	js = iwork[j];
	jsp1 = js + 1;
	je = iwork[j + 1] - 1;
	nb = je - js + 1;
	for (i__ = p; i__ >= 1; --i__) {

	    is = iwork[i__];
	    isp1 = is + 1;
	    ie = iwork[i__ + 1] - 1;
	    mb = ie - is + 1;
	    zdim = mb * nb << 1;

	    if (mb == 1 && nb == 1) {

/*              Build a 2-by-2 system Z * x = RHS. */

		z__[0] = a[is + is * a_dim1];
		z__[1] = -e[js + js * e_dim1];
		z__[8] = -b[js + js * b_dim1];
		z__[9] = d__[is + is * d_dim1];

/*              Set up right hand side(s). */

		rhs[0] = c__[is + js * c_dim1];
		rhs[1] = f[is + js * f_dim1];

/*              Solve Z * x = RHS. */

		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}

		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
		if (scaloc != 1.) {
		    i__2 = *n;
		    for (k = 1; k <= i__2; ++k) {
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
/* L50: */
		    }
		    *scale *= scaloc;
		}

/*              Unpack solution vector(s). */

		c__[is + js * c_dim1] = rhs[0];
		f[is + js * f_dim1] = rhs[1];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

		if (i__ > 1) {
		    i__2 = is - 1;
		    d__1 = -rhs[0];
		    daxpy_(&i__2, &d__1, &a[is * a_dim1 + 1], &c__1, &c__[js *
			     c_dim1 + 1], &c__1);
		    i__2 = is - 1;
		    d__1 = -rhs[1];
		    daxpy_(&i__2, &d__1, &d__[is * d_dim1 + 1], &c__1, &f[js *
			     f_dim1 + 1], &c__1);
		}
		if (j < q) {
		    i__2 = *n - je;
		    daxpy_(&i__2, &rhs[1], &b[js + (je + 1) * b_dim1], ldb, &
			    c__[is + (je + 1) * c_dim1], ldc);
		    i__2 = *n - je;
		    daxpy_(&i__2, rhs, &e[js + (je + 1) * e_dim1], lde, &f[is 
			    + (je + 1) * f_dim1], ldf);
		}

	    } else if (mb == 1 && nb == 2) {

/*              Build a 4-by-4 system Z * x = RHS. */

		z__[0] = a[is + is * a_dim1];
		z__[1] = 0.;
		z__[2] = -e[js + js * e_dim1];
		z__[3] = -e[js + jsp1 * e_dim1];

		z__[8] = 0.;
		z__[9] = a[is + is * a_dim1];
		z__[10] = 0.;
		z__[11] = -e[jsp1 + jsp1 * e_dim1];

		z__[16] = -b[js + js * b_dim1];
		z__[17] = -b[js + jsp1 * b_dim1];
		z__[18] = d__[is + is * d_dim1];
		z__[19] = 0.;

		z__[24] = -b[jsp1 + js * b_dim1];
		z__[25] = -b[jsp1 + jsp1 * b_dim1];
		z__[26] = 0.;
		z__[27] = d__[is + is * d_dim1];

/*              Set up right hand side(s). */

		rhs[0] = c__[is + js * c_dim1];
		rhs[1] = c__[is + jsp1 * c_dim1];
		rhs[2] = f[is + js * f_dim1];
		rhs[3] = f[is + jsp1 * f_dim1];

/*              Solve Z * x = RHS. */

		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}

		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
		if (scaloc != 1.) {
		    i__2 = *n;
		    for (k = 1; k <= i__2; ++k) {
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
/* L60: */
		    }
		    *scale *= scaloc;
		}

/*              Unpack solution vector(s). */

		c__[is + js * c_dim1] = rhs[0];
		c__[is + jsp1 * c_dim1] = rhs[1];
		f[is + js * f_dim1] = rhs[2];
		f[is + jsp1 * f_dim1] = rhs[3];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

		if (i__ > 1) {
		    i__2 = is - 1;
		    dger_(&i__2, &nb, &c_b23, &a[is * a_dim1 + 1], &c__1, rhs,
			     &c__1, &c__[js * c_dim1 + 1], ldc);
		    i__2 = is - 1;
		    dger_(&i__2, &nb, &c_b23, &d__[is * d_dim1 + 1], &c__1, &
			    rhs[2], &c__1, &f[js * f_dim1 + 1], ldf);
		}
		if (j < q) {
		    i__2 = *n - je;
		    daxpy_(&i__2, &rhs[2], &b[js + (je + 1) * b_dim1], ldb, &
			    c__[is + (je + 1) * c_dim1], ldc);
		    i__2 = *n - je;
		    daxpy_(&i__2, rhs, &e[js + (je + 1) * e_dim1], lde, &f[is 
			    + (je + 1) * f_dim1], ldf);
		    i__2 = *n - je;
		    daxpy_(&i__2, &rhs[3], &b[jsp1 + (je + 1) * b_dim1], ldb, 
			    &c__[is + (je + 1) * c_dim1], ldc);
		    i__2 = *n - je;
		    daxpy_(&i__2, &rhs[1], &e[jsp1 + (je + 1) * e_dim1], lde, 
			    &f[is + (je + 1) * f_dim1], ldf);
		}

	    } else if (mb == 2 && nb == 1) {

/*              Build a 4-by-4 system Z * x = RHS. */

		z__[0] = a[is + is * a_dim1];
		z__[1] = a[isp1 + is * a_dim1];
		z__[2] = -e[js + js * e_dim1];
		z__[3] = 0.;

		z__[8] = a[is + isp1 * a_dim1];
		z__[9] = a[isp1 + isp1 * a_dim1];
		z__[10] = 0.;
		z__[11] = -e[js + js * e_dim1];

		z__[16] = -b[js + js * b_dim1];
		z__[17] = 0.;
		z__[18] = d__[is + is * d_dim1];
		z__[19] = 0.;

		z__[24] = 0.;
		z__[25] = -b[js + js * b_dim1];
		z__[26] = d__[is + isp1 * d_dim1];
		z__[27] = d__[isp1 + isp1 * d_dim1];

/*              Set up right hand side(s). */

		rhs[0] = c__[is + js * c_dim1];
		rhs[1] = c__[isp1 + js * c_dim1];
		rhs[2] = f[is + js * f_dim1];
		rhs[3] = f[isp1 + js * f_dim1];

/*              Solve Z * x = RHS. */

		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}

		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
		if (scaloc != 1.) {
		    i__2 = *n;
		    for (k = 1; k <= i__2; ++k) {
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
/* L70: */
		    }
		    *scale *= scaloc;
		}

/*              Unpack solution vector(s). */

		c__[is + js * c_dim1] = rhs[0];
		c__[isp1 + js * c_dim1] = rhs[1];
		f[is + js * f_dim1] = rhs[2];
		f[isp1 + js * f_dim1] = rhs[3];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

		if (i__ > 1) {
		    i__2 = is - 1;
		    dgemv_("N", &i__2, &mb, &c_b23, &a[is * a_dim1 + 1], lda, 
			    rhs, &c__1, &c_b37, &c__[js * c_dim1 + 1], &c__1, 
			    (ftnlen)1);
		    i__2 = is - 1;
		    dgemv_("N", &i__2, &mb, &c_b23, &d__[is * d_dim1 + 1], 
			    ldd, &rhs[2], &c__1, &c_b37, &f[js * f_dim1 + 1], 
			    &c__1, (ftnlen)1);
		}
		if (j < q) {
		    i__2 = *n - je;
		    dger_(&mb, &i__2, &c_b37, &rhs[2], &c__1, &b[js + (je + 1)
			     * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], 
			    ldc);
		    i__2 = *n - je;
		    dger_(&mb, &i__2, &c_b37, rhs, &c__1, &e[js + (je + 1) * 
			    e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
		}

	    } else if (mb == 2 && nb == 2) {

/*              Build an 8-by-8 system Z * x = RHS. */

		dlaset_("All", &c__8, &c__8, &c_b51, &c_b51, z__, &c__8, (
			ftnlen)3);

		z__[0] = a[is + is * a_dim1];
		z__[1] = a[isp1 + is * a_dim1];
		z__[4] = -e[js + js * e_dim1];
		z__[6] = -e[js + jsp1 * e_dim1];

		z__[8] = a[is + isp1 * a_dim1];
		z__[9] = a[isp1 + isp1 * a_dim1];
		z__[13] = -e[js + js * e_dim1];
		z__[15] = -e[js + jsp1 * e_dim1];

		z__[18] = a[is + is * a_dim1];
		z__[19] = a[isp1 + is * a_dim1];
		z__[22] = -e[jsp1 + jsp1 * e_dim1];

		z__[26] = a[is + isp1 * a_dim1];
		z__[27] = a[isp1 + isp1 * a_dim1];
		z__[31] = -e[jsp1 + jsp1 * e_dim1];

		z__[32] = -b[js + js * b_dim1];
		z__[34] = -b[js + jsp1 * b_dim1];
		z__[36] = d__[is + is * d_dim1];

		z__[41] = -b[js + js * b_dim1];
		z__[43] = -b[js + jsp1 * b_dim1];
		z__[44] = d__[is + isp1 * d_dim1];
		z__[45] = d__[isp1 + isp1 * d_dim1];

		z__[48] = -b[jsp1 + js * b_dim1];
		z__[50] = -b[jsp1 + jsp1 * b_dim1];
		z__[54] = d__[is + is * d_dim1];

		z__[57] = -b[jsp1 + js * b_dim1];
		z__[59] = -b[jsp1 + jsp1 * b_dim1];

		z__[62] = d__[is + isp1 * d_dim1];
		z__[63] = d__[isp1 + isp1 * d_dim1];

/*              Set up right hand side(s). */

		k = 1;
		ii = mb * nb + 1;
		i__2 = nb - 1;
		for (jj = 0; jj <= i__2; ++jj) {
		    dcopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &rhs[k 
			    - 1], &c__1);
		    dcopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[ii - 
			    1], &c__1);
		    k += mb;
		    ii += mb;
/* L80: */
		}

/*              Solve Z * x = RHS. */

		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}

		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
		if (scaloc != 1.) {
		    i__2 = *n;
		    for (k = 1; k <= i__2; ++k) {
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
/* L90: */
		    }
		    *scale *= scaloc;
		}

/*              Unpack solution vector(s). */

		k = 1;
		ii = mb * nb + 1;
		i__2 = nb - 1;
		for (jj = 0; jj <= i__2; ++jj) {
		    dcopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
			    c_dim1], &c__1);
		    dcopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
			    f_dim1], &c__1);
		    k += mb;
		    ii += mb;
/* L100: */
		}

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

		k = mb * nb + 1;
		if (i__ > 1) {
		    i__2 = is - 1;
		    dgemm_("N", "N", &i__2, &nb, &mb, &c_b23, &a[is * a_dim1 
			    + 1], lda, rhs, &mb, &c_b37, &c__[js * c_dim1 + 1]
			    , ldc, (ftnlen)1, (ftnlen)1);
		    i__2 = is - 1;
		    dgemm_("N", "N", &i__2, &nb, &mb, &c_b23, &d__[is * 
			    d_dim1 + 1], ldd, &rhs[k - 1], &mb, &c_b37, &f[js 
			    * f_dim1 + 1], ldf, (ftnlen)1, (ftnlen)1);
		}
		if (j < q) {
		    i__2 = *n - je;
		    dgemm_("N", "N", &mb, &i__2, &nb, &c_b37, &rhs[k - 1], &
			    mb, &b[js + (je + 1) * b_dim1], ldb, &c_b37, &c__[
			    is + (je + 1) * c_dim1], ldc, (ftnlen)1, (ftnlen)
			    1);
		    i__2 = *n - je;
		    dgemm_("N", "N", &mb, &i__2, &nb, &c_b37, rhs, &mb, &e[js 
			    + (je + 1) * e_dim1], lde, &c_b37, &f[is + (je + 
			    1) * f_dim1], ldf, (ftnlen)1, (ftnlen)1);
		}

	    }

/* L110: */
	}
/* L120: */
    }
    return 0;
/* *** Last line of SB04OW *** */
} /* sb04ow_ */

