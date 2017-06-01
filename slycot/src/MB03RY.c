/* MB03RY.f -- translated by f2c (version 20100827).
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

static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c_n1 = -1;

/* Subroutine */ int mb03ry_(integer *m, integer *n, doublereal *pmax, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal p[4];
    static integer dk, dl, kk, ll, kk1, lm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static doublereal pnorm;
    extern /* Subroutine */ int dlasy2_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


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

/*     To solve the Sylvester equation -AX + XB = C, where A and B are */
/*     M-by-M and N-by-N matrices, respectively, in real Schur form. */

/*     This routine is intended to be called only by SLICOT Library */
/*     routine MB03RD. For efficiency purposes, the computations are */
/*     aborted when the infinity norm of an elementary submatrix of X is */
/*     greater than a given value PMAX. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A and the number of rows of the */
/*             matrices C and X.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix B and the number of columns of the */
/*             matrices C and X.  N >= 0. */

/*     PMAX    (input) DOUBLE PRECISION */
/*             An upper bound for the infinity norm of an elementary */
/*             submatrix of X (see METHOD). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain the */
/*             matrix A of the Sylvester equation, in real Schur form. */
/*             The elements below the real Schur form are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix B of the Sylvester equation, in real Schur form. */
/*             The elements below the real Schur form are not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix C of the Sylvester equation. */
/*             On exit, if INFO = 0, the leading M-by-N part of this */
/*             array contains the solution matrix X of the Sylvester */
/*             equation, and each elementary submatrix of X (see METHOD) */
/*             has the infinity norm less than or equal to PMAX. */
/*             On exit, if INFO = 1, the solution matrix X has not been */
/*             computed completely, because an elementary submatrix of X */
/*             had the infinity norm greater than PMAX. Part of the */
/*             matrix C has possibly been overwritten with the */
/*             corresponding part of X. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  an elementary submatrix of X had the infinity norm */
/*                   greater than the given value PMAX. */

/*     METHOD */

/*     The routine uses an adaptation of the standard method for solving */
/*     Sylvester equations [1], which controls the magnitude of the */
/*     individual elements of the computed solution [2]. The equation */
/*     -AX + XB = C can be rewritten as */
/*                                 p            l-1 */
/*       -A  X   + X  B   = C   + sum  A  X   - sum  X  B */
/*         kk kl    kl ll    kl  i=k+1  ki il   j=1   kj jl */

/*     for l = 1:q, and k = p:-1:1, where A  , B  , C  , and X  , are */
/*                                         kk   ll   kl       kl */
/*     block submatrices defined by the partitioning induced by the Schur */
/*     form of A and B, and p and q are the numbers of the diagonal */
/*     blocks of A and B, respectively. So, the elementary submatrices of */
/*     X are found block column by block column, starting from the */
/*     bottom. If any such elementary submatrix has the infinity norm */
/*     greater than the given value PMAX, the calculations are ended. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Bavely, C. and Stewart, G.W. */
/*         An Algorithm for Computing Reducing Subspaces by Block */
/*         Diagonalization. */
/*         SIAM J. Numer. Anal., 16, pp. 359-367, 1979. */

/*     NUMERICAL ASPECTS */
/*                               2      2 */
/*     The algorithm requires 0(M N + MN ) operations. */

/*     FURTHER COMMENTS */

/*     Let */

/*            ( A   C )       ( I   X ) */
/*        M = (       ),  Y = (       ). */
/*            ( 0   B )       ( 0   I ) */

/*     Then */

/*         -1      ( A   0 ) */
/*        Y  M Y = (       ), */
/*                 ( 0   B ) */

/*     hence Y is an non-orthogonal transformation matrix which performs */
/*     the reduction of M to a block-diagonal form. Bounding a norm of */
/*     X is equivalent to setting an upper bound to the condition number */
/*     of the transformation matrix Y. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on the RASP routine SYLSM by A. Varga, German Aerospace */
/*     Center, DLR Oberpfaffenhofen. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Diagonalization, real Schur form, Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     For efficiency reasons, this routine does not check the input */
/*     parameters for errors. */

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

    /* Function Body */
    *info = 0;

/*     Column loop indexed by L. */

    l = 1;
/*     WHILE ( L.LE.N ) DO */
L10:
    if (l <= *n) {
	lm1 = l - 1;
	dl = 1;
	if (l < *n) {
	    if (b[l + 1 + l * b_dim1] != 0.) {
		dl = 2;
	    }
	}
	ll = lm1 + dl;
	if (lm1 > 0) {

/*           Update one (or two) column(s) of C. */

	    if (dl == 2) {
		dgemm_("No transpose", "No transpose", m, &dl, &lm1, &c_b5, &
			c__[c_offset], ldc, &b[l * b_dim1 + 1], ldb, &c_b6, &
			c__[l * c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
	    } else {
		dgemv_("No transpose", m, &lm1, &c_b5, &c__[c_offset], ldc, &
			b[l * b_dim1 + 1], &c__1, &c_b6, &c__[l * c_dim1 + 1],
			 &c__1, (ftnlen)12);
	    }
	}

/*        Row loop indexed by KK. */

	kk = *m;
/*        WHILE ( KK.GE.1 ) DO */
L20:
	if (kk >= 1) {
	    kk1 = kk + 1;
	    dk = 1;
	    if (kk > 1) {
		if (a[kk + (kk - 1) * a_dim1] != 0.) {
		    dk = 2;
		}
	    }
	    k = kk1 - dk;
	    if (k < *m) {

/*              Update an elementary submatrix of C. */

		i__1 = ll;
		for (j = l; j <= i__1; ++j) {

		    i__2 = kk;
		    for (i__ = k; i__ <= i__2; ++i__) {
			i__3 = *m - kk;
			c__[i__ + j * c_dim1] += ddot_(&i__3, &a[i__ + kk1 * 
				a_dim1], lda, &c__[kk1 + j * c_dim1], &c__1);
/* L30: */
		    }

/* L40: */
		}

	    }
	    dlasy2_(&c_false, &c_false, &c_n1, &dk, &dl, &a[k + k * a_dim1], 
		    lda, &b[l + l * b_dim1], ldb, &c__[k + l * c_dim1], ldc, &
		    scale, p, &dk, &pnorm, &ierr);
	    if (scale != 1. || pnorm > *pmax) {
		*info = 1;
		return 0;
	    }
	    c__[k + l * c_dim1] = -p[0];
	    if (dl == 1) {
		if (dk == 2) {
		    c__[kk + l * c_dim1] = -p[1];
		}
	    } else {
		if (dk == 1) {
		    c__[k + ll * c_dim1] = -p[1];
		} else {
		    c__[kk + l * c_dim1] = -p[1];
		    c__[k + ll * c_dim1] = -p[2];
		    c__[kk + ll * c_dim1] = -p[3];
		}
	    }
	    kk -= dk;
	    goto L20;
	}
/*        END WHILE 20 */
	l += dl;
	goto L10;
    }
/*     END WHILE 10 */
    return 0;
/* *** Last line of MB03RY *** */
} /* mb03ry_ */

