/* SG03BW.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;
static doublereal c_b18 = 1.;
static doublereal c_b19 = 0.;
static integer c__2 = 2;
static doublereal c_b23 = -1.;

/* Subroutine */ int sg03bw_(char *trans, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *c__, integer *ldc, doublereal *e, 
	integer *lde, doublereal *d__, integer *ldd, doublereal *x, integer *
	ldx, doublereal *scale, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, d_dim1, d_offset, e_dim1, 
	    e_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ma, mb, me;
    static doublereal tm[4]	/* was [2][2] */;
    static integer mai, maj;
    static doublereal mat[16]	/* was [4][4] */, rhs[4];
    static integer piv1[4], piv2[4], info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb02uu_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *), mb02uv_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *);
    static doublereal scale1;
    static integer dimmat;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notrns;


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

/*     To solve for X the generalized Sylvester equation */

/*         T            T */
/*        A  * X * C + E  * X * D  =  SCALE * Y,                      (1) */

/*     or the transposed equation */

/*                 T            T */
/*        A * X * C  + E * X * D   =  SCALE * Y,                      (2) */

/*     where A and E are real M-by-M matrices, C and D are real N-by-N */
/*     matrices, X and Y are real M-by-N matrices. N is either 1 or 2. */
/*     The pencil A - lambda * E must be in generalized real Schur form */
/*     (A upper quasitriangular, E upper triangular). SCALE is an output */
/*     scale factor, set to avoid overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrices A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrices C and D.  N = 1 or N = 2. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper quasitriangular matrix A. The elements below the */
/*             upper Hessenberg part are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper triangular matrix E. The elements below the main */
/*             diagonal are not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,M). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix Y. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution matrix X. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             0 < SCALE <= 1. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the generalized Sylvester equation is (nearly) */
/*                   singular to working precision;  perturbed values */
/*                   were used to solve the equation (but the matrices */
/*                   A, C, D, and E are unchanged). */

/*     METHOD */

/*     The method used by the routine is based on a generalization of the */
/*     algorithm due to Bartels and Stewart [1]. See also [2] and [3] for */
/*     details. */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Gardiner, J.D., Laub, A.J., Amato, J.J., Moler, C.B. */
/*         Solution of the Sylvester Matrix Equation */
/*         A X B**T + C X D**T = E. */
/*         A.C.M. Trans. Math. Soft., vol. 18, no. 2, pp. 223-231, 1992. */

/*     [3] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The routine requires about 2 * N * M**2 flops. Note that we count */
/*     a single floating point arithmetic operation as one flop. */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     FURTHER COMMENTS */

/*     When near singularity is detected, perturbed values are used */
/*     to solve the equation (but the given matrices are unchanged). */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     Decode input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n != 1 && *n != 2) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*ldc < max(1,*n)) {
	*info = -7;
    } else if (*lde < max(1,*m)) {
	*info = -9;
    } else if (*ldd < max(1,*n)) {
	*info = -11;
    } else if (*ldx < max(1,*m)) {
	*info = -13;
    } else {
	*info = 0;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SG03BW", &i__1, (ftnlen)6);
	return 0;
    }

    *scale = 1.;

/*     Quick return if possible. */

    if (*m == 0) {
	return 0;
    }

    if (notrns) {

/*        Solve equation (1). */

/*        Compute block row X(MA:ME,:). MB denotes the number of rows in */
/*        this block row. */

	me = 0;
/*        WHILE ( ME .NE. M ) DO */
L20:
	if (me != *m) {
	    ma = me + 1;
	    if (ma == *m) {
		me = *m;
		mb = 1;
	    } else {
		if (a[ma + 1 + ma * a_dim1] == 0.) {
		    me = ma;
		    mb = 1;
		} else {
		    me = ma + 1;
		    mb = 2;
		}
	    }

/*           Assemble Kronecker product system of linear equations with */
/*           matrix */

/*              MAT = kron(C',A(MA:ME,MA:ME)') + kron(D',E(MA:ME,MA:ME)') */

/*           and right hand side */

/*              RHS = vec(X(MA:ME,:)). */

	    if (*n == 1) {
		dimmat = mb;
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    i__2 = mb;
		    for (j = 1; j <= i__2; ++j) {
			maj = ma + j - 1;
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[maj + 
				mai * a_dim1];
			if (maj <= mai) {
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    maj + mai * e_dim1];
			}
/* L40: */
		    }
		    rhs[i__ - 1] = x[mai + x_dim1];
/* L60: */
		}
	    } else {
		dimmat = mb << 1;
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    i__2 = mb;
		    for (j = 1; j <= i__2; ++j) {
			maj = ma + j - 1;
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[maj + 
				mai * a_dim1];
			mat[mb + i__ + (j << 2) - 5] = c__[(c_dim1 << 1) + 1] 
				* a[maj + mai * a_dim1];
			mat[i__ + (mb + j << 2) - 5] = c__[c_dim1 + 2] * a[
				maj + mai * a_dim1];
			mat[mb + i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) 
				+ 2] * a[maj + mai * a_dim1];
			if (maj <= mai) {
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    maj + mai * e_dim1];
			    mat[mb + i__ + (j << 2) - 5] += d__[(d_dim1 << 1) 
				    + 1] * e[maj + mai * e_dim1];
			    mat[i__ + (mb + j << 2) - 5] += d__[d_dim1 + 2] * 
				    e[maj + mai * e_dim1];
			    mat[mb + i__ + (mb + j << 2) - 5] += d__[(d_dim1 
				    << 1) + 2] * e[maj + mai * e_dim1];
			}
/* L80: */
		    }
		    rhs[i__ - 1] = x[mai + x_dim1];
		    rhs[mb + i__ - 1] = x[mai + (x_dim1 << 1)];
/* L100: */
		}
	    }

/*           Solve the system of linear equations. */

	    mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
	    if (info1 != 0) {
		*info = 1;
	    }
	    mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);
	    if (scale1 != 1.) {
		*scale = scale1 * *scale;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dscal_(m, &scale1, &x[i__ * x_dim1 + 1], &c__1);
/* L120: */
		}
	    }

	    if (*n == 1) {
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    x[mai + x_dim1] = rhs[i__ - 1];
/* L140: */
		}
	    } else {
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    x[mai + x_dim1] = rhs[i__ - 1];
		    x[mai + (x_dim1 << 1)] = rhs[mb + i__ - 1];
/* L160: */
		}
	    }

/*           Update right hand sides. */

/*           X(ME+1:M,:) = X(ME+1:M,:) - A(MA:ME,ME+1:M)'*X(MA:ME,:)*C */

/*           X(ME+1:M,:) = X(ME+1:M,:) - E(MA:ME,ME+1:M)'*X(MA:ME,:)*D */

	    if (me < *m) {
		dgemm_("N", "N", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			c__[c_offset], ldc, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
		i__1 = *m - me;
		dgemm_("T", "N", &i__1, n, &mb, &c_b23, &a[ma + (me + 1) * 
			a_dim1], lda, tm, &c__2, &c_b18, &x[me + 1 + x_dim1], 
			ldx, (ftnlen)1, (ftnlen)1);
		dgemm_("N", "N", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			d__[d_offset], ldd, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
		i__1 = *m - me;
		dgemm_("T", "N", &i__1, n, &mb, &c_b23, &e[ma + (me + 1) * 
			e_dim1], lde, tm, &c__2, &c_b18, &x[me + 1 + x_dim1], 
			ldx, (ftnlen)1, (ftnlen)1);
	    }

	    goto L20;
	}
/*        END WHILE 20 */

    } else {

/*        Solve equation (2). */

/*        Compute block row X(MA:ME,:). MB denotes the number of rows in */
/*        this block row. */

	ma = *m + 1;
/*        WHILE ( MA .NE. 1 ) DO */
L180:
	if (ma != 1) {
	    me = ma - 1;
	    if (me == 1) {
		ma = 1;
		mb = 1;
	    } else {
		if (a[me + (me - 1) * a_dim1] == 0.) {
		    ma = me;
		    mb = 1;
		} else {
		    ma = me - 1;
		    mb = 2;
		}
	    }

/*           Assemble Kronecker product system of linear equations with */
/*           matrix */

/*              MAT = kron(C,A(MA:ME,MA:ME)) + kron(D,E(MA:ME,MA:ME)) */

/*           and right hand side */

/*              RHS = vec(X(MA:ME,:)). */

	    if (*n == 1) {
		dimmat = mb;
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    i__2 = mb;
		    for (j = 1; j <= i__2; ++j) {
			maj = ma + j - 1;
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[mai + 
				maj * a_dim1];
			if (maj >= mai) {
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    mai + maj * e_dim1];
			}
/* L200: */
		    }
		    rhs[i__ - 1] = x[mai + x_dim1];
/* L220: */
		}
	    } else {
		dimmat = mb << 1;
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    i__2 = mb;
		    for (j = 1; j <= i__2; ++j) {
			maj = ma + j - 1;
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[mai + 
				maj * a_dim1];
			mat[mb + i__ + (j << 2) - 5] = c__[c_dim1 + 2] * a[
				mai + maj * a_dim1];
			mat[i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) + 1] 
				* a[mai + maj * a_dim1];
			mat[mb + i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) 
				+ 2] * a[mai + maj * a_dim1];
			if (maj >= mai) {
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    mai + maj * e_dim1];
			    mat[mb + i__ + (j << 2) - 5] += d__[d_dim1 + 2] * 
				    e[mai + maj * e_dim1];
			    mat[i__ + (mb + j << 2) - 5] += d__[(d_dim1 << 1) 
				    + 1] * e[mai + maj * e_dim1];
			    mat[mb + i__ + (mb + j << 2) - 5] += d__[(d_dim1 
				    << 1) + 2] * e[mai + maj * e_dim1];
			}
/* L240: */
		    }
		    rhs[i__ - 1] = x[mai + x_dim1];
		    rhs[mb + i__ - 1] = x[mai + (x_dim1 << 1)];
/* L260: */
		}
	    }

/*           Solve the system of linear equations. */

	    mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
	    if (info1 != 0) {
		*info = 1;
	    }
	    mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);
	    if (scale1 != 1.) {
		*scale = scale1 * *scale;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dscal_(m, &scale1, &x[i__ * x_dim1 + 1], &c__1);
/* L280: */
		}
	    }

	    if (*n == 1) {
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    x[mai + x_dim1] = rhs[i__ - 1];
/* L300: */
		}
	    } else {
		i__1 = mb;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    mai = ma + i__ - 1;
		    x[mai + x_dim1] = rhs[i__ - 1];
		    x[mai + (x_dim1 << 1)] = rhs[mb + i__ - 1];
/* L320: */
		}
	    }

/*           Update right hand sides. */

/*              X(1:MA-1,:) = X(1:MA-1,:) - A(1:MA-1,MA:ME)*X(MA:ME,:)*C' */

/*              X(1:MA-1,:) = X(1:MA-1,:) - E(1:MA-1,MA:ME)*X(MA:ME,:)*D' */

	    if (ma > 1) {
		dgemm_("N", "T", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			c__[c_offset], ldc, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
		i__1 = ma - 1;
		dgemm_("N", "N", &i__1, n, &mb, &c_b23, &a[ma * a_dim1 + 1], 
			lda, tm, &c__2, &c_b18, &x[x_offset], ldx, (ftnlen)1, 
			(ftnlen)1);
		dgemm_("N", "T", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			d__[d_offset], ldd, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
		i__1 = ma - 1;
		dgemm_("N", "N", &i__1, n, &mb, &c_b23, &e[ma * e_dim1 + 1], 
			lde, tm, &c__2, &c_b18, &x[x_offset], ldx, (ftnlen)1, 
			(ftnlen)1);
	    }

	    goto L180;
	}
/*        END WHILE 180 */

    }

    return 0;
/* *** Last line of SG03BW *** */
} /* sg03bw_ */

