/* SG03AY.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 1.;
static doublereal c_b12 = 0.;
static integer c__2 = 2;
static doublereal c_b16 = -1.;
static integer c__4 = 4;

/* Subroutine */ int sg03ay_(char *trans, integer *n, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *x, integer *ldx, 
	doublereal *scale, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, kb, lb, kh, lh, kl, ll;
    static doublereal tm[4]	/* was [2][2] */, ak11, ak12, ak21, ak22, 
	    al11, al12, al21, al22, ek11, ek12, ek22, el11, el12, el22, mat[
	    16]	/* was [4][4] */, rhs[4];
    static integer piv1[4], piv2[4], info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), mb02uu_(integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *), mb02uv_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *);
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

/*     To solve for X either the reduced generalized continuous-time */
/*     Lyapunov equation */

/*         T            T */
/*        A  * X * E + E  * X * A  =  SCALE * Y                       (1) */

/*     or */

/*                 T            T */
/*        A * X * E  + E * X * A   =  SCALE * Y                       (2) */

/*     where the right hand side Y is symmetric. A, E, Y, and the */
/*     solution X are N-by-N matrices. The pencil A - lambda * E must be */
/*     in generalized Schur form (A upper quasitriangular, E upper */
/*     triangular). SCALE is an output scale factor, set to avoid */
/*     overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the quasitriangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the right hand side matrix Y of the equation. Only */
/*             the upper triangular part of this matrix need be given. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the solution matrix X of the equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             (0 < SCALE <= 1) */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  equation is (almost) singular to working precision; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrices A and E are unchanged). */

/*     METHOD */

/*     The solution X of (1) or (2) is computed via block back */
/*     substitution or block forward substitution, respectively. (See */
/*     [1] and [2] for details.) */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     8/3 * N**3 flops are required by the routine. Note that we count a */
/*     single floating point arithmetic operation as one flop. */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */

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
/*     .. Executable Statements .. */

/*     Decode input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*lde < max(1,*n)) {
	*info = -6;
    } else if (*ldx < max(1,*n)) {
	*info = -8;
    } else {
	*info = 0;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SG03AY", &i__1, (ftnlen)6);
	return 0;
    }

    *scale = 1.;

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    if (notrns) {

/*        Solve equation (1). */

/*        Outer Loop. Compute block row X(KL:KH,:). KB denotes the number */
/*        of rows in this block row. */

	kl = 0;
	kb = 1;
/*        WHILE ( KL+KB .LE. N ) DO */
L20:
	if (kl + kb <= *n) {
	    kl += kb;
	    if (kl == *n) {
		kb = 1;
	    } else {
		if (a[kl + 1 + kl * a_dim1] != 0.) {
		    kb = 2;
		} else {
		    kb = 1;
		}
	    }
	    kh = kl + kb - 1;

/*           Copy elements of solution already known by symmetry. */

/*              X(KL:KH,1:KL-1) = X(1:KL-1,KL:KH)' */

	    if (kl > 1) {
		i__1 = kh;
		for (i__ = kl; i__ <= i__1; ++i__) {
		    i__2 = kl - 1;
		    dcopy_(&i__2, &x[i__ * x_dim1 + 1], &c__1, &x[i__ + 
			    x_dim1], ldx);
/* L40: */
		}
	    }

/*           Inner Loop. Compute block X(KL:KH,LL:LH). LB denotes the */
/*           number of columns in this block. */

	    ll = kl - 1;
	    lb = 1;
/*           WHILE ( LL+LB .LE. N ) DO */
L60:
	    if (ll + lb <= *n) {
		ll += lb;
		if (ll == *n) {
		    lb = 1;
		} else {
		    if (a[ll + 1 + ll * a_dim1] != 0.) {
			lb = 2;
		    } else {
			lb = 1;
		    }
		}
		lh = ll + lb - 1;

/*              Update right hand sides (I). */

/*                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) - */
/*                    A(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*E(1:LL-1,LL:LH)) */

/*                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) - */
/*                    E(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*A(1:LL-1,LL:LH)) */

		if (ll > 1) {
		    i__1 = ll - 1;
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &x[kl + x_dim1],
			     ldx, &e[ll * e_dim1 + 1], lde, &c_b12, tm, &c__2,
			     (ftnlen)1, (ftnlen)1);
		    i__1 = lh - kl + 1;
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &a[kl + kl * 
			    a_dim1], lda, tm, &c__2, &c_b11, &x[kl + ll * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		    i__1 = ll - 1;
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &x[kl + x_dim1],
			     ldx, &a[ll * a_dim1 + 1], lda, &c_b12, tm, &c__2,
			     (ftnlen)1, (ftnlen)1);
		    i__1 = lh - kh + 1;
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &e[kl + kh * 
			    e_dim1], lde, tm, &c__2, &c_b11, &x[kh + ll * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		    if (kb == 2) {
			d__1 = -e[kl + kl * e_dim1];
			daxpy_(&lb, &d__1, tm, &c__2, &x[kl + ll * x_dim1], 
				ldx);
		    }
		}

/*              Solve small Sylvester equations of order at most (2,2). */

		if (kb == 1 && lb == 1) {

		    dimmat = 1;

		    mat[0] = e[ll + ll * e_dim1] * a[kl + kl * a_dim1] + a[ll 
			    + ll * a_dim1] * e[kl + kl * e_dim1];

		    rhs[0] = x[kl + ll * x_dim1];

		} else if (kb == 2 && lb == 1) {

		    dimmat = 2;

		    ak11 = a[kl + kl * a_dim1];
		    ak12 = a[kl + kh * a_dim1];
		    ak21 = a[kh + kl * a_dim1];
		    ak22 = a[kh + kh * a_dim1];

		    al11 = a[ll + ll * a_dim1];

		    ek11 = e[kl + kl * e_dim1];
		    ek12 = e[kl + kh * e_dim1];
		    ek22 = e[kh + kh * e_dim1];

		    el11 = e[ll + ll * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = el11 * ak21;
		    mat[1] = el11 * ak12 + al11 * ek12;
		    mat[5] = el11 * ak22 + al11 * ek22;

		    rhs[0] = x[kl + ll * x_dim1];
		    rhs[1] = x[kh + ll * x_dim1];

		} else if (kb == 1 && lb == 2) {

		    dimmat = 2;

		    ak11 = a[kl + kl * a_dim1];

		    al11 = a[ll + ll * a_dim1];
		    al12 = a[ll + lh * a_dim1];
		    al21 = a[lh + ll * a_dim1];
		    al22 = a[lh + lh * a_dim1];

		    ek11 = e[kl + kl * e_dim1];

		    el11 = e[ll + ll * e_dim1];
		    el12 = e[ll + lh * e_dim1];
		    el22 = e[lh + lh * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = al21 * ek11;
		    mat[1] = el12 * ak11 + al12 * ek11;
		    mat[5] = el22 * ak11 + al22 * ek11;

		    rhs[0] = x[kl + ll * x_dim1];
		    rhs[1] = x[kl + lh * x_dim1];

		} else {

		    dimmat = 4;

		    ak11 = a[kl + kl * a_dim1];
		    ak12 = a[kl + kh * a_dim1];
		    ak21 = a[kh + kl * a_dim1];
		    ak22 = a[kh + kh * a_dim1];

		    al11 = a[ll + ll * a_dim1];
		    al12 = a[ll + lh * a_dim1];
		    al21 = a[lh + ll * a_dim1];
		    al22 = a[lh + lh * a_dim1];

		    ek11 = e[kl + kl * e_dim1];
		    ek12 = e[kl + kh * e_dim1];
		    ek22 = e[kh + kh * e_dim1];

		    el11 = e[ll + ll * e_dim1];
		    el12 = e[ll + lh * e_dim1];
		    el22 = e[lh + lh * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = el11 * ak21;
		    mat[8] = al21 * ek11;
		    mat[12] = 0.;

		    mat[1] = el11 * ak12 + al11 * ek12;
		    mat[5] = el11 * ak22 + al11 * ek22;
		    mat[9] = al21 * ek12;
		    mat[13] = al21 * ek22;

		    mat[2] = el12 * ak11 + al12 * ek11;
		    mat[6] = el12 * ak21;
		    mat[10] = el22 * ak11 + al22 * ek11;
		    mat[14] = el22 * ak21;

		    mat[3] = el12 * ak12 + al12 * ek12;
		    mat[7] = el12 * ak22 + al12 * ek22;
		    mat[11] = el22 * ak12 + al22 * ek12;
		    mat[15] = el22 * ak22 + al22 * ek22;

		    rhs[0] = x[kl + ll * x_dim1];
		    if (kl == ll) {
			rhs[1] = x[kl + kh * x_dim1];
		    } else {
			rhs[1] = x[kh + ll * x_dim1];
		    }
		    rhs[2] = x[kl + lh * x_dim1];
		    rhs[3] = x[kh + lh * x_dim1];

		}

		mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
		if (info1 != 0) {
		    *info = 1;
		}
		mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);

/*              Scaling. */

		if (scale1 != 1.) {
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(n, &scale1, &x[i__ * x_dim1 + 1], &c__1);
/* L80: */
		    }
		    *scale *= scale1;
		}

		if (lb == 1 && kb == 1) {
		    x[kl + ll * x_dim1] = rhs[0];
		} else if (lb == 1 && kb == 2) {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kh + ll * x_dim1] = rhs[1];
		} else if (lb == 2 && kb == 1) {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kl + lh * x_dim1] = rhs[1];
		} else {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kh + ll * x_dim1] = rhs[1];
		    x[kl + lh * x_dim1] = rhs[2];
		    x[kh + lh * x_dim1] = rhs[3];
		}

/*              Update right hand sides (II). */

/*                 X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) - */
/*                    A(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*E(LL:LH,LL:LH)) */

/*                 X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) - */
/*                    E(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*A(LL:LH,LL:LH)) */

		if (kl < ll) {
		    if (lb == 2) {
			dgemv_("N", &kb, &c__2, &c_b11, &x[kl + ll * x_dim1], 
				ldx, &e[ll + lh * e_dim1], &c__1, &c_b12, &tm[
				2], &c__1, (ftnlen)1);
		    }
		    dcopy_(&kb, &x[kl + ll * x_dim1], &c__1, tm, &c__1);
		    dscal_(&kb, &e[ll + ll * e_dim1], tm, &c__1);
		    i__1 = lh - kh;
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &a[kl + (kh + 1)
			     * a_dim1], lda, tm, &c__2, &c_b11, &x[kh + 1 + 
			    ll * x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		    dgemm_("N", "N", &kb, &lb, &lb, &c_b11, &x[kl + ll * 
			    x_dim1], ldx, &a[ll + ll * a_dim1], lda, &c_b12, 
			    tm, &c__2, (ftnlen)1, (ftnlen)1);
		    i__1 = lh - kh;
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &e[kl + (kh + 1)
			     * e_dim1], lde, tm, &c__2, &c_b11, &x[kh + 1 + 
			    ll * x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		}

		goto L60;
	    }
/*           END WHILE 60 */

	    goto L20;
	}
/*        END WHILE 20 */

    } else {

/*        Solve equation (2). */

/*        Outer Loop. Compute block column X(:,LL:LH). LB denotes the */
/*        number of columns in this block column. */

	ll = *n + 1;
/*        WHILE ( LL .GT. 1 ) DO */
L100:
	if (ll > 1) {
	    lh = ll - 1;
	    if (lh == 1) {
		lb = 1;
	    } else {
		if (a[ll - 1 + (ll - 2) * a_dim1] != 0.) {
		    lb = 2;
		} else {
		    lb = 1;
		}
	    }
	    ll -= lb;

/*           Copy elements of solution already known by symmetry. */

/*              X(LH+1:N,LL:LH) = X(LL:LH,LH+1:N)' */

	    if (lh < *n) {
		i__1 = lh;
		for (i__ = ll; i__ <= i__1; ++i__) {
		    i__2 = *n - lh;
		    dcopy_(&i__2, &x[i__ + (lh + 1) * x_dim1], ldx, &x[lh + 1 
			    + i__ * x_dim1], &c__1);
/* L120: */
		}
	    }

/*           Inner Loop. Compute block X(KL:KH,LL:LH). KB denotes the */
/*           number of rows in this block. */

	    kl = lh + 1;
/*           WHILE ( KL .GT. 1 ) DO */
L140:
	    if (kl > 1) {
		kh = kl - 1;
		if (kh == 1) {
		    kb = 1;
		} else {
		    if (a[kl - 1 + (kl - 2) * a_dim1] != 0.) {
			kb = 2;
		    } else {
			kb = 1;
		    }
		}
		kl -= kb;

/*              Update right hand sides (I). */

/*                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) - */
/*                    (A(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*E(KL:LH,LL:LH)' */

/*                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) - */
/*                    (E(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*A(KL:LH,LL:LH)' */

		if (kh < *n) {
		    i__1 = *n - kh;
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &a[kl + (kh + 1)
			     * a_dim1], lda, &x[kh + 1 + ll * x_dim1], ldx, &
			    c_b12, tm, &c__2, (ftnlen)1, (ftnlen)1);
		    i__1 = ll - kl + 1;
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &e[
			    kl + ll * e_dim1], lde, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		    if (lb == 2) {
			d__1 = -e[lh + lh * e_dim1];
			daxpy_(&kb, &d__1, &tm[2], &c__1, &x[kl + lh * x_dim1]
				, &c__1);
		    }
		    i__1 = *n - kh;
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &e[kl + (kh + 1)
			     * e_dim1], lde, &x[kh + 1 + ll * x_dim1], ldx, &
			    c_b12, tm, &c__2, (ftnlen)1, (ftnlen)1);
		    i__1 = lh - kl + 1;
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &a[
			    kl + ll * a_dim1], lda, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		}

/*              Solve small Sylvester equations of order at most (2,2). */

		if (kb == 1 && lb == 1) {

		    dimmat = 1;

		    mat[0] = e[ll + ll * e_dim1] * a[kl + kl * a_dim1] + a[ll 
			    + ll * a_dim1] * e[kl + kl * e_dim1];

		    rhs[0] = x[kl + ll * x_dim1];

		} else if (kb == 2 && lb == 1) {

		    dimmat = 2;

		    ak11 = a[kl + kl * a_dim1];
		    ak12 = a[kl + kh * a_dim1];
		    ak21 = a[kh + kl * a_dim1];
		    ak22 = a[kh + kh * a_dim1];

		    al11 = a[ll + ll * a_dim1];

		    ek11 = e[kl + kl * e_dim1];
		    ek12 = e[kl + kh * e_dim1];
		    ek22 = e[kh + kh * e_dim1];

		    el11 = e[ll + ll * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = el11 * ak12 + al11 * ek12;
		    mat[1] = el11 * ak21;
		    mat[5] = el11 * ak22 + al11 * ek22;

		    rhs[0] = x[kl + ll * x_dim1];
		    rhs[1] = x[kh + ll * x_dim1];

		} else if (kb == 1 && lb == 2) {

		    dimmat = 2;

		    ak11 = a[kl + kl * a_dim1];

		    al11 = a[ll + ll * a_dim1];
		    al12 = a[ll + lh * a_dim1];
		    al21 = a[lh + ll * a_dim1];
		    al22 = a[lh + lh * a_dim1];

		    ek11 = e[kl + kl * e_dim1];

		    el11 = e[ll + ll * e_dim1];
		    el12 = e[ll + lh * e_dim1];
		    el22 = e[lh + lh * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = el12 * ak11 + al12 * ek11;
		    mat[1] = al21 * ek11;
		    mat[5] = el22 * ak11 + al22 * ek11;

		    rhs[0] = x[kl + ll * x_dim1];
		    rhs[1] = x[kl + lh * x_dim1];

		} else {

		    dimmat = 4;

		    ak11 = a[kl + kl * a_dim1];
		    ak12 = a[kl + kh * a_dim1];
		    ak21 = a[kh + kl * a_dim1];
		    ak22 = a[kh + kh * a_dim1];

		    al11 = a[ll + ll * a_dim1];
		    al12 = a[ll + lh * a_dim1];
		    al21 = a[lh + ll * a_dim1];
		    al22 = a[lh + lh * a_dim1];

		    ek11 = e[kl + kl * e_dim1];
		    ek12 = e[kl + kh * e_dim1];
		    ek22 = e[kh + kh * e_dim1];

		    el11 = e[ll + ll * e_dim1];
		    el12 = e[ll + lh * e_dim1];
		    el22 = e[lh + lh * e_dim1];

		    mat[0] = el11 * ak11 + al11 * ek11;
		    mat[4] = el11 * ak12 + al11 * ek12;
		    mat[8] = el12 * ak11 + al12 * ek11;
		    mat[12] = el12 * ak12 + al12 * ek12;

		    mat[1] = el11 * ak21;
		    mat[5] = el11 * ak22 + al11 * ek22;
		    mat[9] = el12 * ak21;
		    mat[13] = el12 * ak22 + al12 * ek22;

		    mat[2] = al21 * ek11;
		    mat[6] = al21 * ek12;
		    mat[10] = el22 * ak11 + al22 * ek11;
		    mat[14] = el22 * ak12 + al22 * ek12;

		    mat[3] = 0.;
		    mat[7] = al21 * ek22;
		    mat[11] = el22 * ak21;
		    mat[15] = el22 * ak22 + al22 * ek22;

		    rhs[0] = x[kl + ll * x_dim1];
		    if (kl == ll) {
			rhs[1] = x[kl + kh * x_dim1];
		    } else {
			rhs[1] = x[kh + ll * x_dim1];
		    }
		    rhs[2] = x[kl + lh * x_dim1];
		    rhs[3] = x[kh + lh * x_dim1];

		}

		mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
		if (info1 != 0) {
		    *info = 1;
		}
		mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);

/*              Scaling. */

		if (scale1 != 1.) {
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(n, &scale1, &x[i__ * x_dim1 + 1], &c__1);
/* L160: */
		    }
		    *scale *= scale1;
		}

		if (lb == 1 && kb == 1) {
		    x[kl + ll * x_dim1] = rhs[0];
		} else if (lb == 1 && kb == 2) {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kh + ll * x_dim1] = rhs[1];
		} else if (lb == 2 && kb == 1) {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kl + lh * x_dim1] = rhs[1];
		} else {
		    x[kl + ll * x_dim1] = rhs[0];
		    x[kh + ll * x_dim1] = rhs[1];
		    x[kl + lh * x_dim1] = rhs[2];
		    x[kh + lh * x_dim1] = rhs[3];
		}

/*              Update right hand sides (II). */

/*                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) - */
/*                    (A(KL:KH,KL:KH)*X(KL:KH,LL:LH))*E(KL:LL-1,LL:LH)' */

/*                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) - */
/*                    (E(KL:KH,KL:KH)*X(KL:KH,LL:LH))*A(KL:LL-1,LL:LH)' */

		if (kl < ll) {
		    dgemm_("N", "N", &kb, &lb, &kb, &c_b11, &a[kl + kl * 
			    a_dim1], lda, &x[kl + ll * x_dim1], ldx, &c_b12, 
			    tm, &c__2, (ftnlen)1, (ftnlen)1);
		    i__1 = ll - kl;
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &e[
			    kl + ll * e_dim1], lde, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		    dgemv_("T", &kb, &lb, &c_b11, &x[kl + ll * x_dim1], ldx, &
			    e[kl + kl * e_dim1], lde, &c_b12, tm, &c__2, (
			    ftnlen)1);
		    if (kb == 2) {
			dcopy_(&lb, &x[kh + ll * x_dim1], ldx, &tm[1], &c__2);
			dscal_(&lb, &e[kh + kh * e_dim1], &tm[1], &c__2);
		    }
		    i__1 = ll - kl;
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &a[
			    kl + ll * a_dim1], lda, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
		}

		goto L140;
	    }
/*           END WHILE 140 */

	    goto L100;
	}
/*        END WHILE 100 */

    }

    return 0;
/* *** Last line of SG03AY *** */
} /* sg03ay_ */

