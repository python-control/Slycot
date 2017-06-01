/* SG03BX.f -- translated by f2c (version 20100827).
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

static integer c__2 = 2;
static doublereal c_b8 = 0.;
static doublereal c_b14 = 1.;
static doublereal c_b60 = -1.;
static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int sg03bx_(char *dico, char *trans, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *b, integer *ldb, 
	doublereal *u, integer *ldu, doublereal *scale, doublereal *m1, 
	integer *ldm1, doublereal *m2, integer *ldm2, integer *info, ftnlen 
	dico_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, m1_dim1, 
	    m1_offset, m2_dim1, m2_offset, u_dim1, u_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal l, t, v, w, aa[4]	/* was [2][2] */, b11, bb[4]	/* 
	    was [2][2] */, b22, ai[4]	/* was [2][2] */, bi[4]	/* was [2][2] 
	    */, ci, ee[4]	/* was [2][2] */, ei[4]	/* was [2][2] */, ar[
	    4]	/* was [2][2] */, br[4]	/* was [2][2] */, cr, er[4]	/* 
	    was [2][2] */, qi[4]	/* was [2][2] */, si, ti[4]	/* 
	    was [2][2] */, ui[4]	/* was [2][2] */, xi, yi, qr[4]	/* 
	    was [2][2] */, zi[4]	/* was [2][2] */, sr, tr[4]	/* 
	    was [2][2] */, ur[4]	/* was [2][2] */, xr, yr, zr[4]	/* 
	    was [2][2] */, m1i[4]	/* was [2][2] */, m2i[4]	/* 
	    was [2][2] */, m1r[4]	/* was [2][2] */, m2r[4]	/* 
	    was [2][2] */, b12i, b12r, qbi[4]	/* was [2][2] */, qbr[4]	
	    /* was [2][2] */, eps, qui[4]	/* was [2][2] */, qur[4]	
	    /* was [2][2] */, lami, lamr;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal betai, alpha;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal betar;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sg03by_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal scale1, scale2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum;
    static logical iscont;
    static doublereal smlnum;
    static logical istrns;


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

/*     To solve for X = op(U)**T * op(U) either the generalized c-stable */
/*     continuous-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(E) + op(E)  * X * op(A) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (1) */

/*     or the generalized d-stable discrete-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(A) - op(E)  * X * op(E) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (2) */

/*     where op(K) is either K or K**T for K = A, B, E, U. The Cholesky */
/*     factor U of the solution is computed without first finding X. */

/*     Furthermore, the auxiliary matrices */

/*                                   -1        -1 */
/*        M1 := op(U) * op(A) * op(E)   * op(U) */

/*                           -1        -1 */
/*        M2 := op(B) * op(E)   * op(U) */

/*     are computed in a numerically reliable way. */

/*     The matrices A, B, E, M1, M2, and U are real 2-by-2 matrices. The */
/*     pencil A - lambda * E must have a pair of complex conjugate */
/*     eigenvalues. The eigenvalues must be in the open right half plane */
/*     (in the continuous-time case) or inside the unit circle (in the */
/*     discrete-time case). */

/*     The resulting matrix U is upper triangular. The entries on its */
/*     main diagonal are non-negative. SCALE is an output scale factor */
/*     set to avoid overflow in U. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies whether the continuous-time or the discrete-time */
/*             equation is to be solved: */
/*             = 'C':  Solve continuous-time equation (1); */
/*             = 'D':  Solve discrete-time equation (2). */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  op(K) = K,     K = A, B, E, U; */
/*             = 'T':  op(K) = K**T,  K = A, B, E, U. */

/*     Input/Output Parameters */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,2) */
/*             The leading 2-by-2 part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= 2. */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,2) */
/*             The leading 2-by-2 upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= 2. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
/*             The leading 2-by-2 upper triangular part of this array */
/*             must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= 2. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2) */
/*             The leading 2-by-2 part of this array contains the upper */
/*             triangular matrix U. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U.  LDU >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in U. */
/*             0 < SCALE <= 1. */

/*     M1      (output) DOUBLE PRECISION array, dimension (LDM1,2) */
/*             The leading 2-by-2 part of this array contains the */
/*             matrix M1. */

/*     LDM1    INTEGER */
/*             The leading dimension of the array M1.  LDM1 >= 2. */

/*     M2      (output) DOUBLE PRECISION array, dimension (LDM2,2) */
/*             The leading 2-by-2 part of this array contains the */
/*             matrix M2. */

/*     LDM2    INTEGER */
/*             The leading dimension of the array M2.  LDM2 >= 2. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 2:  the eigenvalues of the pencil A - lambda * E are not */
/*                   a pair of complex conjugate numbers; */
/*             = 3:  the eigenvalues of the pencil A - lambda * E are */
/*                   not in the open right half plane (in the continuous- */
/*                   time case) or inside the unit circle (in the */
/*                   discrete-time case). */

/*     METHOD */

/*     The method used by the routine is based on a generalization of the */
/*     method due to Hammarling ([1], section 6) for Lyapunov equations */
/*     of order 2. A more detailed description is given in [2]. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-323, 1982. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     FURTHER COMMENTS */

/*     If the solution matrix U is singular, the matrices M1 and M2 are */
/*     properly set (see [1], equation (6.21)). */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */
/*     July 2003 (V. Sima; suggested by Klaus Schnepper). */
/*     Oct. 2003 (A. Varga). */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    m1_dim1 = *ldm1;
    m1_offset = 1 + m1_dim1;
    m1 -= m1_offset;
    m2_dim1 = *ldm2;
    m2_offset = 1 + m2_dim1;
    m2 -= m2_offset;

    /* Function Body */
    istrns = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    iscont = lsame_(dico, "C", (ftnlen)1, (ftnlen)1);

/*     Do not check input parameters for errors. */

/*     Set constants to control overflow. */

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1) / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

    *info = 0;
    *scale = 1.;

/*     Make copies of A, E, and B. */

    aa[0] = a[a_dim1 + 1];
    aa[1] = a[a_dim1 + 2];
    aa[2] = a[(a_dim1 << 1) + 1];
    aa[3] = a[(a_dim1 << 1) + 2];
    ee[0] = e[e_dim1 + 1];
    ee[1] = 0.;
    ee[2] = e[(e_dim1 << 1) + 1];
    ee[3] = e[(e_dim1 << 1) + 2];
    bb[0] = b[b_dim1 + 1];
    bb[1] = 0.;
    bb[2] = b[(b_dim1 << 1) + 1];
    bb[3] = b[(b_dim1 << 1) + 2];

/*     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be */
/*     solved, transpose the matrices A, E, B with respect to the */
/*     anti-diagonal. This results in a non-transposed equation. */

    if (istrns) {
	v = aa[0];
	aa[0] = aa[3];
	aa[3] = v;
	v = ee[0];
	ee[0] = ee[3];
	ee[3] = v;
	v = bb[0];
	bb[0] = bb[3];
	bb[3] = v;
    }

/*     Perform QZ-step to transform the pencil A - lambda * E to */
/*     generalized Schur form. The main diagonal of the Schur factor of E */
/*     is real and positive. */

/*     Compute eigenvalues (LAMR + LAMI * I, LAMR - LAMI * I). */

/* Computing MAX */
/* Computing MAX */
    d__2 = abs(ee[0]), d__3 = abs(ee[2]), d__2 = max(d__2,d__3), d__3 = abs(
	    ee[3]);
    d__1 = eps * max(d__2,d__3);
    t = max(d__1,smlnum);
/* Computing MIN */
    d__1 = abs(ee[0]), d__2 = abs(ee[3]);
    if (min(d__1,d__2) < t) {
	*info = 3;
	return 0;
    }
    d__1 = smlnum * eps;
    dlag2_(aa, &c__2, ee, &c__2, &d__1, &scale1, &scale2, &lamr, &w, &lami);
    if (lami <= 0.) {
	*info = 2;
	return 0;
    }

/*     Compute right orthogonal transformation matrix Q. */

    d__1 = scale1 * aa[0] - ee[0] * lamr;
    d__2 = -ee[0] * lami;
    d__3 = scale1 * aa[1];
    sg03by_(&d__1, &d__2, &d__3, &c_b8, &cr, &ci, &sr, &si, &l);
    qr[0] = cr;
    qr[2] = sr;
    qr[1] = -sr;
    qr[3] = cr;
    qi[0] = -ci;
    qi[2] = -si;
    qi[1] = -si;
    qi[3] = ci;

/*     A := Q * A */

    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qr, &c__2, aa, &c__2, &c_b8,
	     ar, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qi, &c__2, aa, &c__2, &c_b8,
	     ai, &c__2, (ftnlen)1, (ftnlen)1);

/*     E := Q * E */

    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qr, &c__2, ee, &c__2, &c_b8,
	     er, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qi, &c__2, ee, &c__2, &c_b8,
	     ei, &c__2, (ftnlen)1, (ftnlen)1);

/*     Compute left orthogonal transformation matrix Z. */

    sg03by_(&er[3], &ei[3], &er[1], &ei[1], &cr, &ci, &sr, &si, &l);
    zr[0] = cr;
    zr[2] = sr;
    zr[1] = -sr;
    zr[3] = cr;
    zi[0] = ci;
    zi[2] = -si;
    zi[1] = -si;
    zi[3] = -ci;

/*     E := E * Z */

    dgemv_("T", &c__2, &c__2, &c_b14, zr, &c__2, er, &c__2, &c_b8, tr, &c__2, 
	    (ftnlen)1);
    dgemv_("T", &c__2, &c__2, &c_b60, zi, &c__2, ei, &c__2, &c_b14, tr, &c__2,
	     (ftnlen)1);
    dgemv_("T", &c__2, &c__2, &c_b14, zi, &c__2, er, &c__2, &c_b8, ti, &c__2, 
	    (ftnlen)1);
    dgemv_("T", &c__2, &c__2, &c_b14, zr, &c__2, ei, &c__2, &c_b14, ti, &c__2,
	     (ftnlen)1);
    dcopy_(&c__2, tr, &c__2, er, &c__2);
    dcopy_(&c__2, ti, &c__2, ei, &c__2);
    er[1] = 0.;
    er[3] = l;
    ei[1] = 0.;
    ei[3] = 0.;

/*     Make main diagonal entries of E real and positive. */
/*     (Note:  Z and E are altered.) */

    v = dlapy2_(er, ei);
    dladiv_(&v, &c_b8, er, ei, &xr, &xi);
    er[0] = v;
    ei[0] = 0.;
    yr = zr[0];
    yi = zi[0];
    zr[0] = xr * yr - xi * yi;
    zi[0] = xr * yi + xi * yr;
    yr = zr[1];
    yi = zi[1];
    zr[1] = xr * yr - xi * yi;
    zi[1] = xr * yi + xi * yr;

/*     A := A * Z */

    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ar, &c__2, zr, &c__2, &c_b8,
	     tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, ai, &c__2, zi, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ar, &c__2, zi, &c__2, &c_b8,
	     ti, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ai, &c__2, zr, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
    dcopy_(&c__4, tr, &c__1, ar, &c__1);
    dcopy_(&c__4, ti, &c__1, ai, &c__1);

/*     End of QZ-step. */

/*     B := B * Z */

    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, bb, &c__2, zr, &c__2, &c_b8,
	     br, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, bb, &c__2, zi, &c__2, &c_b8,
	     bi, &c__2, (ftnlen)1, (ftnlen)1);

/*     Overwrite B with the upper triangular matrix of its */
/*     QR-factorization. The elements on the main diagonal are real */
/*     and non-negative. */

    sg03by_(br, bi, &br[1], &bi[1], &cr, &ci, &sr, &si, &l);
    qbr[0] = cr;
    qbr[2] = sr;
    qbr[1] = -sr;
    qbr[3] = cr;
    qbi[0] = -ci;
    qbi[2] = -si;
    qbi[1] = -si;
    qbi[3] = ci;
    dgemv_("N", &c__2, &c__2, &c_b14, qbr, &c__2, &br[2], &c__1, &c_b8, tr, &
	    c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b60, qbi, &c__2, &bi[2], &c__1, &c_b14, tr, &
	    c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b14, qbi, &c__2, &br[2], &c__1, &c_b8, ti, &
	    c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b14, qbr, &c__2, &bi[2], &c__1, &c_b14, ti, &
	    c__1, (ftnlen)1);
    dcopy_(&c__2, tr, &c__1, &br[2], &c__1);
    dcopy_(&c__2, ti, &c__1, &bi[2], &c__1);
    br[0] = l;
    br[1] = 0.;
    bi[0] = 0.;
    bi[1] = 0.;
    v = dlapy2_(&br[3], &bi[3]);
/* Computing MAX */
/* Computing MAX */
    d__2 = br[0], d__3 = dlapy2_(&br[2], &bi[2]);
    d__1 = eps * max(d__2,d__3);
    if (v >= max(d__1,smlnum)) {
	dladiv_(&v, &c_b8, &br[3], &bi[3], &xr, &xi);
	br[3] = v;
	yr = qbr[1];
	yi = qbi[1];
	qbr[1] = xr * yr - xi * yi;
	qbi[1] = xr * yi + xi * yr;
	yr = qbr[3];
	yi = qbi[3];
	qbr[3] = xr * yr - xi * yi;
	qbi[3] = xr * yi + xi * yr;
    } else {
	br[3] = 0.;
    }
    bi[3] = 0.;

/*     Compute the Cholesky factor of the solution of the reduced */
/*     equation. The solution may be scaled to avoid overflow. */

    if (iscont) {

/*        Continuous-time equation. */

/*        Step I:  Compute U(1,1). Set U(2,1) = 0. */

	v = (ar[0] * er[0] + ai[0] * ei[0]) * -2.;
	if (v <= 0.) {
	    *info = 3;
	    return 0;
	}
	v = sqrt(v);
	t = abs(br[0]) * 2. * smlnum;
	if (t > v) {
	    scale1 = v / t;
	    *scale = scale1 * *scale;
	    br[0] = scale1 * br[0];
	    br[2] = scale1 * br[2];
	    bi[2] = scale1 * bi[2];
	    br[3] = scale1 * br[3];
	}
	ur[0] = br[0] / v;
	ui[0] = 0.;
	ur[1] = 0.;
	ui[1] = 0.;

/*        Step II:  Compute U(1,2). */

/* Computing MAX */
/* Computing MAX */
	d__2 = br[3], d__3 = dlapy2_(&br[2], &bi[2]);
	d__1 = eps * max(d__2,d__3);
	t = max(d__1,smlnum);
	if (abs(br[0]) < t) {
	    ur[2] = 0.;
	    ui[2] = 0.;
	} else {
	    xr = ar[0] * er[2] + ai[0] * ei[2];
	    xi = ai[0] * er[2] - ar[0] * ei[2];
	    xr = xr + ar[2] * er[0] + ai[2] * ei[0];
	    xi = xi - ai[2] * er[0] + ar[2] * ei[0];
	    xr = -br[2] * v - xr * ur[0];
	    xi = bi[2] * v - xi * ur[0];
	    yr = ar[3] * er[0] + ai[3] * ei[0];
	    yi = -ai[3] * er[0] + ar[3] * ei[0];
	    yr = yr + er[3] * ar[0] + ei[3] * ai[0];
	    yi = yi - ei[3] * ar[0] + er[3] * ai[0];
	    t = dlapy2_(&xr, &xi) * 2. * smlnum;
	    if (t > dlapy2_(&yr, &yi)) {
		scale1 = dlapy2_(&yr, &yi) / t;
		*scale = scale1 * *scale;
		br[0] = scale1 * br[0];
		br[2] = scale1 * br[2];
		bi[2] = scale1 * bi[2];
		br[3] = scale1 * br[3];
		ur[0] = scale1 * ur[0];
		xr = scale1 * xr;
		xi = scale1 * xi;
	    }
	    dladiv_(&xr, &xi, &yr, &yi, &ur[2], &ui[2]);
	    ui[2] = -ui[2];
	}

/*        Step III:  Compute U(2,2). */

	xr = (er[2] * ur[0] + er[3] * ur[2] - ei[3] * ui[2]) * v;
	xi = (-ei[2] * ur[0] - er[3] * ui[2] - ei[3] * ur[2]) * v;
	t = dlapy2_(&xr, &xi) * 2. * smlnum;
	if (t > dlapy2_(er, ei)) {
	    scale1 = dlapy2_(er, ei) / t;
	    *scale = scale1 * *scale;
	    ur[0] = scale1 * ur[0];
	    ur[2] = scale1 * ur[2];
	    ui[2] = scale1 * ui[2];
	    br[0] = scale1 * br[0];
	    br[2] = scale1 * br[2];
	    bi[2] = scale1 * bi[2];
	    br[3] = scale1 * br[3];
	    xr = scale1 * xr;
	    xi = scale1 * xi;
	}
	d__1 = -ei[0];
	dladiv_(&xr, &xi, er, &d__1, &yr, &yi);
	yr = br[2] - yr;
	yi = -bi[2] - yi;
	v = (ar[3] * er[3] + ai[3] * ei[3]) * -2.;
	if (v <= 0.) {
	    *info = 3;
	    return 0;
	}
	v = sqrt(v);
	d__1 = dlapy2_(&br[3], &bi[3]);
	d__2 = dlapy2_(&yr, &yi);
	w = dlapy2_(&d__1, &d__2);
	t = w * 2. * smlnum;
	if (t > v) {
	    scale1 = v / t;
	    *scale = scale1 * *scale;
	    ur[0] = scale1 * ur[0];
	    ur[2] = scale1 * ur[2];
	    ui[2] = scale1 * ui[2];
	    br[0] = scale1 * br[0];
	    br[2] = scale1 * br[2];
	    bi[2] = scale1 * bi[2];
	    br[3] = scale1 * br[3];
	    w = scale1 * w;
	}
	ur[3] = w / v;
	ui[3] = 0.;

/*        Compute matrices M1 and M2 for the reduced equation. */

	m1r[1] = 0.;
	m1i[1] = 0.;
	m2r[1] = 0.;
	m2i[1] = 0.;
	dladiv_(ar, ai, er, ei, &betar, &betai);
	m1r[0] = betar;
	m1i[0] = betai;
	m1r[3] = betar;
	m1i[3] = -betai;
	alpha = sqrt(betar * -2.);
	m2r[0] = alpha;
	m2i[0] = 0.;
	v = er[0] * er[3];
	xr = (-br[0] * er[2] + er[0] * br[2]) / v;
	xi = (-br[0] * ei[2] + er[0] * bi[2]) / v;
	yr = xr - alpha * ur[2];
	yi = -xi + alpha * ui[2];
	if (yr != 0. || yi != 0.) {
	    m2r[2] = yr / ur[3];
	    m2i[2] = -yi / ur[3];
	    m2r[3] = br[3] / (er[3] * ur[3]);
	    m2i[3] = 0.;
	    m1r[2] = -alpha * m2r[2];
	    m1i[2] = -alpha * m2i[2];
	} else {
	    m2r[2] = 0.;
	    m2i[2] = 0.;
	    m2r[3] = alpha;
	    m2i[3] = 0.;
	    m1r[2] = 0.;
	    m1i[2] = 0.;
	}
    } else {

/*        Discrete-time equation. */

/*        Step I:  Compute U(1,1). Set U(2,1) = 0. */

/* Computing 2nd power */
	d__1 = er[0];
/* Computing 2nd power */
	d__2 = ei[0];
/* Computing 2nd power */
	d__3 = ar[0];
/* Computing 2nd power */
	d__4 = ai[0];
	v = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
	if (v <= 0.) {
	    *info = 3;
	    return 0;
	}
	v = sqrt(v);
	t = abs(br[0]) * 2. * smlnum;
	if (t > v) {
	    scale1 = v / t;
	    *scale = scale1 * *scale;
	    br[0] = scale1 * br[0];
	    br[2] = scale1 * br[2];
	    bi[2] = scale1 * bi[2];
	    br[3] = scale1 * br[3];
	}
	ur[0] = br[0] / v;
	ui[0] = 0.;
	ur[1] = 0.;
	ui[1] = 0.;

/*        Step II:  Compute U(1,2). */

/* Computing MAX */
/* Computing MAX */
	d__2 = br[3], d__3 = dlapy2_(&br[2], &bi[2]);
	d__1 = eps * max(d__2,d__3);
	t = max(d__1,smlnum);
	if (abs(br[0]) < t) {
	    ur[2] = 0.;
	    ui[2] = 0.;
	} else {
	    xr = ar[0] * ar[2] + ai[0] * ai[2];
	    xi = ai[0] * ar[2] - ar[0] * ai[2];
	    xr = xr - er[2] * er[0] - ei[2] * ei[0];
	    xi = xi + ei[2] * er[0] - er[2] * ei[0];
	    xr = -br[2] * v - xr * ur[0];
	    xi = bi[2] * v - xi * ur[0];
	    yr = ar[3] * ar[0] + ai[3] * ai[0];
	    yi = -ai[3] * ar[0] + ar[3] * ai[0];
	    yr = yr - er[3] * er[0] - ei[3] * ei[0];
	    yi = yi + ei[3] * er[0] - er[3] * ei[0];
	    t = dlapy2_(&xr, &xi) * 2. * smlnum;
	    if (t > dlapy2_(&yr, &yi)) {
		scale1 = dlapy2_(&yr, &yi) / t;
		*scale = scale1 * *scale;
		br[0] = scale1 * br[0];
		br[2] = scale1 * br[2];
		bi[2] = scale1 * bi[2];
		br[3] = scale1 * br[3];
		ur[0] = scale1 * ur[0];
		xr = scale1 * xr;
		xi = scale1 * xi;
	    }
	    dladiv_(&xr, &xi, &yr, &yi, &ur[2], &ui[2]);
	    ui[2] = -ui[2];
	}

/*        Step III:  Compute U(2,2). */

	xr = er[2] * ur[0] + er[3] * ur[2] - ei[3] * ui[2];
	xi = -ei[2] * ur[0] - er[3] * ui[2] - ei[3] * ur[2];
	yr = ar[2] * ur[0] + ar[3] * ur[2] - ai[3] * ui[2];
	yi = -ai[2] * ur[0] - ar[3] * ui[2] - ai[3] * ur[2];
/* Computing 2nd power */
	d__1 = er[3];
/* Computing 2nd power */
	d__2 = ei[3];
/* Computing 2nd power */
	d__3 = ar[3];
/* Computing 2nd power */
	d__4 = ai[3];
	v = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
	if (v <= 0.) {
	    *info = 3;
	    return 0;
	}
	v = sqrt(v);
/* Computing MAX */
	d__1 = abs(br[3]), d__2 = abs(br[2]), d__1 = max(d__1,d__2), d__2 = 
		abs(bi[2]), d__1 = max(d__1,d__2), d__2 = abs(xr), d__1 = max(
		d__1,d__2), d__2 = abs(xi), d__1 = max(d__1,d__2), d__2 = abs(
		yr), d__1 = max(d__1,d__2), d__2 = abs(yi);
	t = max(d__1,d__2);
	if (t <= smlnum) {
	    t = 1.;
	}
/* Computing 2nd power */
	d__1 = br[3] / t;
/* Computing 2nd power */
	d__2 = br[2] / t;
/* Computing 2nd power */
	d__3 = bi[2] / t;
/* Computing 2nd power */
	d__4 = xr / t;
/* Computing 2nd power */
	d__5 = xi / t;
/* Computing 2nd power */
	d__6 = yr / t;
/* Computing 2nd power */
	d__7 = yi / t;
	w = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - d__4 * d__4 - d__5 * 
		d__5 + d__6 * d__6 + d__7 * d__7;
	if (w < 0.) {
	    *info = 3;
	    return 0;
	}
	w = t * sqrt(w);
	t = w * 2. * smlnum;
	if (t > v) {
	    scale1 = v / t;
	    *scale = scale1 * *scale;
	    ur[0] = scale1 * ur[0];
	    ur[2] = scale1 * ur[2];
	    ui[2] = scale1 * ui[2];
	    br[0] = scale1 * br[0];
	    br[2] = scale1 * br[2];
	    bi[2] = scale1 * bi[2];
	    br[3] = scale1 * br[3];
	    w = scale1 * w;
	}
	ur[3] = w / v;
	ui[3] = 0.;

/*        Compute matrices M1 and M2 for the reduced equation. */

	b11 = br[0] / er[0];
	t = er[0] * er[3];
	b12r = (er[0] * br[2] - br[0] * er[2]) / t;
	b12i = (er[0] * bi[2] - br[0] * ei[2]) / t;
	b22 = br[3] / er[3];
	m1r[1] = 0.;
	m1i[1] = 0.;
	m2r[1] = 0.;
	m2i[1] = 0.;
	dladiv_(ar, ai, er, ei, &betar, &betai);
	m1r[0] = betar;
	m1i[0] = betai;
	m1r[3] = betar;
	m1i[3] = -betai;
	v = dlapy2_(&betar, &betai);
	alpha = sqrt((1. - v) * (v + 1.));
	m2r[0] = alpha;
	m2i[0] = 0.;
	xr = (ai[0] * ei[2] - ar[0] * er[2]) / t + ar[2] / er[3];
	xi = (ar[0] * ei[2] + ai[0] * er[2]) / t - ai[2] / er[3];
	xr = betai * -2. * b12i - b11 * xr;
	xi = betai * -2. * b12r - b11 * xi;
	v = (betai - betar) * (betai + betar) + 1.;
	w = betai * -2. * betar;
	dladiv_(&xr, &xi, &v, &w, &yr, &yi);
	if (yr != 0. || yi != 0.) {
	    m2r[2] = (yr * betar - yi * betai) / ur[3];
	    m2i[2] = -(yi * betar + yr * betai) / ur[3];
	    m2r[3] = b22 / ur[3];
	    m2i[3] = 0.;
	    m1r[2] = -alpha * yr / ur[3];
	    m1i[2] = alpha * yi / ur[3];
	} else {
	    m2r[2] = 0.;
	    m2i[2] = 0.;
	    m2r[3] = alpha;
	    m2i[3] = 0.;
	    m1r[2] = 0.;
	    m1i[2] = 0.;
	}
    }

/*     Transform U back:  U := U * Q. */
/*     (Note:  Z is used as workspace.) */

    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ur, &c__2, qr, &c__2, &c_b8,
	     zr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, ui, &c__2, qi, &c__2, &
	    c_b14, zr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ur, &c__2, qi, &c__2, &c_b8,
	     zi, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ui, &c__2, qr, &c__2, &
	    c_b14, zi, &c__2, (ftnlen)1, (ftnlen)1);

/*     Overwrite U with the upper triangular matrix of its */
/*     QR-factorization. The elements on the main diagonal are real */
/*     and non-negative. */

    sg03by_(zr, zi, &zr[1], &zi[1], &cr, &ci, &sr, &si, &l);
    qur[0] = cr;
    qur[2] = sr;
    qur[1] = -sr;
    qur[3] = cr;
    qui[0] = -ci;
    qui[2] = -si;
    qui[1] = -si;
    qui[3] = ci;
    dgemv_("N", &c__2, &c__2, &c_b14, qur, &c__2, &zr[2], &c__1, &c_b8, &u[(
	    u_dim1 << 1) + 1], &c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b60, qui, &c__2, &zi[2], &c__1, &c_b14, &u[(
	    u_dim1 << 1) + 1], &c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b14, qui, &c__2, &zr[2], &c__1, &c_b8, &ui[2]
	    , &c__1, (ftnlen)1);
    dgemv_("N", &c__2, &c__2, &c_b14, qur, &c__2, &zi[2], &c__1, &c_b14, &ui[
	    2], &c__1, (ftnlen)1);
    u[u_dim1 + 1] = l;
    u[u_dim1 + 2] = 0.;
    v = dlapy2_(&u[(u_dim1 << 1) + 2], &ui[3]);
    if (v != 0.) {
	dladiv_(&v, &c_b8, &u[(u_dim1 << 1) + 2], &ui[3], &xr, &xi);
	yr = qur[1];
	yi = qui[1];
	qur[1] = xr * yr - xi * yi;
	qui[1] = xr * yi + xi * yr;
	yr = qur[3];
	yi = qui[3];
	qur[3] = xr * yr - xi * yi;
	qui[3] = xr * yi + xi * yr;
    }
    u[(u_dim1 << 1) + 2] = v;

/*     Transform the matrices M1 and M2 back. */

/*        M1 := QU * M1 * QU**H */
/*        M2 := QB**H * M2 * QU**H */

    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1r, &c__2, qur, &c__2, &
	    c_b8, tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1i, &c__2, qui, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b60, m1r, &c__2, qui, &c__2, &
	    c_b8, ti, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1i, &c__2, qur, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qur, &c__2, tr, &c__2, &
	    c_b8, &m1[m1_offset], ldm1, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, qui, &c__2, ti, &c__2, &
	    c_b14, &m1[m1_offset], ldm1, (ftnlen)1, (ftnlen)1);

    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2r, &c__2, qur, &c__2, &
	    c_b8, tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2i, &c__2, qui, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b60, m2r, &c__2, qui, &c__2, &
	    c_b8, ti, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2i, &c__2, qur, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
    dgemm_("T", "N", &c__2, &c__2, &c__2, &c_b14, qbr, &c__2, tr, &c__2, &
	    c_b8, &m2[m2_offset], ldm2, (ftnlen)1, (ftnlen)1);
    dgemm_("T", "N", &c__2, &c__2, &c__2, &c_b14, qbi, &c__2, ti, &c__2, &
	    c_b14, &m2[m2_offset], ldm2, (ftnlen)1, (ftnlen)1);

/*     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be */
/*     solved, transpose the matrix U with respect to the */
/*     anti-diagonal and the matrices M1, M2 with respect to the diagonal */
/*     and the anti-diagonal. */

    if (istrns) {
	v = u[u_dim1 + 1];
	u[u_dim1 + 1] = u[(u_dim1 << 1) + 2];
	u[(u_dim1 << 1) + 2] = v;
	v = m1[m1_dim1 + 1];
	m1[m1_dim1 + 1] = m1[(m1_dim1 << 1) + 2];
	m1[(m1_dim1 << 1) + 2] = v;
	v = m2[m2_dim1 + 1];
	m2[m2_dim1 + 1] = m2[(m2_dim1 << 1) + 2];
	m2[(m2_dim1 << 1) + 2] = v;
    }

    return 0;
/* *** Last line of SG03BX *** */
} /* sg03bx_ */

