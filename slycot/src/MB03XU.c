/* MB03XU.f -- translated by f2c (version 20100827).
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
static doublereal c_b7 = 1.;
static doublereal c_b9 = 0.;

/* Subroutine */ int mb03xu_(logical *ltra, logical *ltrb, integer *n, 
	integer *k, integer *nb, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *g, integer *ldg, doublereal *q, integer *
	ldq, doublereal *xa, integer *ldxa, doublereal *xb, integer *ldxb, 
	doublereal *xg, integer *ldxg, doublereal *xq, integer *ldxq, 
	doublereal *ya, integer *ldya, doublereal *yb, integer *ldyb, 
	doublereal *yg, integer *ldyg, doublereal *yq, integer *ldyq, 
	doublereal *csl, doublereal *csr, doublereal *taul, doublereal *taur, 
	doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, xa_dim1, xa_offset, xb_dim1, xb_offset, xg_dim1, 
	    xg_offset, xq_dim1, xq_offset, ya_dim1, ya_offset, yb_dim1, 
	    yb_offset, yg_dim1, yg_offset, yq_dim1, yq_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer nb1, nb2, nb3, pdw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, tauq;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), dlarfg_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


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

/*     To reduce 2*nb columns and rows of a real (k+2n)-by-(k+2n) */
/*     matrix H: */

/*             [ op(A)   G   ] */
/*         H = [             ], */
/*             [  Q    op(B) ] */

/*     so that elements in the first nb columns below the k-th */
/*     subdiagonal of the (k+n)-by-n matrix op(A), in the first nb */
/*     columns and rows of the n-by-n matrix Q and in the first nb rows */
/*     above the diagonal of the n-by-(k+n) matrix op(B) are zero. */
/*     The reduction is performed by orthogonal symplectic */
/*     transformations UU'*H*VV and matrices U, V, YA, YB, YG, YQ, XA, */
/*     XB, XG, and XQ are returned so that */

/*                    [ op(Aout)+U*YA'+XA*V'     G+U*YG'+XG*V'    ] */
/*         UU' H VV = [                                           ]. */
/*                    [   Qout+U*YQ'+XQ*V'   op(Bout)+U*YB'+XB*V' ] */

/*     This is an auxiliary routine called by MB04TB. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRA    LOGICAL */
/*             Specifies the form of op( A ) as follows: */
/*             = .FALSE.:  op( A ) = A; */
/*             = .TRUE.:   op( A ) = A'. */

/*     LTRB    LOGICAL */
/*             Specifies the form of op( B ) as follows: */
/*             = .FALSE.:  op( B ) = B; */
/*             = .TRUE.:   op( B ) = B'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix Q. N >= 0. */

/*     K       (input) INTEGER */
/*             The offset of the reduction. Elements below the K-th */
/*             subdiagonal in the first NB columns of op(A) are */
/*             reduced to zero. K >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns/rows to be reduced. N > NB >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDA,N)     if LTRA = .FALSE. */
/*                     (LDA,K+N)   if LTRA = .TRUE. */
/*             On entry with LTRA = .FALSE., the leading (K+N)-by-N part */
/*             of this array must contain the matrix A. */
/*             On entry with LTRA = .TRUE., the leading N-by-(K+N) part */
/*             of this array must contain the matrix A. */
/*             On exit with LTRA = .FALSE., the leading (K+N)-by-N part */
/*             of this array contains the matrix Aout and, in the zero */
/*             parts, information about the elementary reflectors used to */
/*             compute the reduction. */
/*             On exit with LTRA = .TRUE., the leading N-by-(K+N) part of */
/*             this array contains the matrix Aout and in the zero parts */
/*             information about the elementary reflectors. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,K+N),  if LTRA = .FALSE.; */
/*             LDA >= MAX(1,N),    if LTRA = .TRUE.. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDB,K+N)   if LTRB = .FALSE. */
/*                     (LDB,N)     if LTRB = .TRUE. */
/*             On entry with LTRB = .FALSE., the leading N-by-(K+N) part */
/*             of this array must contain the matrix B. */
/*             On entry with LTRB = .TRUE., the leading (K+N)-by-N part */
/*             of this array must contain the matrix B. */
/*             On exit with LTRB = .FALSE., the leading N-by-(K+N) part */
/*             of this array contains the matrix Bout and, in the zero */
/*             parts, information about the elementary reflectors used to */
/*             compute the reduction. */
/*             On exit with LTRB = .TRUE., the leading (K+N)-by-N part of */
/*             this array contains the matrix Bout and in the zero parts */
/*             information about the elementary reflectors. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N),    if LTRB = .FALSE.; */
/*             LDB >= MAX(1,K+N),  if LTRB = .TRUE.. */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix G. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Gout. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix Q. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Qout and in the zero parts information about */
/*             the elementary reflectors used to compute the reduction. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XA. */

/*     LDXA    INTEGER */
/*             The leading dimension of the array XA.  LDXA >= MAX(1,N). */

/*     XB      (output) DOUBLE PRECISION array, dimension (LDXB,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix XB. */

/*     LDXB    INTEGER */
/*             The leading dimension of the array XB. LDXB >= MAX(1,K+N). */

/*     XG      (output) DOUBLE PRECISION array, dimension (LDXG,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix XG. */

/*     LDXG    INTEGER */
/*             The leading dimension of the array XG. LDXG >= MAX(1,K+N). */

/*     XQ      (output) DOUBLE PRECISION array, dimension (LDXQ,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XQ. */

/*     LDXQ    INTEGER */
/*             The leading dimension of the array XQ.  LDXQ >= MAX(1,N). */

/*     YA      (output) DOUBLE PRECISION array, dimension (LDYA,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix YA. */

/*     LDYA    INTEGER */
/*             The leading dimension of the array YA. LDYA >= MAX(1,K+N). */

/*     YB      (output) DOUBLE PRECISION array, dimension (LDYB,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix YB. */

/*     LDYB    INTEGER */
/*             The leading dimension of the array YB.  LDYB >= MAX(1,N). */

/*     YG      (output) DOUBLE PRECISION array, dimension (LDYG,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix YG. */

/*     LDYG    INTEGER */
/*             The leading dimension of the array YG. LDYG >= MAX(1,K+N). */

/*     YQ      (output) DOUBLE PRECISION array, dimension (LDYQ,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix YQ. */

/*     LDYQ    INTEGER */
/*             The leading dimension of the array YQ.  LDYQ >= MAX(1,N). */

/*     CSL     (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2NB elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the left-hand side used to compute the */
/*             reduction. */

/*     CSR     (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2NB-2 elements of this array contain */
/*             the cosines and sines of the symplectic Givens rotations */
/*             applied from the right-hand side used to compute the */
/*             reduction. */

/*     TAUL    (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the left-hand side. */

/*     TAUR    (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the right-hand side. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (5*NB) */

/*     METHOD */

/*     For details regarding the representation of the orthogonal */
/*     symplectic matrices UU and VV within the arrays A, B, CSL, CSR, Q, */
/*     TAUL and TAUR see the description of MB04TB. */

/*     The contents of A, B, G and Q on exit are illustrated by the */
/*     following example with op(A) = A, op(B) = B, n = 5, k = 2 and */
/*     nb = 2: */

/*          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( r  r  r  r  r  )       ( r  r  r  r  r  r  r  ) */
/*      A = ( u2 r  r  r  r  ),  G = ( r  r  r  r  r  r  r  ), */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */

/*          ( t  t  v1 v1 v1 )       ( r  r  r  r  r  v2 v2 ) */
/*          ( u1 t  t  v1 v1 )       ( r  r  r  r  r  r  v2 ) */
/*      Q = ( u1 u1 r  q  q  ),  B = ( b  b  b  r  r  b  b  ). */
/*          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  ) */
/*          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  ) */

/*     where a, b, g and q denote elements of the original matrices, r */
/*     denotes a modified element, t denotes a scalar factor of an */
/*     applied elementary reflector, ui and vi denote elements of the */
/*     matrices U and V, respectively. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires ( 16*K + 32*N + 42 )*N*NB + */
/*     ( 16*K + 112*N - 208/3*NB - 69 )*NB*NB - 29/3*NB floating point */
/*     operations and is numerically backward stable. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. */
/*         Numer. Math., Vol. 78 (3), pp. 329-358, 1998. */

/*     [2] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT Numerical Mathematics, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLASUB). */

/*     KEYWORDS */

/*     Elementary matrix operations, Matrix decompositions, Hamiltonian */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    xa_dim1 = *ldxa;
    xa_offset = 1 + xa_dim1;
    xa -= xa_offset;
    xb_dim1 = *ldxb;
    xb_offset = 1 + xb_dim1;
    xb -= xb_offset;
    xg_dim1 = *ldxg;
    xg_offset = 1 + xg_dim1;
    xg -= xg_offset;
    xq_dim1 = *ldxq;
    xq_offset = 1 + xq_dim1;
    xq -= xq_offset;
    ya_dim1 = *ldya;
    ya_offset = 1 + ya_dim1;
    ya -= ya_offset;
    yb_dim1 = *ldyb;
    yb_offset = 1 + yb_dim1;
    yb -= yb_offset;
    yg_dim1 = *ldyg;
    yg_offset = 1 + yg_dim1;
    yg -= yg_offset;
    yq_dim1 = *ldyq;
    yq_offset = 1 + yq_dim1;
    yq -= yq_offset;
    --csl;
    --csr;
    --taul;
    --taur;
    --dwork;

    /* Function Body */
    if (*n + *k <= 0) {
	return 0;
    }

    nb1 = *nb + 1;
    nb2 = *nb + *nb;
    nb3 = nb2 + *nb;
    pdw = nb3 + *nb + 1;

    if (*ltra && *ltrb) {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first row/column of A and Q. See routine MB04TS. */

	    alpha = q[i__ + i__ * q_dim1];
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
	    q[i__ + i__ * q_dim1] = 1.;
	    i__2 = *n - i__ + 1;
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[i__ 
		    + (*k + i__) * a_dim1], lda);
	    i__2 = *n - i__ + 1;
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[i__ + (*k 
		    + i__) * a_dim1], lda);
	    temp = a[i__ + (*k + i__) * a_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &a[i__ + (*k + i__) * a_dim1]);
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &a[i__ + (*k + i__) * a_dim1], &a[i__ + (*k + i__ 
		    + 1) * a_dim1], lda, &taul[i__]);
	    temp = a[i__ + (*k + i__) * a_dim1];
	    a[i__ + (*k + i__) * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__) * 
		    a_dim1 + 1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[*nb + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &q[
		    i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);

/*           Update XA with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9,
		     &xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[i__ + 1 
		    + (*k + i__) * a_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &c_b7, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &a[i__ + (*k + i__ + 1) * a_dim1], lda, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &a[i__ + (*k + i__) * a_dim1], lda, &
		    c_b9, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1);

/*           Update XG with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + i__ * xg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &g[*k 
		    + i__ + g_dim1], ldg, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &c_b7, 
		    &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)12);

/*           Update XB with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + i__ * xb_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ * b_dim1 + 1], &c__1, 
		    (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &b[i__ 
		    * b_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[*k + i__ + 1 
		    + i__ * b_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &b[*
		    k + i__ + 1 + i__ * b_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ G(k+i,:); B(:,i)' ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &
		    c__1, &c__, &s);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + j * yg_dim1] = 0.;
/* L10: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
/* L20: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + j * ya_dim1] = 0.;
/* L30: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
/* L40: */
	    }

/*           Update XG with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xb[(i__ 
		    + *nb) * xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ * b_dim1 + 1], &c__1);

	    a[i__ + (*k + i__) * a_dim1] = temp;
	    q[i__ + i__ * q_dim1] = tauq;
	    csl[(i__ << 1) - 1] = c__;
	    csl[i__ * 2] = s;

/*           Transform first row/column of Q and B. */

	    alpha = q[i__ + (i__ + 1) * q_dim1];
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    *k + i__ + 1 + i__ * b_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[*k + 
		    i__ + 1 + i__ * b_dim1], &c__1);
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &b[*k + i__ + 1 + i__ * b_dim1]);
	    s = -s;
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &b[*k + i__ + 1 + i__ * b_dim1], &b[*k + i__ + 2 + 
		    i__ * b_dim1], &c__1, &taur[i__]);
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
	    b[*k + i__ + 1 + i__ * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &
		    yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + b_dim1]
		    , ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[*
		    nb + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb,
		     &c_b7, &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &b[
		    *k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);

/*           Update YQ with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &
		    yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &q[
		    i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[*k + i__ 
		    + 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
/* L50: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
/* L60: */
	    }

/*           Update YB with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &b[*k + i__ + 1 + i__ * b_dim1], &c__1,
		     &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 2 + b_dim1]
		    , ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb);

/*           Update YQ with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[i__ * 
		    ya_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[i__ + 1 + (*
		    k + i__ + 1) * a_dim1], lda, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &
		    c_b7, &a[i__ + 1 + (*k + i__ + 1) * a_dim1], lda, (ftnlen)
		    9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[i__ + 1 + 
		    a_dim1], lda, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &a[i__ + 1 
		    + a_dim1], lda, (ftnlen)12);

/*           Update YG with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + i__ * yg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg,
		     &c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1,
		     (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &g[(*k + 
		    i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L70: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L80: */
	    }

/*           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &a[i__ + 1 + a_dim1], lda, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &c_b9, &ya[(
		    i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[
		    i__ + 1 + a_dim1], lda);

/*           Update YG with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

	    b[*k + i__ + 1 + i__ * b_dim1] = temp;
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
	    csr[(i__ << 1) - 1] = c__;
	    csr[i__ * 2] = s;
/* L90: */
	}
    } else if (*ltra) {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first row/column of A and Q. See routine MB04TS. */

	    alpha = q[i__ + i__ * q_dim1];
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
	    q[i__ + i__ * q_dim1] = 1.;
	    i__2 = *n - i__ + 1;
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[i__ 
		    + (*k + i__) * a_dim1], lda);
	    i__2 = *n - i__ + 1;
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[i__ + (*k 
		    + i__) * a_dim1], lda);
	    temp = a[i__ + (*k + i__) * a_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &a[i__ + (*k + i__) * a_dim1]);
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &a[i__ + (*k + i__) * a_dim1], &a[i__ + (*k + i__ 
		    + 1) * a_dim1], lda, &taul[i__]);
	    temp = a[i__ + (*k + i__) * a_dim1];
	    a[i__ + (*k + i__) * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__) * 
		    a_dim1 + 1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[*nb + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)9);

/*           Update XA with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9,
		     &xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[i__ + 1 
		    + (*k + i__) * a_dim1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &
		    c_b7, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)9)
		    ;

/*           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &a[i__ + (*k + i__ + 1) * a_dim1], lda, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__ + 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &a[i__ + (*k + i__) * a_dim1], lda, &
		    c_b9, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1);

/*           Update XG with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &g[*k 
		    + i__ + g_dim1], ldg, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &
		    c_b7, &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (
		    ftnlen)9);

/*           Update XB with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * xb_dim1 + 
		    1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ + b_dim1], ldb, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &b[i__ 
		    + b_dim1], ldb, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[i__ + (*k + 
		    i__ + 1) * b_dim1], ldb, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &
		    b[i__ + (*k + i__ + 1) * b_dim1], ldb, (ftnlen)9);

/*           Apply rotation to [ G(k+i,:); B(i,:) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &
		    c__, &s);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + j * yg_dim1] = 0.;
/* L100: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
/* L110: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + j * ya_dim1] = 0.;
/* L120: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
/* L130: */
	    }

/*           Update XG with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xb[(i__ + *nb) 
		    * xb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ + b_dim1], ldb);

	    a[i__ + (*k + i__) * a_dim1] = temp;
	    q[i__ + i__ * q_dim1] = tauq;
	    csl[(i__ << 1) - 1] = c__;
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

	    alpha = q[i__ + (i__ + 1) * q_dim1];
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    i__ + (*k + i__ + 1) * b_dim1], ldb);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[i__ + (
		    *k + i__ + 1) * b_dim1], ldb);
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (*k + i__ + 1) * b_dim1]
		    );
	    s = -s;
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &b[i__ + (*k + i__ + 1) * b_dim1], &b[i__ + (*k + 
		    i__ + 2) * b_dim1], ldb, &taur[i__]);
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
	    b[i__ + (*k + i__ + 1) * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], 
		    ldq, &c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)
		    12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &
		    yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &dwork[*nb + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb,
		     &c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[i__ + 
		    1 + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);

/*           Update YQ with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &
		    yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12)
		    ;

/*           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, &c__, &s);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
/* L140: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
/* L150: */
	    }

/*           Update YB with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &b[i__ + (*k + i__ + 1) * b_dim1]
		    , ldb, &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1,
		     (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 2) * 
		    b_dim1 + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1);

/*           Update YQ with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[i__ * 
		    ya_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[i__ + 1 + (*
		    k + i__ + 1) * a_dim1], lda, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &
		    c_b7, &a[i__ + 1 + (*k + i__ + 1) * a_dim1], lda, (ftnlen)
		    9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[i__ + 1 + 
		    a_dim1], lda, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &a[
		    i__ + 1 + a_dim1], lda, (ftnlen)12);

/*           Update YG with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + i__ * yg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg,
		     &c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1,
		     (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &g[
		    (*k + i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L160: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L170: */
	    }

/*           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &a[i__ + 1 + a_dim1], lda, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &c_b9, &ya[(
		    i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[
		    i__ + 1 + a_dim1], lda);

/*           Update YG with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

	    b[i__ + (*k + i__ + 1) * b_dim1] = temp;
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
	    csr[(i__ << 1) - 1] = c__;
	    csr[i__ * 2] = s;
/* L180: */
	}

    } else if (*ltrb) {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first columns of A and Q. See routine MB04TS. */

	    alpha = q[i__ + i__ * q_dim1];
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
	    q[i__ + i__ * q_dim1] = 1.;
	    i__2 = *n - i__ + 1;
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[*k 
		    + i__ + i__ * a_dim1], &c__1);
	    i__2 = *n - i__ + 1;
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[*k + i__ + 
		    i__ * a_dim1], &c__1);
	    temp = a[*k + i__ + i__ * a_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &a[*k + i__ + i__ * a_dim1]);
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[*k + i__ + 1 + i__ 
		    * a_dim1], &c__1, &taul[i__]);
	    temp = a[*k + i__ + i__ * a_dim1];
	    a[*k + i__ + i__ * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + a_dim1], 
		    lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[*nb + 1]
		    , &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + a_dim1], lda, &c_b7, &q[i__ 
		    + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &q[
		    i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);

/*           Update XA with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[*k + i__ + (i__ + 
		    1) * a_dim1], lda, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[*k + i__ + a_dim1], lda, &c_b7, &a[*k 
		    + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[*k + 
		    i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &c_b7, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);

/*           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &a[*k + i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + 1 + a_dim1]
		    , lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, 
		    &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda);

/*           Update XG with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + i__ * xg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[*k + i__ + a_dim1], lda, &c_b7, &g[*k + i__ + 
		    g_dim1], ldg, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &c_b7, 
		    &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)12);

/*           Update XB with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__ + 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + i__ * xb_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ * b_dim1 + 1], &c__1, 
		    (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[*k + i__ + a_dim1], lda, &c_b7, &b[i__ * 
		    b_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[*k + i__ + 1 
		    + i__ * b_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &b[*
		    k + i__ + 1 + i__ * b_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ G(k+i,:); B(:,i)' ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &
		    c__1, &c__, &s);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + j * yg_dim1] = 0.;
/* L190: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
/* L200: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + j * ya_dim1] = 0.;
/* L210: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
/* L220: */
	    }

/*           Update XG with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__ + 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xb[(i__ 
		    + *nb) * xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ * b_dim1 + 1], &c__1);

	    a[*k + i__ + i__ * a_dim1] = temp;
	    q[i__ + i__ * q_dim1] = tauq;
	    csl[(i__ << 1) - 1] = c__;
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

	    alpha = q[i__ + (i__ + 1) * q_dim1];
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    *k + i__ + 1 + i__ * b_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[*k + 
		    i__ + 1 + i__ * b_dim1], &c__1);
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &b[*k + i__ + 1 + i__ * b_dim1]);
	    s = -s;
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &b[*k + i__ + 1 + i__ * b_dim1], &b[*k + i__ + 2 + 
		    i__ * b_dim1], &c__1, &taur[i__]);
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
	    b[*k + i__ + 1 + i__ * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + b_dim1]
		    , ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[*
		    nb + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb, &
		    c_b7, &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)
		    12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &b[
		    *k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);

/*           Update YQ with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b7, &
		    q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &q[
		    i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[*k + i__ 
		    + 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
/* L230: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
/* L240: */
	    }

/*           Update YB with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &b[*k + i__ + 1 + i__ * b_dim1], &c__1,
		     &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 2 + b_dim1]
		    , ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb);

/*           Update YQ with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[
		    i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + i__ * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[*k + i__ + 1 
		    + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b7, &
		    a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &a[(i__ + 
		    1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update YG with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + i__ * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &g[(*k + 
		    i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L250: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L260: */
	    }

/*           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(*k + i__ + 1) 
		    * g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &c_b9, 
		    &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);

/*           Update YG with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

	    b[*k + i__ + 1 + i__ * b_dim1] = temp;
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
	    csr[(i__ << 1) - 1] = c__;
	    csr[i__ * 2] = s;
/* L270: */
	}

    } else {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first columns of A and Q. See routine MB04TS. */

	    alpha = q[i__ + i__ * q_dim1];
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
	    q[i__ + i__ * q_dim1] = 1.;
	    i__2 = *n - i__ + 1;
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[*k 
		    + i__ + i__ * a_dim1], &c__1);
	    i__2 = *n - i__ + 1;
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[*k + i__ + 
		    i__ * a_dim1], &c__1);
	    temp = a[*k + i__ + i__ * a_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &a[*k + i__ + i__ * a_dim1]);
	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[*k + i__ + 1 + i__ 
		    * a_dim1], &c__1, &taul[i__]);
	    temp = a[*k + i__ + i__ * a_dim1];
	    a[*k + i__ + i__ * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + a_dim1], 
		    lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[*nb + 1]
		    , &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + a_dim1], lda, &c_b7, &q[i__ 
		    + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)9);

/*           Update XA with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[*k + i__ + (i__ + 
		    1) * a_dim1], lda, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[*k + i__ + a_dim1], lda, &c_b7, &a[*k 
		    + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[*k + 
		    i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &
		    c_b7, &a[*k + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);

/*           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &a[*k + i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + 1 + a_dim1]
		    , lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, 
		    &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)9);
	    i__2 = *n - i__;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda);

/*           Update XG with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[*k + i__ + a_dim1], lda, &c_b7, &g[*k + i__ + 
		    g_dim1], ldg, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &
		    c_b7, &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (
		    ftnlen)9);

/*           Update XB with first Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * xb_dim1 + 
		    1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + i__ * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ + b_dim1], ldb, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[*k + i__ + a_dim1], lda, &c_b7, &b[i__ + 
		    b_dim1], ldb, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[i__ + (*k + 
		    i__ + 1) * b_dim1], ldb, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &
		    b[i__ + (*k + i__ + 1) * b_dim1], ldb, (ftnlen)9);

/*           Apply rotation to [ G(k+i,:); B(i,:) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &
		    c__, &s);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + j * yg_dim1] = 0.;
/* L280: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
/* L290: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + j * ya_dim1] = 0.;
/* L300: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
/* L310: */
	    }

/*           Update XG with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

	    i__2 = *n - i__ + 1;
	    i__3 = *k + *n;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xb[(i__ + *nb) 
		    * xb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
	    i__2 = *k + *n;
	    d__1 = -taul[i__];
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ + b_dim1], ldb);

	    a[*k + i__ + i__ * a_dim1] = temp;
	    q[i__ + i__ * q_dim1] = tauq;
	    csl[(i__ << 1) - 1] = c__;
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

	    alpha = q[i__ + (i__ + 1) * q_dim1];
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    i__ + (*k + i__ + 1) * b_dim1], ldb);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[i__ + (
		    *k + i__ + 1) * b_dim1], ldb);
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
	    dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (*k + i__ + 1) * b_dim1]
		    );
	    s = -s;
	    i__2 = *n - i__;
	    dlarfg_(&i__2, &b[i__ + (*k + i__ + 1) * b_dim1], &b[i__ + (*k + 
		    i__ + 2) * b_dim1], ldb, &taur[i__]);
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
	    b[i__ + (*k + i__ + 1) * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], 
		    ldq, &c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)
		    12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &dwork[*nb + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[i__ + 
		    1 + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);

/*           Update YQ with first Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b7, &
		    q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12)
		    ;

/*           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ]. */

	    i__2 = *n - i__;
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, &c__, &s);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
/* L320: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
/* L330: */
	    }

/*           Update YB with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &b[i__ + (*k + i__ + 1) * b_dim1]
		    , ldb, &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1,
		     (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 2) * 
		    b_dim1 + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1);

/*           Update YQ with second Householder reflection. */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[
		    i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + i__ * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[*k + i__ + 1 
		    + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b7, &
		    a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &a[
		    (i__ + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update YG with first Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + i__ * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (
		    ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &g[
		    (*k + i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L340: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L350: */
	    }

/*           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ]. */

	    i__2 = *k + *n;
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(*k + i__ + 1) 
		    * g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &c_b9, 
		    &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);

/*           Update YG with second Householder reflection. */

	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *n - i__ - 1;
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *k + *n;
	    d__1 = -taur[i__];
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

	    b[i__ + (*k + i__ + 1) * b_dim1] = temp;
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
	    csr[(i__ << 1) - 1] = c__;
	    csr[i__ * 2] = s;
/* L360: */
	}
    }

    return 0;
/* *** Last line of MB03XU *** */
} /* mb03xu_ */

