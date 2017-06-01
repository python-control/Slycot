/* MB04PA.f -- translated by f2c (version 20100827).
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
static doublereal c_b10 = 1.;
static doublereal c_b12 = 0.;
static doublereal c_b586 = -1.;

/* Subroutine */ int mb04pa_(logical *lham, integer *n, integer *k, integer *
	nb, doublereal *a, integer *lda, doublereal *qg, integer *ldqg, 
	doublereal *xa, integer *ldxa, doublereal *xg, integer *ldxg, 
	doublereal *xq, integer *ldxq, doublereal *ya, integer *ldya, 
	doublereal *cs, doublereal *tau, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, xa_dim1, xa_offset, xg_dim1,
	     xg_offset, xq_dim1, xq_offset, ya_dim1, ya_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer nb1, nb2;
    static doublereal aki;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal tauq;
    extern /* Subroutine */ int mb01md_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


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

/*     To reduce a Hamiltonian like matrix */

/*                   [  A   G  ]           T          T */
/*              H =  [       T ] ,    G = G ,    Q = Q, */
/*                   [  Q  -A  ] */

/*     or a skew-Hamiltonian like matrix */

/*                   [  A   G  ]            T          T */
/*              W =  [       T ] ,    G = -G ,   Q = -Q, */
/*                   [  Q   A  ] */

/*     so that elements below the (k+1)-th subdiagonal in the first nb */
/*     columns of the (k+n)-by-n matrix A, and offdiagonal elements */
/*     in the first nb columns and rows of the n-by-n matrix Q are zero. */

/*     The reduction is performed by an orthogonal symplectic */
/*     transformation UU'*H*UU and matrices U, XA, XG, XQ, and YA are */
/*     returned so that */

/*                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'+ XG*U' ] */
/*         UU'*H*UU = [                                           ]. */
/*                    [ Qout + U*XQ'+ XQ*U'  -Aout'- XA*U'- U*YA' ] */

/*     Similarly, */

/*                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'- XG*U' ] */
/*         UU'*W*UU = [                                           ]. */
/*                    [ Qout + U*XQ'- XQ*U'   Aout'+ XA*U'+ U*YA' ] */

/*     This is an auxiliary routine called by MB04PB. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LHAM    LOGICAL */
/*             Specifies the type of matrix to be reduced: */
/*             = .FALSE. :  skew-Hamiltonian like W; */
/*             = .TRUE.  :  Hamiltonian like H. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     K       (input) INTEGER */
/*             The offset of the reduction. Elements below the (K+1)-th */
/*             subdiagonal in the first NB columns of A are reduced */
/*             to zero.  K >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns/rows to be reduced.  N > NB >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading (K+N)-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading (K+N)-by-N part of this array */
/*             contains the matrix Aout and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the reduction. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,K+N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N+K-by-N+1 part of this array must */
/*             contain in the bottom left part the lower triangular part */
/*             of the N-by-N matrix Q and in the remainder the upper */
/*             trapezoidal part of the last N columns of the N+K-by-N+K */
/*             matrix G. */
/*             On exit, the leading N+K-by-N+1 part of this array */
/*             contains parts of the matrices Q and G in the same fashion */
/*             as on entry only that the zero parts of Q contain */
/*             information about the elementary reflectors used to */
/*             compute the reduction. Note that if LHAM = .FALSE. then */
/*             the (K-1)-th and K-th subdiagonals are not referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG. LDQG >= MAX(1,N+K). */

/*     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XA. */

/*     LDXA    INTEGER */
/*             The leading dimension of the array XA.  LDXA >= MAX(1,N). */

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

/*     CS      (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2*NB elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations used */
/*             to compute the reduction. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB elements of this array contain the */
/*             scalar factors of some of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*NB) */

/*     METHOD */

/*     For details regarding the representation of the orthogonal */
/*     symplectic matrix UU within the arrays A, QG, CS, TAU see the */
/*     description of MB04PU. */

/*     The contents of A and QG on exit are illustrated by the following */
/*     example with n = 5, k = 2 and nb = 2: */

/*           ( a  r  r  a  a  )         ( g  g  r  r  g  g  ) */
/*           ( a  r  r  a  a  )         ( g  g  r  r  g  g  ) */
/*           ( a  r  r  a  a  )         ( q  g  r  r  g  g  ) */
/*       A = ( r  r  r  r  r  ),   QG = ( t  r  r  r  r  r  ), */
/*           ( u2 r  r  r  r  )         ( u1 t  r  r  r  r  ) */
/*           ( u2 u2 r  a  a  )         ( u1 u1 r  q  g  g  ) */
/*           ( u2 u2 r  a  a  )         ( u1 u1 r  q  q  g  ) */

/*     where a, g and q denote elements of the original matrices, r */
/*     denotes a modified element, t denotes a scalar factor of an */
/*     applied elementary reflector and ui denote elements of the */
/*     matrix U. */

/*     REFERENCES */

/*     [1] C. F. VAN LOAN: */
/*         A symplectic method for approximating all the eigenvalues of */
/*         a Hamiltonian matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     [2] D. KRESSNER: */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner (Technical Univ. Berlin, Germany) and */
/*     P. Benner (Technical Univ. Chemnitz, Germany), December 2003. */

/*     REVISIONS */

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DLAPVL). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix, */
/*     skew-Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    xa_dim1 = *ldxa;
    xa_offset = 1 + xa_dim1;
    xa -= xa_offset;
    xg_dim1 = *ldxg;
    xg_offset = 1 + xg_dim1;
    xg -= xg_offset;
    xq_dim1 = *ldxq;
    xq_offset = 1 + xq_dim1;
    xq -= xq_offset;
    ya_dim1 = *ldya;
    ya_offset = 1 + ya_dim1;
    ya -= ya_offset;
    --cs;
    --tau;
    --dwork;

    /* Function Body */
    if (*n + *k <= 0) {
	dwork[1] = 1.;
	return 0;
    }

    nb1 = *nb + 1;
    nb2 = *nb + nb1;

    if (*lham) {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform i-th columns of A and Q. See routine MB04PU. */

	    alpha = qg[*k + i__ + 1 + i__ * qg_dim1];
	    i__2 = *n - i__;
/* Computing MIN */
	    i__3 = i__ + 2;
	    dlarfg_(&i__2, &alpha, &qg[*k + min(i__3,*n) + i__ * qg_dim1], &
		    c__1, &tauq);
	    qg[*k + i__ + 1 + i__ * qg_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &a[*k + i__ + 1 + i__ * a_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &a[
		    *k + i__ + 1 + i__ * a_dim1], &c__1);
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
	    dlartg_(&aki, &alpha, &c__, &s, &a[*k + i__ + 1 + i__ * a_dim1]);
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
	    i__2 = *n - i__;
/* Computing MIN */
	    i__3 = i__ + 2;
	    dlarfg_(&i__2, &aki, &a[*k + min(i__3,*n) + i__ * a_dim1], &c__1, 
		    &tau[i__]);
	    a[*k + i__ + 1 + i__ * a_dim1] = 1.;

/*           Update XA with first Householder reflection. */

/*           xa = H(1:n,1:n)'*u1 */
	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + i__ * xa_dim1], &c__1, (
		    ftnlen)9);
/*           w1 = U1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*           xa = xa + XA1*w1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + 
		    xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           w2 = U2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*           xa = xa + XA2*w2 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + nb1 *
		     ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           xa = -tauq*xa */
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update YA with first Householder reflection. */

/*           ya = H(1:n,1:n)*u1 */
	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &ya[i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
/*           temp = XA1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           temp = XA2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ + 
		    1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           ya = ya + YA1*w1 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[ya_offset], ldya,
		     &dwork[1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           ya = ya + YA2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           ya = -tauq*ya */
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);
/*           temp = -tauq*ya'*u1 */
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);
/*           ya = ya + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ + 1 + xa_dim1], ldxa, &c_b10, &a[*
		    k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b10, 
		    &a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)'; */
	    i__2 = *n + *k;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)'; */
	    i__2 = *n + *k;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &a[(i__ 
		    + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update (i+1)-th row of A. */

	    if (*n > i__ + 1) {
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &ya[*k + i__ + 1 + ya_dim1], ldya, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &ya[*k + i__ + 1 + nb1 * ya_dim1], 
			ldya, &c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], 
			lda, (ftnlen)12);
	    }

/*           Annihilate updated parts in YA. */

	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + 1 + j * ya_dim1] = 0.;
/* L10: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + 1 + (*nb + j) * ya_dim1] = 0.;
/* L20: */
	    }

/*           Update XQ with first Householder reflection. */

/*           xq = Q*u1 */
	    i__2 = *n - i__;
	    dsymv_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)5);
/*           xq = xq + XQ1*w1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + 
		    xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = xq + XQ2*w2 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = -tauq*xq */
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
	    i__2 = *n - i__;
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           xq = xq + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

/*           Q(:,i+1) = Q(:,i+1) + U1 * XQ1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ + 1 + xq_dim1], ldxq, &c_b10, &qg[
		    *k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12);
/*           Q(:,i+1) = Q(:,i+1) + U2 * XQ2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b10, 
		    &qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)
		    12);
/*           Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[*k 
		    + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12);
/*           Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &
		    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12)
		    ;

/*           Update XG with first Householder reflection. */

/*           xg = G*u1 */
	    i__2 = *k + i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xg[i__ * xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    dsymv_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1, (ftnlen)
		    5);
/*           xg = xg + XG1*w1 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[xg_offset], ldxg,
		     &dwork[1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           xg = xg + XG2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           temp = XG1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           temp = XG2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + nb1 *
		     xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ + 
		    1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           xg = -tauq*xg */
	    i__2 = *n + *k;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
	    i__2 = *n - i__;
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);
/*           xg = xg + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
	    i__2 = *k + i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[(i__ + 2) *
		     qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
	    i__2 = *k + i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10,
		     &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (ftnlen)
		    12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    nb1 * xg_dim1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &
		    c_b10, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (
		    ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b10,
		     &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (ftnlen)
		    12);
/*           G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b10, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (
		    ftnlen)12);

/*           Annihilate updated parts in XG. */

	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L30: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L40: */
	    }

/*           Apply orthogonal symplectic Givens rotation. */

	    i__2 = *k + i__;
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		    qg_dim1 + 1], &c__1, &c__, &s);
	    if (*n > i__ + 1) {
		i__2 = *n - i__ - 1;
		drot_(&i__2, &a[*k + i__ + 2 + (i__ + 1) * a_dim1], &c__1, &
			qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, &c__, &
			s);
		i__2 = *n - i__ - 1;
		drot_(&i__2, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[*
			k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1, &c__, &s);
	    }
	    temp = a[*k + i__ + 1 + (i__ + 1) * a_dim1];
	    ttemp = qg[*k + i__ + 1 + (i__ + 2) * qg_dim1];
	    a[*k + i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[*k + 
		    i__ + 1 + (i__ + 1) * qg_dim1];
	    qg[*k + i__ + 1 + (i__ + 2) * qg_dim1] = c__ * ttemp - s * temp;
	    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1] = -s * temp + c__ * qg[*k 
		    + i__ + 1 + (i__ + 1) * qg_dim1];
	    ttemp = -s * ttemp - c__ * temp;
	    temp = a[*k + i__ + 1 + (i__ + 1) * a_dim1];
	    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1] = c__ * qg[*k + i__ + 1 + (
		    i__ + 1) * qg_dim1] + s * ttemp;
	    a[*k + i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[*k + 
		    i__ + 1 + (i__ + 2) * qg_dim1];
	    qg[*k + i__ + 1 + (i__ + 2) * qg_dim1] = -s * temp + c__ * qg[*k 
		    + i__ + 1 + (i__ + 2) * qg_dim1];
	    cs[(i__ << 1) - 1] = c__;
	    cs[i__ * 2] = s;
	    qg[*k + i__ + 1 + i__ * qg_dim1] = tauq;

/*           Update XA with second Householder reflection. */

/*           xa = H(1:n,1:n)'*u2 */
	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1,
		     (ftnlen)9);
	    if (*n > i__ + 1) {
/*              w1 = U1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 + 
			qg_dim1], ldqg, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*              xa = xa + XA1*w1 */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 2 
			+ (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              w2 = U2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 + 
			a_dim1], lda, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, 
			&c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*              xa = xa + XA2*w2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[
			i__ + 2 + (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              temp = YA1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &ya[*k + i__ + 2 + 
			ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xa = xa + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xa[(*nb + i__) * xa_dim1 + 1], &
			c__1, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, (ftnlen)12);
/*              temp = YA2'*u1 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 2 + 
			nb1 * ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 (ftnlen)9);
/*              xa = xa + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &c__1, 
			(ftnlen)12);
	    }
/*           xa = -tau*xa */
	    i__2 = *n - i__;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1);

/*           Update YA with second Householder reflection. */

/*           ya = H(1:n,1:n)*u2 */
	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, 
		    &c_b12, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, (ftnlen)12)
		    ;
	    if (*n > i__ + 1) {
/*              temp = XA1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k 
			+ i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XA2'*u1 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + nb1 * 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + 
			i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12);
	    }
/*           ya = ya + YA1*w1 */
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           ya = ya + YA2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[(*nb + i__) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
/*           ya = -tau*ya */
	    i__2 = *k + *n;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &ya[(*nb + i__) * ya_dim1 + 1], &c__1);
/*           temp = -tau*ya'*u2 */
	    i__2 = *n - i__;
	    temp = -tau[i__] * ddot_(&i__2, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);
/*           ya = ya + temp*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &ya[*
		    k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           H(1:n,i+1) = H(1:n,i+1) + ya */
	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);
/*           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + (i__ + 1) * 
		    a_dim1], &c__1);

/*           Update (i+1)-th row of A. */

	    if (*n > i__ + 1) {
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)'; */
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda);
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + 
			(i__ + 2) * a_dim1], lda);
	    }

/*           Annihilate updated parts in YA. */

	    ya[*k + i__ + 1 + (*nb + i__) * ya_dim1] = 0.;

/*           Update XQ with second Householder reflection. */

/*           xq = Q*u2 */
	    i__2 = *n - i__;
	    dsymv_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1, (
		    ftnlen)5);
	    if (*n > i__ + 1) {
/*              xq = xq + XQ1*w1 */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 2 
			+ (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              xq = xq + XQ2*w2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[
			i__ + 2 + (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              temp = XQ1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[(*nb + i__) * xq_dim1 + 1], &
			c__1, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, (ftnlen)12);
/*              temp = XQ2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + nb1 * 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[(*nb + i__) * xq_dim1 + 1], &c__1,
			 &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &c__1, 
			(ftnlen)12);
	    }
/*           xq = -tauq*xq */
	    i__2 = *n - i__;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u2 */
	    i__2 = *n - i__;
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &
		    c__1);
/*           xq = xq + temp*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xq[
		    i__ + 1 + (*nb + i__) * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b10, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1,
		     &qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1);
/*           H(1:n,n+i+1) = H(1:n,n+i+1) + U * XQ(i+1,:)'; */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], &c__1);

/*           Update XG with second Householder reflection. */

/*           xg = G*u2 */
	    i__2 = *k + i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    dsymv_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1, (
		    ftnlen)5);
/*           xg = xg + XG1*w1 */
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           xg = xg + XG2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[(*nb + i__) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    if (*n > i__ + 1) {
/*              temp = XG1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 2 + 
			xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k 
			+ i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XG2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 2 + 
			nb1 * xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + 
			i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12);
	    }
/*           xg = -tauq*xg */
	    i__2 = *n + *k;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xg[(*nb + i__) * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xg'*u1 */
	    i__2 = *n - i__;
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1],
		     &c__1);
/*           xg = xg + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xg[*
		    k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

	    i__2 = *k + i__;
	    daxpy_(&i__2, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &c_b10, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &
		    c__1, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &a[*k + 
		    i__ + 1 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 
		    2) * qg_dim1], ldqg);

/*           Annihilate updated parts in XG. */

	    xg[*k + i__ + 1 + (*nb + i__) * xg_dim1] = 0.;

	    a[*k + i__ + 1 + i__ * a_dim1] = aki;
/* L50: */
	}
    } else {
	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform i-th columns of A and Q. */

	    alpha = qg[*k + i__ + 1 + i__ * qg_dim1];
	    i__2 = *n - i__;
/* Computing MIN */
	    i__3 = i__ + 2;
	    dlarfg_(&i__2, &alpha, &qg[*k + min(i__3,*n) + i__ * qg_dim1], &
		    c__1, &tauq);
	    qg[*k + i__ + 1 + i__ * qg_dim1] = 1.;
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &a[*k + i__ + 1 + i__ * a_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &a[
		    *k + i__ + 1 + i__ * a_dim1], &c__1);
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
	    dlartg_(&aki, &alpha, &c__, &s, &a[*k + i__ + 1 + i__ * a_dim1]);
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
	    i__2 = *n - i__;
/* Computing MIN */
	    i__3 = i__ + 2;
	    dlarfg_(&i__2, &aki, &a[*k + min(i__3,*n) + i__ * a_dim1], &c__1, 
		    &tau[i__]);
	    a[*k + i__ + 1 + i__ * a_dim1] = 1.;

/*           Update XA with first Householder reflection. */

/*           xa = H(1:n,1:n)'*u1 */
	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + i__ * xa_dim1], &c__1, (
		    ftnlen)9);
/*           w1 = U1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*           xa = xa + XA1*w1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + 
		    xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           w2 = U2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*           xa = xa + XA2*w2 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + nb1 *
		     ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           xa = -tauq*xa */
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update YA with first Householder reflection. */

/*           ya = H(1:n,1:n)*u1 */
	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &ya[i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
/*           temp = XA1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           temp = XA2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ + 
		    1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           ya = ya + YA1*w1 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[ya_offset], ldya,
		     &dwork[1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           ya = ya + YA2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           ya = -tauq*ya */
	    i__2 = *k + *n;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);
/*           temp = -tauq*ya'*u1 */
	    i__2 = *n - i__;
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);
/*           ya = ya + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)'; */
	    i__2 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ + 1 + xa_dim1], ldxa, &c_b10, &a[*
		    k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)'; */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b10, 
		    &a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)'; */
	    i__2 = *n + *k;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)'; */
	    i__2 = *n + *k;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &a[(i__ 
		    + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update (i+1)-th row of A. */

	    if (*n > i__ + 1) {
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &ya[*k + i__ + 1 + ya_dim1], ldya, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &ya[*k + i__ + 1 + nb1 * ya_dim1], 
			ldya, &c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], 
			lda, (ftnlen)12);
	    }

/*           Annihilate updated parts in YA. */

	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + 1 + j * ya_dim1] = 0.;
/* L60: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		ya[*k + i__ + 1 + (*nb + j) * ya_dim1] = 0.;
/* L70: */
	    }

/*           Update XQ with first Householder reflection. */

/*           xq = Q*u1 */
	    i__2 = *n - i__;
	    mb01md_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)5);
/*           xq = xq + XQ1*w1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + 
		    xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = xq + XQ2*w2 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq - U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq - U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = -tauq*xq */
	    i__2 = *n - i__;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
	    i__2 = *n - i__;
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           xq = xq + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

	    if (*n > i__ + 1) {
/*              Q(:,i+1) = Q(:,i+1) - U1 * XQ1(i+1,:)'; */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[i__ + 1 + xq_dim1], ldxq, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) - U2 * XQ2(i+1,:)'; */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)'; */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)'; */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &
			c__1, (ftnlen)12);
	    }

/*           Update XG with first Householder reflection. */

/*           xg = G*u1 */
	    i__2 = *k + i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xg[i__ * xg_dim1 + 1], &c__1, (ftnlen)12);
	    i__2 = *n - i__;
	    mb01md_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1, (ftnlen)
		    5);
/*           xg = xg + XG1*w1 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[xg_offset], ldxg,
		     &dwork[1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           xg = xg + XG2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           temp = XG1'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg - U1*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           temp = XG2'*u1 */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + nb1 *
		     xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg - U2*temp */
	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ + 
		    1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           xg = -tauq*xg */
	    i__2 = *n + *k;
	    d__1 = -tauq;
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
	    i__2 = *n - i__;
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);
/*           xg = xg + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
	    i__2 = *k + i__;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[(i__ + 2) *
		     qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
	    i__2 = *k + i__;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1, (ftnlen)12);
	    if (*n > i__ + 1) {
/*              G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b586, &xg[*k + i__ + 2 
			+ xg_dim1], ldxg, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, 
			(ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b586, &xg[*k + i__ + 
			2 + nb1 * xg_dim1], ldxg, &a[*k + i__ + 1 + a_dim1], 
			lda, &c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], 
			ldqg, (ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)'; */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xg[*k + i__ + 1 + xg_dim1], ldxg, &
			c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, 
			(ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)'; */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], 
			ldxg, &c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1],
			 ldqg, (ftnlen)12);
	    }

/*           Annihilate updated parts in XG. */

	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
/* L80: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
/* L90: */
	    }

/*           Apply orthogonal symplectic Givens rotation. */

	    i__2 = *k + i__;
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		    qg_dim1 + 1], &c__1, &c__, &s);
	    if (*n > i__ + 1) {
		i__2 = *n - i__ - 1;
		d__1 = -s;
		drot_(&i__2, &a[*k + i__ + 2 + (i__ + 1) * a_dim1], &c__1, &
			qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, &c__, &
			d__1);
		i__2 = *n - i__ - 1;
		d__1 = -s;
		drot_(&i__2, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[*
			k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1, &c__, &
			d__1);
	    }
	    cs[(i__ << 1) - 1] = c__;
	    cs[i__ * 2] = s;
	    qg[*k + i__ + 1 + i__ * qg_dim1] = tauq;

/*           Update XA with second Householder reflection. */

/*           xa = H(1:n,1:n)'*u2 */
	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1,
		     (ftnlen)9);
	    if (*n > i__ + 1) {
/*              w1 = U1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 + 
			qg_dim1], ldqg, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*              xa = xa + XA1*w1 */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 2 
			+ (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              w2 = U2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 + 
			a_dim1], lda, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, 
			&c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*              xa = xa + XA2*w2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[
			i__ + 2 + (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              temp = YA1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &ya[*k + i__ + 2 + 
			ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xa = xa + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xa[(*nb + i__) * xa_dim1 + 1], &
			c__1, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, (ftnlen)12);
/*              temp = YA2'*u1 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 2 + 
			nb1 * ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 (ftnlen)9);
/*              xa = xa + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &c__1, 
			(ftnlen)12);
	    }
/*           xa = -tau*xa */
	    i__2 = *n - i__;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1);

/*           Update YA with second Householder reflection. */

/*           ya = H(1:n,1:n)*u2 */
	    i__2 = *k + *n;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, 
		    &c_b12, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, (ftnlen)12)
		    ;
	    if (*n > i__ + 1) {
/*              temp = XA1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k 
			+ i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XA2'*u1 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + nb1 * 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + 
			i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12);
	    }
/*           ya = ya + YA1*w1 */
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           ya = ya + YA2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[(*nb + i__) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
/*           ya = -tau*ya */
	    i__2 = *k + *n;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &ya[(*nb + i__) * ya_dim1 + 1], &c__1);
/*           temp = -tau*ya'*u2 */
	    i__2 = *n - i__;
	    temp = -tau[i__] * ddot_(&i__2, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);
/*           ya = ya + temp*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &ya[*
		    k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           H(1:n,i+1) = H(1:n,i+1) + ya */
	    i__2 = *k + *n;
	    daxpy_(&i__2, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);
/*           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + (i__ + 1) * 
		    a_dim1], &c__1);

/*           Update (i+1)-th row of A. */

	    if (*n > i__ + 1) {
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)'; */
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda);
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)' */
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + 
			(i__ + 2) * a_dim1], lda);
	    }

/*           Annihilate updated parts in YA. */

	    ya[*k + i__ + 1 + (*nb + i__) * ya_dim1] = 0.;

/*           Update XQ with second Householder reflection. */

/*           xq = Q*u2 */
	    i__2 = *n - i__;
	    mb01md_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1, (
		    ftnlen)5);
	    if (*n > i__ + 1) {
/*              xq = xq + XQ1*w1 */
		i__2 = *n - i__ - 1;
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 2 
			+ (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              xq = xq + XQ2*w2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[
			i__ + 2 + (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              temp = XQ1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq - U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[(*nb + i__) * xq_dim1 + 1], &
			c__1, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, (ftnlen)12);
/*              temp = XQ2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + nb1 * 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq - U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[(*nb + i__) * xq_dim1 + 1], &c__1,
			 &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &c__1, 
			(ftnlen)12);
	    }
/*           xq = -tauq*xq */
	    i__2 = *n - i__;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u2 */
	    i__2 = *n - i__;
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &
		    c__1);
/*           xq = xq + temp*u2 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xq[
		    i__ + 1 + (*nb + i__) * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

	    if (*n > i__ + 1) {
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1);
/*              H(1:n,n+i+1) = H(1:n,n+i+1) - U * XQ(i+1,:)'; */
		i__2 = *n - i__ - 1;
		d__1 = -xq[i__ + 1 + (*nb + i__) * xq_dim1];
		daxpy_(&i__2, &d__1, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, &
			qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1);
	    }

/*           Update XG with second Householder reflection. */

/*           xg = G*u2 */
	    i__2 = *k + i__;
	    i__3 = *n - i__;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
	    i__2 = *n - i__;
	    mb01md_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1, (
		    ftnlen)5);
/*           xg = xg + XG1*w1 */
	    i__2 = *k + *n;
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           xg = xg + XG2*w2 */
	    i__2 = *k + *n;
	    i__3 = i__ - 1;
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[(*nb + i__) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
	    if (*n > i__ + 1) {
/*              temp = XG1'*u2 */
		i__2 = *n - i__ - 1;
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 2 + 
			xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg - U1*temp */
		i__2 = *n - i__ - 1;
		dgemv_("No Transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k 
			+ i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XG2'*u2 */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 2 + 
			nb1 * xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg - U2*temp */
		i__2 = *n - i__ - 1;
		i__3 = i__ - 1;
		dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + 
			i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12);
	    }
/*           xg = -tauq*xg */
	    i__2 = *n + *k;
	    d__1 = -tau[i__];
	    dscal_(&i__2, &d__1, &xg[(*nb + i__) * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xg'*u1 */
	    i__2 = *n - i__;
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1],
		     &c__1);
/*           xg = xg + temp*u1 */
	    i__2 = *n - i__;
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xg[*
		    k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

	    i__2 = *k + i__;
	    daxpy_(&i__2, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1);
	    if (*n > i__ + 1) {
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &c_b586, &xg[*k + i__ + 2 + (*nb + i__) * 
			xg_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 3) * 
			qg_dim1], ldqg);
		i__2 = *n - i__ - 1;
		daxpy_(&i__2, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 
			+ (i__ + 3) * qg_dim1], ldqg);
	    }

/*           Annihilate updated parts in XG. */

	    xg[*k + i__ + 1 + (*nb + i__) * xg_dim1] = 0.;

	    a[*k + i__ + 1 + i__ * a_dim1] = aki;
/* L100: */
	}
    }

    return 0;
/* *** Last line of MB04PA *** */
} /* mb04pa_ */

