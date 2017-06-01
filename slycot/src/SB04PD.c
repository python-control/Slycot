/* SB04PD.f -- translated by f2c (version 20100827).
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

static doublereal c_b22 = 1.;
static doublereal c_b23 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb04pd_(char *dico, char *facta, char *factb, char *
	trana, char *tranb, integer *isgn, integer *m, integer *n, doublereal 
	*a, integer *lda, doublereal *u, integer *ldu, doublereal *b, integer 
	*ldb, doublereal *v, integer *ldv, doublereal *c__, integer *ldc, 
	doublereal *scale, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen facta_len, ftnlen factb_len, ftnlen trana_len,
	 ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ia, ib, bl, sdim, ierr;
    static logical cont;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), sb04py_(char *,
	     char *, integer *, integer *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static logical bwork[1];
    static integer jwork;
    static logical blas3a, blas3b, nofaca, nofacb, blocka, blockb;
    static integer chunka, chunkb;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern logical select_();
    static integer availw;
    static logical schura, schurb, notrna, notrnb;
    static integer minwrk, maxwrk;
    extern /* Subroutine */ int dtrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen);


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

/*     To solve for X either the real continuous-time Sylvester equation */

/*        op(A)*X + ISGN*X*op(B) = scale*C,                           (1) */

/*     or the real discrete-time Sylvester equation */

/*        op(A)*X*op(B) + ISGN*X = scale*C,                           (2) */

/*     where op(M) = M or M**T, and ISGN = 1 or -1. A is M-by-M and */
/*     B is N-by-N; the right hand side C and the solution X are M-by-N; */
/*     and scale is an output scale factor, set less than or equal to 1 */
/*     to avoid overflow in X. The solution matrix X is overwritten */
/*     onto C. */

/*     If A and/or B are not (upper) quasi-triangular, that is, block */
/*     upper triangular with 1-by-1 and 2-by-2 diagonal blocks, they are */
/*     reduced to Schur canonical form, that is, quasi-triangular with */
/*     each 2-by-2 diagonal block having its diagonal elements equal and */
/*     its off-diagonal elements of opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the equation from which X is to be determined */
/*             as follows: */
/*             = 'C':  Equation (1), continuous-time case; */
/*             = 'D':  Equation (2), discrete-time case. */

/*     FACTA   CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix A is supplied on entry, as follows: */
/*             = 'F':  On entry, A and U contain the factors from the */
/*                     real Schur factorization of the matrix A; */
/*             = 'N':  The Schur factorization of A will be computed */
/*                     and the factors will be stored in A and U; */
/*             = 'S':  The matrix A is quasi-triangular (or Schur). */

/*     FACTB   CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix B is supplied on entry, as follows: */
/*             = 'F':  On entry, B and V contain the factors from the */
/*                     real Schur factorization of the matrix B; */
/*             = 'N':  The Schur factorization of B will be computed */
/*                     and the factors will be stored in B and V; */
/*             = 'S':  The matrix B is quasi-triangular (or Schur). */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op(B) to be used, as follows: */
/*             = 'N':  op(B) = B    (No transpose); */
/*             = 'T':  op(B) = B**T (Transpose); */
/*             = 'C':  op(B) = B**T (Conjugate transpose = Transpose). */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A, and the number of rows in the */
/*             matrices X and C.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix B, and the number of columns in */
/*             the matrices X and C.  N >= 0. */

/*     A       (input or input/output) DOUBLE PRECISION array, */
/*             dimension (LDA,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the matrix A. If FACTA = 'S', then A contains */
/*             a quasi-triangular matrix, and if FACTA = 'F', then A */
/*             is in Schur canonical form; the elements below the upper */
/*             Hessenberg part of the array A are not referenced. */
/*             On exit, if FACTA = 'N', and INFO = 0 or INFO >= M+1, the */
/*             leading M-by-M upper Hessenberg part of this array */
/*             contains the upper quasi-triangular matrix in Schur */
/*             canonical form from the Schur factorization of A. The */
/*             contents of array A is not modified if FACTA = 'F' or 'S'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     U       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDU,M) */
/*             If FACTA = 'F', then U is an input argument and on entry */
/*             the leading M-by-M part of this array must contain the */
/*             orthogonal matrix U of the real Schur factorization of A. */
/*             If FACTA = 'N', then U is an output argument and on exit, */
/*             if INFO = 0 or INFO >= M+1, it contains the orthogonal */
/*             M-by-M matrix from the real Schur factorization of A. */
/*             If FACTA = 'S', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= MAX(1,M), if FACTA = 'F' or 'N'; */
/*             LDU >= 1,        if FACTA = 'S'. */

/*     B       (input or input/output) DOUBLE PRECISION array, */
/*             dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix B. If FACTB = 'S', then B contains */
/*             a quasi-triangular matrix, and if FACTB = 'F', then B */
/*             is in Schur canonical form; the elements below the upper */
/*             Hessenberg part of the array B are not referenced. */
/*             On exit, if FACTB = 'N', and INFO = 0 or INFO = M+N+1, */
/*             the leading N-by-N upper Hessenberg part of this array */
/*             contains the upper quasi-triangular matrix in Schur */
/*             canonical form from the Schur factorization of B. The */
/*             contents of array B is not modified if FACTB = 'F' or 'S'. */

/*     LDB     (input) INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     V       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDV,N) */
/*             If FACTB = 'F', then V is an input argument and on entry */
/*             the leading N-by-N part of this array must contain the */
/*             orthogonal matrix V of the real Schur factorization of B. */
/*             If FACTB = 'N', then V is an output argument and on exit, */
/*             if INFO = 0 or INFO = M+N+1, it contains the orthogonal */
/*             N-by-N matrix from the real Schur factorization of B. */
/*             If FACTB = 'S', the array V is not referenced. */

/*     LDV     INTEGER */
/*             The leading dimension of array V. */
/*             LDV >= MAX(1,N), if FACTB = 'F' or 'N'; */
/*             LDV >= 1,        if FACTB = 'S'. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix C. */
/*             On exit, if INFO = 0 or INFO = M+N+1, the leading M-by-N */
/*             part of this array contains the solution matrix X. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or M+N+1, then: DWORK(1) returns the */
/*             optimal value of LDWORK; if FACTA = 'N', DWORK(1+i) and */
/*             DWORK(1+M+i), i = 1,...,M, contain the real and imaginary */
/*             parts, respectively, of the eigenvalues of A; and, if */
/*             FACTB = 'N', DWORK(1+f+j) and DWORK(1+f+N+j), j = 1,...,N, */
/*             with f = 2*M if FACTA = 'N', and f = 0, otherwise, contain */
/*             the real and imaginary parts, respectively, of the */
/*             eigenvalues of B. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, a+MAX( c, b+d, b+e ) ), */
/*             where a = 1+2*M, if FACTA =  'N', */
/*                   a = 0,     if FACTA <> 'N', */
/*                   b = 2*N,   if FACTB =  'N', FACTA =  'N', */
/*                   b = 1+2*N, if FACTB =  'N', FACTA <> 'N', */
/*                   b = 0,     if FACTB <> 'N', */
/*                   c = 3*M,   if FACTA =  'N', */
/*                   c = M,     if FACTA =  'F', */
/*                   c = 0,     if FACTA =  'S', */
/*                   d = 3*N,   if FACTB =  'N', */
/*                   d = N,     if FACTB =  'F', */
/*                   d = 0,     if FACTB =  'S', */
/*                   e = M,     if DICO  =  'C', FACTA <> 'S', */
/*                   e = 0,     if DICO  =  'C', FACTA =  'S', */
/*                   e = 2*M,   if DICO  =  'D'. */
/*             An upper bound is */
/*             LDWORK = 1+2*M+MAX( 3*M, 5*N, 2*N+2*M ). */
/*             For good performance, LDWORK should be larger, e.g., */
/*             LDWORK = 1+2*M+MAX( 3*M, 5*N, 2*N+2*M, 2*N+M*N ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if INFO = i, i = 1,...,M, the QR algorithm failed */
/*                   to compute all the eigenvalues of the matrix A */
/*                   (see LAPACK Library routine DGEES); the elements */
/*                   2+i:1+M and 2+i+M:1+2*M of DWORK contain the real */
/*                   and imaginary parts, respectively, of the */
/*                   eigenvalues of A which have converged, and the */
/*                   array A contains the partially converged Schur form; */
/*             = M+j:  if INFO = M+j, j = 1,...,N, the QR algorithm */
/*                   failed to compute all the eigenvalues of the matrix */
/*                   B (see LAPACK Library routine DGEES); the elements */
/*                   2+f+j:1+f+N and 2+f+j+N:1+f+2*N of DWORK contain the */
/*                   real and imaginary parts, respectively, of the */
/*                   eigenvalues of B which have converged, and the */
/*                   array B contains the partially converged Schur form; */
/*                   as defined for the parameter DWORK, */
/*                   f = 2*M, if FACTA =  'N', */
/*                   f = 0,   if FACTA <> 'N'; */
/*             = M+N+1:  if DICO = 'C', and the matrices A and -ISGN*B */
/*                   have common or very close eigenvalues, or */
/*                   if DICO = 'D', and the matrices A and -ISGN*B have */
/*                   almost reciprocal eigenvalues (that is, if lambda(i) */
/*                   and mu(j) are eigenvalues of A and -ISGN*B, then */
/*                   lambda(i) = 1/mu(j) for some i and j); */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrices A and B are unchanged). */

/*     METHOD */

/*     An extension and refinement of the algorithms in [1,2] is used. */
/*     If the matrices A and/or B are not quasi-triangular (see PURPOSE), */
/*     they are reduced to Schur canonical form */

/*        A = U*S*U',  B = V*T*V', */

/*     where U, V are orthogonal, and S, T are block upper triangular */
/*     with 1-by-1 and 2-by-2 blocks on their diagonal. The right hand */
/*     side matrix C is updated accordingly, */

/*        C = U'*C*V; */

/*     then, the solution matrix X of the "reduced" Sylvester equation */
/*     (with A and B in (1) or (2) replaced by S and T, respectively), */
/*     is computed column-wise via a back substitution scheme. A set of */
/*     equivalent linear algebraic systems of equations of order at most */
/*     four are formed and solved using Gaussian elimination with */
/*     complete pivoting. Finally, the solution X of the original */
/*     equation is obtained from the updating formula */

/*        X = U*X*V'. */

/*     If A and/or B are already quasi-triangular (or in Schur form), the */
/*     initial factorizations and the corresponding updating steps are */
/*     omitted. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is stable and reliable, since orthogonal */
/*     transformations and Gaussian elimination with complete pivoting */
/*     are used. If INFO = M+N+1, the Sylvester equation is numerically */
/*     singular. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, April 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix algebra, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    cont = lsame_(dico, "C", (ftnlen)1, (ftnlen)1);
    nofaca = lsame_(facta, "N", (ftnlen)1, (ftnlen)1);
    nofacb = lsame_(factb, "N", (ftnlen)1, (ftnlen)1);
    schura = lsame_(facta, "S", (ftnlen)1, (ftnlen)1);
    schurb = lsame_(factb, "S", (ftnlen)1, (ftnlen)1);
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! cont && ! lsame_(dico, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! nofaca && ! lsame_(facta, "F", (ftnlen)1, (ftnlen)1) && ! 
	    schura) {
	*info = -2;
    } else if (! nofacb && ! lsame_(factb, "F", (ftnlen)1, (ftnlen)1) && ! 
	    schurb) {
	*info = -3;
    } else if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -4;
    } else if (! notrnb && ! lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -5;
    } else if (*isgn != 1 && *isgn != -1) {
	*info = -6;
    } else if (*m < 0) {
	*info = -7;
    } else if (*n < 0) {
	*info = -8;
    } else if (*lda < max(1,*m)) {
	*info = -10;
    } else if (*ldu < 1 || ! schura && *ldu < *m) {
	*info = -12;
    } else if (*ldb < max(1,*n)) {
	*info = -14;
    } else if (*ldv < 1 || ! schurb && *ldv < *n) {
	*info = -16;
    } else if (*ldc < max(1,*m)) {
	*info = -18;
    } else {
	if (nofaca) {
	    ia = (*m << 1) + 1;
	    minwrk = *m * 3;
	} else {
	    ia = 0;
	}
	if (schura) {
	    minwrk = 0;
	} else if (! nofaca) {
	    minwrk = *m;
	}
	ib = 0;
	if (nofacb) {
	    ib = *n << 1;
	    if (! nofaca) {
		++ib;
	    }
/* Computing MAX */
	    i__1 = minwrk, i__2 = ib + *n * 3;
	    minwrk = max(i__1,i__2);
	} else if (! schurb) {
	    minwrk = max(minwrk,*n);
	}
	if (cont) {
	    if (! schura) {
/* Computing MAX */
		i__1 = minwrk, i__2 = ib + *m;
		minwrk = max(i__1,i__2);
	    }
	} else {
/* Computing MAX */
	    i__1 = minwrk, i__2 = ib + (*m << 1);
	    minwrk = max(i__1,i__2);
	}
/* Computing MAX */
	i__1 = 1, i__2 = ia + minwrk;
	minwrk = max(i__1,i__2);
	if (*ldwork < minwrk) {
	    *info = -21;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB04PD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	*scale = 1.;
	dwork[1] = 1.;
	return 0;
    }
    maxwrk = minwrk;

    if (nofaca) {

/*        Compute the Schur factorization of A. */
/*        Workspace:  need   1+5*M; */
/*                    prefer larger. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

	jwork = (*m << 1) + 2;
	ia = jwork;
	availw = *ldwork - jwork + 1;
	dgees_("Vectors", "Not ordered", (L_fp)select_, m, &a[a_offset], lda, 
		&sdim, &dwork[2], &dwork[*m + 2], &u[u_offset], ldu, &dwork[
		jwork], &availw, bwork, &ierr, (ftnlen)7, (ftnlen)11);
	if (ierr > 0) {
	    *info = ierr;
	    return 0;
	}
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__2);
    } else {
	jwork = 1;
	ia = 2;
	availw = *ldwork;
    }

    if (! schura) {

/*        Transform the right-hand side:  C <-- U'*C. */
/*        Workspace:  need   a+M, */
/*                    prefer a+M*N, */
/*                    where  a = 1+2*M, if FACTA =  'N', */
/*                           a = 0,     if FACTA <> 'N'. */

	chunka = availw / *m;
	blocka = min(chunka,*n) > 1;
	blas3a = chunka >= *n && blocka;

	if (blas3a) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[jwork], m, (
		    ftnlen)4);
	    dgemm_("Transpose", "NoTranspose", m, n, m, &c_b22, &u[u_offset], 
		    ldu, &dwork[jwork], m, &c_b23, &c__[c_offset], ldc, (
		    ftnlen)9, (ftnlen)11);
	} else if (blocka) {

/*           Use as many columns of C as possible. */

	    i__1 = *n;
	    i__2 = chunka;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
		i__3 = *n - j + 1;
		bl = min(i__3,chunka);
		dlacpy_("Full", m, &bl, &c__[j * c_dim1 + 1], ldc, &dwork[
			jwork], m, (ftnlen)4);
		dgemm_("Transpose", "NoTranspose", m, &bl, m, &c_b22, &u[
			u_offset], ldu, &dwork[jwork], m, &c_b23, &c__[j * 
			c_dim1 + 1], ldc, (ftnlen)9, (ftnlen)11);
/* L10: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &dwork[jwork], &c__1);
		dgemv_("Transpose", m, m, &c_b22, &u[u_offset], ldu, &dwork[
			jwork], &c__1, &c_b23, &c__[j * c_dim1 + 1], &c__1, (
			ftnlen)9);
/* L20: */
	    }

	}
/* Computing MAX */
	i__2 = maxwrk, i__1 = jwork + *m * *n - 1;
	maxwrk = max(i__2,i__1);
    }

    if (nofacb) {

/*        Compute the Schur factorization of B. */
/*        Workspace:  need   1+MAX(a-1,0)+5*N, */
/*                    prefer larger. */

	jwork = ia + (*n << 1);
	availw = *ldwork - jwork + 1;
	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &b[b_offset], ldb, 
		&sdim, &dwork[ia], &dwork[*n + ia], &v[v_offset], ldv, &dwork[
		jwork], &availw, bwork, &ierr, (ftnlen)7, (ftnlen)11);
	if (ierr > 0) {
	    *info = ierr + *m;
	    return 0;
	}
/* Computing MAX */
	i__2 = maxwrk, i__1 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__2,i__1);

	if (! schura) {

/*           Recompute the blocking parameters. */

	    chunka = availw / *m;
	    blocka = min(chunka,*n) > 1;
	    blas3a = chunka >= *n && blocka;
	}
    }

    if (! schurb) {

/*        Transform the right-hand side:  C <-- C*V. */
/*        Workspace:  need   a+b+N, */
/*                    prefer a+b+M*N, */
/*                    where  b = 2*N,   if FACTB =  'N', FACTA =  'N', */
/*                           b = 1+2*N, if FACTB =  'N', FACTA <> 'N', */
/*                           b = 0,     if FACTB <> 'N'. */

	chunkb = availw / *n;
	blockb = min(chunkb,*m) > 1;
	blas3b = chunkb >= *m && blockb;

	if (blas3b) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[jwork], m, (
		    ftnlen)4);
	    dgemm_("NoTranspose", "NoTranspose", m, n, n, &c_b22, &dwork[
		    jwork], m, &v[v_offset], ldv, &c_b23, &c__[c_offset], ldc,
		     (ftnlen)11, (ftnlen)11);
	} else if (blockb) {

/*           Use as many rows of C as possible. */

	    i__2 = *m;
	    i__1 = chunkb;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *m - i__ + 1;
		bl = min(i__3,chunkb);
		dlacpy_("Full", &bl, n, &c__[i__ + c_dim1], ldc, &dwork[jwork]
			, &bl, (ftnlen)4);
		dgemm_("NoTranspose", "NoTranspose", &bl, n, n, &c_b22, &
			dwork[jwork], &bl, &v[v_offset], ldv, &c_b23, &c__[
			i__ + c_dim1], ldc, (ftnlen)11, (ftnlen)11);
/* L30: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &c__[i__ + c_dim1], ldc, &dwork[jwork], &c__1);
		dgemv_("Transpose", n, n, &c_b22, &v[v_offset], ldv, &dwork[
			jwork], &c__1, &c_b23, &c__[i__ + c_dim1], ldc, (
			ftnlen)9);
/* L40: */
	    }

	}
/* Computing MAX */
	i__1 = maxwrk, i__2 = jwork + *m * *n - 1;
	maxwrk = max(i__1,i__2);
    }

/*     Solve the (transformed) equation. */
/*     Workspace for DICO = 'D':  a+b+2*M. */

    if (cont) {
	dtrsyl_(trana, tranb, isgn, m, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, scale, &ierr, (ftnlen)1, (ftnlen)1);
    } else {
	sb04py_(trana, tranb, isgn, m, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, scale, &dwork[jwork], &ierr, (
		ftnlen)1, (ftnlen)1);
/* Computing MAX */
	i__1 = maxwrk, i__2 = jwork + (*m << 1) - 1;
	maxwrk = max(i__1,i__2);
    }
    if (ierr > 0) {
	*info = *m + *n + 1;
    }

/*     Transform back the solution, if needed. */

    if (! schura) {

/*        Transform the right-hand side:  C <-- U*C. */
/*        Workspace:  need   a+b+M; */
/*                    prefer a+b+M*N. */

	if (blas3a) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[jwork], m, (
		    ftnlen)4);
	    dgemm_("NoTranspose", "NoTranspose", m, n, m, &c_b22, &u[u_offset]
		    , ldu, &dwork[jwork], m, &c_b23, &c__[c_offset], ldc, (
		    ftnlen)11, (ftnlen)11);
	} else if (blocka) {

/*           Use as many columns of C as possible. */

	    i__1 = *n;
	    i__2 = chunka;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
		i__3 = *n - j + 1;
		bl = min(i__3,chunka);
		dlacpy_("Full", m, &bl, &c__[j * c_dim1 + 1], ldc, &dwork[
			jwork], m, (ftnlen)4);
		dgemm_("NoTranspose", "NoTranspose", m, &bl, m, &c_b22, &u[
			u_offset], ldu, &dwork[jwork], m, &c_b23, &c__[j * 
			c_dim1 + 1], ldc, (ftnlen)11, (ftnlen)11);
/* L50: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &dwork[jwork], &c__1);
		dgemv_("NoTranspose", m, m, &c_b22, &u[u_offset], ldu, &dwork[
			jwork], &c__1, &c_b23, &c__[j * c_dim1 + 1], &c__1, (
			ftnlen)11);
/* L60: */
	    }

	}
    }

    if (! schurb) {

/*        Transform the right-hand side:  C <-- C*V'. */
/*        Workspace:  need   a+b+N; */
/*                    prefer a+b+M*N. */

	if (blas3b) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[jwork], m, (
		    ftnlen)4);
	    dgemm_("NoTranspose", "Transpose", m, n, n, &c_b22, &dwork[jwork],
		     m, &v[v_offset], ldv, &c_b23, &c__[c_offset], ldc, (
		    ftnlen)11, (ftnlen)9);
	} else if (blockb) {

/*           Use as many rows of C as possible. */

	    i__2 = *m;
	    i__1 = chunkb;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
		i__3 = *m - i__ + 1;
		bl = min(i__3,chunkb);
		dlacpy_("Full", &bl, n, &c__[i__ + c_dim1], ldc, &dwork[jwork]
			, &bl, (ftnlen)4);
		dgemm_("NoTranspose", "Transpose", &bl, n, n, &c_b22, &dwork[
			jwork], &bl, &v[v_offset], ldv, &c_b23, &c__[i__ + 
			c_dim1], ldc, (ftnlen)11, (ftnlen)9);
/* L70: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &c__[i__ + c_dim1], ldc, &dwork[jwork], &c__1);
		dgemv_("NoTranspose", n, n, &c_b22, &v[v_offset], ldv, &dwork[
			jwork], &c__1, &c_b23, &c__[i__ + c_dim1], ldc, (
			ftnlen)11);
/* L80: */
	    }

	}
    }

    dwork[1] = (doublereal) maxwrk;

    return 0;
/* *** Last line of SB04PD *** */
} /* sb04pd_ */

