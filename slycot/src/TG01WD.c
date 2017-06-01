/* TG01WD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 1.;
static doublereal c_b10 = 0.;
static integer c__1 = 1;

/* Subroutine */ int tg01wd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *q, integer *ldq, 
	doublereal *z__, integer *ldz, doublereal *alphar, doublereal *alphai,
	 doublereal *beta, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, bl, sdim;
    static logical blas3;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical block;
    extern /* Subroutine */ int dgges_(char *, char *, char *, L_fp, integer *
	    , doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen, ftnlen), dgemv_(char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical bwork[1];
    extern logical delctg_();
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer maxwrk;


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

/*     To reduce the pair (A,E) to a real generalized Schur form */
/*     by using an orthogonal equivalence transformation */
/*     (A,E) <-- (Q'*A*Z,Q'*E*Z) and to apply the transformation */
/*     to the matrices B and C: B <-- Q'*B and C <-- C*Z. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e., the order of the matrices A and E.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Q' * A * Z in an upper quasi-triangular form. */
/*             The elements below the first subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original descriptor matrix E. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Q' * E * Z in an upper triangular form. */
/*             The elements below the diagonal are set to zero. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix Q' * B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output matrix C * Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             The leading N-by-N part of this array contains the left */
/*             orthogonal transformation matrix used to reduce (A,E) to */
/*             the real generalized Schur form. */
/*             The columns of Q are the left generalized Schur vectors */
/*             of the pair (A,E). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= max(1,N). */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array contains the right */
/*             orthogonal transformation matrix used to reduce (A,E) to */
/*             the real generalized Schur form. */
/*             The columns of Z are the right generalized Schur vectors */
/*             of the pair (A,E). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= max(1,N). */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, if INFO = 0, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), */
/*             j=1,...,N, will be the generalized eigenvalues. */
/*             ALPHAR(j) + ALPHAI(j)*i, and BETA(j), j=1,...,N, are the */
/*             diagonals of the complex Schur form that would result if */
/*             the 2-by-2 diagonal blocks of the real Schur form of */
/*             (A,E) were further reduced to triangular form using */
/*             2-by-2 complex unitary transformations. */
/*             If ALPHAI(j) is zero, then the j-th eigenvalue is real; */
/*             if positive, then the j-th and (j+1)-st eigenvalues are a */
/*             complex conjugate pair, with ALPHAI(j+1) negative. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK.  LDWORK >= 8*N+16. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, the QZ algorithm failed to compute */
/*                   the generalized real Schur form; elements i+1:N of */
/*                   ALPHAR, ALPHAI, and BETA should be correct. */

/*     METHOD */

/*     The pair (A,E) is reduced to a real generalized Schur form using */
/*     an orthogonal equivalence transformation (A,E) <-- (Q'*A*Z,Q'*E*Z) */
/*     and the transformation is applied to the matrices B and C: */
/*     B <-- Q'*B and C <-- C*Z. */

/*     NUMERICAL ASPECTS */
/*                                     3 */
/*     The algorithm requires about 25N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     KEYWORDS */

/*     Orthogonal transformation, generalized real Schur form, similarity */
/*     transformation. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --alphar;
    --alphai;
    --beta;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Check the scalar input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*lde < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*p)) {
	*info = -11;
    } else if (*ldq < max(1,*n)) {
	*info = -13;
    } else if (*ldz < max(1,*n)) {
	*info = -15;
    } else if (*ldwork < (*n << 3) + 16) {
	*info = -20;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TG01WD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Reduce (A,E) to real generalized Schur form using an orthogonal */
/*     equivalence transformation (A,E) <-- (Q'*A*Z,Q'*E*Z), accumulate */
/*     the transformations in Q and Z, and compute the generalized */
/*     eigenvalues of the pair (A,E) in (ALPHAR, ALPHAI, BETA). */

/*     Workspace:  need   8*N+16; */
/*                 prefer larger. */

    dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_, n, &a[a_offset]
	    , lda, &e[e_offset], lde, &sdim, &alphar[1], &alphai[1], &beta[1],
	     &q[q_offset], ldq, &z__[z_offset], ldz, &dwork[1], ldwork, bwork,
	     info, (ftnlen)7, (ftnlen)7, (ftnlen)11);
    if (*info != 0) {
	return 0;
    }
    maxwrk = (integer) dwork[1];

/*     Apply the transformation: B <-- Q'*B. Use BLAS 3, if enough space. */

    chunk = *ldwork / *n;
    block = *m > 1;
    blas3 = chunk >= *m && block;

    if (blas3) {

/*        Enough workspace for a fast BLAS 3 algorithm. */

	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
	dgemm_("Transpose", "No transpose", n, m, n, &c_b9, &q[q_offset], ldq,
		 &dwork[1], n, &c_b10, &b[b_offset], ldb, (ftnlen)9, (ftnlen)
		12);

    } else if (block) {

/*        Use as many columns of B as possible. */

	i__1 = *m;
	i__2 = chunk;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = *m - j + 1;
	    bl = min(i__3,chunk);
	    dlacpy_("Full", n, &bl, &b[j * b_dim1 + 1], ldb, &dwork[1], n, (
		    ftnlen)4);
	    dgemm_("Transpose", "NoTranspose", n, &bl, n, &c_b9, &q[q_offset],
		     ldq, &dwork[1], n, &c_b10, &b[j * b_dim1 + 1], ldb, (
		    ftnlen)9, (ftnlen)11);
/* L10: */
	}

    } else {

/*        Use a BLAS 2 algorithm. Here, M <= 1. */

	if (*m > 0) {
	    dcopy_(n, &b[b_offset], &c__1, &dwork[1], &c__1);
	    dgemv_("Transpose", n, n, &c_b9, &q[q_offset], ldq, &dwork[1], &
		    c__1, &c_b10, &b[b_offset], &c__1, (ftnlen)9);
	}
    }
/* Computing MAX */
    i__2 = maxwrk, i__1 = *n * *m;
    maxwrk = max(i__2,i__1);

/*     Apply the transformation: C <-- C*Z.  Use BLAS 3, if enough space. */

    block = *p > 1;
    blas3 = chunk >= *p && block;

    if (blas3) {
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], p, (ftnlen)4);
	dgemm_("No transpose", "No transpose", p, n, n, &c_b9, &dwork[1], p, &
		z__[z_offset], ldz, &c_b10, &c__[c_offset], ldc, (ftnlen)12, (
		ftnlen)12);

    } else if (block) {

/*        Use as many rows of C as possible. */

	i__2 = *p;
	i__1 = chunk;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
	    i__3 = *p - i__ + 1;
	    bl = min(i__3,chunk);
	    dlacpy_("Full", &bl, n, &c__[i__ + c_dim1], ldc, &dwork[1], &bl, (
		    ftnlen)4);
	    dgemm_("NoTranspose", "NoTranspose", &bl, n, n, &c_b9, &dwork[1], 
		    &bl, &z__[z_offset], ldz, &c_b10, &c__[i__ + c_dim1], ldc,
		     (ftnlen)11, (ftnlen)11);
/* L20: */
	}

    } else {

/*        Use a BLAS 2 algorithm. Here, P <= 1. */

	if (*p > 0) {
	    dcopy_(n, &c__[c_offset], ldc, &dwork[1], &c__1);
	    dgemv_("Transpose", n, n, &c_b9, &z__[z_offset], ldz, &dwork[1], &
		    c__1, &c_b10, &c__[c_offset], ldc, (ftnlen)9);
	}

    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = *p * *n;
    maxwrk = max(i__1,i__2);

    dwork[1] = (doublereal) maxwrk;

    return 0;
/* *** Last line of TG01WD *** */
} /* tg01wd_ */

