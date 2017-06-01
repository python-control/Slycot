/* TB01KD.f -- translated by f2c (version 20100827).
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

static integer c_n1 = -1;
static doublereal c_b14 = 1.;
static doublereal c_b22 = 0.;

/* Subroutine */ int tb01kd_(char *dico, char *stdom, char *joba, integer *n, 
	integer *m, integer *p, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	integer *ndim, doublereal *u, integer *ldu, doublereal *wr, 
	doublereal *wi, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen stdom_len, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u_dim1, 
	    u_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer nr, ndim1;
    static doublereal scale;
    extern /* Subroutine */ int tb01ld_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static logical ljobg;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrsyl_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen);


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

/*     To compute an additive spectral decomposition of the transfer- */
/*     function matrix of the system (A,B,C) by reducing the system */
/*     state-matrix A to a block-diagonal form. */
/*     The system matrices are transformed as */
/*     A <-- inv(U)*A*U, B <--inv(U)*B and C <-- C*U. */
/*     The leading diagonal block of the resulting A has eigenvalues */
/*     in a suitably defined domain of interest. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     STDOM   CHARACTER*1 */
/*             Specifies whether the domain of interest is of stability */
/*             type (left part of complex plane or inside of a circle) */
/*             or of instability type (right part of complex plane or */
/*             outside of a circle) as follows: */
/*             = 'S':  stability type domain; */
/*             = 'U':  instability type domain. */

/*     JOBA    CHARACTER*1 */
/*             Specifies the shape of the state dynamics matrix on entry */
/*             as follows: */
/*             = 'S':  A is in an upper real Schur form; */
/*             = 'G':  A is a general square dense matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION. */
/*             Specifies the boundary of the domain of interest for the */
/*             eigenvalues of A. For a continuous-time system */
/*             (DICO = 'C'), ALPHA is the boundary value for the real */
/*             parts of eigenvalues, while for a discrete-time system */
/*             (DICO = 'D'), ALPHA >= 0 represents the boundary value for */
/*             the moduli of eigenvalues. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the unreduced state dynamics matrix A. */
/*             If JOBA = 'S' then A must be a matrix in real Schur form. */
/*             On exit, the leading N-by-N part of this array contains a */
/*             block diagonal matrix inv(U) * A * U with two diagonal */
/*             blocks in real Schur form with the elements below the */
/*             first subdiagonal set to zero. */
/*             The leading NDIM-by-NDIM block of A has eigenvalues in the */
/*             domain of interest and the trailing (N-NDIM)-by-(N-NDIM) */
/*             block has eigenvalues outside the domain of interest. */
/*             The domain of interest for lambda(A), the eigenvalues */
/*             of A, is defined by the parameters ALPHA, DICO and STDOM */
/*             as follows: */
/*             For a continuous-time system (DICO = 'C'): */
/*               Real(lambda(A)) < ALPHA if STDOM = 'S'; */
/*               Real(lambda(A)) > ALPHA if STDOM = 'U'; */
/*             For a discrete-time system (DICO = 'D'): */
/*               Abs(lambda(A)) < ALPHA if STDOM = 'S'; */
/*               Abs(lambda(A)) > ALPHA if STDOM = 'U'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix inv(U) * B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output matrix C * U. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NDIM    (output) INTEGER */
/*             The number of eigenvalues of A lying inside the domain of */
/*             interest for eigenvalues. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array contains the */
/*             transformation matrix used to reduce A to the block- */
/*             diagonal form. The first NDIM columns of U span the */
/*             invariant subspace of A corresponding to the eigenvalues */
/*             of its leading diagonal block. The last N-NDIM columns */
/*             of U span the reducing subspace of A corresponding to */
/*             the eigenvalues of the trailing diagonal block of A. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= max(1,N). */

/*     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N) */
/*             WR and WI contain the real and imaginary parts, */
/*             respectively, of the computed eigenvalues of A. The */
/*             eigenvalues will be in the same order that they appear on */
/*             the diagonal of the output real Schur form of A. Complex */
/*             conjugate pairs of eigenvalues will appear consecutively */
/*             with the eigenvalue having the positive imaginary part */
/*             first. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX(1,N)   if JOBA = 'S'; */
/*             LDWORK >= MAX(1,3*N) if JOBA = 'G'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0: successful exit; */
/*             < 0: if INFO = -i, the i-th argument had an illegal */
/*                  value; */
/*             = 1: the QR algorithm failed to compute all the */
/*                  eigenvalues of A; */
/*             = 2: a failure occured during the ordering of the real */
/*                  Schur form of A; */
/*             = 3: the separation of the two diagonal blocks failed */
/*                  because of very close eigenvalues. */

/*     METHOD */

/*     A similarity transformation U is determined that reduces the */
/*     system state-matrix A to a block-diagonal form (with two diagonal */
/*     blocks), so that the leading diagonal block of the resulting A has */
/*     eigenvalues in a specified domain of the complex plane. The */
/*     determined transformation is applied to the system (A,B,C) as */
/*       A <-- inv(U)*A*U, B <-- inv(U)*B and C <-- C*U. */

/*     REFERENCES */

/*     [1] Safonov, M.G., Jonckheere, E.A., Verma, M., Limebeer, D.J.N. */
/*         Synthesis of positive real multivariable feedback systems. */
/*         Int. J. Control, pp. 817-842, 1987. */

/*     NUMERICAL ASPECTS */
/*                                     3 */
/*     The algorithm requires about 14N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SADSDC. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Invariant subspace, real Schur form, similarity transformation, */
/*     spectral factorization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
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
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    ljobg = lsame_(joba, "G", (ftnlen)1, (ftnlen)1);

/*     Check input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (lsame_(stdom, "S", (ftnlen)1, (ftnlen)1) || lsame_(stdom, 
	    "U", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lsame_(joba, "S", (ftnlen)1, (ftnlen)1) || ljobg)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (discr && *alpha < 0.) {
	*info = -7;
    } else if (*lda < max(1,*n)) {
	*info = -9;
    } else if (*ldb < max(1,*n)) {
	*info = -11;
    } else if (*ldc < max(1,*p)) {
	*info = -13;
    } else if (*ldu < max(1,*n)) {
	*info = -16;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3;
	if (*ldwork < max(1,*n) || *ldwork < max(i__1,i__2) && ljobg) {
	    *info = -20;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB01KD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *ndim = 0;
    if (*n == 0) {
	return 0;
    }

/*     Reduce A to an ordered real Schur form using an orthogonal */
/*     similarity transformation A <- U'*A*U and accumulate the */
/*     transformations in U. The reordering of the real Schur form of A */
/*     is performed in accordance with the values of the parameters DICO, */
/*     STDOM and ALPHA. Apply the transformation to B and C: B <- U'*B */
/*     and C <- C*U. The eigenvalues of A are computed in (WR,WI). */

/*     Workspace:  need   3*N (if JOBA = 'G'), or N (if JOBA = 'S'); */
/*                 prefer larger. */

    tb01ld_(dico, stdom, joba, n, m, p, alpha, &a[a_offset], lda, &b[b_offset]
	    , ldb, &c__[c_offset], ldc, ndim, &u[u_offset], ldu, &wr[1], &wi[
	    1], &dwork[1], ldwork, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

    if (*info != 0) {
	return 0;
    }

    if (*ndim > 0 && *ndim < *n) {

/*        Reduce A to a block-diagonal form by a similarity */
/*        transformation of the form */
/*               -1                  ( I -X ) */
/*         A <- T  AT,  where    T = (      )  and X satisfies the */
/*                                   ( 0  I ) */
/*        Sylvester equation */

/*          A11*X - X*A22 = A12. */

	nr = *n - *ndim;
	ndim1 = *ndim + 1;
	dtrsyl_("N", "N", &c_n1, ndim, &nr, &a[a_offset], lda, &a[ndim1 + 
		ndim1 * a_dim1], lda, &a[ndim1 * a_dim1 + 1], lda, &scale, 
		info, (ftnlen)1, (ftnlen)1);
	if (*info != 0) {
	    *info = 3;
	    return 0;
	}
/*                      -1 */
/*        Compute B <- T  B,  C <- CT,  U <- UT. */

	scale = 1. / scale;
	dgemm_("N", "N", ndim, m, &nr, &scale, &a[ndim1 * a_dim1 + 1], lda, &
		b[ndim1 + b_dim1], ldb, &c_b14, &b[b_offset], ldb, (ftnlen)1, 
		(ftnlen)1);
	d__1 = -scale;
	dgemm_("N", "N", p, &nr, ndim, &d__1, &c__[c_offset], ldc, &a[ndim1 * 
		a_dim1 + 1], lda, &c_b14, &c__[ndim1 * c_dim1 + 1], ldc, (
		ftnlen)1, (ftnlen)1);
	d__1 = -scale;
	dgemm_("N", "N", n, &nr, ndim, &d__1, &u[u_offset], ldu, &a[ndim1 * 
		a_dim1 + 1], lda, &c_b14, &u[ndim1 * u_dim1 + 1], ldu, (
		ftnlen)1, (ftnlen)1);

/*        Set A12 to zero. */

	dlaset_("Full", ndim, &nr, &c_b22, &c_b22, &a[ndim1 * a_dim1 + 1], 
		lda, (ftnlen)4);
    }

/*     Set to zero the lower triangular part under the first subdiagonal */
/*     of A. */

    if (*n > 2) {
	i__1 = *n - 2;
	i__2 = *n - 2;
	dlaset_("L", &i__1, &i__2, &c_b22, &c_b22, &a[a_dim1 + 3], lda, (
		ftnlen)1);
    }
    return 0;
/* *** Last line of TB01KD *** */
} /* tb01kd_ */

