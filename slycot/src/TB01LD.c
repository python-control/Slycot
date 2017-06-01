/* TB01LD.f -- translated by f2c (version 20100827).
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

static doublereal c_b13 = 0.;
static doublereal c_b14 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tb01ld_(char *dico, char *stdom, char *joba, integer *n, 
	integer *m, integer *p, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	integer *ndim, doublereal *u, integer *ldu, doublereal *wr, 
	doublereal *wi, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen stdom_len, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u_dim1, 
	    u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, sdim, ierr, ldwp;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), mb03qd_(char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical ljobg;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    extern /* Subroutine */ int mb03qx_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical bwork[1];
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static doublereal wrkopt;


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

/*     To reduce the system state matrix A to an ordered upper real */
/*     Schur form by using an orthogonal similarity transformation */
/*     A <-- U'*A*U and to apply the transformation to the matrices */
/*     B and C: B <-- U'*B and C <-- C*U. */
/*     The leading block of the resulting A has eigenvalues in a */
/*     suitably defined domain of interest. */

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
/*             (DICO = 'D'), ALPHA >= 0 represents the boundary value */
/*             for the moduli of eigenvalues. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the unreduced state dynamics matrix A. */
/*             If JOBA = 'S' then A must be a matrix in real Schur form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the ordered real Schur matrix U' * A * U with the elements */
/*             below the first subdiagonal set to zero. */
/*             The leading NDIM-by-NDIM part of A has eigenvalues in the */
/*             domain of interest and the trailing (N-NDIM)-by-(N-NDIM) */
/*             part has eigenvalues outside the domain of interest. */
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
/*             the transformed input matrix U' * B. */

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
/*             orthogonal transformation matrix used to reduce A to the */
/*             real Schur form and/or to reorder the diagonal blocks of */
/*             real Schur form of A. The first NDIM columns of U form */
/*             an orthogonal basis for the invariant subspace of A */
/*             corresponding to the first NDIM eigenvalues. */

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
/*                  Schur form of A. */

/*     METHOD */

/*     Matrix A is reduced to an ordered upper real Schur form using an */
/*     orthogonal similarity transformation A <-- U'*A*U. This */
/*     transformation is determined so that the leading block of the */
/*     resulting A has eigenvalues in a suitably defined domain of */
/*     interest. Then, the transformation is applied to the matrices B */
/*     and C: B <-- U'*B and C <-- C*U. */

/*     NUMERICAL ASPECTS */
/*                                     3 */
/*     The algorithm requires about 14N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRSFOD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001. */

/*     KEYWORDS */

/*     Invariant subspace, orthogonal transformation, real Schur form, */
/*     similarity transformation. */

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
	xerbla_("TB01LD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *ndim = 0;
    if (*n == 0) {
	return 0;
    }

    if (lsame_(joba, "G", (ftnlen)1, (ftnlen)1)) {

/*        Reduce A to real Schur form using an orthogonal similarity */
/*        transformation A <- U'*A*U, accumulate the transformation in U */
/*        and compute the eigenvalues of A in (WR,WI). */

/*        Workspace:  need   3*N; */
/*                    prefer larger. */

	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, 
		&sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)7, (ftnlen)11);
	wrkopt = dwork[1];
	if (*info != 0) {
	    *info = 1;
	    return 0;
	}
    } else {

/*        Initialize U with an identity matrix. */

	dlaset_("Full", n, n, &c_b13, &c_b14, &u[u_offset], ldu, (ftnlen)4);
	wrkopt = 0.;
    }

/*     Separate the spectrum of A. The leading NDIM-by-NDIM submatrix of */
/*     A corresponds to the eigenvalues of interest. */
/*     Workspace:  need   N. */

    mb03qd_(dico, stdom, "Update", n, &c__1, n, alpha, &a[a_offset], lda, &u[
	    u_offset], ldu, ndim, &dwork[1], info, (ftnlen)1, (ftnlen)1, (
	    ftnlen)6);
    if (*info != 0) {
	return 0;
    }

/*     Compute the eigenvalues. */

    mb03qx_(n, &a[a_offset], lda, &wr[1], &wi[1], &ierr);

/*     Apply the transformation: B <-- U'*B. */

    if (*ldwork < *n * *m) {

/*        Not enough working space for using DGEMM. */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[1], &c__1);
	    dgemv_("Transpose", n, n, &c_b14, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b13, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)9);
/* L10: */
	}

    } else {
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
	dgemm_("Transpose", "No transpose", n, m, n, &c_b14, &u[u_offset], 
		ldu, &dwork[1], n, &c_b13, &b[b_offset], ldb, (ftnlen)9, (
		ftnlen)12);
/* Computing MAX */
	d__1 = wrkopt, d__2 = (doublereal) (*n * *m);
	wrkopt = max(d__1,d__2);
    }

/*     Apply the transformation: C <-- C*U. */

    if (*ldwork < *n * *p) {

/*        Not enough working space for using DGEMM. */

	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dcopy_(n, &c__[i__ + c_dim1], ldc, &dwork[1], &c__1);
	    dgemv_("Transpose", n, n, &c_b14, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b13, &c__[i__ + c_dim1], ldc, (ftnlen)9);
/* L20: */
	}

    } else {
	ldwp = max(1,*p);
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], &ldwp, (ftnlen)
		4);
	dgemm_("No transpose", "No transpose", p, n, n, &c_b14, &dwork[1], &
		ldwp, &u[u_offset], ldu, &c_b13, &c__[c_offset], ldc, (ftnlen)
		12, (ftnlen)12);
/* Computing MAX */
	d__1 = wrkopt, d__2 = (doublereal) (*n * *p);
	wrkopt = max(d__1,d__2);
    }

    dwork[1] = wrkopt;

    return 0;
/* *** Last line of TB01LD *** */
} /* tb01ld_ */

