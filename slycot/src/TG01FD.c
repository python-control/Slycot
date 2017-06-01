/* TG01FD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static doublereal c_b42 = 0.;
static doublereal c_b43 = 1.;

/* Subroutine */ int tg01fd_(char *compq, char *compz, char *joba, integer *l,
	 integer *n, integer *m, integer *p, doublereal *a, integer *lda, 
	doublereal *e, integer *lde, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *q, integer *ldq, doublereal *z__, 
	integer *ldz, integer *ranke, integer *rnka22, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen compq_len, ftnlen compz_len, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, nb, lh, ln, kw, ir1, la22, na22;
    static logical ilq, ilz;
    static integer lwr, ire1;
    static logical reda;
    static doublereal sval[3];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static logical withb, withc, redtr;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icompq, icompz;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal svlmax;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer wrkopt;


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

/*     To compute for the descriptor system (A-lambda E,B,C) */
/*     the orthogonal transformation matrices Q and Z such that the */
/*     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is */
/*     in a SVD-like coordinate form with */

/*                  ( A11  A12 )             ( Er  0 ) */
/*         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) , */
/*                  ( A21  A22 )             (  0  0 ) */

/*     where Er is an upper triangular invertible matrix. */
/*     Optionally, the A22 matrix can be further reduced to the form */

/*                  ( Ar  X ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix, and X either a full */
/*     or a zero matrix. */
/*     The left and/or right orthogonal transformations performed */
/*     to reduce E and A22 can be optionally accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPQ   CHARACTER*1 */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     orthogonal matrix Q is returned; */
/*             = 'U':  Q must contain an orthogonal matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     orthogonal matrix Z is returned; */
/*             = 'U':  Z must contain an orthogonal matrix Z1 on entry, */
/*                     and the product Z1*Z is returned. */

/*     JOBA    CHARACTER*1 */
/*             = 'N':  do not reduce A22; */
/*             = 'R':  reduce A22 to a SVD-like upper triangular form. */
/*             = 'T':  reduce A22 to an upper trapezoidal form. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*A*Z. If JOBA = 'T', this matrix */
/*             is in the form */

/*                           ( A11  *   *  ) */
/*                  Q'*A*Z = (  *   Ar  X  ) , */
/*                           (  *   0   0  ) */

/*             where A11 is a RANKE-by-RANKE matrix and Ar is a */
/*             RNKA22-by-RNKA22 invertible upper triangular matrix. */
/*             If JOBA = 'R' then A has the above form with X = 0. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*E*Z. */

/*                      ( Er  0 ) */
/*             Q'*E*Z = (       ) , */
/*                      (  0  0 ) */

/*             where Er is a RANKE-by-RANKE upper triangular invertible */
/*             matrix. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the transformed matrix Q'*B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,L) */
/*             If COMPQ = 'N':  Q is not referenced. */
/*             If COMPQ = 'I':  on entry, Q need not be set; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the orthogonal matrix Q, */
/*                              where Q' is the product of Householder */
/*                              transformations which are applied to A, */
/*                              E, and B on the left. */
/*             If COMPQ = 'U':  on entry, the leading L-by-L part of this */
/*                              array must contain an orthogonal matrix */
/*                              Q1; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the orthogonal matrix */
/*                              Q1*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If COMPZ = 'N':  Z is not referenced. */
/*             If COMPZ = 'I':  on entry, Z need not be set; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the orthogonal matrix Z, */
/*                              which is the product of Householder */
/*                              transformations applied to A, E, and C */
/*                              on the right. */
/*             If COMPZ = 'U':  on entry, the leading N-by-N part of this */
/*                              array must contain an orthogonal matrix */
/*                              Z1; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the orthogonal matrix */
/*                              Z1*Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. */
/*             LDZ >= 1,        if COMPZ = 'N'; */
/*             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'. */

/*     RANKE   (output) INTEGER */
/*             The estimated rank of matrix E, and thus also the order */
/*             of the invertible upper triangular submatrix Er. */

/*     RNKA22  (output) INTEGER */
/*             If JOBA = 'R' or 'T', then RNKA22 is the estimated rank of */
/*             matrix A22, and thus also the order of the invertible */
/*             upper triangular submatrix Ar. */
/*             If JOBA = 'N', then RNKA22 is not referenced. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the rank of E */
/*             and of A22. If the user sets TOL > 0, then the given */
/*             value of TOL is used as a lower bound for the */
/*             reciprocal condition numbers of leading submatrices */
/*             of R or R22 in the QR decompositions E * P = Q * R of E */
/*             or A22 * P22 = Q22 * R22 of A22. */
/*             A submatrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = L*N*EPS,  is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH). TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, N+P, MIN(L,N)+MAX(3*N-1,M,L) ). */
/*             For optimal performance, LDWORK should be larger. */

/*             If LDWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             DWORK array, returns this value as the first entry of */
/*             the DWORK array, and no error message related to LDWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated QR factorization with column */
/*     pivoting of E, in the form */

/*                       ( E11 E12 ) */
/*           E * P = Q * (         ) */
/*                       (  0  E22 ) */

/*     and finds the largest RANKE-by-RANKE leading submatrix E11 whose */
/*     estimated condition number is less than 1/TOL. RANKE defines thus */
/*     the rank of matrix E. Further E22, being negligible, is set to */
/*     zero, and an orthogonal matrix Y is determined such that */

/*           ( E11 E12 ) = ( Er  0 ) * Y . */

/*     The overal transformation matrix Z results as Z = P * Y' and the */
/*     resulting transformed matrices Q'*A*Z and Q'*E*Z have the form */

/*                          ( Er  0 )                      ( A11  A12 ) */
/*         E <- Q'* E * Z = (       ) ,  A <- Q' * A * Z = (          ) , */
/*                          (  0  0 )                      ( A21  A22 ) */

/*     where Er is an upper triangular invertible matrix. */
/*     If JOBA = 'R' the same reduction is performed on A22 to obtain it */
/*     in the form */

/*                  ( Ar  0 ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix. */
/*     If JOBA = 'T' then A22 is row compressed using the QR */
/*     factorization with column pivoting to the form */

/*                  ( Ar  X ) */
/*            A22 = (       ) */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix. */

/*     The transformations are also applied to the rest of system */
/*     matrices */

/*          B <- Q' * B, C <- C * Z. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( L*L*N )  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSSV. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     May 2003, Jan. 2009. */

/*     KEYWORDS */

/*     Descriptor system, matrix algebra, matrix operations, */
/*     orthogonal transformation. */

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

/*     Decode COMPQ. */

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
    --iwork;
    --dwork;

    /* Function Body */
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
	ilq = FALSE_;
	icompq = 1;
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
	ilq = TRUE_;
	icompq = 2;
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
	ilq = TRUE_;
	icompq = 3;
    } else {
	icompq = 0;
    }

/*     Decode COMPZ. */

    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
	ilz = FALSE_;
	icompz = 1;
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
	ilz = TRUE_;
	icompz = 2;
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
	ilz = TRUE_;
	icompz = 3;
    } else {
	icompz = 0;
    }
    reda = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
    redtr = lsame_(joba, "T", (ftnlen)1, (ftnlen)1);
    withb = *m > 0;
    withc = *p > 0;
    lquery = *ldwork == -1;

/*     Test the input parameters. */

    ln = min(*l,*n);
    *info = 0;
/* Computing MAX */
/* Computing MAX */
    i__3 = *n * 3 - 1, i__3 = max(i__3,*m);
    i__1 = 1, i__2 = *n + *p, i__1 = max(i__1,i__2), i__2 = ln + max(i__3,*l);
    wrkopt = max(i__1,i__2);
    if (icompq <= 0) {
	*info = -1;
    } else if (icompz <= 0) {
	*info = -2;
    } else if (! lsame_(joba, "N", (ftnlen)1, (ftnlen)1) && ! reda && ! redtr)
	     {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*p < 0) {
	*info = -7;
    } else if (*lda < max(1,*l)) {
	*info = -9;
    } else if (*lde < max(1,*l)) {
	*info = -11;
    } else if (*ldb < 1 || withb && *ldb < *l) {
	*info = -13;
    } else if (*ldc < max(1,*p)) {
	*info = -15;
    } else if (ilq && *ldq < *l || *ldq < 1) {
	*info = -17;
    } else if (ilz && *ldz < *n || *ldz < 1) {
	*info = -19;
    } else if (*tol >= 1.) {
	*info = -22;
    } else {
	if (lquery) {
/* Computing MIN */
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LC", l, n, &ln, &c_n1,
		     (ftnlen)6, (ftnlen)2);
	    nb = min(i__1,i__2);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = ln + *n * nb;
	    wrkopt = max(i__1,i__2);
	    if (withb) {
/* Computing MIN */
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LC", l, m, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = ln + *m * nb;
		wrkopt = max(i__1,i__2);
	    }
	    if (ilq) {
/* Computing MIN */
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "RN", l, l, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = ln + *l * nb;
		wrkopt = max(i__1,i__2);
	    }
	    nb = ilaenv_(&c__1, "DGERQF", " ", l, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = ln + *n * nb;
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", l, n, n, &c_n1, (
		    ftnlen)6, (ftnlen)2);
	    nb = min(i__1,i__2);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = *n + max(1,*l) * nb;
	    wrkopt = max(i__1,i__2);
	    if (withc) {
/* Computing MIN */
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", p, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = *n + max(1,*p) * nb;
		wrkopt = max(i__1,i__2);
	    }
	    if (ilz) {
/* Computing MIN */
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", n, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = *n + max(1,*n) * nb;
		wrkopt = max(i__1,i__2);
	    }
	} else if (*ldwork < wrkopt) {
	    *info = -25;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TG01FD", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

/*     Initialize Q and Z if necessary. */

    if (icompq == 3) {
	dlaset_("Full", l, l, &c_b42, &c_b43, &q[q_offset], ldq, (ftnlen)4);
    }
    if (icompz == 3) {
	dlaset_("Full", n, n, &c_b42, &c_b43, &z__[z_offset], ldz, (ftnlen)4);
    }

/*     Quick return if possible. */

    if (*l == 0 || *n == 0) {
	dwork[1] = 1.;
	*ranke = 0;
	if (reda || redtr) {
	    *rnka22 = 0;
	}
	return 0;
    }

    toldef = *tol;
    if (toldef <= 0.) {

/*        Use the default tolerance for rank determination. */

	toldef = (doublereal) (*l * *n) * dlamch_("EPSILON", (ftnlen)7);
    }

/*     Set the estimate of maximum singular value of E to */
/*     max(||E||,||A||) to detect negligible A or E matrices. */

/* Computing MAX */
    d__1 = dlange_("F", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)1), d__2 =
	     dlange_("F", l, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
    svlmax = max(d__1,d__2);

/*     Compute the rank-revealing QR decomposition of E, */

/*                        ( E11 E12 ) */
/*           E * P = Qr * (         ) , */
/*                        (  0  E22 ) */

/*     and determine the rank of E using incremental condition */
/*     estimation. */
/*     Workspace: MIN(L,N) + 3*N - 1. */

    lwr = *ldwork - ln;
    kw = ln + 1;

    mb03oy_(l, n, &e[e_offset], lde, &toldef, &svlmax, ranke, sval, &iwork[1],
	     &dwork[1], &dwork[kw], info);

/*     Apply transformation on the rest of matrices. */

    if (*ranke > 0) {

/*        A <-- Qr' * A. */
/*        Workspace: need   MIN(L,N) + N; */
/*                   prefer MIN(L,N) + N*NB. */

	dormqr_("Left", "Transpose", l, n, ranke, &e[e_offset], lde, &dwork[1]
		, &a[a_offset], lda, &dwork[kw], &lwr, info, (ftnlen)4, (
		ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
	wrkopt = max(i__1,i__2);

/*        B <-- Qr' * B. */
/*        Workspace: need   MIN(L,N) + M; */
/*                   prefer MIN(L,N) + M*NB. */

	if (withb) {
	    dormqr_("Left", "Transpose", l, m, ranke, &e[e_offset], lde, &
		    dwork[1], &b[b_offset], ldb, &dwork[kw], &lwr, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
	    wrkopt = max(i__1,i__2);
	}

/*        Q <-- Q * Qr. */
/*        Workspace: need   MIN(L,N) + L; */
/*                   prefer MIN(L,N) + L*NB. */

	if (ilq) {
	    dormqr_("Right", "No Transpose", l, l, ranke, &e[e_offset], lde, &
		    dwork[1], &q[q_offset], ldq, &dwork[kw], &lwr, info, (
		    ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
	    wrkopt = max(i__1,i__2);
	}

/*        Set lower triangle of E to zero. */

	if (*l >= 2) {
	    i__1 = *l - 1;
	    dlaset_("Lower", &i__1, ranke, &c_b42, &c_b42, &e[e_dim1 + 2], 
		    lde, (ftnlen)5);
	}

/*        Compute A*P, C*P and Z*P by forward permuting the columns of */
/*        A, C and Z based on information in IWORK. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    iwork[j] = -iwork[j];
/* L10: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iwork[i__] < 0) {
		j = i__;
		iwork[j] = -iwork[j];
L20:
		k = iwork[j];
		if (iwork[k] < 0) {
		    dswap_(l, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			    c__1);
		    if (withc) {
			dswap_(p, &c__[j * c_dim1 + 1], &c__1, &c__[k * 
				c_dim1 + 1], &c__1);
		    }
		    if (ilz) {
			dswap_(n, &z__[j * z_dim1 + 1], &c__1, &z__[k * 
				z_dim1 + 1], &c__1);
		    }
		    iwork[k] = -iwork[k];
		    j = k;
		    goto L20;
		}
	    }
/* L30: */
	}

/*        Determine an orthogonal matrix Y such that */

/*           ( E11 E12 ) = ( Er  0 ) * Y . */

/*        Compute E <-- E*Y', A <-- A*Y', C <-- C*Y', Z <-- Z*Y'. */

	if (*ranke < *n) {

/*           Workspace: need   2*N; */
/*                      prefer N + N*NB. */

	    kw = *ranke + 1;
	    i__1 = *ldwork - kw + 1;
	    dtzrzf_(ranke, n, &e[e_offset], lde, &dwork[1], &dwork[kw], &i__1,
		     info);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);

/*           Workspace: need   N + MAX(L,P,N); */
/*                      prefer N + MAX(L,P,N)*NB. */

	    lh = *n - *ranke;
	    i__1 = *ldwork - kw + 1;
	    dormrz_("Right", "Transpose", l, n, ranke, &lh, &e[e_offset], lde,
		     &dwork[1], &a[a_offset], lda, &dwork[kw], &i__1, info, (
		    ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	    if (withc) {
		i__1 = *ldwork - kw + 1;
		dormrz_("Right", "Transpose", p, n, ranke, &lh, &e[e_offset], 
			lde, &dwork[1], &c__[c_offset], ldc, &dwork[kw], &
			i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
		i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
		wrkopt = max(i__1,i__2);
	    }
	    if (ilz) {
		i__1 = *ldwork - kw + 1;
		dormrz_("Right", "Transpose", n, n, ranke, &lh, &e[e_offset], 
			lde, &dwork[1], &z__[z_offset], ldz, &dwork[kw], &
			i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
		i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
		wrkopt = max(i__1,i__2);
	    }

/*           Set E12 and E22 to zero. */

	    dlaset_("Full", l, &lh, &c_b42, &c_b42, &e[kw * e_dim1 + 1], lde, 
		    (ftnlen)4);
	}
    } else {
	dlaset_("Full", l, n, &c_b42, &c_b42, &e[e_offset], lde, (ftnlen)4);
    }

/*     Reduce A22 if necessary. */

    if (reda || redtr) {
	la22 = *l - *ranke;
	na22 = *n - *ranke;
	if (min(la22,na22) == 0) {
	    *rnka22 = 0;
	} else {

/*           Compute the rank-revealing QR decomposition of A22, */

/*                              ( R11 R12 ) */
/*              A22 * P2 = Q2 * (         ) , */
/*                              (  0  R22 ) */

/*           and determine the rank of A22 using incremental */
/*           condition estimation. */
/*           Workspace: MIN(L,N) + 3*N - 1. */

	    ir1 = *ranke + 1;
	    mb03oy_(&la22, &na22, &a[ir1 + ir1 * a_dim1], lda, &toldef, &
		    svlmax, rnka22, sval, &iwork[1], &dwork[1], &dwork[kw], 
		    info);

/*           Apply transformation on the rest of matrices. */

	    if (*rnka22 > 0) {

/*              A <-- diag(I, Q2') * A */
/*              Workspace: need   MIN(L,N) + N; */
/*                         prefer MIN(L,N) + N*NB. */

		dormqr_("Left", "Transpose", &la22, ranke, rnka22, &a[ir1 + 
			ir1 * a_dim1], lda, &dwork[1], &a[ir1 + a_dim1], lda, 
			&dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9);

/*              B <-- diag(I, Q2') * B */
/*              Workspace: need   MIN(L,N) + M; */
/*                         prefer MIN(L,N) + M*NB. */

		if (withb) {
		    dormqr_("Left", "Transpose", &la22, m, rnka22, &a[ir1 + 
			    ir1 * a_dim1], lda, &dwork[1], &b[ir1 + b_dim1], 
			    ldb, &dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9)
			    ;
		}

/*              Q <-- Q * diag(I, Q2) */
/*              Workspace: need   MIN(L,N) + L; */
/*                         prefer MIN(L,N) + L*NB. */

		if (ilq) {
		    dormqr_("Right", "No transpose", l, &la22, rnka22, &a[ir1 
			    + ir1 * a_dim1], lda, &dwork[1], &q[ir1 * q_dim1 
			    + 1], ldq, &dwork[kw], &lwr, info, (ftnlen)5, (
			    ftnlen)12);
		}

/*              Set lower triangle of A22 to zero. */

		if (la22 >= 2) {
		    i__1 = la22 - 1;
		    dlaset_("Lower", &i__1, rnka22, &c_b42, &c_b42, &a[ir1 + 
			    1 + ir1 * a_dim1], lda, (ftnlen)5);
		}

/*              Compute A*diag(I,P2), C*diag(I,P2) and Z*diag(I,P2) */
/*              by forward permuting the columns of A, C and Z based */
/*              on information in IWORK. */

		i__1 = na22;
		for (j = 1; j <= i__1; ++j) {
		    iwork[j] = -iwork[j];
/* L40: */
		}
		i__1 = na22;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (iwork[i__] < 0) {
			j = i__;
			iwork[j] = -iwork[j];
L50:
			k = iwork[j];
			if (iwork[k] < 0) {
			    dswap_(ranke, &a[(*ranke + j) * a_dim1 + 1], &
				    c__1, &a[(*ranke + k) * a_dim1 + 1], &
				    c__1);
			    if (withc) {
				dswap_(p, &c__[(*ranke + j) * c_dim1 + 1], &
					c__1, &c__[(*ranke + k) * c_dim1 + 1],
					 &c__1);
			    }
			    if (ilz) {
				dswap_(n, &z__[(*ranke + j) * z_dim1 + 1], &
					c__1, &z__[(*ranke + k) * z_dim1 + 1],
					 &c__1);
			    }
			    iwork[k] = -iwork[k];
			    j = k;
			    goto L50;
			}
		    }
/* L60: */
		}

		if (reda && *rnka22 < na22) {

/*                 Determine an orthogonal matrix Y2 such that */

/*                 ( R11 R12 ) = ( Ar  0 ) * Y2 . */

/*                 Compute A <-- A*diag(I, Y2'), C <-- C*diag(I, Y2'), */
/*                         Z <-- Z*diag(I, Y2'). */
/*                 Workspace: need   2*N. */
/*                            prefer N + N*NB. */

		    kw = *ranke + 1;
		    i__1 = *ldwork - kw + 1;
		    dtzrzf_(rnka22, &na22, &a[ir1 + ir1 * a_dim1], lda, &
			    dwork[1], &dwork[kw], &i__1, info);
/* Computing MAX */
		    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
		    wrkopt = max(i__1,i__2);

/*                 Workspace: need   N + MAX(P,N); */
/*                            prefer N + MAX(P,N)*NB. */

		    lh = na22 - *rnka22;
		    if (withc) {
			i__1 = *ldwork - kw + 1;
			dormrz_("Right", "Transpose", p, n, rnka22, &lh, &a[
				ir1 + ir1 * a_dim1], lda, &dwork[1], &c__[
				c_offset], ldc, &dwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)9);
/* Computing MAX */
			i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
			wrkopt = max(i__1,i__2);
		    }
		    if (ilz) {
			i__1 = *ldwork - kw + 1;
			dormrz_("Right", "Transpose", n, n, rnka22, &lh, &a[
				ir1 + ir1 * a_dim1], lda, &dwork[1], &z__[
				z_offset], ldz, &dwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)9);
/* Computing MAX */
			i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
			wrkopt = max(i__1,i__2);
		    }
		    ire1 = *ranke + *rnka22 + 1;

/*                 Set R12 and R22 to zero. */

		    dlaset_("Full", &la22, &lh, &c_b42, &c_b42, &a[ir1 + ire1 
			    * a_dim1], lda, (ftnlen)4);
		}
	    } else {
		dlaset_("Full", &la22, &na22, &c_b42, &c_b42, &a[ir1 + ir1 * 
			a_dim1], lda, (ftnlen)4);
	    }
	}
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of TG01FD *** */
} /* tg01fd_ */

