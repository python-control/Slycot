/* TG01CD.f -- translated by f2c (version 20100827).
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

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;

/* Subroutine */ int tg01cd_(char *compq, integer *l, integer *n, integer *m, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, doublereal *
	b, integer *ldb, doublereal *q, integer *ldq, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen compq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer ln;
    static logical ilq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer icompq;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

/*     To reduce the descriptor system pair (A-lambda E,B) to the */
/*     QR-coordinate form by computing an orthogonal transformation */
/*     matrix Q such that the transformed descriptor system pair */
/*     (Q'*A-lambda Q'*E, Q'*B) has the descriptor matrix Q'*E */
/*     in an upper trapezoidal form. */
/*     The left orthogonal transformations performed to reduce E */
/*     can be optionally accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPQ   CHARACTER*1 */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     orthogonal matrix Q is returned; */
/*             = 'U':  Q must contain an orthogonal matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A and E.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*E in upper trapezoidal form, */
/*             i.e. */

/*                      ( E11 ) */
/*               Q'*E = (     ) ,     if L >= N , */
/*                      (  0  ) */
/*             or */

/*               Q'*E = ( E11 E12 ),  if L < N , */

/*             where E11 is an MIN(L,N)-by-MIN(L,N) upper triangular */
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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, MIN(L,N) + MAX(L,N,M)). */
/*             For optimum performance */
/*             LWORK >= MAX(1, MIN(L,N) + MAX(L,N,M)*NB), */
/*             where NB is the optimal blocksize. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes the QR factorization of E to reduce it */
/*     to the upper trapezoidal form. */

/*     The transformations are also applied to the rest of system */
/*     matrices */

/*         A <- Q' * A ,  B <- Q' * B. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( L*L*N )  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara, University "Politehnica" Bucharest. */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSQR. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     May 2003. */

/*     KEYWORDS */

/*     Descriptor system, matrix algebra, matrix operations, */
/*     orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
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
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
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

/*     Test the input parameters. */

    *info = 0;
/* Computing MAX */
/* Computing MAX */
    i__3 = max(*l,*n);
    i__1 = 1, i__2 = min(*l,*n) + max(i__3,*m);
    wrkopt = max(i__1,i__2);
    if (icompq == 0) {
	*info = -1;
    } else if (*l < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*lda < max(1,*l)) {
	*info = -6;
    } else if (*lde < max(1,*l)) {
	*info = -8;
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
	*info = -10;
    } else if (ilq && *ldq < *l || *ldq < 1) {
	*info = -12;
    } else if (*ldwork < wrkopt) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TG01CD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize Q if necessary. */

    if (icompq == 3) {
	dlaset_("Full", l, l, &c_b7, &c_b8, &q[q_offset], ldq, (ftnlen)4);
    }

/*     Quick return if possible. */

    if (*l == 0 || *n == 0) {
	dwork[1] = 1.;
	return 0;
    }

    ln = min(*l,*n);

/*     Compute the QR decomposition of E. */

/*     Workspace: need   MIN(L,N) + N; */
/*                prefer MIN(L,N) + N*NB. */

    i__1 = *ldwork - ln;
    dgeqrf_(l, n, &e[e_offset], lde, &dwork[1], &dwork[ln + 1], &i__1, info);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
    wrkopt = max(i__1,i__2);

/*     Apply transformation on the rest of matrices. */

/*     A <-- Q' * A. */
/*     Workspace: need   MIN(L,N) + N; */
/*                prefer MIN(L,N) + N*NB. */

    i__1 = *ldwork - ln;
    dormqr_("Left", "Transpose", l, n, &ln, &e[e_offset], lde, &dwork[1], &a[
	    a_offset], lda, &dwork[ln + 1], &i__1, info, (ftnlen)4, (ftnlen)9)
	    ;
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
    wrkopt = max(i__1,i__2);

/*     B <-- Q' * B. */
/*     Workspace: need   MIN(L,N) + M; */
/*                prefer MIN(L,N) + M*NB. */

    if (*m > 0) {
	i__1 = *ldwork - ln;
	dormqr_("Left", "Transpose", l, m, &ln, &e[e_offset], lde, &dwork[1], 
		&b[b_offset], ldb, &dwork[ln + 1], &i__1, info, (ftnlen)4, (
		ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
	wrkopt = max(i__1,i__2);
    }

/*     Q <-- Q1 * Q. */
/*     Workspace: need   MIN(L,N) + L; */
/*                prefer MIN(L,N) + L*NB. */

    if (ilq) {
	i__1 = *ldwork - ln;
	dormqr_("Right", "No Transpose", l, l, &ln, &e[e_offset], lde, &dwork[
		1], &q[q_offset], ldq, &dwork[ln + 1], &i__1, info, (ftnlen)5,
		 (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
	wrkopt = max(i__1,i__2);
    }

/*     Set lower triangle of E to zero. */

    if (*l >= 2) {
	i__1 = *l - 1;
	dlaset_("Lower", &i__1, &ln, &c_b7, &c_b7, &e[e_dim1 + 2], lde, (
		ftnlen)5);
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of TG01CD *** */
} /* tg01cd_ */

