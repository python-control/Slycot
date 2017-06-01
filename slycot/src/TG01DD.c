/* TG01DD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tg01dd_(char *compz, integer *l, integer *n, integer *p, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, doublereal *
	c__, integer *ldc, doublereal *z__, integer *ldz, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, e_dim1, e_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer ln;
    static logical ilz;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgerqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer icompz;
    extern /* Subroutine */ int dormrq_(char *, char *, integer *, integer *, 
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

/*     To reduce the descriptor system pair (C,A-lambda E) to the */
/*     RQ-coordinate form by computing an orthogonal transformation */
/*     matrix Z such that the transformed descriptor system pair */
/*     (C*Z,A*Z-lambda E*Z) has the descriptor matrix E*Z in an upper */
/*     trapezoidal form. */
/*     The right orthogonal transformations performed to reduce E can */
/*     be optionally accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPZ   CHARACTER*1 */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     orthogonal matrix Z is returned; */
/*             = 'U':  Z must contain an orthogonal matrix Z1 on entry, */
/*                     and the product Z1*Z is returned. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix A*Z. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix E*Z in upper trapezoidal form, */
/*             i.e. */

/*                      ( E11 ) */
/*                E*Z = (     ) ,  if L >= N , */
/*                      (  R  ) */
/*             or */

/*                E*Z = ( 0  R ),  if L < N , */

/*             where R is an MIN(L,N)-by-MIN(L,N) upper triangular */
/*             matrix. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, MIN(L,N) + MAX(L,N,P)). */
/*             For optimum performance */
/*             LWORK >= MAX(1, MIN(L,N) + MAX(L,N,P)*NB), */
/*             where NB is the optimal blocksize. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes the RQ factorization of E to reduce it */
/*     the upper trapezoidal form. */

/*     The transformations are also applied to the rest of system */
/*     matrices */

/*         A <- A * Z,  C <- C * Z. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( L*N*N )  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara, University "Politehnica" Bucharest. */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSRQ. */

/*     REVISIONS */

/*     July 1999, V. Sima, Research Institute for Informatics, Bucharest. */

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

/*     Decode COMPZ. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --dwork;

    /* Function Body */
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

/*     Test the input parameters. */

    *info = 0;
/* Computing MAX */
/* Computing MAX */
    i__3 = max(*l,*n);
    i__1 = 1, i__2 = min(*l,*n) + max(i__3,*p);
    wrkopt = max(i__1,i__2);
    if (icompz == 0) {
	*info = -1;
    } else if (*l < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*lda < max(1,*l)) {
	*info = -6;
    } else if (*lde < max(1,*l)) {
	*info = -8;
    } else if (*ldc < max(1,*p)) {
	*info = -10;
    } else if (ilz && *ldz < *n || *ldz < 1) {
	*info = -12;
    } else if (*ldwork < wrkopt) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TG01DD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize Q if necessary. */

    if (icompz == 3) {
	dlaset_("Full", n, n, &c_b7, &c_b8, &z__[z_offset], ldz, (ftnlen)4);
    }

/*     Quick return if possible. */

    if (*l == 0 || *n == 0) {
	dwork[1] = 1.;
	return 0;
    }

    ln = min(*l,*n);

/*     Compute the RQ decomposition of E, E = R*Z. */

/*     Workspace: need   MIN(L,N) + L; */
/*                prefer MIN(L,N) + L*NB. */

    i__1 = *ldwork - ln;
    dgerqf_(l, n, &e[e_offset], lde, &dwork[1], &dwork[ln + 1], &i__1, info);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
    wrkopt = max(i__1,i__2);

/*     Apply transformation on the rest of matrices. */

/*     A <--  A * Z'. */
/*     Workspace: need   MIN(L,N) + L; */
/*                prefer MIN(L,N) + L*NB. */

    i__1 = *ldwork - ln;
    dormrq_("Right", "Transpose", l, n, &ln, &e[*l - ln + 1 + e_dim1], lde, &
	    dwork[1], &a[a_offset], lda, &dwork[ln + 1], &i__1, info, (ftnlen)
	    5, (ftnlen)9);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
    wrkopt = max(i__1,i__2);

/*     C <-- C * Z'. */
/*     Workspace: need   MIN(L,N) + P; */
/*                prefer MIN(L,N) + P*NB. */

    i__1 = *ldwork - ln;
    dormrq_("Right", "Transpose", p, n, &ln, &e[*l - ln + 1 + e_dim1], lde, &
	    dwork[1], &c__[c_offset], ldc, &dwork[ln + 1], &i__1, info, (
	    ftnlen)5, (ftnlen)9);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
    wrkopt = max(i__1,i__2);

/*     Z <-- Z1 * Z'. */
/*     Workspace: need   MIN(L,N) + N; */
/*                prefer MIN(L,N) + N*NB. */

    if (ilz) {
	i__1 = *ldwork - ln;
	dormrq_("Right", "Transpose", n, n, &ln, &e[*l - ln + 1 + e_dim1], 
		lde, &dwork[1], &z__[z_offset], ldz, &dwork[ln + 1], &i__1, 
		info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[ln + 1] + ln;
	wrkopt = max(i__1,i__2);
    }

/*     Set lower triangle of E to zero. */

    if (*l < *n) {
	i__1 = *n - *l;
	dlaset_("Full", l, &i__1, &c_b7, &c_b7, &e[e_offset], lde, (ftnlen)4);
	if (*l >= 2) {
	    i__1 = *l - 1;
	    dlaset_("Lower", &i__1, l, &c_b7, &c_b7, &e[(*n - *l + 1) * 
		    e_dim1 + 2], lde, (ftnlen)5);
	}
    } else {
	if (*n >= 2) {
	    i__1 = *n - 1;
	    dlaset_("Lower", &i__1, n, &c_b7, &c_b7, &e[*l - *n + 2 + e_dim1],
		     lde, (ftnlen)5);
	}
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of TG01DD *** */
} /* tg01dd_ */

