/* AG07BD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static doublereal c_b24 = -1.;

/* Subroutine */ int ag07bd_(char *jobe, integer *n, integer *m, doublereal *
	a, integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *ai, integer *ldai, doublereal *ei, integer *ldei, 
	doublereal *bi, integer *ldbi, doublereal *ci, integer *ldci, 
	doublereal *di, integer *lddi, integer *info, ftnlen jobe_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, ai_dim1, ai_offset, b_dim1, b_offset, bi_dim1, 
	    bi_offset, c_dim1, c_offset, ci_dim1, ci_offset, d_dim1, d_offset,
	     di_dim1, di_offset, e_dim1, e_offset, ei_dim1, ei_offset, i__1, 
	    i__2;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical unite;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute the inverse (Ai-lambda*Ei,Bi,Ci,Di) of a given */
/*     descriptor system (A-lambda*E,B,C,D). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether E is a general square or an identity */
/*             matrix as follows: */
/*             = 'G':  E is a general square matrix; */
/*             = 'I':  E is the identity matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the square matrices A and E; */
/*             also the number of rows of matrix B and the number of */
/*             columns of matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs and outputs, i.e., the number */
/*             of columns of matrices B and D and the number of rows of */
/*             matrices C and D.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the original system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             If JOBE = 'G', the leading N-by-N part of this array must */
/*             contain the descriptor matrix E of the original system. */
/*             If JOBE = 'I', then E is assumed to be the identity */
/*             matrix and is not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E. */
/*             LDE >= MAX(1,N), if JOBE = 'G'; */
/*             LDE >= 1,        if JOBE = 'I'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input matrix B of the original system. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading M-by-N part of this array must contain the */
/*             output matrix C of the original system. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,M). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading M-by-M part of this array must contain the */
/*             feedthrough matrix D of the original system. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,M). */

/*     AI      (output) DOUBLE PRECISION array, dimension (LDAI,N+M) */
/*             The leading (N+M)-by-(N+M) part of this array contains */
/*             the state matrix Ai of the inverse system. */
/*             If LDAI = LDA >= N+M, then AI and A can share the same */
/*             storage locations. */

/*     LDAI    INTEGER */
/*             The leading dimension of the array AI. */
/*             LDAI >= MAX(1,N+M). */

/*     EI      (output) DOUBLE PRECISION array, dimension (LDEI,N+M) */
/*             The leading (N+M)-by-(N+M) part of this array contains */
/*             the descriptor matrix Ei of the inverse system. */
/*             If LDEI = LDE >= N+M, then EI and E can share the same */
/*             storage locations. */

/*     LDEI    INTEGER */
/*             The leading dimension of the array EI. */
/*             LDEI >= MAX(1,N+M). */

/*     BI      (output) DOUBLE PRECISION array, dimension (LDBI,M) */
/*             The leading (N+M)-by-M part of this array contains */
/*             the input matrix Bi of the inverse system. */
/*             If LDBI = LDB >= N+M, then BI and B can share the same */
/*             storage locations. */

/*     LDBI    INTEGER */
/*             The leading dimension of the array BI. */
/*             LDBI >= MAX(1,N+M). */

/*     CI      (output) DOUBLE PRECISION array, dimension (LDCI,N+M) */
/*             The leading M-by-(N+M) part of this array contains */
/*             the output matrix Ci of the inverse system. */
/*             If LDCI = LDC, CI and C can share the same storage */
/*             locations. */

/*     LDCI    INTEGER */
/*             The leading dimension of the array CI.  LDCI >= MAX(1,M). */

/*     DI      (output) DOUBLE PRECISION array, dimension (LDDI,M) */
/*             The leading M-by-M part of this array contains */
/*             the feedthrough matrix Di = 0 of the inverse system. */
/*             DI and D can share the same storage locations. */

/*     LDDI    INTEGER */
/*             The leading dimension of the array DI.  LDDI >= MAX(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrices of the inverse system are computed with the formulas */

/*                ( E  0 )        ( A  B )         (  0 ) */
/*           Ei = (      ) , Ai = (      ) ,  Bi = (    ), */
/*                ( 0  0 )        ( C  D )         ( -I ) */

/*           Ci = ( 0  I ),  Di = 0. */

/*     FURTHER COMMENTS */

/*     The routine does not perform an invertibility test. This check can */
/*     be performed by using the SLICOT routines AB08NX or AG08BY. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     KEYWORDS */

/*     Descriptor system, inverse system, state-space representation. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    ai_dim1 = *ldai;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ei_dim1 = *ldei;
    ei_offset = 1 + ei_dim1;
    ei -= ei_offset;
    bi_dim1 = *ldbi;
    bi_offset = 1 + bi_dim1;
    bi -= bi_offset;
    ci_dim1 = *ldci;
    ci_offset = 1 + ci_dim1;
    ci -= ci_offset;
    di_dim1 = *lddi;
    di_offset = 1 + di_dim1;
    di -= di_offset;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    unite = lsame_(jobe, "I", (ftnlen)1, (ftnlen)1);
    if (! (lsame_(jobe, "G", (ftnlen)1, (ftnlen)1) || unite)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*lde < 1 || ! unite && *lde < *n) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*m)) {
	*info = -11;
    } else if (*ldd < max(1,*m)) {
	*info = -13;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n + *m;
	if (*ldai < max(i__1,i__2)) {
	    *info = -15;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 1, i__2 = *n + *m;
	    if (*ldei < max(i__1,i__2)) {
		*info = -17;
	    } else /* if(complicated condition) */ {
/* Computing MAX */
		i__1 = 1, i__2 = *n + *m;
		if (*ldbi < max(i__1,i__2)) {
		    *info = -19;
		} else if (*ldci < max(1,*m)) {
		    *info = -21;
		} else if (*lddi < max(1,*m)) {
		    *info = -23;
		}
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AG07BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0) {
	return 0;
    }

/*     Form Ai. */

    dlacpy_("Full", n, n, &a[a_offset], lda, &ai[ai_offset], ldai, (ftnlen)4);
    dlacpy_("Full", m, n, &c__[c_offset], ldc, &ai[*n + 1 + ai_dim1], ldai, (
	    ftnlen)4);
    dlacpy_("Full", n, m, &b[b_offset], ldb, &ai[(*n + 1) * ai_dim1 + 1], 
	    ldai, (ftnlen)4);
    dlacpy_("Full", m, m, &d__[d_offset], ldd, &ai[*n + 1 + (*n + 1) * 
	    ai_dim1], ldai, (ftnlen)4);

/*     Form Ei. */

    if (unite) {
	i__1 = *n + *m;
	dlaset_("Full", &i__1, n, &c_b10, &c_b11, &ei[ei_offset], ldei, (
		ftnlen)4);
    } else {
	dlacpy_("Full", n, n, &e[e_offset], lde, &ei[ei_offset], ldei, (
		ftnlen)4);
	dlaset_("Full", m, n, &c_b10, &c_b10, &ei[*n + 1 + ei_dim1], ldei, (
		ftnlen)4);
    }
    i__1 = *n + *m;
    dlaset_("Full", &i__1, m, &c_b10, &c_b10, &ei[(*n + 1) * ei_dim1 + 1], 
	    ldei, (ftnlen)4);

/*     Form Bi. */

    dlaset_("Full", n, m, &c_b10, &c_b10, &bi[bi_offset], ldbi, (ftnlen)4);
    dlaset_("Full", m, m, &c_b10, &c_b24, &bi[*n + 1 + bi_dim1], ldbi, (
	    ftnlen)4);

/*     Form Ci. */

    dlaset_("Full", m, n, &c_b10, &c_b10, &ci[ci_offset], ldci, (ftnlen)4);
    dlaset_("Full", m, m, &c_b10, &c_b11, &ci[(*n + 1) * ci_dim1 + 1], ldci, (
	    ftnlen)4);

/*     Set Di. */

    dlaset_("Full", m, m, &c_b10, &c_b10, &di[di_offset], lddi, (ftnlen)4);

    return 0;
/* *** Last line of AG07BD *** */
} /* ag07bd_ */

