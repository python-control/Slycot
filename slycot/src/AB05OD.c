/* AB05OD.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;
static doublereal c_b32 = 1.;

/* Subroutine */ int ab05od_(char *over, integer *n1, integer *m1, integer *
	p1, integer *n2, integer *m2, doublereal *alpha, doublereal *a1, 
	integer *lda1, doublereal *b1, integer *ldb1, doublereal *c1, integer 
	*ldc1, doublereal *d1, integer *ldd1, doublereal *a2, integer *lda2, 
	doublereal *b2, integer *ldb2, doublereal *c2, integer *ldc2, 
	doublereal *d2, integer *ldd2, integer *n, integer *m, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, integer *info, ftnlen over_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, 
	    b_offset, b1_dim1, b1_offset, b2_dim1, b2_offset, c_dim1, 
	    c_offset, c1_dim1, c1_offset, c2_dim1, c2_offset, d_dim1, 
	    d_offset, d1_dim1, d1_offset, d2_dim1, d2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lover;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
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

/*     To obtain the state-space model (A,B,C,D) for rowwise */
/*     concatenation (parallel inter-connection on outputs, with separate */
/*     inputs) of two systems, each given in state-space form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     OVER    CHARACTER*1 */
/*             Indicates whether the user wishes to overlap pairs of */
/*             arrays, as follows: */
/*             = 'N':  Do not overlap; */
/*             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B, */
/*                     C1 and C, and D1 and D, i.e. the same name is */
/*                     effectively used for each pair (for all pairs) */
/*                     in the routine call.  In this case, setting */
/*                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD */
/*                     will give maximum efficiency. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The number of state variables in the first system, i.e. */
/*             the order of the matrix A1.  N1 >= 0. */

/*     M1      (input) INTEGER */
/*             The number of input variables for the first system. */
/*             M1 >= 0. */

/*     P1      (input) INTEGER */
/*             The number of output variables from each system.  P1 >= 0. */

/*     N2      (input) INTEGER */
/*             The number of state variables in the second system, i.e. */
/*             the order of the matrix A2.  N2 >= 0. */

/*     M2      (input) INTEGER */
/*             The number of input variables for the second system. */
/*             M2 >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             A coefficient multiplying the transfer-function matrix */
/*             (or the output equation) of the second system. */

/*     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1) */
/*             The leading N1-by-N1 part of this array must contain the */
/*             state transition matrix A1 for the first system. */

/*     LDA1    INTEGER */
/*             The leading dimension of array A1.  LDA1 >= MAX(1,N1). */

/*     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1) */
/*             The leading N1-by-M1 part of this array must contain the */
/*             input/state matrix B1 for the first system. */

/*     LDB1    INTEGER */
/*             The leading dimension of array B1.  LDB1 >= MAX(1,N1). */

/*     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1) */
/*             The leading P1-by-N1 part of this array must contain the */
/*             state/output matrix C1 for the first system. */

/*     LDC1    INTEGER */
/*             The leading dimension of array C1. */
/*             LDC1 >= MAX(1,P1) if N1 > 0. */
/*             LDC1 >= 1 if N1 = 0. */

/*     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1) */
/*             The leading P1-by-M1 part of this array must contain the */
/*             input/output matrix D1 for the first system. */

/*     LDD1    INTEGER */
/*             The leading dimension of array D1.  LDD1 >= MAX(1,P1). */

/*     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2) */
/*             The leading N2-by-N2 part of this array must contain the */
/*             state transition matrix A2 for the second system. */

/*     LDA2    INTEGER */
/*             The leading dimension of array A2.  LDA2 >= MAX(1,N2). */

/*     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M2) */
/*             The leading N2-by-M2 part of this array must contain the */
/*             input/state matrix B2 for the second system. */

/*     LDB2    INTEGER */
/*             The leading dimension of array B2.  LDB2 >= MAX(1,N2). */

/*     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2) */
/*             The leading P1-by-N2 part of this array must contain the */
/*             state/output matrix C2 for the second system. */

/*     LDC2    INTEGER */
/*             The leading dimension of array C2. */
/*             LDC2 >= MAX(1,P1) if N2 > 0. */
/*             LDC2 >= 1 if N2 = 0. */

/*     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M2) */
/*             The leading P1-by-M2 part of this array must contain the */
/*             input/output matrix D2 for the second system. */

/*     LDD2    INTEGER */
/*             The leading dimension of array D2.  LDD2 >= MAX(1,P1). */

/*     N       (output) INTEGER */
/*             The number of state variables (N1 + N2) in the connected */
/*             system, i.e. the order of the matrix A, the number of rows */
/*             of B and the number of columns of C. */

/*     M       (output) INTEGER */
/*             The number of input variables (M1 + M2) for the connected */
/*             system, i.e. the number of columns of B and D. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2) */
/*             The leading N-by-N part of this array contains the state */
/*             transition matrix A for the connected system. */
/*             The array A can overlap A1 if OVER = 'O'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N1+N2). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M1+M2) */
/*             The leading N-by-M part of this array contains the */
/*             input/state matrix B for the connected system. */
/*             The array B can overlap B1 if OVER = 'O'. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N1+N2). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2) */
/*             The leading P1-by-N part of this array contains the */
/*             state/output matrix C for the connected system. */
/*             The array C can overlap C1 if OVER = 'O'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P1) if N1+N2 > 0. */
/*             LDC >= 1 if N1+N2 = 0. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M1+M2) */
/*             The leading P1-by-M part of this array contains the */
/*             input/output matrix D for the connected system. */
/*             The array D can overlap D1 if OVER = 'O'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     After rowwise concatenation (parallel inter-connection with */
/*     separate inputs) of the two systems, */

/*     X1'     = A1*X1 + B1*U */
/*     Y1      = C1*X1 + D1*U */

/*     X2'     = A2*X2 + B2*V */
/*     Y2      = C2*X2 + D2*V */

/*     (where  '  denotes differentiation with respect to time), */

/*     with the output equation for the second system multiplied by a */
/*     scalar alpha, the following state-space model will be obtained: */

/*     X'      = A*X + B*(U) */
/*                       (V) */

/*     Y       = C*X + D*(U) */
/*                       (V) */

/*     where matrix  A  has the form    ( A1   0  ), */
/*                                      ( 0    A2 ) */

/*           matrix  B  has the form    ( B1   0  ), */
/*                                      ( 0    B2 ) */

/*           matrix  C  has the form    ( C1   alpha*C2 ) and */

/*           matrix  D  has the form    ( D1   alpha*D2 ). */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1996. */
/*     Supersedes Release 2.0 routine AB05CD by C.J.Benson, Kingston */
/*     Polytechnic, United Kingdom, January 1982. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003, */
/*     Feb. 2004. */

/*     KEYWORDS */

/*     Continuous-time system, multivariable system, state-space model, */
/*     state-space representation. */

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
    a1_dim1 = *lda1;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    b1_dim1 = *ldb1;
    b1_offset = 1 + b1_dim1;
    b1 -= b1_offset;
    c1_dim1 = *ldc1;
    c1_offset = 1 + c1_dim1;
    c1 -= c1_offset;
    d1_dim1 = *ldd1;
    d1_offset = 1 + d1_dim1;
    d1 -= d1_offset;
    a2_dim1 = *lda2;
    a2_offset = 1 + a2_dim1;
    a2 -= a2_offset;
    b2_dim1 = *ldb2;
    b2_offset = 1 + b2_dim1;
    b2 -= b2_offset;
    c2_dim1 = *ldc2;
    c2_offset = 1 + c2_dim1;
    c2 -= c2_offset;
    d2_dim1 = *ldd2;
    d2_offset = 1 + d2_dim1;
    d2 -= d2_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;

    /* Function Body */
    lover = lsame_(over, "O", (ftnlen)1, (ftnlen)1);
    *n = *n1 + *n2;
    *m = *m1 + *m2;
    *info = 0;

/*     Test the input scalar arguments. */

    if (! lover && ! lsame_(over, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n1 < 0) {
	*info = -2;
    } else if (*m1 < 0) {
	*info = -3;
    } else if (*p1 < 0) {
	*info = -4;
    } else if (*n2 < 0) {
	*info = -5;
    } else if (*m2 < 0) {
	*info = -6;
    } else if (*lda1 < max(1,*n1)) {
	*info = -9;
    } else if (*ldb1 < max(1,*n1)) {
	*info = -11;
    } else if (*n1 > 0 && *ldc1 < max(1,*p1) || *n1 == 0 && *ldc1 < 1) {
	*info = -13;
    } else if (*ldd1 < max(1,*p1)) {
	*info = -15;
    } else if (*lda2 < max(1,*n2)) {
	*info = -17;
    } else if (*ldb2 < max(1,*n2)) {
	*info = -19;
    } else if (*n2 > 0 && *ldc2 < max(1,*p1) || *n2 == 0 && *ldc2 < 1) {
	*info = -21;
    } else if (*ldd2 < max(1,*p1)) {
	*info = -23;
    } else if (*lda < max(1,*n)) {
	*info = -27;
    } else if (*ldb < max(1,*n)) {
	*info = -29;
    } else if (*n > 0 && *ldc < max(1,*p1) || *n == 0 && *ldc < 1) {
	*info = -31;
    } else if (*ldd < max(1,*p1)) {
	*info = -33;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB05OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = *n, i__2 = min(*m,*p1);
    if (max(i__1,i__2) == 0) {
	return 0;
    }

/*     First form the matrix  A. */

    if (lover && *lda1 <= *lda) {
	if (*lda1 < *lda) {

	    for (j = *n1; j >= 1; --j) {
		for (i__ = *n1; i__ >= 1; --i__) {
		    a[i__ + j * a_dim1] = a1[i__ + j * a1_dim1];
/* L10: */
		}
/* L20: */
	    }

	}
    } else {
	dlacpy_("F", n1, n1, &a1[a1_offset], lda1, &a[a_offset], lda, (ftnlen)
		1);
    }

    if (*n2 > 0) {
	dlacpy_("F", n2, n2, &a2[a2_offset], lda2, &a[*n1 + 1 + (*n1 + 1) * 
		a_dim1], lda, (ftnlen)1);
	dlaset_("F", n1, n2, &c_b10, &c_b10, &a[(*n1 + 1) * a_dim1 + 1], lda, 
		(ftnlen)1);
	dlaset_("F", n2, n1, &c_b10, &c_b10, &a[*n1 + 1 + a_dim1], lda, (
		ftnlen)1);
    }

/*     Now form the matrix  B. */

    if (lover && *ldb1 <= *ldb) {
	if (*ldb1 < *ldb) {

	    for (j = *m1; j >= 1; --j) {
		for (i__ = *n1; i__ >= 1; --i__) {
		    b[i__ + j * b_dim1] = b1[i__ + j * b1_dim1];
/* L30: */
		}
/* L40: */
	    }

	}
    } else {
	dlacpy_("F", n1, m1, &b1[b1_offset], ldb1, &b[b_offset], ldb, (ftnlen)
		1);
    }

    if (*m2 > 0) {
	if (*n2 > 0) {
	    dlacpy_("F", n2, m2, &b2[b2_offset], ldb2, &b[*n1 + 1 + (*m1 + 1) 
		    * b_dim1], ldb, (ftnlen)1);
	}
	dlaset_("F", n1, m2, &c_b10, &c_b10, &b[(*m1 + 1) * b_dim1 + 1], ldb, 
		(ftnlen)1);
    }
    if (*n2 > 0) {
	dlaset_("F", n2, m1, &c_b10, &c_b10, &b[*n1 + 1 + b_dim1], ldb, (
		ftnlen)1);
    }

/*     Now form the matrix  C. */

    if (lover && *ldc1 <= *ldc) {
	if (*ldc1 < *ldc) {

	    for (j = *n1; j >= 1; --j) {
		for (i__ = *p1; i__ >= 1; --i__) {
		    c__[i__ + j * c_dim1] = c1[i__ + j * c1_dim1];
/* L50: */
		}
/* L60: */
	    }

	}
    } else {
	dlacpy_("F", p1, n1, &c1[c1_offset], ldc1, &c__[c_offset], ldc, (
		ftnlen)1);
    }

    if (*n2 > 0) {
	dlacpy_("F", p1, n2, &c2[c2_offset], ldc2, &c__[(*n1 + 1) * c_dim1 + 
		1], ldc, (ftnlen)1);
	if (*alpha != 1.) {
	    dlascl_("G", &c__0, &c__0, &c_b32, alpha, p1, n2, &c__[(*n1 + 1) *
		     c_dim1 + 1], ldc, info, (ftnlen)1);
	}
    }

/*     Now form the matrix  D. */

    if (lover && *ldd1 <= *ldd) {
	if (*ldd1 < *ldd) {

	    for (j = *m1; j >= 1; --j) {
		for (i__ = *p1; i__ >= 1; --i__) {
		    d__[i__ + j * d_dim1] = d1[i__ + j * d1_dim1];
/* L70: */
		}
/* L80: */
	    }

	}
    } else {
	dlacpy_("F", p1, m1, &d1[d1_offset], ldd1, &d__[d_offset], ldd, (
		ftnlen)1);
    }

    if (*m2 > 0) {
	dlacpy_("F", p1, m2, &d2[d2_offset], ldd2, &d__[(*m1 + 1) * d_dim1 + 
		1], ldd, (ftnlen)1);
	if (*alpha != 1.) {
	    dlascl_("G", &c__0, &c__0, &c_b32, alpha, p1, m2, &d__[(*m1 + 1) *
		     d_dim1 + 1], ldd, info, (ftnlen)1);
	}
    }

    return 0;
/* *** Last line of AB05OD *** */
} /* ab05od_ */

