/* AB09DD.f -- translated by f2c (version 20100827).
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

static doublereal c_b14 = 1.;

/* Subroutine */ int ab09dd_(char *dico, integer *n, integer *m, integer *p, 
	integer *nr, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	 doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *rcond, integer *iwork, doublereal *dwork, integer *info, 
	ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, ns;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal a22nrm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen), dgetrs_(char *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);


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

/*     To compute a reduced order model by using singular perturbation */
/*     approximation formulas. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the state vector, i.e. the order of the */
/*             matrix A; also the number of rows of matrix B and the */
/*             number of columns of the matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector, i.e. the number of columns */
/*             of matrices B and D.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector, i.e. the number of rows of */
/*             matrices C and D.  P >= 0. */

/*     NR      (input) INTEGER */
/*             The order of the reduced order system.  N >= NR >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix of the original system. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the state dynamics matrix Ar of the reduced order system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix of the original system. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the input/state matrix Br of the reduced order system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix of the original system. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the state/output matrix Cr of the reduced order system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix of the original system. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the input/output matrix Dr of the reduced order system. */
/*             If NR = 0 and the given system is stable, then D contains */
/*             the steady state gain of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal condition number of the matrix A22-g*I */
/*             (see METHOD). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*(N-NR) */

/*     DWORK   DOUBLE PRECISION array, dimension 4*(N-NR) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1: if the matrix A22-g*I (see METHOD) is numerically */
/*                  singular. */

/*     METHOD */

/*     Given the system (A,B,C,D), partition the system matrices as */

/*            ( A11 A12 )        ( B1 ) */
/*        A = (         ) ,  B = (    ) ,  C = ( C1  C2 ), */
/*            ( A21 A22 )        ( B2 ) */

/*     where A11 is NR-by-NR, B1 is NR-by-M, C1 is P-by-NR, and the other */
/*     submatrices have appropriate dimensions. */

/*     The matrices of the reduced order system (Ar,Br,Cr,Dr) are */
/*     computed according to the following residualization formulas: */
/*                                -1                               -1 */
/*        Ar = A11 + A12*(g*I-A22)  *A21 ,  Br = B1 + A12*(g*I-A22)  *B2 */
/*                              -1                               -1 */
/*        Cr = C1 + C2*(g*I-A22)  *A21   ,  Dr = D + C2*(g*I-A22)  *B2 */

/*     where g = 0 if DICO = 'C' and g = 1 if DICO = 'D'. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRESID. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Model reduction, multivariable system, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the input scalar arguments. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*nr < 0 || *nr > *n) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*p)) {
	*info = -11;
    } else if (*ldd < max(1,*p)) {
	*info = -13;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09DD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*nr == *n) {
	*rcond = 1.;
	return 0;
    }

    k = *nr + 1;
    ns = *n - *nr;

/*     Compute: T = -A22   if  DICO = 'C' and */
/*              T = -A22+I if  DICO = 'D'. */

    i__1 = *n;
    for (j = k; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = k; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = -a[i__ + j * a_dim1];
/* L10: */
	}
	if (discr) {
	    a[j + j * a_dim1] += 1.;
	}
/* L20: */
    }

/*     Compute the LU decomposition of T. */

    a22nrm = dlange_("1-norm", &ns, &ns, &a[k + k * a_dim1], lda, &dwork[1], (
	    ftnlen)6);
    dgetrf_(&ns, &ns, &a[k + k * a_dim1], lda, &iwork[1], info);
    if (*info > 0) {

/*        Error return. */

	*rcond = 0.;
	*info = 1;
	return 0;
    }
    dgecon_("1-norm", &ns, &a[k + k * a_dim1], lda, &a22nrm, rcond, &dwork[1],
	     &iwork[ns + 1], info, (ftnlen)6);
    if (*rcond <= dlamch_("E", (ftnlen)1)) {

/*        Error return. */

	*info = 1;
	return 0;
    }

/*     Compute A21 <- INV(T)*A21. */

    dgetrs_("NoTranspose", &ns, nr, &a[k + k * a_dim1], lda, &iwork[1], &a[k 
	    + a_dim1], lda, info, (ftnlen)11);

/*     Compute B2 <- INV(T)*B2. */

    dgetrs_("NoTranspose", &ns, m, &a[k + k * a_dim1], lda, &iwork[1], &b[k + 
	    b_dim1], ldb, info, (ftnlen)11);

/*     Compute the residualized systems matrices. */
/*     Ar = A11 + A12*INV(T)*A21. */

    dgemm_("NoTranspose", "NoTranspose", nr, nr, &ns, &c_b14, &a[k * a_dim1 + 
	    1], lda, &a[k + a_dim1], lda, &c_b14, &a[a_offset], lda, (ftnlen)
	    11, (ftnlen)11);

/*     Br = B1 + A12*INV(T)*B2. */

    dgemm_("NoTranspose", "NoTranspose", nr, m, &ns, &c_b14, &a[k * a_dim1 + 
	    1], lda, &b[k + b_dim1], ldb, &c_b14, &b[b_offset], ldb, (ftnlen)
	    11, (ftnlen)11);

/*     Cr = C1 + C2*INV(T)*A21. */

    dgemm_("NoTranspose", "NoTranspose", p, nr, &ns, &c_b14, &c__[k * c_dim1 
	    + 1], ldc, &a[k + a_dim1], lda, &c_b14, &c__[c_offset], ldc, (
	    ftnlen)11, (ftnlen)11);

/*     Dr = D + C2*INV(T)*B2. */

    dgemm_("NoTranspose", "NoTranspose", p, m, &ns, &c_b14, &c__[k * c_dim1 + 
	    1], ldc, &b[k + b_dim1], ldb, &c_b14, &d__[d_offset], ldd, (
	    ftnlen)11, (ftnlen)11);

    return 0;
/* *** Last line of AB09DD *** */
} /* ab09dd_ */

