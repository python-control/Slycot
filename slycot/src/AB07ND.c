/* AB07ND.f -- translated by f2c (version 20100827).
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
static doublereal c_b15 = -1.;
static doublereal c_b16 = 0.;
static doublereal c_b31 = 1.;

/* Subroutine */ int ab07nd_(integer *n, integer *m, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *rcond, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, bl, ierr;
    static logical blas3;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical block;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal dnorm;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgetri_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *);
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

/*     To compute the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs and outputs.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state matrix A of the original system. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the state matrix Ai of the inverse system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the original system. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the input matrix Bi of the inverse system. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the output matrix C of the original system. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the output matrix Ci of the inverse system. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,M). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the feedthrough matrix D of the original system. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the feedthrough matrix Di of the inverse system. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,M). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The estimated reciprocal condition number of the */
/*             feedthrough matrix D of the original system. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or M+1, DWORK(1) returns the optimal */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,4*M). */
/*             For good performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  the matrix D is exactly singular; the (i,i) diagonal */
/*                   element is zero, i <= M; RCOND was set to zero; */
/*             = M+1:  the matrix D is numerically singular, i.e., RCOND */
/*                   is less than the relative machine precision, EPS */
/*                   (see LAPACK Library routine DLAMCH). The */
/*                   calculations have been completed, but the results */
/*                   could be very inaccurate. */

/*     METHOD */

/*     The matrices of the inverse system are computed with the formulas: */
/*                   -1              -1         -1           -1 */
/*       Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  . */

/*     NUMERICAL ASPECTS */

/*     The accuracy depends mainly on the condition number of the matrix */
/*     D to be inverted. The estimated reciprocal condition number is */
/*     returned in RCOND. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000. */
/*     D. Sima, University of Bucharest, April 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */
/*     Based on the routine SYSINV, A. Varga, 1992. */

/*     REVISIONS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000. */

/*     KEYWORDS */

/*     Inverse system, state-space model, state-space representation. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -6;
    } else if (*ldc < max(1,*m)) {
	*info = -8;
    } else if (*ldd < max(1,*m)) {
	*info = -10;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m << 2;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -14;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB07ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0) {
	*rcond = 1.;
	dwork[1] = 1.;
	return 0;
    }

/*     Factorize D. */

    dgetrf_(m, m, &d__[d_offset], ldd, &iwork[1], info);
    if (*info != 0) {
	*rcond = 0.;
	return 0;
    }

/*     Compute the reciprocal condition number of the matrix D. */
/*     Workspace: need   4*M. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*      minimal amount of workspace needed at that point in the code, */
/*      as well as the preferred amount for good performance. */
/*      NB refers to the optimal block size for the immediately */
/*      following subroutine, as returned by ILAENV.) */

    dnorm = dlange_("1-norm", m, m, &d__[d_offset], ldd, &dwork[1], (ftnlen)6)
	    ;
    dgecon_("1-norm", m, &d__[d_offset], ldd, &dnorm, rcond, &dwork[1], &
	    iwork[*m + 1], &ierr, (ftnlen)6);
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
	*info = *m + 1;
    }
/*                   -1 */
/*     Compute Di = D  . */
/*     Workspace: need   M; */
/*                prefer M*NB. */

/* Computing MAX */
    i__1 = *m << 2, i__2 = *m * ilaenv_(&c__1, "DGETRI", " ", m, &c_n1, &c_n1,
	     &c_n1, (ftnlen)6, (ftnlen)1);
    maxwrk = max(i__1,i__2);
    dgetri_(m, &d__[d_offset], ldd, &iwork[1], &dwork[1], ldwork, &ierr);
    if (*n > 0) {
	chunk = *ldwork / *m;
	blas3 = chunk >= *n && *m > 1;
	block = min(chunk,*m) > 1;
/*                          -1 */
/*        Compute  Bi = -B*D  . */

	if (blas3) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
	    dgemm_("NoTranspose", "NoTranspose", n, m, m, &c_b15, &dwork[1], 
		    n, &d__[d_offset], ldd, &c_b16, &b[b_offset], ldb, (
		    ftnlen)11, (ftnlen)11);

	} else if (block) {

/*           Use as many rows of B as possible. */

	    i__1 = *n;
	    i__2 = chunk;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
		i__3 = *n - i__ + 1;
		bl = min(i__3,chunk);
		dlacpy_("Full", &bl, m, &b[i__ + b_dim1], ldb, &dwork[1], &bl,
			 (ftnlen)4);
		dgemm_("NoTranspose", "NoTranspose", &bl, m, m, &c_b15, &
			dwork[1], &bl, &d__[d_offset], ldd, &c_b16, &b[i__ + 
			b_dim1], ldb, (ftnlen)11, (ftnlen)11);
/* L10: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dcopy_(m, &b[i__ + b_dim1], ldb, &dwork[1], &c__1);
		dgemv_("Transpose", m, m, &c_b15, &d__[d_offset], ldd, &dwork[
			1], &c__1, &c_b16, &b[i__ + b_dim1], ldb, (ftnlen)9);
/* L20: */
	    }

	}

/*        Compute  Ai = A + Bi*C. */

	dgemm_("NoTranspose", "NoTranspose", n, n, m, &c_b31, &b[b_offset], 
		ldb, &c__[c_offset], ldc, &c_b31, &a[a_offset], lda, (ftnlen)
		11, (ftnlen)11);
/*                        -1 */
/*        Compute  C <-- D  *C. */

	if (blas3) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (ftnlen)
		    4);
	    dgemm_("NoTranspose", "NoTranspose", m, n, m, &c_b31, &d__[
		    d_offset], ldd, &dwork[1], m, &c_b16, &c__[c_offset], ldc,
		     (ftnlen)11, (ftnlen)11);

	} else if (block) {

/*           Use as many columns of C as possible. */

	    i__2 = *n;
	    i__1 = chunk;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* Computing MIN */
		i__3 = *n - j + 1;
		bl = min(i__3,chunk);
		dlacpy_("Full", m, &bl, &c__[j * c_dim1 + 1], ldc, &dwork[1], 
			m, (ftnlen)4);
		dgemm_("NoTranspose", "NoTranspose", m, &bl, m, &c_b31, &d__[
			d_offset], ldd, &dwork[1], m, &c_b16, &c__[j * c_dim1 
			+ 1], ldc, (ftnlen)11, (ftnlen)11);
/* L30: */
	    }

	} else {

/*           Use a BLAS 2 algorithm. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &dwork[1], &c__1);
		dgemv_("NoTranspose", m, m, &c_b31, &d__[d_offset], ldd, &
			dwork[1], &c__1, &c_b16, &c__[j * c_dim1 + 1], &c__1, 
			(ftnlen)11);
/* L40: */
	    }

	}
    }

/*     Return optimal workspace in DWORK(1). */

/* Computing MAX */
    i__1 = maxwrk, i__2 = *n * *m;
    dwork[1] = (doublereal) max(i__1,i__2);
    return 0;

/* *** Last line of AB07ND *** */
} /* ab07nd_ */

