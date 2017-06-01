/* TB01MD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tb01md_(char *jobu, char *uplo, integer *n, integer *m, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	u, integer *ldu, doublereal *dwork, integer *info, ftnlen jobu_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, u_dim1, u_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer j, m1, n1, ii, nj;
    static doublereal dz;
    static integer par1, par2, par3, par4, par5, par6;
    static logical ljoba, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical luplo;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlatzm_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);


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

/*     To reduce the pair (B,A) to upper or lower controller Hessenberg */
/*     form using (and optionally accumulating) unitary state-space */
/*     transformations. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the unitary state-space transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is initialized to the unit matrix and the */
/*                     unitary transformation matrix U is returned; */
/*             = 'U':  The given matrix U is updated by the unitary */
/*                     transformations used in the reduction. */

/*     UPLO    CHARACTER*1 */
/*             Indicates whether the user wishes the pair (B,A) to be */
/*             reduced to upper or lower controller Hessenberg form as */
/*             follows: */
/*             = 'U':  Upper controller Hessenberg form; */
/*             = 'L':  Lower controller Hessenberg form. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the */
/*             matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e. the number of columns of */
/*             the matrix B.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A to be transformed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state transition matrix U' * A * U. */
/*             The annihilated elements are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B to be transformed. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix U' * B. */
/*             The annihilated elements are set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, if JOBU = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix U (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix U and the state-space transformation */
/*             matrix which reduces the given pair to controller */
/*             Hessenberg form. */
/*             On exit, if JOBU = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated unitary */
/*             similarity transformations which reduces the given pair */
/*             to controller Hessenberg form. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. If JOBU = 'U' or */
/*             JOBU = 'I', LDU >= MAX(1,N); if JOBU = 'N', LDU >= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(N,M-1)) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a unitary state-space transformation U, which */
/*     reduces the pair (B,A) to one of the following controller */
/*     Hessenberg forms: */

/*                    |*  . . .  *|*  . . . . . .  *| */
/*                    |   .      .|.               .| */
/*                    |     .    .|.               .| */
/*                    |       .  .|.               .| */
/*       [U'B|U'AU] = |          *|.               .| N */
/*                    |           |*               .| */
/*                    |           |   .            .| */
/*                    |           |     .          .| */
/*                    |           |       .        .| */
/*                    |           |         * . .  *| */
/*                         M               N */

/*     if UPLO = 'U', or */

/*                    |*  . . *         |           | */
/*                    |.        .       |           | */
/*                    |.          .     |           | */
/*                    |.            .   |           | */
/*       [U'AU|U'B] = |.               *|           | N */
/*                    |.               .|*          | */
/*                    |.               .|.  .       | */
/*                    |.               .|.    .     | */
/*                    |.               .|.      .   | */
/*                    |*  . . . . . .  *|*  . . .  *| */
/*                            N               M */
/*     if UPLO = 'L'. */

/*     IF M >= N, then the matrix U'B is trapezoidal and U'AU is full. */

/*     REFERENCES */

/*     [1] Van Dooren, P. and Verhaegen, M.H.G. */
/*         On the use of unitary state-space transformations. */
/*         In : Contemporary Mathematics on Linear Algebra and its Role */
/*         in Systems Theory, 47, AMS, Providence, 1985. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O((N + M) x N**2) operations and is */
/*     backward stable (see [1]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01AD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     February 1997. */

/*     KEYWORDS */

/*     Controllability, controller Hessenberg form, orthogonal */
/*     transformation, unitary transformation. */

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
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ljobi = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
    ljoba = ljobi || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljoba && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (! ljoba && *ldu < 1 || ljoba && *ldu < max(1,*n)) {
	*info = -10;
    }

    if (*info != 0) {

/*        Error return */

	i__1 = -(*info);
	xerbla_("TB01MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0) {
	return 0;
    }

    m1 = *m + 1;
    n1 = *n - 1;

    if (ljobi) {

/*        Initialize U to the identity matrix. */

	dlaset_("Full", n, n, &c_b9, &c_b10, &u[u_offset], ldu, (ftnlen)4);
    }

/*     Perform transformations involving both B and A. */

    i__1 = min(*m,n1);
    for (j = 1; j <= i__1; ++j) {
	nj = *n - j;
	if (luplo) {
	    par1 = j;
	    par2 = j;
	    par3 = j + 1;
	    par4 = *m;
	    par5 = *n;
	} else {
	    par1 = *m - j + 1;
	    par2 = nj + 1;
	    par3 = 1;
	    par4 = *m - j;
	    par5 = nj;
	}

	i__2 = nj + 1;
	dlarfg_(&i__2, &b[par2 + par1 * b_dim1], &b[par3 + par1 * b_dim1], &
		c__1, &dz);

/*        Update A. */

	i__2 = nj + 1;
	dlatzm_("Left", &i__2, n, &b[par3 + par1 * b_dim1], &c__1, &dz, &a[
		par2 + a_dim1], &a[par3 + a_dim1], lda, &dwork[1], (ftnlen)4);
	i__2 = nj + 1;
	dlatzm_("Right", n, &i__2, &b[par3 + par1 * b_dim1], &c__1, &dz, &a[
		par2 * a_dim1 + 1], &a[par3 * a_dim1 + 1], lda, &dwork[1], (
		ftnlen)5);

	if (ljoba) {

/*           Update U. */

	    i__2 = nj + 1;
	    dlatzm_("Right", n, &i__2, &b[par3 + par1 * b_dim1], &c__1, &dz, &
		    u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &dwork[
		    1], (ftnlen)5);
	}

	if (j != *m) {

/*           Update B */

	    i__2 = nj + 1;
	    i__3 = par4 - par3 + 1;
	    dlatzm_("Left", &i__2, &i__3, &b[par3 + par1 * b_dim1], &c__1, &
		    dz, &b[par2 + par3 * b_dim1], &b[par3 + par3 * b_dim1], 
		    ldb, &dwork[1], (ftnlen)4);
	}

	i__2 = par5;
	for (ii = par3; ii <= i__2; ++ii) {
	    b[ii + par1 * b_dim1] = 0.;
/* L10: */
	}

/* L20: */
    }

    i__1 = n1;
    for (j = m1; j <= i__1; ++j) {

/*        Perform next transformations only involving A. */

	nj = *n - j;
	if (luplo) {
	    par1 = j - *m;
	    par2 = j;
	    par3 = j + 1;
	    par4 = *n;
	    par5 = j - *m + 1;
	    par6 = *n;
	} else {
	    par1 = *n + m1 - j;
	    par2 = nj + 1;
	    par3 = 1;
	    par4 = nj;
	    par5 = 1;
	    par6 = *n + *m - j;
	}

	i__2 = nj + 1;
	dlarfg_(&i__2, &a[par2 + par1 * a_dim1], &a[par3 + par1 * a_dim1], &
		c__1, &dz);

/*        Update A. */

	i__2 = nj + 1;
	i__3 = par6 - par5 + 1;
	dlatzm_("Left", &i__2, &i__3, &a[par3 + par1 * a_dim1], &c__1, &dz, &
		a[par2 + par5 * a_dim1], &a[par3 + par5 * a_dim1], lda, &
		dwork[1], (ftnlen)4);
	i__2 = nj + 1;
	dlatzm_("Right", n, &i__2, &a[par3 + par1 * a_dim1], &c__1, &dz, &a[
		par2 * a_dim1 + 1], &a[par3 * a_dim1 + 1], lda, &dwork[1], (
		ftnlen)5);

	if (ljoba) {

/*           Update U. */

	    i__2 = nj + 1;
	    dlatzm_("Right", n, &i__2, &a[par3 + par1 * a_dim1], &c__1, &dz, &
		    u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &dwork[
		    1], (ftnlen)5);
	}

	i__2 = par4;
	for (ii = par3; ii <= i__2; ++ii) {
	    a[ii + par1 * a_dim1] = 0.;
/* L30: */
	}

/* L40: */
    }

    return 0;
/* *** Last line of TB01MD *** */
} /* tb01md_ */

